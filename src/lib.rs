use std::io::{BufRead, Write};

use rsomics_common::{Result, RsomicsError};

pub mod dm;
pub mod mantel;
mod rng;

pub use dm::DistanceMatrix;
pub use mantel::{Alternative, MantelResult, Method, mantel};

/// A distance matrix paired with the label it carries into the results table.
pub struct LabeledMatrix {
    pub label: String,
    pub matrix: DistanceMatrix,
    pub n: usize,
}

impl LabeledMatrix {
    pub fn read<R: BufRead>(reader: R, label: String) -> Result<Self> {
        let matrix = DistanceMatrix::read(reader, &label)?;
        let n = matrix.n;
        Ok(LabeledMatrix { label, matrix, n })
    }
}

pub struct PairResult {
    pub dm1: String,
    pub dm2: String,
    pub statistic: f64,
    pub p_value: f64,
    pub n: usize,
}

/// One Mantel test for every (i<j) pair of `dms`, in input order
/// (`itertools.combinations` semantics, i outer / j inner). Each pair gets an
/// independent permutation draw by folding its index into the base seed,
/// mirroring how skbio threads one rng through the pairwise loop.
pub fn pwmantel(
    dms: &[LabeledMatrix],
    method: Method,
    permutations: usize,
    alternative: Alternative,
    seed: u64,
) -> Result<Vec<PairResult>> {
    if dms.len() < 2 {
        return Err(RsomicsError::InvalidInput(
            "need at least two distance matrices".into(),
        ));
    }

    let mut out = Vec::with_capacity(dms.len() * (dms.len() - 1) / 2);
    let mut pair = 0u64;
    for i in 0..dms.len() {
        for j in (i + 1)..dms.len() {
            let (x, y) = (&dms[i], &dms[j]);
            if x.n != y.n {
                return Err(RsomicsError::InvalidInput(format!(
                    "{} ({}x{}) and {} ({}x{}) differ in size",
                    x.label, x.n, x.n, y.label, y.n, y.n
                )));
            }
            let y_data = y.matrix.reorder_like(&x.matrix.ids, &y.label)?;
            let res = mantel(
                &x.matrix.data,
                &y_data,
                x.n,
                method,
                permutations,
                alternative,
                seed.wrapping_add(pair.wrapping_mul(0x9E37_79B9_7F4A_7C15)),
            );
            out.push(PairResult {
                dm1: x.label.clone(),
                dm2: y.label.clone(),
                statistic: res.r,
                p_value: res.p_value,
                n: res.n,
            });
            pair += 1;
        }
    }
    Ok(out)
}

pub fn write_results<W: Write>(
    out: &mut W,
    results: &[PairResult],
    method: Method,
    permutations: usize,
    alternative: Alternative,
) -> Result<()> {
    writeln!(
        out,
        "dm1\tdm2\tstatistic\tp-value\tn\tmethod\tpermutations\talternative"
    )
    .map_err(RsomicsError::Io)?;
    for r in results {
        writeln!(
            out,
            "{}\t{}\t{:.12}\t{}\t{}\t{}\t{}\t{}",
            r.dm1,
            r.dm2,
            r.statistic,
            fmt_p(r.p_value),
            r.n,
            method.name(),
            permutations,
            alternative.name(),
        )
        .map_err(RsomicsError::Io)?;
    }
    Ok(())
}

fn fmt_p(p: f64) -> String {
    if p.is_nan() {
        "nan".to_string()
    } else {
        format!("{p:.12}")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn matrix(label: &str, body: &str) -> LabeledMatrix {
        LabeledMatrix::read(Cursor::new(body.to_string()), label.to_string()).unwrap()
    }

    fn doctest_dms() -> Vec<LabeledMatrix> {
        vec![
            matrix("x", "\ta\tb\tc\na\t0\t1\t2\nb\t1\t0\t3\nc\t2\t3\t0\n"),
            matrix("y", "\ta\tb\tc\na\t0\t2\t7\nb\t2\t0\t6\nc\t7\t6\t0\n"),
            matrix("z", "\ta\tb\tc\na\t0\t5\t6\nb\t5\t0\t1\nc\t6\t1\t0\n"),
        ]
    }

    #[test]
    fn pairwise_matches_skbio_doctest() {
        let dms = doctest_dms();
        let res = pwmantel(&dms, Method::Pearson, 0, Alternative::TwoSided, 1).unwrap();
        assert_eq!(res.len(), 3);
        let labels: Vec<_> = res
            .iter()
            .map(|r| (r.dm1.as_str(), r.dm2.as_str()))
            .collect();
        assert_eq!(labels, vec![("x", "y"), ("x", "z"), ("y", "z")]);
        assert!((res[0].statistic - 0.7559289460184544).abs() < 1e-12);
        assert!((res[1].statistic + 0.7559289460184544).abs() < 1e-12);
        assert!((res[2].statistic + 0.14285714285714285).abs() < 1e-12);
        assert!(res.iter().all(|r| r.p_value.is_nan()));
    }

    #[test]
    fn single_matrix_is_an_error() {
        let dms = vec![doctest_dms().pop().unwrap()];
        assert!(pwmantel(&dms, Method::Pearson, 0, Alternative::TwoSided, 1).is_err());
    }
}
