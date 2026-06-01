use rayon::prelude::*;

use crate::dm::DistanceMatrix;
use crate::rng;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Method {
    Pearson,
    Spearman,
}

impl Method {
    pub fn parse(s: &str) -> Option<Self> {
        match s {
            "pearson" => Some(Method::Pearson),
            "spearman" => Some(Method::Spearman),
            _ => None,
        }
    }
    pub fn name(self) -> &'static str {
        match self {
            Method::Pearson => "pearson",
            Method::Spearman => "spearman",
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Alternative {
    TwoSided,
    Greater,
    Less,
}

impl Alternative {
    pub fn parse(s: &str) -> Option<Self> {
        match s {
            "two-sided" => Some(Alternative::TwoSided),
            "greater" => Some(Alternative::Greater),
            "less" => Some(Alternative::Less),
            _ => None,
        }
    }
    pub fn name(self) -> &'static str {
        match self {
            Alternative::TwoSided => "two-sided",
            Alternative::Greater => "greater",
            Alternative::Less => "less",
        }
    }
}

pub struct MantelResult {
    pub r: f64,
    pub p_value: f64,
    pub n: usize,
}

/// One Mantel test. `x_data`/`y_data` are full square matrices on a shared id
/// order. The statistic reproduces skbio's `mantel()` to floating-point
/// tolerance; the p-value is a seeded permutation estimate.
pub fn mantel(
    x_data: &[f64],
    y_data: &[f64],
    n: usize,
    method: Method,
    permutations: usize,
    alternative: Alternative,
    seed: u64,
) -> MantelResult {
    let (x_flat, y_flat) = match method {
        Method::Pearson => (
            DistanceMatrix::condensed(x_data, n),
            DistanceMatrix::condensed(y_data, n),
        ),
        Method::Spearman => (
            rankdata(&DistanceMatrix::condensed(x_data, n)),
            rankdata(&DistanceMatrix::condensed(y_data, n)),
        ),
    };

    // Permutation acts on the full matrix; for Spearman that is the
    // rank-transformed matrix, rebuilt from the ranked condensed form.
    let x_full = match method {
        Method::Pearson => x_data.to_vec(),
        Method::Spearman => square_from_condensed(&x_flat, n),
    };

    let xmean = mean(&x_flat);
    let normx = norm_centered(&x_flat, xmean);
    let ym = normalize(&y_flat);
    let r = match (&ym, normx) {
        (Some(ymn), Some(nx)) => dot_centered(&x_flat, nx, ymn).clamp(-1.0, 1.0),
        _ => f64::NAN,
    };

    let p_value = if permutations == 0 || r.is_nan() {
        f64::NAN
    } else {
        let ymn = ym.unwrap();
        let nx = normx.unwrap();
        let count_extreme: usize = (0..permutations)
            .into_par_iter()
            .map(|k| {
                let perm = rng::permutation(n, seed, k as u64);
                let stat = permuted_stat(&x_full, n, &perm, nx, &ymn).clamp(-1.0, 1.0);
                match alternative {
                    Alternative::TwoSided => usize::from(stat.abs() >= r.abs()),
                    Alternative::Greater => usize::from(stat >= r),
                    Alternative::Less => usize::from(stat <= r),
                }
            })
            .sum();
        (count_extreme + 1) as f64 / (permutations + 1) as f64
    };

    MantelResult { r, p_value, n }
}

fn normalize(v: &[f64]) -> Option<Vec<f64>> {
    let m = mean(v);
    let mut out: Vec<f64> = v.iter().map(|&x| x - m).collect();
    let norm = out.iter().map(|&x| x * x).sum::<f64>().sqrt();
    if norm == 0.0 {
        return None;
    }
    for x in &mut out {
        *x /= norm;
    }
    Some(out)
}

fn mean(v: &[f64]) -> f64 {
    v.iter().sum::<f64>() / v.len() as f64
}

fn norm_centered(v: &[f64], m: f64) -> Option<f64> {
    let s = v.iter().map(|&x| (x - m) * (x - m)).sum::<f64>().sqrt();
    (s != 0.0).then_some(s)
}

fn dot_centered(x: &[f64], normx: f64, ym_norm: &[f64]) -> f64 {
    x.iter().zip(ym_norm).map(|(&xv, &yv)| xv * yv).sum::<f64>() / normx
}

/// ym_norm is mean-centered, so it sums to zero and the xmean term drops out of
/// the dot product. The inner loop is then a plain gather-multiply over the
/// upper triangle.
fn permuted_stat(x_full: &[f64], n: usize, perm: &[usize], normx: f64, ym_norm: &[f64]) -> f64 {
    let mut acc = 0.0;
    let mut k = 0;
    for i in 0..n {
        let base = perm[i] * n;
        for j in (i + 1)..n {
            acc += x_full[base + perm[j]] * ym_norm[k];
            k += 1;
        }
    }
    acc / normx
}

/// Average-rank of each element, scipy `rankdata` default (ties averaged).
fn rankdata(v: &[f64]) -> Vec<f64> {
    let mut order: Vec<usize> = (0..v.len()).collect();
    order.sort_by(|&a, &b| v[a].partial_cmp(&v[b]).unwrap());
    let mut ranks = vec![0.0f64; v.len()];
    let mut i = 0;
    while i < order.len() {
        let mut j = i + 1;
        while j < order.len() && v[order[j]] == v[order[i]] {
            j += 1;
        }
        let avg = ((i + 1 + j) as f64) / 2.0;
        for &idx in &order[i..j] {
            ranks[idx] = avg;
        }
        i = j;
    }
    ranks
}

fn square_from_condensed(cond: &[f64], n: usize) -> Vec<f64> {
    let mut out = vec![0.0f64; n * n];
    let mut k = 0;
    for i in 0..n {
        for j in (i + 1)..n {
            out[i * n + j] = cond[k];
            out[j * n + i] = cond[k];
            k += 1;
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn square(rows: &[&[f64]]) -> (Vec<f64>, usize) {
        let n = rows.len();
        let mut d = vec![0.0; n * n];
        for (i, r) in rows.iter().enumerate() {
            for (j, &v) in r.iter().enumerate() {
                d[i * n + j] = v;
            }
        }
        (d, n)
    }

    #[test]
    fn skbio_doc_example_pearson() {
        let (x, n) = square(&[&[0.0, 1.0, 2.0], &[1.0, 0.0, 3.0], &[2.0, 3.0, 0.0]]);
        let (y, _) = square(&[&[0.0, 2.0, 7.0], &[2.0, 0.0, 6.0], &[7.0, 6.0, 0.0]]);
        let res = mantel(&x, &y, n, Method::Pearson, 0, Alternative::TwoSided, 1);
        assert!((res.r - 0.7559289460184544).abs() < 1e-12, "r={}", res.r);
        assert!(res.p_value.is_nan());
    }

    #[test]
    fn rankdata_ties_averaged() {
        assert_eq!(rankdata(&[1.0, 2.0, 2.0, 3.0]), vec![1.0, 2.5, 2.5, 4.0]);
    }
}
