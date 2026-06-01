use std::collections::HashMap;
use std::io::BufReader;
use std::process::Command;

use rsomics_pwmantel::{Alternative, LabeledMatrix, Method, pwmantel};

const GOLDEN: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/golden");

fn load(name: &str) -> LabeledMatrix {
    let path = format!("{GOLDEN}/{name}");
    let f = std::fs::File::open(&path).unwrap();
    let label = std::path::Path::new(name)
        .file_stem()
        .unwrap()
        .to_string_lossy()
        .into_owned();
    LabeledMatrix::read(BufReader::new(f), label).unwrap()
}

fn ours(method: Method) -> HashMap<(String, String), f64> {
    let dms = vec![load("dm1.tsv"), load("dm2.tsv"), load("dm3.tsv")];
    let res = pwmantel(&dms, method, 0, Alternative::TwoSided, 1).unwrap();
    res.into_iter()
        .map(|r| ((r.dm1, r.dm2), r.statistic))
        .collect()
}

/// Always runs: every pairwise statistic must match the committed skbio
/// pwmantel capture to ~1e-9, for both pearson and spearman.
#[test]
fn statistic_matches_skbio_golden() {
    let pearson = ours(Method::Pearson);
    let spearman = ours(Method::Spearman);
    let golden = std::fs::read_to_string(format!("{GOLDEN}/golden.tsv")).unwrap();
    for line in golden.lines().skip(1) {
        let f: Vec<&str> = line.split('\t').collect();
        let method = Method::parse(f[2]).unwrap();
        let want: f64 = f[3].parse().unwrap();
        let table = match method {
            Method::Pearson => &pearson,
            Method::Spearman => &spearman,
        };
        let got = table[&(f[0].to_string(), f[1].to_string())];
        assert!(
            (got - want).abs() < 1e-9,
            "{} {}-{}: ours {got} vs skbio golden {want}",
            method.name(),
            f[0],
            f[1]
        );
    }
}

fn skbio_python() -> Option<String> {
    let candidates = [
        std::env::var("SKBIO_PYTHON").ok(),
        Some(format!(
            "{}/oracle-venvs/skbio/bin/python",
            std::env::var("HOME").unwrap_or_default()
        )),
    ];
    candidates.into_iter().flatten().find(|p| {
        Command::new(p)
            .args(["-c", "import skbio.stats.distance"])
            .status()
            .map(|s| s.success())
            .unwrap_or(false)
    })
}

/// Live differential against scikit-bio pwmantel. Loud-skips when the venv is
/// absent so CI stays green via the committed golden.
#[test]
fn live_skbio_pwmantel() {
    let Some(py) = skbio_python() else {
        eprintln!("SKIP live_skbio_pwmantel: scikit-bio venv not found");
        return;
    };

    let dir = std::env::temp_dir().join("rsomics-pwmantel-compat");
    std::fs::create_dir_all(&dir).unwrap();
    let script = dir.join("oracle.py");
    std::fs::write(
        &script,
        format!(
            r#"
from skbio import DistanceMatrix
from skbio.stats.distance import pwmantel
paths = ["{GOLDEN}/dm1.tsv","{GOLDEN}/dm2.tsv","{GOLDEN}/dm3.tsv"]
labels = ["dm1","dm2","dm3"]
dms = [DistanceMatrix.read(p) for p in paths]
for method in ("pearson","spearman"):
    res = pwmantel(dms, labels=labels, method=method, permutations=0)
    for (l1,l2), row in res.iterrows():
        print(f"{{method}}\t{{l1}}\t{{l2}}\t{{float(row['statistic'])!r}}")
"#
        ),
    )
    .unwrap();

    let out = Command::new(&py).arg(&script).output().unwrap();
    assert!(
        out.status.success(),
        "skbio oracle failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let text = String::from_utf8(out.stdout).unwrap();
    for line in text.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        let method = Method::parse(f[0]).unwrap();
        let want: f64 = f[3].trim().parse().unwrap();
        let got = ours(method)[&(f[1].to_string(), f[2].to_string())];
        assert!(
            (got - want).abs() < 1e-9,
            "{} {}-{}: ours {got} vs live skbio {want}",
            method.name(),
            f[1],
            f[2]
        );
    }
}

/// The permutation p-value is a Monte-Carlo estimate: the strongly correlated
/// dm1-dm2 pair must land significant with 999 permutations.
#[test]
fn p_value_in_expected_range() {
    let dms = vec![load("dm1.tsv"), load("dm2.tsv"), load("dm3.tsv")];
    let res = pwmantel(&dms, Method::Pearson, 999, Alternative::Greater, 42).unwrap();
    let strong = res
        .iter()
        .find(|r| r.dm1 == "dm1" && r.dm2 == "dm2")
        .unwrap();
    assert!(
        strong.p_value <= 0.05,
        "strong correlation should be significant, got p={}",
        strong.p_value
    );
}
