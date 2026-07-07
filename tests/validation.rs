use std::io::Write;
use std::process::Command;

/// A malformed distance matrix must be rejected before any computation, exactly
/// as skbio's `DistanceMatrix` constructor rejects it. The parser's earlier
/// checks (raggedness, labels, numerics) pass these, so without the extra
/// validation each would silently yield a confident wrong correlation.
fn write(dir: &std::path::Path, name: &str, body: &str) -> std::path::PathBuf {
    let p = dir.join(name);
    let mut f = std::fs::File::create(&p).unwrap();
    f.write_all(body.as_bytes()).unwrap();
    p
}

const VALID: &str = "\ta\tb\tc\na\t0.0\t1.0\t2.0\nb\t1.0\t0.0\t3.0\nc\t2.0\t3.0\t0.0\n";

fn run(a: &std::path::Path, b: &std::path::Path) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_rsomics-pwmantel"))
        .args([a, b])
        .args(["-p", "0"])
        .output()
        .unwrap()
}

#[test]
fn rejects_asymmetric_matrix() {
    let dir = tempfile::tempdir().unwrap();
    let valid = write(dir.path(), "valid.tsv", VALID);
    let asym = write(
        dir.path(),
        "asym.tsv",
        "\ta\tb\tc\na\t0.0\t1.0\t2.0\nb\t9.0\t0.0\t3.0\nc\t2.0\t3.0\t0.0\n",
    );
    let out = run(&valid, &asym);
    assert!(
        !out.status.success(),
        "asymmetric matrix must exit non-zero"
    );
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(err.contains("symmetric"), "stderr was: {err}");
}

#[test]
fn rejects_nonhollow_matrix() {
    let dir = tempfile::tempdir().unwrap();
    let valid = write(dir.path(), "valid.tsv", VALID);
    let nonhollow = write(
        dir.path(),
        "nonhollow.tsv",
        "\ta\tb\tc\na\t0.0\t1.0\t2.0\nb\t1.0\t5.0\t3.0\nc\t2.0\t3.0\t0.0\n",
    );
    let out = run(&valid, &nonhollow);
    assert!(
        !out.status.success(),
        "non-hollow matrix must exit non-zero"
    );
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(err.contains("diagonal"), "stderr was: {err}");
}

#[test]
fn rejects_nan_matrix() {
    let dir = tempfile::tempdir().unwrap();
    let valid = write(dir.path(), "valid.tsv", VALID);
    let nan = write(
        dir.path(),
        "nan.tsv",
        "\ta\tb\tc\na\t0.0\tnan\t2.0\nb\tnan\t0.0\t3.0\nc\t2.0\t3.0\t0.0\n",
    );
    let out = run(&valid, &nan);
    assert!(!out.status.success(), "NaN matrix must exit non-zero");
}

#[test]
fn rejects_duplicate_ids() {
    let dir = tempfile::tempdir().unwrap();
    let valid = write(dir.path(), "valid.tsv", VALID);
    let dup = write(
        dir.path(),
        "dup.tsv",
        "\ta\ta\tc\na\t0.0\t1.0\t2.0\na\t1.0\t0.0\t3.0\nc\t2.0\t3.0\t0.0\n",
    );
    let out = run(&valid, &dup);
    assert!(!out.status.success(), "duplicate ids must exit non-zero");
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(err.contains("unique"), "stderr was: {err}");
}

#[test]
fn negatives_are_allowed() {
    let dir = tempfile::tempdir().unwrap();
    let a = write(dir.path(), "a.tsv", VALID);
    let neg = write(
        dir.path(),
        "neg.tsv",
        "\ta\tb\tc\na\t0.0\t-1.0\t2.0\nb\t-1.0\t0.0\t-3.0\nc\t2.0\t-3.0\t0.0\n",
    );
    let out = run(&a, &neg);
    assert!(
        out.status.success(),
        "negative distances are valid: {}",
        String::from_utf8_lossy(&out.stderr)
    );
}
