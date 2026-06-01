use std::io::BufReader;

use criterion::{Criterion, criterion_group, criterion_main};
use rsomics_pwmantel::{Alternative, LabeledMatrix, Method, pwmantel};

fn fixture() -> Vec<LabeledMatrix> {
    let dir = std::env::var("PWMANTEL_BENCH_DIR")
        .unwrap_or_else(|_| concat!(env!("CARGO_MANIFEST_DIR"), "/tests/golden").to_string());
    let mut dms = Vec::new();
    for entry in std::fs::read_dir(&dir).unwrap() {
        let path = entry.unwrap().path();
        let name = path.file_name().unwrap().to_string_lossy().into_owned();
        if name.starts_with("dm") && name.ends_with(".tsv") {
            let f = std::fs::File::open(&path).unwrap();
            dms.push(LabeledMatrix::read(BufReader::new(f), name).unwrap());
        }
    }
    dms.sort_by(|a, b| a.label.cmp(&b.label));
    dms
}

fn bench(c: &mut Criterion) {
    let dms = fixture();
    c.bench_function("pwmantel_pearson_999", |b| {
        b.iter(|| pwmantel(&dms, Method::Pearson, 999, Alternative::TwoSided, 42))
    });
}

criterion_group!(benches, bench);
criterion_main!(benches);
