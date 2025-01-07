use criterion::{criterion_group, criterion_main, Criterion};
use hgvs::sequences::{translate_cds, TranslationTable};
use std::sync::LazyLock;

/// TTN FASTA string from https://www.ncbi.nlm.nih.gov/nuccore/NM_001126114.1
static TTN_FASTA: &str = include_str!("TTN.fasta");

/// Raw TTN sequence.
static SEQ_TTN: LazyLock<String> = LazyLock::new(|| {
    let mut seq = String::new();
    for line in TTN_FASTA.lines() {
        if !line.starts_with('>') {
            seq.push_str(line);
        }
    }
    seq
});

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("translate_cds TTN", |b| {
        b.iter(|| translate_cds(&SEQ_TTN, true, "*", TranslationTable::Standard).unwrap())
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
