[![Crates.io](https://img.shields.io/crates/d/hgvs.svg)](https://crates.io/crates/hgvs)
[![Crates.io](https://img.shields.io/crates/v/hgvs.svg)](https://crates.io/crates/hgvs)
[![Crates.io](https://img.shields.io/crates/l/hgvs.svg)](https://crates.io/crates/hgvs)
[![CI](https://github.com/bihealth/hgvs-rs/actions/workflows/rust.yml/badge.svg)](https://github.com/bihealth/hgvs-rs/actions/workflows/rust.yml)
[![codecov](https://codecov.io/gh/bihealth/hgvs-rs/branch/main/graph/badge.svg?token=aZchhLWdzt)](https://codecov.io/gh/bihealth/hgvs-rs)
[![DOI](https://zenodo.org/badge/601272076.svg)](https://zenodo.org/badge/latestdoi/601272076)

# hgvs-rs

This is a port of [biocommons/hgvs](https://github.com/biocommons/hgvs) to the Rust programming language.
The `data::cdot::*` code is based on a port of  [SACGF/cdot](https://github.com/SACGF/cdot) to Rust.

## Running Tests

The tests need an instance of UTA to run.
Either you setup a local copy (with minimal dataset in `tests/data/data/*.pgd.gz`) or use the public one.
You will have to set the environment variables `TEST_UTA_DATABASE_URL` and `TEST_UTA_DATABASE_SCHEMA` appropriately.
To use the public database:

```
export TEST_UTA_DATABASE_URL=postgres://anonymous:anonymous@uta.biocommons.org:/uta
export TEST_UTA_DATABASE_SCHEMA=uta_20210129
```

Note that [seqrepo-rs](https://github.com/bihealth/seqrepo-rs) is used for access to the genome contig sequence.
It is inconvenient to provide sub sets of sequences in SeqRepo format.
Instead, we use a build-cache/read-cache approach that is also used by `biocommons/hgvs`.

To build the cache, you will first need a download of the seqrepo [as described in biocommons/biocommons.seqrepo Quickstart](https://github.com/biocommons/biocommons.seqrepo#quick-start).
Then, you configure the running of tests for `hgvs-rs` as follows:

```
export TEST_SEQREPO_CACHE_MODE=write
export TEST_SEQREPO_PATH=path/to/seqrepo/instance
export TEST_SEQREPO_CACHE_PATH=tests/data/seqrepo_cache.fasta
```

When running the tests with `cargo test`, the cache file will be (re-)written.
Note that you have to use `cargo test --release -- --test-threads 1 --include-ignored` when writing the cache for enforcing a single test writing to the cache at any time.
If you don't want to regenerate the cache then you can use the following settings.
With these settings, the cache will only be read.

```
export TEST_SEQREPO_CACHE_MODE=read
export TEST_SEQREPO_CACHE_PATH=tests/data/seqrepo_cache.fasta
```

After either this, you can run the tests.

```
cargo test
```

## Creating Reduced UTA Databases

The script `tests/data/data/bootstrap.sh` allows to easily build a reduced set of the UTA database given a list of genes.
The process is as follows:

1. You edit `bootstrap.sh` to include the HGNC gene symbols of the transcripts that you want to use.
2. You run the bootstrap script.
   This will download the given UTA dump and reduce it to the information related to these transcripts.

```
$ bootstrap.sh http://dl.biocommons.org/uta uta_20210129
```

The `*.pgd.gz` file is added to the Git repository via `git-lfs` and in CI, this minimal database will be used.

## Some Timing Results

(I don't want to call it "benchmarks" yet.)

### Deserialization of large cdot JSON files.

Host:

- CPU: Intel(R) Xeon(R) E-2174G CPU @ 3.80GHz
- Disk: NVME (WDC CL SN720 SDAQNTW-1T00-2000)

Single Running Time Results (no repetitions/warm start etc.)

- ENSEMBL: 37s
- RefSeq: 67s

This includes loading and deserialization of the records only.
