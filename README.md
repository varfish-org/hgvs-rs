[![CI](https://github.com/bihealth/hgvs-rs/actions/workflows/rust.yml/badge.svg)](https://github.com/bihealth/hgvs-rs/actions/workflows/rust.yml)
[![codecov](https://codecov.io/gh/bihealth/hgvs-rs/branch/main/graph/badge.svg?token=aZchhLWdzt)](https://codecov.io/gh/bihealth/hgvs-rs)

# hgvs-rs

This is a port of [biocommons/hgvs](https://github.com/biocommons/hgvs) to the Rust programming language.

## Running Tests

The tests need an instance of UTA to run.
Either you setup a local copy (with minimal dataset in `tests/data/data/*.pgd.gz`) or use the public one.
You will have to set the environment variables `TEST_UTA_DATABASE_URL` and `TEST_UTA_DATABASE_SCHEMA` appropriately.
To use the public database:

```
export TEST_UTA_DATABASE_URL=postgres://anonymous:anonymous@uta.biocommons.org:/uta
export TEST_UTA_DATABASE_SCHEMA=uta_20210129
$ cargo test
```
