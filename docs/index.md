# SNARF_hashkey

## What this repository is

`SNARF_hashkey` is a C++ code repository for **SNARF** (Sparse Numerical Array Based Range Filter), described in the repository README as a *fast and compact updatable range filter* and linked to a SIGMOD paper.

Based on the checked-in source, this repository is primarily a research/prototype-style implementation rather than a production service. It provides:

- a templated SNARF implementation in `include/snarf_hash.cpp`
- supporting model, bitset, and bloom filter code under `include/`
- a runnable example in `example.cpp`
- a workload-oriented test/demo in `workload_tests.cpp`

## What SNARF does in this repo

From the example and implementation, the main supported operations are:

- **range query**: test whether any key may exist in a numeric interval
- **insert**: add a key
- **delete**: remove a key
- **size reporting**: inspect the space used by the data structure

The example demonstrates the usual filter semantics:

- true negatives are possible
- false positives are possible
- false negatives are not expected

## Who should use these docs

These docs are aimed at maintainers or catalog users who need to understand:

- what the repository contains
- how the implementation is structured
- how to build and run the example code
- what caveats apply before operational or product use

## Source-grounded notes

These docs were generated from the repository contents present in the workspace, especially:

- `README.md`
- `example.cpp`
- `workload_tests.cpp`
- `include/snarf_hash.cpp`
- `include/snarf_model.cpp`
- `include/snarf_bitset.cpp`
- `include/bloom_filter.cpp`
- `makefile`

Where the repository is incomplete or ambiguous, the docs say so explicitly instead of guessing.
