# Repository layout

## Top level

### `README.md`

Very brief project overview with:
- the SNARF name expansion
- a short claim about being a fast and compact updatable range filter
- a paper link
- a simple compile/run example
- a warning about integer-overflow behavior near large numeric limits

### `makefile`

Provides a minimal build target for `example.cpp` and a `clean` target.

Notable details:
- uses `g++`
- compiles with `-std=c++17 -O3 -w -fpermissive`
- includes a machine-specific Boost include path:
  - `/Users/lucas/C++_lib/1.79.0/include`

That hard-coded path is not portable and should be adjusted for other environments.

### `example.cpp`

Primary example and interactive/demo-oriented driver.

Contains:
- synthetic key-generation helpers for normal, uniform, and exponential distributions
- false-positive checking helper logic
- SNARF initialization and experiment code
- an interface demonstrating range-query behavior and update operations

### `workload_tests.cpp`

A heavier workload example using a large generated dataset (`N = 10,000,000`) and demonstrating range query, insert, and delete operations.

## `include/`

### `include/snarf_hash.cpp`

Main implementation file for the updatable hashed/Golomb-coded SNARF structure.

### `include/snarf_model.cpp`

Model-building and inference logic used by the filter implementation.

### `include/snarf_bitset.cpp`

Bitset abstraction built on Boost `dynamic_bitset`, with serialization helpers.

### `include/bloom_filter.cpp`

Auxiliary bloom filter implementation.

## Generated or local-only artifacts

The repository also contains local/editor artifacts that are not core documentation targets:

- `.vscode/`
- `.DS_Store`
- compiled outputs such as `example.out` and `workload_tests.out`

These should generally be treated as local development artifacts rather than source of truth.
