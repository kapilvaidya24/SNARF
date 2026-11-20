# Architecture

## High-level design

The core implementation lives in `include/snarf_hash.cpp` as the templated type:

- `snarf_updatable_gcs_hash<T>`

This structure combines several internal pieces:

1. **Learned model (`snarf_model<T>`)**
   - stored in the `rmi` field
   - used to infer a CDF-like position for a key
   - maps keys toward bit locations used by the filter

2. **Compressed bit-array blocks (`vector<snarf_bitset> bb_bitset_vec`)**
   - store Golomb-coded data in block form
   - each block is backed by `boost::dynamic_bitset<>`
   - construction is handled through helper methods such as `create_new_gcs_block` and `build_bb`

3. **Bloom filter (`BloomFilter bf`)**
   - present as an auxiliary probabilistic structure
   - initialized via `BloomFilter_init` and used with hashed values

4. **Hash tracking (`unordered_map<uint64_t, uint64_t> map_hash`)**
   - tracks hashed keys or related state for update behavior

## Construction flow

From `snarf_init(...)` and its helper calls, the initialization path is:

1. validate or derive target false-positive-related parameters from `bits_per_key`
2. set internal parameters such as:
   - `N`
   - `P`
   - `bit_size`
   - `block_size`
   - `total_blocks`
3. build the learned model with `rmi.snarf_model_builder(keys)`
4. derive bit locations for all keys with `get_locations(...)`
5. build Golomb-coded bit blocks with `build_bb(...)`
6. initialize additional hashed/filter state

## Query and update model

From the example program, the main data-structure interface includes:

- `snarf_init(...)`
- `range_query(left_end, right_end)`
- `insert_key(key)`
- `delete_key(key)`
- `return_size()`

This means the repository implements an **updatable probabilistic range filter** rather than a point-lookup hash table or a networked service.

## Supporting components

### `snarf_bitset`

Defined in `include/snarf_bitset.cpp`.

Responsibilities:
- allocate and reset bit storage
- write binary and unary-coded fragments into the bitset
- read bits back during inference/query operations
- serialize and deserialize block contents
- report approximate byte size

It depends on Boost's `dynamic_bitset`.

### `BloomFilter`

Defined in `include/bloom_filter.cpp`.

Responsibilities:
- allocate bit storage
- hash values into multiple bit positions
- support `add(...)` and `possiblyContains(...)`
- report size in bytes

### `snarf_model<T>`

Defined in `include/snarf_model.cpp`.

The full model internals are not documented in the README, but the code clearly shows it is responsible for inference and model building used to place keys into the filter structure.

## Runtime shape

This repository does **not** define a long-running daemon, HTTP API, CLI packaging, or deployment artifact. Its runtime shape is instead:

- C++ source compiled into local binaries
- example-style execution from `example.cpp`
- experiment/workload execution from `workload_tests.cpp`

So in RHDH terms, this is best documented as a **library/research component repository**, not a production service.
