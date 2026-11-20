# API and operations

## Main data-structure interface

From the example and implementation, the principal API surface of the SNARF structure includes:

- `snarf_init(keys, bits_per_key, batch_size, num_hash_bits)`
- `range_query(left_end, right_end)`
- `insert_key(key)`
- `delete_key(key)`
- `return_size()`

## Semantics

### `snarf_init(...)`

Builds the filter from a key set and configuration parameters.

Key inputs inferred from the source:
- a vector of numeric keys
- `bits_per_key`
- `batch_size`
- `num_hash_bits`

### `range_query(left_end, right_end)`

Checks whether a key may exist in the closed or implementation-defined numeric interval bounded by `left_end` and `right_end`.

Operationally, callers should treat this as a probabilistic query:
- `false` means no matching key was found by the filter
- `true` may be a true positive or a false positive

### `insert_key(key)` / `delete_key(key)`

Support updates after initialization, which is a key property called out by the repository.

### `return_size()`

Returns the current memory footprint metric exposed by the implementation.

## Data assumptions

The repository examples use:
- `uint64_t` keys
- synthetic numeric distributions
- sorted and unsorted vectors depending on workflow stage

The README explicitly warns about using values near integer representation limits (greater than roughly `2^60`), where overflow may lead to incorrect behavior.

## Operational model

This repository does not define service operations like deployment, health checks, metrics endpoints, or authentication.

Instead, “operations” here means maintainers running local binaries or experiments.

## External dependencies observed in source

- standard C++ library
- Boost `dynamic_bitset`

No package manager manifest beyond the handwritten build instructions is present.
