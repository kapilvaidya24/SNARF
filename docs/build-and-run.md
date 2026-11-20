# Build and run

## Toolchain expectations

From the repository contents, the expected build environment is:

- `g++`
- C++17 support
- Boost headers, specifically for `boost/dynamic_bitset.hpp`

## Minimal compile example

The README shows compilation of `example.cpp` using `g++` and C++17:

```bash
g++ -std=c++17 -O3 -w -fpermissive example.cpp -o example.out
./example.out
```

The checked-in `makefile` uses this variant:

```bash
g++ -std=c++17 -O3 -w -fpermissive -I /Users/lucas/C++_lib/1.79.0/include example.cpp -o example.out
```

## Portability note

The include path in `makefile` is hard-coded to a specific developer machine:

- `/Users/lucas/C++_lib/1.79.0/include`

In another environment, maintainers should replace that with the local Boost include location or remove it if Boost is already in the compiler's include search path.

## Example behavior

The example code demonstrates:

1. building a vector of numeric keys
2. initializing `snarf_updatable_gcs_hash<uint64_t>`
3. querying a numeric range
4. inserting a key and querying again
5. deleting a key and querying again
6. showing a false-positive example

## Workload test behavior

`workload_tests.cpp` creates a much larger synthetic dataset and exercises the range filter over it. This is useful as a demonstration or performance experiment, but it implies significantly higher memory and runtime requirements than the small example.

## Verification limits in this workspace

In this documentation run, the repository was inspected directly, but an actual compile test could not be completed because `g++` was not available in the current execution environment.

That means the documented build commands are source-grounded from the repo, but not compiler-verified here.
