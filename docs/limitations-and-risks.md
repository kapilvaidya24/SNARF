# Limitations and risks

## Not a production service repository

There is no evidence in the checked-in files of:
- CI configuration
- packaging or release automation
- deployment manifests
- container build files
- service endpoints
- observability instrumentation
- formal test framework setup

This should be documented in RHDH as a code repository/component rather than a deployed service unless other repositories provide the operational layer.

## Integer-overflow caveat

The README warns maintainers to be careful with values near the integer representation limit (noted there as greater than `2^60`). It says overflow may cause inaccurate results.

That caveat should remain visible in any published TechDocs.

## Build portability issues

The checked-in `makefile` contains a developer-specific Boost include path. That means builds are not currently environment-neutral without adjustment.

## Include style

The repository includes `.cpp` implementation files directly from other `.cpp` files (for example in `include/snarf_hash.cpp`). This works as a repository style choice here, but it is less conventional than separating declarations into headers and linking translation units normally.

Maintainers should keep that structure in mind before refactoring or integrating it into a larger build system.

## Verification limits for this documentation pass

This documentation was generated from direct source inspection.

A compile verification step was attempted, but the current workspace lacked `g++`, so binary execution could not be confirmed here.
