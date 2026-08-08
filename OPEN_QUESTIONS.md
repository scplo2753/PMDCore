# Open Questions

## Semantics of `--range 0`

Status: Open

The original PMDtools implementation does not clearly establish
whether a range value of zero is valid.

The current C++ implementation allows `FLAGS_range == 0`. In this
case:

- no forward Platypus statistics are collected;
- no backward Platypus statistics are collected;
- denominator tables have zero length;
- PMD likelihood calculation still runs normally.

It remains undecided whether `--range 0` should:

1. be accepted as a way to disable terminal pattern statistics;
2. be rejected as invalid input; or
3. have another compatibility behavior matching the original tool.

Do not treat the current behavior as a stable interface until this
question is resolved.

Related code:

- `include/args.list`
- `src/calPMD.cpp`
- `include/statics_types.hpp`
- `tests/calPMD_test.cpp`