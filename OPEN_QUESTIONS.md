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
- `include/statistics/statistics_types.hpp`
- `tests/calPMD_test.cpp`

## Reference or read context for 5-prime CpG deamination

Status: Resolved — use the reconstructed reference

In the original `pmdtools.0.60.py`, the 5-prime branch of
`options.deamination` determines CpG context using the observed read:

```python
if real_read[i + 1] != 'G':
    continue
```

The corresponding 3-prime branch uses the reconstructed reference:

```python
if real_ref_seq[i - 1] != 'C':
    continue
```

Other CpG-sensitive paths, including Platypus statistics and PMD
likelihood scoring, use the reconstructed reference for context on
both ends.

The C++ implementation intentionally differs from the original script here:
the 5-prime deamination branch now uses `real_ref_seq[i + 1]`. This makes
`--CpG` and `--noCpG` classify 5-prime C sites from the reconstructed
reference, consistently with the 3-prime branch, Platypus statistics, and PMD
likelihood scoring.

Consequently, when the read and reconstructed reference disagree at the base
following a 5-prime C, the reference determines whether that C is in CpG
context. This is a deliberate behavior change rather than strict compatibility
with `pmdtools.0.60.py`.

Related code:

- `pmdtools.0.60.py`, `options.deamination`
- `src/calPMD.cpp`, `calPMD::deamination`
- `src/calPMD.cpp`, `calPMD::computeDegradationScore`
- `src/calPMD.cpp`, `calPMD::platypus_forward`
- `tests/deamination_tests/unit/deamination_CpG_test.cpp`

## Simultaneous use of `--deamination` and `--platypus`

Status: Resolved for the current C++ implementation; retained here as a
compatibility note.

The original `pmdtools.0.60.py` allows `options.deamination` and
`options.platypus` to be enabled at the same time. Both branches update
some of the same mismatch dictionaries.

The Platypus branch is executed first:

```python
if options.platypus:
    # Updates mismatch_dict, mismatch_dict_rev,
    # mismatch_dict_CpG, and mismatch_dict_CpG_rev.
```

The deamination branch is then executed independently:

```python
if options.deamination:
    # Updates mismatch_dict and mismatch_dict_rev.
    continue
```

When both options are enabled, a C- or G-reference position may
therefore be counted twice in `mismatch_dict` or `mismatch_dict_rev`:

1. once by the Platypus branch;
2. once by the deamination branch.

The two branches do not always update identical dictionaries. In CpG
mode, Platypus stores CpG-context observations in:

```text
mismatch_dict_CpG
mismatch_dict_CpG_rev
```

The original deamination branch filters positions by CpG context but
still stores the resulting observations in:

```text
mismatch_dict
mismatch_dict_rev
```

Consequently, the original program's Platypus output may include
observations added by the deamination branch, and some observations may
be counted twice when both options are active. The source does not make
clear whether this coupling was intentional or an accidental consequence
of reusing global dictionaries.

### Current C++ behavior

PMDCore allows both options to be enabled, but deliberately keeps their
statistics independent:

- Platypus updates its mismatch dictionaries and nucleotide denominator
  tables;
- deamination updates a separate `deamination_statics_t`, containing
  forward and reverse terminal counts;
- thread-local results for the two modes are merged independently;
- deamination and Platypus results are printed from their respective
  statistics, so neither mode contributes counts to the other mode's
  output.

Each mode can therefore observe the same input position, but this is not
double counting within a shared result. It is one observation in each of
two independent outputs.

This is an intentional compatibility difference. When both options are
enabled, PMDCore does not reproduce the original program's shared-dictionary
side effects or its possible double counting. PMD likelihood scoring remains
skipped after the deamination branch, matching the original control flow.

`inputParams_validator()` emits a `Caution` message when both options are
present so that users comparing results with the original implementation are
made aware of this difference. The combination remains valid and does not
terminate the program.

Dedicated compatibility tests should cover:

- non-CpG C-reference positions;
- non-CpG G-reference positions;
- CpG C-reference positions;
- CpG G-reference positions;
- positions inside and outside `FLAGS_range`;
- independent Platypus and deamination counts when both options are enabled;
- whether PMD likelihood scoring is skipped after deamination, as in
  the original implementation.

Related code:

- `pmdtools.0.60.py`, `options.platypus`
- `pmdtools.0.60.py`, `options.deamination`
- `src/calPMD.cpp`, `calPMD::calPMD_loop`
- `src/calPMD.cpp`, `calPMD::platypus_forward`
- `src/calPMD.cpp`, `calPMD::platypus_backward`
- `src/calPMD.cpp`, `calPMD::deamination`
- `include/statistics/statistics_types.hpp`, `platypus_statics_dicts_t`
- `include/statistics/deam_types.hpp`, `deamination_statics_t`
- `main_multithreaded.cpp`, independent thread-result merging and output
- `src/arguments.cpp`, simultaneous-option caution

## Reference-dependent checks when an existing DS field skips reconstruction

Status: Open for future DS reuse; current PMDCore behavior is defined below.

In the original `pmdtools.0.60.py`, reconstruction of `real_ref_seq` is
conditional. A record with an existing `DS:Z:` field can skip reconstruction
when no option requiring the alignment or a new PMD calculation is enabled:

```python
if (DSfield == False) or options.writesamfield or (options.basic > 0) or \
        options.terminal or (options.perc_identity > 0.01) or \
        options.printalignments or options.adjustss or options.adjustbaseq or \
        options.adjustbaseq_all or options.deamination or options.dry or \
        options.estimate or options.first:
    # Reconstruct real_read and real_ref_seq.
```

However, the GC-content and reconstructed-reference validity checks are
executed unconditionally after this block:

```python
GCcontent = 1.0 * (real_ref_seq.count('G') +
                   real_ref_seq.count('C')) / readlen
if GCcontent > options.maxGC: continue
elif GCcontent < options.minGC: continue

if ('G' not in real_ref_seq and 'C' not in real_ref_seq and
        'T' not in real_ref_seq and 'A' not in real_ref_seq):
    continue
```

Because `real_ref_seq` is assigned only inside the reconstruction block, the
original program may encounter an undefined variable when the first eligible
record skips reconstruction. If an earlier record did reconstruct a reference,
the later record may instead reuse that earlier record's stale `real_ref_seq`
because Python loop bodies do not introduce a new variable scope.

### Current C++ behavior

PMDCore currently disables reuse of an existing `DS:Z:` value. The input tag
does not suppress reference reconstruction or PMD calculation. Every record
that reaches this part of the pipeline must therefore have its reference
reconstructed from the MD tag, after which the normal GC-content and
reference-validity checks are applied.

This temporary policy means that:

- records with and without an existing `DS:Z:` tag follow the same
  reconstruction and calculation path;
- an existing DS value is not used for threshold filtering;
- Platypus, deamination, masking, and PMD likelihood calculation are not
  skipped merely because an input DS tag exists;
- records still require a usable MD tag even when they already contain a DS
  value; and
- PMDCore avoids both reuse of a stale reference and checks against an empty
  reference sequence.

This behavior intentionally prioritizes consistent processing of currently
implemented features over the original optimization of reusing a stored PMD
score.

If DS reuse is implemented in the future, it remains undecided whether such
records should:

1. skip GC-content and reference-validity checks when reconstruction is not
   otherwise required;
2. reconstruct the reference solely to apply those checks; or
3. follow another explicitly documented compatibility policy.

Until that decision is made, the unconditional-reconstruction policy above is
the defined PMDCore behavior. The original undefined/stale-variable behavior
should not be treated as a stable compatibility requirement.

Related code:

- `pmdtools.0.60.py`, conditional reference reconstruction and subsequent
  GC/reference checks
- `main_multithreaded.cpp`, disabled DS reuse, unconditional reconstruction,
  and reference-dependent checks
- `src/seqProcedures.cpp`, `isGCcontentInRange`
- `src/Filters.cpp`, `badRefSeq_Vailder`
