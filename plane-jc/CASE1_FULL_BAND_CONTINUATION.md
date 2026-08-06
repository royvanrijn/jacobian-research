# Exact continuation of the `(72,108)` Case-1 lower bands

> **Status.** This is an exact characteristic-zero computation over the
> pinned degree-35 coefficient field. It continues the archived necessary
> Case-1 recurrence from bracket layer `-3` through layer `-11`, reconstructing
> all three formerly omitted `P` bands and all eight formerly omitted `Q`
> bands. Every new linear block has full column rank and introduces no new
> parameter. The calculation adds 66 compatibility equations. This closes the
> explicit band-recovery prerequisite for the alternate-chart residue; it does
> not construct a Case-1 point, since the earlier thirteen equations already
> define the empty scheme, and it does not yet perform the alternate-chart
> transport or polynomial right-component sieve.

## 1. The stopped recurrence

Write the audited Laurent pair as

\[
 P=\sum_i P_i(t)z^i,\qquad Q=\sum_j Q_j(t)z^j.
\]

In the chart `t=xy^2,z=y^-1`, the bracket layer contributed by `P_i,Q_j`
is

\[
 z^{i+j-1}\left(iP_iQ_j'-jP_i'Q_j\right).       \tag{1.1}
\]

The exact archive solved layers `1,0,-1,-2,-3`. Its last reconstructed bands
were `P_-5` and `Q_-4`; this was sufficient for the existing unit-ideal
contradiction, so the replay stopped there. The full Case-1 Newton polygons
instead require

\[
 P_{-6},P_{-7},P_{-8},\qquad
 Q_{-5},Q_{-6},\ldots,Q_{-12}.                   \tag{1.2}
\]

At a new bracket layer `k`, the only new unknowns are

\[
 P_{k-2},\qquad Q_{k-1}.                          \tag{1.3}
\]

All other pairs in (1.1) are known from preceding layers. The same triangular
linear solve used by the archive therefore continues without changing the
coefficient field or the six existing parameters. Once `P_-8` is reached,
the published polygon has no further `P` band and the recurrence determines
only the shrinking `Q` tail.

## 2. Exact ledger

The continuation gives:

| bracket layer | new `P` band | new `Q` band | columns | rank | nullity | new compatibility equations | maximum parameter degree |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `-4` | `-6` | `-5` | 11 | 11 | 0 | 6 | 8 |
| `-5` | `-7` | `-6` | 9 | 9 | 0 | 7 | 9 |
| `-6` | `-8` | `-7` | 7 | 7 | 0 | 8 | 10 |
| `-7` | -- | `-8` | 5 | 5 | 0 | 9 | 11 |
| `-8` | -- | `-9` | 4 | 4 | 0 | 9 | 12 |
| `-9` | -- | `-10` | 3 | 3 | 0 | 9 | 13 |
| `-10` | -- | `-11` | 2 | 2 | 0 | 9 | 14 |
| `-11` | -- | `-12` | 1 | 1 | 0 | 9 | 15 |

Thus all eight blocks have full column rank, the parameter count remains
exactly six, and the new compatibility count is

\[
 6+7+8+9+9+9+9+9=66.                            \tag{2.1}
\]

The deterministic
[`case1_full_band_continuation.json`](../artifacts/generated-results/case1_full_band_continuation.json)
records, for every band and compatibility block, its support, coefficient-slot
count, nonzero count, parameter-term count, maximum parameter degree, and
SHA-256 hash. It also pins the two immutable inputs:

```text
case1_checkpoint.pkl
  2dcf13d924530cdc9a8728e943efdc73d003ce1c187d5cec273f6f701e0240ba
exact_core.py
  3ba2d44e52a8028044dd73a9394449f6a26c638aebed6cc9f09288e49d77ff82
```

## 3. Mathematical interpretation

This computation does not use the conductor/contact-loss theorem. The
alternate residue is a polynomial map along the entire dicritical
`P^1`, while the conductor theorem controls classes in a finite conductor
quotient. A finite supported jet ledger cannot by itself decide a global
polynomial right component. Direct continuation of the Newton recurrence is
therefore the correct route here.

The result must also not be read as a new family of Laurent solutions. The
six parameters still satisfy the original thirteen compatibility equations,
whose ideal is already the unit ideal. Equations (1.2) are universal formal
necessary expressions: if a Case-1 solution existed, its lower bands would
be these expressions and would satisfy the additional 66 equations. Lower
layers cannot resurrect the empty upper scheme.

What has changed is the data boundary. The former statement

```text
the exact Case-1 residue cannot be constructed because eleven bands are absent
```

is no longer current. The complete band expressions are now algorithmically
derived and hash-pinned. The remaining tasks are:

1. serialize a reusable full-band checkpoint;
2. transport the reconstructed pair through the certified alternate chart
   `Y2=Y1-beta/X+delta` and restrict it to `X=0`;
3. extract the exact degree `(8,12)` residue coefficient vectors on both sign
   branches; and
4. run the existing general degree-two and degree-four right-component
   remainder sieve on those vectors.

Because the coefficient scheme is independently empty, steps 2--4 are a
conceptual compression and reuse project, not a prerequisite for the proved
degree-125 frontier.

## 4. Reproduction

The quick pinned audit checks the full deterministic ledger and exactly
replays the first omitted layer:

```bash
.venv/bin/python scripts/verify_case1_full_band_continuation.py
```

The complete specialized replay is:

```bash
.venv/bin/python scripts/continue_case1_full_bands.py \
  --stop-layer -11 \
  --ledger-output artifacts/generated-results/case1_full_band_continuation.json
```

It uses `python-flint` exact arithmetic and is intentionally not part of the
routine regression suite. On the reference run, the bottom exact
convolutions dominated the runtime; the final three assemblies took about
8.4, 13.0, and 19.8 minutes. The script never mutates the external replay.
