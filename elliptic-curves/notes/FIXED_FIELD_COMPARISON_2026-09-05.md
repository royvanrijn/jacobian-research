# Fixed-field comparison: five deformations leave one candidate each

Status: **the frozen comparison is complete; its success criterion was not
met**. No tested deformation improves the two-dimensional restricted radical
at `u=-1`, and no nonzero inherited class was realized on a new curve.
Mathematical authority is `EC-FIXED-FIELD-COMPARISON` in
[`MATH_STATUS.json`](../../MATH_STATUS.json).

The family, labelled identification of the cubic algebra, twenty independent
anchor Kummer classes, and local-completeness argument are in the
[family proof](FIXED_CUBIC_FIELD_VARYING_CURVE_EXPERIMENT_2026-09-04.md).
The [u=-1 pairing proof](FIXED_CUBIC_U_MINUS1_CASSELS_TATE_2026-09-05.md)
is the retained baseline. Its separately known rational point is outside the
inherited span and contributes nothing to the fourth column here.

## Results

| Parameter u | dim W_u | rank CT restricted to W_u | Restricted radical dimension | Certified realized dimension |
|---:|---:|---:|---:|---:|
| -3 | 17 | 16 | 1 | 0 |
| -2 | 13 | 12 | 1 | 0 |
| -1, retained baseline | 18 | 16 | 2 | 0 |
| 0, correctness control | 20 | 0 | 20 | 20 |
| 1 | 13 | 12 | 1 | 0 |
| 2 | 13 | 12 | 1 | 0 |
| 3 | 15 | 14 | 1 | 0 |

“Certified realized dimension” means the dimension of the explicitly certified
realized subspace. Zero is a lower bound, not a proof that the actual rational
Kummer intersection is zero. At `u=0`, all twenty independent anchor points
are already certified, and all 190 computed matrix entries vanish.

The five new matrices all have the maximum possible alternating rank in
odd dimension. Each radical has exactly one nonzero element. The result is
uniform across local dimensions 13, 15 and 17: in this panel, larger local
spaces account for more obstructed directions, without more surviving
point candidates. This is a statement about these five curves only.

| u | Nonzero radical class, as an anchor mask | Bounded search result |
|---:|---:|---|
| -3 | 370313 | No point returned |
| -2 | 702532 | No point returned |
| 1 | 154245 | No point returned |
| 2 | 513585 | No point returned |
| 3 | 602844 | No point returned |

Every mask refers to the original ordered twenty-element anchor basis. Each
search used a reduced quartic with its rational square content removed,
PARI `hyperellratpoints` height argument 100,000, a twenty-second cover limit,
and an explicit check at parameter infinity. All five searches completed
without a point. Their point-or-Sha status remains **UNKNOWN**.

## Frozen scope and exact evidence

The [initial protocol](../data/fixed_field_comparison_v1.json) fixes
`u=0,-2,1,2`, sixty seconds per pairing entry, 1,200 seconds per curve for
CT, 2 GiB worker RSS and a 256 MB PARI stack. Only after all four matrices and
their point phases replayed was the
[extension selection](../data/fixed_field_comparison_extension_freeze_v1.json)
frozen at `u=-3,3`, the next two integers. It pins the first-stage evidence
hashes and fixes the final parameter set. No later score or point outcome
was used to enlarge it. The [extension protocol](../data/fixed_field_comparison_extension_v1.json)
uses the same CT and point limits.

The local extension had a 600-second limit and retained a complete
[seven-parameter local certificate](../../artifacts/generated-results/elliptic-curves/fixed_cubic_field_local_kummer_u3_v1.json).
Its replay certifies every finite local-image dimension, the real condition,
complete bad-prime support and the simultaneous kernel on all twenty anchor
classes. The old five-parameter certificate is unchanged. The new dimensions
are 17 at `u=-3` and 15 at `u=3`.

The [comparison certificate](../../artifacts/generated-results/elliptic-curves/fixed_field_comparison_v1.json)
indexes six portable compressed witness files. They retain **665 independent
matrix entries, 772 exact covering maps, 17 extra symmetry/bilinearity
checks, and five point-search observations**. The baseline computation is
not included in those counts.

The shared runtime constructs the exact labelled minimal-model transport,
parametrizes the norm conic, and optionally minimizes and reduces the
quartic. Each retained map satisfies the original two-quadric equations.
A separate square identity binds the quartic's cubic invariant to its
specified anchor class. All quartics on one curve have the same exact
invariants. Fisher's Theorem 3.1 is then replayed from cubic square roots,
primitive pairing quadratics, rational local square witnesses, and complete
prime factorizations. The elementary Hilbert symbols were also compared
with PARI during discovery. Reverse-order and sum-class checks agree with
the computed matrices; each nonzero radical generator pairs trivially in
its extra check.

No incomplete matrix is assigned a radical. No point search is allowed
outside a replayed restricted radical. The frozen cap is fifteen nonzero
combinations in increasing radical-coordinate weight, with a 400-second
point-phase limit per curve; each new curve needed only one combination.
The twenty-dimensional control uses its certified points without new point
enumeration. A separate small positive control tests recovery of actual
rational points through the new reduced-cover transport.

One initial `u=0` attempt failed during quartic reduction because completing
the square had introduced rational coefficients. The implementation now
clears their denominators by a rational square before every PARI reduction,
and checks the complete resulting map. Failed-attempt logs are retained
under the local checkpoint directory. No failed entry was treated as zero.

## Meaning and limits

For each curve let R_u be the radical of the pairing on W_u. Rational Kummer
classes in W_u lie in R_u. Thus every class in W_u outside R_u has nonzero
image in Sha[2], and the rational Kummer intersection has dimension at most
one on each new deformation. This gives no bound on rational classes outside
W_u or on the rank of the whole curve.

The remaining class may pair nontrivially with Selmer classes outside W_u,
or may remain a nontrivial Sha class despite all pairings tested here. A
bounded point-search miss decides neither possibility. The experiment
computes no full Selmer group or class group.

This panel supplies no successful rank-transfer deformation. It also gives
no evidence for spending a larger point budget on these five curves: each
has fewer remaining inherited candidates than `u=-1`. Any subsequent panel
should declare a new mathematical selection rule and bounds before running;
this frozen experiment does not continue by increasing |u|.

## Replay and regeneration

Cheap replay checks exact maps, pairings, profiles, the frozen extension gate
and positive witnesses without conic solving or point enumeration. Local
completeness is the separate source theorem and is replayed explicitly:

```sh
sage -python elliptic-curves/cas/run_fixed_cubic_field_curve_family.sage --check
sage -python elliptic-curves/cas/run_fixed_cubic_field_curve_family.sage \
  --parameter-bound 3 --check \
  --output artifacts/generated-results/elliptic-curves/fixed_cubic_field_local_kummer_u3_v1.json
sage -python elliptic-curves/cas/run_fixed_field_comparison.py --check
sage -python -m unittest elliptic-curves/tests/test_fixed_field_comparison.py \
  elliptic-curves/tests/test_sage_subspace_runtime.py
```

The global status audit still reports ten checker-source hash mismatches
already present in `HEAD` from the shared-runtime migration. Historical
bindings are preserved. The new entry and all other registry/consumer
constraints pass after isolating that existing provenance drift; `STATUS.md`
was generated from the unchanged authority using `scripts.render_status.render`.
This does not turn the full status audit into a passing check.

The eighteen targeted tests pass. They include the zero control, an empty
arithmetic cache, tampered cover maps, missing entries, wrong-parameter
binding, an out-of-radical search, a fabricated point, and a positive
reduced-cover point recovery. Negative search completeness is not asserted
by cheap replay; the original bounded CAS outputs are retained observations.

For fresh arithmetic discovery, use a fresh work directory and run both
phases for each initial parameter, in this order:

```sh
python3 elliptic-curves/cas/run_fixed_field_comparison.py --u 0 \
  --workdir artifacts/local/fixed-field-comparison-fresh
python3 elliptic-curves/cas/run_fixed_field_comparison.py --u 0 --phase points \
  --workdir artifacts/local/fixed-field-comparison-fresh
```

Repeat those two commands for `--u -2`, `--u 1`, and `--u 2`. For `--u -3`
and `--u 3`, also pass
`--protocol elliptic-curves/data/fixed_field_comparison_extension_v1.json`;
the runner first replays the pinned first-stage gate. Every cover and pair
is checkpointed through the shared fact store; interrupted work can resume
in the same directory. A completed result is preserved, so use a new
directory for a new complete run. Raw logs and supervision records remain
under `artifacts/local/fixed-field-comparison-v1/`.
