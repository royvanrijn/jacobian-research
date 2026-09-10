# Autonomous search over constructed fibrations

The active programme constructs elliptic fibrations on both X948 and X1092,
then evaluates modest, fixed specialization panels. It replaces the stopped
[parameter foundry](HIGH_RANK_SEARCH_FOUNDRY_2026-09-09.md). Operational counts
are in the [live report](../../artifacts/generated-results/elliptic-curves/parent-foundry-v2/REPORT.md);
this note records the mathematical gates and frozen protocol.

Two cold-reconstructed commissioning examples establish that the new lane is
executable. They are genuine changes of elliptic fibration on the declared
existing K3 surfaces, with explicit equations and sixteen rational sections.

| Constructed fibration | Source NS determinant | Generic rank | Height determinant | Compact coefficient bits | Commissioning fibre |
|---|---:|---:|---:|---:|---|
| `x948-11952-a1-00414` | 948 | exactly 16 | 474 | 158 | t=1/2: bound 16 after 8 calls |
| `x1092-det1092-a1-121385` | 1092 | exactly 16 | 546 | 129 | t=-1/5: 16→17 at call 1, →18 at call 4; bound 18 after 8 calls |

The complete [equations, sections and point certificates](../../artifacts/generated-results/elliptic-curves/parent-foundry-v2/commissioning/)
are retained. These are new fibration classes relative to the declared
repository roster, not claims of unpublished K3 surfaces, new world records,
or superior specialization tails. The two commissioning tests are not included
in the production panel statistics. A production slot can repeat one of their
parameters, but receives the unchanged frozen generic-only detector; no
commissioning exceptional point is supplied to it.

## Why A1 is still an admissible search lane

The [104-fibre sensitivity experiment](MW16_SENSITIVITY_RECOVERY_2026-09-05.md)
completed 856 boxes with no prospective gains. That was a real detector null,
despite recovering all 55 control directions. Later compact searches on the
same five A1 families produced, for example, three new bound-26 curves in the
[extended twenty-fibre campaign](MILLION_HEIGHT_AND_MW16_EXTENSION_2026-09-06.md).
The [broader higher-parameter trial](BROAD_MW16_HIGHER_POPULATION_2026-09-06.md)
produced ten seeded fibres among sixty and two bound-24 curves; the subsequent
[corrected-score trial](CORRECTED_MW16_HIGHER_POPULATION_2026-09-06.md) seeded
five of sixty, with one bound-23 curve. These exposures differ and cannot be
pooled as one homogeneous discovery rate.

The pool is therefore usable but not established as superior. The new search
does not interpret those failures as rank upper bounds or repeat the same
five parents indefinitely. It tests newly constructed fibrations on both
surfaces, alongside the existing six R17, five MW16 and det1092 MW17 baselines.

The recent [norm-eight degeneration obstruction](DET1092_NORM8_DEGENERATION_AND_SMOOTH_ATLAS_OBSTRUCTION_2026-09-09.md)
concerns rational degenerations of one particular `5I2+14I1` pencil and the
use of its rational curves to supply seeds on an original fibre. That
fibration has generic MW12. The present construction explicitly rejects
multi-I2 pencils and directly searches fibres of new `I2+22I1` fibrations.
It does not reopen that obstructed rational-degeneration mechanism.

## The difference from the blinded curve302 MW16 cores

| Object | Reducible-fibre root rank | Full generic rank | Meaning of the supplied sixteen |
|---|---:|---:|---|
| Genuine A1 parent on X948 | 1 | 16 | Full-rank generic subgroup; no seventeenth independent section exists |
| Genuine A1 parent on X1092 | 1 | 16 | Same generic upper bound, on a different NS lattice |
| Recovered curve302 parent on X1092 | 0 | 17 | In the blinded experiment, a deliberately incomplete subgroup of known MW17 |

Both surfaces have geometric and rational Picard rank 19. Shioda–Tate gives
`19−2−1=16` for an A1 pencil and `19−2=17` for a rootless pencil. The blinded
experiment did not change the fibration or create a root component. Its
seventeenth direction existed before a section was hidden.

There is also an [exact lattice distinction](../../artifacts/generated-results/elliptic-curves/parent-foundry-v2/commissioning/lattice-comparison.json). The blinded sixteen-section
principal Grams are integral even forms; their determinants range from 1445
to 8237. The actual new A1 bases have half-integral height pairings, including
height 7/2, because the I2 component contributes the Shioda correction. Their
determinants are respectively 474 and 546. The full rootless curve302 parent
has saturated height determinant 1092 and 24 I1 fibres.

These differences explain why generic rank completion is possible for the
blinded cores and impossible for an exact A1/MW16 pencil. They do **not**
explain the fourteen exceptional directions of curve302. Its own rootless
parent had a [48-fibre prospective null](DET1092_LOW_SHELL_CASCADE_2026-09-07.md).
The programme measures specialization tails instead of assuming that the
determinant, generic rank or known extreme fibre predicts them.

No supplied input in this launch is a genuinely incomplete arithmetic MW16
fibration with upper bound at least 17. Consequently the RR section-completion
stage is theorem-blocked for these A1 inputs. New rank-17 parents would require
a different fibration construction or a new admissible incomplete input;
artificially deleting a section is not used as a substitute.

## Construction and proof gates

The X948 lane starts with the exact 1,266-class A1 stratum in the committed
11952 norm-eight table. The X1092 lane starts with 63,922 norm-eight proposals
from its committed degree-two census. A deterministic 4,096-proposal lattice
preflight finds 133 singleton unsigned minimum classes; these receive the
equation gate first. The remaining unexamined proposals stay in the queue.
The count 133/4096 is a lattice observation, not a specialization success rate.

For a norm-eight section P_w, construct `D=O+P_w`, the residual-chord quartic,
its Jacobian, and degree-one source sections. Require coefficient degrees
8,12, finite discriminant degree 22 and squarefree finite discriminant. This
establishes `I2 at infinity + 22 I1`. Require sixteen explicit sections with
positive exact Shioda Gram, including the A1 correction. Compactification
uses bounded auxiliary minimization and exact weighted coefficient identities.
Every section is transported and checked as a rational-function identity.

An exact separating invariant prevents relabelled fibrations entering as new
parents. Any equivalence must preserve the unique I2 fibre, so in the raw
coordinates it is affine. Depress the monic finite discriminant and compare
ratios of equal-weight coefficients. Different invariant vectors prove
inequivalence even over the algebraic closure. Equal fingerprints are treated
conservatively as duplicates or unresolved collisions. The nine known A1
presentations reproduce exactly five keys, including the required curve398
collision regression. A new I2 pencil also differs from every 24-I1 source.

The commissioning replay reconstructs both complete parent packets from their
generic inputs in fresh temporary directories. It then replays independent
finite certificates for their labelled specialization points. Sixteen exact
generic sections with independent images at a smooth rational fibre are
independent over Q(t); the Picard/root upper bound supplies equality 16.
The geometry reconstruction shares its construction implementation, while
the final point independence is checked by two finite-group implementations.
No fully independent symbolic implementation or external review is claimed.

## Autonomous allocation and evidence

Each new parent receives twelve frozen parameter slots: eight at projective
height 2..8 and four at height 17..64. These are coordinate-dependent sampling
bands, not intrinsic fibration invariants. No parameter score discards a slot.
Each slot receives an allowance of 64 point calls, height 125000, ten seconds
per call and five seconds per map. Singular fibres and failed finite generic
certificates remain unresolved slots, with no replacement and no rank claim.
Point or map timeouts are recorded as incomplete exposure.

The generic parent banks work in rank 16 or 17. Each bank samples 64 fresh
full-space parity classes, checks their exact minima using two CVP solvers,
retains all minimum representatives, and selects twelve deep and four
exploratory classes. After each gain the existing adaptive engine rebuilds
the subgroup. Every returned cloud is reconciled before a result is accepted.

The first twelve slots give separate observed detection fractions for
certified jumps at least 3,5,8. Missing or incomplete exposures are shown,
including a conservative upper fraction if all unresolved slots succeeded.
These are finite protocol measurements, not estimates of exact-rank density
with established confidence guarantees. Promoted or deeper fibre results
are kept separate from this first-panel table.

A completed panel promotes its parent if some fibre gains at least five
directions or at least three fibres seed. Further panels grow in bounded
increments. Approximately 30% of evaluation worker time goes to productive
individual cascades, weighted by gain recency, gain rate and proximity to 32.
Stalled fibres cool after 128 no-gain calls below 28, or 256 at 28 and above.
Individual exploitation allowances grow 64,96,128,... up to 256 calls per job.
The old six-parent parameter search is not inherited as a stalled work queue.

Four workers run from a frozen source snapshot. Construction jobs have
600-second bounds; point evaluation jobs have 7200-second bounds and 3-GiB
process-tree RSS caps. A guardian renews each 64-job controller run, including
clean exits, and restarts failed controllers with backoff. One interrupted
point job is retried from its saved calls; repeated failures are quarantined.
There is no daily budget and no AI/model call in the controller or workers.
An explicit stop drains the active jobs. A 20-GiB disk reserve causes automatic
waiting and resumption when space returns. On eventual proposal-pool exhaustion,
the controller opens larger fresh panels on accepted parents.

Exact equations, independent points, trajectories, maps, full clouds and
replay receipts are preserved. Novelty comparisons occur after evaluation,
against the pinned public/repository catalogues and earlier foundry fibres.
Possible conductor payoffs receive a cheap discriminant-based flag; this
launch does not add an unbounded factorization lane. The separate existing
conductor campaign remains untouched. Curated inventory additions remain
separate from machine result publication.

## Operation and replay

From the repository root:

```sh
python3 research/elliptic-curves/cas/run_parent_foundry.py status
python3 research/elliptic-curves/cas/run_parent_foundry.py stop
python3 research/elliptic-curves/cas/run_parent_foundry.py launch

sage -python research/elliptic-curves/cas/replay_parent_foundry_constructions.sage
python3 research/elliptic-curves/cas/verify_parent_foundry_certificate.py \
  research/artifacts/generated-results/elliptic-curves/parent-foundry-v2/commissioning/point-1.json \
  --parent research/artifacts/generated-results/elliptic-curves/parent-foundry-v2/commissioning/parent-1.json
```

The [configuration](../../artifacts/generated-results/elliptic-curves/parent-foundry-v2/config.json),
[source manifest](../../artifacts/generated-results/elliptic-curves/parent-foundry-v2/manifest.json)
and frozen runtime archive bind this launch. The initial prepared parent-foundry
v1 snapshot was never launched; it is retained as development evidence.
The old high-rank-foundry v3 drained before the replacement launch.

The live v2 snapshot contains an opt-in generic-seed gate in its frozen V3
preparer. In the working sources the equivalent gate is extracted to
`parent_foundry_seed.py`, preserving the byte-pinned historical V3 checker.
Production always runs the manifest-bound snapshot, not changing checkout files.
