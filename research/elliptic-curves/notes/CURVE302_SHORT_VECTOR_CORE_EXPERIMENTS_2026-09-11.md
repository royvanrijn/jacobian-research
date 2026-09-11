# Curve302 short-vector/core experiments

These experiments test the strongest interpretation left after the closure work:
the fourteen seeded V3 trajectories may be repeatedly recovering a distinguished
low-height **integral sublattice** of the displayed quotient
`D/M17 ~= Z^14`, rather than following a predictive adaptive gradient.

They consume only the sealed `curve302-closure-structure-v1` outputs
(`quotient-relations.json`, `trajectories.json`, `REPORT.json`). No point search,
Mestre scoring, V3 execution, rank search or new curve construction is launched.

**Completed result:** the intrinsic norm-shell filtration recovers the common
rank-one and rank-three local cores, but **does not reproduce the entire V3
common-core chain**. All 180 acquisitions have exact ranks in the complete
1,288,441-direction vocabulary. The full deterministic check passes.
See [completed results](#completed-real-results) and the
[compact certificate bundle](../../artifacts/generated-results/elliptic-curves/curve302_short_vector_core_v1/manifest.json).
Status authority: `EC-CURVE302-SHORT-VECTOR-CORE-20260911`.

Run from `research/`:

```sh
python3 elliptic-curves/cas/run_curve302_short_vector_core.py run
python3 elliptic-curves/cas/run_curve302_short_vector_core.py status
python3 elliptic-curves/cas/run_curve302_short_vector_core.py check
```

Use `--source` if the passed closure-structure evidence is not in the canonical
local folder. The default output is
`artifacts/local/elliptic-curves/curve302-short-vector-core-v1/`.

## Exact performance fixes after the first real run

The first attached run stopped without a receipt and is retained as UNKNOWN.
The detached `curve302-short-vector-core-v2` run completed and sealed exact
enumeration: **1,288,441 primitive directions**, 5,789,129 nodes, covering all
180 acquisitions. Filtration then reached its 3,600-second cap; its timeout
receipt and the complete enumeration remain preserved. No basin result or
completed full check came from that run.

The optimized implementation retains the same mathematical definitions and
limits:

- Incremental rational echelon reduction replaces repeated general matrix
  rank calls. Filtration consumes the complete final norm shell, then stops
  after reaching rank 14: saturation is Z^14 at every later shell.
- A single shell scan computes all observed rank intervals with exact ties.
- Basin annihilators are cleared to primitive integer rows. Optional NumPy
  int64 dot products run only when the exact bound
  `sum(abs(a[j])*max(abs(v[j]))) <= 2^63-1` proves every intermediate fits.
  Unsafe cases and environments without NumPy use arbitrary-size Python
  integers. No floating-point containment or cutoff is introduced.
- Enumeration obtains each leaf's norm from the exact accumulated LDL sum,
  instead of evaluating all 196 matrix terms again. Full deterministic `check`
  still re-enumerates the entire ball and compares the TSV byte-for-byte.

Seventeen tests pass, including comparison with exact matrix ranks and rational
basin tests, overflow fallbacks, exact ties and boundary vectors, a full
synthetic pipeline/check, successful enumeration reuse and rejection of a
corrupt donor before output creation.

From the repository root, the fresh optimized run is:

```sh
sage -python research/elliptic-curves/cas/run_curve302_short_vector_core.py run \
  --reuse-enumeration research/artifacts/local/elliptic-curves/curve302-short-vector-core-v2 \
  --folder research/artifacts/local/elliptic-curves/curve302-short-vector-core-v3
sage -python research/elliptic-curves/cas/run_curve302_short_vector_core.py check \
  --folder research/artifacts/local/elliptic-curves/curve302-short-vector-core-v3
```

`--reuse-enumeration` is for a fresh folder only. It verifies the donor plan,
stage seal, output hashes, identical source bindings, exact observed norm bound
and the current node/vector limits. Imported files are copied and rehashed;
the new plan and import receipt preserve their producer's code hashes. The
failed donor is never resumed, edited or relabelled. Sage Python is used here
because system Python lacks SymPy; the installed SymPy is 1.14.0.

## Completed real results

Optimization commit `5195d945` completed the fresh v3 run in 52.17 seconds:
filtration 20.81 seconds, ranks/basins 29.43 seconds, plus verified enumeration
import. The preceding filtration attempt reached its 3,600-second cap. The
full deterministic check took 272.06 seconds, re-enumerated the entire ball,
and matched enumeration metadata, all TSV bytes, filtration and ranks/basins
outputs exactly. It did not merely trust the imported enumeration.

The [export](../../artifacts/generated-results/elliptic-curves/curve302_short_vector_core_v1/manifest.json)
preserves output hashes, the frozen plan and import provenance, execution
receipts and check log. The full TSV stays under local artifacts, with hash
`173e1b5f4ff9a250e451957ef200c4e457f99b0d8fec56da392a006fe7e33b96`.
A separate Sage integer-module comparison verifies that all 14 independently
reconstructed common intersections equal the earlier exact-core certificate.

Write L1,L2,L4 for recovered-local-01,02,04. In units Q(v)/10^6, the first
three primitive directions and saturation steps are:

| Primitive rank | Direction | Norm | Intrinsic saturated lattice |
| --- | --- | --- | --- |
| 1 | L2 | 16.4776247118 | <L2> |
| 2 | L2-L4 | 17.6213023391 | <L2,L4> |
| 3 | L1 | 19.8433942551 | <L1,L2,L4> |

The observed common rank-two lattice is <L1,L2>. Since L4 enters the intrinsic
filtration before L1, **no norm threshold produces that rank-two common core**.
The first and third local cores agree exactly. The later rank-five, rank-six
and rank-nine common cores are also absent as exact intrinsic steps:

| Trajectory dimension | Runs | Common rank | Intrinsic rank when first contained | Exact intrinsic step? |
| --- | --- | --- | --- | --- |
| 6 | 14 | 1 | 1 | Yes |
| 8 | 14 | 2 | 3 | No |
| 9 | 14 | 3 | 3 | Yes |
| 11 | 14 | 5 | 10 | No |
| 12 | 14 | 6 | 11 | No |
| 13 | 13 | 9 | 11 | No |
| 14 | 13 | 14 | 14 | Yes, full ambient lattice |

The intrinsic filtration reaches full rank 14 at norm 27.4396586377, within the
first 29 primitive directions. Its dependence only on the rounded quotient
metric does not make it the same filtration as the common V3 histories.

The complete rank census has no missing-vocabulary cases:

| Exact static rank cutoff | Acquisitions |
| --- | --- |
| Top 10 | 59/180 |
| Top 100 | 83/180 |
| Top 1000 | 136/180 |

The two middle ranks are 168 and 171, giving the usual median 169.5. The output's
`median_rank_worst` retains the upper middle rank, 171; the maximum is 1,288,441.
The five most frequent acquired primitive directions have global
static ranks 1,8,3,4,2, respectively. This supports recurrence of very short
directions without establishing the complete intrinsic-filtration explanation.

The 96 basin rows must be read with their `deficit_before` fields. At the
landmarks 6,8,9, respectively 12,13,13 of 14 prefixes already contain the target
core; their basin fraction is forced to one. Among the two prefixes still
missing L2 at dimension 6, the actual-norm counts are 1/1 and 104/12,580.
At the rank-two and rank-three landmarks the one remaining deficient prefix
has counts 1/2 and 1/1. Full-rank endpoint basin fractions are also forced to one.
Thus the large aggregate basin means do not by themselves establish a broad
propagation basin. The complete per-prefix counts are retained, with no
replacement scoring rule or prospective probability interpretation.

## Experiment 1: complete primitive-vector vocabulary

Let `Q` be the exact rational 14-dimensional Schur form stored by the previous
closure experiment (exact for its entrywise-rounded height metric). Compute

```
H = max Q(v)
```

over the 180 already acquired primitive quotient vectors, then enumerate **all**
primitive integer directions modulo sign satisfying `Q(v) <= H`.

Enumeration is exact Fincke--Pohst recursion using a rational LDL decomposition.
Integer interval endpoints are obtained with integer square roots; floating point
cannot exclude a boundary vector. `--max-directions` and `--max-nodes` are
fail-closed resource guards: hitting either makes the experiment UNKNOWN rather
than returning a truncated vocabulary.

This removes the 49/180 and 70/180 vocabulary-coverage limitation of the previous
next-move experiment. Every acquisition must occur in the complete vocabulary.

## Experiment 2: intrinsic saturation filtration

Order the complete primitive directions by exact `Q`-norm. At every norm shell
where the rational span grows, form

```
Lambda(h) = Sat_Z < v primitive : Q(v) <= h >.
```

Smith/Hermite arithmetic returns a canonical saturated integral basis. This gives
an intrinsic short-vector filtration depending only on the quotient height form.

Independently reconstruct every one of the **194** actual trajectory prefix
lattices. Each is required to have saturation index one. At every attained
quotient dimension, intersect the prefix lattices of all runs present there to
recover the exact common integral core. The experiment then asks when the
intrinsic `Lambda(h)` first contains, or exactly equals, each observed common
core.

A close match would support the statement that V3 is exposing the successive
short-vector saturation filtration of the exceptional quotient. A mismatch is a
clean falsification of that stronger explanation.

## Experiment 3: exact acquisition ranks and basin multiplicity

Every one of the 180 acquisitions receives its exact static-height rank interval
inside the complete primitive vocabulary, including exact norm ties. There are no
out-of-vocabulary events.

For each dimension at which the common integral core grows, inspect every actual
prefix one step earlier. For every complete-vocabulary short direction `v` not
already in that prefix, test whether

```
Sat(prefix + Z v)
```

contains the next common core. Since the inputs are saturated, this is exactly a
rational-span containment test. Report basin fractions both:

- through the actual next acquisition's norm; and
- through the experiment's complete global bound `H`.

This distinguishes a simple **energy mechanism** (the core is represented by a
few exceptionally short vectors) from a **multiplicity/basin mechanism** (many
different short extensions saturate into the same underlying core).

## Interpretation boundary

The quotient form is the already published rounded-height Schur form, not an exact
canonical-height matrix. The experiment is retrospective and uses the 180 known
acquisitions only to choose the complete enumeration bound and to measure the
observed common cores. It does not infer a V3 transition probability, prove that
a short lattice vector has a rational point/chart of bounded complexity, or prove
a self-propagation theorem.

The original conditional target below remains unproved; the completed
experiment does not establish it:

> A distinguished low-height saturated sublattice of `D/M17` has a large basin of
> short primitive extensions under the half-lattice geometry; once enough of that
> sublattice is acquired, later exceptional directions admit substantially cheaper
> half-lattice representatives.
