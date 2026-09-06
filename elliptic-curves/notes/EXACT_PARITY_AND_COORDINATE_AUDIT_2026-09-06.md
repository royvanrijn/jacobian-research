# Exact parity geometry and coordinate visibility

The six displayed rank-17 Gram lattices now have an exact certificate for
**all maximum parity classes**. The 48 previously omitted chart attempts
and all 301 unfinished boxes on the new rank-at-least-26 curve have also
been completed and replayed. Neither experiment increased a certified
rank. The inventory remains 39 curves with lower bounds 22–26; new
rank-at-least-28/32 curves remain open.

## Exact discrete lattice statement

For each integral Gram matrix `G` in the compact six-family atlas, define

\[
 m_G(\bar a)=\min_{z\in\mathbf Z^{17},\ z\equiv a\pmod2}z^TGz.
\]

Then `max m_G = 12` in all six frames. In atlas order the complete numbers
of classes attaining 12 are:

| Family | Classes with exact minimum 12 | Outside the initial 43 |
|---|---:|---:|
| `103b2` | 43 | 0 |
| `11952` | 49 | 6 |
| `074d9` | 43 | 0 |
| `07ca9` | 43 | 0 |
| `08234` | 43 | 0 |
| `08f72` | 49 | 6 |

The [complete dataset](../../artifacts/generated-results/elliptic-curves/r17_exact_maximum_parity_classes_v1.json)
contains the matrices, all 270 maximum-class masks and witnesses. This is a
statement about the discrete quotient of each **displayed integral lattice**.
It is not a continuous covering-radius theorem, a Mordell–Weil saturation
statement or a specialized rank bound.

For the proof, all 786,432 census vectors are checked with exact integer
arithmetic: each has the specified parity and norm at most 12. Every class
outside the displayed maximum-candidate lists already has a witness of norm
at most 10. For each remaining class, an integral unimodular change of basis
is certified by its inverse and the identity `H = U^T G U`. Rational LDL
factorization proves positive definiteness. Descending enumeration exhausts
every vector in the transported parity coset of norm at most 10, using
rational inequalities and conservative integer coordinate intervals.
Every such list is empty. Since the integral Gram lattice is even and a
norm-12 witness exists, its minimum is exactly 12.

Sage LLL supplies only an integer basis change. The proof does not rely on
LLL optimality or a floating CVP solver. The
[Sage-free checker](../cas/export_r17_exact_maximum_parity_classes.py) replays
all six censuses and all 270 empty ellipsoids in about 16 seconds locally.
Two targeted enumeration tests pass, including exhaustive comparison on a
skew two-dimensional lattice.

The 49-versus-43 counts were already reported in the
[historical calibration](HALF_LATTICE_FAKE_DESCENT_REPLAY_2026-09-04.md).
The 43 cap was an intentional equal-budget choice. The additions here are
exact minimum certification and execution of the omitted-class experiment.

## Completed finite search gaps

The eight retained `11952` and `08f72` candidates each received their six
omitted classes: **48 charts**, at height 100,000 and four seconds per chart,
under 90-second / 1-GiB worker caps with at most two workers. All chart and
admission histories, eight complete point clouds and eight standalone rank
checks pass. Lower bounds remain respectively 19,20,17,17,20,18,17,17 in the
retained ledger order. Individual boxes may be partial: this experiment
completed the declared attempts, not every denominator interval.

Exact pointed-involution compression reduced 68 raw admission events to
34 retained representatives. All raw points and the relations `P+Q=C`
remain available. Compression preserves the generated subgroup; it does
not prove a rank upper bound or a wall-time speedup.

On the new rank-at-least-26 curve, a separate fixed-tail experiment searched
exactly the unvisited denominator intervals of its 301 adaptive charts.
It held the original rational maps and height 100,000 fixed, with two
seconds per tail and a 900-second / 1.5-GiB worker cap. All **301 interval
unions now reach 100,000**. Runtime was 468.36 seconds. The
[coverage proof v2](../../artifacts/generated-results/elliptic-curves/new_rank26_tail_coverage_v2.json)
replays every map, point admission and contiguous interval union. Coverage
completeness relies on the pinned sieve executions, not an independent
second enumeration. The separate
[whole retained cloud proof](../../artifacts/generated-results/elliptic-curves/new_rank26_all_retained_mod2_v1.json)
checks 2,338 unique points up to sign through prime 997 and still proves
rank at least 26. Its input explicitly concatenates point clouds, without
pretending the 645 source charts form a single chronological state history.

**Verifier correction:** v1 reused the parent-chart index variable in its
inner admission loop. This mislabeled 27 coverage rows, although the chart,
interval, membership, rank and total-count checks used the correct inputs.
The frozen v1 source and output are preserved. The
[versioned checker](../cas/replay_new_rank26_tails_v2.py) uses a distinct
admission index; fresh replay passed, and v2 labels are exactly 0 through 300.
The old unfinished-tail ledger is historical evidence, superseded by v2
coverage rather than overwritten.

## Native `11952` coordinate control

The current metric/GMP policy was calibrated on the known ICARM curve 12,
using only its 17 reconstructed generic points. Public exceptional points
were excluded from the search. This is a known-curve control, not an
unknown-curve experiment. The input is exactly transported by scaling six
from the previously frozen generic-only boundary. All 49 centres were fixed
from the explicit 384-bit generic-point height Gram before searching.

The first 43 and additional six classes both end at lower bound **25**.
After completing all 49 missing denominator tails, the bound is still 25.
All initial and tail histories, coverage unions and complete clouds replay.
Thus incomplete denominator coverage no longer explains this control's
shortfall at the declared metric boxes.

The older blind PARI generic43 discoveries independently certify **27**.
The old and new generic43 masks agree, and all 43 centre words agree up to
sign. A retrospective exact audit of both signs of the 28 retained old
candidates on all 49 current boxes makes 2,744 comparisons: 20 observations
are visible and recorded; 2,724 are outside the finite boxes. Twenty of the
28 candidates are visible somewhere; eight are outside every metric box.
There are no visible-but-missing discrepancies. Combining the old and new
clouds still certifies 27, not an additional direction.

This isolates a concrete point-visibility limitation of the current boxes.
It does not prove that those eight representatives exhaust missing quotient
directions, or that coordinate reduction will improve unseen curves. The
[visibility proof](../../artifacts/generated-results/elliptic-curves/native11952_metric_visibility_v1.json)
is explicitly retrospective and kept out of prospective centre selection.

## Controlled coordinate replacement

All 49 native control centres were then held fixed while PARI supplied only
`hyperellminimalmodel`/`hyperellred` horizontal transformations. The current
GMP engine performed every point search. Composition, square ordinate scale
and binary-quartic identities are replayed over the rationals; the checker
does not call PARI. All 49 boxes completed to height 100,000 in 264.29 seconds
under a 500-second / 1.5-GiB cap. The original 43 recover lower bound **27**;
the additional six add no further certified direction. The complete-cloud
proof also certifies 27. Thus the coordinate change recovers the two
directions missed by the metric policy on this control with the same engine.
This is a calibration result, not a new curve or universal policy comparison.

The first map-history checker failed before replay because it compared
`MWState` rational strings directly with parsed fractions. Its failure log
and source remain intact. Version 2 parses the stored coordinates before
equality; its full exact replay passes. Independent cloud build/check also
pass. No failed check was converted into a pass or overwritten.

The same coordinate policy was then applied to the new rank-at-least-26
curve's original 43 generic centres, beginning from its certified 26-point
subgroup. All 43 boxes completed to height 100,000 in 301.20 seconds under
a 500-second / 1.5-GiB cap. All exact map/admission histories and the complete
573-point cloud replay, still with lower bound 26. The extra coordinate
policy recovered more retained points but no further certified independent
direction. This finite miss does not bound the curve's rank.

## Fixed-box backend comparison

A separate bounded comparison ran PARI `hyperellratpoints` on the exact
49 reduced `P,Q` models already searched by GMP. Each call had a three-second
cap, with a 180-second / 1-GiB supervisor. All calls completed. The complete
sets of **72 affine square coordinates** agree exactly across all 49 boxes;
every returned coordinate is checked against the final integer quartic.
The standalone comparison replay also passes. Infinity is handled separately
by the shared pointed engine and is outside this affine-set comparison.

PARI reports 17.314 seconds of search CPU time; the retained GMP workers
report 239.474 seconds. Aggregate PARI call wall time was 18.030 seconds,
versus 243.696 seconds for the shared GMP search calls, which also perform
point transport. These clocks and surrounding operations differ, so this
is an observed efficiency gap on one control, not a universal speed ratio.
The [comparison certificate](../../artifacts/generated-results/elliptic-curves/native11952_engine_comparison_v1.json)
retains every exact hit set. This supports testing a separately versioned
PARI backend while retaining exact shared chart maps, point witnesses and
independence checks; it does not change any frozen worker or default.

## Reproduction

```sh
python3 elliptic-curves/cas/export_r17_exact_maximum_parity_classes.py \
  --check artifacts/generated-results/elliptic-curves/r17_exact_maximum_parity_classes_v1.json
python3 elliptic-curves/cas/audit_recorded_point_mod2_rank_v2.py \
  --check artifacts/generated-results/elliptic-curves/new_rank26_all_retained_mod2_v1.json
python3 elliptic-curves/cas/audit_native11952_visibility.py --check
```

Raw protocols, checkpoints and verification limits are retained under
`artifacts/local/elliptic-curves/`: `r17-norm12-exact-minima-v1`,
`r17-norm12-exact-minima-remaining4-v1`, `r17-omitted-generic-classes-v1`,
`new-rank26-fixed-tails-v1`, `native11952-metric49-control-v1`,
`native11952-metric49-tails-v1`, and `native11952-control-union-v1`.

The [evidence supplement](../../artifacts/generated-results/elliptic-curves/exact_parity_coordinate_evidence_v1.json) records its two pinned base archives and all additional sources, raw records and proofs. Extract the bases in listed order, then the supplement. The [isolated replayer](../cas/verify_exact_parity_coordinate_bundle.py) verifies member hashes and reruns the new mathematical checks without launching point searches.
