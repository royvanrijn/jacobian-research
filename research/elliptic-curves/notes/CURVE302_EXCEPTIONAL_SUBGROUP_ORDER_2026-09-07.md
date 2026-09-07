# Curve 302: exhaustive exceptional-subgroup order diagnostic

This is a retrospective finite experiment on the completed 302 data. It does
not alter the certified lower bound `rank E_302(Q) >= 31`, select a future
centre, or establish a point-search runtime law.

## Fixed diagnostic universe

The 14 known exceptional directions are used only as a fixed displayed
complement to the certified generic rank-17 subgroup. The prior visibility
artifact proves that this is a direct integral `M17 -> M24 -> D` presentation.
Consequently there are exactly `2^14 = 16,384` displayed subgroups

\[
 M_{17+S}=M_{17}+\sum_{i\in S}\mathbb ZR_i.
\]

Every state was measured in one entrywise-`10^6`-rounded, 384-bit canonical
height Gram metric. For each remaining direction the stage-local score is a
nearest-plane candidate followed by deterministic exact coordinate descent for

\[
 4\min_{q\in M}\widehat h_{\rm rounded}(R_i-q/2).
\]

It is an upper bound for the rounded-metric CVP problem, not an exact CVP
minimum, canonical-height theorem, or pointed-quartic coordinate minimum. A
separate continuous projection lower bound is retained for calibration.

The cumulative score is the minimum over every candidate obtained in every
contained subgroup. It is therefore a persistent finite atlas: it cannot rise
under subgroup inclusion even if a fresh nearest-plane calculation gets worse.

## Order results

The global minimax path over all fourteen axes has maximum score
`123,778,523 / (4·10^6) = 30.94463075`. The minimum-total path has total
`1,375,488,055 / (4·10^6) = 343.87201375`.

The first seven historical gains are certified only as an `M17 -> M24`
seven-dimensional block. No source attests an internal order on the fixed
diagnostic coordinates, so the experiment does **not** fabricate one. It
checks all `7! = 5,040` compatible orders before appending the attested
`M24 -> M31` tail

\[
06,04,02,01,05,03,07.
\]

For that fixed start, the historical tail has bottleneck
`135,591,028 / (4·10^6) = 33.897757`; the tail minimax optimum is
`102,542,623 / (4·10^6) = 25.63565575`. Its total falls from
`669,602,723 / (4·10^6) = 167.40068075` to
`629,602,582 / (4·10^6) = 157.4006455` under the tail minimum-total order.
Thus the historical tail was viable but not close to this surrogate's best
bottleneck order. This diagnoses available half-lattice geometry, not the
actual time required by the historical quartic searches.

At the root, the largest retrospective unlock values are obtained by
`recovered-strict-03`, `residual-strict-04`, `residual-strict-05`, and
`residual-strict-03` (in that order). Here unlock uses the declared finite
regularization `log(1+C)`, so a genuine half-lattice endpoint with `C=0`
remains meaningful.

## Local versus cumulative check

Across all 745,472 comparable one-direction enlargement edges, 88,832 fresh
nearest-plane scores rise. The retained finite-atlas score rises on zero
edges. This is exactly the requested separation: the raw rises are reducer or
scheduling instability, while the retained curve obeys the inclusion law.

The autonomous V1 subgroups were also expressed in the fixed 31-point basis.
Each inferred integral word was checked with the exact elliptic group law;
they are not relabelled as a permutation of the fourteen coordinate axes. By
its ranks `19,22,23,26,28`, the number of fixed directions that become exact
half-lattice endpoints is `1,3,4,5,7`; at rank 28 these include
`residual-strict-04`. This is strong
structural evidence for the half-lattice interpretation, but not evidence
that V1 prospectively knew those labels.

## Reproducibility and boundary

- [Complete 16,384-state ledger](../../artifacts/generated-results/elliptic-curves/curve302_exceptional_subgroup_landscape_v1.json)
- [Order report](../../artifacts/generated-results/elliptic-curves/curve302_exceptional_subgroup_landscape_report_v1.json)
- [Path table](../../artifacts/generated-results/elliptic-curves/curve302_exceptional_subgroup_landscape_v1.tsv)
- [Raw-versus-cumulative figure](../../artifacts/generated-results/elliptic-curves/curve302_exceptional_subgroup_landscape_v1.svg)
- [Core enumerator](../cas/audit_curve302_exceptional_subgroup_landscape.sage)
- [Independent report/replay](../cas/report_curve302_exceptional_subgroup_landscape.py)

Replay with:

```bash
cd research
sage -python elliptic-curves/cas/audit_curve302_exceptional_subgroup_landscape.sage --check
python3 elliptic-curves/cas/report_curve302_exceptional_subgroup_landscape.py --check
```

The targets, their directions, and the autonomous point words are all known
data confined to this retrospective analysis. A prospective beam search must
rank geometric orbit landscapes without reading this diagnostic basis,
exceptional labels, point coordinates, or catalogue ranks.
