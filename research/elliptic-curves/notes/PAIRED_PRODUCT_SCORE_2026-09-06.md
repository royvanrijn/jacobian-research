# Paired product-score experiment and a visibility countercheck

The alternative product score produced no stronger curve in this fixed trial.
Both twelve-curve arms have eleven certified lower bounds of 17 and one of 18;
the rank-at-least-18 curve, `11952-005`, belongs to both arms. The two arms share
two curves, leaving 22 distinct curves and 988 completed point-search boxes.
All retained-point clouds give the same bounds modulo 2, 3 and 5. These are
lower bounds from finite searches, not measurements of the curves' true ranks.
The [aggregate certificate](../../artifacts/generated-results/elliptic-curves/product22_comparison_v1.json)
binds the selection, exact proofs, complete clouds and visibility audits.

## Fixed comparison

The input is the already saved 6,144-address higher-parameter population and
its finite-field traces. No new parameter scan or trace computation was used.
The product score sums `log((p+1-a_p)/p)` over nonsingular displayed reductions;
each binary logarithm contribution is rounded to units of 10^-12 before summing.
This numerical ordering is reproducible, not an exact transcendental comparison.
The original S1 score is recomputed from the same data.

The 24 previously exposed addresses are excluded uniformly. For each of six
families, each score selects its top two remaining addresses using primes
through 32,749. Primes from 32,771 through 65,521 are validation only and cannot
break ties. Overlaps are merged without refill. The original 24 outcomes were
known to the researcher, but their points and ranks are not read by selection.
The entire input pool was already truncated by short-prime S1: this experiment
cannot recover candidates discarded from the original 122,368,792 addresses.

Every selected curve receives the same generic-17-only policy: all 43 or 49
declared generic maximum-parity charts, height 125,000, ten seconds per chart.
All 988 boxes completed. Exact history, rational-map, raw-point, full-cloud and
independence replays passed. Post-search Q-isomorphism checks find no matches
among the 22 curves, the pinned 593-curve ICARM snapshot, or 426 previous
address-equations. These low bounds do not qualify for the high-rank inventory;
the inventory remains 100 curves, including six with certified bound 27.

This is a null result for improved bounded detection in this paired sample.
It is not evidence that product scoring is generally ineffective or that the
two policies select identical true-rank distributions. The product formula is
established prior art, already tested on a smaller ordinary-fibre panel here;
the new contribution is this higher-parameter comparison and its recorded
outcomes. See [Elkies--Klagsbrun, sections 2 and 7](https://arxiv.org/pdf/2003.00077)
and [Watkins et al., section 5](https://www.dpmms.cam.ac.uk/~taf1000/papers/rankcongr.pdf).

## Direct representative visibility is not direction recovery

An exact retrospective audit checks both signs of all 17 generic sections on
the initial completed charts of the old and higher 24-curve cohorts: 73,440
point/chart observations. No known point lying within completed coverage was
omitted from the recorded output. Seventeen higher-cohort curves have no
directly visible generic representative; ten old-cohort curves do likewise.
These points already entered centre construction, so this is not a masked
control or an independent test of withheld-direction recovery.

The essential countercheck is the known curve where the same 49-chart,
height-125,000 method independently recovered 28 directions from 17 generic
seeds. Its 1,666 exact generic-point observations also have **zero visible
original generic representatives**, with no within-box omission. Thus an
original basis representative outside every box does not imply failure to
detect exceptional directions. Other representatives and combinations matter.
The [control certificate](../../artifacts/generated-results/elliptic-curves/native28_generic_visibility_v1.json)
preserves this negative result for the proposed diagnostic. It does not prove
uniform sensitivity, saturation or absence of missing directions on new curves.

## Replay

`report_product22_comparison.py --check` checks the aggregate and its bindings.
`compare_higher_r17_product_score.py replay` recomputes all saved scores and
selection. `certify_product22_r17_results.py --check` checks the 22 point proofs
and Q-isomorphism comparisons. `replay_product22_geometry.py` checks all 988
maps and raw-point provenance. The per-curve full clouds have independent
modulo-2 and modulo-3/5 checkers. `audit_generic_point_box_visibility.py replay`
and `audit_native28_generic_visibility.py --check` reproduce the visibility
observations without a new point search.
