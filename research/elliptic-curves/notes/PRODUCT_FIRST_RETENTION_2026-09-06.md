# Product-first retention reaches discarded candidates but adds no high-rank curve

Changing the initial selector retains 4,338 addresses that the earlier S1
screen discarded. A fixed 24-curve point cohort from these newly retained
addresses completes all 1,080 boxes. Twenty-three curves have certified lower
bound 17 and one, `103b2-543`, has bound 18. All full clouds give the same
bounds modulo 2, 3 and 5. There is no addition to the 100-curve high-rank
inventory and no new near-record result.

The [aggregate certificate](../../artifacts/generated-results/elliptic-curves/product_first_experiment_v1.json)
binds the complete retention experiment, score extensions, point proofs and
local-scaling diagnostic. It also verifies agreement between two product-score
implementations on every overlapping address.

## The retention-stage gap is tested separately

The [earlier paired experiment](PAIRED_PRODUCT_SCORE_2026-09-06.md) could only
reorder 6,144 addresses already retained by short-prime S1. It could not test
whether product scoring would retain useful discarded addresses. This new
experiment uses precisely the same twelve signed denominator slices and
122,368,792 primitive addresses at parameter height 32,768. It introduces no
larger height or new slice.

Only the cached table contributions change, to
`round(log((p+1-a_p)/p)*10^12)` at good displayed reductions. The binary,
primitive-address enumeration and tie order remain fixed. Before the large
calls, all 6,912 primitive scores and complete ordered outputs in twelve
31-by-29 boxes pass exact checks against canonical residue traces. The
returned 6,144 large-box scores also replay. The retained pool overlaps the
old one in 1,806 addresses, leaving 4,338 newly retained addresses.

The extension reuses the existing 10,796,268 trace values for the overlaps
and computes 25,932,564 for the newly retained addresses. A six-curve
benchmark first passes 48 independent character sums. Selection uses the
product score through 32,749; validation primes through 65,521 do not break
ties. The two implementations agree exactly on selection scores, validation
scores and selection-band good counts at all 1,806 overlapping addresses.

Four candidates per family are selected only from addresses absent from the
old 6,144 pool, without reading point results or catalogue labels. All receive
the same generic-17-only 43/49-chart policy, height 125,000 and ten seconds per
chart. Every box, admission history, exact map and full-cloud proof completes
and replays. Post-search Q-isomorphism checks find no matches within the 24
curves, the pinned 593-equation catalogue, or 448 previous address-equations.
The low certified bounds do not meet the inventory threshold of 22.

Retention and replay took 52.996 and 13.919 seconds. Extended traces and replay
took 433.927 and 38.027 seconds, reusing cached overlaps. The point cohort and
its history/cloud verification took 553.909 and 62.123 seconds. No automatic
point refill or increased point budget followed the null result.

## A local-model audit and a descriptive difference

Both score policies omit singular displayed reductions. On all 10,482
distinct address-models in the two retained pools, an exact bounded audit
checks p-power scalings for score primes from 5 through 65,521. If p^4 divides
A and p^6 divides B in `y^2=x^3+Ax+B`, dividing by those powers gives an
explicitly Q-isomorphic model. The audit checks its discriminant modulo p.

There are 766 such removable scalings, all by 13 once. Every scaled display
remains singular at 13. Thus this audit restores no omitted good-reduction
contribution. It does not classify all bad reductions, audit primes 2 and 3,
or establish global minimality. The [scaling certificate](../../artifacts/generated-results/elliptic-curves/higher_displayed_reduction_scalings_v1.json)
contains the exact coefficient identities and residues.

The original higher-parameter 24-curve cohort selected eight of these
13-scaled models, the paired 22-curve cohort selected one, and the new
product-first cohort selected none. Scaling by 13 divides the displayed
coefficients A and B by 13^4 and 13^6. This is a concrete model-size difference
between the cohorts, not a causal explanation of the observed point bounds.
Model cost and bad-prime information remain possible diagnostic directions;
no predictor or new selection policy is promoted by this comparison.

[Elkies--Klagsbrun, sections 6--7](https://arxiv.org/pdf/2003.00077) discuss the
tradeoff between initial cutoffs and prime bounds, skewed parameter regions,
and a heuristic bonus at split multiplicative primes. That bonus is not the
unmodified product score tested here. The source does not give a general
calibration theorem for these curves.

## Replay

`replay_product_first_population_portable.py` checks the retained population
with recorded command paths bound to the original checkout and data read from
the extracted workspace. `extend_product_first_higher.py replay` checks all
trace rosters, cached inputs and final ordering. `certify_productfirst24_r17_results.py
--check` checks the point proofs and Q-isomorphism comparisons;
`replay_productfirst24_geometry.py` checks all 1,080 chart maps and point-cloud
provenance. The aggregate links the per-curve modulo-2 and modulo-3/5 proofs.
`audit_displayed_reduction_scalings.py replay` checks the local-model audit.
