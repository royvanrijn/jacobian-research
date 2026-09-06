# A new rank-at-least-25 curve from a model-scaling condition

The compact `11952` family at `t = 102/1525` gives inventory curve
`new-20260906-101`, with 25 exactly certified independent rational points.
Its globally minimal integral equation is

```
y^2 + x*y = x^3
 - 1786694850006932466619529831652324354433346080*x
 + 28897710108108857599038628725248039072984480888228845021215611110400.
```

The [minimal-model proof](../../artifacts/generated-results/elliptic-curves/scaled13_24_rank25_minimal_proof_v1.json)
contains the exact model transport, minimality witness and independent-point
certificate. The largest coefficient has 225 bits. The
[Sage export](../../artifacts/generated-results/elliptic-curves/new_scaled13_rank25_curve_11952.sage)
loads the equation and all 25 points. No rational-isomorphism match occurs in
the pinned 593-equation catalogue or the 472 earlier measured address-equations.
This establishes catalogue absence, not universal novelty or exact rank.

The [V12 inventory](../../artifacts/generated-results/elliptic-curves/new_high_rank_curve_index_v12.json)
now contains 101 distinct curves: six lower bounds 27, eleven 26, twenty-two
25, twenty-four 24, twenty 23 and eighteen 22. Every point proof and the CSV
replay, with previous curve IDs preserved. The new curve has only its own
rational preimage in the twelve recorded presentations. The
[incidence aggregate](../../artifacts/generated-results/elliptic-curves/inventory101_incidence_v1.json)
now covers 1,212 pairs: 1,090 exclusions and 122 certified preimages. The
21 extra presentations remain the previously proved duplicate R17 subgroup.
The strongest new lower bound is still 27; no new near-record is established.

## Two distinct local-feature tests

The preceding [product-first experiment](PRODUCT_FIRST_RETENTION_2026-09-06.md)
identified 766 retained models with a removable scaling by 13. The original
higher-parameter cohort contained eight, while the two product-score cohorts
contained one and zero. This observation motivated a separate finite test;
it did not establish a causal effect or a predictor.

The first new protocol conditions the original 6,144 saved S1 candidates on
that exact scaling identity, excludes all 70 addresses in the three frozen
prior point cohorts, and retains four per family by unchanged S1 selection
scores. There are 454 eligible addresses. No new traces, parameter population,
public rank labels or outcome-dependent refill enter selection. Validation
primes remain excluded from ties. All 1,080 generic-17-only boxes complete at
height 125,000 with ten seconds per chart. The certified lower bounds are
13 at 17, four at 18, four at 19, two at 20, and the new curve at 25.

The second protocol tests the split-multiplicative-prime heuristic discussed
by [Elkies--Klagsbrun, section 7, question 2](https://arxiv.org/pdf/2003.00077).
The earlier product experiments lacked this term. A fixed `c = 7/5` adds
`round(log(7*(p-1)/(5*p))*10^12)` at split nodes to the saved product score.
Primes through 32,749 select; primes through 65,521 above that bound only
validate. The 10,482-address retained union and existing traces are reused.
All 94 addresses in four already frozen point protocols are excluded; the
concurrent scaling cohort's outcomes are not read. Four candidates per family
receive the same point policy. All 1,080 boxes complete: 23 lower bounds 17
and one 18. This test adds no high-rank inventory curve.

At each score prime at least 5 dividing the short-model discriminant, the
local calculation repeatedly divides integral powers `p^4` and `p^6` out of
A and B. A nonsingular resulting display would fail the protocol rather
than silently omit its trace. None occurs. For a node, the double root
`r = -3B/(2A)` satisfies `A = -3r^2`, `B = 2r^3`, and its tangents split
exactly when `3r` is a nonzero square. Every residue and Euler-criterion
calculation is retained. Cusps receive no bonus. No Tamagawa number, conductor
or analysis at primes 2 and 3 is asserted.

All 48 measured curves are internally distinct within their cohorts and
unmatched in the catalogue and respective prior equation sets. All admission
histories, chart geometry and complete-cloud proofs modulo 2, 3 and 5 replay.
The [aggregate](../../artifacts/generated-results/elliptic-curves/local_feature_experiments_v1.json)
binds the two policies and their different outcomes. These lower bounds are
not true-rank distributions or a demonstration that either selector is better.

## The complete 301-box follow-up

After the new 25-point proof, a separate frozen protocol pairs 301 generic
census labels with nonzero combinations of the curve's eight discovered
directions beyond its generic seventeen. Every box completes at the same
height and per-chart limit. The 350 initial-plus-adaptive chart inputs yield
3,260 points up to sign; their full cloud still has certified lower bound 25
modulo 2, 3 and 5. The original minimal-model proof remains valid.

This is a finite search of chosen representatives. In particular, pairing
generic labels cyclically with quotient words does not enumerate all parity
classes in the specialized 25-point subgroup, prove specialized covering
optimality, or exclude a 26th point. No further automatic height escalation
follows this bounded null.

## Reproduction

The [evidence manifest](../../artifacts/generated-results/elliptic-curves/local_feature_evidence_v1.json)
contains the new inputs and names nineteen pinned base archives. Extract the
bases in listed order and then the supplement. `verify_local_feature_bundle.py`
replays the two selectors, both cohorts' exact maps and point proofs, all
48 complete clouds modulo 2, 3 and 5, the adaptive union and geometry,
minimality, Sage export, V12 inventory and twelve new incidence pairs.
It launches no new point search. The archive's isolated replay report is
linked from the artifact README after completion.
