# Search changes supported by the recent parallel results

Follow-up: the [fixed own27 point-search control](INVENTORY188_EXCEPTIONAL_DIRECTION_RECOVERY_2026-09-07.md) recovers the known28 direction on chart5. The representative misses below remain correct, but do not imply failure to recover their quotient direction. The [outer/native exposure follow-up](RETAINED_OUTER_AND_NATIVE_EXPOSURE_2026-09-07.md) completes196 boxes without a rank gain and removes a global-minimalization preparation bottleneck through a separately calibrated factor-free mapping policy.

The review covers the completed corrected/stratified/near-finalist searches,
the rank-jump commits through `e16d3dc1`, and the curve302 and class-span
results included in `1845b0f4`. It produces a reusable exact norm preflight,
tested on retained inputs. It supplies **no new arithmetic-score term, curve
exclusion, parameter sweep, or rank-record claim**. The supplement below
extends the review through `73602fb5`.

## Newer commits: useful construction gates, no new score terms

`1f7f8632` proves a [fixed-incidence six-direction solubility switch](../rank-jump/FIXED_INCIDENCE_SIX_DIRECTION_SOLUBILITY_SWITCH.md).
The same strict incidence block can contain rational directions on one curve
and a nondegenerate obstructed six-dimensional block on its twist. This
supports separating incidence from solubility. It does not supply an
equation-only extractor for the missing additional block in our candidate
families, and its conditional exact-rank statements remain conditional.

`55b3ed8d` proves that [keeping the complete cubic field constant](../rank-jump/FIXED_CUBIC_TRANSFER_REQUIRES_HIGH_GENUS.md)
requires a degree-six carrier of genus31 in the three tested R17 families
and genus28 in the two tested MW16 families. All eight retained high/low
pairs have different discriminant squareclasses. A rational or elliptic
base change preserving the full cubic therefore cannot be the proposed
transfer mechanism in these families. This does not exclude a carrier
transporting fewer data or allowing the cubic field to vary.

`2751b3d2` bounds the [everywhere geometrically locally soluble generic pool](../rank-jump/LARGE_JUMPS_EXCEED_THE_GENERIC_SELMER_POOL.md)
by20−c, where c is the reducible-fibre root-lattice rank. In the seven
verified presentations, the generic17/16-point subgroup leaves capacity
at most three for additional classes in this particular pool. A proposed
global block explaining specialized rank R needs geometric local obstruction
rank at least R−(20−c). If its obstruction support consists of g good and
b bad geometric base places, it must satisfy2g+b ≥ R−(20−c).
These are parameter-direction conditions, not rational-prime score features
or upper bounds on specialized Mordell–Weil rank.

Both narrow portable verifiers pass in this review: five constant-field
geometry rows, and seven Selmer-capacity presentations with sixteen exact
specialization/generic checks. The cohomological statements still depend on
their written mathematical proofs. The production consequence is to apply
these necessary capacity checks to proposed class/transfer constructions
before expensive searches. They do not justify filtering individual retained
fibres or changing the completed score-stratified comparison. The missing
knowledge remains a specialization-dependent additional-class construction
with a rational-solubility mechanism; a new score term is not established.

The subsequent `73602fb5` makes that restricted class carrier explicit as
[rational2-torsion on the cubic root Jacobian](../rank-jump/ROOT_CURVE_TORSION_AND_REAL_CAPACITY.md),
with one possible node-parity bit for MW16. Exact real topology sharpens the
additional generic-pool cap from three to **two** in its six verified panel
presentations: total capacity19 for R17 and18 for MW16. These sharper totals
replace20/19 in the necessary obstruction-rank test for those presentations.
The three-prime Frobenius parity test gives no improvement; its unipotent
ambiguity is certified, so more of the same trace sampling is not justified
by this result. The independent finite-root, Sturm-isolation and topology
replay passes. This supplies a concrete object for future class construction,
but neither an additional class basis nor rational points on the elliptic
covers. The real-topology bound is constant within a family and cannot rank
its specializations. No score term or candidate exclusion follows.

## Applied improvement: reject impossible norm words before expensive work

The [ramification lemma](../rank-jump/FRESH_NORM_PROJECTION_RAMIFICATION_GATE.md)
in `d602fd60` and the [complete retained dictionary audit](../rank-jump/RETAINED_NORM_RELATIONS_DO_NOT_YET_SUPPLY_THE_BLOCK.md)
in `e16d3dc1` can be used earlier in class construction.

[`norm_ramification.py`](../cas/research_runtime/norm_ramification.py) accepts
a separable monic rational cubic and a finite dictionary of nonzero elements
of degree less than three, each with nonzero norm. It recomputes exact norms
as determinants. For each generator alpha_i, it removes from its norm every
prime factor shared with the forbidden support or another active norm.
Forbidden support includes 2, the cubic discriminant, coefficient
denominators and polynomial contents. A nonsquare remaining integer proves
that coefficient i must be zero in **every unramified product** of the norm
projections N(alpha_j)alpha_j. Such coordinates are removed and the procedure
repeats until stable.

The reason is exact: at some odd good prime the selected norm has odd
valuation and all other active norms are units. The selected polynomial
cannot vanish in every component of the separable cubic residue algebra.
On a component where it is a unit, multiplication by its scalar norm
introduces odd valuation. Other generators cannot cancel it. Induction
justifies subsequent peeling rounds. No prime factorization is needed.

The filter preserves possible cancellation between generators sharing
support. Two identical ramified generators, for example, must remain
unresolved because their product is a square. Separately rejecting every
ramified basis vector would give a false conclusion. A nonempty residual
dictionary is **UNKNOWN**, not a certificate of unramifiedness, class
independence, Selmer incidence or rational solubility.

For elliptic Kummer work, generic correction cannot cancel the obstruction
provided all bad elliptic places are in the excluded support. The retained
integration audit checks that inclusion for every case below. An arbitrary
caller must establish the same condition before making that interpretation;
the bare cubic API only proves its stated unramifiedness restriction.

The [frozen protocol](../../artifacts/generated-results/elliptic-curves/retained_norm_preflight_protocol_v1.json)
selects all 132 original panel generators and all 296 fixed-box reference
generators. The projected worker input contains only cubics and element
coefficients. The dictionaries were already studied: this is a retrospective
equivalence/calibration test, not a blind discrimination trial.

| Retained dictionary | Generators | Forced-zero coefficients | Remaining coefficient cap |
|---|---:|---:|---:|
| Eleven original panel dictionaries | 132 | 132 | 0 |
| Fixed-box MW16-05, t=3/17 reference | 296 | 296 | 0 |
| Total | 428 | 428 | 0 |

All cases finish in a single peeling round. The supervised complete build
takes **0.1273 seconds**, its replay **0.1271 seconds**, and the independent
Sage verification **0.6413 seconds**, including process startup. These are
single local measurements, not a universal runtime or whole-search speedup.
The frozen cap is 60 seconds per invocation. A display-only wrapper lookup
used an incorrect timing-key name after the successful build; the completed
worker log and supervisor record remain intact and were read correctly
afterward. No arithmetic rerun was needed to repair that display.

The [result](../../artifacts/generated-results/elliptic-curves/retained_norm_preflight_v1.json)
and [independent replay](../../artifacts/generated-results/elliptic-curves/retained_norm_preflight_sage_v1.json)
agree. Sage recomputes all 428 norms by polynomial resultants and checks all
428 isolation witnesses with a separate product-gcd calculation and exact
square-root inequalities. Six focused tests cover shared-support
cancellation, dependent peeling rounds, repeated prime powers, rational
scaling, unresolved results and invalid arithmetic.

**Use this before local-support factorization or CT work on future norm
dictionaries.** Apply it to the whole proposed dictionary; adding generators
can create new cancellations, so old forced-zero conclusions cannot simply
be retained after expanding that dictionary. The filter does not require
maximal-order initialization. It rejects constructions, never the underlying
elliptic curve. The known-soluble +6 reference demonstrates the distinction:
its genuine additional directions exist despite this dictionary's zero
capacity. No frozen worker, score or point budget is edited.

## Other findings and their search consequences

| Recent result | Supported use | Still missing |
|---|---|---|
| `61a5d666`: conditionally exact +6 reference | Calibrate an additional-class extractor against six genuine extra strict directions; retain the GRH qualification on the upper bound | A point-independent extractor reaching those six directions and a prospective solubility criterion |
| `69163e44`, `85bbb451`: strict boundary bounds | Track additional class capacity relative to the generic subgroup; only completed, applicable upper certificates can exclude rank targets | Most localized class dimensions; an observed-zero point search is still censored |
| `874c875f`: symbolic discriminants | Reuse the verified family factors in local arithmetic; stop repeating this failed split on the five frozen cofactors | Integer factorization and missing boundary data; irreducible polynomial values need not be prime |
| General class-span GRH machinery | Use audited principal relations and interval character exclusion for conditional upper bounds without demanding a complete class presentation | New curves' own relation/prime coverage and local Selmer corrections; conditional bounds cannot become unconditional production vetoes |
| Curve302 full MW9 in `1845b0f4` | Stop looking for generic MW17–20 by changing fibrations on this constructed K3; use the saturated basis for any work on it | A different high-rank parent; the same rank31 fibre on this MW9 surface is not a new curve |
| Completed score-strata experiment | Retain strong-score selection as the supported baseline; its top arm supplied 10 directions versus 1 and 0 with matched exposure | A stable production optimum; one +9 curve dominates the result |
| Completed near-finalist60 | Retained ranks7–13 are productive: 16 directions and one new rank-at-least23 curve, without a new scan | A matched comparison proving this cutoff better than the top six; the two trials are not interchangeable replicates |
| Known28 chart control | Treat exceptional-direction coverage as an unresolved point-search bottleneck before increasing population size | A chart policy with demonstrated improvement on held-out directions; own27 geometry worsened both tested public representatives |

The [class-span machinery](CLASS_SPAN_GRH_MACHINERY.md) is useful proof
infrastructure. It must not be confused with constructing new classes: many
independent principal relations can improve an upper bound while their norm
projections supply no unramified excess. The new preflight makes that
distinction executable. Generic governing-field degree and inherited
CT-switch size did not distinguish the completed controls; neither is added
to the production score.

## A stronger calibration target from the public28 reproduction

The parallel strict-boundary panel retained rank27 for curve188. Its exact
equation is identical to the [independently reproduced public28 model](CURRENT_CATALOGUE_AND_PUBLIC28_2026-09-07.md).
The [integration audit](../../artifacts/generated-results/elliptic-curves/search_result_integration_v1.json)
checks that equality and replays the finite28-point independence proof.
The existing boundary certificate gives m=17, k=0 and a=4, hence

```
additional strict rational dimension >= 28 - 17 - 4 = 7,
localized class dimension c_S >= 7.
```

This strengthens the previous necessary lower bound six to **seven**.
It is a deduction after joining public points, not an independently measured
class-group feature. Its class upper bound and exact curve rank remain
UNKNOWN. The historical panel and its frozen labels are preserved.

## Replay and next decision

```sh
python3 -m unittest discover -s elliptic-curves/tests -p test_norm_ramification.py
python3 elliptic-curves/cas/audit_retained_norm_preflight.py check
sage -python elliptic-curves/cas/verify_retained_norm_preflight.sage --check
python3 elliptic-curves/cas/audit_search_result_integration.py --check
```

The next useful point-search policy change needs a frozen visibility test
with identical completed exposure. The next class-construction change needs
a dictionary surviving the complete parity and generic-dependence gates.
Neither result calls for another parameter scan. Validation primes remain
separate, and the completed retained-score comparison is unchanged.
