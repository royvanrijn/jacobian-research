# Search changes supported by the recent parallel results

Follow-up: the [fixed own27 point-search control](INVENTORY188_EXCEPTIONAL_DIRECTION_RECOVERY_2026-09-07.md) recovers the known28 direction on chart5. The representative misses below remain correct, but do not imply failure to recover their quotient direction. The [outer/native exposure follow-up](RETAINED_OUTER_AND_NATIVE_EXPOSURE_2026-09-07.md) completes196 boxes without a rank gain and removes a global-minimalization preparation bottleneck through a separately calibrated factor-free mapping policy.

The review covers the completed corrected/stratified/near-finalist searches,
the rank-jump commits through `e16d3dc1`, and the curve302 and class-span
results included in `1845b0f4`. It produces a reusable exact norm preflight,
tested on retained inputs. It supplies **no new arithmetic-score term, curve
exclusion, parameter sweep, or rank-record claim**. The supplement below
extends the review through `e575d0e5`.

## Newer commits: useful construction gates, no new score terms

The follow-up review through `e575d0e5` independently replays three further
restrictions. These refine constructor requirements while preserving the
completed score experiments and their separate validation primes.

- `b6bcf2bd`: the [standard S3 class block](../rank-jump/STRICT_CLASS_CREATION_IS_A_STANDARD_S3_BLOCK.md)
  carries the strict elliptic classes. The trivial quadratic-resolvent genus
  component contributes none. The subgroup cochain checks, 1,933 Hilbert
  bits and four small certified class groups replay. Increasing resolvent
  genus rank alone is therefore not a constructor for the missing elliptic
  directions; the cubic standard component remains to be found.
- `085690b6`: the [relation-root construction](../rank-jump/RELATION_ROOTS_CREATE_RAMIFICATION_NOT_THE_MISSING_BLOCK.md)
  leaves the old span but has seven odd ideal valuations. The entire
  4,134-element retained dictionary plus generic corrections cannot remove
  them; independent replay checks 33,864 lattice valuations. Apply its
  necessary modulo-four valuation condition before treating a square-root
  norm relation as a new unramified class. This excludes the retained
  construction, not its elliptic curve or every possible dictionary.
- `e575d0e5`: [branch-fibre divisibility](../rank-jump/BRANCH_FIBRES_BOUND_NATIVE_CLASS_CREATION.md)
  bounds each of the 37 retained native quadratic twist supports by three
  new generic directions, for every nonzero rational scalar. Their distinct
  branch kernels bound any product involving two or more supports by two;
  the full fibre product of a pair can add at most eight directions. All
  1,604 finite blocks and 81,804 character bits replay independently. A
  proposed generic gain of nine or more through this specific construction
  requires at least three supports and cover genus at least five. New
  supports with different branch-fibre two-torsion or divisibility require
  their own proofs; none of these bounds excludes specialization gains.

The [six additional parents](MESTRE_DETERMINANT468_PARENTS_2026-09-07.md)
now give an actual change of geometric NS determinant, from production948
to468. Their rational Picard rank18 rules out arithmetic MW17 on those
surfaces and focuses their next construction gate on an explicit useful
MW16 fibration. No such fibration is yet constructed. This parent result
does not use the native-cover capacity bounds outside their proved scope,
and no larger parameter or point sweep follows either result.

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

## Subsequent review through 6b832b2b

The independent verifiers for `fc264181`, `34c83f03` and `6b832b2b` pass.
These sharpen construction requirements without changing the completed
score-strata experiment or launching another parameter scan:

- [The native global pool](../rank-jump/ONE_COMMON_CLASS_FOR_A_LARGE_NATIVE_BLOCK.md)
  has dimension 17 or 18. New generic-rank capacities for one through four
  supports are at most 2,5,10,19, at genera 0,1,5,17. A four-support gain
  at least fourteen requires a common extra class rational on at least
  ten character twists. That class and simultaneous solubility remain
  unknown. The producer's direct capacity check has a tuple-versus-JSON-list
  comparison defect; the read-only
  [`verify_native_common_class_json.py`](../cas/verify_native_common_class_json.py)
  replays its exact object after JSON normalization. Frozen producer files
  and certificates are preserved; no arithmetic values are changed.
- [A degree-four Galois scheme](../rank-jump/A_QUARTIC_GOVERNS_THE_LAST_GLOBAL_CLASS.md)
  reduces the remaining global-class question to an S4 versus S3 action.
  Its defining coefficient polynomial and any extra class remain
  unconstructed. The independent replay checks 136 pairings, the Arf
  invariant and the permutation representation. This is a fixed-root-curve
  incidence result, not a detector for new specialized classes.
- [Fixing a resolvent polynomial](../rank-jump/FIXING_THE_RESOLVENT_POLYNOMIAL_TESTS_SOLUBILITY.md)
  can silently impose point solubility. On the retained MW16-05 control,
  all 63 nonzero classes in a strict six-dimensional block are admitted
  with an arbitrary cubic generator, none with the fixed polynomial
  Z^3+AZ+B, and all with Z^3+AZ-B, although the two polynomials define
  the same cubic field. The independent root-sum quartic, prime-witness
  and CT replay passes. Future incidence constructors must allow arbitrary
  cubic generators and certify field identity, or explicitly label the
  stronger fixed-polynomial condition as a solubility restriction.

These classes are already present in the generic control. No new
specialization direction follows. The concrete parent expansion and seed
repair are recorded in the
[Kihara experiment](KIHARA_PARENT_EXPANSION_AND_SEED_INDEX_2026-09-07.md):
four further Q-distinct parents, exact index-six seed enlargement, smaller
models, and 196 completed boxes with zero rank gains. These bounded outcomes
supply no parent-superiority or specialized-rank upper bound.

The subsequent `0dd1e731` result also passes its independent parity-rank,
norm-one-relation and inclusion replay. It
[exhausts all in-field radical operations on the retained pool](../rank-jump/ALL_IN_FIELD_ROOTS_OF_THE_RETAINED_POOL_ARE_EXHAUSTED.md):
after its one required square root, the 4150-generator group is
2-saturated and its Selmer intersection is precisely the generic16
subgroup. Repeated multiplication, inversion, rational rescaling, norm
projection and integer-root extraction cannot create an extra class from
that fixed dictionary. A future constructor needs an element outside
that radical closure; the theorem does not bound the whole cubic class
or elliptic Selmer group, or exclude the underlying fibre. No production
score changes follow.

## Subsequent review through 9c5c4d3b

The independent checks for `ade6aebf` and `9c5c4d3b` pass. The former
[central-governing-field theorem](../rank-jump/CENTRAL_GOVERNING_FIELDS_CANNOT_CREATE_THE_JUMP_CLASSES.md)
shows that each frozen degree192 pair field contains only its two input
Kummer directions; all three nonzero combinations fail strictness on all
sixteen fibres, with48 explicit local witnesses. Central extensions and
composita of generic pair-governing fields cannot add a Kummer direction:
centrality forces every cocycle into the zero invariant subspace of the
standard S3 module. Three independent finite-prime witnesses also separate
the third generic class field from the first-pair governing field. These
are capacity restrictions on the construction, not fibre-rank bounds.

The latter [four-division result](../rank-jump/FOUR_DIVISION_RAMIFICATION_CANNOT_SUPPLY_THE_JUMP.md)
identifies the sole possible Kummer direction in Q(E[4]) as the derivative
class -Delta f'(theta). On all sixteen frozen fibres it is nonzero but
locally outside the2-Selmer group. Adjoining this field to the generic
governing compositum therefore supplies no additional Selmer class on
the panel. The independent replay covers all67 sign-module subspaces and
the exact local exclusions. It does not assume a full mod4 Galois image.

Future class constructors need a new standard S3 module and independently
verified local admissibility before a solubility test. More central bits,
division-field degree or operations within the exhausted dictionary cannot
stand in for those inputs. Neither commit supplies a production score or
a new rational point, and neither changes the frozen score comparison.

The [retained26 source audit](RETAINED26_SOURCE_GAPS_AND_COMPLETED_EXPOSURE_2026-09-07.md)
separately finds seven source-prefix exposure gaps among eighteen curves,
avoids eleven already-exposed curves and completes98 new boxes on two
fixed training-score selections with zero gains. Five gaps remain
unscheduled. Parent expansion remains a separate need: the twelve labels
of the high-rank inventory all belong to X948, whereas the six Mestre and
four Kihara additions have only bounded pilot exposure so far.
## Subsequent review through 3abeefeb

All three new independent arithmetic replays pass; their
[source bindings and transcripts](../../artifacts/generated-results/elliptic-curves/parallel_parent_review_through_3abeefeb_v1.json)
record20.969733168 supervised seconds in total.

- `3c3a5e73` proves [full private ramification](../rank-jump/ADDITIVE_COLLISION_BLOCKS_HAVE_FULL_PRIVATE_RAMIFICATION.md)
  for every additive three-radical block on the frozen sixteen fibres.
  The10520 independent private odd-prime witnesses exclude the entire
  respective block from Selmer, including combinations within the block.
  Adding the generic subgroup gives no extra Selmer dimension. This is a
  construction exclusion; the fibre ranks remain unaffected.
- `4e4ef90a` identifies the [ten-section node factor](../rank-jump/ADDITIVE_NORM_CARRIERS_AND_THE_TEN_SECTION_NODE.md).
  On MW16 at t=-2,120 triple norms share(t+2)^4. Removing that square factor
  reduces some necessary carrier genera from11 to9. The remaining3280
  branch polynomials are squarefree and pairwise coprime within their
  families. The carriers still impose necessary incidence conditions,
  not global point solubility. The generic node factor is not a new
  specialized rank selector.
- `3abeefeb` closes [the remaining global class](../rank-jump/THE_LAST_GLOBAL_CLASS_IS_ABSENT.md)
  on the fixed published R17 root curve: its rational Jacobian2-torsion
  is exactly the seventeen generic classes. The residual quartic has
  S4 Galois group and no rational point. On any k of the37 retained native
  supports with their original scalars, the multiquadratic pullback has
  exact generic rank17+k; all mixed-support character ranks vanish. This
  removes the proposed hidden common class from this model, while leaving
  exceptional rational specializations open. The new replay checks the
  height, finite characters and residual trace arithmetic; the previously
  retained independent surface count and the companion note's Artin–Tate
  and local descent arguments remain explicit dependencies.

None supplies a prospective condition constructing new specialized strict
classes. The production scores and completed frozen comparisons remain as
recorded. More generic collisions on these fixed inputs cannot substitute
for a new parent or a certified new specialized class.

The parent-side consequences are now concrete. The
[full554-pencil Mestre audit](MESTRE_DEGREE_THREE_PENCILS_2026-09-07.md)
finds no generic rank improvement in that retained portfolio. The
[Kihara section and coverage proof](KIHARA_QUADRATIC_SECTION_AND_PARENT_COVERAGE_2026-09-07.md)
closes the first parent's full geometric lattice and exhibits a quadratic
section over the unrestricted parent-ratio line. It also proves that the
existing rank14 path covers only negative ratios between-1/2 and0.
Broader intake should test the unrestricted parent coordinate, check
surface equivalence and rational section independence, and compare equal
completed exposure. Neither a fresh coordinate nor another fibration label
alone counts as a new parent.
## Subsequent review through 9f59caad

The independent arithmetic checks for `79a3487a` and `9f59caad` pass;
the [retained review](../../artifacts/generated-results/elliptic-curves/parallel_review_through_9f59caad_v1.json)
records6.364229609 supervised seconds. These results require a distinction
between S-units alone and S-units combined with generic classes.

The [bad-prime principalization experiment](../rank-jump/STRICT_BLOCKS_NEED_IDEALS_OUTSIDE_BAD_SUPPORT.md)
finds no generator among60 fixed squared bad-prime ideals on the fresh
matched pair, while recovering all six positive controls. The replay
checks3234 norms and exact ideal transports. Those bounded misses prove
neither nonprincipality nor a small class group. A separate retrospective
Artin-kernel bound limits the S-unit portion of three earlier known strict
blocks to dimensions1,0,0. That bound does not apply to the larger
generic-plus-S-unit constructor.

The [enlarged constructor](../rank-jump/GENERIC_IDEALS_AND_SUNIT_CORRECTIONS_CARRY_INCIDENCE.md)
has a positive incidence criterion. Write G for the generic Kummer group,
E_S for norm-square S-unit classes, U for the strict subgroup, and Phi_S
for the localized half-ideal map. Then

```
dim((G+E_S) intersect U) >= dim Phi_S(G).
```

The complete generic half-ideal images detect elementary S-class factors
of dimensions10,8,6 on the earlier controls. Relative to generic strict
dimensions1,2,6, this forces at least9,6,0 additional strict directions
in the enlarged constructor. The arithmetic replay verifies the Artin
evaluations and explicit dual generic-ideal words. These independence
witnesses still use old exceptional-point-derived characters; they do not
provide a prospective selector or explicit new rational elliptic points.

This corrects a possible overreading of the preceding exclusions: new
strict Kummer classes need not have half ideals outside the span of the
generic half ideals. A usable next arithmetic input is a certified
sufficient S-unit span and its joint localization matrix, together with
point-independent independence witnesses for the generic ideal image.
Explicit corrections and their rationality versus Sha remain unresolved.
No production score or completed experiment changes follow from these
retrospective matrices.

Separately, the [positive-parent fibre pilot](KIHARA_POSITIVE_PARENT_POINT_PILOT_2026-09-07.md)
now completes294 fixed boxes, with three certified directions on two
curves and strongest lower bound14. It also closes two specialization
losses and uses odd-prime certification for the new rational2-torsion
fibres. Those are concrete visibility and admission improvements in the
broader parent portfolio; they do not test the new S-unit criterion.
