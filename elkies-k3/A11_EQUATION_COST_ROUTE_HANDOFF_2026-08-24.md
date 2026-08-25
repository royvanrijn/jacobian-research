# A11 equation-cost route search handoff (2026-08-24)

## Status

**Promoted physical q4/orbit208 replacement (2026-08-25).** The physical
component-chamber search has superseded the q10/RR15 lifting target by the
fully certified splice

```text
component-9-zero 2A5/MW7 --q4/orbit208, degree 2--> canonical current 3A3/MW8.
```

Its primitive nef fibre is
`[2,2,1,1,1,1,1,1,1,1,1,1,0,0,-1,0,0,0,1]`.  The complete
all-section and finite-horizontal-wall gates pass.  More importantly for the
compiler, its special member is the literal equation-explicit I4 cycle

```text
F_q4 = old_zero + P1229 + old_A11_component_10 + old_A11_component_8.
```

Thus `P.O=0`, and P1229 is already exact over QQ with polynomial degrees
`(4,6,0)`.  The characteristic-zero resolved RR computation now passes
exactly with dimensions `4 -> 2 -> 2` (ambient, local-condition rank, kernel),
and maximum kernel rational height 420,710 bits.  Exact removal of the two old
I6 squared factors gives a degree-four quartic.  Its globally minimal Jacobian
has degrees `(8,12,24)` and fibres `3I4 + 12I1`, smooth infinity and Euler
number 24, hence `3A3/MW8` conditional on rho 19.  This took 355.70 seconds
and 782,296 KB maximum RSS, without Groebner basis or full discriminant
factorization.  Eleven known curves
have degree zero and four have degree one.  Reframing all four physical
degree-one curves selects effective C5 by minimum canonical-3A3 transport
growth.  The exact marked U, full 3A3 roots, both determinant-one NS
transports, and the pinned-R17 Gram identification pass.  The operational
score is `-1412` (`1388` gross positive burden before explicit-curve credits),
down by 5,883 from q10's 4,471.  The full q sequence from A11 is
`(8,4,4,4,4,4,4,4,6)`.
An exact direct-landing closure over all 181 `3A3/MW8` candidates at
q=4,6,8,10 leaves 56 fully nef candidates and confirms q4/orbit208 as the
unique equation-cost leader.

The lifting agent should switch from q10 to q4/orbit208.  Canonical artifacts
are
[`../artifacts/generated-results/elkies-k3-h3-a5a5-physical-q4o208-to-pinned-r17-certificate.json`](../artifacts/generated-results/elkies-k3-h3-a5a5-physical-q4o208-to-pinned-r17-certificate.json)
and
[`../artifacts/generated-results/elkies-k3-h3-a5a5-physical-q4o208-promoted-route-certificate.json`](../artifacts/generated-results/elkies-k3-h3-a5a5-physical-q4o208-promoted-route-certificate.json).
The exact equation artifact is
[`../artifacts/local/elkies-k3/q24-2a5-physical-q4o208-rr-qq.json`](../artifacts/local/elkies-k3/q24-2a5-physical-q4o208-rr-qq.json),
with status `PASS_EXACT_QQ_PHYSICAL_Q4O208_3A3_RR_AND_JACOBIAN`.
Only the equation-effective C5 point/sign and full old-curve equation marking
remain to be attached.
The q10 certificates remain exact and reusable, but q10 is no longer the
active lifting target.

**Equation-effective-zero correction (2026-08-25 10:38 CEST).** The former
4,199-point q4/orbit230--q6/orbit1315 suffix is still an exact marked-lattice
path, but it is **not an equation-realizable cost promotion as scored**. The
zero in the stored q4-return chamber is the class

```text
[27,1,-2,3,3,-8,-9,1,-8,-8,-6,-4,-2,0,1,-1,-1,-1,-1]
```

in the parent `2A5` coordinates, hence has `P.O=26`. Its MW quotient is P230,
but it differs from the exact effective P230 section by vertical roots. The
q6/orbit1315 fibre is intrinsically an exact degree-two self-zero neighbour,
yet its intersections with the effective P230 section and the original zero
are respectively 54 and 58. Therefore the low-pole term used in the 4,199
score came from a Weyl/chamber pseudo-zero. The exact q4 equation return and
its effective changed zero (the nonidentity component of the forward I2 at
infinity) remain valid; the returned frontier is being reranked using only
equation-effective curves. Do not use q6/orbit1315 as a lifting target on the
strength of the old score.

The earlier q6/orbit1307 score 10,334 is also withdrawn as an equation target.
In the exact physical component-9-zero `2A5` Gram, its stored fibre has
intersection `-1` with the first physical I6 affine component.  Reflections in
that affine component and old components 0 and 3 produce an exact primitive,
isotropic, physical-nef fibre, but its degrees are `C3=C5=C9=1` and `C10=0`.
Thus component 10 is not a section of the corrected pencil, invalidating the
certified component-10-zero return and landing.  The q6 horizontal itself
survives, with an improved expected RR profile `9 -> 3 -> 2`; reranking from
physical zeros C3, C5, and C9 remains useful as a lateral search.  The former
q104/13,518 comparator is also withdrawn: it has negative degree on physical
components C0 and C6.  Exact reduction against both I6 cycles takes 61
reflections and produces the promoted q10 target described below.
The durable replay is
[`../artifacts/local/elkies-k3/q24-2a5-q6o1307-physical-nef-audit.json`](../artifacts/local/elkies-k3/q24-2a5-q6o1307-physical-nef-audit.json),
with status `PASS_EXACT_Q6O1307_PHYSICAL_WEYL_REPAIR_REJECT_C10_ZERO`.

The first-q8 and D13 lattice improvements remain separately certified, but
their combined cost totals that included 4,199 are withdrawn.  The physical
q4/orbit208 target above is now the safe lifting suffix from the exact orbit12
equation; q10 below is retained as the superseded physical repair.

## Superseded physical q10 target

The canonical historical q104 class is not nef in the equation chamber.  Its
complete physical Weyl reduction yields

```text
component-9-zero 2A5/MW7 --q10, degree 2--> canonical current 3A3/MW8
```

The reduced fibre has `P.O=5`, three connected vertical layers, expected RR
ambient 15, eight already-explicit degree-zero curves, and four degree-one
curves.  It passes exact component, all-section, and complete finite
horizontal-wall gates.  The current-3A3 transport has determinant `-1`; its
composition with the stored suffix identifies pinned R17 by a full
determinant-one basis.  The operational score is 4,471, or 5,071 under the
older convention that omitted the two current-I6 affine curves.  This saves
9,047 (66.9%) against the withdrawn 13,518 presentation.

The canonical current-3A3 chamber zero is not effective in the q10 equation
marking.  All four physical degree-one curves were reframed exactly; old
component C5 minimizes child-frame coefficient growth (`58`) and is the
selected equation-effective zero.  Its full current-3A3 and pinned-R17 basis
transports are recorded in
[`../artifacts/local/elkies-k3/q24-2a5-direct-physical-q10-effective-c5-zero-certificate.json`](../artifacts/local/elkies-k3/q24-2a5-direct-physical-q10-effective-c5-zero-certificate.json).

The lifting agent should no longer switch to the q10 RR pencil; use physical
q4/orbit208 above.  The exact q10 artifacts retained for provenance are
[`../artifacts/local/elkies-k3/q24-2a5-direct-physical-q10-certificate.json`](../artifacts/local/elkies-k3/q24-2a5-direct-physical-q10-certificate.json)
and
[`../artifacts/generated-results/elkies-k3-h3-a5a5-direct-physical-q10-promoted-route-certificate.json`](../artifacts/generated-results/elkies-k3-h3-a5a5-direct-physical-q10-promoted-route-certificate.json).

The exact replayable obstruction is
[`../artifacts/generated-results/elkies-k3-h3-a5a5-q4o230-effective-return-zero-audit.json`](../artifacts/generated-results/elkies-k3-h3-a5a5-q4o230-effective-return-zero-audit.json).

The compact machine-readable handoff is
[`../artifacts/generated-results/elkies-k3-h3-a11-route-optimization-handoff.json`](../artifacts/generated-results/elkies-k3-h3-a11-route-optimization-handoff.json).
It records input hashes, exact marked degrees, certified transitions, rejected
branches, the exact endpoint identification, and the promotion decision.

## Promoted first-q8 improvement (2026-08-25)

After the initial H3 q6 reaches the equation-explicit `E8+E6/MW3` model, use

```text
E8+E6/MW3
 --q4 orbit11, zero=old_E8E6_component_1--> A2+D5+E7/MW3
 --q4--> E8+E6/MW3 (changed zero)
 --q4--> equation D13/MW4.
```

All three old-fibre degrees are two. Every fibre is primitive, nef, and
isotropic in its exact marked chamber; the component, affine-component,
all-section, and finite horizontal-wall gates pass. Every full 19-dimensional
NS transport and inverse is integral unimodular. The D13 landing is identified
with the current equation-D13 marking by a determinant-minus-one full basis,
and its canonical pinned-R17 basis has determinant one and Gram exactly `U`
plus the negative stored rank-17 frame.

The inherited-explicit operational cost falls from 5,802 for the direct q8 to
3,961 for q4,q4,q4, saving 1,841 or 31.730438%. With the unchanged promoted
D13 and A11 continuation, the recorded combined score falls from 35,324 to
33,483. A second beam layer from the returned E8+E6 marking certifies 38 exact
q4/q6/q8 changed-zero presentations; none beats the direct 1,952-point D13
exit (the best costs 4,291).

The widened first-edge boundary also has no winner: q10 degree two gives 16
exact presentations with best cost 10,699; q6/q9/q12 degree three gives 96
with best 7,515; and q8/q12/q16 degree four gives 167 with best 12,792. Thus
the promotion is not an artifact of optimizing only for small q or degree two.

Carrying the newly landed D13 zero directly into the next stage does not help:
the exact current-D12 target becomes a degree-two q60 fibre with direct score
42,833. The q4/q6/q8 crossover scan certifies 73 presentations; its best raw
credit score is 38,034 but its operational score is 44,934, above both that
direct comparator and the canonical 25,323-point D13 splice.

The full promoted H3 q sequence is

```text
6,4,4,4,4,4,24,6,8,4,4,6,4,4,4,4,4,4,4,4,6.
```

This remains an exact lattice detour, but the lifting agent should not switch
on the planning score alone until its q4 return identifies the effective
equation zero and reproduces the exit cost in the physical chamber. The
canonical lattice certificate is
[`../artifacts/generated-results/elkies-k3-h3-first-q8-q4o11-promoted-route-certificate.json`](../artifacts/generated-results/elkies-k3-h3-first-q8-q4o11-promoted-route-certificate.json),
and the exact negative second-loop boundary is
[`../artifacts/generated-results/elkies-k3-h3-first-q8-q4o11-c1-second-zero-changing-d13-presentations.json`](../artifacts/generated-results/elkies-k3-h3-first-q8-q4o11-c1-second-zero-changing-d13-presentations.json).

## Withdrawn-cost double-zero A11 lattice path (2026-08-25)

The following records the exact lattice path and historical score. Its cost
promotion is superseded by the equation-effective-zero correction above.

After A11 q8/orbit12 reaches the equation-explicit `2A5/MW7` model, use

```text
2A5/MW7
 --q4 orbit230, zero=old_A11_component_10--> A1+A4+A5/MW7
 --q4--> 2A5/MW7 (changed zero)
 --q6 orbit1315, zero=old_A5A5_component_1--> 3A2+A3/MW8
 --q4--> 2A5/MW7 (second changed zero)
 --q4--> current 3A3/MW8.
```

All five old-fibre degrees are two. Every fibre passes exact component,
affine, all-section, and finite horizontal-wall nef gates. All full NS maps
and inverses are integral unimodular. The current-3A3 landing and pinned R17
endpoint are full-basis identifications with determinant one and exact Gram,
not ADE/MW matches.

The formerly reported inherited-explicit operational score, with every degree-two edge floored
at its unavoidable horizontal cost 500, is 4,199. This improves the fully
certified q3372/q2052 route at 4,504, the q1307/q1581 route at 8,545, and the
direct q104 presentation at 13,518. The unfloored credit score is also
recorded (1,958), but is not used to promote over q3372 because negative curve
credits can make individual edges spuriously negative.

Together with the D13 q4/orbit11 splice, the conservative two-bottleneck score
is 29,522 instead of 41,403, saving 11,881 (28.695988%). The full H3 q
sequence is

```text
6,8,4,4,24,6,8,4,4,6,4,4,4,4,4,4,4,4,6.
```

The lifting agent must **not** switch to q6/orbit1315 on this score. The
historical marked-lattice certificate is
[`../artifacts/generated-results/elkies-k3-h3-a5a5-q4o230-q6o1315-promoted-route-certificate.json`](../artifacts/generated-results/elkies-k3-h3-a5a5-q4o230-q6o1315-promoted-route-certificate.json).
Its 1,000-candidate second-loop search is
[`../artifacts/generated-results/elkies-k3-h3-a5a5-q4o230-c10-second-zero-changing-3a3-presentations.json`](../artifacts/generated-results/elkies-k3-h3-a5a5-q4o230-c10-second-zero-changing-3a3-presentations.json).

## Promoted D13-prefix improvement (2026-08-25)

The marked lattice path below remains exact. Its compiler-cost promotion is
provisional after the q4/orbit230 physical-chamber counterexample: the D13
q4-return zero must be identified as an effective equation curve before the
25,323 score can be used for a lifting switch.

The new route prefix is

```text
D13/MW4
 --q4 orbit11--> D5+D9/MW3
 --q4, zero=old_D13_component_5--> D13/MW4 (changed zero)
 --q24--> current D12/MW5
 --q6--> equation A11/MW6.
```

All three replacement fibres pass exact component, affine, all-section, and
finite horizontal-wall gates.  Their marked U splittings and full NS maps are
integral and unimodular in both directions.  The q24 landing is identified by
a determinant-minus-one full basis with the stored current D12 frame, and the
composite endpoint basis has determinant one and Gram exactly `U` plus the
negative pinned `rank17_gram.txt`; this is not an ADE/MW-only match.

The new splice scores 25,323.  The direct q24 comparator scores 28,485 before
calibration and 27,885 after replacing its estimated RR ambient dimension 61
by the measured dimension 56.  Thus the certified planning reduction is
2,562, or 9.187735%. Combining this with the promoted double-zero A11 splice
reduces the conservative two-bottleneck score from 41,403 to 29,522, a saving
of 11,881. The full q sequence from H3 is

```text
6,8,4,4,24,6,8,4,4,6,4,4,4,4,4,4,4,4,6.
```

The standalone certificate is
[`../artifacts/generated-results/elkies-k3-h3-d13-q4o11-promoted-route-certificate.json`](../artifacts/generated-results/elkies-k3-h3-d13-q4o11-promoted-route-certificate.json).
Its compiler scores remain estimates until the replacement equations are
actually lifted.

## Critical zero-frame audit

The newly documented quintic AJ bridge is not composable as written.  Its
selected profile uses the D12 `A0` zero, while the orbit64 transition uses the
distinct `R3` zero.  In the selected equation-A11 Gram the stored `close_P24`
class has square `-3210`, not `-2`; the replayed word is
`(-23,73,-40,34,-4,-1)`, not `(1,0,0,0,0,1)`.

The exact rejection and the corrected R3-zero transports are in
[`../artifacts/generated-results/elkies-k3-h3-a11-quintic-bridge-zero-mismatch.json`](../artifacts/generated-results/elkies-k3-h3-a11-quintic-bridge-zero-mismatch.json).
Do not use the mismatched bridge as a certified equation-side carrier.

## A11 decision

The exhaustive 2,333-orbit q8 ranking uses exact marked intersections with
pinned R17 and the q25/MW reverse hubs together with compiler-facing cost
data.  Equation-side orbit12 is the unique marked-distance winner:

| orbit | pinned-R17 fibre degree | cost | `P.O` | RR estimate |
|---:|---:|---:|---:|---:|
| 12 | 29,900,919 | 17,869 | 6 | 14 |
| 2162 | 17,326,326,081 | 17,861 | 6 | 14 |

Orbit2162 has now also been tested as an EC branch rather than rejected only
by marked distance.  Its exact first-edge score is eight points smaller, but
its direct degree against orbit12 is 23.  The complete q4/q6 layer has 4,829
nef survivors and best degree 48.  All four root-adapted states tied at degree
48 have no q4/q6 continuation below 48, and all 26,477 q8 candidates have no
continuation below the direct degree 23.  Thus orbit2162 does not improve the
route in the certified local shells.  See
[`../artifacts/generated-results/elkies-k3-h3-a11-marked-target-neighbor-ranking.json`](../artifacts/generated-results/elkies-k3-h3-a11-marked-target-neighbor-ranking.json).

The first q8 was also tested by a dedicated zero-changing search.  All 1,119
q8 candidates passing the declared curve-nef gate were scanned; 897 pass the
independent exact nef/C2 gate and 1,482 complete return-and-exit presentations
were certified.  Raw abstract-component scoring produces two apparent wins,
but both disappear when degree-zero/one credits are restricted to curves
inherited explicitly from the equation A11 model and the first-edge target
coset penalty is restored.  The best retained presentation is
q8/orbit2013, q4 return, q10 exit at operational score 8,115 (raw 7,784),
versus 7,024 for direct
orbit12, so q8/orbit12 remains the lifting target.

The 24 re-zeroings with nonprimitive full root span are now handled by a
saturated unimodular frame retaining the actual simple-root lattice as an
embedded sublattice. The complete audit has 1,506 exact return/exit
presentations and no operational-cost winner. No fully certified q8
replacement was found. The machine audit is
[`../artifacts/generated-results/elkies-k3-h3-a11-zero-changing-q8-presentations.json`](../artifacts/generated-results/elkies-k3-h3-a11-zero-changing-q8-presentations.json).

Historically, the abstract q4/orbit230 returned chamber was additionally
widened beyond its complete q4/q6 frontier.  The q8 shell has 26,450 primitive neighbours and
12,775 nef survivors; its exact minimum inherited first-edge operational
score is 1,182.  The q10 shell has 88,755 primitive neighbours and 28,352 nef
survivors; its minimum is 3,459.  Adding the unavoidable 500-point return and
exit floors gives loop floors 2,182 and 4,459.  Both exceed every residual
strict bound 1,852--1,911 arising from the four abstract returned states that
could have beaten 4,199.  These are now retained only as marked-chamber search
closures: the effective-zero obstruction invalidates 4,199 as an equation
benchmark, so the bounds do not exclude an equation-effective replacement.

The combined hashed audit also records 31 expanded D13 two- and three-zero
beam states.  Their best combined raw score is 26,254, still above the
certified 25,323 D13 splice.  It also closes the compiler-relevant
current-D12 return topology around the direct 3,979-point D12-to-A11 edge.
The complete degree-two q4/q6/q8 and q10 shells have minima 5,167 and 9,446;
the degree-three q6/q9/q12 and degree-four q8/q12/q16 shells have minima
7,164 and 9,654.  All 12 twice-returned states with operational prefix 2,000
that could still beat 3,979 were expanded.  The best accumulated result is
6,666, and any deeper return has prefix at least 3,000 plus a loop-and-exit
floor of 1,500.  Thus deeper zero loops cannot beat the direct edge.  No new
lifting target is promoted.  Reproduce the compact handoff with

```text
python3 elkies-k3/scripts/build_h92_route_ec_extension_audit.py
```

which writes
[`../artifacts/generated-results/elkies-k3-h3-route-ec-extension-audit-2026-08-25.json`](../artifacts/generated-results/elkies-k3-h3-route-ec-extension-audit-2026-08-25.json).

The compiler-cheap orbit1991 branch is fully certified through its explicit
zero, but all 2,131 root-adapted q4/q6/q8/q10 continuations move farther from
pinned R17.  It is retained as a certified rejected detour, not a lifting
target.

A separate target-directed scan checked whether the earlier general-degree
search had simply stopped below the right shell scale.  In the exact
equation-A11 frame, the q8 degree-two probe reproduces orbit12 and its pinned
degree 29,900,919.  The adjacent q10 and q12 minima are 1,406,711,848 and
1,475,487,300; q14 has no nef hit, and q16 is 2,370,778,283.  For horizontal
degree three the real target slope predicts q near 18, beyond the earlier
q6/q9/q12 audit, but the new q15 minimum is already 1,446,758,396 and q18,
q21, and q24 have no nef candidates in the fixed A11 chamber.  The q24--q40
degree-four scan likewise has no nef candidate.  These are bounded
target-directed searches, not exhaustive non-existence theorems, but they
close the most plausible larger-q scale gap around orbit12.  See
[`../artifacts/generated-results/elkies-k3-h3-equation-a11-pinned-q8q16-targeted-shell-cvp.json`](../artifacts/generated-results/elkies-k3-h3-equation-a11-pinned-q8q16-targeted-shell-cvp.json),
[`../artifacts/generated-results/elkies-k3-h3-equation-a11-pinned-degree3-q15q24-targeted-shell-cvp.json`](../artifacts/generated-results/elkies-k3-h3-equation-a11-pinned-degree3-q15q24-targeted-shell-cvp.json), and
[`../artifacts/generated-results/elkies-k3-h3-equation-a11-pinned-degree4-q24q40-targeted-shell-cvp.json`](../artifacts/generated-results/elkies-k3-h3-equation-a11-pinned-degree4-q24q40-targeted-shell-cvp.json).

## Explicit-zero orbit12 branches

The abstract orbit12 child zero badly misprices the next edge: the historical
q4 fibre has q110 in that zero, and the historical zero is not one of the
known explicit curves.  Reframing with the already-explicit old A11 component
9 gives an `A5+A5/MW7` frame with maximum coefficient 58 and exact unimodular
transport in both directions.

In that selected explicit zero, the stored historical `2A5 -> 3A3` class was
formerly represented as q104, with `P.O=5` and estimated RR ambient 30.  That
representative is physically non-nef and is no longer a compiler-cost
baseline.  Its exact 61-reflection movable reduction is the promoted q10/RR15
presentation above.  The superseded comparison and its correction are in
[`../artifacts/generated-results/elkies-k3-h3-a5a5-current-route-equation-cost-audit.json`](../artifacts/generated-results/elkies-k3-h3-a5a5-current-route-equation-cost-audit.json).

Two compiler-facing exits pass the complete lattice gate, including the
all-section closest-vector nef test, marked `U`, full root data, and integral
forward/inverse NS transports:

| edge from explicit-zero `A5+A5/MW7` | child | `P.O` | RR estimate | explicit degree 0 / 1 | pinned degree |
|---|---|---:|---:|---:|---:|
| q4 orbit32 | `A1+D5+A5/MW6` | 0 | 5 | 7 / 2 | 1,562,674,900 |
| q6 orbit3372 | `2A1+A3+A5/MW7` | 3 | 11 | 6 / 4 | 1,461,309,907 |

The certificates are
[`../artifacts/generated-results/elkies-k3-h3-a5a5-explicit-zero-q4-orbit32-lattice-certificate.json`](../artifacts/generated-results/elkies-k3-h3-a5a5-explicit-zero-q4-orbit32-lattice-certificate.json)
and
[`../artifacts/generated-results/elkies-k3-h3-a5a5-explicit-zero-q6-orbit3372-lattice-certificate.json`](../artifacts/generated-results/elkies-k3-h3-a5a5-explicit-zero-q6-orbit3372-lattice-certificate.json).

Orbit32 still has no certified continuation to pinned R17.  Orbit3372 now has
the exact zero-changing continuation described below, but it is not yet a
strict compiler-cost win.  Neither branch is promoted.

The larger-q possibility was checked independently of the low-q search.
In the explicit-zero `A5+A5` frame, exhaustive dominant-Weyl enumeration
gave 88,755 q10 orbits and 245,617 q12 orbits.  After exact gates against all
18 identity-shell sections, the physical A11 components/affine component,
`A0`, `P24`, and both parent affine components, 1,109 q10 and 1,742 q12
candidates remain.  Their best compiler-cost profiles are nevertheless
dominated by the certified q4/q6 branches:

| shell | best cost | `P.O` | RR estimate | explicit degree 0 / 1 | pinned degree |
|---|---:|---:|---:|---:|---:|
| q10 | 6,650 | 6 | 17 | 6 / 2 | 1,469,032,088 |
| q12 | 8,931 | 8 | 21 | 6 / 2 | 1,485,985,660 |

The exact gates and scores are
[`../artifacts/generated-results/elkies-k3-h3-a5a5-explicit-zero-q10-gate.json`](../artifacts/generated-results/elkies-k3-h3-a5a5-explicit-zero-q10-gate.json),
[`../artifacts/generated-results/elkies-k3-h3-a5a5-explicit-zero-q10-equation-cost.json`](../artifacts/generated-results/elkies-k3-h3-a5a5-explicit-zero-q10-equation-cost.json),
[`../artifacts/generated-results/elkies-k3-h3-a5a5-explicit-zero-q12-gate.json`](../artifacts/generated-results/elkies-k3-h3-a5a5-explicit-zero-q12-gate.json), and
[`../artifacts/generated-results/elkies-k3-h3-a5a5-explicit-zero-q12-equation-cost.json`](../artifacts/generated-results/elkies-k3-h3-a5a5-explicit-zero-q12-equation-cost.json).
Because they are already cost-dominated, these survivors were not promoted
to child-root/transport certification.

An exact marked crossover audit from both certified branches to every stage
of the current suffix also found no splice.  The common `2A5/MW7` parent is
the only nef hit.  The first genuine target, `3A3/MW8`, already has degree
102 and q10,404 from q4 orbit32 (degree 96 and q14,016 from q6 orbit3372),
and has minimum section intersection -9 in both candidate chambers.  Later
targets grow rapidly and fail still more strongly.  See
[`../artifacts/generated-results/elkies-k3-h3-candidate-current-suffix-crossovers.json`](../artifacts/generated-results/elkies-k3-h3-candidate-current-suffix-crossovers.json).

### Superseded q6/orbit1307 zero-changing splice

An exhaustive equation-cost search now reframes every full-gate-passing q4/q6
child by every already-explicit old-A11 component of degree one.  It then
tests the inverse return to the exact marked current `2A5` fibre and the exit
to the exact marked current `3A3` fibre.  All retained returns and exits pass
the finite horizontal-wall test of Proposition C2 in addition to component,
affine, and all-section gates.

The search covers all 283 full declared-nef first edges and 558 available
explicit zeros.  It finds 109 exact zero loops and four strict cost winners.
The unique leader is

```text
2A5(explicit orbit12 zero)
 --q6 orbit1307--> A1+A3+A5/MW8
 --zero old_A11_component_10--
 --q4--> 2A5(new zero) --q6--> current 3A3/MW8.
```

All three horizontal degrees are two.  Its deterministic equation-cost score
is 10,334, versus 13,518 for the direct q104 presentation: an improvement of
3,184, or 23.553780%.  The direct comparator also has two named explicit
curves of negative degree, while every new edge passes the exact nef gate.

The landing is certified by a full determinant-minus-one basis change to the
stored current `3A3` stage.  Composing that marking with the unchanged suffix
gives a determinant-minus-one full-basis identification with canonical pinned
R17.  Thus the promoted route from A11 has q sequence

```text
8, 6, 4, 6, 4, 4, 4, 4, 4, 4, 6
```

and resumes the existing route immediately after `3A3`. This remains an exact
historical certificate, but the lifting agent should now use the cheaper
q4/orbit230 and q6/orbit1315 double-zero splice above after completing the
already-active A11 q8/orbit12 lift.
The canonical machine certificate is
[`../artifacts/generated-results/elkies-k3-h3-a5a5-q6o1307-promoted-route-certificate.json`](../artifacts/generated-results/elkies-k3-h3-a5a5-q6o1307-promoted-route-certificate.json),
and the exhaustive ranking is
[`../artifacts/generated-results/elkies-k3-h3-a5a5-zero-changing-loop-search.json`](../artifacts/generated-results/elkies-k3-h3-a5a5-zero-changing-loop-search.json).

### Exact zero-changing loop through q6 orbit3372

The negative crossover conclusion above applies to a direct edge in the
selected orbit3372 zero.  A target-directed second layer found a different
phenomenon: the unique q6 candidate nearest both pinned R17 and `3A3` is the
original `2A5` fibre itself, now equipped with a new zero.  Exact
certification gives the loop and exit

```text
2A5(explicit orbit12 zero) --q6--> 2A1+A3+A5/MW7
 --q6--> 2A5(new zero) --q12--> 3A3/MW8.
```

All three old-fibre degrees are two.  The first edge is the previously
certified orbit3372 exit.  On the return edge, nine source components have
degree zero and two have degree one; on the q12 exit, eight have degree zero
and four have degree one.  Both have exact all-section minimum zero.  The
return is identified with the exact marked current `2A5` fibre, not merely by
its ADE type.

Following the corresponding marked fibres produces the fully certified route

```text
A11 --q8--> 2A5 --q6--> 2A1+A3+A5 --q6--> 2A5
 --q12--> 3A3 --q18--> A3+2A2 --q12--> 5A1
 --q16--> 4A1 --q30--> 3A1 --q30--> 2A1
 --q20--> A1 --q30--> pinned rootless MW17.
```

Every edge is primitive, nef, isotropic, and has exact component/all-section
gates, root data, marked `U`, and determinant-one transports in both
directions.  The final full basis has determinant one and identifies its Gram
exactly with `U` plus the negative pinned `rank17_gram.txt`; this is therefore
a complete pinned route, not an ADE/MW match.

It is not promoted yet.  Under the current deterministic compiler weights,
the three-edge zero-changing replacement scores 14,641, compared with 13,518
for the formal direct q104 presentation.  The loop has a real gate advantage:
all its edges pass exact nefness, whereas the q104 presentation has two named
explicit curves of negative degree and no named degree-one curve.  That gate
advantage is not enough to claim a strict measured cost improvement without
an equation replay.  The complete certificate and no-promotion flag are in
[`../artifacts/generated-results/elkies-k3-h3-a5a5-explicit-zero-detour-route-certificate.json`](../artifacts/generated-results/elkies-k3-h3-a5a5-explicit-zero-detour-route-certificate.json).

Correct-scale bounded probes from orbit3372 found no degree-three or
degree-four continuation through q32, and its degree-two q4--q24 optimum is
exactly the q6 return above.  The analogous orbit32 probe found no pinned-nef
candidate in those shells.  Thus the loop is a reusable exact zero-translation
mechanism, but not yet the requested cheaper lifting target.

## New certified reverse hub

The semistable `A5+A4+2A3/MW2` frame is now a valid pinned reverse hub.  Its
lost bounded-beam q14 witness was recovered by an exact dominant-Weyl search:
22,971 orbits were tested, only one all-section-nef orbit had root data
`(15,74,480)` and an integral isometry to the stored frame.  The recovered
reverse q14 fibre has factorization `7*2`, minimum section intersection zero,
and a unique exact D6-hub match.

The remaining inverse q4,q4,q25 transports initially lay outside their source
component chambers.  Exact affine-Weyl reduction followed by rebuilding and
integral endpoint identification gives the fully nef suffix

```text
A5+A4+2A3/MW2 --q14--> D6+A5+A3/MW3
 --q4--> D4+A3+2A2+2A1/MW4
 --q4--> A3+7A1/MW7 --q25--> pinned rootless/MW17.
```

The old-fibre degrees are `2,2,2,5`; the all-section minima are `0,1,0,1`.
Every edge has an exact marked U and full determinant-one transports in both
directions.  The certificates are
[`../artifacts/generated-results/elkies-k3-h3-semistable-mw2-pinned-transport.json`](../artifacts/generated-results/elkies-k3-h3-semistable-mw2-pinned-transport.json)
and
[`../artifacts/generated-results/elkies-k3-h3-semistable-mw2-reverse-suffix-nef.json`](../artifacts/generated-results/elkies-k3-h3-semistable-mw2-reverse-suffix-nef.json).

This is a substantial bidirectional-search improvement, but not yet a new
lifting route.  Orbit12 remains the best of all 2,333 A11 q8 neighbours even
when ranked against this hub.  Its direct hub degree is 4,244,273,479,663;
the best of the 306 compiler-gated q4/q6 exits still has degree
205,669,321,746,575.  Thus no cheap certified A11-to-hub meeting path is yet
known.  See
[`../artifacts/generated-results/elkies-k3-h3-a5a5-semistable-hub-survivor-ranking.json`](../artifacts/generated-results/elkies-k3-h3-a5a5-semistable-hub-survivor-ranking.json).

## Reverse-hub beam and general-degree audit

The reverse search has now been extended beyond the stored suffix rather than
stopping at matching ADE/MW labels. Every state was root-adapted, marked in
equation-A11 coordinates, and filtered by the exact component and all-section
closest-vector nef gates.

From the semistable MW2 hub, the complete degree-two q4/q6 and q8/q10/q12
frontiers contain 305 and 12,789 candidates, with 150 and 485 nef survivors.
The only one-step improvement is q12 orbit7798,
`A1+A2+A5+D5/MW4`, which lowers the orbit12 degree only from
4,244,273,479,663 to 4,133,250,357,677 while costing 9,270 (`P.O=9`, RR 26).
Its complete q4/q6 second layer moves back out to degree 6,809,438,300,451.
The next two one-step Pareto states, q10 orbits 2380 and 3063, also move out in
their complete second layers. Degree one has no nef representative through
q12; degree three has 3,752 nef survivors but closest degree
10,545,541,040,102; the tested degree-four and degree-five shells have none.

The nearer certified q25 reverse hubs were also searched. Their low-q
frontiers behave as follows:

| certified source hub | q4/q6 candidates | nef | source orbit12 degree | closest new degree |
|---|---:|---:|---:|---:|
| `A3+7A1/MW7` | 49,494 | 35,818 | 397,811,472 | 2,222,985,364 |
| `D4+A3+2A2+2A1/MW4` | 3,633 | 2,040 | 3,277,287,981 | 26,783,292,667 |
| `A5+D4+2A2+A1/MW3` | 993 | 520 | 62,185,731,764 | 933,617,155,681 |

Following the best MW7 candidate through a complete q4 second layer returns
only to the original marked MW7 hub at degree 397,811,472; its best genuine
non-backtrack has degree 2,385,065,530. Thus the q25,4,4,4 reverse corridor
has no low-q meeting with orbit12.

The principal artifacts are
[`../artifacts/generated-results/elkies-k3-h3-semistable-mw2-q8q10q12-equation-cost.json`](../artifacts/generated-results/elkies-k3-h3-semistable-mw2-q8q10q12-equation-cost.json),
[`../artifacts/generated-results/elkies-k3-h3-semistable-mw2-q12o7798-q4q6-equation-cost.json`](../artifacts/generated-results/elkies-k3-h3-semistable-mw2-q12o7798-q4q6-equation-cost.json),
[`../artifacts/generated-results/elkies-k3-h3-q25-mw7-q4q6-marked-frontier-compact.json`](../artifacts/generated-results/elkies-k3-h3-q25-mw7-q4q6-marked-frontier-compact.json), and
[`../artifacts/generated-results/elkies-k3-h3-q25-mw7-q6o36810-q4-marked-frontier-compact.json`](../artifacts/generated-results/elkies-k3-h3-q25-mw7-q6o36810-q4-marked-frontier-compact.json).

## Non-degree-two A11 exits

The equation-side A11 degree-three q6/q9/q12 shells contain 18,152 exact
orbits. After the known equation-curve, old-affine, and all-section gates,
6,431 survive. The apparent cost leader q6 orbit385 (`A11/MW6`, score -536)
is rejected by an undeclared section wall with intersection -1.

The real compiler-cost leader is q9 orbit1802,
`A3+A3+D6/MW5`, with score 914, `P.O=1`, RR estimate 8, exact-section tier
zero, nine explicit degree-zero components, and three degree-one components.
It is fully certified, but its pinned degree is 11,068,422,340 and its
orbit12 degree is 15. Five cost/distance Pareto states (q9 orbits 1802, 2800,
2793, 1542, and 956) were given exact marked transports and complete q4/q6
second layers; every non-backtracking continuation moves farther from
orbit12. The closest degree-three pinned candidate is q12 orbit474 at
1,343,863,080, still far worse than orbit12's 29,900,919. Candidates meeting
orbit12 at degree two start at cost 17,640 and therefore add an extra edge to
only a 229-point first-edge saving over the direct orbit12 cost 17,869.

See
[`../artifacts/generated-results/elkies-k3-h3-a11-q6q9q12-degree3-full-nef-equation-cost.json`](../artifacts/generated-results/elkies-k3-h3-a11-q6q9q12-degree3-full-nef-equation-cost.json)
and
[`../artifacts/generated-results/elkies-k3-h3-a11-q9d3o1802-lattice-certificate.json`](../artifacts/generated-results/elkies-k3-h3-a11-q9d3o1802-lattice-certificate.json).

These negative A11 degree-three branches did not change the decision at that
stage.  The later q6/orbit1307 zero-loop search above supplies the strict-cost
fully certified route and supersedes that earlier no-promotion boundary.

## Earlier equation-D13 branch

The search was also moved back to the current physical equation-D13 marking,
without reopening the certified q24/orbit85 equation lift.  Exact marked
degree-two, degree-three, and degree-four presentations through q20 produced
9,559 candidates and 753 full-component/full-section-nef survivors.  The
closest candidate in each shell is still much farther from pinned R17 than
the current q24 D12 fibre, whose exact pinned degree is 19,775,425,087:

| equation-D13 shell | candidates | nef | closest pinned degree |
|---|---:|---:|---:|
| q4/q6/q8, degree 2 | 194 | 62 | 11,554,571,751,627 |
| q10--q20, degree 2 | 4,427 | 258 | 11,562,968,221,003 |
| q9--q18, degree 3 | 2,185 | 140 | 17,125,829,545,574 |
| q8--q20, degree 4 | 2,753 | 293 | 23,531,813,928,125 |

The best compiler-facing first hop is q6 orbit42 to `E6+D8/MW3`: it has
`P.O=1`, thirteen component degrees zero, one component degree one, minimum
section intersection zero, exact marked U, and determinant-one transports.
This is nevertheless a rejected detour.  Its complete q4/q6 second layer
moves still farther out, from pinned degree 11,554,571,751,627 to
952,659,391,354,350.

The compact exact artifacts are
[`../artifacts/generated-results/elkies-k3-h3-equation-d13-marking.json`](../artifacts/generated-results/elkies-k3-h3-equation-d13-marking.json),
[`../artifacts/generated-results/elkies-k3-h3-equation-d13-q4q6q8-marked-frontier-adapted-compact.json`](../artifacts/generated-results/elkies-k3-h3-equation-d13-q4q6q8-marked-frontier-adapted-compact.json),
[`../artifacts/generated-results/elkies-k3-h3-equation-d13-q10to20-degree2-marked-frontier-compact.json`](../artifacts/generated-results/elkies-k3-h3-equation-d13-q10to20-degree2-marked-frontier-compact.json),
[`../artifacts/generated-results/elkies-k3-h3-equation-d13-q9to18-degree3-marked-frontier-compact.json`](../artifacts/generated-results/elkies-k3-h3-equation-d13-q9to18-degree3-marked-frontier-compact.json), and
[`../artifacts/generated-results/elkies-k3-h3-equation-d13-q8to20-degree4-marked-frontier-compact.json`](../artifacts/generated-results/elkies-k3-h3-equation-d13-q8to20-degree4-marked-frontier-compact.json).

## Target-directed reverse search from pinned R17

A target-directed closest-vector probe now searches the correct shell scale
rather than treating small q as the primary objective.  The bounded method is
calibrated by reproducing the exact exhaustive q4 optimum.  At q8 it finds a
primitive-nef `A1/MW16` neighbour whose orbit12 marked degree is 4,595,988,
down from 29,900,919.  Exact comparison with the stored current-route
transport shows that this q8 fibre is not a new neighbour: it is precisely the
current route's penultimate `A1/MW16` fibre.  The fibre classes agree in pinned
coordinates and the complete A1 bases are related by a determinant-one
isometry.  The identification is recorded in
[`../artifacts/generated-results/elkies-k3-h3-pinned-r17-q8-current-route-a1-identification.json`](../artifacts/generated-results/elkies-k3-h3-pinned-r17-q8-current-route-a1-identification.json).

Repeating the same slope calculation from that existing suffix state gives
the exact certified lateral chain

```text
pinned R17 --q8--> A1/MW16 --q18--> 2A1/MW15
 --q10--> (root data 5,12,24)/MW12
 --q14--> 7A1/MW10.
```

The successive orbit12 marked degrees are `4,595,988`, `2,380,028`,
`1,543,114`, and `1,461,358`.  Every edge has a primitive-nef audit, exact
child roots, marked U, and determinant-one transports composed back to both
equation A11 and pinned R17.  The final state is nevertheless a local diamond:
the searched degree-two and degree-three continuations return to earlier
states.  This is a marked-distance improvement along a detour from the known
suffix, not a certified meeting with orbit12 and not a lifting target.  The
first certificate is
[`../artifacts/generated-results/elkies-k3-h3-pinned-r17-q8-orbit12-cvp-lattice-certificate.json`](../artifacts/generated-results/elkies-k3-h3-pinned-r17-q8-orbit12-cvp-lattice-certificate.json),
and the current beam endpoint is
[`../artifacts/generated-results/elkies-k3-h3-pinned-r17-q8q18q10q14-orbit12-cvp-lattice-certificate.json`](../artifacts/generated-results/elkies-k3-h3-pinned-r17-q8q18q10q14-orbit12-cvp-lattice-certificate.json).

Every current-route suffix basis was then transported into pinned coordinates,
with exact checks that `current_A11=equation_A11`,
`current_A5A5=orbit12`, and `current_rootless=pinned_R17`.  Target-directed
q4/q6/q8/q10 probes against the middle and late suffix again choose the
existing current-route A1 fibre as the q8 winner.  Its degrees against
`3A3`, `A3+2A2`, `5A1`, `4A1`, `3A1`, and `2A1` are respectively
365,373, 55,374, 2,761, 149, 8, and 2.  No tested state gives a new
two-edge degree-two crossover that skips a material part of the suffix.  The
full suffix marking is
[`../artifacts/generated-results/elkies-k3-h3-pinned-r17-current-suffix-marking.json`](../artifacts/generated-results/elkies-k3-h3-pinned-r17-current-suffix-marking.json).

Two exact q25/MW7 corridors were also isolated.  The better all-degree-two
profile is

```text
q25/MW7 --q8,d2,P.O2--> 6A1/MW11
 --q16,d2,P.O6--> pinned R17.
```

The alternative degree-three detour is

```text
q25/MW7 --q6,d2,P.O1--> 2A2+2A1/MW11
 --q39,d3,P.O10--> pinned R17.
```

The latter is identified with the canonical physical q25 basis by a full
unimodular change of basis, not merely by its `(10,26,512)` root data.  It
also improves the marked contact with the physical MW3 hub from 302 to 250.
Its complete certificate is
[`../artifacts/generated-results/elkies-k3-h3-q25mw7-pinned-r17-degree2-degree3-detour.json`](../artifacts/generated-results/elkies-k3-h3-q25mw7-pinned-r17-degree2-degree3-detour.json).
Neither corridor currently has a certified cheap A11/orbit12 prefix, and the
degree-three detour's `P.O=10` tradeoff is not clearly preferable to the
all-degree-two corridor.  The physical MW3 frame was separately probed at its
correct q24--q40 scale; the bounded fixed-chamber search found no nef
degree-two candidate.  These reverse-corridor results do not affect the later
q6/orbit1307 promotion.

## Requested exchange with the lifting agent

The lifting agent supplied an exact characteristic-zero bridge section during
this search. The artifact
[`../artifacts/local/elkies-k3/q24-a11-bridge-m-section-qq.json`](../artifacts/local/elkies-k3/q24-a11-bridge-m-section-qq.json)
has status `PASS_EXACT_Q24_A11_BRIDGE_M_SECTION_QQ`, equation `P.O=8`, pinned
`P.O=5`, degrees `(20,30,8)`, and maximum coefficient size 909,707 bits. This
closes the requested bridge-class gap but confirms that the current lift is
extremely expensive.

Please add a note or machine artifact with either of the following when it
becomes available:

- the actual zero chosen by the equation compiler after orbit12;
- measured resolved-RR dimensions or coefficient growth for the current
  `2A5 -> 3A3` edge and later suffix.

Those data can materially change the equation-cost ordering without changing
the lattice certificates above.

The requested equation data are now partly available.  The q4/orbit230
forward edge is exact over QQ with divisor `O+P230`, `P230.O=2`, section
degrees `(8,12,2)`, resolved-RR `ambient=6`, collision rank `4`, and `h0=2`.
Its quartic has bidegree `(4,4)` and maximum rational coefficient height
2,123,040 bits; the Jacobian has degrees `(8,12,22)`, height 5,931,683 bits,
and fibres `I2+I5+I6+11I1`.  The exact component-10 marking and q4 return
identify the physical effective changed zero as the P230 branch of the
infinity `I2`; no large Groebner calculation was used.  The artifacts are
[`../artifacts/local/elkies-k3/q24-2a5-to-a1a4a5-q4o230-resolved-rr-qq.json`](../artifacts/local/elkies-k3/q24-2a5-to-a1a4a5-q4o230-resolved-rr-qq.json),
[`../artifacts/local/elkies-k3/q24-2a5-to-a1a4a5-q4o230-equation-marking-qq.json`](../artifacts/local/elkies-k3/q24-2a5-to-a1a4a5-q4o230-equation-marking-qq.json),
and
[`../artifacts/local/elkies-k3/q24-a1a4a5-to-2a5-q4-return-resolved-rr-qq.json`](../artifacts/local/elkies-k3/q24-a1a4a5-to-2a5-q4-return-resolved-rr-qq.json).

After the effective-zero pivot, P1229 for the q6/orbit1307 route was certified
exactly over QQ with degrees `(4,6,0)`, 1,000-digit branch agreement, and
maximum rational coefficient height 1,259,550 bits.  P146 is being rebuilt as
the exact short word `P_affine+P1+P32` modulo the trivial lattice; the regular
P1 and P32 Hensel branches are active.  The later physical-nef audit preserves
this horizontal and improves the expected RR profile, but withdraws its
component-10 continuation and 10,334 score.  See
[`../artifacts/local/elkies-k3/q24-2a5-p1229-scaled-x-qq.json`](../artifacts/local/elkies-k3/q24-2a5-p1229-scaled-x-qq.json).

On the reverse side, the complete current-`3A3` degree-two q4/q6/q8 box has
169,725 primitive inputs and 94,403 exact-nef root-adapted candidates.  Its
best q9 contact has degree 4,091 versus the direct 735, and its best contact
with the forward compiler-friendly MW2 state has degree 46,707 versus direct
8,391.  Therefore this entire box cannot improve either meeting.  The compact
rankings are
[`../artifacts/generated-results/elkies-k3-h3-current_3A3-d2-q4q6q8-q9-frontier.json`](../artifacts/generated-results/elkies-k3-h3-current_3A3-d2-q4q6q8-q9-frontier.json)
and
[`../artifacts/generated-results/elkies-k3-h3-current_3A3-d2-q4q6q8-mw2-frontier.json`](../artifacts/generated-results/elkies-k3-h3-current_3A3-d2-q4q6q8-mw2-frontier.json).
The complete degree-three q6/q9/q12 dominant-orbit box contains respectively
25,844, 286,748, and 1,659,547 orbits.  An exact pre-adaptation marked-degree
filter retains none below the direct q9 contact degree 735.  Thus this entire
degree-three box also cannot improve the q9 meeting.  Its replay artifact is
[`../artifacts/generated-results/elkies-k3-h3-current_3A3-d3-q6q9q12-q9lt735-neighbors.json`](../artifacts/generated-results/elkies-k3-h3-current_3A3-d3-q6q9q12-q9lt735-neighbors.json).

## Reproduction

Use the documented Sage Python launcher:

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/rank_h92_a11_marked_target_neighbors.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_h92_a11_o12_explicit_zero_frames.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_h92_a5a5_explicit_zero_q4_orbit32.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/search_h92_a5a5_zero_changing_loops.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_h92_a5a5_q6o1307_promoted_route.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/search_h92_d13_zero_changing_d12_presentations.sage --mode a11

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/search_h92_d13_zero_changing_d12_presentations.sage --mode d13

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_h92_d13_q4o11_promoted_route.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/gate_h92_a5a5_explicit_zero_large_q.sage --q 10 \
  --output artifacts/generated-results/elkies-k3-h3-a5a5-explicit-zero-q10-gate.json

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/audit_h92_candidate_current_suffix_crossovers.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/export_h92_d13_equation_marking.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/search_root_adapted_weyl_neighbors.sage \
  --frame artifacts/generated-results/elkies-k3-h3-equation-d13-root-adapted-frame.txt \
  --root-rank 13 --q 4 --q 6 --q 8 --degree 2 --adapt-mw-at-least 0 \
  --output artifacts/generated-results/elkies-k3-h3-equation-d13-q4q6q8-degree2-all-adapted.json \
  --frames-dir artifacts/generated-results/elkies-k3-h3-equation-d13-q4q6q8-degree2-all-adapted-frames

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/rank_h92_marked_root_adapted_frontier.sage \
  --neighbors artifacts/generated-results/elkies-k3-h3-equation-d13-q4q6q8-degree2-all-adapted.json \
  --marking artifacts/generated-results/elkies-k3-h3-equation-d13-marking.json \
  --target pinned_R17 --retain 100 \
  --output artifacts/generated-results/elkies-k3-h3-equation-d13-q4q6q8-marked-frontier-adapted-compact.json

# Physical 2A5 component chamber and promoted q4/orbit208.
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/export_h92_a5a5_physical_component_chamber_marking.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/search_root_adapted_weyl_neighbors.sage \
  --frame artifacts/generated-results/elkies-k3-h3-a5a5-physical-component-chamber-frame.txt \
  --root-rank 10 --q 4 --q 6 --q 8 --q 10 --degree 2 \
  --filter-marking artifacts/generated-results/elkies-k3-h3-a5a5-physical-component-chamber-marking.json \
  --filter-target physical_q10_current_3A3 --filter-max-degree 20 \
  --adapt-mw-at-least 8 --rank-growth-only \
  --output artifacts/generated-results/elkies-k3-h3-a5a5-physical-component-chamber-d2-q4q6q8q10-physical3a3le20-neighbors.json

python3 elkies-k3/scripts/extract_h92_root_data_frontier.py \
  --input artifacts/generated-results/elkies-k3-h3-a5a5-physical-component-chamber-d2-q4q6q8q10-physical3a3le20-neighbors.json \
  --root-data 9,36,64 \
  --output artifacts/generated-results/elkies-k3-h3-a5a5-physical-component-chamber-d2-q4q6q8q10-3a3-only-neighbors.json

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/rank_h92_marked_root_adapted_frontier.sage \
  --neighbors artifacts/generated-results/elkies-k3-h3-a5a5-physical-component-chamber-d2-q4q6q8q10-3a3-only-neighbors.json \
  --marking artifacts/generated-results/elkies-k3-h3-a5a5-physical-component-chamber-marking.json \
  --target physical_q10_current_3A3 --retain 181 \
  --output artifacts/generated-results/elkies-k3-h3-a5a5-physical-component-chamber-d2-q4q6q8q10-3a3-only-frontier.json

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/score_h92_marked_frontier_equation_cost.sage \
  --neighbors artifacts/generated-results/elkies-k3-h3-a5a5-physical-component-chamber-d2-q4q6q8q10-3a3-only-neighbors.json \
  --marking artifacts/generated-results/elkies-k3-h3-a5a5-physical-component-chamber-marking.json \
  --nef-frontier artifacts/generated-results/elkies-k3-h3-a5a5-physical-component-chamber-d2-q4q6q8q10-3a3-only-frontier.json \
  --retain 181 \
  --output artifacts/generated-results/elkies-k3-h3-a5a5-physical-component-chamber-d2-q4q6q8q10-3a3-only-equation-cost.json

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_h92_marked_degree2_candidate.sage \
  --source-marking artifacts/generated-results/elkies-k3-h3-a5a5-physical-component-chamber-marking.json \
  --source-frame artifacts/generated-results/elkies-k3-h3-a5a5-physical-component-chamber-frame.txt \
  --fibre 2,2,1,1,1,1,1,1,1,1,1,1,0,0,-1,0,0,0,1 \
  --candidate-label physical-q4-orbit208-3A3 \
  --target physical_q10_current_3A3 \
  --frame-output artifacts/generated-results/elkies-k3-h3-a5a5-physical-q4o208-3a3-frame.txt \
  --output artifacts/generated-results/elkies-k3-h3-a5a5-physical-q4o208-3a3-lattice-certificate.json

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_h92_a5a5_physical_q4o208_to_pinned_r17.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_h92_a5a5_physical_q4o208_rr_qq.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_h92_a5a5_physical_q4o208_promoted_route.sage

.venv/bin/python \
  elkies-k3/scripts/build_h92_a11_route_optimization_handoff.py
```

The exhaustive shell dumps are intentionally not retained after their
survivors are copied into compact gate artifacts.  Regenerate a shell with
`search_root_adapted_weyl_neighbors.sage` before replaying a new branch.
