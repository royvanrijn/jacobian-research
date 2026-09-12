# Exceptional ancestry, equation-only arithmetic, and strict descent classes

The subsequent [completed-wide-population arithmetic profile](WIDE_ARITHMETIC_PROFILE_2026-09-12.md)
extends the equation-only lane to all 2,080 fibres under a separately authorized
fixed budget. Its BASE+LOCAL census is running; CLASS is off. The finite
three-lane results below are unchanged, not superseded by an unfinished census.

## Result and boundary

The first frozen three-lane comparison is complete. It does **not** yet explain
the difference between the known jumps 17→31 and 17→25.

| Exact finite comparison | 302, original X1092 fibre at 0 | 11952 at 921/653 |
|---|---:|---:|
| Chosen exceptional directions | 14 | 8 |
| Minimum multisection degree through each chosen point | 2 | 2 |
| Best genus constructed through every chosen point | 1 | 1 |
| Global minimum genus | UNKNOWN | UNKNOWN |
| Tested degree-two carriers | 504 | 288 |
| Distinct tested covers over the fixed base | 504 | 288 |
| Different targets sharing a tested cover | 0 | 0 |
| Equation-only Brumer–Kramer offset | 12 | 9 |
| Independent class-group 2-rank upper bound | UNKNOWN | UNKNOWN |

Genus-one ancestry is available on **both** fibres through their existing
alternate fibrations. It is not a distinguishing property of 302. One
alternate-fibration pencil containing many carrier curves is not one quadratic
cover carrying many independent sections.

The fixed-target [records](../../artifacts/generated-results/elliptic-curves/rank_triangle_v1/targets.json)
contain every point, exact model transport, best carrier, branch polynomial,
cover lift, rational map, and the 302 local-class labels. The
[combined replay](../../artifacts/generated-results/elliptic-curves/rank_triangle_v1/verified.json)
checks 792 carrier identities, 792 altered-carrier negative controls, eleven
complete local-arithmetic rows, and the strict dictionary. The
[manifest](../../artifacts/generated-results/elliptic-curves/rank_triangle_v1/manifest.json)
binds the inputs, scripts, outputs and retained failure logs.

No new point-search campaign was launched. The running broad search was not
changed. Neither fibre has a new exact-rank certificate.

## 1. What the ancestry experiment actually tests

### Fixed representatives and parent surfaces

Use the prior sealed accessibility input's integral bases `M1,...,M17` followed
by `E1,...,E14` or `E1,...,E8`. The 302 exceptional representatives are public
points 1,3,4,6,9,13,14,17,19,22,26,29,30,31, with one-based indexing. They are
not the different fourteen words in the historical seed-universality packet.
Exact Weierstrass isomorphisms are selected by matching all seventeen generic
sections up to sign, and every chosen point is checked after transport.

The 302 parent is the recovered MW17 fibration on X1092. The compact 11952
parent is on X948. Its existing direct old-R17↔11952 map supplies the matched
alternate-fibration control; the historical degree-11511 route is not used.

Let the short parent equation be `y²=x³+A(t)x+B(t)`, and the target over `t0`
be `(px,py)`. The common atlas comprises:

1. The vertical-x curve `v²=px³+A(t)px+B(t)`.
2. Constant-slope chords through each of the 34 signed generic sections.
   For a signed section `(tx,ty)`, put
   `m=(py+ty(t0))/(px-tx(t0))`. The other two intersections of
   `y=m(x-tx)-ty` with the cubic form a quadratic multisection with radical

   `D=m⁴-6 tx m²-8 ty m-3 tx²-4A`,

   and map `x=(m²-tx+v)/2`, `y=m(x-tx)-ty`, where `v²=D`.

3. One carrier in a known alternate-fibration pencil: class 1 on 302 and
   old R17 on 11952. Its parameter is obtained by evaluating the exact map at
   the target. Fixing that parameter gives the emitted genus-one bisection.

The first atlas alone has best genus 3 on both fibres. The matched alternate
fibrations lower this to 1 on both. Counts by normalization genus are:

| Genus | 302 | 11952 |
|---|---:|---:|
| 1 | 14 | 8 |
| 3 | 420 | 80 |
| 5 | 70 | 168 |
| 7 | 0 | 32 |

The exact [original atlas](../../artifacts/generated-results/elliptic-curves/rank_triangle_v1/geometry.json)
and [11952 matched supplement](../../artifacts/generated-results/elliptic-curves/rank_triangle_v1/native-carriers.json)
retain all maps, not merely winning costs. Every genus-one member has a
different parameter and a different fixed-base cover within its fibre.

### Degree, genus, branch divisor, and equivalence

Each radical is normalized exactly as `D=c q(t) s(t)²`, with `q` monic and
squarefree. The normalized cover is `v²=c q(t)`. Its branch divisor is the
reduced zero divisor of `q`, plus infinity exactly when `deg(q)` is odd.
Riemann–Hurwitz gives genus `floor((deg(q)-1)/2)`. The linear-in-v point map
and the rational lift over `t0` are recorded and replayed. At least one map
coefficient of v is nonzero, so the curve really retains degree two.

Two such quadratic extensions over the **fixed** t-line are equal exactly
when their monic branch polynomials agree and the ratio of their constants
is a rational square. This is the equivalence used here. We do not identify
covers under arbitrary PGL2 changes of the base, or claim a canonical
integer representative of every constant square class.

Degree one is excluded: a horizontal degree-one curve gives a rational
generic section. Since M17 spans the full rank-17 generic Mordell–Weil group
over Q, its specialization lies in the rational span of the specialized
M17. Each certified exceptional target lies outside that span. The displayed
degree-two constructions therefore prove **minimum degree two**, not minimum
genus one. Rational bisections through these targets remain UNKNOWN.

The genus-one curves have the known alternate-fibre family labels. They do
not receive labels from the rational-bisection orbit census. No new complete
Néron–Severi orbit classification is asserted.

There is also a useful trace restriction: if one degree-two multisection has
two distinct rational points over the same smooth fibre, their sum is the
specialization of its generic trace section. Thus their classes modulo
`M17 ⊗ Q` are opposite. Two independent chosen exceptional directions cannot
be those two intersections of a single bisection. A genuinely rich shared
quadratic cover would instead have to support **different multisections /
independent sections after base change**. The present atlas finds none shared
by two chosen targets; it does not exclude other common covers.

Finally, the question “which curves pass through these fourteen points?” still
depends on the chosen exceptional representatives. Canonical cover equations
do not make that choice invariant under arbitrary rebasing or translation by
generic points. A subgroup-level ancestry invariant needs an additional
definition and proof; this report does not silently supply one.

## 2. Frozen equation-only panel

The [population snapshot](../../artifacts/generated-results/elliptic-curves/rank_triangle_v1/population.json)
contains all 2,080 accepted lower-bound states. Their known jumps 0 through 8
have counts `1241,333,209,143,79,45,25,4,1`. Selection is frozen by
SHA256(`triangle-v1:` + case id): at most two per stratum, preferring different
parent presentations, plus 302. The single jump-8 case is the named 11952
control. Packet hashes and independent lower-bound replays are bound before
arithmetic starts. Failed arithmetic cases are **not replaced**.

Workers receive equations, not points or jump labels. Each field attempt has
15 seconds, one worker, a 2 GiB PARI stack and a 3 GiB RSS cap. Each completed
field receives a separate 3-second provisional BNF probe. All ten probes
failed to finish; no provisional class invariants were obtained. The initial
unhinted 302 field attempt also timed out. A separately recorded 15-second
replay with the already certified equation-support primes completed its local
data. This caching advantage is explicit, not treated as matched runtime.

| Known jump | Presentation | t | BK offset | Field-discriminant digits | Root number |
|---:|---|---|---:|---:|---:|
| 0 | 08f72 | -665/944 | 8 | 120 | -1 |
| 0 | 074d9 | -915/562 | UNKNOWN | — | — |
| 1 | 11952 | 1086/299 | 6 | 117 | +1 |
| 1 | 08f72 | 504/113 | UNKNOWN | — | — |
| 2 | 103b2 | -896/529 | 8 | 111 | -1 |
| 2 | 07ca9 | 293/750 | UNKNOWN | — | — |
| 3 | 103b2 | 808/1277 | UNKNOWN | — | — |
| 3 | x1092-original | 372/1003 | 8 | 91 | +1 |
| 4 | 07ca9 | -108/311 | 6 | 105 | -1 |
| 4 | 11952 | 276/787 | UNKNOWN | — | — |
| 5 | 08f72 | -1079/246 | 5 | 108 | +1 |
| 5 | 07ca9 | 493/706 | 5 | 113 | +1 |
| 6 | 074d9 | -601/716 | 5 | 101 | -1 |
| 6 | 08f72 | -1102/775 | UNKNOWN | — | — |
| 7 | 11952 | -1670/647 | 9 | 116 | +1 |
| 7 | 07ca9 | 371/601 | UNKNOWN | — | — |
| 8 | 11952 | 921/653 | 9 | 107 | -1 |
| 14 | 302 reference | 0 | 12 | 150 | -1 |

The [machine-readable panel](../../artifacts/generated-results/elliptic-curves/rank_triangle_v1/panel.csv)
and per-case arithmetic records retain exact discriminants, conductor,
minimal model and transport, bad-prime factorization, Kodaira data, signature,
prime decomposition `(e,f)`, and local Kummer-image dimensions. At an odd
finite prime the latter is `#primes(K above p)-1`; at 2 it is
`#primes(K above 2)`. The real image has dimension one for the totally real
cubic and zero otherwise. These dimensions do not compute the compatibility
of local conditions inside the global Selmer group.

For the irreducible cubic `f=X³+b2 X²+8 b4 X+16 b6`, the calculation proves
the discriminant factors prime, certifies the maximal order, and checks
`disc(f)=256 disc(E)=[O_K:Z[X]]² disc(K)`. It computes the Brumer–Kramer bound

`rank E(Q) ≤ dim Sel₂(E) ≤ g + offset`,

with ordinary class-group 2-rank g. The offset consists of the signature term
(2 for positive discriminant, 1 otherwise), one for each multiplicative
bad prime with even minimal discriminant valuation, and `#primes(K above p)-1`
for each additive bad prime. This is the equation-derived bound discussed in
[Klagsbrun–Sherman–Weigandt, §3.1](https://arxiv.org/html/1606.07178#S3.SS1).

302's offset 12 is larger than the offsets 5–9 of the ten completed search
cases. It is **not** established to be extreme in the full panel or population:
seven selected search cases remain unresolved. There is no monotone pattern
in the displayed complete offsets. Without g, this is not a comparison of
numerical Selmer capacities, and it supports no enrichment claim. Search
exposure is adaptive, parent presentations differ, and jumps are lower bounds.
No p-values or censored-rank regressions are reported.

The earlier `g≤16 ⇒ rank(11952)=25` criterion remains valid and unresolved.
Incomplete BNF relations are not a certified class-group overgroup, and a
known-point-derived lower bound on g is not independent arithmetic evidence.

## 3. Exact known Kummer group → strict classes → ideals

The [strict dictionary](../../artifacts/generated-results/elliptic-curves/rank_triangle_v1/strict.json)
freshly recomputes the 31 Kummer generators and their finite/local signatures.
In the existing cubic field, the public points give
`γ_i=a_i-d_i² z`, with norm a square. Complete-split good-prime characters
certify their rank 31. Let W be this known Kummer image, G the generic
17-dimensional image, and loc full localization at the recorded bad/dyadic
and real places. The exact ranks are

`dim W=31`, `dim G=17`, `rank loc(W)=21`, `rank loc(G)=17`.

Consequently `V=ker(loc|W)` has dimension ten and `V∩G=0`. Its image in W/G
is canonical relative to these specified subgroups and local conditions:

`0 → V₁₀ → W/G → loc(W)/loc(G) ≅ F₂⁴ → 0`.

The four chosen coordinates and a splitting are not canonical. Neither V nor
W is asserted to exhaust the strict Selmer kernel or full Selmer group.
The retained equation is `dim Sel₂=21+c_S`, with `c_S≥10` and no upper bound.

The table gives all ten strict words in the fixed geometry basis. Brackets
list the indices to add modulo 2. Each row's full Kummer product β has
principal ideal `(β)=J_i²`; the exact HNF ideals are in the dictionary.

| Half-ideal | Exceptional indices | Generic correction indices |
|---|---|---|
| J0 | 1,6,7,11,14 | 2,3,5,6,7,9,10 |
| J1 | 7,9,11,13 | 1,2,3,4,7,8,9,11,12,13,15,16,17 |
| J2 | 2,8,12,14 | 1,3,4,5,7,9,11,12,13,14,15,16,17 |
| J3 | 3 | 3,7,15,16 |
| J4 | 6,7,8,9,10,11,12,14 | 1,4,6,8,11,12,13,14,16,17 |
| J5 | 4,6,7,11,12,13 | 4,5,6,7,8,9,10,14 |
| J6 | 9,10,12,13,14 | 2,3,4,5,6,8,10,16,17 |
| J7 | 6,7,8,10,12,13,14 | 1,4,6,7,11,13,15 |
| J8 | 5,6,7,10,11,12 | 5,6,7,10,13,14,15,16 |
| J9 | 6,7,8,9,11,12 | 4,5,6,9,11,14,15 |

These are **mod-2 Kummer equalities**, not integral relations among points.
The ten ordinary half-ideal classes are independent and the strict unit
kernel is zero, as certified by the retained
[descent anatomy](CURVE302_DESCENT_ANATOMY_2026-09-10.md). Its checker was rerun:
all ten ideal-square identities, 100 Artin entries, the remaining-class
argument, and two deliberately altered ideals pass/reject as expected.

In this chosen complement, only E3 and E7 individually admit a generic
correction into V. Explicitly:

* `E3 + M3 + M7 + M15 + M16` gives J3.
* `E7 + M1 + M2 + M3 + M4 + M9 + M11 + M13 + M14` gives J1+J4+J7.

All fourteen targets nevertheless have genus-one carriers, including both
zero and nonzero local quotient classes. Thus this genus-one test does not
separate strict from non-strict exceptional axes. With every tested cover
singleton on the chosen target list, there is no demonstrated shared-cover /
strict-class clustering. Mixed strict words have not all been given separate
minimal-ancestry searches.

## Replay and next mathematical gates

From the repository root, with the recorded Sage environment:

```sh
sage -python research/elliptic-curves/cas/verify_rank_triangle.sage
sage -python research/elliptic-curves/cas/verify_curve302_descent_anatomy.sage --negative-controls
python3 research/elliptic-curves/cas/report_rank_triangle.py
python3 research/scripts/render_status.py
python3 research/scripts/audit_status.py
```

The original producers and frozen protocols are retained. The panel producer's
`freeze` deliberately refuses to replace the snapshot; `run` resumes only
unattempted cases and does not retry failed ones. Development serialization
failures remain in separate logs. The first combined verification timed out;
the diagnosed repeat was stopped, and the final replay supplied already
proved discriminant primes to PARI's factor cache, avoiding unnecessary
refactorization. All 792 identities and eleven local rows then passed within
the original two-minute verifier cap. The large running searches were untouched.

The next discriminating geometric gate is **genus-zero incidence or multiple
independent anti-invariant sections on one common quadratic cover**, not the
mere existence of genus-one carriers. Any enlarged census needs a frozen
budget and an explicitly defined equivalence on subgroups/targets. The
arithmetic gate is completion of missing equation support and genuinely
independent class-group/Selmer upper information. The exact strict dictionary
is now available to test a positive geometric construction when one exists.

The canonical-height-ball benchmark stays paused. A further long BNF run,
broader multisection sweep, or new point-search policy is not scheduled by
this report. Rank ≥31 and rank ≥25 remain the respective fibre statements.
