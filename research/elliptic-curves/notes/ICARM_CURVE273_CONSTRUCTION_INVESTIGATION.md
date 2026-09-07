# ICARM curve 273: construction and family investigation

Status: **exactly excluded from the published R17 chart; broader source audit
ongoing**.  This note does not identify a construction for curve 273.

## Bottom line

The curve itself and 30 independent points are reproducible: see
[`ICARM_CURVE273_RANK30.md`](ICARM_CURVE273_RANK30.md).  The construction is not
yet reproducible.  The current public ICARM entry attributes the find to
Claude, with Levent Alpöge and Ava Howell, but neither that entry nor the public
leaderboard discussion gives a family equation, search parameter, or
specialization certificate.

The hypothesis that curve 273 is a direct rational specialization of the
published rank-17 chart is now exactly false: the primitive degree-24 equation
`j_R17(t)=j_273` has an irreducible degree-24 reduction modulo `367`, and hence
no rational root.  A broader common construction lineage remains possible.
Establishing it requires one of the following:

1. a discoverer-supplied construction record;
2. a different explicit family and a rational parameter whose specialization
   is isomorphic over `Q` to curve 273; or
3. an equivalent exact specialization certificate, including the transport of
   the generic sections.

The exact published-chart exclusion is replayed in
[`ELKIES_BISECTION_VISIBILITY_AND_RECORD_CURVES.md`](ELKIES_BISECTION_VISIBILITY_AND_RECORD_CURVES.md).

## Public provenance recovered on 2026-08-20

- The [ICARM curve 273 page](https://elliptic-rank.icarm.cloud/curve/273)
  records the equation, 30 points, the submitter name `ranksunbounded`, and the
  attribution "Claude, with Levent Alpöge and Ava Howell".
- The [public ICARM Zulip topic](https://icarm.zulipchat.com/#narrow/channel/519875-general/topic/Elliptic.20Curve.20Rank.20Leaderboard/near/603443505)
  records that the initial 29-point display was a parser omission, after which
  the thirtieth submitted point was restored.  It also records the conditional
  analytic-rank discussion, but no construction recipe.
- Dujella's maintained [rank-30 page](https://web.math.pmf.unizg.hr/~duje/tors/rk30.html)
  reproduces the curve and its 30 points.  As rechecked on 2026-08-23, that
  page and Dujella's rank-history page now attribute the record to Levent
  Alpöge and Ava Howell; this updates the public attribution but supplies no
  construction data.

The earlier pinned source audit
[`elliptic_rank30_public_source_audit.json`](../../archive/elliptic-curves/artifacts/generated-results/elliptic_rank30_public_source_audit.json)
predates curve 273.  It remains useful for the rank-29 construction history,
but its conclusion that no public rank-30 curve was available has been
superseded by the 2026-08-20 record.

## What is known about the preceding record family

The public rank-29 announcement says that its search used a rank-17 elliptic
fibration on the same K3 surface used for the rank-28 record, sieved rational
specializations, and searched outside the generic `Z^17`.  The announcement
did not publish the family equation or the 17 generic sections.

The independent reconstruction programme in [`../../elkies-k3/`](../../elkies-k3/)
has identified the recovered rank-17 Mordell--Weil lattice with the Shimura
datum

```text
quaternion discriminant D = 6
level M = 79
rank-17 lattice determinant = 948 = 2^2 * 3 * 79.
```

The most economical current explicit route is the `E6/MW3` neighbour with

```text
ADE = E6 + A3^2 + A1^2
MW rank = 3
fibres = IV* + I4 + I4 + I2 + I2 + 4 I1
reduced MW determinant = 79/16.
```

An exact genus check plus Nikulin's rank-versus-discriminant-length uniqueness
theorem first proved that its `U`-extended frame is in the same integral
Neron--Severi isometry class as the recovered rank-17 fibration.  The stronger
explicit lattice transport is now also pinned.  The actual path is

```text
rank17 --q=90--> MW7 --q=4--> MW4 --q=4--> E6/MW3,
```

and every raw child frame plus the composite determinant-one isometry is
checked by
[`../../elkies-k3/scripts/verify_e6_neighbor_chain.sage`](../../elkies-k3/scripts/verify_e6_neighbor_chain.sage).
See [`../../elkies-k3/E6_NEIGHBOR_CHAIN.md`](../../elkies-k3/E6_NEIGHBOR_CHAIN.md).

What remains missing is no longer an abstract or integral-lattice
identification.  The source audit, distinct from the reverse MW17 backtracks,
now identifies the correct upstream class of genuine models: the
Dolgachev--Kumar `E7+E8` fibration of
MW rank two and regulator `474`.  Exact binary-form/glue classification leaves
three compatible determinant-948 Kumar frames, with height Grams
`[[5/2,1],[1,190]]`, `[[4,0],[0,237/2]]`, and
`[[21/2,3],[3,46]]`.  The middle one is uniquely involutive.  Exact elliptic
quotient labels and its action on the discriminant group identify that
involution as `w2=w237`, the hyperelliptic involution of Elkies's genus-two
curve, so the middle frame is the canonical involutive H2 anchor.  The
corrected H3 source marking is described next.  See
[`../../elkies-k3/KUMAR_E7E8_BACKTRACK.md`](../../elkies-k3/KUMAR_E7E8_BACKTRACK.md).
The route directions and the qualified meanings of endpoint, source, lattice
child, equation child, and specialization endpoint are fixed in
[`../../elkies-k3/CONSTRUCTION_ROUTES.md`](../../elkies-k3/CONSTRUCTION_ROUTES.md).

### Corrected H3 source and low-degree route

The third frame, not the involutive middle frame, is now the preferred source
marking.  Pulling `H21` back to the Elkies--Kumar `H92` chart identifies the
level-474 component exactly and normalizes it over `QQ` to

```text
y^2 = -27*x^6 + 198*x^4 - 171*x^2 + 576.
```

At the published non-CM point `(x,y)=(13/7,12048/343)`, exact weighted-Igusa
matching gives rational H92 and H21 presentations whose unmarked short
Weierstrass K3 models are isomorphic over `QQ`.  Their oriented Hilbert-cover
values have the same nonsquare class `-52203427`, and their oriented-coordinate
ratio is rational.  The marked section descent is now unconditional: the
split source fibers have squarefree trivial-lattice determinant 21, so the
height-`21/2` generator and both signs are individually rational.  Moreover,
modular nonflex conversion followed by structured CRT/LLL has recovered exact
H92 coordinates for that section.  The reconstruction uses 204 good primes
and a 1,945-bit modulus; `x` has degrees `(10,12)`, `y` has degrees `(15,18)`,
and exact substitution proves the characteristic-zero Weierstrass identity.
An exact marked-fiber incidence fixes the sign of `y`.  The certificate is
[`../../artifacts/generated-results/elkies-k3-h92-p1-lift.json`](../../artifacts/generated-results/elkies-k3-h92-p1-lift.json),
SHA-256
`c323bf6346bb239934a5a2d8b1a3f4067e70e993d2e4eb32aaa30f469fca6397`.
This supersedes historical SHA-256
`0602c3b199629c6f460c9b7c728e048822418ecf85bf54807852be3d97b66616`,
which lacked only the orientation-incidence block and used the older H21
status label.  The remaining source gate is more specific: construct the
chord and an explicit basis of `H0(O+(-P1)-F)` giving the desired
`E7+E8 -> E8+E6` q=6 pencil.

The complementary H92 source direction is rational as well.  Its pinned
split `D6+A8+A1` fibration has a section of height `23/18`; adjoining the
section gives determinant `92` and Smith factors `(2,46)`.  An exhaustive
check of all 92 discriminant classes finds no nonzero isotropic class, so the
source lattice is saturated.  The rational two- and three-neighbor sequence
then makes the height-`46` generator `P2` on the final `E7+E8` fibration
individually rational.  This is certified by
[`../../elkies-k3/scripts/verify_h92_section_descent.sage`](../../elkies-k3/scripts/verify_h92_section_descent.sage)
and
[`../../artifacts/generated-results/elkies-k3-h92-section-descent.json`](../../artifacts/generated-results/elkies-k3-h92-section-descent.json),
SHA-256
`fe525f75fa87c31afb34755fe63fc778349d2843010eb5c9b17ce6d8b8712e40`.
Explicit `P2` coordinates on the short H92 model remain unrecovered; they are
not needed for the first `P1` q=6 chord, but may be needed for later section
transport.

The defining H3 source family is now equation-level, not only a Humbert-locus
identification.  Over

```text
Y^2=-27*X^6+198*X^4-171*X^2+576,
```

the exact normalization functions recover the H92 chart by
`r=(a+Y0)/2`, `s=2/(Y0-a)`, after inverting the published linear-fractional
`X`-map and applying the pinned `Y` multiplier.  Function-field substitution
proves the degree-21 H21/H92 component equation.  Composing the five H92
coefficients gives

```text
v^2=u^3+(A1*tau^3+A*tau^4)*u+(B1*tau^5+B*tau^6+B2*tau^7),
```

the exact `E7+E8/MW2` source family.  The verifier
[`../../elkies-k3/scripts/export_h3_level474_source_family.sage`](../../elkies-k3/scripts/export_h3_level474_source_family.sage)
also specializes `(X,Y)=(13/7,12048/343)` to the pinned H92 point and short
model.  Its artifact is
[`../../artifacts/generated-results/elkies-k3-h3-level474-source-family.json`](../../artifacts/generated-results/elkies-k3-h3-level474-source-family.json),
SHA-256
`8f5afd11e1d8979d57cb1a569833309f9664c19cd47194af0581a5cbbf8f1d59`.
This identifies the source family of the proposed construction.  The later
published rootless `MW17` chart can now be tested directly, and its exact
degree-24 `j`-recognition equation excludes curve 273 at every rational
parameter.  Identifying a broader shared construction would therefore require
a genuinely different family or an isogeny-level explanation.
<!-- status-consumer: EC-K3-H3-SOURCE a4bb40c9c9d0ff09 -->

The rational points of the level-474 base are now determined globally.  On
the Q-isomorphic smaller model

```text
y^2=-3*x^6+22*x^4-19*x^2+64,
```

an exact two-cover descent has three locally soluble classes.  Factoring the
sextic over `L=Q(a)`, with `-3*a^3+22*a^2-19*a+64=0`, gives two elliptic
covers over `L`; their twisting factors are

```text
(263403*a^2-123771*a+818724)/4,
(45*a^2-21*a+144)/16.
```

The exact Magma certificate
[`../../elkies-k3/scripts/prove_h3_level474_rational_points.m`](../../elkies-k3/scripts/prove_h3_level474_rational_points.m)
computes finite-odd-index pseudo-Mordell--Weil groups for both covers and
runs elliptic Chabauty at `p=41`.  In each case the Chabauty upper bound equals
the known number of cover points (`2` and `8`) and its residual index is `1`.
Thus their complete rational `x`-image sets are `{0}` and
`{-13/7,-1,1,13/7}`.  Exact substitution supplies both signs of `y`; the
leading coefficient `-3` is nonsquare, so there are no rational points at
infinity.  Returning to the published coordinates gives exactly

```text
(X,Y) = (0,+-24), (+-1,+-24),
        (+-13/7,+-12048/343).
```

The pinned terminal output is
[`../../artifacts/generated-results/elkies-k3-h3-level474-rational-points.txt`](../../artifacts/generated-results/elkies-k3-h3-level474-rational-points.txt),
SHA-256 `7008b9536d82c03ad4b568324192b954a67bebe9120880d5a5327b36d080df02`.
The historical
[`../../elkies-k3/scripts/sieve_h3_level474_rational_points.sage`](../../elkies-k3/scripts/sieve_h3_level474_rational_points.sage)
remains a separate bounded quotient-sieve cross-check.
<!-- status-consumer: EC-K3-H3-PTS 8f0a27c947843b4a -->

On the lattice and chamber side, the corrected source now has the exact path

```text
H3 E8+E7/MW2 --q6--> E8+E6/MW3 --q8--> D13/MW4
 --q24--> D12/MW5 --q6--> A11/MW6 --q8--> 2A5/MW7
 --q4--> 3A3/MW8 --q4--> A3+2A2/MW10
 --q4--> 5A1/MW12 --q4--> 4A1/MW13 --q4--> 3A1/MW14
 --q4--> 2A1/MW15 --q4--> A1/MW16 --q6--> rootless/MW17.
```

Every displayed factor presentation is a certified nef old-fiber-degree-two
pencil.  The selected D13-to-MW17 composite transport is integral and
unimodular.  The authoritative source, chamber, and transport replays are
[`../../elkies-k3/scripts/verify_h3_noncm_q6_source_anchor.sage`](../../elkies-k3/scripts/verify_h3_noncm_q6_source_anchor.sage),
[`../../elkies-k3/scripts/verify_h21_q6_section_descent.sage`](../../elkies-k3/scripts/verify_h21_q6_section_descent.sage),
[`../../elkies-k3/scripts/lift_h21_p1_modular.sage`](../../elkies-k3/scripts/lift_h21_p1_modular.sage),
[`../../elkies-k3/scripts/analyze_h3_rank_growing_degree2_chain.sage`](../../elkies-k3/scripts/analyze_h3_rank_growing_degree2_chain.sage),
[`../../elkies-k3/scripts/analyze_h3_mw10_to_rootless_chambers.sage`](../../elkies-k3/scripts/analyze_h3_mw10_to_rootless_chambers.sage),
[`../../elkies-k3/scripts/verify_h3_d13_to_mw17_path.sage`](../../elkies-k3/scripts/verify_h3_d13_to_mw17_path.sage),
and
[`../../elkies-k3/scripts/verify_rank17_to_h3_reverse_transport.sage`](../../elkies-k3/scripts/verify_rank17_to_h3_reverse_transport.sage).
The last checker supplies the previously missing determinant-one isometry to
the pinned recovered rank-17 Gram and a lossless inverse NS transport to H3.
It also proves that the dominant D13 lattice record and component-nef D13
equation record require a non-`U`-preserving marking bridge rather than an
unqualified identification.
This supersedes the H2/q80 deformation as the preferred source transport;
the latter remains a valuable independent downstream comparison.
<!-- status-consumer: EC-K3-H3-D13-MW17-LATTICE-CHAIN 2c6a2a36699933ab -->

The older H2 route remains a separate downstream comparison.  Its height-4
and height-`237/2` directions define the Humbert intersection
`H8 cap H237`, and Elkies--Kumar's discriminant-8 ancillary data give the
complete two-parameter Kumar coefficients, Clebsch--Igusa invariants, and
oriented double cover in coordinates `(r,s,z)`.  Recovering the
discriminant-237 equation and its marked function-field map would complete
that H2 comparison, but it is no longer the missing source datum: the source
curve has instead been normalized exactly on `H21 cap H92` through H3.  See
[`../../elkies-k3/KUMAR_E7E8_BACKTRACK.md`](../../elkies-k3/KUMAR_E7E8_BACKTRACK.md).

The same split is now visible intrinsically in the recovered rootless
Neron--Severi lattice.  A short primitive square-`-4`, divisibility-four class
has the exact `w2` discriminant action `475 mod 948`; it splits off with a
rank-18 determinant-237 orthogonal complement in the Kumar fixed genus.  The
canonical proof and verifier are in
[`../../elkies-k3/KUMAR_E7E8_BACKTRACK.md`](../../elkies-k3/KUMAR_E7E8_BACKTRACK.md)
and
[`../../elkies-k3/scripts/verify_rank17_h8_split.sage`](../../elkies-k3/scripts/verify_rank17_h8_split.sage).
This confirms the intrinsic H2/`H8 cap H237` backtrack independently of the
long q80 path; it does not identify that backtrack with the corrected H3
source polarization.
It also rules out a direct q8 fiber in the distinguished orthogonal
complement. The first such fibers occur at q9: exactly thirteen classes up to
sign, with child root ranks 8--10 and MW ranks 9--7; none is rootless. The
best root data are `(8,20,144)` and `(8,22,108)`, certified by
[`../../elkies-k3/scripts/analyze_rank17_h8_q9_fibers.sage`](../../elkies-k3/scripts/analyze_rank17_h8_q9_fibers.sage).

The parallel H2 Humbert-8 construction has also been replayed one step farther
back.  Elkies--Kumar first use the exact `D9+E7` family

```text
Y^2 = X^3 + T*(r+(2*r+1)*T)*X^2
            + 2*r*s*T^4*(T+1)*X + r*s^2*T^7
```

and then take the two-neighbor parameter
`U=(X+s*T^3)/T^4`.  The resulting pointed quartic has rational point
`(T,V)=(0,s)` and binary-quartic Jacobian equal exactly to the displayed
`E7+E8` Kumar model.  See
[`../../elkies-k3/scripts/verify_humbert8_d9e7_two_neighbor.sage`](../../elkies-k3/scripts/verify_humbert8_d9e7_two_neighbor.sage).
Thus the q=8 signal has a precise structural meaning: `det(D9+E7)=8`, and the
height-four Humbert direction is absorbed into this pre-neighbor fiber frame.
On the desired determinant-948 curve, its remaining generic MW section is
torsion-free and has forced height `948/8=237/2`.  The missing `H237` equation
is equivalently the rank-one MW-jump locus of this explicit two-parameter
family.

The second quartic point `(0,-s)` gives the generic height-four Kumar section
directly.  At the CM-`43` point, if `P1,P2` are the two exact height-`5/2`
sections and `P3` is this height-four section, then the level-79 generator is

```text
Q79 = 4*P1 - 5*P2 + P3.
```

The exact height pairing gives `height(Q79)=237/2` and `P3.Q79=0`; its
Weierstrass coordinates have denominator pattern `h^2,h^3` with
`deg(h)=58`.  Thus the original large-section obstruction is now reproduced
from a completely marked CM boundary point, and the q60/q80 transport has a
specific horizontal class rather than only an isometry class.

The direct local deformation of this pole-58 section has now been audited far
enough to reject an initially tempting shortcut.  Linearizing its cleared
identity over a good finite field gives a `361 x 362` matrix: both the section
block and full matrix have rank `358`.  The quadratic obstruction kills two
nilpotent section directions but leaves the whole `(dr,ds)` tangent plane.  A
canonical right inverse at the next order produces the binary cubic

```text
 446402445291973586382160910173306149076992 * dr^3
-2642585311483808090087377295677498775040000 * dr^2*ds
-6535448784280617898603638836147191706250000 * dr*ds^2
+16726681288628536079919155997729849853515625 * ds^3.
```

Twenty good primes near `10^6` reconstruct this canonical-slice cubic uniquely
with a 399-bit CRT modulus.  It has the rational factor

```text
dr/ds = 223593125/30934224,
```

but this is **not** an `H237` tangent.  The second correction is defined only
modulo the four-dimensional first-order kernel.  Polarizing the quadratic
residual shows that these kernel additions map to the third-order cokernel
with rank one; augmenting by the displayed cubic still has rank one.  Hence
the entire cubic, including its rational factor, is absorbed by the correction
gauge for generic `(dr,ds)`.  The high-section deformation remains
nontransverse through third order.  The later full chamber audit also rejects
the apparent marked q8/q60 CM-43 shortcut: both classes have extra fixed CM
sections and reduce to the old fiber.  Replay the modular audit and the
canonical-slice CRT
certificate with
[`../../elkies-k3/scripts/compute_cm43_h237_tangent.sage`](../../elkies-k3/scripts/compute_cm43_h237_tangent.sage)
and
[`../../elkies-k3/scripts/verify_cm43_h237_tangent_crt.sage`](../../elkies-k3/scripts/verify_cm43_h237_tangent_crt.sage).

The H2 boundary is exact rather than conjectural.  The CM points at
`t=infinity` and `t=+/-2` give primitive Kumar-frame closures with respective
root/MW data `E8+E8+A2, MW=0` and `E8+E8, MW=diag(4,6)`.  Their standard Inose
equations are

```text
Y^2 = X^3 + T^5*(T-1)^2,
Y^2 = X^3 - 51*T^4*X + T^5*(T^2-92*T+1),
```

with the second fixed up to rational quadratic twist.  Thus the next family
calculation is a local deformation/interpolation from known boundary models,
not a search for an unspecified K3 surface.  The high Kumar generator should
also not be imposed directly: its height is `237/2`, its intersection with the
zero section is `58`, and `u` records its sign/descent.  The height-4 generator
is the small direction that descends to the rational `t`-line.

A third boundary anchor is now exact in the Humbert-8 chart.  The added class
of the rational Picard-rank-20 surface is the Gross norm-`43` vector

```text
beta = -11*i-j-5*k.
```

It transports to `(169,167,-128)` in the generic transcendental lattice, with
square `-40764`, divisibility `948`, and complement
`[[22,1],[1,2]]`.  Two height-`5/2` section loci in the Humbert-8 plane meet at

```text
r = -1225/722,
s = -93312/442225,
z^2 = -43*(11664/6859)^2.
```

The two short sections have height Gram
`[[5/2,-1/2],[-1/2,5/2]]`.  Direct point counts at fourteen good primes match
the weight-three CM-`43` coefficients, including zero at the tested inert
primes.  This makes the discriminant-43 surface a concrete interpolation
anchor for the H2/`H8 cap H237` comparison.  The finite Frobenius fingerprint
is not itself an
explicit characteristic-zero isomorphism between the rational MW3 equation
and the Kumar presentation, so the inverse neighbor/normalization gate remains.
Replay all of these checks with
[`../../elkies-k3/scripts/verify_cm43_humbert8_anchor.sage`](../../elkies-k3/scripts/verify_cm43_humbert8_anchor.sage).

An earlier exact H2 neighbor forced to use the level-79 direction is:

```text
Kumar E7+E8/MW2 --q=60--> E8+E6/MW3.
```

A constrained norm-120 enumeration is complete with only 56 sign-pairs.  The
selected `(a,b)=(5,12)` frame has reducible root type `E8+E6`, so the expected
minimal generic fiber configuration is `II*+IV*+6 I1`; its reduced MW height
Gram is `[[4,0,0],[0,20/3,1],[0,1,12]]`, of determinant 316.  This is a direct
lattice bridge from the H2 frame.  The corrected H3 source instead has the
smaller q6 degree-two entrance above, so this q60 bridge is retained as a
comparison and marking check rather than the preferred source transport.

Shioda's formula reduces the three zero-section intersections to `0,2,4`.
After placing `IV*` at zero and `II*` at infinity, the whole ambient equation
has only the six coefficients in

```text
Y^2 = X^3 + T^3*(a0+a1*T)*X
            + T^4*(b0+b1*T+b2*T^2+b3*T^3).
```

This is a concrete small reconstruction system, not a request for another
large blind fibration scan.

The CM anchors now give exact boundary conditions in this chart.  Transport
of their primitive lattice closures through the `q=60` neighbor gives
`E8+E7+A2` with MW height 4 at discriminant 24 and `E8+E8+A2` with MW rank
zero at discriminant 3.  Consequently the former lies on `b0=0`, where the
zero fiber has orders `(3,5,9)`, and the latter lies on `a0=b0=0`, where the
orders are `(4,5,10)`.  The discriminant-3 coefficient point is explicitly
`(a0,a1,b0,b1,b2,b3)=(0,0,0,1,-2,1)`.  See
[`../../elkies-k3/scripts/verify_e8e6_cm_boundaries.sage`](../../elkies-k3/scripts/verify_e8e6_cm_boundaries.sage).

The section complexity has now been reduced further without a search.  The
`0,2,4` pole profile is optimal, but the height-4 polynomial section admits an
exact three-parameter normal form: its endpoint jets determine all six
surface coefficients linearly from `y(P1)^2-x(P1)^3`.  Thus recovering the
rank-19 locus requires imposing only the remaining `P.O=2` and `P.O=4`
sections.  The same chart exposes the CM enhancement through the factor
`ell=c3-epsilon*c1`, with divisibilities `ell|a0`, `ell^2|b0`, and
`ell^3|b3`.  See
[`../../elkies-k3/scripts/verify_q60_height4_normal_form.sage`](../../elkies-k3/scripts/verify_q60_height4_normal_form.sage).

The discriminant-24 CM seed in this chart has also been recovered exactly:

```text
Y^2 = X^3 - (4096/19683)*T^3*(7*T+2)*X
            - (262144/14348907)*T^5*(T^2+34*T+19),
Delta = -(2^40/3^27)*T^9*(T-1)^3*(T^2+71*T+32).
```

It has fibers `II*+III*+I3+2 I1`, a polynomial height-4 section, and the
exact transported discriminant-24 frame.  This is now a concrete local seed
and an independent marking check, rather than only an abstract Inose boundary
condition.  See
[`../../elkies-k3/scripts/verify_q60_delta24_anchor.sage`](../../elkies-k3/scripts/verify_q60_delta24_anchor.sage).

Backtracking over the neighbor presentation materially improves this route.
The five proper `q=60` factorizations show that `(4,15)` retains two MW
directions at discriminant 24, whereas `(5,12)` retains only one.  More
importantly, an exact dominant-weight enumeration tested all 2,869 proper
level-79 presentations through `q=80`.  There was no root-stable chart in
this bounded range, but the best compact additive frame is

```text
q=80, (a,b)=(4,20)
generic:    E6+D5+A3,     MW rank 3, optimal P.O=0,0,3
Delta=-24: E6+D5+A3+A1,  MW rank 3, optimal P.O=0,0,0.
```

Only one root is added and no generic MW direction is lost.  Its exact frame
is
[`../../elkies-k3/data/fibrations/kumar_q80_e6_d5_a3_mw3_frame.txt`](../../elkies-k3/data/fibrations/kumar_q80_e6_d5_a3_mw3_frame.txt),
verified by
[`../../elkies-k3/scripts/classify_kumar_cm_frame_extensions.sage`](../../elkies-k3/scripts/classify_kumar_cm_frame_extensions.sage).

This chart has a four-parameter Weierstrass ambient with `I1*`, `I4`, and
`IV*` at `0,1,infinity`:

```text
A=T^2*(-3+p*T+q*T^2+(-3*d^2+3-p-q)*T^3),
B=T^3*(2+b1*T+b2*T^2+b3*T^3+b4*T^4+e*T^5).
```

The `I4` jets determine `b1,...,b4` linearly, and
`Delta=T^7*(T-1)^4*R5(T)` up to scale.  The discriminant-24 enhancement is
just `disc_T(R5)=0`, an extra `I2`.  These formulas are exact in
[`../../elkies-k3/scripts/verify_q80_ambient_normal_form.sage`](../../elkies-k3/scripts/verify_q80_ambient_normal_form.sage).
This was the leading noncollapsed H2 chart before the H3 branch was
identified.  It remains a useful downstream cross-check, but it is no longer
the preferred source reconstruction.  The claimed smooth `GF(7)`
rank-19 seed was too strong:
after resolving the forced component coordinates, its Jacobian has rank `37`
in `39` variables.  The exact quadratic obstruction is

```text
tau0^2 - tau0*tau1 + tau1^2
  = (tau0+2*tau1)(tau0-3*tau1)  over GF(7).
```

Thus the marked CM point is a crossing of two tangent branches, not a smooth
one-dimensional point in that chart.  Fixing `p=4+h`, both directions lift
uniquely through the bounded formal order `h^20`.  A coordinatewise Pade test
with numerator and denominator degrees at most ten recovers only 24/39 and
22/39 active coordinates, respectively, so neither branch has yet been
algebraized.  The exact cone and bounded extension are checked by
[`../../elkies-k3/scripts/verify_q80_rank19_deformation_gf7.sage`](../../elkies-k3/scripts/verify_q80_rank19_deformation_gf7.sage)
and
[`../../elkies-k3/scripts/extend_q80_rank19_branches_gf7.sage`](../../elkies-k3/scripts/extend_q80_rank19_branches_gf7.sage).
Neither finite-field branch has yet been algebraized or identified with the
characteristic-zero Shimura curve.

The CM boundary point itself is exact over `QQ`:

```text
A=T^2*(-3+9/4*T-9/4*T^2+9/4*T^3),
B=T^3*(2-315/32*T+9*T^2-9/16*T^3-27/32*T^5).
```

Three selected marked sections are defined over
`QQ(sqrt(-3),sqrt(-6))`.  The characteristic-zero resolved tangent cone is
`(tau0-8/87*tau1)(tau0-1/12*tau1)`; its two lines reduce to the GF(7) slopes
`5` and `3`.  This is an exact CM boundary and tangent-direction certificate,
not an algebraized one-parameter branch.  See
[`../../elkies-k3/scripts/verify_q80_cm24_rational_model.sage`](../../elkies-k3/scripts/verify_q80_cm24_rational_model.sage).

With `p=9/4+h`, their normalized surface tangents are

```text
(d',p',q',e')=(8/87,1,-24/29,-45/116),
(d',p',q',e')=(1/12,1,-45/52,-261/832).
```

Both full marked systems lift exactly over the CM compositum through order
14; see
[`../../elkies-k3/scripts/extend_q80_cm24_branches_qq.sage`](../../elkies-k3/scripts/extend_q80_cm24_branches_qq.sage).
This remains a formal-series certificate, not an algebraic family.

The CM-43 q8/q60 comparison is now closed exactly.  In the marked Kumar
divisor basis, the q60 class first reduces to `Q79+4O-43F`, then has fixed
parts `4(P1-P2)` and `P3-P2`, leaving `F`.  The q8 class likewise has fixed
sections `-P2` and `P1-P2`, again leaving `F`.  Thus

```text
D60_raw = F + initial_fixed + 4*(P1-P2) + (P3-P2),
D8_raw  = F + initial_fixed + (-P2) + (P1-P2).
```

This is a boundary automorphism/fixed-part certificate, not an equation seed.
See
[`../../elkies-k3/scripts/verify_cm43_marked_divisor_transport.sage`](../../elkies-k3/scripts/verify_cm43_marked_divisor_transport.sage).

### Lower-q backtrack and rational discriminant-43 anchor

A rerun from the recovered rank-17 lattice found the substantially shorter
exact ancestry

```text
MW17
  -- q=25 --> MW7
  -- q=4  --> MW4
  -- q=4  --> MW3
  -- q=4  --> MW2, E6+D4+2A2+A1.
```

Every transition and the determinant-one composite transport are exact.  The
terminal fibration has fibers

```text
IV* + I0* + 2 I3 + I2 + 2 I1
```

and two displayed sections with height Gram
`(1/6)[9,2;2,18]`, of determinant `79/18`.  More importantly, this endpoint
has been reconstructed over `QQ`; literal substitution verifies its equation,
sections, fiber valuations, component profiles, and pair intersection.  The
model and the full transport replay are checked by
[`../../elkies-k3/scripts/verify_mw2_e6d4a2a2a1_qq.sage`](../../elkies-k3/scripts/verify_mw2_e6d4a2a2a1_qq.sage)
and
[`../../elkies-k3/scripts/verify_mw2_rank17_transport.sage`](../../elkies-k3/scripts/verify_mw2_rank17_transport.sage).

The rational endpoint has a third section.  Adding it gives a saturated
Picard-rank-20 Neron--Severi lattice of discriminant `43`; see
[`../../elkies-k3/scripts/verify_picard20_ns_extension.sage`](../../elkies-k3/scripts/verify_picard20_ns_extension.sage).
This explains the repeated small-neighbor signal: on this enlarged
discriminant-43 lattice, the exact path `q=8,q=9` reaches an MW-rank-one
fibration.  It is not a factorization of the determinant-948 Kumar `q=80`
neighbor; the two computations live on different Neron--Severi lattices.

The first `q=8` move now has an exact geometric interpretation.  Translating
its isotropic vector back to the displayed divisor basis and replaying 25
Weyl reflections gives a primitive isotropic class nonnegative on every
displayed effective curve.  Its old-fiber degree is two and its horizontal
restriction is

```text
3 O - S,
```

where `S` is the third rational section.  Thus the neighbor is the line pencil
through `S`, generated fiberwise by `x-Sx` and `y-Sy`; it is not an
undirected large-coefficient section search.  The exact raw class, reflection
certificate, and reduced class are replayed in
[`../../elkies-k3/scripts/verify_picard20_ns_extension.sage`](../../elkies-k3/scripts/verify_picard20_ns_extension.sage).
Determining the required base-dependent normalization of that pencil remains
a specialization-only experiment, not the next generic construction gate.

This `q=8` move is intrinsically specialization-only: its lattice vector has
coefficient `-1` in the added `S` direction, and the reduced divisor literally
contains `-S`.  The determinant-948 generic rank-19 lattice has only `P1,P2`
in this fibration.  Indeed the determinant ratio between the three- and
two-section height lattices is exactly

```text
(43/216) / (79/18) = 43/948,
```

the orthogonal specialization height contributed by `S`.  Consequently the
`q=8` pencil is a useful Noether--Lefschetz boundary model, not by itself the
generic rank-17 inverse path.  To inform the construction it must either be
deformed off the `S` divisor or replaced by a separate small-neighbor witness
already present in the determinant-948 lattice.

The conspicuous generic `q=8` witnesses already present in the bounded
`A13+A1` continuation do not supply that replacement.  Exact root and glue
analysis of its retained frame 007 gives `A10+A2+2A1` and reduced height Gram

```text
(1/66) [ 79 -17   1]
       [-17 106  19]
       [  1  19 259].
```

Negating the first basis vector recovers the original A10 branch Gram exactly.
The two retained raw q8 witnesses need 39 and 40 old-chamber reflections but
reduce to the same primitive nef class
`(2,2,-32,23,22,15,-8,-4,-2,-23,6,-4,11,11,-18,17,-15,15,23)`.
An exact rank-three section bound and bisection norm argument prove nefness.
Thus this generic `q=8` is a genuine fibration but an A10 return loop, while
the discriminant-43 `q=8` is a special boundary move that collapses to the old
fiber; their common norm does not identify a missing generic path. Replay with

```bash
sage elkies-k3/scripts/analyze_mw3_branch.sage \
  --frame artifacts/generated-results/elkies-k3-lower-q-path-search/a13-mw3-child-search-frames/frame-007.txt \
  --name generic-q8-from-a13
sage elkies-k3/scripts/analyze_a13_q8_neighbors.sage
```

The first geometric inverse step from the rational MW2 model is also explicit.
The chord pencil through `O+P1`, followed by the binary-quartic Jacobian,
produces over `QQ` the fibration

```text
I1* + I4 + 2 I3 + 2 I2 + 3 I1,
root system D5+A3+2A2+2A1.
```

Three rational sections are independently certified by specialization.  The
obvious fourth polynomial point is only a saturation correction: it satisfies
`2*P4=P1`.  Since the source K3 has geometric Picard rank 20 and the new root
rank is 14, Shioda--Tate predicts geometric MW rank four, but a fourth
independent `QQ(u)` section has not been proved.  The verifier must therefore
distinguish `three independent rational sections` from `geometric MW rank
four`.  Exact component orientations and Shioda corrections give the displayed
three-section height Gram

```text
[ 2/3   -1/6   1/6 ]
[-1/6  13/12   1/4 ]
[ 1/6    1/4  7/12 ],       det = 23/72.
```

Replacing `P1` by its exact half divides the determinant by four, giving a
certified rational rank-three sublattice of determinant `23/288`; together
with the reducible fibers it defines a rank-19 sublattice of discriminant
`46`.  This visible basis is therefore neither the full geometric MW lattice
nor a determinant-948 transported generic basis.  See
[`../../elkies-k3/scripts/derive_mw3_from_mw2_chord.sage`](../../elkies-k3/scripts/derive_mw3_from_mw2_chord.sage)
and
[`../../elkies-k3/scripts/verify_mw3_d5a3a2a2a1a1_twist_qq.sage`](../../elkies-k3/scripts/verify_mw3_d5a3a2a2a1a1_twist_qq.sage).

This rational lower-q anchor is an exact Picard-rank-20 boundary diagnostic,
not the preferred generic construction frontier.  Its additional section is
Noether--Lefschetz special, and the corresponding small-neighbor classes do
not automatically deform in the determinant-948 rank-19 family.  The leading
route is instead the noncollapsed CM-24 q80 chart; a rootless rank-17 or
rank-18 Weierstrass model has not yet been constructed.

Curve 273 is not a direct rational fiber of this particular discriminant-43
anchor.  Clearing denominators in `j_anchor(t)=j_273` gives `t^6` times an
irreducible polynomial of degree ten; `t=0` is the singular `I0*` fiber, so
there are no smooth rational solutions.  This exact exclusion is checked by
[`verify_curve273_not_disc43_anchor_fiber.sage`](../../archive/elliptic-curves/cas/verify_curve273_not_disc43_anchor_fiber.sage).
It does **not** exclude curve 273 from the determinant-948 family containing
the anchor as a special K3: that question has an additional family parameter
which has not yet been reconstructed.

See [`../../elkies-k3/RECONSTRUCTION_PROGRESS.md`](../../elkies-k3/RECONSTRUCTION_PROGRESS.md)
and [`../../elkies-k3/E6_P2_REDUCTION_2026-08-20.md`](../../elkies-k3/E6_P2_REDUCTION_2026-08-20.md).
Recovering this family remains more valuable than repeating large blind
specialization searches: it would expose the actual parameter geometry and
make searches beyond rank 30 reproducible.

A secondary equation-level shortcut is now excluded.  The exact CM embedding
gives an abstract inherited root frame `E8+A2^3`, but promoting all three
`A2` factors to `IV` gives a constant-`j=0` family.  Its order-three CM
automorphism forces even geometric Mordell--Weil rank, contradicting the
target Shioda--Tate rank `19-2-14=3`.  Thus at least one `A2` must be an `I3`
fiber.  See
[`../../elkies-k3/E8_A2_KODAIRA_CORRECTION.md`](../../elkies-k3/E8_A2_KODAIRA_CORRECTION.md).
Any use of this secondary CM chart must therefore be non-isotrivial and must
first be matched to one of the rational Kumar polarizations before neighbor
transport.
One such ambient family is now exact: with `D=t(t-1)`, take
`A=-3r^2D^2` and `B=D^2((t-lambda)^3-2r^3D)`.  Its fibers are
`II*+2 IV+I3+3 I1`, its `j` is nonconstant, and `(r,lambda)=(0,0)` is the CM
endpoint.  This does not yet prove that its one-dimensional rank-three locus
is `X(6,79)`; see
[`../../elkies-k3/E8_A2_MIXED_FAMILY.md`](../../elkies-k3/E8_A2_MIXED_FAMILY.md).

## Corrected modular reconstruction frontier

Two semistable rank-3 neighbors have now been tested over `GF(31)`, with all
promotion conditions kept separate from the closed polynomial equations.

For the `E6/MW3` neighbor, exact frame-glue recovery now proves that the
canonical component-label orbit used by the search is the correct one; the
second height-compatible orbit is excluded.  The reconstruction error was
instead a missing off-diagonal height gate.  For the certified P1/P2 profiles,
the target pairing `-5/6` is equivalent to `(P1+P2).O=1`.

A compiled exhaustive scan of the declared rational chart over
`GF(5),GF(7),GF(11),GF(13),GF(17)` tested 406,655,040 exact cores.  It found 69
P2 sections passing the earlier fiber/component gates and none passing the
target height gate.  One characteristic-11 surface even has three independent
sections, but its height determinant is `13/48`, not `79/16`.  These are nearby
wrong K3 surfaces, not modular seeds for the recovered rank-17 frame.  See
[`../../elkies-k3/E6_MW3_ATTACK.md`](../../elkies-k3/E6_MW3_ATTACK.md) and
[`../../elkies-k3/scripts/recover_e6_mw3_component_glue.sage`](../../elkies-k3/scripts/recover_e6_mw3_component_glue.sage).
Because this split chart was inferred from the Kodaira configuration rather
than derived from the recovered neighbor chain, the bounded emptiness result
does not reject the E6 frame.  It diagnoses the need to backtrack the genuine
geometric construction before any further modular search.

For the `A10/MW3` neighbor, two exact P1+P2 surfaces survive.  A fast
meet-in-the-middle search exhausts the canonical P3 system on both.  One has
no raw solutions; the other has 58 numerator solutions, all rejected because
the alleged pole cancels.  Those cancellations reduce to one polynomial
section `R` (up to sign), and an exact finite-quotient certificate proves
`P1,P2,R` independent.  This gives a genuine modular rank-at-least-three
surface, but `R.O=0` rather than the target `P3.O=1`, so it lies on the wrong
height-lattice component.  See
[`../../elkies-k3/MW3_A10_REDUCTION_2026-08-20.md`](../../elkies-k3/MW3_A10_REDUCTION_2026-08-20.md).

These failures are structurally useful.  They identify the required search
gates: exact Kodaira multiplicity, residual squarefreeness, denominator
noncancellation, deeper component labels, full height Gram, and only then a
Jacobian/lift.  No current modular point passes all of them, so the explicit
`X(6,79)` family is still unrecovered.

## Structural comparison with the rank-28 and rank-29 curves

Exact invariant calculations give the following small-prime valuations in the
displayed integral discriminants.  These are fingerprints, not family
certificates; minimalization and the distinction between discriminant and
conductor remain essential.

| curve | selected discriminant valuations |
|---|---|
| rank 28 | `2^15 3^6 5^6 7^4 11^2 13^4 17^5` |
| rank 29 | `2^19 3^7 5^7 7^4 11^5 13^3 17^4 31^3 41^2` |
| curve 273 | `2^16 3^12 5^8 7^5 13^5 31^2 41^2 47^4 53^3 67^3` |

Thus all three share the bad-prime pattern `2,3,5,7,13`, and the last two also
share `31,41`.  For the rank-29 and rank-30 displayed models,

```text
gcd(|Delta_29|, |Delta_30|)
  = 95418385098324986880000000
  = 2^16 * 3^7 * 5^7 * 7^4 * 13^3 * 31^2 * 41^2.
```

On the other hand, the three `j`-invariants are pairwise different.  For the
rank-29 and rank-30 models, `gcd(c4_29,c4_30)=1` and
`gcd(c6_29,c6_30)=1`.  In particular curve 273 is not merely the same curve in
a scaled Weierstrass presentation, and the invariant comparison gives no
simple twist explanation.

Interpretation: the shared discriminant support is consistent with a common
CRT-shaped search lineage, but it does not establish membership in `X(6,79)`.

An independent exact fixed-root recognition screen now closes one simpler
alternative.  The exact `j`-equation for curve 273 has no rational-square
parameter in any of the 2,329 normalized nonsingular six-root Mestre families
of diameter at most 300, nor in the larger Fermigier control tuple.  A scan of
the generated artifact tree likewise finds no previously unrecognized model
with the same exact `j`.  Since curve 273 has trivial rational torsion, it is
also not a direct specialization of the implemented Elkies--Klagsbrun model,
which has the rational 2-torsion point `(0,0)`.  These are bounded/direct-model
exclusions only: larger or generalized Mestre tuples, isogenous K3 models and
the proposed rootless rank-17 descendant remain open.  The replay is
[`../cas/analyze_icarm_construction_fingerprints.py`](../cas/analyze_icarm_construction_fingerprints.py),
with pinned output
[`../../artifacts/generated-results/elliptic-curves/icarm_construction_fingerprints_v1.json`](../../artifacts/generated-results/elliptic-curves/icarm_construction_fingerprints_v1.json).

A calibrated bounded run of
[`../../elkies-k3/scripts/search_rank17_embedding_graph_v2.py`](../../elkies-k3/scripts/search_rank17_embedding_graph_v2.py)
used 12,000 ambient short-vector lines and the first eight eligible norm shells.
Both the rank-29 control and curve 273 reached only partial graph depth `2`.
The result is therefore non-discriminating and does not justify a larger
heuristic search.  The local checkpoints are
`artifacts/local/elkies-k3/rank17-E29-control-bounded-20260820-best-partial.txt`
and `artifacts/local/elkies-k3/rank17-E30-bounded-20260820-best-partial.txt`.
The common options were

```text
--limit 12000 --max-shells 8 --max-shell-lines 1000
--node-limit 100000 --seconds-per-shell 4
```

## Exact rank and conditional upper bound

The repository certificate proves unconditionally

```text
rank E(Q) >= 30.
```

It does not prove exact rank 30.  The ICARM discussion reports an analytic
upper bound of 31 under GRH and uses root number `+1` with BSD parity to obtain
conditional exact rank 30.  The hypotheses are indispensable; the repository
therefore retains only the unconditional lower-bound claim.

## Residual 2-descent: exact bounded progress

The relation search began with a six-large-prime target and has now produced
the exact support-size chain

```text
6 -> 4 -> 2 -> 3 -> 2 -> 4 -> 2 -> 4 -> 1.
```

The successful stage-two target was

```text
505724623:356162826
84664160213:21346805921
541738517197:261717997519
28691731813798755604363789:17957201189903465826327159
```

Forcing the first three ideals with `--top 2000` gave 399 exactly factored
candidates (the remaining shortlisted candidates exceeded the 160-bit
factorization guard), 268 improvements, and the two-ideal residual

```text
28691731813798755604363789:17957201189903465826327159
7159638381133483906634203654283170391:12780381281373253031851035100853459
```

The new relation is represented by

```text
m = 14332057143548341066300258343667194241.
```

Raw log:
`artifacts/local/elliptic-curves/crt-cycle-stage2.log`.

Two bounded continuations explain why the same one-dimensional search should
not simply be enlarged:

| forced ideal | shortlist | exactly tested | improvement |
|---|---:|---:|---:|
| old 85-bit ideal | 2000 | 2000 | none |
| new 123-bit ideal | 2000 | 6 | none |

The 123-bit modulus is already larger than the real-root scale.  Almost all
points in the integral arithmetic progression then have residual cofactors
above the factorization guard.  Raising the guard would mostly buy harder
factorizations, not better geometry.

The clean continuation is to leave the restricted elements `m-theta` and work
in the full cubic field.  Construct the product of the two target prime ideals,
reduce its three-dimensional Minkowski lattice (with the bad-prime `S` support
handled explicitly), enumerate short elements, and factor their exact norms.
This searches all three coefficient directions and aligns the enumeration with
algebraic norm.  A naive unweighted two-coordinate `a-b*theta` lattice was
tested and rejected: the large curve coefficients make its norm residuals
worse, not better.

### Full-ideal descent implemented

That continuation is now implemented in
[`../cas/search_curve273_ideal_lattice_relations.sage`](../cas/search_curve273_ideal_lattice_relations.sage).
For declared degree-one ideals `P_i=(q_i,theta-r_i)`, it forms the exact product
`I=prod P_i`, constructs its ternary norm form, and LLL-reduces several
determinant-one twists of the `(1,1)` Minkowski embedding.  It then enumerates
the resulting three-dimensional coefficient boxes.  Candidates are ranked
after a small-prime primorial presieve rather than by raw norm; exact rational
factorization and degree-one root labels are used for every reported relation.

This is a bounded relation collector, not a class-group or Selmer completeness
proof.  It materially outperforms the one-dimensional slice once the target
modulus is too large.  Starting with the `85+123`-bit two-ideal residual above,
the first bounded pass (`radius=18`, five embedding shapes, 248,747 distinct
elements) found an exact principal relation replacing it by

```text
8108401645961:3241943091229
13757921887007:10290585023712
1998874339580503775477:182072033289697848802
```

of approximate sizes `43,44,71` bits.  Forcing the first two in the existing
CRT sieve uses only an `87`-bit modulus; 284 of 1,000 candidates improved the
three-ideal target and gave the two-ideal residual

```text
1998874339580503775477:182072033289697848802
210549181078644643738293358684427:12055357087325822035141506894195
```

of sizes `71+108` bits.  A primorial-presieved full-ideal pass on this pair
then gave the four-ideal residual

```text
6757889:3837804
9828720251573:3497176632472
48922801174561:47639959748846
50777509904197:8431789759716
```

of sizes `23,44,46,46` bits.  Forcing the first three uses a `112`-bit modulus;
38 of the 56 exactly factored shortlisted candidates improved the target and
gave

```text
50777509904197:8431789759716
28667248277432793199773894520960218073:12725678201649134529227295218218488637
```

of sizes `46+125` bits.  One more presieved ideal pass returned the four-ideal
residual

```text
110693621:3604894
118153954113061:26623425927990
222635866003829:168812229950784
241667315076287:75987985552523
```

of sizes about `27,47,48,48`.  The first report treated this as a local
plateau because the final candidate ranking minimized the largest residual
prime before residual support count.  Correcting that objective exposed a
strictly better relation already present in the same `radius=18` bounded pass:

```text
coordinates = (1876028045451,-1977200499447,-8769733)
remaining =
110082763798456567967312309869978698968187793:12357613814041997795032595200900142542684867
```

The remaining degree-one ideal is `147` bits.  A first report from a recursive
`radius=18` pass again preferred support count and therefore displayed a larger
one-ideal continuation of `163` bits.  Ranking nonclosures by maximum residual
prime and then total residual size instead exposed the much better five-ideal
continuation already present in the same bounded pass:

```text
3069949:1911346
63217513:22678367
168510611:8392763
7754996948873:5321926736053
38169724803043:26317327110158
```

Its prime sizes are only `22,26,28,43,46` bits.  Two further exact full-ideal
steps gave

```text
19767179:6801961
2024658179:1400637434
12736098241512293:9617972020917348
231573912732876911:136784656449821585
```

of sizes `25,31,54,58`, and then

```text
7636709:6045380
56740549:8493651
104379329:22584173
413582930291:366762253550
555442705507:488627836116
```

of sizes `23,26,27,39,40`.  The next identical bounded pass did not improve
the `40`-bit ceiling (its best new ceiling was `55` bits), so the recursive
walk stops here rather than silently enlarging its radius.  A separately
certified side branch `147 -> (49,74) -> (25,96)` is retained in the data but
is arithmetically inferior.  A 2,000-candidate CRT pass forcing its 25-bit
ideal found no factor-base-smooth quotient and no replacement below 96 bits.
Five separate `radius=18` full-ideal searches, each forcing just one of the
current endpoint ideals, also found no exact closure.  Their best genuine new
residual ceilings, in endpoint order, were `64,52,56,66,64` bits.  The
apparently smooth rational basis vector is not a closure: exact ideal
factorization exposes an undeclared conjugate prime ideal above the same
rational `q`.  The search now reports that obstruction explicitly and removes
nonprimitive scalar multiples before preselection, so they no longer consume
the exact-factorization budget.

All nine main-chain relations, the two side-branch relations, and every
support transition are independently pinned by
[`../cas/verify_curve273_full_ideal_descent_chain.sage`](../cas/verify_curve273_full_ideal_descent_chain.sage).
It reconstructs each arbitrary ideal element from its integral ideal basis,
checks valuation exactly one at every declared degree-one prime ideal, and
proves that the remaining norm is supported on the ordinary factor base plus
`S`.  Its certified support-count chain is

```text
2 -> 3 -> 2 -> 4 -> 2 -> 4 -> 1 -> 5 -> 4 -> 5.
```

The principal replay command and raw bounded continuation logs are:

```bash
sage -python elliptic-curves/cas/search_curve273_ideal_lattice_relations.sage \
  --target 28691731813798755604363789:17957201189903465826327159 \
  --target 7159638381133483906634203654283170391:12780381281373253031851035100853459 \
  --radius 18 --shape-shifts=-48,-24,0,24,48 \
  --preselect 2000 --factor-top 400 --factor-base-bound 1000000

sage -python elliptic-curves/cas/verify_curve273_full_ideal_descent_chain.sage

sage -python elliptic-curves/cas/analyze_curve273_relation_pool.py \
  --include-full-ideal-chain --include-crt-cycle-logs

artifacts/local/elliptic-curves/crt-cycle-fullideal-stage1.log
artifacts/local/elliptic-curves/ideal-lattice-stage2-presieved.log
artifacts/local/elliptic-curves/crt-cycle-fullideal-stage3.log
artifacts/local/elliptic-curves/ideal-lattice-stage3-presieved.log
artifacts/local/elliptic-curves/ideal-lattice-stage4-one-lp.log
```

The exact pooled audit now includes 115 earlier special-`q`/multi-`q` rows,
eleven certified chain/side-branch rows, and 228 distinct CRT-cycle rows after
removing three duplicates.  On 676 large-prime ideal columns, all 354 rows are
independent: the large-prime rank is 354 and the nullity is zero.  Hence none
of the logged candidates hides a cycle.  Solving the same sparse system against
the final five-ideal endpoint succeeds with a thirteen-row witness and 84
factor-base ideals.  This independently certifies that the original relation
component reaches the `23,26,27,39,40`-bit endpoint, while making clear that a
new independent incidence at one of those five endpoint ideals is still
required for closure.  The `IN_SPAN` endpoint diagnostic is deliberately only
a reachability statement: the original six-ideal target was itself the XOR of
pool rows 0 and 5, so telescoping that seed pair with the certified chain must
put the endpoint in the same row span.  A second representation would instead
produce an LP-column dependency; the computed nullity zero rules that out for
the retained pool.

### Bounded conjugate/twisted endpoint extension

A subsequent bounded pass replaced the selected ideals above each of the two
split endpoint primes by the product of their conjugates and applied all
individual factor-base twists above primes at most `50`.  Neither branch
closed; their best new residual ceilings were `53` and `65` bits.  Pairwise
factor-base twists of the complete five-ideal endpoint found exact branches

```text
5 ideals -> 4 ideals, maximum size 46 bits
5 ideals -> 3 ideals, maximum size 48 bits.
```

Standard `radius=18` continuations returned the first branch to the known
40-bit endpoint and worsened the three-ideal branch to a 56-bit ceiling.  These
are bounded search outcomes, not class-group or Selmer completeness results.

`analyze_curve273_relation_pool.py --include-ideal-lattice-logs` reconstructs
each retained `R30IDEAL|best` element from its exact, possibly
factor-base-twisted ideal basis, checks the declared prime-ideal valuations,
and includes it in the mod-two sparse audit.  Combining the pinned chain, all
retained CRT-cycle rows, and the new bounded logs gives

```text
exact rows = 444
large-prime ideal columns = 962
large-prime rank = 444
nullity = 0
```

Thus none of the 90 new nonduplicate ideal-lattice rows closes a hidden cycle
with the previous pool.  The raw searches and combined audit are local
checkpoints in `artifacts/local/elliptic-curves/` with date suffix `20260821`.
Each LP-free dependency now also records the product of its actual principal
generators in power-basis coordinates.  The portable quotient/certification
layer is documented in
[`BNF_FREE_RESIDUAL_2SELMER.md`](BNF_FREE_RESIDUAL_2SELMER.md); it makes no
class-group, Selmer, or rank claim from this bounded pool.
The 30 known Kummer generators now also have a reproducible 59-local/
54-witness signature map from
[`../cas/analyze_curve273_kummer_fingerprint.py`](../cas/analyze_curve273_kummer_fingerprint.py),
so any candidate global squareclass can be reduced against the known image
before more relation collection is attempted.

## Efficient next gates

1. **Execute the next equation-level neighbour.** The first q6 and corrected
   q8 pencils are exact over `QQ`, reaching the parent `D13/MW4` equation.
   Execute the selected orbit-85 q24 divisor on that parent to construct the
   immediate equation child `D12/MW5`, with the same binary-quartic multiplier
   and fully cleared-residue regressions.  This local D12 child is not the
   recovered R17 route endpoint.
2. **Execute the remaining certified degree-two chain.** Realize the ten
   later q6/q8/q4-series/final-q6 neighbours in characteristic zero. The
   complete D13-to-rootless lattice transport and chamber/nef certification
   are already exact; only the equation-level execution remains.
3. **Retain the exact specialization exclusion.**  The published rootless
   family is explicit, and its degree-24 `j`-recognition equation for curve 273
   is irreducible over `QQ`.  Do not spend further work looking for a rational
   parameter in this chart.  Test only genuinely different proposed families
   or isogeny constructions.
4. **Keep exact-rank descent independent.**  Closing the remaining
   `23,26,27,39,40`-bit ideal support can prove an upper bound or expose a
   new point, but it does not reconstruct the family.

## Historical and parallel gates

1. **Family certificate first.**  Distinguish and algebraize the two exact
   formal branches with rational tangents `8/87` and `1/12` through the
   noncollapsed CM-24 `q=80,(4,20)` chart, which retains all three generic MW
   directions.  Both are certified over the CM compositum through order 18;
   their GF(7) reductions lift through `h^20`.  The incomplete degree-ten
   coordinatewise Pade test is only a bounded negative recognition result.
   The `8/87` branch has an exact two-dimensional cubic space matching the
   order-145 slope-`5` space modulo seven, and the
   modular filtered quotient dimensions `1,5,15,33,48,63` fit `15*n-12` from
   degree three. A modular `(D,Q)` plane projection subsequently gives
   normalization genus zero, so the genus-two Shimura curve must be sought on
   the marked-section cover rather than in the unmarked coefficient curve.
   The first marked section gives the exact mod-seven cover
   `u^2=(t+3)(t+4)(t+6)(t^2+6t+4)`. It has genus two, but its absolute Igusa
   invariants `(4,2,5)` differ from `(4,4,3)` for the known source model
   modulo seven. The second marked section needs a distinct genus-three cover;
   together they give a genus-six biquadratic marked curve with genus-one
   third quotient. This two-short-section descent is therefore a high-genus
   detour, not the desired source marking. It does not reject the unmarked
   slope-5 surface branch, whose true level-79 marking may involve `P3`. Its
   ideal does pass five withheld orders
   through order 90, but its P1 branch set has trivial projective stabilizer
   over `GF(49)` and its genus-one third quotient has `j=0`, not the known
   quotient values `1,3`. The other slope has no relation through degree five,
   but its generic split-prime degree-six space has fifteen generators. At
   both `p=73` and `p=79`, the irreducible degree-32 `(P,D)` projection has
   total delta 464 and normalization genus one, the expected elliptic quotient
   shape. Point counts nevertheless give traces `15,22` at `p=73,127`, which
   match neither known source elliptic factor nor their twists. Conditional on
   the fifteen-sextic ideal being the global branch, slope `1/12` is therefore
   rejected. Return to the rational unmarked slope-`8/87` surface and recover
   the true `P3/Q79` quadratic marking; the failed P1/P2 covers do not exclude
   that different genus-two cover. A bounded order-230 recovery of the
   selected `P3` pole and numerator instead puts both on the same genus-one
   `j=0` field as the third P1/P2 quotient,
   `w^2=t^4+3t^3+3t^2+1`. Thus this particular three-section orientation also
   fails to expose the known genus-two source; try a different CM orientation
   or descend `Q79` directly. This is not an unmarked-surface rejection.
   Its homogeneous normalization now gives an explicit modular parameter,
   with surface-function degrees `5/4,10/8,8/6,15/12`; all fourteen modular
   ideal generators vanish identically. The 26-term modular plane support has
   no exact order-30 relation on that same support, so lift this finite rational
   ansatz (or enlarge the characteristic-zero plane support), then recover the
   marked cover and verify it in
   characteristic zero. The exact q80-to-rootless lattice continuation is
   now the short chain `q=4,4,12,12,4,6`, with MW ranks
   `3,4,5,6,13,16,17`; execute those pencils after the family is algebraized.
   The first q4 class is already reduced to a degree-two class with zero MW
   projection in `L(2O+4F)`, and an exact height-shell argument excludes
   every section wall. Root-lattice primitivity also excludes the only
   possible negative bisection, so this first class is fully nef. The exact
   `D5+E6` local resolution now gives `L(D)=<T^2,x-T>` and
   `U=(x-T)/T^2`. The ambient child has `D9+A3`; the marked rank-19 collision
   gives the required `D9+A4`, while CM24 enhances to `D9+A5+2A1`. The second
   q4 class is also fully nef, by an exact rank-four MW height-shell
   calculation and the same primitive-root bisection argument. Its exact
   coordinate on the first child is `W=(X-3v^3-x1v-x0)/v^2`. At CM24 its
   finite special fiber has valuations `(3,4,8)`, hence type `IV*`, and the
   full signature is `D7+E6+3A1/MW2`, exactly matching the transported pinned
   frame. Lift the third pencil below to characteristic zero, then execute the
   remaining three downstream pencils. The third q12 divisor reduces
   to `S+2O+2F+root_correction`, where `S` is an integral height-eight section
   with `S.O=2`. Thus its old-degree-three pencil lies in `L(S+2O+2F)` and
   needs marked-section data rather than another zero-MW compensated
   coordinate. At CM24 the exact specialized marking is the polynomial
   height-three section `Q_CM=P1+3P2` over `QQ(sqrt(-6))`, reducing the
   generic-fiber calculation to
   `<1,X,(Y+y(Q_CM))/(X-x(Q_CM))>` plus vertical gates. Exact transport makes
   the CM target `2A6+3A1/MW3`. The unique effective lift is
   `D=Q_CM+2O+4F+R`, with old-component coefficients
   `A1:(0),(1),(0)`, `E6:(1,2,3,2,1,2)`, and
   `D7:(2,4,3,3,5,6,3)`. The two pure chord-translation completions
   have the wrong roots (`D8+E7+A1` and `D5+E6+A3+A1`), proving that a
   nonzero `X` coefficient is required. A bounded scan of all 2,401 naive
   linear-X/quadratic-translation tuples over `GF(7)` has branch degree 17 or
   23 whenever the X coefficient is nonzero, so the next ansatz must be the
   nine-dimensional ambient `a(W)+bX+c(W)z_Q`. Its two selected-I2 rows and
   one E6 row combine with four exact D7 rows from the complete local ideal
   `(Y,U^2,ZU,Z^3)`. The rank-seven kernel gives `Vnew=N1/N0`; clearing the
   chord leaves `X-Qx` times a W-degree-nine cubic. At split `p=73`, generic
   values `Vnew=1` and `Vnew=7` are irreducible and normalize to genus one.
   The weighted infinity form is a scalar multiple of
   `(xi-3)^2*(xi+6)`, independent of `Vnew`. At `Vnew=7`, the rational
   simple branch `xi=-6` gives the bounded certificate
   `y^2+12xy+27y=x^3+51x^2+40x+26` over `GF(73)`.
   Exact interpolation from 49 canonical-origin fibers, followed by seven
   withheld checks, now gives the full CM24 Jacobian over `GF(73)(V)`. Its
   discriminant has multiplicities `7,7,2,2,2,1,1,1,1`, hence fibers
   `2I7+3I2+4I1` and the required `2A6+3A1/MW3` marking. The reconstructed
   `A,B` have degrees `8,12`; see
   [`../../elkies-k3/scripts/reconstruct_q80_third_q12_jacobian_gf73.sage`](../../elkies-k3/scripts/reconstruct_q80_third_q12_jacobian_gf73.sage).
   The earlier unresolved-cusp rows gave genera four and three and are
   rejected. The CM24 third q12 pencil and its mod-73 Jacobian marking are
   therefore exact; lifting the compact model to characteristic zero is the
   next gate. Do not impose both quadratic roots; that
   is the ineffective lift. The
   alternate q8 move is genuinely
   noncollapsed (`E6+A7/MW4`, CM24 `E7+A7+2A1/MW2`) but loses a generic MW
   direction, so it is only a possible coordinate normalization rather than
   a replacement for q12.
   Downstream, the fifth q4 step now has two distinct equation strategies.
   The pinned `A1/MW16` child still requires the saturated local module for
   `L(O+(-R))`. A bounded alternate q4 class is easier over `GF(73)`: after
   normalizing every independent fiber by the marked section pair `(1,4)`,
   all twelve coefficients reconstruct from fourteen samples, pass four
   withheld samples, and give a degree-four branch squareclass, hence a
   genus-one fifth-child gate. The artifact is
   [`../../artifacts/generated-results/q80-fifth-q4-marked-projection-pair14-gf73.json`](../../artifacts/generated-results/q80-fifth-q4-marked-projection-pair14-gf73.json),
   SHA256
   `e46c9925c6870a6f9185f36994a5aef682382bba7a9bf8d2adc3d897420988fa`.
   Restoring the omitted factorization unit gives polynomial degrees `(8,12)`
   and exact fibers `I6+I6+I5+I2+5I1`, hence CM24 root rank 15/MW3; the
   twisted-Jacobian artifact SHA256 is
   `78c654d35acccb907a3b019bc309c84d7d7b705d8d6a521e17f3f169fad67ca9`.
   The pair-14 equation is now matched exactly to a CM24 lattice class. The
   oriented norm-eight affine shell for `(a,b)=(2,2)` and horizontal `(1,0)`
   has an exact first retained class with root data `(15,82,360)`:
   `(2,2,-98,42,47,21,8,66,-179,29,81,219,30,65,48,23,32,53,-89,-8)`.
   It is already reduced, has old-fiber degree 49, and gives
   `2A5+A4+A1/MW3`. The exact marking artifact SHA256 is
   `c3949e37638b138184bf9591e127aa19c93d8d927a4b224b865df7c20dcf6cac`.
   The bounded run stopped at this first full hit after 193 primitive
   presentations; it is not a uniqueness certificate over all 797,472 signed
   affine classes.
   This is a CM24 boundary class, not a generic determinant-948 divisor.
   In the three retained productive generic windows, all 48 `A1/MW16`
   candidates specialize instead to `(16,66,2048)` =
   `D4+3A3+3A1/MW2`; they form CM nef classes of old-fiber degrees 47
   (32 candidates) and 43 (16 candidates). Thus the pair-14 equation cannot
   yet be attached to the generic rootless suffix. Saturating the CM24 root
   lattice resolves the order-three torsion ambiguity: the degree-43 class is
   the `(0,1)` lift, while the productive degree-47 class is the `(2,3)` lift.
   The latter is exactly
   `O+P_(2,3)+2*Theta0+Theta1+Theta2` on one `I6` fiber, with no residual
   whole fibers. Thus its equation problem is a single-I6 compensated
   projection; the previously computed genus-three **uncompensated** pair-23
   cover is not the target divisor. The exact decomposition is certified by
   [`../../elkies-k3/scripts/analyze_q80_deforming_fifth_vertical_compensation.sage`](../../elkies-k3/scripts/analyze_q80_deforming_fifth_vertical_compensation.sage).
   Pair `(0,1)` gives the
   distinct CM24 fibration `D4+A5+2A3/MW3`; it too has an exact first CM24
   lattice hit, with old-fiber degree 47, but bounded q6 sampling improves
   only to root rank 14/MW4.
   Pair14's wider bounded low-q and p-neighbor continuation likewise bottoms
   out at root rank 14/MW4. Separately, the selected generic lattice node was continued by a q6
   witness with `(a,b)=(2,3)` to a different rootless/MW17 frame, and the full
   q80-to-selected-q4-to-q6 NS composite has determinant one. Its artifact SHA256 is
   `48381d91e288b2cefb85b1484d351d659748f801ea57d190453bd2db0a56eaab`.
   The selected fifth equation is now exact over `GF(73)`: the pair-23 secant
   gauge has fibers `I0*+3I4+3I2` and root data `(16,66,2048)`, matching the
   productive degree-47 lattice class. Its generated artifact SHA256 is
   `23fc49bce2618a6d3c5f5e18ded34b4ffbee220be83523ae250bf7774a91db14`.
   The raw final class has degree three, but one zero and one `A1` reflection
   reduce it to old degree two, `D.O=1`, MW norm `23/2`; hence the final
   equation will be another degree-two marked problem once the correct fifth
   vertical class is reconstructed. This reduced q6 class is now fully nef:
   its exact section CVP has minimum pairing one, and parity excludes a
   negative bisection. Its horizontal section is explicitly
   `S=(5,1,-1,-2,4,2,-1,2,1,-1,1,0,1,-1,1,0,0,0,0)`, with `S.O=4`, and
   `D=O+S-F` generically. At CM24 the section is identity at all three `I2`
   fibers, hits endpoints of two `I4` fibers, and is nonidentity at `I0*`.
   In deterministic special-fifth roots its vertical class is
   `D-O-S=2F-(R1+R5+R8+R11+R13+R16)`. Locally the nonzero cycles are
   `F-R1`, `F-R5-R8`, and `F-R11-R13-R16`; three affine copies are required
   but only two global fibers are available. Hence the constant is not in the
   specialized local module, explaining the degree-12 raw chord failure.
   The corrected q6 section artifact (no I2 hits, I0* correction one) has
   SHA256 `fcd61f89daab0a68785a006e6b10dc3829b1c30c24243b67a4e1b80c7d6e6e09`.
   The CM24 q6 child has `A1+2A3+2A4`, root data `(15,66,800)`, and MW rank
   three; generic rootlessness/MW17 remains exact. The saturated three-fiber
   transform is now exact over `GF(73)`.  The cleared chord module has Smith
   diagonal `(1,s-27)`, its compensated generator is
   `q_sat=(q0+63)/(s-27)`, and in the basis
   `(1,s,s^2,q_sat,s*q_sat)` the local gates are
   `I2:(72,6)`, `I4:(64,3)`, and `I0*:a2=0`.  The kernel rows are
   `(1,0,0,41,48)` and `(0,1,0,6,72)`.  Its unit-preserving Jacobian has
   fibers `2I5+2I4+I2+4I1`, exactly the expected root data `(15,66,800)`.
   The artifact SHA256 is
   `7d3866855b9995b733193de9c5d5e3ba1cea6aa1292141e06d0d5011f28975e3`;
   the artifact includes the explicit `q0(s,R)` and globally minimal
   polynomial `A`, `B`, and discriminant coefficients.
   The live downstream gate is now characteristic-zero lifting, not another
   finite-field local-module search.

   This construction remains a downstream `H2=diag(4,237/2)` chain.  The
   exact source branch has `H3=[[21/2,3],[3,46]]`; the doubled forms have the
   same determinant `1896` but minima eight and 21 and are not integrally
   isometric.  The recurrent q=8 signal is therefore the height-four
   Humbert-8 entrance or a special Noether--Lefschetz return/collapse, not the
   local Smith factor in the final q6 module and not a direct `H3` fiber.
   Compare this with a compact
   `q=60,(4,15)` chart, which retains two directions; keep the exact
   `q=60,(5,12)` equation as an anchor and marking check only.  Do not resume
   CM-43 Riemann--Roch construction: its marked q8 and q60 divisors both
   collapse to the old fiber after their CM-only fixed sections are removed.
   Only after obtaining a rootless characteristic-zero model should curve 273
   be tested by exact specialization and `Q`-isomorphism.
2. **Exact specialization certificate second.**  Solve for a rational family
   parameter and a `Q`-isomorphism to curve 273, then substitute that parameter
   into every generic section and match the resulting points.  Canonical
   heights are not preserved isometrically under specialization, so an
   isometric sublattice search in the numerical height Gram cannot supply this
   certificate.  Height growth and short-vector graphs may prioritize
   candidates, but literal section coordinates must finish the check.
3. **Close the certified small-ideal residual.**  The global sparse audit is now
   complete for every retained special-`q`, multi-`q`, CRT-cycle, and chain
   row, and has nullity zero.  Target a second incidence at one of the final
   `23,26,27,39,40`-bit degree-one ideals with a marked CRT batch, an improved
   full-ideal presieve, or a class-group/S-unit computation.  For recursive
   descent rank nonclosures by maximum residual-prime size, then total residual
   size, and only then support count; the opposite order hid every productive
   continuation of the 147-bit ideal.  Keep every factorization bound and
   enumeration radius explicit so a negative run remains a bounded experiment.
4. **Only then search for improvements.**  Determine the rational generic rank
   of the reconstructed determinant-948 rootless fibration; geometric Picard
   rank at a CM specialization alone is insufficient.  With its parameter and
   rational sections exposed, optimize specializations by exact local
   conditions and finite-quotient escape of new points.  Re-minimize every
   specialization before interpreting discriminant valuations or conductor.

These gates distinguish two separate objectives: reconstructing how the curve
was found, and proving whether its algebraic rank is exactly 30.  Progress on
one does not automatically settle the other.
