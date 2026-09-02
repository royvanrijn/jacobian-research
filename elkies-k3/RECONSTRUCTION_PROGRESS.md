# Elkies rank-17 K3 reconstruction progress

Status checkpoint: 2026-08-21, consolidated after recovery of the canonical
Kumar `E7+E8` entrance.

> **Historical snapshot.** The reconstruction goal below has been completed:
> the selected H3 corridor reaches q12/orbit5867, the rootless Jacobian has
> `24I1`, and its full saturated determinant-948 Mordell--Weil lattice is the
> pinned R17 lattice.  This file remains active only as process evidence for
> the early numerical-collision and route-selection lessons.  Use
> [`README.md`](README.md) for current status and
> [`ELKIES_K3_PROCESS_ATLAS.md`](ELKIES_K3_PROCESS_ATLAS.md) for chronology.

## Historical goal

Recover an explicit elliptic K3 fibration realizing the already-recovered
rank-17 Mordell--Weil lattice, together with explicit generic sections.  Once
that model exists, feed real specializations into the ignition/cascade search in
`RANK_GROWTH_SEARCH.md`.

This document separates three evidence levels:

- **exact**: certified by integer/lattice/algebraic computation in the repo;
- **numerical**: high-accuracy solutions of the reconstruction equations;
- **pending**: not yet exactified or proved.

## Exact starting point

The recovered rank-17 lattice is stored in `data/lattice/rank17_gram.txt` and
`data/lattice/short_vector_basis_gram.txt`.  The current identification work
matches the Elkies `X(6,79)` example: quaternion discriminant `D=6`, level
`M=79`; the corrected local criterion at the odd ramified prime `p=3` is
satisfied.  See `results/IDENTIFICATION-NOTES.md`.

The norm-4 shell is unusually rich:

- 1311 unsigned `+/-` pairs;
- 2622 signed minimal vectors;
- 184242 exact additive identities `a+b=c` among signed minimal vectors;
- a 17-vector seed generates the complete signed shell under these additive
  relations.

These files are the principal reconstruction data:

```text
data/lattice/rank17_gram.txt
data/lattice/short_vector_basis_gram.txt
data/relations/all_2622_signed_short_basis.npy
data/relations/minimal_additive_triples.npy
```

Public Elkies lecture material gives the working model expected for the
reconstruction: after a suitable change of variables one may work with

```text
y^2 = x^3 + A(T) x + B(T)
```

with `deg A <= 8`, `deg B <= 12`; minimal height-4 sections have quartic
`x(T)` and sextic `y(T)`.  The canonical final representation is still to be
stored under `data/k3-model/` once exact coefficients are recovered.

## Exact E6 lattice route and invalid split chart

The direct Coxeter continuation described below remains a useful diagnostic,
but it is no longer the primary reconstruction route. Exact neighbor searches
from the 17-by-17 positive frame have produced lower-Mordell--Weil-rank
fibrations with larger reducible root systems. The most promising node has

```text
ADE = E6 + A3^2 + A1^2
MW rank = 3
reduced MW determinant = 79/16
fibers = IV* + I4 + I4 + I2 + I2 + 4 I1.
```

This frame is certified to lie on the same Neron--Severi lattice in two
independent ways.  First, the two `U`-extended lattices have signature
`(1,18)`, matching local genus, and cyclic discriminant group of order 948, so
Nikulin's indefinite-genus uniqueness theorem applies.  Second, the actual
three-arrow discovery path has been recovered exactly:

```text
rank17 --q=90--> MW7 --q=4--> MW4 --q=4--> E6/MW3.
```

[`scripts/verify_e6_neighbor_chain.sage`](scripts/verify_e6_neighbor_chain.sage)
reconstructs every stored child frame and checks the composite integral
isometry in
[`data/fibrations/e6_ns_transport_from_rank17.txt`](data/fibrations/e6_ns_transport_from_rank17.txt).
See [`E6_NEIGHBOR_CHAIN.md`](E6_NEIGHBOR_CHAIN.md).

This closes the lattice-transport gap, not the geometric one.  The repository
still lacks the explicit rational elliptic parameters and Weierstrass models
obtained by executing those three neighbor operations on a genuine K3 model.

Putting `IV*` at infinity reduces the short Weierstrass degrees to
`deg A <= 5`, `deg B <= 8`. Exact triangular elimination and a rational
parametrization reduce the one-section locus to eight fiber equations in
eleven variables, of expected dimension three.

The first complete point of the closed equation scheme has been reconstructed
exactly over `GF(31)`. Its full section identity and closed fiber equations
vanish, and the reduced `8 x 8` Jacobian is nonsingular.  An exact open-stratum
audit nevertheless gives `ord_1(Delta)=5`, so this is an `I5` boundary point,
not a promoted target-fiber seed.  Its polynomial companion also satisfies
`P1=-2*P3`.

The initial reduced-core scan was finite and concrete, but its P2 promotion
criteria were incomplete. Exact frame-glue recovery confirms the chosen
component orbit; the omitted condition was the off-diagonal Shioda pairing.
For this basis, `<P1,P2>=-5/6` is equivalent to `(P1+P2).O=1`. Complete
rational-chart scans over `GF(5),GF(7),GF(11),GF(13),GF(17)` tested
406,655,040 cores and found 69 component-valid P2 sections, all with the wrong
pairing.  The declared split E6 chart was inferred from the Kodaira
configuration rather than derived by geometric transport along the recovered
neighbor chain.  It must not be searched further.  The E6 route must instead
be rebuilt from that chain with the intersection condition in its defining
equations. See
[`E6_MW3_ATTACK.md`](E6_MW3_ATTACK.md).

## Correct upstream anchor: the H3 Kumar `E7+E8/MW2` frame

The primary-source construction does not begin by deforming a chosen Kodaira
realization of the CM endpoint.  It begins with the principally polarized QM
abelian surface attached to the non-CM point on `X(6,79)`, forms its
Dolgachev--Kumar K3 surface, and uses the canonical fibration

```text
E7+E8, MW rank 2, MW regulator 474.
```

## Reusable exact degree-two neighbor engine

The repeated lattice operations for a degree-two neighbor are now centralized
in [`scripts/exact_neighbor_engine.sage`](scripts/exact_neighbor_engine.sage).
Given an isotropic divisor, old fiber, and explicitly supplied section/
component walls, it performs deterministic fixed-component reduction, exact
primitive `U`-splitting, and root/MW minimization of the child.  Its q80 first
`q=4` regression checker is
[`scripts/verify_exact_neighbor_engine.sage`](scripts/verify_exact_neighbor_engine.sage);
it reproduces the known `D9+A4/MW4` child.  See
[`EXACT_NEIGHBOR_ENGINE.md`](EXACT_NEIGHBOR_ENGINE.md) for the explicit
interface and the deliberate boundary between supplied-wall nonnegativity and
a full global-nef proof.  This is the reusable lattice layer for executing
the certified H3 degree-two chain; it does not itself produce a
characteristic-zero pencil.

An exact finite binary-form/glue classification now finds precisely three
Kumar frames in the recovered determinant-948 genus.  Their height Grams are

```text
[[5/2,1],[1,190]],  [[4,0],[0,237/2]],  [[21/2,3],[3,46]].
```

All three have exactly 366 roots, root lattice `E7+E8`, cyclic discriminant
group `Z/948`, and the same local genus as the recovered rank-17 frame.  The
middle height lattice is uniquely distinguished by an automorphism other than
`-1`.  Exact elliptic quotient labels identify its extra involution as
`w2=w237`, the hyperelliptic involution of Elkies's genus-two model.  This
identifies the canonical involutive H2 comparison frame, but a later exact
Humbert calculation corrects the source marking to the third frame H3.  See
[`KUMAR_E7E8_BACKTRACK.md`](KUMAR_E7E8_BACKTRACK.md) and
[`scripts/deconstruct_x0679_quotients.sage`](scripts/deconstruct_x0679_quotients.sage).

Pulling the primitive Humbert-21 equation back to the H92 chart now isolates
the H3 component exactly over `QQ`.  Its normalization is birational to the
published level-474 curve

```text
y^2 = -27*x^6 + 198*x^4 - 171*x^2 + 576.
```

At the published rational non-CM point, the H21 and H92 unmarked K3 models
are isomorphic over `QQ`.  The remaining source datum is therefore the marked
descent: identify the rational H21 entrance with the signed H3 q6 divisor
`O+(-P1)-F` and determine the field of definition of the height-`21/2`
section.  The exact normalization and non-CM anchor are replayed by
[`scripts/normalize_h21_h92_level474_qq.sage`](scripts/normalize_h21_h92_level474_qq.sage)
and
[`scripts/verify_h3_noncm_q6_source_anchor.sage`](scripts/verify_h3_noncm_q6_source_anchor.sage).

The lattice/chamber continuation from H3 is exact through rootless MW rank
seventeen:

```text
H3 E8+E7/MW2 --q6--> E8+E6/MW3 --q8--> D13/MW4
 --q24--> D12/MW5 --q6--> A11/MW6 --q8--> 2A5/MW7
 --q4--> 3A3/MW8 --q4--> A3+2A2/MW10
 --q4--> 5A1/MW12 --q4--> 4A1/MW13 --q4--> 3A1/MW14
 --q4--> 2A1/MW15 --q4--> A1/MW16 --q6--> rootless/MW17.
```

Every selected arrow is a certified nef degree-two pencil. Equation-level
execution in characteristic zero beyond the exact D13/MW4 child remains open;
the MW10-to-rootless suffix is no longer an open lattice or chamber gate.

For the parallel H2 comparison, there are three exact CM interpolation
anchors.  The
`w3`-fixed points at `t=infinity` transport to the singular K3 with
transcendental lattice `[[2,1],[1,2]]`, determinant `3`.  The free CM orbit at
`t=+/-2` transports to `diag(4,6)`, determinant `24`; its primitive K3 vector
is `(70,86,-3)` with square `-158` and divisibility `79`.  Primitive closure
of the `H2` frame gives `E8+E8+A2` with MW zero at discriminant `-3`, and
`E8+E8` with MW lattice `diag(4,6)` at discriminant `-24`.  The resulting
standard Inose equations are

```text
Delta=-3:   Y^2 = X^3 + T^5*(T-1)^2,
Delta=-24:  Y^2 = X^3 - 51*T^4*X + T^5*(T^2-92*T+1),
```

the second up to rational quadratic twist.  These exact boundary equations
replace an unanchored one-parameter search by a local-deformation and
interpolation problem.  See
[`scripts/classify_kumar_cm_frame_extensions.sage`](scripts/classify_kumar_cm_frame_extensions.sage)
and
[`scripts/verify_kumar_cm_inose_anchors.sage`](scripts/verify_kumar_cm_inose_anchors.sage).

The third anchor is the rational discriminant-43 K3.  Its Gross vector
`-11*i-j-5*k` transports to `(169,167,-128)`, with square `-40764`,
divisibility `948`, and complement `[[22,1],[1,2]]`.  In the Humbert-8 plane
it is located at

```text
r=-1225/722,  s=-93312/442225,
z^2=-43*(11664/6859)^2.
```

Two exact short sections have height Gram
`[[5/2,-1/2],[-1/2,5/2]]`; fourteen good-prime counts give the expected
CM-`43` weight-three coefficients.  This is replayed in
[`scripts/verify_cm43_humbert8_anchor.sage`](scripts/verify_cm43_humbert8_anchor.sage).
The finite Frobenius fingerprint does not replace the still-missing explicit
inverse neighbor between the two characteristic-zero Weierstrass models.

The two Kumar sections should not be reconstructed symmetrically.  The
height-4 generator meets the identity `E7` component and has `P.O=0`.  The
height-`237/2` generator has local correction `3/2`, hence `Q.O=58`.  Up to
fiberwise negation, the hyperelliptic Atkin--Lehner involution fixes the small
generator and negates the large one.  Thus the surface and small section
descend to `t`, while `u` carries the level-79 section descent.  This explains
why a direct two-section Kumar system is unnecessarily large.

For the parallel H2 route, the height-4 summand identifies the ambient moduli
surface exactly.  Since the `E7` root determinant is two, it gives the
Humbert discriminant-8 family;
the height-`237/2` summand alone gives discriminant 237.  Hence this comparison
locus is the intersection `H8 cap H237`, not an unspecified curve in the four
Kumar coefficients.  It is not the corrected H3 source curve.
Elkies--Kumar's discriminant-8 ancillary data give explicit rational functions
for all Kumar and Clebsch--Igusa coefficients
in two parameters `(r,s)`, as well as the oriented double cover `z`.  They are
recorded and verified in
[`KUMAR_E7E8_BACKTRACK.md`](KUMAR_E7E8_BACKTRACK.md) and
[`scripts/verify_humbert8_kumar_entrance.sage`](scripts/verify_humbert8_kumar_entrance.sage).

The H2 ancillary construction also resolves its preceding geometric step.
Before the Kumar equation it uses the two-parameter fibration

```text
Y^2 = X^3 + T*(r+(2*r+1)*T)*X^2
            + 2*r*s*T^4*(T+1)*X + r*s^2*T^7,
```

with fibers `D9+E7+4I1`, and performs the explicit two-neighbor
`U=(X+s*T^3)/T^4`.  The resulting pointed quartic has rational point
`(T,V)=(0,s)` and Jacobian exactly equal to the published Kumar model.  This
symbolic characteristic-zero replay is
[`scripts/verify_humbert8_d9e7_two_neighbor.sage`](scripts/verify_humbert8_d9e7_two_neighbor.sage).
It shows that q=8 is structural: it is the determinant of the pre-neighbor
root frame `D9+E7`, not evidence that every observed q=8 lattice move is the
same neighbor.  The generic determinant-948 specialization adds one
torsion-free MW direction of height `237/2` to this frame.

The second quartic point also supplies the generic height-four Kumar section
without a search.  At the CM-`43` anchor, write the two exact height-`5/2`
sections as `P1,P2` and this height-four section as `P3`.  Then the generic
level-79 direction is exactly

```text
Q79 = 4*P1 - 5*P2 + P3.
```

Its height is `237/2`, it is orthogonal to `P3`, and exact group law gives
denominators `h^2,h^3` with `deg(h)=58`.  This both reproduces the pole-58
obstruction and supplies the missing horizontal marking for the compact
q60/q80 neighbor.

The direct pole-58 deformation at CM-43 is now known to remain nontransverse
through third order.  The first-order level-79 section system has rank `358`,
and its quadratic obstruction has a two-dimensional reduced radical.  A
canonical second-correction slice produces a binary cubic whose rational
factor is

```text
dr/ds = 223593125/30934224.
```

but the second correction is only defined modulo the four-dimensional
first-order kernel.  Its polarized gauge map to the third-order cokernel has
rank one, and augmenting by the slice cubic still has rank one.  Thus the
whole cubic is absorbed: the displayed ratio is **not** an `H237` tangent.
The exact diagnostic explains why another pole-58 Taylor order is inferior to
the exact CM-24 q80 formal branches.  The later chamber audit shows that the
marked CM-43 q8/q60 classes themselves collapse to the old fiber.  See
[`scripts/compute_cm43_h237_tangent.sage`](scripts/compute_cm43_h237_tangent.sage)
and
[`scripts/verify_cm43_h237_tangent_crt.sage`](scripts/verify_cm43_h237_tangent_crt.sage).

The remaining H2 comparison equation is therefore the discriminant-237 divisor
inside this explicit `(r,s)` chart.  Its normalization is no longer the source
gate, because the H3 `H21 cap H92` component is already normalized exactly.

The first neighbor that is forced to use that high section is now exact.  A
complete constrained norm-120 enumeration has 56 sign-pairs and yields a
direct `q=60` transition

```text
Kumar E7+E8/MW2  --q=60-->  E8+E6/MW3.
```

For the selected `(a,b)=(5,12)` presentation, the new MW height Gram reduces
to `[[4,0,0],[0,20/3,1],[0,1,12]]`, of determinant 316.  The pinned child
frame and exact witness are recorded in
[`KUMAR_E7E8_BACKTRACK.md`](KUMAR_E7E8_BACKTRACK.md).  Its reducible fibers
are expected to be `II*+IV*`, with Euler budget six left for irreducible
singular fibers.  This is the first explicit compact neighbor, but the later
stability audit below finds a better deformation chart.

More concretely, the three reduced generators have zero-section intersections
`0,2,4`, and the corresponding short Weierstrass ambient space is

```text
Y^2 = X^3 + T^3*(a0+a1*T)*X
            + T^4*(b0+b1*T+b2*T^2+b3*T^3).
```

This is the decisive complexity reduction: six ambient coefficients and
largest section pole four replace the direct pole-58 Kumar system.

The two CM closures also survive this exact neighbor and locate boundary
divisors in the six-coefficient chart.  At `Delta=-24` the child frame is
`E8+E7+A2` with one MW generator of height 4; hence the generic `E6` fiber at
zero enhances from orders `(3,4,8)` to `(3,5,9)`, exactly the divisor `b0=0`.
At `Delta=-3` the child frame is `E8+E8+A2` with MW rank zero; the enhancement
to `(4,5,10)` is `a0=b0=0`.  The standard equation
`Y^2=X^3+T^5(T-1)^2` is the coefficient point `(0,0,0,1,-2,1)` and realizes
the extra `A2` as an `IV` fiber at `T=1`.  These statements are checked by
[`scripts/verify_e8e6_cm_boundaries.sage`](scripts/verify_e8e6_cm_boundaries.sage)
and [`scripts/classify_kumar_cm_frame_extensions.sage`](scripts/classify_kumar_cm_frame_extensions.sage).

The `0,2,4` pole profile is optimal among integral MW bases: vectors with
`P.O<=2` span rank only two.  More importantly, the height-4 section can be
solved in closed form.  After normalizing
`x(P1)=1+c1*T+c2*T^2+c3*T^3+T^4`, its endpoint jets determine `y(P1)` and all
six surface coefficients, leaving a three-parameter chart and only the two
rational sections with `P.O=2,4` to impose.  In the two sign charts the
linear form `ell=c3-epsilon*c1` divides `a0,b0,b3` to orders `1,2,3`, exposing
the CM boundary directly.  See
[`scripts/verify_q60_height4_normal_form.sage`](scripts/verify_q60_height4_normal_form.sage)
and
[`scripts/verify_q60_pole_profile_optimal.sage`](scripts/verify_q60_pole_profile_optimal.sage).

The discriminant-24 endpoint is now explicit in this compact chart.  A
finite-field seed modulo 31, followed by exact Hensel lift and rational
reconstruction, gives

```text
Y^2 = X^3 - (4096/19683)*T^3*(7*T+2)*X
            - (262144/14348907)*T^5*(T^2+34*T+19),
Delta = -(2^40/3^27)*T^9*(T-1)^3*(T^2+71*T+32).
```

Thus its fibers are `II*+III*+I3+2 I1`, and the displayed polynomial section
has height four.  Its `E8+E7+A2+<4>` frame is exactly isometric to the
transported `Delta=-24` closure.  See
[`scripts/verify_q60_delta24_anchor.sage`](scripts/verify_q60_delta24_anchor.sage).

## Leading noncollapsed chart: `q=80`, `E6+D5+A3/MW3`

The five proper `q=60` factor presentations were compared at both CM
closures.  The original `(5,12)` chart has optimal generic pole profile
`0,2,4` but retains only one MW direction at discriminant 24.  The `(4,15)`
chart retains two directions, but its `E6` still enhances to `E8`.

An exact Weyl-orbit search then tested every proper presentation using the
level-79 section through `q=80`.  This is a complete bounded search, not a
sample: for `q<237` that section has coefficient `+/-1`, so dominant-weight
enumeration in `E7+E8` replaces the enormous raw root shells.  Among 2,869
presentations there were no root-stable hits, but 313 had root jump one.

The best compact additive result is

```text
q=80, (a,b)=(4,20)
generic frame:    E6+D5+A3,     MW rank 3, optimal P.O=0,0,3
Delta=-24 frame: E6+D5+A3+A1,  MW rank 3, optimal P.O=0,0,0.
```

Thus the CM point adds only an `A1`; no generic section is lost.  The frame is
pinned in
[`data/fibrations/kumar_q80_e6_d5_a3_mw3_frame.txt`](data/fibrations/kumar_q80_e6_d5_a3_mw3_frame.txt),
and the exact orbit, CM transport, height lattices, component corrections,
and pole optimality are verified by
[`scripts/classify_kumar_cm_frame_extensions.sage`](scripts/classify_kumar_cm_frame_extensions.sage).

An exact `QQ` boundary surface in this chart is

```text
A=T^2*(-3+9/4*T-9/4*T^2+9/4*T^3),
B=T^3*(2-315/32*T+9*T^2-9/16*T^3-27/32*T^5).
```

Its discriminant is exactly

```text
(19683/1024)*T^7*(T-1)^4*(T+3)^2
              *(T^3+10/27*T^2+67/27*T-128/27),
```

so its fibers are `I1*+I4+I2+IV*+3 I1`.  Three short polynomial sections
with profiles `D1=(1,1,1,1)`, `D2=(0,0,1,2)`, and `D3=(1,1,0,2)` in
`(A1,A3,D5,E6)` order live over `QQ(sqrt(-3),sqrt(-6))`.  The exact lattice
basis change

```text
G1=D1, G2=D1-D2-D3, G3=-4*D1-D3
```

gives the transported generic basis with visible CM pole profile `(0,0,1)`
and resolved pair intersections `(2,3,1)`.  Inflating the cancelled quadratic
factor in `G3` restores the generic pole profile `(0,0,3)`.  The resolved
characteristic-zero tangent cone is
`(tau0-8/87*tau1)(tau0-1/12*tau1)`.  All surface identities, section
identities, pair intersections, and tangent calculations are certified by
[`scripts/verify_q80_cm24_rational_model.sage`](scripts/verify_q80_cm24_rational_model.sage).
After normalizing `p=9/4+h`, the corresponding surface tangents are

```text
(d',p',q',e')=(8/87,1,-24/29,-45/116),
(d',p',q',e')=(1/12,1,-45/52,-261/832).
```

Both full marked systems lift exactly over the CM compositum through order
18.  The certificate
[`scripts/extend_q80_cm24_branches_qq.sage`](scripts/extend_q80_cm24_branches_qq.sage)
checks every formal residual coefficient.  The order-18 jets also rule out
linear and quadratic implicit relations among `(d,p,q,e)` on either branch.
At order 36 the `8/87` branch has a two-dimensional exact kernel of centered
surface cubics; the `1/12` branch has no cubic relation modulo seven.  The
characteristic-zero cubic computation is reproduced by

```bash
sage scripts/extend_q80_cm24_branches_qq.sage \
  --slope=8/87 --order=36 --relation-degree=3 \
  --pair-max-degree=1 --quiet-series
```

These are exact finite-jet relations, not yet proved global equations: they
still require higher-order validation or direct substitution into an
algebraized family.  The asymmetry makes `8/87` the current recognition
target, but does not by itself prove that it is the `X(6,79)` branch.

The exact cubic plane reduces to the same two-dimensional cubic space as the
selected slope-`5` branch modulo seven, independently continued through order
145.  For that modular jet the filtered quotient dimensions through degrees
zero to five are

```text
1, 5, 15, 33, 48, 63.
```

Thus degrees three through five agree with the Hilbert polynomial `15*n-12`,
giving bounded formal evidence for a degree-15 projective coefficient image
of arithmetic genus 13. A subsequent modular plane projection has generic
degree five over `D` and its local singularity calculation gives normalization
genus zero. Therefore the known genus-two Shimura curve should not be
identified with the normalization of the **unmarked surface-coefficient**
image; it must enter through the marked-section cover, subject to the pending
characteristic-zero identification. These remain modular/formal computations,
not a proved global characteristic-zero model.

The order-85 computation separates the two cubics, their ten affine-linear
multiples, and twelve new quartic generators.  Its centered affine and
homogeneous `(z,D,P,Q,E)` forms are exported in
[`../artifacts/generated-results/q80-cm24-slope-8-87-gf7-ideal.json`](../artifacts/generated-results/q80-cm24-slope-8-87-gf7-ideal.json),
SHA-256 `d620f3443551ca11c510de3be5776ea9a2180f4a9a21b40ca7a941223eba08b8`.
This artifact is explicitly labelled bounded `GF(7)` formal evidence. The
candidate ideal is analyzed reproducibly by
[`scripts/analyze_q80_rank19_branch_ideal.sage`](scripts/analyze_q80_rank19_branch_ideal.sage).
Its `(D,Q)` eliminant is irreducible, has degree five over `D`, and the full
four-coordinate candidate also has generic `D`-fiber length five. Thus this
projection is birational at the level of the modular candidate function
field. Local delta invariants total `21` on its degree-eight plane closure,
so that normalization has genus zero; the affine normalization alone has
delta `11`. Normalizing the homogeneous coordinate ring gives a degree-eight
rational normal curve with Hilbert polynomial `8*n+1`. Its osculating flag at
CM-24 supplies an explicit parameter `t=0` there, and exact substitution gives
rational functions of numerator/denominator degrees

```text
D: 5/4,  P: 10/8,  Q: 8/6,  E: 15/12.
```

All fourteen modular ideal generators vanish identically after substitution.
The deterministic parameter artifact is
[`../artifacts/generated-results/q80-cm24-slope-8-87-gf7-parameter.json`](../artifacts/generated-results/q80-cm24-slope-8-87-gf7-parameter.json),
SHA-256 `ed35aab6a6be653985cc342fe79b067971e27c84babf44446bae92c37db4a166`.
This is still a certificate for the bounded modular candidate ideal, not yet a
characteristic-zero family.

The marking test now sharply separates this modular candidate from the desired
`X(6,79)` source.  The component gates force the first polynomial section to
have `X1=T+(d-1)T^2`.  After substituting the rational surface parameter, its
full Weierstrass identity is an exact square over the quadratic cover

```text
v^2 = t^5 + 5*t^4 + 3*t^3 + 4*t + 1                 (mod 7),
```

which has genus two and absolute Igusa invariants `(4,2,5)`.  Elkies's known
model `u^2=16*t^6-19*t^4+88*t^2-48` has good reduction with invariants
`(4,4,3)`, so the two covers are not geometrically isomorphic by a base Mobius
change.  The five bounded rational candidates for the second polynomial
section's `X`-coordinates also pass their complete degree-twelve square
identity exactly.  That section requires a distinct genus-three quadratic
cover.  The two branch divisors overlap in degree five; their combined marked
curve is therefore a genus-six biquadratic cover of the rational surface line.
Its third quadratic quotient is

```text
w^2 = t^4 + 3*t^3 + 3*t^2 + 1,
```

with `j=0`, whereas the two elliptic quotients of `X(6,79)` reduce to
`j=1,3`.  Thus neither individual short-section cover, nor their biquadratic
compositum, is the missing source hyperelliptic coordinate.  Exact lattice
transport now identifies the horizontal Kumar markings in the generic q80
optimal basis `(G1,G2,G3)` as

```text
height-4 = -G2,
Q79      = -3*G1 - 2*G2 + 4*G3.
```

In the q80 fibration these classes have heights `4,120`; `Q79` has `P.O=59`
and necessarily uses the pole-bearing `G3`.  Thus direct construction of
`Q79` would merely recreate the old large-pole bottleneck.  Conversely, the
height-four section must descend to the rational Shimura quotient, while the
modular slope-`8/87` candidate makes `G2` live on the nontrivial genus-three
cover above.  This rules out that two-short-section descent, but not the
unmarked `8/87` surface because the true level-79 class uses `G3`.  The later
split-prime point counts below conditionally reject `1/12`; therefore `8/87`
remains the live unmarked candidate, but only with a different `P3/Q79`
marking or orientation.  The exact horizontal map and q80 height/pole
profiles are checked by
[`scripts/classify_kumar_cm_frame_extensions.sage`](scripts/classify_kumar_cm_frame_extensions.sage).
The modular cover checker is
[`scripts/analyze_q80_rank19_marked_cover.sage`](scripts/analyze_q80_rank19_marked_cover.sage).
Its bounded input candidates are exported in
[`../artifacts/generated-results/q80-cm24-slope-8-87-gf7-partial-marked-parameter.json`](../artifacts/generated-results/q80-cm24-slope-8-87-gf7-partial-marked-parameter.json),
SHA-256 `83e72d38223561381c29a766a26ec66b667355b6391318fee5f76dfdeb76fb7f`;
the global square identities, rather than the finite jets, promote precisely
the displayed `X`-coordinate candidates.

The first characteristic-zero support test also supplies a useful negative
result.  The modular `(D,Q)` eliminant has 26 nonzero monomials, but the exact
order-30 branch has kernel dimension zero on precisely that support (four
extra jet coefficients are retained for validation). Hence some coefficients
vanish accidentally modulo seven, and the modular 26-term equation cannot be
lifted coefficient-for-coefficient. The next gate is to recover a rational
parameter or a slightly larger characteristic-zero plane support, then lift
the marked-section extension and identify that cover with the known genus-two
function field.

The first marked-section extension is a genus-two double cover modulo seven,
but its absolute Igusa invariants `(4,2,5)` differ from `(4,4,3)` for the
known source model. The second marked section lies on a distinct genus-three
cover; the combined marked curve has genus six and its third quadratic
quotient has genus one. Hence the two-short-section slope-5 descent is a
high-genus detour, not the desired source marking. This does not reject the
unmarked slope-5 surface branch: the true level-79 coordinate may involve
`P3`. Its fourteen surface generators pass five withheld orders through
order 90, so this is not merely an
order-85 fitting failure. Its P1 branch set has trivial stabilizer over
`GF(49)` despite the known source's bielliptic involution, and its genus-one
third quotient has `j=0` rather than either known value `1,3`.

The slope-3 branch now has a cross-prime bounded global candidate. No centered
relation occurs through degree five. Over `GF(7)`, degree six gave seventeen
relations at order 214 and all passed sixteen withheld coefficients through
order 230, but split-prime replay shows that two are exceptional: ordinary
split primes have rank 195/kernel 15, `p=7` has kernel 17, and `p=97,103`
have kernel 16. The generic fifteen-sextic modular ideal has affine dimension
one, projective degree 48 and arithmetic genus 94, and an irreducible `(P,D)`
projection of total degree 32 and bidegree `(15,32)`. Thus the full
seventeen-generator `GF(7)` ideal must not be lifted wholesale. At both
`p=73` and `p=79`, exact local delta calculations give affine/infinity
contributions `223+241=464`; the irreducible degree-32 plane curve therefore
has normalization genus one. This is the expected elliptic quotient shape.
However its traces are `15,22` at `p=73,127`, while the known source elliptic
factors have traces `(-9,7)` and `(8,0)` there (including twists still gives
no match). Conditional on the multi-prime fifteen-sextic ideal being the true
global branch, slope `1/12` is therefore not the desired component. Return to
slope `8/87` and recover the actual `P3/Q79` marked quadratic cover; the
high-genus P1/P2 compositum does not rule out that marking.
The first bounded order-230 `P3` reconstruction supplies an additional
negative gate: both its pole and one numerator coordinate lie on the same
genus-one `j=0` field as the third P1/P2 quotient, with squarefree model
`w^2=t^4+3t^3+3t^2+1`.  Thus the selected three-section orientation still
does not yield the known genus-two source.  This is not an unmarked-surface
rejection; it requires either a different CM marking/orientation or direct
descent of the combination `Q79`.
See the canonical route audit
[`Q80_TO_ROOTLESS_PATH_2026-08-21.md`](Q80_TO_ROOTLESS_PATH_2026-08-21.md)
and analyzers
[`scripts/analyze_q80_rank19_marked_cover.sage`](scripts/analyze_q80_rank19_marked_cover.sage)
and
[`scripts/analyze_q80_rank19_branch_ideal.sage`](scripts/analyze_q80_rank19_branch_ideal.sage).

Independently, the missing lattice continuation from this q80 frame to the
rootless rank-17 fibration is now exact. The six-step path has neighbor norms
`4,4,12,12,4,6` and MW ranks `3,4,5,6,13,16,17`; its terminal frame is
integrally isometric to the pinned determinant-948 rootless lattice. The
complete witness table and composite Neron--Severi transport are certified in
[`Q80_TO_ROOTLESS_PATH_2026-08-21.md`](Q80_TO_ROOTLESS_PATH_2026-08-21.md)
and
[`scripts/verify_q80_to_rootless_path.sage`](scripts/verify_q80_to_rootless_path.sage).
This closes the lattice route.  A second downstream suffix has now also been
executed through the final equation gate over `GF(73)`: the productive
pair-23 fifth model has root data `(16,66,2048)`, and Smith saturation of its
final q6 chord module gives `diag(1,s-27)`,
`q_sat=(q0+63)/(s-27)`, kernel rows `(1,0,0,41,48)` and
`(0,1,0,6,72)`, and a unit-preserving Jacobian with fibers
`2I5+2I4+I2+4I1` and root data `(15,66,800)`.  The certificate
[`../artifacts/generated-results/q80-final-q6-saturated-module-gf73.json`](../artifacts/generated-results/q80-final-q6-saturated-module-gf73.json)
has SHA256
`7d3866855b9995b733193de9c5d5e3ba1cea6aa1292141e06d0d5011f28975e3`;
it includes `q0(s,R)` and globally minimal polynomial `A`, `B`, and
discriminant coefficients.
Characteristic-zero lifting and field-of-definition tracking remain open.

This q80 construction uses the `H2=diag(4,237/2)` polarization; it must not
be identified with the exact source polarization
`H3=[[21/2,3],[3,46]]`.  Their doubled forms have equal determinant but
different minima and are not integrally isometric.  The recurring q=8 signal
comes from the height-four Humbert-8 direction or special
Noether--Lefschetz return/collapse behavior, not from the degree-one Smith
factor in the final q6 local module.

The equation ambient is also exact.  Put `I1*`, `I4`, and `IV*` at
`0,1,infinity`.  After Weierstrass scaling it has the four-parameter form

```text
A=T^2*(-3+p*T+q*T^2+(-3*d^2+3-p-q)*T^3),
B=T^3*(2+b1*T+b2*T^2+b3*T^3+b4*T^4+e*T^5),
```

where the four `I4` jets determine `b1,...,b4` linearly.  Its discriminant is
`T^7*(T-1)^4*R5(T)` up to a nonzero scalar.  The discriminant-24 boundary is
simply `disc_T(R5)=0`, the extra `I2`.  The formulas are derived and verified
by
[`scripts/verify_q80_ambient_normal_form.sage`](scripts/verify_q80_ambient_normal_form.sage).
The CM sections and both rank-19 formal branches are now recovered.  The
remaining equation-level task is to select the `X(6,79)` branch using the
global marking and recognize it by a common modular function or higher-degree
implicit relations.

The section marking is now exact.  In an optimal generic basis the height
Gram, component profiles in `(A3,D5,E6)` order, zero intersections, and pair
intersections are

```text
H = [[2/3,0,3/4],[0,4,2],[3/4,2,37/4]]
P1=(1,1,1), P2=(0,0,0), P3=(3,0,0)
P.O=(0,0,3), (P1.P2,P1.P3,P2.P3)=(2,4,3).
```

At discriminant 24 this same transported basis has

```text
H = [[1/6,0,-3/4],[0,4,2],[-3/4,2,19/4]]
profiles in (A1,A3,D5,E6):
P1=(1,1,1,1), P2=(0,0,0,0), P3=(1,3,0,0)
P.O=(0,0,1), (P1.P2,P1.P3,P2.P3)=(2,3,1).
```

Thus the smooth deformation chart must retain a cubic denominator for `P3`:
its visible linear denominator at the CM point is caused by cancellation of a
quadratic factor.  A different CM basis has three polynomial sections, but it
deforms with pole profile `(0,16,12)` and is therefore the wrong lifting
basis.  These claims, including the exact generic-to-CM basis matrix, are
checked by [`scripts/classify_kumar_cm_frame_extensions.sage`](scripts/classify_kumar_cm_frame_extensions.sage).

The earlier standard-`P1` modular seeds at `p=7,11,29,53` are retracted.  The
compiled `P2` gate omitted avoidance of the residual `I2` node: the proposed
`P2` passed through that node, giving nonzero `A1` component label when the
target label is zero.  After adding this condition and the exact pair gates,
all of those candidates disappear.

The standard branch itself has nevertheless been solved exactly:

```text
q=18-2*p,  e=(p-42)^2/36,
P2 eliminant=(p+105/8)*(p^2-144*p+7371).
```

The rational root gives an explicit rank-19 surface, but every modular `P3`
seed tested at `11,29,53` has a full-rank first Hensel obstruction.  A second
mod-11 seed avoiding the residual node lifts as a one-dimensional family but
fails the required `P1.P2` and `P2.P3` gates.  Pair intersections are now
checked on the resolved surface by `P.Q=(P-Q).O`: the chord numerator and
denominator are square-cancelled, and the remaining denominator degree is the
exact gate.  This avoids raw-gcd overcount at singular Weierstrass nodes.

With the node and resolved-pair gates in place, complete searches at
`p=11,13,17,19` have no marked hit.  At `p=7` there is exactly one fully
gated hit.  It is the reduction of the exact rational CM-24 surface and its
marked short-section basis:

```text
(d,p,q,e)=(3,4,3,2), rho=4, node=1,
P1: X=[0,1,2,0], Y=[0,0,5,6,3,0],
P2: X=[1,2,2,6,4], Y=[6,4,4,5,4,4,1],
P3: Z=T-6,
    X=[1,4,3,0,6,1,4], Y=[1,6,2,5,3,4,1,2,3,1].
```

The former 88-equation pair-cancellation lift obscured this point with six
vertical cancellation directions.  The minimal pole-one chart removes those
directions.  Its exact characteristic-zero Jacobian has rank 37 and nullity
two, and its Kuranishi equation is

```text
tau0^2-(61/348)*tau0*tau1+(2/261)*tau1^2
 = (tau0-8/87*tau1)*(tau0-1/12*tau1).
```

Reducing the two rational slopes modulo seven gives exactly `5` and `3`, the
two finite-field branches that lift uniquely through `h^20`.  The finite-field
Pade search remains a useful bounded negative result, but the existence and
splitting of the branches are now characteristic-zero facts.  See
[`scripts/verify_q80_rank19_deformation_gf7.sage`](scripts/verify_q80_rank19_deformation_gf7.sage)
and
[`scripts/extend_q80_rank19_branches_gf7.sage`](scripts/extend_q80_rank19_branches_gf7.sage).

## Secondary CM chart: non-isotrivial `E8+A2^3`

The discriminant-3 CM endpoint supplies a useful independent geometric
cross-check.  The determinant-948 rank-19 lattice embeds exactly in its
Neron--Severi lattice, and the inherited abstract root frame is
`E8+A2^3`, with exact MW rank three and reduced height Gram

```text
(1/3) * [[8,-1,0],[-1,10,0],[0,0,12]].
```

However, the former Kodaira lift `II*+3 IV+II` is impossible.  Its equation
has `j=0`, so the order-three CM automorphism forces the geometric MW rank to
be even; the target requires `19-2-14=3`.  An `A2` root factor may instead be
an `I3` fiber.  At least one target factor must be `I3`, and the correct short
Weierstrass family must allow nonzero `A(t)`.  This exact obstruction is
recorded in
[`E8_A2_KODAIRA_CORRECTION.md`](E8_A2_KODAIRA_CORRECTION.md).

This identifies a viable secondary chart, but it is not by itself the
source-level deconstruction: an abstract CM embedding does not identify the
rational polarization, Kodaira realization, or Galois descent used by Elkies.

The first viable ambient deformation is now explicit.  With `D=t(t-1)`,

```text
A = -3r^2 D^2,
B = D^2((t-lambda)^3-2r^3D)
```

has exact fibers `II*+2 IV+I3+3 I1`, nonconstant `j`, and specializes at
`(r,lambda)=(0,0)` to the CM model.  Its discriminant factors without any
elimination.  This is an ambient candidate, not yet an identification of the
target one-dimensional locus; see
[`E8_A2_MIXED_FAMILY.md`](E8_A2_MIXED_FAMILY.md).

Exact glue recovery also fixes a reduced generator basis with `A2^3`
component profiles

```text
P1=(1,1,0),  P2=(0,2,0),  P3=(0,0,0),
```

up to the natural symmetries of the first two factors.  Their nonzero counts
`(2,1,0)` recover heights `(8/3,10/3,4)` and force `P_i.O=0` for all three.
Thus the corrected equation search can use three polynomial sections with
pinned component gates and pairwise intersections all equal to `2`.

## Why direct reconstruction was reduced

The naive symbolic system (`build_rank17_system.sage`) introduces `A`, `B`, and
17 quartic/sextic sections simultaneously and is too large for direct Groebner
work.

The recovered shell contains a 9-vector clique with Gram

```text
4 on the diagonal, 2 off the diagonal.
```

Write its oriented sections as `V_0,...,V_8`.  Since
`<V_i,V_j>=2`, every difference

```text
D_ij = V_i - V_j
```

is again minimal.  Thus the clique immediately gives 45 known minimal
x-coordinate classes: 9 vertices plus 36 differences.

For every minimal additive triple `P+Q=R`, the chord slope is quadratic and

```text
x_P + x_Q + x_R = m_PQR^2.
```

For the Coxeter subsystem this gives exactly 120 square relations on 45
quartic x-polynomials.  `build_coxeter9_x_reconstruction.py` constructs the
120x45 incidence matrix and verifies every relation against the global 184242
triple catalog.

Observed exact/numerical rank data:

```text
incidence_shape = 120 x 45
numeric_rank = 45
rank mod 5,7,11,... = 45
exact full-column rank over Q = certified
left nullity = 75
```

Hence all 225 x-coordinate coefficients can be eliminated linearly once the
quadratic slopes are known.

## Coherent-slope reduction

The 84 triangle slopes are not independent.  If `m_ij` is the slope through
`V_i` and `-V_j`, then

```text
slope(D_ij,D_jk) = m_ik - m_ij - m_jk.
```

Therefore the nonlinear reconstruction can be expressed using only the 36 pair
slopes, i.e. 108 scalar coefficients before gauge fixing.  The numerical gauge
fixes four coefficients, leaving 104 free Coxeter variables.

Given the slopes, the implementation reconstructs in sequence:

```text
36 pair slopes
    -> 9 quartic x_i
    -> 9 sextic y_i
    -> one common degree-8 A(T)
    -> one common degree-12 B(T).
```

`solve_coxeter9_slopes_numeric.py` repeatedly finds roots with algebraic
residuals around `1e-13` to `1e-16`.  Many roots lie on or close to the
cuspidal/isotrivial component, so residual size alone is not a useful ranking
criterion.

## Coxeter-root diagnostics

`analyze_coxeter9_numeric_roots.py` ranks roots using:

- scale-aware `A/B` and discriminant strength;
- effective degrees of `A`, `B`, and the discriminant;
- the polynomial non-isotriviality invariant

  ```text
  3 A' B - 2 A B';
  ```

- numerical Jacobian nullity.

Two roots became especially important:

### `root-000001`

```text
raw ~ 6.97e-15
non-isotrivial
(deg A, deg B, deg Delta) = (8,12,24)
nullity at 1e-8 ~ 9
```

This is an excellent generic Coxeter-scaffold point.  However, attempts to
impose the first rank-10 extension on its initial lattice labeling drove the
surface strongly toward the isotrivial boundary.

### `root-000029`

```text
raw ~ 8.14e-13
non-isotrivial
(deg A, deg B, deg Delta) = (8,12,24)
discriminant diagnostic ~ 0.18
```

This root became the best testbed for the first extension, but the first
apparently machine-precision rank-10 hit on it was later found to be a
near-collision and is not yet an independent-section certificate.

## Exact rank-17 extension chain

Starting with the raw 9-vector Coxeter clique, the combinatorial extension
search selected the following signed minimal vectors:

```text
2313 2525 307 1303 1859 1441 683 2351 2143
961 2402 1642 1300 1023 2216 2392 2610
```

The final matrix has:

```text
rank = 17
coordinate determinant = 2
77 pairings with absolute value 2
additive closure = 1396/2622 before explicit saturation bridging
```

The determinant `2` is **optimal**, not a failure.  The raw Coxeter clique has
Gram determinant 5120 while its saturated rank-9 lattice has determinant 1280,
so the raw clique itself has saturation index

```text
sqrt(5120/1280) = 2.
```

Any 17-vector coordinate matrix retaining those nine raw generators must
therefore have determinant divisible by 2.  The selected chain attains this
minimum, and its saturation is the complete recovered rank-17 lattice.

`audit_rank17_extension_chain.py` exposes the index-2 parity class, searches a
minimal saturation bridge, and exports exact `|pairing|=2` continuation anchors.

## Rank-10 continuation: first attempts

The first continuation used section `961` with only its `|pairing|=2` links to
the nine raw Coxeter generators.  It found excellent algebraic fits, often
`1e-9` or better, but these collapsed toward the isotrivial boundary.  Adding
all nine available anchors against the full 45-section Coxeter shell improved
the algebraic fit further but showed a clear Pareto split:

- algebraic residual near `1e-10`/`1e-11` -> discriminant and j-variation
  collapse by many orders of magnitude;
- healthy non-isotrivial branch -> algebraic residual only around `1e-5`.

The line/group-law formulas were checked; the remaining issue was not a sign
bug but a symmetry/labeling ambiguity and, as found later, a collision escape.

## Coxeter S9 ambiguity

The raw clique Gram is invariant under permutation of its nine generators.
Thus a numerical Coxeter solution has no intrinsic labeling telling us which
numerical `V_i` corresponds to which full rank-17 lattice generator.

For lattice section `961`, the original generator-pairing signature is

```text
-1 0 -1 0 -1 -2 -2 -2 -1
```

Its multiset has 1260 distinct `S9` permutations; including the opposite
orientation of the new section gives 2520 distinct fingerprints.  These were
scanned exhaustively on fixed Coxeter roots by
`scan_rank10_coxeter_fingerprints.py`.

### Root-000001 scan

Best fixed-surface result:

```text
score ~ 4.412e-7
```

No compelling rank-10 hit appeared.

### Root-000029 scan: apparent hit, later downgraded

The exhaustive scan found:

```text
fingerprint index = 2452
V pairings = 0 0 2 2 2 1 1 1 1
anchors = 9
curve residual = 5.220e-14
line residual = 5.850e-13
score = 5.850e-13
```

The winning fingerprint is a permutation of the **negative** of the original
section-961 signature.  This left 288 compatible full-lattice Coxeter mappings.

A joint refinement then reached

```text
base residual  ~ 3.36e-13
curve residual ~ 7.73e-16
line residual  ~ 4.11e-14
healthy discriminant / non-isotrivial diagnostics
```

However, `validate_rank10_independence.py` subsequently found

```text
nearest_x = V3
x_distance ~ 2.52e-7
```

and the intended quadratic-line identities were not numerically separated from
many accidental non-target line fits.  Six intended edges had absolute errors
around `3e-14` to `4e-14`, but several non-target sections had comparable
absolute errors.  The original validator's strict relative threshold therefore
marked only 3 of 9 target edges; simply relaxing that threshold would create
many false extra edges.

The key conclusion is that the original `5.85e-13` fixed-surface hit and its
joint refinement are **near-collision candidates, not yet evidence of an
independent tenth section**.  Their small residuals can be explained by the
ill-conditioned gauge together with `Q` approaching an existing Coxeter
section.

The rank-10 tangent/Jacobian experiments are also not being used as evidence at
this checkpoint.  Several versions showed that strongest-SVD-gap nullities are
parameterization- and scaling-sensitive, and an incremental-codimension test
was initially invalid because numerical tangent leakage was amplified by
column normalization.  These diagnostics remain useful for conditioning work,
but they are secondary to direct section/incidence validation.

## Current experiment: nondegenerate rank-10 refinement

`refine_rank10_nondegenerate.py` now repeats the winning-fingerprint refinement
with explicit structural guards:

1. `Q` must stay a scale-free positive distance from all 45 known Coxeter
   sections, using the median pairwise Coxeter x-distance as normalization;
2. the discriminant and non-isotriviality branch must remain healthy;
3. all nine required target line identities are imposed;
4. after solving, required edges must be numerically separated from all
   non-target and wrong-sign quadratic-line fits.

The principal threshold-free validation statistic is

```text
edge_gap = min(non-target line error) / max(required-edge line error).
```

The previous near-collision solution has edge gap of order one.  The new solver
requires, by default, both target/non-target and correct/wrong-sign gaps of at
least `1e2`, in addition to the collision and branch guards.

A successful result should therefore simultaneously have:

```text
small algebraic residual
healthy discriminant / j-variation
distinct_ratio >= 0.02
edge_gap >= 1e2
wrong_sign_gap >= 1e2
```

If no such solution exists for the current fingerprint/root, the next step is
to rescan the 2520 Coxeter fingerprints with the same non-collision and
incidence-separation criteria rather than interpreting the old machine-small
residual as rank growth.

## Next steps

Primary H3 geometric route:

1. Identify the rational H21 entrance data with the signed H3 q6 divisor
   `O+(-P1)-F`, and determine the actual field of definition of the
   height-`21/2` section.  The equality of the H21 and H92 oriented square
   classes is compatibility evidence, not a section-descent proof.
2. Execute the certified degree-two pencils in characteristic zero through
   `E8+E6/MW3`, `D13/MW4`, `D12/MW5`, `A11/MW6`, `2A5/MW7`,
   `3A3/MW8`, and `A3+2A2/MW10`.
3. Continue the Weyl-quotient search from that single MW10 frame to a rootless
   MW17 frame, then construct the remaining divisor functions and track all
   section and component fields.
4. Only after the rootless family is explicit, solve for a rational parameter
   and a `QQ`-isomorphism to the record curve and literally transport all
   generic sections.

Parallel H2/q80 route:

1. Compute the discriminant-237 Noether--Lefschetz divisor in the explicit
   Humbert-8 `(r,s)` Kumar chart.  Normalize it with the exact Inose anchors at
   `t=infinity` (discriminant `-3`) and `t=+/-2` (discriminant `-24`), and
   constrain it with the CM-`43` point
   `(-1225/722,-93312/442225)`.  Then recover `(r(t),s(t))` together with the
   double-cover coordinate `z/u`.
2. Do not construct q=60 from CM-43.  Full fixed-component reduction gives
   `D60_raw=F+initial_fixed+4(P1-P2)+(P3-P2)` and
   `D8_raw=F+initial_fixed+(-P2)+(P1-P2)`, so both marked classes collapse to
   the old fiber on that CM surface.  The q=8 nonmembership and the q8/q9
   factorization remain useful marking certificates, not deformation seeds.
3. Continue from the exact q80 `(4,20)` CM-24 model, which retains all three
   generic MW directions.  Its characteristic-zero tangent slopes are
   `8/87` and `1/12`, and both marked systems lift exactly through order 18;
   modulo seven they lift uniquely through `h^20`.  First select the H237
   branch using its global marking. Do not use the two short slope-5 sections
   alone as the source coordinate: their fields combine to genus six, while
   the true marking may involve `P3`. The selected P3 coordinates still lie
   on the genus-one `j=0` quotient, so change the CM orientation or descend
   `Q79` directly. The generic fifteen-sextic slope-3 candidate has wrong
   elliptic traces at two good primes and is conditionally rejected; do not
   invest further in that finite-jet ideal without a different global marking.
   Impose the discriminant/component equations directly. Derive the compact
   q60 `(4,15)` chart only if its
   smaller system offsets the loss of one MW direction.
4. Execute the exact `q=4,4,12,12,4,6` q80-to-rootless path sequentially on
   the algebraized model. Reduce each raw fiber to the nef chamber and track
   all section and component fields; the exact lattice endpoint is already
   the pinned MW17 frame. The first q4 class has old degree two, zero MW
   projection, and exact class `4F+2O` minus displayed `D5+E6` components. A
   rank-three height argument and root-primitivity proof make it fully nef.
   The local component resolution now gives the four linear gates on
   `L(2O+4F)` and the exact pencil `U=(x-T)/T^2`. Its ambient child has
   `D9+A3`; the marked rank-19 collision gives the pinned `D9+A4`, and CM24
   specializes further to `D9+A5+2A1`. The second q4 class is also fully nef:
   its exact saturated MW height shell checks every section and root
   primitivity excludes a negative bisection. It is again old-degree two with
   zero MW projection, namely `4F+2O` minus
   `(1,2,2,1)_A4+(1,4,4,4,2,2,3,4,2)_D9`. Its exact coordinate on the first
   child is `W=(X-3v^3-x1v-x0)/v^2`, where `v=U-d+1`,
   `x0=-3B1(d-1)/(2A1(d-1))`, and `x1=-A1'(d-1)/(6x0)`.
   At CM24 it gives `D7+E6+3A1/MW2`; the finite `(3,4,8)` valuation is
   `IV*`, resolving the earlier false `I2*` classification. Continue with the
   third q12 step. Its chamber-reduced divisor has old degree three and
   nonzero integral MW projection of norm eight. This gives a height-eight
   section `S` with `S.O=2` and exact decomposition
   `D3=S+2O+2F+root_correction`, so the next equation is an actual MW-marked
   trisection in `L(S+2O+2F)` rather than another compensated zero-MW pencil.
   At CM24 the transported class is the polynomial section
   `Q_CM=P1+3P2` over `QQ(sqrt(-6))`, of height three and `P.O=0`; the earlier
   `2P1-P2` identification was a basis-order error. Thus the CM generic-fiber
   space is exactly `<1,X,(Y+y(Q_CM))/(X-x(Q_CM))>`, and only the vertical
   gates remain. The unique effective norm-four lift gives
   `D=Q_CM+2O+4F+R`, with integral component coefficients
   `A1:(0),(1),(0)`, `E6:(1,2,3,2,1,2)`, and
   `D7:(2,4,3,3,5,6,3)`; the extra `2F` compensates the CM pole drop. Exact
   CM transport fixes the target as
   `2A6+3A1/MW3`. Exhausting the pure translations
   `z_Q+k2*W^2+k1*W+k0` gives only `D8+E7+A1` and
   `D5+E6+A3+A1`, so a nonzero `X` coefficient is required. The alternate
   polynomial ansatz `(a1*W+a0)*X+b2*W^2+b1*W+z_Q` is also excluded by all
   2,401 coefficient tuples over `GF(7)`: its nonzero-X branch degree is 17
   or 23, never genus one. The correct nine-dimensional ambient module is
   `a(W)+bX+c(W)z_Q`. Its two selected-I2 rows and one E6 row combine with
   four exact D7 rows from the complete local ideal `(Y,U^2,ZU,Z^3)`. The
   resulting rank-seven kernel gives `Vnew=N1/N0`; clearing the chord leaves
   `X-Qx` times a cubic of old-W degree nine. At `p=73`, generic fibers
   `Vnew=1` and `Vnew=7` are irreducible and normalize to genus one. The
   earlier unresolved-cusp rows gave genera four and three and are rejected.
   Canonical-origin interpolation from 49 fibers plus seven withheld checks
   recovers its short Jacobian over `GF(73)(V)`, with
   `deg(A),deg(B)=(8,12)` and discriminant `2I7+3I2+4I1`. Thus the CM24
   third q12 pencil and the transported `2A6+3A1/MW3` fiber marking are exact
   modulo 73; lifting this compact model to characteristic zero is the next
   gate. Do not impose vanishing at both quadratic roots; that came from
   an ineffective lift. The alternate
   second-child q8 class is genuine and noncollapsed (`E6+A7/MW4`,
   specializing to `E7+A7+2A1/MW2`), but it loses a generic MW direction
   rather than gaining one; retain it only as a possible degree-two
   coordinate normalization, not as a replacement for the q12 rank-growing
   step.
5. Determine how the double-cover coordinate `u` selects the rational
   quadratic twist and the level-79 section of the involutive middle
   `E7+E8/MW2` frame.  Track
   the zero section, components, fields of definition, and section classes.
6. Retain the retraction of the earlier q80 modular seeds, which were false
   residual-`I2` profile hits.  Execute only a branch that passes the corrected
   component and pair gates on the genuine Kumar Weierstrass model.  Only then
   choose an E6 base normalization and, if still useful, build its
   P1/P2 system with `(P1+P2).O=1` imposed from the start.
7. Verify the characteristic-zero Weierstrass identities, Kodaira fibers,
   full height lattice, and exact specialization/Q-isomorphism to the record
   curve.

Secondary CM-chart route:

8. Use the exact mixed/all-`I3` `E8+A2^3` families only to cross-check the
   Kumar model or if the Shimura-to-Igusa map must be reconstructed from CM
   interpolation.  Do not enlarge their finite-field searches without first
   matching a polarization and rational descent.

Fallback direct-Coxeter route:

9. If the corrected upstream route fails its exact lattice gates, rescan the 2520 Coxeter fingerprints
   with collision guards and edge-gap scoring rather than reviving the old
   near-collision.
10. Use the index-2 saturation bridge as a gluing constraint for any genuine
   continuation through sections
   `961,2402,1642,1300,1023,2216,2392,2610`.

Final certification remains unchanged: exactify all coefficients, reconstruct
the complete 1311-pair shell, store the canonical model under `data/k3-model/`,
and only then resume specialization searches.

## Current claim boundary

What is exact today:

- the rank-17 lattice and its 1311-pair minimal shell;
- the 184242 additive relations;
- the Coxeter-9 decomposition and incidence-rank reduction;
- the coherent-slope algebra;
- the determinant-2 optimal rank-17 extension chain and its saturation
  interpretation.

What is numerical today:

- reconstruction of many Coxeter-9 elliptic K3 points;
- the non-isotrivial `root-000029` scaffold;
- a very small-residual **near-collision** candidate for the first extension,
  now explicitly downgraded pending nondegenerate validation.

What is exact over a finite field today:

- a complete E6/P1 closed-equation reconstruction over `GF(31)`, subsequently
  identified as an `I5` boundary point;
- nonsingularity of its reduced `8 x 8` closed-equation Jacobian;
- dependency `P1=-2*P3` on that boundary branch;
- completeness of the seven common-factor candidates on the fixed coordinate
  slice, their exact reconstruction audit, and the degree reduction
  `X2=C+q0*F`;
- two exact target-fiber P1 surfaces on that slice, with zero canonical P2
  hits after the deeper component test.

What remains unproved:

- an independent generic tenth section realizing the intended fingerprint;
- an exact Weierstrass model for the intended rank-17 fibration;
- exact generic coordinates for sections 10 through 17;
- any new rank-21 specialization arising from this family.
