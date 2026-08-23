# Correct upstream backtrack for the rank-17 K3

## Status

This note corrects the order of the geometric reconstruction.  The exact
rank-17 Mordell--Weil lattice and the previously recovered neighbor frames
remain valid lattice computations.  What was not justified was treating the
CM-derived `E8+A2^3` frame as the primary equation-level entrance.

The primary construction instead supplies a canonical Dolgachev--Kumar
fibration with root lattice `E7+E8` and Mordell--Weil rank two.  Run the exact
finite classification with

```text
sage elkies-k3/scripts/classify_kumar_e7e8_anchors.sage
```

## Source-level construction order

For the non-CM rational point on the genus-two Shimura curve

```text
X(6,79)/<w_(6*79)>,
u^2 = 16 t^6 - 19 t^4 + 88 t^2 - 48,
|t| = 14/13,
|u| = 2^6*251/13^3,
```

the construction described in the primary sources runs as follows:

```text
principally polarized QM abelian surface A
  -> Dolgachev--Kumar K3 surface S_A
  -> canonical E7+E8 elliptic fibration, MW rank 2
  -> elliptic-neighbor transformations
  -> rootless elliptic fibration, MW rank 17.
```

The canonical Kumar equation is

```text
Y^2 = X^3 + (a*t^4 + a'*t^3)*X
            + (b''*t^7 + b*t^6 + b'*t^5),
```

with the `E7` and `E8` fibers at zero and infinity.  After normalizing
`a'=-1`, its parameters recover the Clebsch--Igusa invariants by

```text
(I2,I4,I6,I10) = (-24*b'/a', -12*a,
                  96*a*b'/a' - 36*b, -4*a'*b'').
```

For a non-CM surface with QM of level datum `N`, this fibration has MW rank
two and regulator `N`.  Here `N=6*79=474`; hence Shioda's determinant formula
gives

```text
det(NS) = det(E7+E8) * regulator = 2*474 = 948.
```

The formula and construction order are recorded in the local primary-source
audit summarized by
[`../elliptic-curves/ELKIES_RANK18_SOURCE_RECOVERY_AUDIT.md`](../elliptic-curves/ELKIES_RANK18_SOURCE_RECOVERY_AUDIT.md).

## Exact finite classification of Kumar anchors

Write the rank-two height Gram as

```text
H = (1/2) * [a b; b c].
```

Then `det(H)=474` is equivalent to `a*c-b^2=1896`.  Every positive binary
form has a reduced representative satisfying

```text
0 <= 2*b <= a <= c,
a <= sqrt(4*1896/3).
```

There are only 36 such integral reduced forms.  Gluing a section to the
nonzero discriminant class of `E7` adds the minuscule-weight correction
`3/2`.  Requiring the resulting `E7+E8` frame to be integral and even leaves
12 candidates.  Comparing exact local genus symbols with the recovered
rank-17 frame leaves exactly three:

```text
H1 = [5/2   1]       H2 = [4       0]       H3 = [21/2  3]
     [  1 190],           [0   237/2],           [   3 46].
```

For all three, the verifier checks:

```text
det(H)             = 474
root system        = E7+E8
root rank/count    = 15 / 366
root determinant   = 2
frame determinant  = 948
discriminant group = Z/948
local genus        = the recovered rank-17 genus.
```

Their exact frames are pinned in
[`data/fibrations/`](data/fibrations/) as
`kumar_e7e8_mw2_frame_1.txt`, `..._2.txt`, and `..._3.txt`.

The integral automorphism groups of the three height lattices have orders

```text
|Aut(H1)| = 2,   |Aut(H2)| = 4,   |Aut(H3)| = 2.
```

Thus `H2` is the unique compatible anchor with an involution other than
`-1`.  The quotient labels now identify that involution exactly.  On Elkies's
curve let

```text
alpha:(t,u) -> (-t,u),
beta: (t,u) -> (-t,-u),
h=alpha*beta:(t,u) -> (t,-u).
```

The two elliptic quotients have minimal models

```text
C/<alpha>: y^2+xy = x^3+x^2+81x-27       (474a1),
C/<beta>:  y^2+xy+y = x^3-7x+14          (474b1).
```

The exact Atkin--Lehner quotient table identifies these as

```text
474a1 = X_0^6(79)/<w6,w79>,
474b1 = X_0^6(79)/<w3,w158>.
```

Their subgroup intersection is `<w474>`, the genus-two quotient defining
`C`.  Consequently

```text
alpha = w6=w79,  beta = w3=w158,  h = w2=w237.
```

The label is also visible intrinsically in the lattice.  The extra
automorphism of

```text
H2 = diag(4,237/2)
```

acts on the cyclic discriminant group `Z/948` as multiplication by `475`:
it is `-1 mod 4` and `+1 mod 3,79`.  This is exactly the `w2` local sign
pattern; composing with global negation gives the equivalent `w237` pattern.
It also matches the polarization criterion: `474=2*237`, and both quaternion
algebras `(-474,2)` and `(-474,237)` are ramified exactly at `2,3`.

Thus the extra symmetry of `H2` is the hyperelliptic Atkin--Lehner involution
`w2=w237`.  The Kumar/abelian-surface moduli therefore factor through the
rational `t`-line.  The square root `u` is still essential: it supplies the
rational descent of the Mordell--Weil sections at a rational point of the
double cover.

The height formula makes this division of labor concrete.  In the basis
`H2=diag(4,237/2)`, the height-4 generator meets the identity component of the
`E7` fiber and has `P.O=0`.  The other generator meets the nonidentity
component, whose local correction is `3/2`, so

```text
237/2 = 4 + 2*(Q.O) - 3/2,   hence Q.O=58.
```

Solving simultaneously for both sections in the Kumar chart is therefore the
wrong-sized equation system.  After composing the Atkin--Lehner action with
fiberwise negation if necessary, the height-4 direction is fixed and the
height-`237/2` direction is negated.  This is consistent with the surface and
the low section descending to `Q(t)`, while `u` chooses the sign/descent of the
level-79 section.

These claims are checked exactly by
[`scripts/deconstruct_x0679_quotients.sage`](scripts/deconstruct_x0679_quotients.sage).
The Atkin--Lehner labels and the alternative equation for this same quotient
are independently recorded in the
[Padurariu--Saia data repository](https://github.com/fsaia/GenusAtMost2)
for their complete low-genus quotient tables.

That independent model is not merely abstractly isomorphic.  At source commit
`6cc368fe37aa67187783118f18d149b2b1fd6230`, the entry
`[6,79,{474}]` is

```text
y^2 = -27*x^6 + 198*x^4 - 171*x^2 + 576.
```

It is related to Elkies's coordinates by the exact birational change

```text
x = 2/t,       y = 6*u/t^3,
t = 2/x,       u = 4*y/(3*x^3).
```

Indeed, clearing denominators gives

```text
t^6*(-27*(2/t)^6 + 198*(2/t)^4 - 171*(2/t)^2 + 576)
  = 36*(16*t^6 - 19*t^4 + 88*t^2 - 48).
```

Thus the non-CM point maps to `(x,y)=(13/7,12048/343)`, the two
points at infinity map to `(0,+/-24)`, and the `t=+/-2` CM orbit maps to
`(x,y)=(+/-1,+/-24)`.  The involutions `(alpha,beta,h)` become
`((-x,-y),(-x,y),(x,-y))`.  This pins the source coordinate and all three
Atkin--Lehner sign conventions independently of a numerical genus-two
isomorphism; the exact replay is part of
[`scripts/deconstruct_x0679_quotients.sage`](scripts/deconstruct_x0679_quotients.sage).

## Correct H21/H92 branch at the CM-24 cusp

The third compatible Kumar polarization must be treated as a different
`E7+E8` fibration, not as an extra section imposed on the `H2` model.  Its
doubled height form is

```text
2*H3 = [21  6]
       [ 6 92].
```

It therefore lies on `H21 intersect H92`; the difference of its two marked
vectors has square `21+92-2*6=101`, so the same branch also lies on `H101`.
The determinant is `1896=4*474`, as required for the level-474 Shimura curve.

At the CM-24 point `(r,s)=(w,-w)`, where `w^2-w+1=0`, the pullback of the
primitive Humbert-21 Satake equation to the Elkies--Kumar `H92` chart has six
distinct tangent components.  In order, their slopes are

```text
(3*w+5)/7, (-21*w+5)/19, w-1,
(11*w-35)/31, (19*w-99)/91, (-7*w-33)/37.
```

An exact calculation modulo 11 strips the chart-base factor `r^88*s^28`
from the degree-700 pullback and factors the degree-584 residual.  The six
components through CM-24 have the following plane degrees:

| branch | tangent slope | degree modulo 11 | role |
|---:|---|---:|---|
| 1 | `(3*w+5)/7` | 25 | extra component |
| 2 | `(-21*w+5)/19` | 30 | `b=42`, level 42 / `H29` branch |
| 3 | `w-1` | 2 | decomposable/boundary branch |
| 4 | `(11*w-35)/31` | 8 | `b=30`, level 258 / `H53` branch |
| 5 | `(19*w-99)/91` | 20 | `b=18`, level 402 / `H77` branch |
| **6** | **`(-7*w-33)/37`** | **21** | **`b=6`, level 474 / `H101` target** |

The branch-6 plane component has arithmetic genus 190 and total delta
invariant 188, hence normalized genus two.  Its normalized point counts are

```text
#C(F_11) = 22,  #C(F_121) = 116,
L_11(T) = 1 + 10*T + 47*T^2 + 110*T^3 + 121*T^4.
```

These two counts agree with the independently published level-474 model

```text
y^2 = -27*x^6 + 198*x^4 - 171*x^2 + 576,
```

whereas the nearby level-402 control has counts `(12,166)`.  This pins the
correct backtrack component to branch 6 at this good prime.  It also explains
why the earlier `H2`/q80 reconstruction was aimed at the wrong polarization.

The exact replay is
[`scripts/verify_h21_h92_level474_branch.sage`](scripts/verify_h21_h92_level474_branch.sage),
with generated result
[`../artifacts/generated-results/elkies-k3-h21-h92-level474-branch-mod11.json`](../artifacts/generated-results/elkies-k3-h21-h92-level474-branch-mod11.json).
It verifies the pinned source hashes, sparse pullback, complete factor-degree
list, tangent labels, target normalization genus, and point counts.

The one-prime calculation first pins which modular factor must be lifted.  The
lift itself is now exact.  The worker
[`scripts/factor_h21_h92_level474_modp.sage`](scripts/factor_h21_h92_level474_modp.sage)
extracts the normalized degree-21 factor at a caller-selected good prime.
Eight images at

```text
17, 19, 23, 29, 31, 41, 43, 47
```

have CRT modulus `553401357731`.  Centering their common 133 coefficients
gives a primitive polynomial `C(r,s)` of total degree 21, with
`deg_r(C)=deg_s(C)=13` and maximum coefficient `26378808832`.  Its leading
`r^13` row is exactly `r^13*(s-1)^8`.

This reconstruction is not promoted merely from an unused-prime match.
[`scripts/reconstruct_h21_h92_level474_qq.sage`](scripts/reconstruct_h21_h92_level474_qq.sage)
proves exact divisibility over `QQ`.  The structural change

```text
r = x/(s-1),
F(x,s) = (s-1)^5*C(x/(s-1),s)
```

makes `F` polynomial and monic of degree 13 in `x`.  Give `x` weight 3 and
`s` weight 2.  Then `F` has weight at most 39.  The Satake coordinate of
weight `k` becomes polynomial after multiplication by `(s-1)^(4*k)` and has
new weight `16*k`; since the Humbert-21 equation is homogeneous of weight
120, its transformed pullback has weight at most 1920.  Reduction by the
monic `F` therefore has `x`-degree at most 12 and coefficient `s`-degree at
most 960.  The verifier evaluates that exact remainder at all 961 rational
values `s=-480,...,480`; every value is zero, proving the remainder polynomial
is identically zero.

The generated characteristic-zero certificate is
[`../artifacts/generated-results/elkies-k3-h21-h92-level474-factor-qq.json`](../artifacts/generated-results/elkies-k3-h21-h92-level474-factor-qq.json).
It records every coefficient and the full finite degree-bound argument.

The plane equation has the exact involution

```text
(r,s) -> (-1/s,-1/r),
r^13*s^13*C(-1/s,-1/r) = -C(r,s).
```

Its invariant quotient is now normalized exactly.  Put

```text
t = r/s,
m = s - 1/r,
a = t*m = r - 1/s,
b = t - a.
```

Eliminating `s` gives a 73-term degree-20 equation in `(t,m)`.  The change to
`a` removes a spurious factor `t^7`, and the further change to `b` produces a
37-term model of bidegree `(5,8)` and total degree 13.  Its two infinity
singularities have multiplicities eight and five.  Together with the rational
node `(-1,46)`, they define an exact quadratic Cremona transformation whose
strict transform has degree 11.

On that reduced curve, Singular's local ideal-quotient calculation returns a
ten-dimensional degree-nine adjoint space.  Saturating its integer coefficient
lattice and applying LLL lowers the maximum coefficient size from 154 bits to
29 bits without changing its `QQ`-span.  There is then no need for the slow
generic rational-normal-curve inversion.  At the smooth rational point
`[-54:1:-54]`, the adjoint jet matrix has rank ten; the spaces of sections
vanishing to orders at least eight and nine have dimensions two and one.
Their ratio is therefore a degree-one parameter on the normalization.  A
32-term Newton expansion and degree-11 Padé reconstruction give the inverse,
and exact substitution verifies both the reduced equation and the parameter.

The double cover is intrinsic.  Since

```text
Y = r + 1/s,
Y^2 = a^2 + 4*t,
```

its squarefree branch polynomial in the recovered parameter has degree six.
Writing that parameter as `x0`, the exact base change

```text
x = 164590478411323*x0
    /(18542755847723*x0 + 52351461697088)
```

together with the recorded rational multiplier on `Y`, satisfies identically

```text
y^2 = -27*x^6 + 198*x^4 - 171*x^2 + 576.
```

Thus the characteristic-zero normalization and the birational identification
with the published level-474 curve are proved, not inferred from matching
point counts.  The exact replay is
[`scripts/normalize_h21_h92_level474_qq.sage`](scripts/normalize_h21_h92_level474_qq.sage),
and its generated certificate is
[`../artifacts/generated-results/elkies-k3-h21-h92-level474-normalization.json`](../artifacts/generated-results/elkies-k3-h21-h92-level474-normalization.json).

### Rational H21 entrance at the non-CM point

The published rational point `(x,y)=(13/7,12048/343)` maps exactly to

```text
(r92,s92)=(-3621005/690947, 158286/143585)
```

on the H92 chart.  Weighted Igusa matching gives the rational H21
presentation

```text
r21 = 271621954946208883/51863976688786080,
s21 = 33855626850015548642165023481946578233811
      /37279290235251051531892004373889536000.
```

All four weighted Igusa identities and all five short-Weierstrass coefficient
identities hold exactly.  More strongly, the resulting twist parameter is a
rational square, so the two unmarked models are isomorphic over `QQ`, rather
than only after a quadratic extension.  Replaying the pinned Elkies--Kumar
`21/21.txt` construction produces a ten-term plane cubic of bidegree `(3,3)`
over `QQ(u)` with an exact rational nonflex point at old base coordinate zero.
Consequently the H21 entrance genus-one curve has a rational origin and is
its own Jacobian over `QQ(u)` at the target non-CM point.

The oriented Hilbert-cover coordinates are separate.  At these rational H21
and H92 points, their squares satisfy

```text
z21^2 = -52203427 * (rational square),
z92^2 = -52203427 * (rational square).
```

The ratio `z21/z92` is rational, so both orientations use the same field
`QQ(sqrt(-52203427))`; no biquadratic marking field appears.  This is strong
exact compatibility evidence for the H21/H92 intersection.  It does not by
itself determine the field of definition of the K3 Mordell--Weil section: an
oriented Hilbert cover and a rational K3 section need not have the same
descent field.  The independent lattice and coordinate certificates below
settle that section question.

This statement has a deliberate boundary.  The cubic in `21/21.txt` is the
`A2+A6+E8 -> E7+E8` entrance to H21.  It has not been identified with the
oppositely directed `E7+E8 -> E8+E6` q=6 pencil below.  A separate exact
integral-lattice argument now certifies the section's field of definition;
the explicit divisor function remains the equation-level step.

The exact checker is
[`scripts/verify_h3_noncm_q6_source_anchor.sage`](scripts/verify_h3_noncm_q6_source_anchor.sage),
with generated certificate
[`../artifacts/generated-results/elkies-k3-h3-noncm-q6-source-anchor.json`](../artifacts/generated-results/elkies-k3-h3-noncm-q6-source-anchor.json)
(SHA-256 `0560b1921c87ad2d8db6c293ce070cb30aa75626315801e3e4a71cad59573ea5`).

The marked section descent is unconditional.  In the pinned ancillary model,
the starting `A2+A6+E8` fibers and 3-neighbor parameter are defined over
`QQ(r,s)`.  Their trivial lattice has determinant `3*7=21`, hence is already
integrally saturated because 21 is squarefree.  The transported `E7+E8` MW
generator is consequently individually `QQ(r,s)`-rational and has height
`21/2`.  Thus `P1`, `-P1`, and `D=O+(-P1)-F` are Galois fixed at the rational
H21/H92 point; the common oriented-cover field
`QQ(sqrt(-52203427))` does not exchange the signs.  Replay with
[`scripts/verify_h21_q6_section_descent.sage`](scripts/verify_h21_q6_section_descent.sage);
the artifact is
[`../artifacts/generated-results/elkies-k3-h21-q6-section-descent.json`](../artifacts/generated-results/elkies-k3-h21-q6-section-descent.json)
(SHA-256 `9ccdfc7b7a1ca79d549c161e9922051e9f90d7b89ddc81057e17188eedc2a4d2`).

The second H3 direction is independently rational and saturated.  The pinned
H92 ancillary construction begins with split `D6+A8+A1` fibers and an
explicit section whose local corrections are `1/2` and `20/9`, hence whose
height is `23/18`.  The resulting rank-18 intersection lattice has
determinant `72*(23/18)=92` and Smith factors `(2,46)`.  Exhaustion of all 92
discriminant classes finds no nonzero isotropic class, so it admits no proper
even overlattice.  The rational two- then three-neighbor sequence therefore
transports this saturated source lattice to the `E7+E8` fibration over
`QQ(r,s)`, where its rank-one quotient generator has height `46`.  On the
H21/H92 intersection this is `P2` in the H3 Gram
`[[21/2,3],[3,46]]`; it is individually rational, not merely rational modulo
the span of `P1`.

Replay this with
[`scripts/verify_h92_section_descent.sage`](scripts/verify_h92_section_descent.sage).
The generated certificate is
[`../artifacts/generated-results/elkies-k3-h92-section-descent.json`](../artifacts/generated-results/elkies-k3-h92-section-descent.json)
(SHA-256 `fe525f75fa87c31afb34755fe63fc778349d2843010eb5c9b17ce6d8b8712e40`).
The proof boundary is coordinate-level: explicit height-`46` coordinates on
the final short H92 model have not been recovered.  They are unnecessary for
the first q=6 chord, which uses `P1`, but may be required to track later
sections through the equation chain.

### Pole-profile reduction for the second direction

The marked H3 height matrix is

```text
H = [[21/2, 3], [3, 46]].
```

Here `P1` has the nontrivial `E7` component class (correction `3/2`), while
`P2` has the identity class.  Thus for `Q_m=P2+m*P1`, its correction is
`3/2` precisely when `m` is odd, and Shioda's formula gives

```text
Q_m.O = (height(Q_m) - 4 + correction(Q_m))/2.
```

The exact small-translate table has its unique minimum at `m=0`:

| `m` | `height(Q_m)` | correction | `Q_m.O` |
|---:|---:|---:|---:|
| -2 | 76 | 0 | 36 |
| -1 | 101/2 | 3/2 | 24 |
| 0 | 46 | 0 | **21** |
| 1 | 125/2 | 3/2 | 30 |
| 2 | 100 | 0 | 48 |

Indeed, relative to `P2` the pole increase is
`3*m*(7*m+4)/4` for even `m`, and
`(21*m^2+12*m+3)/4` for odd `m`; the latter has negative discriminant and the
former is positive at every nonzero even integer.  Therefore a small
unimodular basis change cannot lower the denominator degree: the correct
target remains `P2` itself, with pole order 21.  The exact replay is
[`scripts/analyze_h3_p2_pole_profile.sage`](scripts/analyze_h3_p2_pole_profile.sage).

The source involutions do not yet come with a coordinate-level action on this
marked section.  Since the lattice/descent certificate already makes `P2`
individually rational over the H92 function field, a quadratic trace/norm
ansatz would be redundant unless an explicit involution sends it to a distinct
section.  The next economical equation-level step is consequently the marked
degree-29 divisor transport, not an unrestricted eigenspace interpolation.

The smaller first-neighbor marked point has now been lifted exactly before
that final divisor calculation.  The source degree-three divisor gives a
Jacobian point with coordinate degrees `(22,18)` and `(33,27)`.  A fast worker
uses only split marked fibers, rather than the broader diagnostics in the
transport notebook; 80 additional 61-bit good-prime records raise the CRT
modulus to 9,147 bits.  Every coefficient of both intermediate coordinates
then rationally reconstructs, and direct substitution in the exact
first-neighbor Weierstrass equation holds.  Reproduce this stage with

```bash
H92P2_PRIME=2305843009213693967 \
H92P2_MODULAR_OUTPUT=artifacts/generated-results/h92-p2-modular/intermediate-2305843009213693967.json \
sage elkies-k3/scripts/sample_h92_p2_intermediate_modp.sage
sage elkies-k3/scripts/crt_h92_p2_intermediate.sage
H92P2_PROBE=0 sage elkies-k3/scripts/probe_h92_p2_final_divisor.sage
```

The resulting intermediate CRT record has SHA-256
`dfe18860793ad2fae6012716e07d5289d94a74aea27401aafdb4b9ef3c0ac60a`.

The final transport is now exact.  For the degree-29 divisor `D`, plane class
`H`, and degree-10 residual point `R` with `D+R ~ 10H`, the intrinsic
degree-zero class is `3D-29H = H-3R`, not the nonzero-degree expression
`D-7H`.  Modular interpolation of this canonical class gives the doubled
section with `x` degrees `(184,180)`.  The degree-two division polynomial has
a unique rational linear factor whose half has precisely the lattice-predicted
profile

```text
x(P2): (46,42),    y(P2): (69,63).
```

Rather than CRT-lifting the doubled coordinates, the normalized half solves a
138-variable cleared Weierstrass system in `X,Y,Z`, where
`x=X/Z^2`, `y=Y/Z^3`, and `deg(Z,X,Y)=(21,46,69)`.  Its Jacobian has full rank
138 modulo `100003`; p-adic Newton lifting to 1,024 digits rationally
reconstructs every coefficient.  Exact substitution verifies the H92
Weierstrass identity, and the doubled section agrees with the canonical
degree-29 class on 100 fibers modulo `100003`.  The resulting certificate is
[`../artifacts/generated-results/elkies-k3-h92-p2-lift.json`](../artifacts/generated-results/elkies-k3-h92-p2-lift.json),
with status `PASS_EXACT_H92_P2` and SHA-256
`fd3eccccd94ba59e0e0472ac13f635c209b4fbfa3df6fe6ff6a7847a38c365e4`.
The y-sign is fixed by the retained modular half; replacing it gives the
opposite, equally valid MW generator.

The defining H3 family itself is now explicit.  Let `(X,Y)` lie on

```text
Y^2=-27*X^6+198*X^4-171*X^2+576.
```

The exact normalization supplies rational functions `t(x)`, `a(x)` and the
multiplier relating its double-cover coordinate to the published `Y`.
After inverting the linear-fractional `X`-map and writing `Y0=Y/m(x)`, the H92
chart is

```text
r=(a+Y0)/2,   s=2/(Y0-a).
```

Direct function-field substitution proves the degree-21 H21/H92 component
equation.  Evaluating the pinned H92 coefficients `A1,A,B1,B,B2` at these
functions gives the exact source K3

```text
v^2=u^3+(A1*tau^3+A*tau^4)*u+(B1*tau^5+B*tau^6+B2*tau^7).
```

Thus the H3 source is an equation over the function field of the published
genus-two level-474 curve, with `E7` at `tau=0`, `E8` at infinity, and the two
rational MW directions above.  Exact specialization at
`(X,Y)=(13/7,12048/343)` recovers the pinned non-CM H92 point and all five
short-Weierstrass coefficients.  Replay with
[`scripts/export_h3_level474_source_family.sage`](scripts/export_h3_level474_source_family.sage);
the generated certificate is
[`../artifacts/generated-results/elkies-k3-h3-level474-source-family.json`](../artifacts/generated-results/elkies-k3-h3-level474-source-family.json)
(SHA-256 `8f5afd11e1d8979d57cb1a569833309f9664c19cd47194af0581a5cbbf8f1d59`).
This closes the defining-source-family gate only: the q=6 and later neighbors
still have to be executed before the rootless MW17 equation exists.
<!-- status-consumer: EC-K3-H3-SOURCE a4bb40c9c9d0ff09 -->

The simplest alternative-specialization search can also be reduced exactly
before the neighbor chain is available.  Quotienting the even genus-two model
by `X -> -X`, with `xE=-27*X^2` and `yE=-27*Y`, gives

```text
yE^2=xE^3+198*xE^2+4617*xE+419904.
```

The base is now determined globally by the exact Magma certificate
[`scripts/prove_h3_level474_rational_points.m`](scripts/prove_h3_level474_rational_points.m).
Its two-cover descent has three locally soluble classes; their quartic-factor
elliptic quotients give two cubic-field covers.  Elliptic Chabauty at `p=41`
has upper bounds equal to the found cover points and residual index `1`, giving
the complete rational `X`-set

```text
0, -1, 1, -13/7, 13/7.
```

There are no rational points at infinity, and exact substitution gives both
signs of `Y` above each value.  Hence the ten points are precisely
`(0,+-24)`, `(+-1,+-24)`, and `(+-13/7,+-12048/343)`.  The old local quotient
sieve remains a bounded cross-check only; the canonical proof and output are
recorded in
[`../elliptic-curves/notes/ICARM_CURVE273_CONSTRUCTION_INVESTIGATION.md`](../elliptic-curves/notes/ICARM_CURVE273_CONSTRUCTION_INVESTIGATION.md)
and
[`../artifacts/generated-results/elkies-k3-h3-level474-rational-points.txt`](../artifacts/generated-results/elkies-k3-h3-level474-rational-points.txt).
<!-- status-consumer: EC-K3-H3-PTS 8f0a27c947843b4a -->

The section is also explicit, not merely an integral-lattice class.  Modular
nonflex conversion of the pinned H21 entrance cubic, transported to the
smaller H92 short model, determines its `x`-coordinate with entrance-base
degrees `(10,12)`.  Nine disjoint windows retain 204 good primes and give a
1,945-bit CRT modulus.  Structured simultaneous LLL exploits
`D(u)=d4*u^4*Z4(u)^2`, `Z4(0)=1`; the recovered `y`-coordinate has degrees
`(15,18)`, direct substitution proves the H92 Weierstrass identity over
`QQ`, and an exact marked-fiber incidence fixes the square-root sign.  Replay
the exact lift from the pinned windows with

```bash
sage -python elkies-k3/scripts/lift_h21_p1_modular.sage \
  artifacts/generated-results/elkies-k3-h92-p1-mod-window-100.json \
  artifacts/generated-results/elkies-k3-h92-p1-mod-window-199.json \
  artifacts/generated-results/elkies-k3-h92-p1-mod-window-383.json \
  artifacts/generated-results/elkies-k3-h92-p1-mod-window-557.json \
  artifacts/generated-results/elkies-k3-h92-p1-mod-window-733.json \
  artifacts/generated-results/elkies-k3-h92-p1-mod-window-887.json \
  artifacts/generated-results/elkies-k3-h92-p1-mod-window-1069.json \
  artifacts/generated-results/elkies-k3-h92-p1-mod-window-1300.json \
  artifacts/generated-results/elkies-k3-h92-p1-mod-window-1600.json \
  --output artifacts/generated-results/elkies-k3-h92-p1-lift.json
```

The certificate is
[`../artifacts/generated-results/elkies-k3-h92-p1-lift.json`](../artifacts/generated-results/elkies-k3-h92-p1-lift.json),
with 204 primes, 1,945 CRT bits, status `PASS_EXACT_H92_P1`, and SHA-256
`c323bf6346bb239934a5a2d8b1a3f4067e70e993d2e4eb32aaa30f469fca6397`.
It supersedes historical SHA-256
`0602c3b199629c6f460c9b7c728e048822418ecf85bf54807852be3d97b66616`,
which lacked only the exact orientation-incidence block and retained the old
H21 status label.  The proof boundary is exact and narrow: the marked
height-`21/2` H92 section is certified, while the chord and explicit basis of
`H0(O+(-P1)-F)` remain to be constructed before the q=6 pencil exists at the
equation level.

There is now also an exact lattice-level gate for the tempting descent by
pairing the two signed q=6 pencils.  The inverse section is not obtained by
only changing the MW sign while retaining the E7 correction: the two exact
section classes have different root corrections.  For
`D-=O+(-P1)-F` and `D+=O+(+P1)-F`, the checker proves

```text
D-^2 = D+^2 = 0,  D-.F = D+.F = 2,  D-.D+ = 21.
```

Thus `H=D-+D+` is primitive with `H^2=42`, `H.F=4`, `h0(H)=23` if nef,
and arithmetic genus 22.  Conditional on the still-unproved equation-level
Galois rule exchanging the signs, the natural tensor module has rank four
and the unordered trace/norm module rank three; neither is an elliptic
pencil.  The preceding certificate proves that this hypothetical sign
exchange is not the actual H21 marking: both signed sections are fixed.  The
calculation still rejects that shortcut in any sign-exchanged setting and
does not rule out an accidental elliptic factor of the genus-22 Jacobian.  Replay
with
[`scripts/verify_h3_q6_signed_descent_gate.sage`](scripts/verify_h3_q6_signed_descent_gate.sage);
the artifact is
[`../artifacts/generated-results/elkies-k3-h3-q6-signed-descent-gate.json`](../artifacts/generated-results/elkies-k3-h3-q6-signed-descent-gate.json)
(SHA-256 `3be2de6e2f7c722bc04dde0ad5eba81924b130b93d6850009cd266398b4b60d7`).

## Preferred first transport from H3

The first neighbor out of the now labeled `H3` family need not use the much
larger q80 class.  Start from
[`data/fibrations/kumar_e7e8_mw2_frame_3.txt`](data/fibrations/kumar_e7e8_mw2_frame_3.txt)
and enumerate the complete q=6 shell subject to the section marking
`coordinate 15 = 1`.  Taking only the proper factor presentation `(a,b)=(2,3)`
and only `a<=b`, the exact Fincke--Pohst traversal has 441 nodes and returns 56
sign-pair vectors.  The traversal reports `exhaustive=true`, and all 56
vectors give

```text
root rank = 14,  number of roots = 312,  root determinant = 3,
MW rank = 3.
```

The ADE root lattice with these invariants is `E8+E6`: its two factors have
240 and 72 roots and determinant `1*3`.  Thus this constrained q=6 shell gives
an exact, uniform transport

```text
H3 (E8+E7, MW 2)  --q=6, (a,b)=(2,3)-->  E8+E6, MW 3.
```

A compact witness, in the coordinates of the source frame, is

```text
v = (0,0,-1,-1,-1,-1,-1,0,0,0,0,0,0,0,0,1,0).
```

The four frames archived under
`artifacts/local/elkies-k3/h3-q6-mw1-frames/` are merely the representatives
retained by `--per-root-data-cap 4`.  There is no separate isometry or Weyl
orbit certificate proving that the 56 vectors have only four orbits.  The
proved statement is instead stronger in the direction needed here: every one
of the 56 vectors in the specified constrained shell has the same root data
and MW rank.

For the displayed compact witness, the saturated reduced MW height Gram is

```text
[8/3  1/3  -1]
[1/3  8/3   1]
[ -1    1  46].
```

An exhaustive short-vector calculation proves that the optimal section-pole
profile is `(0,0,21)`: two height-`8/3` directions are polynomial, while every
integral MW basis needs a third direction with pole at least 21.  This child is
not the older H2 q=60 `E8+E6/MW3` frame.  That frame has reduced height Gram
`[[4,0,0],[0,20/3,1],[0,1,12]]`, and exact `qfisom` on the two scaled ternary
forms returns no isometry.  Thus the shared root type does not provide a q=60
shortcut; equation-level continuation should exploit the clean degree-two
source chord without trying to realize the high-pole third child generator.
The exact chamber, height, pole, and non-isometry replay is

```bash
sage elkies-k3/scripts/analyze_h3_first_q6_chamber.sage
```

This is the preferred first transport from the normalized `H3` family.  Small
bounded q=4 and q=6 searches from a root-adapted copy of the first child found
only MW-rank-three frames (including root data `(14,184,16)` and returns to
`(14,312,3)`).  Those second-step searches are experiments, not exhaustive
obstructions.  They point instead to the horizontal q=8 shell treated next.

## Rank-growing q=8 second transport

Blind enumeration is especially misleading here.  The exact norm-16 shell
of the q=6 child contains 219,758,670 signed vectors; 139,006,800 of them are
outside the `E8+E6` root span.  The primitive root lattice nevertheless makes
the quotient classification small.  Modulo `W(E6) x W(E8)`, the projected MW
norm bound leaves ten nonzero sign representatives.  Combining these with
the 45 and 12 bounded dominant weights of the two root components and imposing
the exact discriminant congruence gives only 63 horizontal q=8 orbits.

Of these 63 orbits, two are nonprimitive and 61 give primitive U-neighbors.
Exactly two have child root data

```text
root rank = 13,  number of roots = 312,  root determinant = 4,
MW rank = 4.
```

The connected simply-laced rank-13 root system with 312 roots and determinant
four is `D13`.  In the pinned simple-root/MW coordinates, one dominant witness
has

```text
MW projection = (0,-2,0),
dominant labels = (0,0,0,0,0,0,0,0,0,1,1,0,0),
v = (-20,-10,2,4,6,8,5,10,7,4,-4,-9,-14,1,0,-2,0).
```

This is not only a lattice neighbor.  In the deterministic q=6 child
coordinates, the equivalent chamber witness is

```text
v = (-5,-4,-3,4,5,7,10,8,6,4,2,-4,2,-2,-4,0,-2).
```

The raw `(a,b)=(2,4)` class has old-fiber degree four and intersection `-2`
with the old zero section.  One zero-section reflection gives the `(4,2)`
presentation, with old-fiber degree two and zero-section intersection two.
Its finite-component pairings are

```text
(0,0,1,0,0,0,0,0,0,0,0,1,0,0),
```

and its affine `E6,E8` pairings are `(1,0)`.  The MW projection is `-2e2`.
The height minimum `8/3`, together with the nonintegrality of `v/2`, excludes
a negative section; the norm identity for a degree-two curve would force
`w=v` and simultaneously `v^2=4k+2`, contradicting `v^2=16`.  This proves
nefness rather than inferring it from the child root data.

The exact chain is therefore

```text
H3: E8+E7/MW2  --q=6, degree 2-->  E8+E6/MW3
                  --q=8, degree 2-->  D13/MW4.
```

[`scripts/analyze_h3_first_q6_chamber.sage`](scripts/analyze_h3_first_q6_chamber.sage)
checks both chamber reductions and the nef proof.  The exhaustive Weyl-orbit
classification is
[`scripts/classify_h3_q6_child_q8_orbits.sage`](scripts/classify_h3_q6_child_q8_orbits.sage),
with generated basis/transport certificate
[`../artifacts/generated-results/elkies-k3-h3-q6-q8-orbits.json`](../artifacts/generated-results/elkies-k3-h3-q6-q8-orbits.json).
For continuation searches it pins the simple-root/LLL-quotient frame
[`data/fibrations/h3_q6_q8_d13_mw4_root_adapted_frame.txt`](data/fibrations/h3_q6_q8_d13_mw4_root_adapted_frame.txt),
whose rank-four MW height Gram is

```text
[ 3/4   1/4  -1/4   0]
[ 1/4  11/4   1/4   1]
[-1/4   1/4  11/4  -1]
[   0     1    -1  46].
```

The remaining equation-level gate is to realize these two compensated
degree-two chords on the characteristic-zero H3 family and continue from the
pinned D13/MW4 frame toward the rootless MW-rank-17 fibration.

A deterministic q=4 witness from the pinned D13 frame is

```text
(a,b) = (2,2),
v = (3,-1,1,4,8,5,7,6,5,4,-1,-3,3,-1,1,-1,0).
```

It is already in the effective chamber: `D.F=2`, `D.O=0`, the only nonzero
D13 simple-component pairing is one at node six, and the affine-component
pairing is one.  Full nefness is exact.  In the rank-four MW quotient only
eight shifted cosets can have projected distance below two; exact D13 CVP
gives total squared distances `2,2,3,3,3,3,3,3`, so no section is negative.
The standard bisection norm identity and parity exclude a negative degree-two
curve.  The child has root data `(13,158,26)`, hence `A12+A1/MW4`.  Therefore
this is a genuine degree-two lateral presentation, but it does not raise the
MW rank.  The replay is
[`scripts/analyze_h3_d13_q4_chamber.sage`](scripts/analyze_h3_d13_q4_chamber.sage).

The exhaustive Weyl-quotient barrier from D13 now covers every proper
presentation through q=20.  The q=20 factor orders `(a,b)=(10,2)` and `(5,4)`
each have 1,567 dominant orbits.  The former has 1,533 primitive neighbors
and the latter 1,567; both have maximum MW rank four.  The exact endpoint
certificates are
[`../artifacts/generated-results/elkies-k3-h3-q6-q8-d13-q20-degree2.json`](../artifacts/generated-results/elkies-k3-h3-q6-q8-d13-q20-degree2.json)
and
[`../artifacts/generated-results/elkies-k3-h3-q6-q8-d13-q20-degree4.json`](../artifacts/generated-results/elkies-k3-h3-q6-q8-d13-q20-degree4.json).
This is a bounded no-growth theorem for q at most 20, not a global barrier.
The q=21 degree-three and q=22 degree-two proper presentations likewise have
maximum MW rank four; q=23 has no proper factor presentation.

The first rank growth occurs at q=24 and, importantly, already has degree two.
The exact Weyl quotient for `(a,b)=(12,2)` has 2,709 dominant orbits, of which
2,653 are primitive.  Exactly three have root data

```text
(root rank, root count, root determinant) = (12,264,4),
```

so their child is `D12/MW5`.  The preferred orbit 85 has MW projection
`(0,-1,1,1)`, its sole nonzero D13 label is node three, and its witness is

```text
v = (0,5,0,1,2,1,2,2,2,2,4,8,2,0,-1,1,1).
```

The raw `(12,2,v)` class is already in the effective chamber, with `D.F=2`,
`D.O=10`, simple-component pairings `(0,0,1,0,...,0)`, and affine pairing
one.  The complete shifted rank-four MW quotient ball of squared radius two
is empty, excluding a negative section before any root completion.  The
degree-two bisection norm/parity identity excludes the remaining horizontal
case.  Hence the class is nef and its geometric pencil again has degree two.
The classification artifact is
[`../artifacts/generated-results/elkies-k3-h3-q6-q8-d13-q24-degree2.json`](../artifacts/generated-results/elkies-k3-h3-q6-q8-d13-q24-degree2.json)
(SHA-256 `66d5a7ff6ec26f8aa8344cdbd779a6c96707b041ba4f89d7dbfe460c95485a93`),
and [`scripts/analyze_h3_d13_q4_chamber.sage`](scripts/analyze_h3_d13_q4_chamber.sage)
certifies the chamber.  The exact low-degree chain is now

```text
H3 MW2 --q6/degree2--> MW3 --q8/degree2--> MW4
       --q24/degree2--> D12/MW5.
```

## Degree-two continuation to MW rank eight

The root-adapted quotient search continues rank growth at every next chosen
step:

```text
D12/MW5 --q6--> A11/MW6 --q8--> A5+A5/MW7
          --q4--> 3A3/MW8.
```

All three arrows again have old-fiber degree two.  The selected D12 q6 orbit
42 has witness

```text
(1,4,-1,0,7,4,3,7,7,7,7,7,-1,0,-1,-1,0).
```

Its component labels are at D12 nodes three and six, its affine pairing is
zero, and the eight potentially dangerous shifted MW cosets all have exact
total squared distance three.  The selected A11 q8 orbit 922 has witness

```text
(0,-1,-5,2,-2,-5,-4,-3,-3,-3,-3,-1,0,0,1,0,-1).
```

Its shifted MW quotient ball below squared radius two is empty.  Finally, the
selected A5+A5 q4 orbit 472 has the compact witness

```text
(1,1,1,1,1,1,1,1,1,1,0,0,-1,0,0,0,1).
```

The corresponding exact section distances are `2,2,3,3,3,3,3,3`.  In all
three cases the degree-two norm identity and parity exclude a negative
bisection.  Thus these are nef geometric pencils, not merely abstract child
frames.  The exact replay is
[`scripts/analyze_h3_rank_growing_degree2_chain.sage`](scripts/analyze_h3_rank_growing_degree2_chain.sage).

The selected search artifacts are
[`../artifacts/generated-results/elkies-k3-h3-d12-o85-q6-degree2.json`](../artifacts/generated-results/elkies-k3-h3-d12-o85-q6-degree2.json),
[`../artifacts/generated-results/elkies-k3-h3-a11-middle-q8-degree2.json`](../artifacts/generated-results/elkies-k3-h3-a11-middle-q8-degree2.json),
and
[`../artifacts/generated-results/elkies-k3-h3-a5a5-c2-q4-degree2.json`](../artifacts/generated-results/elkies-k3-h3-a5a5-c2-q4-degree2.json).
The resulting root-adapted `3A3/MW8` frame is
[`../artifacts/generated-results/elkies-k3-h3-a5a5-c2-q4-degree2-frames/q4-o0472-r9-n36-d64-4841c34fa442.txt`](../artifacts/generated-results/elkies-k3-h3-a5a5-c2-q4-degree2-frames/q4-o0472-r9-n36-d64-4841c34fa442.txt).

Because the stable Neron--Severi lattices have signature `(1,18)`, rank 19,
and cyclic discriminant group of length one, Nikulin's rank-versus-length
uniqueness theorem identifies each `U+(-frame)` with the recovered target NS
inside this genus.

### Exact lattice suffix to rootless MW rank seventeen

The selected `3A3/MW8` frame has another productive q=4 degree-two shell.
Its Weyl quotient has 2,481 dominant orbits, all primitive.  Exactly 127
reduce the root rank below nine, and the maximum child MW rank is ten.  A
compact selected orbit is 323, with

```text
MW projection = (-1,0,0,1,0,1,0,-1),
dominant labels = (1,0,1,0,0,0,1,1,0),
v = (-3,0,2,2,3,1,2,3,1,-1,0,0,1,0,1,0,-1).
```

Its child root data are `(7,24,36)`, identifying `A3+2A2`, and hence its
MW rank is ten.  Continuing from that single root-adapted frame, deterministic
degree-two first hits give

```text
D13/MW4 --q24--> D12/MW5 --q6--> A11/MW6
         --q8--> 2A5/MW7 --q4--> 3A3/MW8
         --q4--> A3+2A2/MW10 --q4--> 5A1/MW12
         --q4--> 4A1/MW13 --q4--> 3A1/MW14
         --q4--> 2A1/MW15 --q4--> A1/MW16
         --q6--> rootless/MW17.
```

[`scripts/verify_h3_d13_to_mw17_path.sage`](scripts/verify_h3_d13_to_mw17_path.sage)
independently reconstructs all eleven primitive U-neighbors from their
witnesses, verifies every raw and root-adapted child frame, and checks a
determinant-one composite transport.  Its transport SHA-256 is
`6542f74b2780b4143999e346d519bb72690fac2eeeb99293a97192f305d24c40`;
the generated replay is
[`../artifacts/generated-results/elkies-k3-h3-d13-to-mw17-path.json`](../artifacts/generated-results/elkies-k3-h3-d13-to-mw17-path.json),
with SHA-256
`f6eac2339c86de84b79a0ddfec3229df9b9c1617110bdd9c474443e7e39fd484`.

All eleven arrows in this suffix are certified degree-two nef pencils.  For
the q=4 orbit into `A3+2A2/MW10`, the source component pairings are
`(1,0,1,0,0,0,1,1,0)`, the three
affine pairings are `(0,1,1)`, and the complete shifted root/MW closest-vector
profile is two distances two, sixteen distances three, and two distances
four.  The bisection norm/parity identity excludes the remaining horizontal
case.  The exact chamber replay is
[`scripts/analyze_h3_rank_growing_degree2_chain.sage`](scripts/analyze_h3_rank_growing_degree2_chain.sage).
The six subsequent raw divisors are also already in their respective
chambers, with no reflections.  Exact full-frame closest-vector calculations
and bisection parity prove nefness in
[`scripts/analyze_h3_mw10_to_rootless_chambers.sage`](scripts/analyze_h3_mw10_to_rootless_chambers.sage),
whose terminal status is `PASS_H3_MW10_TO_ROOTLESS_NEF`.  Equation-level
divisor functions remain to be constructed.
<!-- status-consumer: EC-K3-H3-D13-MW17-LATTICE-CHAIN 30c9f060a5da7ed5 -->

At the final `A1/MW16` marking, q=4 was tested only as a bounded obstruction:
the exact 9,000-orbit prefix and seven 1,000-orbit strata across the 160,308
dominant-orbit quotient contained no rootless child.  The q=6 degree-two
search therefore used a PARI storage cap of 10,000 MW quotient vectors and
streamed exact children.  Witness 2,247 is rootless, with

```text
dominant A1 label = (1),
MW projection = (-5,-1,2,3,-2,1,1,-1,0,2,1,0,-1,1,-1,3),
v = (-2,-5,-1,2,3,-2,1,1,-1,0,2,1,0,-1,1,-1,3).
```

The full q=6 shell has 7,187,438 vectors according to PARI, so the search is
explicitly bounded; the selected witness, integral neighbor, and rootless
child are nevertheless exact.

## Intrinsic recovery inside the rootless frame

The Humbert split can also be read backwards from the pinned rootless
rank-17 lattice, without using the q80 neighbor chain.  In the coordinates
`U+(-rank17_gram)`, set

```text
h = (4,4,-1,0,-3,0,2,-2,1,-2,1,1,0,1,0,0,-2,-2,2).
```

Exact integral arithmetic gives

```text
h^2 = -4,                    div(h) = 4,
disc(reflection in h) = 475 mod 948,
h^perp: rank 18, det -237, Smith (1^17,237).
```

The multiplier `475` is `-1 mod 4` and `+1 mod 3,79`, exactly the `w2`
local sign pattern found independently in the Kumar frame.  Moreover
`div(h)=|h^2|` makes the splitting integral:

```text
NS = <h> orthogonal_sum h^perp.
```

Deleting the orthogonal height-four row from the Kumar `H2` frame gives
`U+(-H2_fixed)`.  It has the same indefinite genus as `h^perp`; both have
rank 18 and cyclic discriminant group of length one.  Nikulin uniqueness
therefore identifies the two complements.  This is an independent intrinsic
certificate that the correct rootless backtrack is the height-four
`H8` direction plus its determinant-237 fixed complement, rather than an
arbitrary low-MW neighbor.

There is also an exact warning about the recurring `q=8` signal.  Relative to
the pinned rootless `U`, no isotropic vector `(a,b,v)` with `ab=8` is
orthogonal to this `h`.  Orthogonality and isotropy would require

```text
v*R*v = 16,       (x*R/4).v = a+b,
```

where `x` is the last 17 coordinates of `h` and
`(x*R/4)*R^-1*(R*x/4)=9/4`.  For the factor sums `a+b=9,6`,
Cauchy--Schwarz either exceeds norm 16 or forces the nonintegral equality
vector `v=2*x/3`.  Thus the generic and CM `q=8` occurrences discussed below
are ambient alternate fibrations, return loops, or boundary factorizations;
they are not the direct `H2` fiber in this short rootless presentation.

Replay all of these checks with

```text
sage elkies-k3/scripts/verify_rank17_h8_split.sage
```

The verifier proves the splitting and the genus identification.  It does not
yet pin a small coordinate isometry from `h^perp` to the Kumar fixed frame;
the existing q80-to-rootless transport remains the explicit basis-map
certificate.

## Exact Humbert-8 ambient family

The height-4 direction identifies the correct ambient surface, not merely a
convenient section.  Keeping only that identity-component section in the
`E7+E8` frame gives determinant

```text
det(E7)*height(P) = 2*4 = 8.
```

It is therefore the Humbert discriminant-8 locus.  Conversely, keeping only
the nonidentity height-`237/2` section gives determinant

```text
2*(237/2) = 237 = 3*79.
```

Thus the rank-19 Shimura curve is the intersection of the discriminant-8 and
discriminant-237 Noether--Lefschetz/Humbert conditions; together their height
Gram has determinant 474 and the full K3 Neron--Severi determinant is 948.

Elkies--Kumar's ancillary calculation for discriminant 8 gives the whole first
condition explicitly.  To avoid confusing coordinates, let `T` denote the
elliptic base and `(r,s)` the Humbert-surface coordinates.  The Kumar equation
is

```text
Y^2 = X^3 + (A1*T^3 + A*T^4)*X
            + B1*T^5 + B*T^6 + B2*T^7,

A1 = 2*r*s^2,
A  = -(9*r*s + 4*r^2 + 4*r + 1)/3,
B1 = r*s^2*(3*s + 8*r - 2)/3,
B  = -(54*r^2*s + 81*r*s - 16*r^3
       - 24*r^2 - 12*r - 2)/27,
B2 = r^2.
```

The inverse Kumar map simplifies to

```text
I2  = -4*(8*r + 3*s - 2),
I4  =  4*(4*r^2 + 9*r*s + 4*r + 1),
I6  = -4*(48*r^3 + 94*r^2*s + 36*r*s^2
          + 40*r^2 - 35*r*s + 4*r + 4*s - 2),
I10 = -8*r^3*s^2.
```

The oriented Hilbert modular double cover is

```text
z^2 = 2*(16*r*s^2 + 32*r^2*s - 40*r*s - s
         + 16*r^3 + 24*r^2 + 12*r + 2).
```

The same ancillary file supplies a birational rationalization of this
*oriented* surface.  After completing the square in `s`, introduce `(m,n)` by

```text
r = (m^2-1)/(16*(2*n^2-1)),
s = m*(16*r-1)/(32*r) - r + 5/4 + 1/(32*r),
z = (16*r-1)*n.
```

The inverse on the generic open set is

```text
m = 32*r*(s+r-5/4-1/(32*r))/(16*r-1),
n = z/(16*r-1).
```

This separates the unmarked Humbert coordinates from the orientation much
more cleanly than `(r,s,z)`.  The ramification divisor has the rational
parametrization

```text
r = (1-v^2)/16,
s = (v+3)^3/(16*(v+1)).
```

It is the extra-`I2` branch of `Y_-(8) -> H_8`, **not** the missing
discriminant-237 curve.  Confusing these two curves would manufacture a
spurious short reconstruction.

These identities, the `E7`/`E8` valuations, and the simple branch divisor are
checked by
[`scripts/verify_humbert8_kumar_entrance.sage`](scripts/verify_humbert8_kumar_entrance.sage).
The source is the discriminant-8 ancillary file for
[Elkies--Kumar](https://arxiv.org/abs/1209.3527).

### The skipped `D9+E7` two-neighbor

The ancillary calculation contains a simpler construction *before* the Kumar
equation.  It starts from

```text
Y^2 = X^3 + T*(r+(2*r+1)*T)*X^2
            + 2*r*s*T^4*(T+1)*X + r*s^2*T^7.
```

This has fibers `I5*=D9` at zero and `III*=E7` at infinity, together with
four generic `I1` fibers.  Its root determinant is `4*2=8`, so the
height-four Humbert direction has been absorbed into the reducible fibers.
The exact two-neighbor parameter is

```text
U = (X+s*T^3)/T^4.
```

Putting `X=U*T^4-s*T^3` and `Y=V*T^4` gives a binary quartic in `T`.  It has
the rational point `(T,V)=(0,s)`, and its classical binary-quartic Jacobian is
coefficient-for-coefficient the `E7+E8` Kumar equation displayed above.
Thus this is a characteristic-zero construction of the Humbert-8 entrance,
not merely a matching of invariants.  The full symbolic replay is
[`scripts/verify_humbert8_d9e7_two_neighbor.sage`](scripts/verify_humbert8_d9e7_two_neighbor.sage).
Taking `(0,s)` as the quartic origin, the second rational point `(0,-s)` maps
to the generic polynomial height-four Kumar section.  Its `x` and `y`
coordinates have degrees four and six and are printed implicitly by the same
compact quartic-to-Jacobian formulas; no section search is required.

This corrects another piece of the backtracking order.  The clean source
chain is

```text
D9+E7 Humbert-8 model --explicit 2-neighbor--> E7+E8 Kumar model
                       --level-79 neighbor--> compact MW3 chart.
```

On the determinant-948 Shimura curve, the pre-neighbor model has one generic
MW direction.  Generic torsion is trivial and Shioda's determinant formula
forces its height to be `948/8=237/2`.  Hence the still-missing `H237`
condition can equivalently be phrased as the rank-one Mordell--Weil jump of
this exact two-parameter `D9+E7` family.  That formulation does not by itself
make the high section small, but it removes the previously implicit
Humbert-8 construction and gives a second exact marking for any q60/q80
reconstruction.

This changes the missing-family problem materially.  We no longer seek four
unrelated Kumar coefficients: we seek the single discriminant-237 condition
inside the explicit `(r,s)` plane, together with its lift to `z`.  The expected
map has `(r,s)` rational in Elkies's quotient coordinate `t`; identifying `z`
with a rational multiple of `u` is the natural descent prediction, but that
last identification still requires an exact function-field certificate.

## Exact CM anchors on the correct component

Elkies's rational-point description supplies two CM loci before the non-CM
orbit at `t=+/-14/13`.  In the infinity chart

```text
s=1/t,  v=u/t^3,
v^2 = 16 - 19*s^2 + 88*s^4 - 48*s^6,
```

the two points at infinity are `(s,v)=(0,+/-4)`.  The involution
`beta=w3` acts as `(s,v)->(-s,v)`, so it fixes both.  In contrast, the four
points `(t,u)=(+/-2,+/-32)` form a free orbit under the visible Atkin--Lehner
Klein four-group.  This geometry is now checked in
[`scripts/deconstruct_x0679_quotients.sage`](scripts/deconstruct_x0679_quotients.sage).

The exact Gross-lattice candidates are distinguished by precisely this
stabilizer behavior: the norm-3 vector is `w3`-fixed and the norm-24 vector is
not.  Transport through the inverse-Clifford isometry gives

```text
source coordinate   Gross norm   K3 vector       T_CM
infinity             3           (81,95,-52)     [[2,1],[1,2]]
t=+/-2              24           (70,86,-3)      [[4,0],[0,6]]
CM-43 orbit          43           (169,167,-128)  [[22,1],[1,2]]
```

For the second row the exact checks are

```text
v^2=-158,  div(v)=79,
v^perp = [[4,0],[0,6]],  det(v^perp)=24,
948*158/79^2 = 24.
```

Thus the missing one-variable Kumar map has concrete singular-K3 anchors at
`t=infinity,+2,-2`.  The transport is reproducible with

```text
sage elkies-k3/scripts/transport_cm_delta3_to_k3.sage --target 3
sage elkies-k3/scripts/transport_cm_delta3_to_k3.sage --target 24 \
  --out artifacts/local/elkies-k3/cm-delta24-k3-vector.txt

sage elkies-k3/scripts/transport_cm_delta3_to_k3.sage --target 43 \
  --out artifacts/local/elkies-k3/cm-delta43-k3-vector.txt
```

For the third row, the exact Gross vector is

```text
beta = -11*i-j-5*k,       Nr(beta)=43.
```

It transports to a primitive vector of square `-40764` and divisibility
`948`; the determinant identity is

```text
948*40764/948^2 = 43.
```

This point can now also be located in the explicit Humbert-8 plane.  Two
height-`5/2` section loci meet at

```text
r = -1225/722,
s = -93312/442225.
```

The oriented double-cover coordinate satisfies

```text
z^2 = -43*(11664/6859)^2.
```

In the rational oriented chart this becomes

```text
m = -2468019/407569,
n = -1296/21451 * sqrt(-43)
```

for the displayed choice `z=(11664/6859)*sqrt(-43)` (simultaneously reverse
the signs of `z,n` for the other sheet).  Thus `m` remains rational and the
entire orientation field `Q(sqrt(-43))` is carried by `n`.  This is an exact
normalization constraint for the eventual map from Elkies's `(t,u)` curve:
`m` should descend through `t`, while `n/u` should be rational in `t`.
The two displayed
sections have height Gram

```text
[[5/2,-1/2],[-1/2,5/2]].
```

At fourteen good primes, direct point counts agree with the weight-three
CM-`43` coefficients, including the zeros at inert primes.  This supplies a
strong independent identification fingerprint and makes the discriminant-43
surface a third exact interpolation anchor in the Humbert-8 chart.  The
finite Frobenius list is not, by itself, a characteristic-zero isomorphism
certificate between two Weierstrass presentations; an explicit inverse
neighbor or the missing `H237` normalization is still required for that.
The complete two-second replay is
[`scripts/verify_cm43_humbert8_anchor.sage`](scripts/verify_cm43_humbert8_anchor.sage).

Its primitive closure in the `H2` Kumar frame has index `948`, retains the
`E7+E8` root system, and has MW rank three with reduced height Gram

```text
[[ 5/2,-1/2,-1],
 [-1/2, 5/2, 0],
 [  -1,   0, 4]].
```

The basis can now be marked geometrically.  Let `P1,P2` be the two exact
height-`5/2` sections and let `P3` be the height-four section recovered from
the `D9+E7` quartic.  The two generic `H2` directions specialize as

```text
height-four:  P3,
level-79:     Q79 = 4*P1 - 5*P2 + P3.
```

The Gram matrix gives `P3.Q79=0` and `height(Q79)=237/2`.  Exact group law on
the CM surface gives

```text
x(Q79)=N120/h58^2,    y(Q79)=N180/h58^3,
```

where subscripts denote polynomial degrees.  This recovers `Q79.O=58`
directly and fixes the horizontal marking that any q60/q80 divisor transport
must preserve.

Unlike the two earlier Inose anchors, this point also has a fully explicit
rational semistable q=8 equation and saturated MW basis; see
[`PICARD20_Q8_CHORD_2026-08-21.md`](PICARD20_Q8_CHORD_2026-08-21.md).  The q=8
fiber class itself uses the extra CM section and is not generic.  For reverse
construction, the useful stable presentation is instead the Kumar
`E7+E8/MW3` closure, especially its two height-`5/2` directions.

The corresponding primitive closures of the `H2` Kumar frame determine the
enhanced elliptic frames, not just the transcendental lattices:

```text
Delta   glue index   determinant   enhanced root frame   MW lattice
-3      316          3             E8+E8+A2              0
-24     79           24            E8+E8                 diag(4,6)
-43     948          43            E7+E8                 rank 3, Gram above
```

The finite discriminant-form calculation is checked by
[`scripts/classify_kumar_cm_frame_extensions.sage`](scripts/classify_kumar_cm_frame_extensions.sage).
It identifies the two anchors with standard Inose fibrations.  The Hilbert
class polynomials give the following rational equations:

```text
Delta=-3:   Y^2 = X^3 + T^5*(T-1)^2,
Delta=-24:  Y^2 = X^3 - 51*T^4*X + T^5*(T^2-92*T+1).
```

The first is exactly the `E8+E8+A2` Utsumi No.1 model.  The second is the
`E8+E8` Inose model with MW lattice `diag(4,6)`, up to a rational quadratic
twist.  Their class-polynomial identities, discriminants, and the two `II*`
fibers are checked by
[`scripts/verify_kumar_cm_inose_anchors.sage`](scripts/verify_kumar_cm_inose_anchors.sage).

## First neighbor that actually uses the level-79 direction

The small norm shells of the middle Kumar frame do not leave the
`E7+E8/MW2` presentation in a useful way: through norm 24 they can be formed
inside the roots and the height-4 section direction.  The first basis vector
that is forced to use the height-`237/2` section has integral frame norm 120,
so its quadratic norm is `q=60`.

Fixing that section coordinate to one turns the apparently large norm-120
enumeration into a tiny exact calculation.  It visits 441 recursion nodes and
finds the complete constrained shell of 56 sign-pairs.  The five proper factor
presentations

```text
(a,b) = (2,30), (3,20), (4,15), (5,12), (6,10)
```

all expose a rank-14 root frame with component invariants

```text
(rank, roots, determinant) = (8,240,1) + (6,72,3),
```

hence `E8+E6` and Mordell--Weil rank three.  The balanced `(5,12)` witness is

```text
v=(0,0,-1,-1,-1,-1,-1,0,0,0,0,0,0,0,0,0,1),
```

and its pinned child frame is
[`data/fibrations/kumar_q60_e8_e6_mw3_frame.txt`](data/fibrations/kumar_q60_e8_e6_mw3_frame.txt).
The saturated Mordell--Weil height lattice has determinant 316 and a reduced
Gram

```text
[[4,0,0],[0,20/3,1],[0,1,12]].
```

The denominator identifies the middle generator as the only nonzero `E6`
component class; its local correction is `4/3`.  Shioda's height formula then
gives zero-section intersections `0,2,4` for the three displayed generators.
Their pairwise section intersections are respectively `4,6,7`.  This replaces
the pole-58 Kumar generator by a three-section system whose largest pole is
only four.

Putting `IV*` at `T=0` and `II*` at infinity forces the short Weierstrass
ambient form

```text
Y^2 = X^3 + T^3*(a0+a1*T)*X
            + T^4*(b0+b1*T+b2*T^2+b3*T^3).
```

Thus there are only six ambient coefficients before the base and Weierstrass
scalings.  The exact one-dimensional target is cut out by sections with pole
orders `0,2,4`, rather than by the original pole-58 section.

This pole profile cannot be improved by another integral MW basis.  In the
reduced height lattice, every section with `P.O<=2` lies in the rank-two span
of the first two generators; hence every basis contains a section with
`P.O>=4`, and the displayed basis attains the sharp profile `0,2,4`.  The
finite exact enumeration is
[`scripts/verify_q60_pole_profile_optimal.sage`](scripts/verify_q60_pole_profile_optimal.sage).

The height-4 section nevertheless removes half of the ambient unknowns
without elimination.  Away from the two endpoint-degeneracy divisors, use the
base and Weierstrass scalings to write

```text
x(P1)=1+c1*T+c2*T^2+c3*T^3+T^4.
```

The first three coefficients of `y(P1)` are forced by the equation at zero,
and its last four are forced at infinity.  For a relative leading sign
`epsilon=+/-1`, the residual `y(P1)^2-x(P1)^3` is then supported only in
degrees `3,...,8`.  Its degree-3 and degree-8 coefficients give `a0,a1`, and
the four intervening coefficients give `b0,...,b3` linearly.  Thus the entire
`E8+E6` surface together with `P1` has an explicit three-parameter normal
form; only the `P.O=2,4` sections remain to be imposed.  The derivation and
the full coefficient formulas are checked and printed by
[`scripts/verify_q60_height4_normal_form.sage`](scripts/verify_q60_height4_normal_form.sage).

In this normal form let `ell=c3-epsilon*c1`.  Exact factorization gives

```text
ell | a0,       ell^2 | b0,       ell^3 | b3.
```

Consequently the discriminant-3 enhancement is already a visible boundary
of the normalized section chart.  The other factor of `b0` supplies the
generic discriminant-24 branch where `b0=0` but `a0` remains nonzero.  This is
the correct small starting system for a finite-field seed and local lift; a
classical equation for the full Humbert-237 surface is no longer required.

That boundary calculation can now be completed exactly.  Retain the constant
normalization of `x(P1)`, move the extra `A2` fiber to `T=1`, and write the
leading coefficient of `x(P1)` as `d^2`.  On the nondegenerate factor of
`b0=0`, requiring a triple residual-discriminant root at `T=1` gives a
zero-dimensional three-equation system.  Its nonsingular solution
`(c1,c3,d)=(16,12,8) mod 31` Hensel-lifts and rationally reconstructs to

```text
c1=20/9,  c2=334/729,  c3=68/729,  d=-1/27.
```

The resulting equation and height-4 section are

```text
Y^2 = X^3 - (4096/19683)*T^3*(7*T+2)*X
            - (262144/14348907)*T^5*(T^2+34*T+19),

x(P1) = (T^4+68*T^3+334*T^2+1620*T+729)/729,
y(P1) = (-T^6-102*T^5-2235*T^4+188*T^3
         +49977*T^2+65610*T+19683)/19683.
```

Its discriminant is

```text
-(2^40/3^27)*T^9*(T-1)^3*(T^2+71*T+32),
```

so the fibers are exactly `III*+I3+II*+2 I1`.  The root determinant
`det(E7+A2)=6` times the section height four is 24.  Moreover, the resulting
frame is exactly isometric to the discriminant-24 frame transported through
the `q=60` neighbor, not merely a lattice with the same determinant.  The
equation, section identity, valuations, and residual discriminant are checked
by
[`scripts/verify_q60_delta24_anchor.sage`](scripts/verify_q60_delta24_anchor.sage).
This supplies the first explicit equation-level point on a correct compact
neighbor.  The broader CM-stability audit below subsequently found a better
deformation frame, so the `q=60` point is now a cross-check rather than the
preferred endpoint chart.

Transporting the two primitive CM closures through this same `q=60` neighbor
gives sharper boundary data:

```text
Delta   q=60 enhanced frame   MW data
-3      E8+E8+A2              rank 0
-24     E8+E7+A2              rank 1, height 4
```

In the six-coefficient equation these force explicit valuation divisors.  The
generic `E6` fiber at `T=0` has orders `(3,4,8)`.  At the discriminant-24
anchor it enhances to `E7`, so `b0=0` and the orders become `(3,5,9)`.  At the
discriminant-3 anchor it enhances to `E8`, so `a0=b0=0` and the orders become
`(4,5,10)`.  In particular the standard discriminant-3 equation is the point

```text
(a0,a1,b0,b1,b2,b3)=(0,0,0,1,-2,1),
```

for which `a6=T^5*(T-1)^2`; the additional `A2` is the `IV` fiber at `T=1`.
The coefficient and valuation assertions are checked by
[`scripts/verify_e8e6_cm_boundaries.sage`](scripts/verify_e8e6_cm_boundaries.sage),
while the enhanced root frames are certified independently by
[`scripts/classify_kumar_cm_frame_extensions.sage`](scripts/classify_kumar_cm_frame_extensions.sage).
Thus the missing `H237` curve must pass through a codimension-two point of the
`a0=b0=0` boundary and through the `b0=0` divisor at the other CM orbit.

The exact discovery and verification commands are

```text
sage elkies-k3/scripts/search_alternate_fibrations.sage \
  --frame elkies-k3/data/fibrations/kumar_e7e8_mw2_frame_2.txt \
  --min-qnorm 60 --max-qnorm 60 --fixed-coordinate 16:1 \
  --enum-baseline-cap 1000 --enum-restarts 1 --enum-cap 1 \
  --proper-factors-only --max-candidates 1000 --report 5 \
  --quiet-candidates \
  --frames-dir artifacts/local/elkies-k3/kumar-h2-q60-frames \
  --out artifacts/local/elkies-k3/kumar-h2-q60-level79.txt

sage elkies-k3/scripts/verify_fibration_neighbor.sage \
  --parent elkies-k3/data/fibrations/kumar_e7e8_mw2_frame_2.txt \
  --child elkies-k3/data/fibrations/kumar_q60_e8_e6_mw3_frame.txt \
  --q 60 --a 5 --b 12 \
  --v 0,0,-1,-1,-1,-1,-1,0,0,0,0,0,0,0,0,0,1

sage elkies-k3/scripts/analyze_mw3_branch.sage \
  --frame elkies-k3/data/fibrations/kumar_q60_e8_e6_mw3_frame.txt \
  --name kumar-q60-e8-e6-mw3

sage elkies-k3/scripts/verify_e8e6_cm_boundaries.sage

sage elkies-k3/scripts/verify_q60_pole_profile_optimal.sage

sage elkies-k3/scripts/verify_q60_height4_normal_form.sage

sage elkies-k3/scripts/verify_q60_delta24_anchor.sage

sage elkies-k3/scripts/verify_humbert8_d9e7_two_neighbor.sage

sage elkies-k3/scripts/verify_kumar_cm43_q8_q9_factor.sage

sage elkies-k3/scripts/verify_cm43_q8_short_section.sage
```

## Backtracked CM-stable neighbor search

The choice `(a,b)=(5,12)` is not intrinsic to the norm-120 class.  The five
proper presentations of `60` have exact generic pole profiles and
discriminant-24 behavior

```text
(a,b)   optimal generic P.O   frame at Delta=-24       surviving MW rank
(2,30)  0,0,13                E8+E7+A2                  1
(3,20)  0,1,6                 E8+E7+A2                  1
(4,15)  0,3,4                 E8+E8                     2
(5,12)  0,2,4                 E8+E7+A2                  1
(6,10)  0,2,5                 E8+E8                     2
```

Thus `(4,15)` is the best `q=60` presentation for deforming from the
discriminant-24 anchor, but none is root-stable.  This prompted an exact
search beyond `q=60`.

The search does not enumerate the enormous raw root shells.  In the Kumar
frame a vector involving the level-79 section has orthogonal contribution
`237/2`; for `60<=q<237` its coefficient is necessarily `+/-1`.  Weyl
reduction in `E7+E8`, together with the independent sign of the height-four
section, gives a complete finite list of orbit representatives.  Through
`q=80`, all 2,869 proper factor presentations were transported to the
discriminant-24 closure.  No root-stable presentation occurred, but 313 had
root jump one.

The best compact additive chart in this bounded search is

```text
q=80,  (a,b)=(4,20),
witness=(2,3,4,6,5,4,3,6,9,12,18,15,12,8,4,2,1).

generic:    E6+D5+A3,       MW rank 3, optimal P.O = 0,0,3
Delta=-24: E6+D5+A3+A1,    MW rank 3, optimal P.O = 0,0,0.
```

Only one new root appears at the CM point and all three generic MW directions
survive.  The generic frame is pinned in
[`data/fibrations/kumar_q80_e6_d5_a3_mw3_frame.txt`](data/fibrations/kumar_q80_e6_d5_a3_mw3_frame.txt).
The orbit enumeration, CM transport, saturated height lattices, component
corrections, and optimal-pole claims are all checked by
[`scripts/classify_kumar_cm_frame_extensions.sage`](scripts/classify_kumar_cm_frame_extensions.sage).
The bounded audit is reproduced by

```text
sage elkies-k3/scripts/classify_kumar_cm_frame_extensions.sage \
  --search-max-q 80
```

This frame also has a four-parameter equation chart.  Put `I1*`, `I4`, and
`IV*` at `T=0,1,infinity` and normalize the leading nodal cubic at zero:

```text
A = T^2*(-3 + p*T + q*T^2 + r*T^3),
B = T^3*( 2 + b1*T + b2*T^2 + b3*T^3 + b4*T^4 + e*T^5),
r = -3*d^2 + 3 - p - q.
```

The four `I4` discriminant jets at `T=1` determine `b1,...,b4` linearly,
leaving exactly `(d,p,q,e)`.  The discriminant has the certified form

```text
Delta = constant * T^7*(T-1)^4*R5(T).
```

At the discriminant-24 anchor the added `A1` is simply a double root of the
residual quintic, so its boundary is `Res_T(R5,R5')=0`.  The full coefficient
formulas and an exact nonvanishing specialization are checked in
[`scripts/verify_q80_ambient_normal_form.sage`](scripts/verify_q80_ambient_normal_form.sage).
At the discriminant-24 anchor this remains the least-degenerate local chart.
The complete discriminant-43 chamber calculation below shows that the
apparent q=60 shortcut there is entirely fixed.  Thus the q=80/CM-24 chart is
the leading noncollapsed deformation chart; q=60 at CM-24 remains an exact
comparison route.

There is now an exact characteristic-zero boundary model in this chart:

```text
A = T^2*(-3 + 9/4*T - 9/4*T^2 + 9/4*T^3),
B = T^3*(2 - 315/32*T + 9*T^2 - 9/16*T^3 - 27/32*T^5).
```

It is defined over `QQ`; the three selected marked sections are defined over
the CM compositum `QQ(sqrt(-3),sqrt(-6))`.  Its resolved characteristic-zero
tangent cone is

```text
(tau0-8/87*tau1)*(tau0-1/12*tau1).
```

Modulo seven these are exactly the two slopes `5` and `3` below.  This
certifies the CM boundary point, its marking, and two Q-rational tangent
directions; it does **not** algebraize either one-parameter branch.  The exact
certificate is
[`scripts/verify_q80_cm24_rational_model.sage`](scripts/verify_q80_cm24_rational_model.sage).

With the local parameter normalized by `p=9/4+h`, the two surface tangents in
`(d,p,q,e)` are exactly

```text
(8/87, 1, -24/29, -45/116),
(1/12, 1, -45/52, -261/832).
```

Both branches solve the full marked section system over
`QQ(sqrt(-3),sqrt(-6))[[h]]` through order 18, with every residual coefficient
through that order checked exactly.  This is a formal certificate, not an
algebraization.  At degree at most seven only the normalized coordinate
`p=9/4+h` has a unique Pade recognition.  Replay it with
[`scripts/extend_q80_cm24_branches_qq.sage`](scripts/extend_q80_cm24_branches_qq.sage).

The corrected local calculation at the fully marked GF(7) CM-24 seed has
Jacobian rank 37 in 39 resolved variables.  Its quadratic tangent cone splits
as

```text
tau0^2-tau0*tau1+tau1^2=(tau0+2*tau1)*(tau0-3*tau1).
```

Both reduced lines lift uniquely through order `h^20` over GF(7).  A bounded
Pade test with numerator and denominator degrees at most ten recovers only
24/39 and 22/39 active coordinates, respectively, so neither finite-field jet
has been algebraized.  This is only a bounded negative recognition result.
Replay the cone and the bounded recognition test with
[`scripts/verify_q80_rank19_deformation_gf7.sage`](scripts/verify_q80_rank19_deformation_gf7.sage)
and
[`scripts/extend_q80_rank19_branches_gf7.sage`](scripts/extend_q80_rank19_branches_gf7.sage).

The two branches are no longer equally opaque.  On the slope-`5` reduction of
the characteristic-zero `8/87` branch, the filtered quotient dimensions in
degrees zero through five are `1,5,15,33,48,63`, agreeing from degree three
onward with `15*n-12`.  Its order-85 ideal has two cubic generators and twelve
new quartic generators modulo cubic multiples.  A deterministic centered and
homogeneous export is
[`../artifacts/generated-results/q80-cm24-slope-8-87-gf7-ideal.json`](../artifacts/generated-results/q80-cm24-slope-8-87-gf7-ideal.json).
The apparent degree 15 and arithmetic genus 13 describe the singular
surface-coefficient image. A later modular `(D,Q)` plane projection has
generic degree five over `D` and normalization genus zero, so the known
genus-two Shimura curve must instead occur on the marked-section cover if
this branch is the target. The plane relation and marked cover still require
characteristic-zero substitution certificates. The homogeneous normalization
is now explicitly a rational normal degree-eight curve, and its CM osculating
flag yields a modular parameter with surface-function degrees
`D,P,Q,E = 5/4,10/8,8/6,15/12`; all fourteen ideal generators vanish after
exact substitution. This gives a finite ansatz for the characteristic-zero
lift rather than a reason to extend unstructured jets.

The first marked section lies on a genus-two cover of the rational coefficient
line whose absolute Igusa invariants `(4,2,5)` do not match `(4,4,3)` for the
known source model modulo seven. The second marked section lies on a distinct
genus-three cover; together the markings give a genus-six biquadratic curve,
with genus-one third quotient. Thus descent through these two short sections
is a high-genus detour, not the desired source marking. This does not reject
the unmarked slope-5 surface branch by itself, because `Q79` involves the
third pole section. Exact lattice transport now makes that statement precise
in the generic q80 optimal basis:

```text
height-4 = -G2,
Q79      = -3*G1 - 2*G2 + 4*G3.
```

Their q80 heights are `4,120`; the latter has `P.O=59`. Since the height-four
class must descend to the rational Shimura quotient, the slope-5 `G2` cover
rules out the two-short-section descent but not the unmarked surface: `Q79`
uses `G3`. The later split-prime trace comparison conditionally rejects the
other exact tangent `1/12`, so the live unmarked candidate returns to `8/87`
with a different `P3/Q79` marking or orientation. This remains conditional on
global algebraization of the finite modular ideals, rather than a
characteristic-zero exclusion theorem.
The horizontal map is certified by
[`scripts/classify_kumar_cm_frame_extensions.sage`](scripts/classify_kumar_cm_frame_extensions.sage).
The slope-5 ideal's fourteen generators
nevertheless pass five withheld jet orders through order 90. This exact audit is
certified by
[`scripts/analyze_q80_rank19_marked_cover.sage`](scripts/analyze_q80_rank19_marked_cover.sage)
and interpreted in
[`Q80_TO_ROOTLESS_PATH_2026-08-21.md`](Q80_TO_ROOTLESS_PATH_2026-08-21.md).
It rules out further all-coordinate Padé extension of the chosen P1/P2
marking coordinates. The
slope-5 P1 branch set has trivial projective stabilizer over `GF(49)`, unlike
the known bielliptic source, and its genus-one third quotient has `j=0` rather
than either known quotient value `1,3`.

The slope-3 branch has now been continued uniquely through order 230. Its
first centered relations occur in degree six. The initial `GF(7)` kernel of
dimension 17 passes sixteen withheld coefficients, but split-prime replay
shows it is exceptional: ordinary split primes give rank 195/kernel 15,
`p=97,103` give kernel 16, and only `p=7` gives 17. The generic fifteen-sextic
bounded ideal is a curve of projective degree 48 and arithmetic genus 94; its
irreducible `(P,D)` plane projection still has total degree 32 and bidegree
`(15,32)`. Thus the two extra `GF(7)` sextics must not be lifted. Exact
local-delta calculations at both `p=73` and `p=79` give total delta 464 on
the irreducible degree-32 plane curve, hence normalization genus one. This is
the expected elliptic quotient shape, but exact point counts distinguish it
from both known source factors: the candidate traces are `15,22` at
`p=73,127`, versus `(-9,7)` and `(8,0)` for the two quotients (and neither
twist matches). Conditional on the fifteen-sextic ideal being the true global
branch, `1/12` is rejected as the source component. The live branch choice
returns to `8/87`, where the rational unmarked surface may still acquire the
correct genus-two cover through the third pole section `P3`/`Q79` rather than
through the failed P1/P2 marking.

The first bounded order-230 `P3` recovery adds a marking-level obstruction:
its pole and an independent numerator coordinate lie on the same genus-one
`j=0` field as the third P1/P2 quotient, with squarefree model
`w^2=t^4+3t^3+3t^2+1`. Thus this selected three-section orientation is still
not the known genus-two source. This does not reject the unmarked rational
surface branch; it calls for a different CM orientation or direct descent of
the combination `Q79`.

That direct combination has now been tested on the same formal branch.  At a
good elliptic-base value, exact group law forms
`Q79=-3*G1-2*G2+4*G3`; all four sign classes of the selected triple (up to
simultaneous negation) fail both rational Padé bounds `100/100` and quadratic
relations of parameter degree at most `70`, using an order-230 jet with sixteen
withheld coefficients.  Hence a sign-only repair of this selected marking is
not the missing source cover.  This remains bounded modular evidence: it
narrows the live option to a genuinely different CM basis/marking or a direct
construction in the rational Humbert-8 `(m,n)` chart; it does not reject the
unmarked slope-`8/87` coefficient curve in characteristic zero.  The exact
commands and four artifacts are recorded in
[`Q80_TO_ROOTLESS_PATH_2026-08-21.md`](Q80_TO_ROOTLESS_PATH_2026-08-21.md).

The other CM24 third-section hit is not such a different marking.  Exact group
law gives `G2'=-G2, G3'=G3-G2`, under which
`-3*G1-2*G2'+4*G3'` is the same `Q79` point.  The finite-pole and
infinity-pole seed representatives therefore lead to the same failed direct
cover test.

The q80-to-rootless lattice continuation is now explicit. The exact path uses
neighbor norms `4,4,12,12,4,6` and passes through MW ranks
`3,4,5,6,13,16,17`, ending in a rootless frame integrally isometric to the
pinned determinant-948 rank17 lattice. See
[`Q80_TO_ROOTLESS_PATH_2026-08-21.md`](Q80_TO_ROOTLESS_PATH_2026-08-21.md)
and
[`scripts/verify_q80_to_rootless_path.sage`](scripts/verify_q80_to_rootless_path.sage).
There is now a second exact lattice suffix after the first four steps. A
bounded q2--q6 continuation from a selected alternate fifth q4 class finds a
q6 rootless child. The latter frame is not
integrally isometric to the previously pinned rootless frame, but the full
q80-to-alternate-q4-to-q6 NS composite is unimodular and it is a genuine
MW17 fibration on the same K3. See
[`scripts/verify_q80_alternate_fifth_q6_rootless.sage`](scripts/verify_q80_alternate_fifth_q6_rootless.sage).

This is an exact geometric-lattice route. The first two effective pencils
have been constructed on the q80 equation, and the two q12 steps have exact
CM24/mod-73 equation gates. A separate marked pair-14 fifth equation has
fibers `I6+I6+I5+I2+5I1`, but its CM root data `(15,82,360)` do not match the
selected q4 lattice witness `(16,66,2048)`. The productive pair-23 fifth
equation does match that witness.  The final degree-two q6 pencil is now also
exact over `GF(73)`: Smith saturation gives
`q_sat=(q0+63)/(s-27)`, the kernel rows are
`(1,0,0,41,48)` and `(0,1,0,6,72)`, and the unit-preserving Jacobian has
fibers `2I5+2I4+I2+4I1`, root data `(15,66,800)`.  The generated certificate
[`q80-final-q6-saturated-module-gf73.json`](../artifacts/generated-results/q80-final-q6-saturated-module-gf73.json)
has SHA256
`7d3866855b9995b733193de9c5d5e3ba1cea6aa1292141e06d0d5011f28975e3`;
it includes the explicit `q0(s,R)` and globally minimal polynomial
Weierstrass coefficients.
Characteristic-zero lifting and field-of-definition tracking remain open.

This q80 chain uses the `H2=diag(4,237/2)` polarization and is downstream
from the source problem.  It is not the recovered `H3=[[21/2,3],[3,46]]`
polarization: after doubling, the forms have the same determinant and Smith
invariants but minima eight and 21 respectively, and are not integrally
isometric.  Likewise, the final-q6 Smith factor `s-27` is a local chord-module
saturation, not an explanation of the recurrent norm-eight neighbor.  The
latter remains the Humbert-8/height-four `H2` entrance or a
Noether--Lefschetz return/collapse at special fibers.

This downstream chain does not simplify the missing source marking.  Exact
transport of the level-79 frame direction makes its MW quotient height jump
from `120` in q80 to `148039` after the first q4 step and then grow at every
later step.  A direct q80-to-q60 conversion is also large: among all five q60
presentations its minimum neighbor norm is `q=400`.  Consequently neither
route is an efficient substitute for the marked H3 source descent.  The
source curve itself is now normalized exactly as the `H21 cap H92` level-474
component; the remaining source gate is to identify the rational H21 entrance
data with the signed H3 q6 divisor and track its field of definition.  The
rational `(m,n)` recovery of `H8 cap H237` remains an H2 comparison problem,
not a prerequisite for the H3 construction.

At CM-43, the q=80 frame has `D4+D5+E6` roots and MW rank three.  A CM-only
basis is polynomial, but the three generic MW directions specialize with
pole profile `(0,0,3)`, exactly the generic bound.  In the q=60 `(5,12)`
presentation the generic directions instead specialize with profile
`(0,1,1)`, and the child returns to an `E7+E8/MW3` frame integrally isometric
to the original Kumar closure.  Its fiber in the Kumar CM Neron--Severi basis
is

```text
(5,12,0,0,-1,-1,-1,-1,-1,0,0,0,0,0,0,0,0,0,1,0).
```

This root-frame return is not a movable q=60 pencil at CM-43.  In the explicit
divisor basis `[F,O,E7(7),E8(8),P1,P2,P3]`, fixed-component reduction first
gives

```text
D60_raw  ->  D1 = Q79 + 4*O - 43*F.
```

The CM-only section `S=P1-P2` then has `D1.S=-4`.  Subtracting `4*S` gives
`F+T`, where `T=P3-P2` is again fixed, and subtracting `T` leaves exactly the
old fiber `F`.  Hence

```text
D60_raw = F + (initial O/E7 fixed part) + 4*(P1-P2) + (P3-P2).
```

The exact certificate is
[`scripts/verify_cm43_marked_divisor_transport.sage`](scripts/verify_cm43_marked_divisor_transport.sage).
CM-43 therefore cannot seed the q=60 equation; it only certifies the marking
and the boundary collapse.

The semistable Picard-20 q=8 fibration is not directly adjacent to the Kumar
closure.  The exact norm-16 shell has 1,421,331,656 signed vectors but only
303 `W(E7)xW(E8)` dominant orbits; all have been classified, and none of the
292 primitive q=8 neighbors is the `A7+A4+A3+A2` frame.  Replay this complete
orbit calculation with

```text
sage elkies-k3/scripts/classify_kumar_cm43_q8_orbits.sage
```

The same complete classification does find a shorter decomposition of the
q=60 CM automorphism:

```text
E7+E8/MW3 --q=8--> E7+E8/MW3 --q=9--> E7+E8/MW3.
```

For the first step the exact witness is

```text
(156,-78,0,0,-78,0,-78,0,0,0,0,0,0,0,0,-1,-155,-32).
```

Thus in the pinned marked NS coordinates `[F,O,frame18]` its full isotropic
fiber is `(2,4,witness)`.  The marked closure Gram is
[`data/fibrations/kumar_cm43_marked_e7e8_mw3_frame.txt`](data/fibrations/kumar_cm43_marked_e7e8_mw3_frame.txt).
Projection to the explicit section basis `(P1,P2,P3)` is exactly `(1,-2,0)`.
Exact function-field group law gives

```text
R=P1-2*P2,  height(R)=29/2,  R.O=6,
x(R)=N16/h6^2,  y(R)=N24/h6^3.
```

This pole-six section gives a compact formula for the q=8 fixed-component
factorization at the CM-43 boundary; it is dramatically smaller than the
pole-58 `Q79` section.
The exact certificate is
[`scripts/verify_cm43_q8_short_section.sage`](scripts/verify_cm43_q8_short_section.sage).

It does **not** define a fibration class on the generic rank-19 family.  In
the marked MW quotient, with `P3=(0,0,1)` and `Q79=(4,-5,1)`, the primitive
height-`40764` vector orthogonal to the generic plane is
`W=(116,92,29)`, and

```text
P1-2*P2 = -1/4*P3 + 27/79*Q79 - 1/316*W.
```

The residual has height `129/316 = 9*(43/948)`.  Equivalently, in the full
glue-211 primitive closure the q=8 frame witness has nonzero coordinate
`1/316` along the added CM rank-one summand, while the q=60 witness has
coordinate zero.  Thus q=8 is only a low-pole CM boundary factorization;
q=60, whose horizontal part is `Q79`, remains an integral class of the generic
rank-19 lattice, but its CM-43 specialization is fixed as described above.
This full membership check is certified by
[`scripts/verify_cm43_q8_generic_membership.sage`](scripts/verify_cm43_q8_generic_membership.sage).

The q=8 class itself also collapses after completing the chamber reduction.
The first O/E7 cascade gives `D8'=O+(P1-2*P2)-2F`; the sections `-P2` and
`P1-P2` are then fixed successively, leaving `F`.  Equivalently,

```text
D8_raw = F + (initial O/E7 fixed part) + (-P2) + (P1-P2).
```

Thus the q=8 then q=9 return records a boundary fixed-part decomposition, not
a second elliptic pencil on the CM-43 equation.

Its marked horizontal projection is `P1-2*P2`, of height `29/2`.  In the q=8
child the q=60 fiber begins `(3,3,...)`, proving that the second move has q=9.
The characteristic-zero lattice certificate, including the full second
witness, ADE components, and return isometry, is
[`scripts/verify_kumar_cm43_q8_q9_factor.sage`](scripts/verify_kumar_cm43_q8_q9_factor.sage).

This arithmetic-optimal intermediate is not the H2 Humbert-8 `D9+E7` chart.
Exactly ten q=8 orbits have roots `D9+E7`; the cheapest root-equivalent one
expresses the q=60 fiber with `(a,b)=(17,17)`, hence second q=289.  The actual
H2 two-neighbor is selected by horizontal projection `+P3`, not by
root invariants alone.  Its q=60 coordinates begin `(56056,44,...)`, hence
second q=2,466,464.  The explicit `D9+E7` equation can still make its rational
functions shorter, but its marked q=8 class is not the optimal q=8 step.
One exact source-marked witness is

```text
(0,0,0,0,0,0,0,1,1,1,1,1,1,2,3,4,0,-1),
```

and subtracting the marked height-four frame vector from it lies in the
`E7+E8` root span.  This certifies that it is the `P3` class rather than merely
a determinant-equivalent `D9+E7` frame.

Finally, the explicit CM equation fixes the level marking.  With the two
height-`5/2` sections and the height-four section ordered as `(P1,P2,P3)`,

```text
Q79 = 4*P1 - 5*P2 + P3.
```

The eight normalized closure signs have now been audited.  Glues `211` and
`737` carry the ordered pair `(P3,Q79)` to the explicit equation; glue `53`
does not.  Its old horizontal coordinates `(5,-4,2)` were therefore an
unmarked CM isometry.  The equation-level q=60 divisor has marked horizontal
coordinates `(4,-5,1)`.

The generic `q=60` construction remains the first direct lattice bridge out of
the correct Kumar anchor that cannot ignore the level-79 generator.
Geometrically it predicts the compact generic fiber configuration
`II*+IV*+6 I1`.  The CM-43 equation cannot be used to execute that pencil,
because the added CM sections make the specialization fixed.  At CM-24, the
known `(5,12)` equation retains only the height-four MW direction, while the
`(4,15)` presentation retains two of the three generic directions.  The
`q=80,(4,20)` CM-24 chart retains all three and already has two certified
formal branches, so it is the current leading reconstruction route.

## What this proves and what it does not

This proves that there are three exact Kumar-type lattice entrances compatible
with the recovered determinant-948 K3 and identifies `H2` as the entrance for
Elkies's quotient construction.  It also explains why a direct search inside
an arbitrarily normalized `E8+A2^3` Kodaira chart was poorly anchored.

It does **not** yet give the rational functions from the `t`-line to the four
Kumar/Clebsch--Igusa coefficients, nor the birational maps realizing the now
certified lattice path on the Kumar equation. The CM boundary equations and
all six lattice neighbors to the rootless frame are known, but the local
branch of the `H2` family and its rational twist/marking have not yet been
recovered. The next exact tasks are therefore:

1. use the three exact CM anchors to recover the `H237` normalization inside
   the explicit Humbert-8 plane and determine how `u` selects the quadratic
   twist and level-79 section;
2. distinguish and algebraize the two exact characteristic-zero formal
   branches with tangents `tau0/tau1=8/87,1/12` in the `q=80,(4,20)` CM-24
   chart, which retains all three generic MW directions.  They are certified
   over the CM compositum through order 18; their reductions lift through
   `h^20`, but neither has complete bounded coordinatewise Pade recognition.
   The `8/87` line has an exact modular rational surface image, but its P1 and
   P2 marking fields combine to genus six and its P1 genus-two quotient is not
   the known source modulo seven. Treat this **two-short-section descent** as a
   high-genus detour, while retaining slope 5 as an unmarked candidate whose
   true level-79 marking may involve `P3`. The selected P3 coordinates still
   lie on the same genus-one `j=0` field, so change the CM orientation or
   descend `Q79` directly. The generic fifteen-sextic model of the `1/12`
   line has the wrong elliptic traces at two good primes and is conditionally
   rejected. Use the global H237 marking or direct component equations to
   decide the remaining orientation; do not merely increase every
   coordinatewise Pade bound;
3. in parallel derive the compact `q=60,(4,15)` chart, which retains two MW
   directions, and compare its equation complexity against q=80.  Do not use
   CM-43 as a q=60 equation seed: both its marked q=60 and q=8 classes reduce
   to the old fiber after all CM-only fixed sections are removed;
4. execute the now-certified `q=4,4,12,12,4,6` neighbor chain on the generic
   q80 equation, reducing each fiber to its nef chamber and tracking sections,
   components, and fields of definition. The first q4 class is already an
   exact degree-two, zero-MW-projection class in `L(2O+4F)` and is
   fully nef: the height shell handles all sections and root primitivity
   excludes a negative bisection. Its `D5+E6` component resolution gives
   `L(D)=<T^2,x-T>` and the exact first coordinate `U=(x-T)/T^2`. The ambient
   child has `D9+A3`, the marked rank-19 collision upgrades it to `D9+A4`,
   and CM24 has `D9+A5+2A1`. The second q4 class is likewise fully nef: its
   saturated rank-four MW height shell handles every section and root
   primitivity excludes a negative bisection. Its exact first-child coordinate
   is `W=(X-3v^3-x1v-x0)/v^2`, with `v=U-d+1`, the nodal center `x0`, and its
   first derivative `x1`. The identity clears exactly to a cubic. At CM24 it
   has `I3*+IV*+3I2+I1`, hence `D7+E6+3A1/MW2`; the finite valuations
   `(3,4,8)` are `IV*`, not `I2*`, and agree with the transported pinned
   frame. The third q12 divisor reduces to old degree three with an integral
   MW projection of norm eight. It is exactly
   `S+2O+2F+root_correction`, where the height-eight section has `S.O=2`, so
   the next pencil is a marked trisection in `L(S+2O+2F)`. At CM24 the
   saturated-basis profile gives `Q_CM=P1+3P2`, a polynomial height-three
   section over `QQ(sqrt(-6))`; hence the generic-fiber space is already
   `<1,X,(Y+y(Q_CM))/(X-x(Q_CM))>`. The unique effective integral lift has
   `D=Q_CM+2O+4F+R` and old-component coefficients
   `A1:(0),(1),(0)`, `E6:(1,2,3,2,1,2)`,
   `D7:(2,4,3,3,5,6,3)`. The transported CM target is
   `2A6+3A1/MW3`. Pure translations of the chord have only two genus-one
   completions, with `D8+E7+A1` and `D5+E6+A3+A1`; both are wrong, so the
   next vertical ansatz must use a nonzero `X` coefficient. All 2,401 naive
   linear-X/quadratic-translation tuples over `GF(7)` still have branch
   degree 17 or 23. The correct nine-dimensional ambient is
   `a(W)+bX+c(W)z_Q`. Two selected-I2 rows and one E6 row combine with the
   four resolved D7 rows from the complete ideal `(Y,U^2,ZU,Z^3)`. The exact
   matrix has rank seven, its kernel gives `Vnew=N1/N0`, and clearing the
   chord leaves `X-Qx` times a cubic of W-degree nine. At `p=73`, both
   `Vnew=1` and `Vnew=7` normalize to genus one. The earlier unresolved-cusp
   rows instead gave genera four and three and remain rejected. Do not impose
   both quadratic roots, which corresponds to an ineffective lift.
   Canonical-origin interpolation at `p=73` now recovers the full short
   Jacobian with `deg(A),deg(B)=(8,12)` and discriminant
   `2I7+3I2+4I1`, exactly the transported `2A6+3A1/MW3` CM24 frame. The next
   gate is its characteristic-zero lift and fourth-neighbor transport.
   The genuine
   alternate q8 move gives `E6+A7/MW4` (and `E7+A7+2A1/MW2` at CM24), but
   loses one generic MW direction and therefore does not replace the q12
   rank-growing step. The
   lattice endpoint is already the recovered rootless frame; after lifting
   this third step, three downstream geometric pencils remain.

This order also matches Elkies's published
[K3 reconstruction method](https://arxiv.org/abs/0802.1301): choose a simpler
elliptic fibration containing a primitive sublattice, find one exact point
modulo a suitable prime, lift the one-dimensional branch, and recognize two
low-degree modular functions.  The source explicitly says that the
`X_0^6(79)` calculation was completed but its full data were omitted from the
paper.  The exact CM-24 q=60 anchor, the q=80 two-branch formal certificate,
and the CM-43 fixed-part audit now provide the boundary data needed to replay
that method without mistaking a special-fiber collapse for a generic pencil.

The `E8+A2^3` and E6 frames remain useful cross-checks after this upstream
model is recovered; they are no longer the starting assumptions.
