# Affine sections on the six-root Mestre moduli space

Status: exact symbolic computation and explicit identity families.  The
universal elimination alone does not prove that an extra section is
independent of the displayed Mestre points; for the centers
`(0,25,95,143,168,205)`, the later exact specialization certificate proves
one companion supplies a thirteenth generic direction after the
split-infinity base change.  There is no generic-rank equality or
rank-21/rank-30 claim here.

## Affine normalization

Fix two labelled roots at `0` and `1` and write

\[
 q(X)=X(X-1)h(X),\qquad
 h(X)=X^4+c_1X^3+c_2X^2+c_3X+c_4.
\]

The rational six-root locus is obtained by substituting

\[
 h(X)=\prod_{i=2}^5(X-r_i).
\]

The open distinct-root locus additionally requires
\(c_4h(1)\operatorname{disc}(h)\ne0\).  For

\[
 q(X-T)q(X+T)=g(X,T)^2-T^2R(X,T),
\]

where the monic degree-six `g` matches degrees 12 through 6, the condition
\(\deg_XR\le4\) is

\[
\begin{aligned}
M={}&c_1^5+c_1^4-6c_1^3c_2-5c_1^2c_2+8c_1c_2^2+7c_1^2c_3
 +4c_2^2+6c_1c_3\\
 &-12c_2c_3-8c_1c_4-c_1-c_2-c_3-16c_4-1=0.
\end{aligned}
\]

This is the three-dimensional affine-normalized Mestre moduli hypersurface.

## Smallest affine extra-section ansatz

Use

\[
x(T)=x_0+x_1T,
\qquad
y(T)=z_0+z_1T+z_2T^2+z_3T^3.
\]

The cubic ordinate is necessary: after substitution, `R(x(T),T)` has degree
six.  Write

\[
R(x_0+x_1T,T)=\sum_{i=0}^6 f_iT^i.
\]

The universal leading coefficient identity, verified modulo `M`, is

\[
4f_6=(1-x_1^2)^2D,
\]

where

\[
\begin{aligned}
D={}&5c_1^4+4c_1^3-24c_1^2c_2-2c_1^2-16c_1c_2+16c_2^2
 +32c_1c_3\\
 &+4c_1+8c_2+32c_3-64c_4+5.
\end{aligned}
\]

Thus every nonvisible solution on the open chart
\((1-x_1^2)D\ne0\) lies over the double cover `D=w^2`.  Put

\[
s=z_3=\frac{(1-x_1^2)w}{2}.
\]

The top three coefficient equations give a triangular elimination:

\[
z_2=\frac{f_5}{2s},\qquad
N_1=4s^2f_4-f_5^2,
\qquad z_1=\frac{N_1}{8s^3},
\]

\[
N_0=8s^4f_3-N_1f_5,
\qquad z_0=\frac{N_0}{16s^5}.
\]

The remaining conditions, with denominators cleared only by powers of `s`,
are

\[
\begin{aligned}
E_2&=64s^6f_2-N_1^2-4N_0f_5,\\
E_1&=64s^8f_1-N_0N_1,\\
E_0&=256s^{10}f_0-N_0^2.
\end{aligned}
\]

All three are even in `w`.  Reducing by `w^2-D` therefore removes `w` and
all four ordinate coefficients.  On the double cover `w^2=D`, the resulting
abscissa conditions are compactly represented by

\[
(M,\overline E_2,\overline E_1,\overline E_0)
 \subset \mathbf Q[c_1,c_2,c_3,c_4,x_0,x_1].
\]

Its fully expanded residuals have total degrees `40, 50, 60` and respectively
`21707, 63285, 164807` terms.  The recursive representation is therefore the
canonical executable form.  Over `Q`, a rational solution of these projected
equations lifts to a rational ordinate only when `D` is a rational square;
the equation `w^2=D` must therefore remain part of any arithmetic search.

For two extra affine sections, introduce two pairs `(x0_j,x1_j)` and impose a
copy of the three residual equations for each, sharing `M` and `D`.  Away from
coincident sections and the excluded leading locus, the naive dimension count
is one: eight variables and seven equations.  No explicit two-extra
positive-dimensional component is asserted here.

## Explicit rational two-dimensional solution locus

Fermigier's six roots `alpha_i(u,v)` and his extra section

\[
x=A(u,v)+B(u,v)T,
\qquad
B=\frac{-u^2-v^2+2u+2v+1}{u^2+v^2+1}
\]

give an exact rational two-parameter solution locus.  The full source formulas
for all six `alpha_i` and `A` are transcribed in
`cas/verify_fermigier_affine_section_component.sing`.  That verifier sends
`alpha_2` to zero and `alpha_4` to one, constructs `q,g,R` independently, and
factors the specialized value exactly as

\[
R(x(T),T)=U(u,v)C(u,v,T)^2,
\]

where `C` is cubic and `U=S^2`, with

\[
S=\frac{
(u-1)(v-1)(v-u)(u+v+1)
(-u^2+uv+u+2v)(-uv+v^2-2u-v)
}{
(u^2+v^2+1)^2(2uv+u+v-1)^4
(-u^3+u^2v+u^2+uv+v^2-2u-v)^4
}.
\]

Thus `y=S*C` is an explicit cubic ordinate over `Q(u,v)`.  At `(u,v)=(3,5)`
the normalized roots other than zero and one are

\[
\left(\frac{157}{518},\frac{27}{74},
      \frac{1007}{1036},\frac{55}{1036}\right),
\]

which are the affine normalization of `(0,55,314,378,1007,1036)`.  The extra
section is

\[
x(T)=\frac{314}{1295}-\frac{17}{35}T
\]

with raw normalized ordinate coefficients

\[
(z_0,z_1,z_2,z_3)=
\left(
\frac{49747491}{41141582272},
-\frac{989134281}{97294282400},
-\frac{3107583}{328696900},
\frac{80028}{2220925}
\right).
\]

All eight original equations vanish exactly there.  Their Jacobian has rank
eight in the ten genuine variables, so this point is smooth of local dimension
two.  A separate `2x2` minor of the `(u,v)` parameterization is

\[
-\frac{21645237}{657393800}\ne0.
\]

This proves that the displayed formulas give a genuine positive-dimensional
solution locus, not an isolated specialization.

## Exact fiber at `(0,23,93,128,133,175)`

The newly found six-root tuple is a point of the root moduli space, not by
itself a positive-dimensional moduli component.  Normalize the first two
roots to zero and one.  Its exact coefficient point is

\[
(c_1,c_2,c_3,c_4)=
\left(-23,\frac{4489}{23},-\frac{8810207}{12167},
\frac{277065600}{279841}\right).
\]

It lies on `M=0`, and the leading invariant is already a rational square:

\[
D=\frac{207360000}{279841}=\left(\frac{14400}{529}\right)^2.
\]

Specializing the three ordinate-eliminated equations at this moduli point,
removing only the excluded visible factors `x1=+/-1`, and computing over
`Q[x1,x0]` gives the reduced degree-six fiber

\[
\begin{aligned}
14283x_0^3-214866x_0^2+882963x_0-940880&=0,\\
9261x_1^2-2760x_0+4579&=0.
\end{aligned}
\]

Equivalently, the intercept polynomial is

\[
(69x_0-619)(69x_0-304)(69x_0-115),
\]

and eliminating `x0` gives six simple rational slopes

\[
(21x_1-31)(21x_1+31)(21x_1-19)(21x_1+19)
(21x_1-1)(21x_1+1).
\]

Thus these are exactly the six nonvisible affine abscissae in the open chart,
not merely six points found by inspection.  The independent Singular verifier
reconstructs a rational cubic ordinate for each one from the triangular
recursion; it imports no ordinate formula from the rank-family code.  At all
six lifts, the Jacobian of the original eight equations has rank eight in ten
variables, so every lift is a smooth point of local dimension two.  This fiber
calculation is moduli geometry only and makes no independence claim.

### Distinctness from Fermigier's component

The six-root point above is not in the direct image or the Zariski closure of
Fermigier's two-parameter root family.  This is stronger than failing to solve
for one finite parameter pair: the exact audit also resolves the base points
of the rational parameterization.

For a centered monic sextic write its coefficients in degrees four through
zero as (K_4,K_3,K_2,K_1,K_0).  Comparing the Fermigier sextic with the
centered target

\[
 Z^6-11546Z^4+139944Z^3+30087529Z^2
 -807896712Z+777680784
\]

gives three scale-invariant equations from (K_2,K_3,K_0).  Factoring the
first and taking resultants of its two factors with the other two equations
shows that every possible finite parameter has

\[
 u^{24}(u-1)^{24}(u^2+u+1)^{48}=0.
\]

At (u=0,1), the common `v`-factor is (v^2(v-1)^2), and at
(u^2+u+1=0) it is ((v+u+1)^4).  In every case that common factor divides
(K_4): these are collapsed base points, not nondegenerate direct matches.

The exceptional directions are also exact.  Above the four rational finite
base points, every first nonzero centered sextuple consists of repeated pairs
and is never totally collapsed.  Above the conjugate quadratic base points,
the three target-invariant residuals have gcd one in the finite tangent
parameter; the missing vertical direction has residuals

\[
(73616305923,-111258226205248,904342225400),
\]

so it is not the target either.  On the projective parameter line at infinity,
the generic residual gcd is (L^4(L-1)^4).  Blowing up the two exceptional
directions (L=0,1), together with (L=\infty), again gives noncollapsed
sextuples with three repeated pairs.  This exhausts the boundary of the
parameter plane.  Since the six target roots are distinct, none of those
boundary sextuples is the target.

Consequently each smooth local two-dimensional incidence component through
the six lifts computed above is distinct from Fermigier's known
two-parameter component.  This is a component-membership statement in the
affine-section incidence geometry, not a new rank claim.

## Rank-at-least-13 family from `(0,25,95,143,168,205)`

The same corrected elimination produces a particularly useful second exact
fiber.  After normalizing the roots by `25`, its moduli point is

\[
(c_1,c_2,c_3,c_4)=
\left(-\frac{611}{25},\frac{136799}{625},-\frac{530557}{625},
\frac{18714696}{15625}\right),
\]

with

\[
D=\frac{278784}{625}=\left(\frac{528}{25}\right)^2.
\]

After removing only the visible slopes (x_1=\pm1), the open section fiber is
the reduced degree-six scheme

\[
\begin{aligned}
(23x_0+95)(115x_0-583)(575x_0-3444)&=0,\\
279841x_1^2+69000x_0-439201&=0.
\end{aligned}
\]

Thus all six nonvisible slopes are rational:

\[
x_1=\pm\frac{37}{23},\quad
\pm\frac{13}{23},\quad
\pm\frac7{23}.
\]

In the original integral-root coordinate their abscissae are

\[
\frac{-2375\pm37T}{23},\qquad
\frac{2915\pm13T}{23},\qquad
\frac{3444\pm7T}{23}.
\]

The verifier pins all six cubic ordinates.  For example, the first is

\[
\begin{aligned}
x&=\frac{-2375+37T}{23},\\
y&=8295400-\frac{71029947}{529}T
 +\frac{356162}{529}T^2-\frac{840}{529}T^3,
\end{aligned}
\]

and its square is the primitive Mestre quartic identically in `T`.

Now make the split-infinity base change

\[
T=\frac{39146-u^2}{2u}.
\]

At `u=197`, hence `T=337/394`, the thirteen displayed
visible-plus-infinity points have exact mod-3 span 12.  Adding the first
companion raises the span to 13, with pivot indices
`1,...,11,13,14`.  Any relation over `Q(u)` would specialize to a relation
here, so the base-changed family has generic Mordell--Weil rank at least 13.
Adding all six companions at this specialization still certifies only 13; no
generic or specialized upper bound is claimed.  Independently, the same
specialization has a pinned exact rank lower bound 17 and
`ln(N)=173.5948911449...` in
[the D-square family note](MESTRE_DSQUARE_FOUR_SCREEN.md).

The conductor frontier is compact.  In the direct parameter `T`, the
primitive discriminant core is squarefree and irreducible of degree 20.  Its
pullback to `u` is squarefree and irreducible of degree 40.  With coefficients
encoded as compact JSON integer lists in ascending order, their SHA-256 hashes
are respectively

```
fc36f00ad71a6b30126402aae310cdd2c9d35553e9f22910334c2ba4b9a05590
876a5e46a21c20cf531eb63469b55fe2cecf58d4fd5fdedfedacb3950a0e3a41
```

The degree-40 squarefree frontier is far leaner than the degree-398 frontier
in Kihara's rank-at-least-14 family.  This is a precise polynomial-degree
comparison, not by itself a conductor bound.

## Local two-section continuation at the rank-13 seed

The two labelled nonvisible sections

\[
x_1(T)=-\frac{95}{23}+\frac{37}{23}T,
\qquad
x_2(T)=\frac{583}{115}+\frac{13}{23}T
\]

at the preceding normalized six-root point give a useful local test case for
the two-section incidence.  A first-order evaluator reconstructs `R` by the
same monic-square recursion as above and propagates values and derivatives
only.  Thus it does not form any expanded `E_i` residual.

For the seven equations `M,E2_1,E1_1,E0_1,E2_2,E1_2,E0_2` in
`(c1,c2,c3,c4,x01,x11,x02,x12)`, the exact Jacobian at this labelled pair has
rank `6`, hence tangent-space dimension `2`, rather than the naively expected
one.  The same rank occurs at the clean primes `17`, `29`, and `37`.  This is
a tangent degeneracy, not by itself a proof of a two-dimensional component:
the seventh residual must still be retained at higher order.

A bounded lift does retain it.  At each of those primes, choose a nonsingular
six-by-six minor, enumerate the two free first-order digits, solve the six
selected equations, and reject a candidate unless all seven residuals vanish.
All `p^2` first-order choices survive to the next prime power; a nonzero
choice is then continued through `17^4`, `29^4`, and `37^4`.  These are finite
local computations, not an infinite Hensel theorem or a rational component
construction.  In characteristic zero, fixing
`x02=583/115+t` and `x12=13/23` and recursively solving the same six rows
gives a formal line through order `t^12`; the retained seventh row `E0_2`
vanishes at every checked coefficient.  Its exact quadratic obstruction after
solving the other six rows is identically zero on the two-dimensional tangent
plane, and four cubic samples spanning the cubic forms vanish as well.
More strongly, letting `x02=583/115+u` and `x12=13/23+v` and recursively
solving those same six rows gives the full bivariate formal germ through total
degree five: all twenty nonconstant coefficients of the retained `E0_2` row
are zero.  This remains a finite local calculation, not a proof that `E0_2`
is in the localized six-row ideal.
Conventional low-height rational reconstruction of the recorded nonzero
branches does not produce an exact rational solution; nor does a checked
Padé scan through that order find any nonconstant coordinate of numerator and
denominator degrees at most five.

The first intersection audit is deliberately limited.  The two selected
abscissae coincide only at `T=-529/60`; there the two reconstructed ordinates
are neither equal nor opposite, so their finite affine intersection number at
this seed is zero.  This neither computes an intersection at infinity nor any
reducible-fiber correction, and hence is not a Shioda height or Gram-matrix
calculation.

The excess is not repaired by relabelling the pair.  At this same root point,
the exact seven-by-eight Jacobian has rank `6` for each of the `15` unordered
pairs of the six nonvisible companions.  Thus the whole six-companion fiber
belongs to the same rank-six local phenomenon; it is not a missed
rank-seven/one-dimensional pair at this seed.  A promising continuation must
therefore move to another smooth multi-companion seed or identify this
two-dimensional germ globally before any rank-promotion attempt.
The independent six-companion seed `(0,23,93,128,133,175)` has the same
exact `15`-of-`15` rank-six pair audit.  Thus this is not an artefact of the
rank-13 D-square normalization; the naive seven-equation dimension count is
not a reliable guide at the presently known multi-companion fibers.

The existing finite-reduction certificate at `u=197` is also an explicit
non-promotion gate for this seed.  It proves rank at least `13` after adding
the first selected companion, and its run with all six companions (hence with
the selected pair) still exhibits only a rank-`13` independent subset.  This
does not prove a relation among the two companions or rule out rank `14`; it
does show that this fixed specialization cannot certify the desired new
direction by the present quotient set.  Any rational reconstruction of the
local germ needs a fresh specialization and a new finite-reduction escape
certificate.

For this particular D-square branch, the generic group-law calculation now
sharpens that last point.  After the split-infinity base change

\[
T=\frac{39146-u^2}{2u},
\]

write `P1` for the previously certified companion
`(-2375+37T)/23`, and let `V(r,+/-)` be the primitive visible section at
`x=r+/-T`.  Exact arithmetic over \(\mathbb Q(u)\), using the recursive
Mestre square root and covariant map, gives

\[
\begin{aligned}
P_2&=V(168,-)-V(168,+)-V(205,-)+V(205,+)-P_1,\\
P_3&=V(0,+)+V(143,-)+V(168,+)+V(205,-)+P_1,\\
P_4&=-V(0,-)-V(143,+)-V(168,+)-V(205,-)-P_1,\\
P_5&=-V(25,+)-V(95,-)-V(168,-)-V(205,+)+P_1,\\
P_6&=V(25,-)+V(95,+)+V(168,-)+V(205,+)-P_1.
\end{aligned}
\]

Thus all six displayed companions on the fixed-root D-square family lie in
the generic subgroup generated by the twelve visible sections, split
infinity, and `P1`; the other five cannot contribute a second generic
direction there.  This resolves the companion-independence question for this
branch, but does not identify the rank-six two-section *moduli* germ globally,
or compute saturation, heights, intersections, a Shioda Gram matrix, or a
rank upper bound.  The exact certificate is
`elliptic_mestre_dsquare_all_companion_generic_relations.json`.

In particular, this investigation does **not** establish a rational parameter
or a plane model, pair intersections, a Shioda Gram matrix, saturation,
independence from the existing generic rank-13 subgroup, or generic rank at
least 14.  The rank-six tangent degeneracy makes those later checks more
important, not less.

## A transverse rational two-section component through `(0,1,7,8,9,11)`

There is, however, a different six-root seed where the expected local
dimension does occur.  At

\[
(c_1,c_2,c_3,c_4)=(-35,455,-2605,5544),
\]

whose roots are `(0,1,7,8,9,11)`, the leading invariant is
`D=2304=48^2`.  The two affine abscissae

\[
x_1(T)=\frac{61+7T}{5},\qquad x_2(T)=\frac{33}{5}
\]

satisfy the recursive residuals.  Here the exact seven-by-eight Jacobian has
rank `7`, with a characteristic-zero transverse coordinate `c4`.  The same
rank persists modulo `17`, `19`, and `23`; regular lifts moving that free
coordinate by the prime have all seven recursive residuals zero through the
fourth prime power.  Solving the seven equations formally with
`c4=5544+t` recognizes the following component.  Put `z=c1+35` and
`d=30-z`:

\[
\begin{aligned}
c_1&=-35+z,\\
c_2&=455-\frac{77}{3}z+\frac{13}{36}z^2,\\
c_3&=-2605+\frac{652}{3}z-\frac{217}{36}z^2+\frac1{18}z^3,\\
c_4&=5544-608z+\frac{449}{18}z^2-\frac{49}{108}z^3+\frac1{324}z^4,\\
x_{01}&=\frac{z^2-66z+1098}{3d},\qquad
x_{11}=\frac{42-z}{d},\\
x_{02}&=\frac{z^2-69z+1188}{6d},\qquad x_{12}=0.
\end{aligned}
\]

The recursive checker evaluates `M,E2,E1,E0` for both sections at 301
admissible rational values of `z`.  After multiplication by `d^60`, the
respective numerator-degree bounds are `20,200,250,300`; the values therefore
certify all seven identities over `Q(z)`.  In particular,

\[
D=\frac{16}{9}(z-36)^2,
\]

and the same triangular recursion gives two rational cubic ordinates on the
open set `z != 30,36`.  No giant expanded residual is used in this check.

The component really lies in the split six-root base after one elementary
conic parametrization.  Its quartic factor is

\[
h(X)=\frac{(3X+z-33)(6X+z-42)
 (18X^2+9Xz-306X+z^2-72z+1296)}{324},
\]

and its discriminant has the factorization

\[
\operatorname{disc}(h)=
\frac{(z-36)^2(z-24)^2(z^2-36z+36)}{11664}.
\]

Writing `w^2=z^2-36z+36`, take

\[
z=\frac{12(r+3)}{1-r^2},\qquad w=6+rz.
\]

The six roots are then

\[
0,\ 1,\ \frac{33-z}{3},\ \frac{42-z}{6},\
\frac{102-3z-w}{12},\ \frac{102-3z+w}{12}.
\]

The exact verifier compares the product with `X(X-1)h(X)` at 38 rational
`r` values.  Clearing `(1-r^2)^6` gives coefficient numerators of degree at
most 12, so this is also an identity proof.  Thus this is an explicit
rational one-parameter two-section component on the split-six-root open
moduli locus, rather than merely a Hensel branch.

There is one exact finite intersection check, with its scope kept narrow.  At
the seed `z=0`, take the `D`-square root `+48` in the triangular recurrence.
The resulting primitive quartic ordinates are

\[
\begin{aligned}
y_1&=-196-\frac{2811}{25}T-\frac{434}{25}T^2-\frac{24}{25}T^3,\\
y_2&=-\frac{634}{25}T+T^3.
\end{aligned}
\]

Their squares are checked against the respective substituted quartics at
seven `T` values, hence identically.  The two abscissae coincide at `T=-4`,
where both signed sections pass through the finite point
`(x,y)=(33/5,936/25)`.  This records a finite meeting of this selected pair;
it is not the full section intersection number and says nothing about an
infinite-fibre contribution, reducible fibres, or a Shioda correction.

The rank audit currently gives a negative gate, not a promotion.  At
`r=8`, `T=3`, the twelve visible Jacobian images have combined exact mod-3
finite-reduction rank `9` using the displayed primes through `139`; adjoining
both affine points leaves the same rank, the same nine pivot indices, and the
same independent-subset hash.  This says only that these two points add no
direction to that recorded finite quotient span.  It is not a Mordell--Weil
relation, saturation calculation, intersection computation at infinity, or
Shioda Gram-matrix calculation.  In particular, this new component does
**not** yet prove generic rank at least `14`, and its height/intersection and
independence work remains open.

The first two-parameter quotient screen is also a non-promotion result.  For
every reduced `r` and `T` of numerator/denominator height at most five on the
rational component, 1,324 of the 1,330 parameter pairs admitted an exact
mod-3 certificate using primes at most 101.  The twelve visible points and
the augmented fourteen-point set had the same quotient rank in every such
case; the largest recorded rank was `9`.  Thus this explicit low-height grid
contains no specialization witness that the two affine sections are
independent of the visible subgroup.  The six uncertified pairs merely lacked
a rational-3-torsion exclusion within that deliberately short prime range.
This is a bounded negative screen, not a generic dependence theorem or a
Mordell--Weil upper bound.

One clean point on that grid gives a stronger exact non-promotion audit.  At
`r=2`, `T=1`, the short model is

\[
y^2=x^3-436481168886243x+1186790811178179337758.
\]

Let `V(a,+/-)` denote the Jacobian image of the visible quartic point at
`x=a+/-T`, retaining the displayed split-root labels rather than sorting the
roots.  Exact rational group-law computation gives

\[
\begin{aligned}
P_1&=V(0,-)-V(0,+)+V(1,+),\\
P_2&=-V(0,-)-V(1,+)-V((42-z)/6,-).
\end{aligned}
\]

So both affine points are literally in the visible subgroup at this fibre;
this is stronger than their images merely failing to enlarge a finite
quotient.  The visible and augmented exact mod-3 spans both have rank `9`.
At 64 and 128 decimal digits PARI/GP's height matrix has stable numerical
rank `9`, and `ellrank` returns bounds `[9,9]`.  These diagnostics support the
same rank-nine fibre picture, but do not prove saturation of the chosen
visible sublattice.

The generic relation is now proved.  Retaining the triangular-recursion signs
of the two cubic ordinates, exact covariant-map and group-law arithmetic over
\(\mathbb Q(r,T)\) gives

\[
P_1=-V(0,-)+V(0,+)-V(1,+),\qquad
P_2=V(0,-)+V(1,+)+V(r_3,-).
\]

Thus both selected affine sections are generically visible.  This first
rational component, like the later conic-rational component, is a genuine
positive-dimensional two-section locus but not a generic rank-jump locus.
The certificate is
`elliptic_mestre_transverse_component_generic_relations.json`; it does not
determine saturation, full intersections, or a Shioda Gram matrix.

### Bounded seed census beyond the initial families

To avoid extrapolating from the handful of hand-selected roots, the frozen
max-200 exact panel was screened with the same recursive projected equations.
Of its 803 new root tuples, 167 have square affine-normalized leading
invariant.  For each such seed the calculation exhausts the two abscissa
coordinates modulo the first usable prime in `7,11,13,17`, lifts every
nonsingular solution to precision eight, reconstructs it over `Q`, and checks
the three recursive equations exactly.  It then computes the exact
seven-by-eight Jacobian rank of every recovered pair.  This is a bounded
recovery procedure, not an absence proof for seeds with no recovered section.

Only one new seed passes that recovery gate:

\[
(0,7,79,81,128,137),
\]

with normalized moduli

\[
\left(-\frac{425}{7},\frac{66335}{49},-
\frac{4501495}{343},\frac{112212864}{2401}\right)
\quad\hbox{and}\quad D=\left(\frac{576}{7}\right)^2.
\]

It has the five reconstructed affine abscissae

\[
\left(\frac{137}{65},\pm\frac{47}{65}\right),
\quad\left(\frac{4932}{455},0\right),
\quad\left(\frac{9409}{455},\pm\frac{79}{65}\right).
\]

The central slope-zero section pairs transversely with each member of the
first or second signed pair: these are four exact rank-seven pairs.  The raw
`c4` and `x02` formal charts did not initially show a low-degree Padé model,
but their relation exposes the correct conic parameter.

### Second conic-rational component

Put `u=c1+425/7`.  The moduli coordinates are polynomial in `u`:

\[
\begin{aligned}
c_2&=\frac{66335}{49}-\frac{929}{21}u+\frac{13}{36}u^2,\\
c_3&=-\frac{4501495}{343}+\frac{93718}{147}u
      -\frac{2599}{252}u^2+\frac1{18}u^3,\\
c_4&=\frac{112212864}{2401}-\frac{1029264}{343}u
      +\frac{63671}{882}u^2-\frac{583}{756}u^3+\frac1{324}u^4.
\end{aligned}
\]

The residual root-splitting conic is

\[
w^2=49u^2-4284u+79524.
\]

On its rational parametrization

\[
u=-\frac{12(47s+357)}{s^2-49},\qquad w=282+su,
\]

the selected pair is rational:

\[
\begin{aligned}
x_{01}&=\frac{137s^2+1316s+3283}{65s^2+658s+1813},\\
x_{11}&=-\frac{47s^2+714s+2303}{65s^2+658s+1813},\\
x_{02}&=\frac{(4s+21)(9s+35)(137s^2+1316s+3283)}
 {7(s-7)(s+7)(65s^2+658s+1813)},\qquad x_{12}=0.
\end{aligned}
\]

Moreover

\[
D=\frac{16(7u-432)^2}{441},
\]

and the six split roots are

\[
0,\ 1,\ \frac{411-7u}{21},\ \frac{474-7u}{42},\
\frac{1254-21u-w}{84},\ \frac{1254-21u+w}{84}.
\]

The recursive verifier checks all seven residuals at 1,922 rational `s`
values.  A common-denominator numerator degree bound of 1,920 makes this an
identity certificate, and a separate 39-value degree-12 check certifies the
root product.  Thus this seed is now a second explicit rational
two-section component on the split six-root base, not merely a local branch.

At the seed `s=-357/47` the selected signed pair has a finite common point at
`T=3973/329`, namely

\[
\left(x,y\right)=
\left(\frac{4932}{455},-
\frac{107740485691272}{438652175}\right).
\]

As for the first component, this is only a finite-fibre intersection audit;
the full intersection number, Shioda corrections, saturation, and any
independence from a rank-13 subgroup are still open.  The bounded census and
the new exact component certificate are pinned in
`elliptic_mestre_two_section_seed_screen_max200.json` and
`elliptic_mestre_transverse_two_section_conic_component.json`, respectively.

The first-good-prime condition in that census is deliberately not a negative
test: a rational affine section can be singular in its two projected
abscissa coordinates at that prime.  A separate replay taking the union of
exact reconstructions at `7,11,13,17` recovers `38` rank-seven pairs at eight
base points.  Exact moduli comparison puts seven of those base points on the
previous rational two-section base curve.  The remaining eight-companion
point `(0,8,58,77,85,102)` is the exact affine normalization of Fermigier's
two-parameter roots at `(u,v)=(-3,-8/3)`, with `alpha_1 -> 0` and
`alpha_2 -> 1`.  Thus this stronger bounded replay finds no unclassified
rank-seven *base point* in the frozen panel.  It does not prove a generic
relation for every labelled pair at those points, or a height, saturation, or
independence statement.  The reproducible artifacts are
`elliptic_mestre_two_section_seed_screen_max200_all_primes.json` and
`elliptic_mestre_two_section_seed_screen_max200_all_primes_classification.json`.

At the Fermigier base point the first recovered affine line is exactly the
generic Fermigier extra line after that affine normalization.  Appending the
second line to the visible-plus-extra baseline at each `T=1,ldots,40` gives
no positive marginal dimension in exact mod-2 or mod-3 finite quotients
through prime `151`.  This is bounded non-escape evidence only: it neither
finds a dependence relation nor supplies saturation, heights, or a rank upper
bound.

### Fermigier local two-section branch

The remaining Fermigier base point has a smooth local continuation in the
two-parameter root surface.  Retaining the first (Fermigier) section and
writing the second normalized abscissa as `a_2+b_2 T`, the three recursive
residuals in `(u,v,a_2,b_2)` have rank three at
`(u,v)=(-3,-8/3)`.  The nonzero transverse determinant permits `u=-3+t` as
the formal parameter.  A compact square-root-recursion lift through `t^18`
recognizes

\[
 v=\frac{u^2+u+2}{u},\qquad
 b_2=\frac{u}{u^2+2},\qquad
 a_2=-\frac{2u^6-u^5+4u^4-u^3-8u^2-4u-16}
 {2(u+2)(u^2+2)(u^2-u+4)}.
\]

The degree-eight Padé search has independent holdout coefficients, and the
displayed formulas reproduce the compact recursive formal lift through order
18.  The follow-up exact verifier compares the source Fermigier chart and
then evaluates the recursive residuals at 1,141 admissible rational values.
Each coordinate denominator divides
`((u+2)(u^2+2)(u^2-u+4))^4`; with the established residual total-degree bound
60, the cleared residual numerator has degree at most 1,140.  Thus this is an
exact rational two-section component identity, without materializing an
expanded residual.  Its leading invariant is the square

\[
 D=\left(
 \frac{8u(u-2)(u-1)(u+1)(u^2+2)}
 {(u+2)(u^2-u+4)^2}
 \right)^2,
\]

so the triangular recursion supplies rational cubic ordinates on the stated
open locus.  The pinned discovery and identity artifacts are
`elliptic_mestre_fermigier_two_section_local_branch.json` and
`elliptic_mestre_fermigier_two_section_component.json`.

There is now an exact generic rank lower bound on this particular component.
Take the first eleven primitive visible points (discarding only the twelfth
visible point) and the two selected affine sections.  Thirteen fixed
finite-reduction quotients, distributed over the smooth specializations
`(u,T)=(-5,1),(-5,2),(-3,1),(-1/2,1)`, have full stacked column rank 13 over
`F_3`.  At `u=-5,T=1`, good reduction modulo 19 has group order 28, excluding
rational 3-torsion; smooth specialization also excludes generic rational
3-torsion.  Infinite 3-divisibility of a hypothetical generic relation then
proves these thirteen sections independent over `Q(u)(T)`.  Thus the
one-parameter Fermigier two-section curve has generic Mordell--Weil rank at
least 13.  The certificate uses the compact root/triangular-ordinate
representation throughout, not an expanded residual.

This is a genuine independent second direction on this component, but it is
not a rank-14 claim: it does not compare against a separate pre-existing
rank-13 family.  Saturation, the full pair intersections, and the Shioda Gram
matrix remain open.  The pinned certificate is
`elliptic_mestre_fermigier_two_section_generic_rank13.json`.

One signed finite intersection is now fixed before those lattice calculations:
at the continuation seed `u=-3`, the two cubic ordinates meet at

\[
 T=-\frac{479}{56},\qquad
 (x,y)=\left(\frac{445}{56},-\frac{1141635}{28}\right).
\]

Each cubic square identity is checked at seven distinct Mestre parameters,
which is exact for the degree-six specialized quartic.  At `T=1`, the exact
mod-3 finite-reduction quotient through prime 499 has rank nine both for the
twelve visible points and after adjoining this signed pair.  This records a
finite intersection and rejects this seed as an independence witness; it is
not the complete pair intersection number, a height calculation, saturation,
or a generic dependence proof.  The pinned audit is
`elliptic_mestre_fermigier_two_section_intersection.json`.

The component-level mod-2/mod-3 quotient grid now covers all 136 admissible
pairs with `H(u)<=5` and `T=1,...,4`, using reduction primes through 151.
The first Fermigier section raises the visible mod-3 quotient on 57 of those
fibres, but the reconstructed second section never raises the
visible-plus-first baseline for either modulus.  This is bounded non-escape
data only: the later stacked multi-specialization certificate shows that it
does not refute a generic independent second direction.  It remains neither a
generic relation, saturation, nor Shioda calculation.  The pinned record is
`elliptic_mestre_fermigier_two_section_component_escape_h5_t4.json`.

A separate numerical-height triage covers all 366 smooth specializations in
the rational panel `H(u)<=10`, `T=1,2,3`, always with the twelve visible and
two affine points.  Its maximum canonical-height-matrix rank is 13, reached
at five listed fibres (for example `u=-7/8,T=3`); no fibre has numerical rank
14.  This is only an efficiently reproducible target-selection result: a
canonical-height rank is neither an algebraic rank bound nor a generic
dependence statement.  In particular, it does not give a rank-14 witness;
the separate exact quotient certificate gives rank at least 13.  The pinned
record is
`elliptic_mestre_fermigier_two_section_height_triage_h10_t123.json`.

### A non-transverse conjugate-slope germ at diameter 233

The first-prime scout of the disjoint diameter band `231`--`235` has one
two-section seed even though it has no rank-seven pair.  Its roots are

\[
(0,7,127,128,225,233),
\]

and the two normalized affine abscissae are conjugate:

\[
x_1=\frac{233}{113}-\frac{97}{113}T,\qquad
x_2=\frac{233}{113}+\frac{97}{113}T.
\]

The full labelled recursive incidence Jacobian has exact rank six over
`Q`, and the same rank at the usable small primes `17,23,29`; it is therefore
not a naively transverse curve point.  Taking the second intercept and slope
as formal parameters, the two tangent directions preserve the conjugate form
`x_1=a-bT`, `x_2=a+bT`.  Solving six recursive rows in the bivariate formal
ring gives zero for the seventh row through total order four.  On the slice
varying the common intercept and fixing the slope, that residual remains zero
through order twelve.  This is local evidence for a two-dimensional formal
germ in the conjugate-slope locus.

One root-motion chart is already more informative than the original
intercept/slope coordinates.  Fix the largest moving root at `r_6=233/7` and
put `r_3=127/7+t`.  The six-row recursive implicit solve continues through
order twelve, keeps the slopes conjugate, and recognizes both affine
intercepts exactly as

\[
a=\frac{233}{113-7t}=\frac{233}{240-7r_3}.
\]

The other two roots and the common slope do not yet have a low-degree rational
recognition in this chart, so this is a coordinate-level local identity, not
a global reconstruction.

The pattern is not confined to that seed.  The disjoint next diameter band
contains `(0,21,151,169,200,239)` with the conjugate pair
`239/109-31T/109`, `239/109+31T/109`.  Its exact tangent rank is again six,
the rank persists at `11,17,19,23,29`, and its bivariate recursive lift has
zero remaining row through total order three.  This corroborates a recurring
conjugate-slope local locus, but still does not identify a common global
surface or establish any Mordell--Weil property.

There is now a useful seed-level Mordell--Weil audit, at `T=1`.  On the first
seed the two affine points are exactly

\[
P_1=-V_7+V_8-V_{12},\qquad P_2=-V_7+V_8+V_{11},
\]

and on the second they are exactly

\[
P_1=-V_7+V_8-V_{12},\qquad P_2=V_7-V_8-V_{11}.
\]

Here the `V_i` are the ordered primitive visible sections after mapping the
quartic to its short Jacobian.  The exact mod-3 quotient certificate through
prime `251` has rank nine both before and after adjoining the affine pair at
each seed.  Thus these two specializations provide no additional direction;
they are evidence against, not proof of, a rank jump on the formal germ.  In
particular, they neither supply a generic relation nor determine saturation,
heights, intersections, or Shioda data.  The reproducible audit is
`elliptic_mestre_conjugate_two_section_seed_relations.json`.

It has not yet been recognized as a rational surface or low-degree plane
model, nor reconstructed to a global identity.  No Mordell--Weil
independence, intersection, saturation, height, Shioda, or rank assertion is
made for it.  The reproducible local record is
`elliptic_mestre_conjugate_two_section_germ.json`.

The same all-prime method begins the post-census search without modifying the
frozen panel.  The diameter band `201`--`205` has no rank-seven pair; its sole
six-companion recovery is the already-known D-square germ, whose fifteen
pairs all have rank six.  The disjoint band `206`--`210` has 224
nonreflection square-leading candidates and no reconstructed affine section
at the declared primes and precision.  Both statements are bounded screens,
not an absence theorem for the full two-section locus.

The following band `211`--`215` does recover six rank-seven pairs at two
multi-companion bases.  Their normalized moduli lie on the rational component
above, and the exact kernel tangent of every pair agrees with the derivative
of that known component.  This rejects a new tangent direction at these two
seeds, but does not supply a generic relation for each additionally observed
specialized affine line.

The disjoint all-prime panel `226`--`230` has 26,556 obstruction-zero tuples
and 126 nonreflection square-leading candidates.  No candidate reconstructs
an affine section at precision eight over the union of `7,11,13,17`, hence it
contains no rank-seven pair in this declared method.  This is again a bounded
chart/primes/precision result, not an absence theorem for the two-section
locus.

### A smooth eight-companion seed at diameter 235

The next all-prime panel discovers a new local target at

\[
(0,17,136,161,207,235).
\]

It has eight reconstructed affine companions.  For the signed pair

\[
x_1=-\frac{68}{27}-\frac{37}{27}T,\qquad
x_2=\frac{3067}{459}-\frac{5}{27}T,
\]

the compact recursive incidence Jacobian has exact rank seven.  Taking `c1`
as a local parameter gives non-singular all-seven-row Hensel continuation
modulo `29^4` and `41^4`.  Over \(\mathbb Q[[t]]\), with
`c1=-739/17+t`, the recursive solve vanishes through order 12 and recognizes

\[
c_2=\frac{201815}{289}-\frac{1621}{51}t+\frac{13}{36}t^2.
\]

The moduli-coordinate projection was misleadingly complicated, but the
split-root chart is not.  Hensel-lifting the four nonfixed roots recognizes

\[
r_6=\frac{235}{17}-\frac{t}{3},\qquad c_1+3r_6=-\frac{34}{17}.
\]

Writing \(U=r_6-235/17\) and \(V=r_3-8\), the order-32 root jet gives a
sparse plane quartic.  Its three rational ordinary nodes (two affine and one
at infinity) make its normalization rational.  Projection from a node gives
a parameter `p`, with the local seed at `p=-294`.  The resulting expressions
for all four moduli and both affine lines are rational over \(\mathbb Q(p)\).
Their denominators divide

\[
\bigl((p-66)(p+54)(p^2+18p+456)\bigr)^4
  (3p^2+4p+1068).
\]

The six roots split rationally on this parameter line, and the leading
invariant is the square

\[
\left(
 \frac{15(p-26)(p-6)(p+6)(p+14)(3p^2+4p+1068)}
 {(p-66)(p+54)(p^2+18p+456)^2}
\right)^2.
\]

Most importantly, the exact verifier evaluates the compact recursive
residuals at 1,099 admissible rational `p` values.  The established residual
degree bound is 60, so after clearing the displayed degree-18 common
denominator, every numerator has degree at most 1,080.  This proves all seven
residual identities over \(\mathbb Q(p)\) without materializing an expanded
residual.  Thus this is an exact rational split-six-root two-section
component, not merely a local branch.  The certificate is
`elliptic_mestre_diameter235_eight_companion_component.json`; the earlier
order-32 low-bidegree record remains useful provenance for its discovery.

The two selected affine **abscissae** collide at
`T=-4223/544`, `x=4417/544`.  The compact triangular calculation is more
informative than the earlier square-only check: with the common leading
ordinate normalization, the two points have opposite raw-quartic ordinates.
The generic \(\mathbb Q(p)\) recursive calculation now verifies the same
opposite-ordinate identity at the unique affine-line crossing.  Thus this is
a hyperelliptic-conjugate collision, not a same-point affine intersection,
on the displayed component.  It neither supplies the full pair intersection
number nor rules out contributions at infinity or reducible fibres.

The finite-reduction audit is materially different from the previously known
multi-companion seeds: at `T=1,-1` the visible and every one-companion mod-3
quotient has rank 10, while at `T=2,3` every companion separately raises it
to 11 (and adjoining all eight still has rank 11).  At the regular smooth
specialization `p=-294,T=2`, exact finite-reduction certificates make the
first ten visible points together with either selected affine point
independent modulo 3.  A generic visible-subgroup equality would specialize,
so each selected affine section is outside the generic subgroup generated by
the twelve visible sections; the component therefore has generic rank at
least 11.  This does not separate the two new points from each other, prove
generic rank at least 12 or 14, establish saturation, or determine a Shioda
Gram matrix.  Full intersection/height work remains necessary.  The
reproducible local record is
`elliptic_mestre_diameter235_eight_companion_local.json`.

The same identity is now exact over the full field \(\mathbb Q(p,T)\):

\[
 P_1+P_2=V(0,-)+V(1,-)+V(1,+)+V(r_3,-)+V(r_3,+)
          +V(r_4,-)+V(r_5,+)+V(r_6,-),
\]

where the component ordering is `(0,1,r3,r4,r5,r6)` and the checker uses the
corresponding `(-,+)` visible ordering.  Thus \(P_2\) lies in the subgroup
generated by the visible sections and \(P_1\), so the displayed pair supplies
at most one quotient generator over the visible subgroup.  The same exact
calculation proves the two independent visible relations

\[
 V(1,-)+V(1,+)+V(r_3,+)+V(r_4,+)+V(r_5,+)+V(r_6,-)=0,
\]
\[
 V(0,-)+V(0,+)+V(r_3,-)+V(r_4,-)+V(r_5,-)+V(r_6,+)=0.
\]

Thus the twelve visible sections have rank at most 10, while the regular
finite-reduction certificate above makes the first ten visible sections and
either affine section independent.  Consequently the subgroup generated by
the twelve visible sections and the displayed pair has **exact rank 11**
over \(\mathbb Q(p,T)\).  This is a subgroup statement, not an upper bound
for the full Mordell--Weil group: it does not prove generic rank 12 or 14,
saturation, a full intersection calculation, or a Shioda Gram matrix.  The
generic certificate is
`elliptic_mestre_diameter235_eight_companion_generic_relation.json`; the
seed-only checkpoint remains
`elliptic_mestre_diameter235_eight_companion_seed_pair_relation.json`.

The first ten visible sections plus (P_1) have also been audited at the
regular seed fibre `p=-294,T=2`.  The real height matrix has numerical rank
11 at both 72 and 120 decimal digits (smallest eigenvalue about `2.045`).  A
single PARI `ellsaturation` call through primes below 20 returns eleven exact
points and reports a height-determinant ratio

\[
1048576=2^{20}.
\]

Thus the displayed independent seed basis is visibly not a saturated basis
in this bounded computational audit (the corresponding determinant heuristic
is index (2^{10})).  This is precisely the kind of warning that prevents a
rank statement from being mistaken for a Mordell--Weil lattice statement.  It
does **not** establish a saturated basis—at the seed or generically—because
the PARI routine is documented under a finite-index hypothesis and this audit
does not determine the full fibre rank.  Nor does it supply the missing full
intersection calculation, reducible-fibre contributions, or Shioda Gram
matrix.  The replayable record is
`elliptic_mestre_diameter235_displayed_lattice_seed_audit.json`.

At the seed specialization `s=-357/47,T=1`, exact group law puts both affine
points in the visible subgroup:

\[
P_1=-V_1-V_4-V_5-V_7-V_{10},\qquad P_2=-V_1-V_4-V_5.
\]

The `V_i` denote the fixed ordered visible points of that audit, so these are
not labelled generic relations.  The mod-3 finite-reduction quotient has
rank nine before and after adjoining the pair.  A separate height-five screen
in the conic parameter `s` and Mestre parameter `T` covers all 1,482
candidate pairs; each is admissible and none raises that rank-nine quotient.
Both are negative bounded-search results only, not a saturation calculation,
Shioda Gram matrix, or generic rank upper bound.  The reproducible artifacts
are `elliptic_mestre_transverse_two_section_conic_component_seed_relations.json`
and `elliptic_mestre_transverse_two_section_conic_component_independence_h5.json`.

There is one stronger generic result on this component.  Keeping the
triangular-recursion sign of the first selected cubic ordinate, exact
covariant-map and short-Weierstrass arithmetic over \(\mathbb Q(s,T)\) proves

\[
P_1=V(0,-)+V(1,+)+V(r_3,-)+V(r_4,-)+V(r_5,+).
\]

Thus the first selected affine section cannot supply a new generic direction.
This uses only the compact recursive square-root representation, not an
expanded two-section residual.  The same exact calculation gives the shorter
second relation

\[
P_2=-V(0,-)-V(1,+)-V(r_3,-).
\]

Consequently both selected affine sections are already in the generic visible
subgroup on this component: it is a genuine positive-dimensional two-section
locus, but it gives no new Mordell--Weil direction.  This does not determine
the visible subgroup's rank or saturation, nor its Shioda data.  The combined
certificate is
`elliptic_mestre_transverse_two_section_conic_component_generic_relations.json`;
the earlier first-section-only artifact remains a smaller replay checkpoint.

## Reproduction

From the repository root:

```bash
Singular -q archive/elliptic-curves/cas/mestre_affine_section_elimination.sing
Singular -q archive/elliptic-curves/cas/verify_fermigier_affine_section_component.sing
Singular -q archive/elliptic-curves/cas/verify_fermigier_affine_section_jacobian.sing
Singular -q archive/elliptic-curves/cas/verify_mestre_02393128133175_moduli_fiber.sing
Singular -q archive/elliptic-curves/cas/audit_fermigier_component_mestre_02393128133175.sing
Singular -q archive/elliptic-curves/cas/verify_mestre_02595143168205_moduli_fiber.sing
PYTHONPATH=elliptic-curves/cas python3 elliptic-curves/cas/verify_mestre_02595143168205_rank13_section.py
PYTHONPATH=elliptic-curves/cas python3 elliptic-curves/cas/verify_mestre_02595143168205_discriminants.py
python3 elliptic-curves/cas/probe_mestre_two_section_local_continuation.py \
  --precision 4 --bivariate-order 5 \
  --output archive/elliptic-curves/artifacts/generated-results/elliptic_mestre_two_section_local_continuation.json
python3 archive/elliptic-curves/cas/probe_mestre_transverse_two_section.py --precision 4
PYTHONPATH=elliptic-curves/cas python3 \
  archive/elliptic-curves/cas/verify_mestre_transverse_two_section_component.py
PYTHONPATH=elliptic-curves/cas python3 \
  archive/elliptic-curves/cas/audit_mestre_transverse_two_section_specialization.py
PYTHONPATH=elliptic-curves/cas python3 \
  archive/elliptic-curves/cas/verify_mestre_transverse_two_section_conic_component.py
PYTHONPATH=elliptic-curves/cas python3 \
  archive/elliptic-curves/cas/audit_mestre_transverse_conic_component_relations.py
PYTHONPATH=elliptic-curves/cas .venv/bin/python \
  archive/elliptic-curves/cas/verify_mestre_transverse_conic_component_generic_relations.py \
  --include-second
PYTHONPATH=elliptic-curves/cas python3 \
  archive/elliptic-curves/cas/screen_mestre_transverse_conic_component_independence.py \
  --root-height 5 --parameter-height 5 --prime-bound 101
PYTHONPATH=elliptic-curves/cas python3 \
  elliptic-curves/cas/screen_mestre_two_section_transverse_seeds.py --output \
  archive/elliptic-curves/artifacts/generated-results/elliptic_mestre_two_section_seed_screen_max200.json
PYTHONPATH=elliptic-curves/cas python3 \
  elliptic-curves/cas/screen_mestre_two_section_transverse_seeds.py --all-primes --output \
  archive/elliptic-curves/artifacts/generated-results/elliptic_mestre_two_section_seed_screen_max200_all_primes.json
PYTHONPATH=elliptic-curves/cas python3 \
  archive/elliptic-curves/cas/audit_mestre_multprime_seed_classification.py --output \
  archive/elliptic-curves/artifacts/generated-results/elliptic_mestre_two_section_seed_screen_max200_all_primes_classification.json
python3 -m unittest elliptic-curves/tests/test_mestre_02595143168205_rank13.py
```

The first command derives and eliminates symbolically.  The second verifies
the exact two-parameter square identity and its rational square unit.  The
third independently rebuilds the original coefficient ideal and proves the
rank-eight Jacobian statement at `(u,v)=(3,5)`.  The fourth reconstructs the
new six-root point independently and proves its exact rational degree-six
affine-section fiber.  The fifth performs the finite-image and compactified
base-point audit proving distinctness from Fermigier's component.  The last
four commands verify the second rational fiber, its generic rank-13
specialization certificate, the degree-20/40 discriminant geometry, and their
focused test replay.
