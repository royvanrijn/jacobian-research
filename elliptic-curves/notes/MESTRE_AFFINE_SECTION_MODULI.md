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
\left(-\frac{636}{25},\frac{152074}{625},-\frac{667356}{625},
\frac{31978621}{15625}\right),
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

## Reproduction

From the repository root:

```bash
Singular -q elliptic-curves/cas/mestre_affine_section_elimination.sing
Singular -q elliptic-curves/cas/verify_fermigier_affine_section_component.sing
Singular -q elliptic-curves/cas/verify_fermigier_affine_section_jacobian.sing
Singular -q elliptic-curves/cas/verify_mestre_02393128133175_moduli_fiber.sing
Singular -q elliptic-curves/cas/audit_fermigier_component_mestre_02393128133175.sing
Singular -q elliptic-curves/cas/verify_mestre_02595143168205_moduli_fiber.sing
PYTHONPATH=elliptic-curves/cas python3 elliptic-curves/cas/verify_mestre_02595143168205_rank13_section.py
PYTHONPATH=elliptic-curves/cas python3 elliptic-curves/cas/verify_mestre_02595143168205_discriminants.py
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
