# Rank-three collision-framed descent audit

The rank-three collision sheet removes the finite labeling ambiguity in the
universal cubic Keller fiber.  More precisely, the off-diagonal ordered-pair
cover of a rank-three finite-etale algebra is already its full `S_3` frame
torsor.  On that torsor, two primitive cubic presentations determine a unique
projective change of root coordinate, and uniqueness makes its descent
cocycle exact.

For the foundational cubic Keller map, this projective transition lifts
canonically to the normalized linear--quadratic factorization map after
localizing the target at one explicit linear denominator.  The resulting
target-localized, or germ-level, Keller incidence therefore descends through
arbitrary cubic primitive-coordinate changes.

This does **not** prove a choice-free global polynomial map
\(BS_3\to\mathscr I_3\).  Within the canonical factorization transport, the
only projective transitions that act polynomially on the whole affine chart
are the known scaling torus.  Arbitrary translations and quadratic
Tschirnhaus changes meet a genuine chart-boundary denominator.  No claim is
made here that every possible nonlinear polynomial self-equivalence of the
foundational map has been classified.

Work over a characteristic-zero base.  All assertions remain valid over a
base on which the displayed discriminants and determinants are units.

## 1. The collision sheet is the cubic frame torsor

Let `E -> S` be finite etale of rank three.  Following the
collision-algebra interface adopted in the
[external audit](COLLISION_IDEALS_EXTERNAL_AUDIT.md), define

\[
 \operatorname{Off}_2(E/S)
 =(E\times_SE)\setminus\Delta_E.                       \tag{1.1}
\]

The tensor-product realization of the ordered self-collision fiber and its
diagonal/off-diagonal splitting are credited to Chloe van der Vlugt's
*Collision Ideals and Off-Diagonal Sheets*.  The statement below that the
rank-three factor is the frame torsor, and its use for Keller-presentation
descent, are deductions in this repository.

A geometric point of (1.1) is an ordered pair `(r_i,r_j)` with `i!=j`.
There is a unique missing third point `r_k`.  Scheme-theoretically, two
disjoint sections of a rank-three finite-etale cover have a rank-one
open-and-closed complement, hence a third section.  This construction is
functorial under base change and gives

\[
\boxed{
 \operatorname{Off}_2(E/S)
 \simeq
 \operatorname{Isom}_S(\{1,2,3\}_S,E).}               \tag{1.2}
\]

Consequently `Off_2(E/S)` is not merely a rank-six finite-etale cover.  It is
canonically the right `S_3`-torsor of complete orderings of `E`.  Equivalently,
for rank three,

\[
 \operatorname{Conf}_2(E/S)\simeq\operatorname{Conf}_3(E/S). \tag{1.3}
\]

This strengthens the generic normal-closure statement in the
[universal relative theorem](UNIVERSAL_RELATIVE_KELLER_MAP.md): (1.2) is
valid for every rank-three finite-etale cover, including split and
disconnected fibers, without a generic connectedness hypothesis.

## 2. Exact projective transition between framed cubics

On the frame torsor, let `(r_1,r_2,r_3)` and `(u_1,u_2,u_3)` be the ordered
roots of two primitive presentations of the same algebra.  Form the matrix

\[
 M(r,u)=
 \begin{pmatrix}
 r_1&1&-u_1r_1&-u_1\\
 r_2&1&-u_2r_2&-u_2\\
 r_3&1&-u_3r_3&-u_3
 \end{pmatrix}.                                      \tag{2.1}
\]

Let `m_j=(-1)^j det(M_{\widehat j})`, with columns numbered from zero, and
put

\[
 g(r,u)=
 \begin{pmatrix}m_0&m_1\\m_2&m_3\end{pmatrix}.        \tag{2.2}
\]

The maximal-minor identity `M(m_0,m_1,m_2,m_3)^t=0` gives

\[
 \frac{m_0r_i+m_1}{m_2r_i+m_3}=u_i
 \qquad(i=1,2,3).                                    \tag{2.3}
\]

Moreover,

\[
\boxed{
 \det g(r,u)
 =
 \prod_{i<j}(r_i-r_j)\prod_{i<j}(u_i-u_j).}           \tag{2.4}
\]

Both Vandermonde factors are units on the ordered squarefree cover, so (2.2)
defines an element of `PGL_2`.  It is the unique projective transformation
sending the first ordered triple to the second.

A simultaneous permutation of the two triples multiplies all four maximal
minors by the sign of the permutation.  It therefore leaves the class of
`g(r,u)` in `PGL_2` unchanged.  Effective descent along the `S_3`-torsor
(1.2) gives a projective transition already on the presentation overlap.

For a third presentation `v`, uniqueness gives

\[
\boxed{g(u,v)g(r,u)=g(r,v)\quad\text{in }PGL_2.}       \tag{2.5}
\]

Thus the finite collision-framing ambiguity and the projective-coordinate
cocycle both close exactly in rank three.

## 3. Quadratic Tschirnhaus formulas and their two boundaries

The projective description can be compared directly with primitive-element
coordinates.  Write

\[
 E_{\pi,b,c}(R)=\pi R^3+bR^2+R-\frac c2
 =\pi\prod_{i=1}^3(R-r_i),                            \tag{3.1}
\]

and let `e_1,e_2,e_3` be the elementary symmetric functions of the roots.
The normalization of the linear coefficient says

\[
 \pi=\frac1{e_2},\qquad
 b=-\frac{e_1}{e_2},\qquad
 c=\frac{2e_3}{e_2}.                                  \tag{3.2}
\]

Every second generator has a unique expression

\[
 u=q_0+q_1r+q_2r^2.                                   \tag{3.3}
\]

First put

\[
\begin{aligned}
 A_1={}&q_1e_1+q_2(e_1^2-2e_2),\\
 A_2={}&q_1^2e_2+q_1q_2(e_1e_2-3e_3)
          +q_2^2(e_2^2-2e_1e_3),\\
 A_3={}&e_3(q_1^3+q_1^2q_2e_1+q_1q_2^2e_2+q_2^3e_3).
                                                               \tag{3.4}
\end{aligned}
\]

The elementary symmetric functions of the three conjugates
`u_i=q_0+q_1r_i+q_2r_i^2` are

\[
\begin{aligned}
 f_1&=A_1+3q_0,\\
 f_2&=A_2+2q_0A_1+3q_0^2,\\
 f_3&=A_3+q_0A_2+q_0^2A_1+q_0^3.                    \tag{3.5}
\end{aligned}
\]

Whenever `f_2` is a unit, the normalized target of the second presentation is

\[
\boxed{
 \pi'=\frac1{f_2},\qquad
 b'=-\frac{f_1}{f_2},\qquad
 c'=\frac{2f_3}{f_2}.}                                \tag{3.6}
\]

There are exactly two transition boundaries.  The basis-change determinant
from `(1,r,r^2)` to `(1,u,u^2)` is

\[
\begin{aligned}
 \Theta={}&q_1^3+2q_1^2q_2e_1
 +q_1q_2^2(e_1^2+e_2)
 +q_2^3(e_1e_2-e_3),                                 \tag{3.7}
\end{aligned}
\]

because

\[
 \prod_{i<j}(u_i-u_j)
 =\Theta\prod_{i<j}(r_i-r_j).                         \tag{3.8}
\]

Thus `Theta` is a unit exactly when `u` remains primitive.  The second
factor `f_2` is a unit exactly when its relation can be normalized to have
linear coefficient one.  On the overlap

\[
\boxed{f_2\Theta\ne0,}                                 \tag{3.9}
\]

and the discriminants satisfy

\[
\boxed{
 \Delta(E_{\pi',b',c'})
 =
 \Delta(E_{\pi,b,c})
 \left(\frac{\pi'}{\pi}\right)^4\Theta^2.}             \tag{3.10}
\]

The coefficient cocycle is also explicit.  If
`v=p_0+p_1u+p_2u^2`, reduce `u^2` modulo the original cubic and put

\[
\begin{aligned}
 K_0={}&q_0^2+2q_1q_2e_3+q_2^2e_1e_3,\\
 K_1={}&2q_0q_1-2q_1q_2e_2+q_2^2(-e_1e_2+e_3),\\
 K_2={}&2q_0q_2+q_1^2+2q_1q_2e_1+q_2^2(e_1^2-e_2).
                                                               \tag{3.11}
\end{aligned}
\]

Then the composite is represented in the original basis by

\[
\boxed{
 (p_0+p_1q_0+p_2K_0,\;
  p_1q_1+p_2K_1,\;
  p_1q_2+p_2K_2).}                                    \tag{3.12}
\]

Substitution of (3.12) in (3.4)--(3.6) equals the iterated target
transformation.  Hence the Tschirnhaus presentation groupoid itself has no
remaining rank-three cocycle defect.

## 4. Application to the foundational Keller map

Let `ell` be the coefficient functional

\[
 \ell(C)=[UV^2]C
\]

on binary cubics, and define

\[
\begin{aligned}
 H_\ell&=\{C\in\operatorname{Sym}^3(k^2):\ell(C)=1\},\\
 X_\ell&=\{(L,Q):
   \operatorname{Res}(L,Q)=1,\ \ell(LQ)=1\}.
                                                               \tag{4.1}
\end{aligned}
\]

Multiplication

\[
 \mu_\ell:X_\ell\longrightarrow H_\ell,\qquad(L,Q)\longmapsto LQ
                                                               \tag{4.2}
\]

is the normalized linear--quadratic factorization model.  After the explicit
polynomial source isomorphism and fixed linear target normalization in the
[factorization theorem](NORMALIZED_FACTORIZATION_MODEL.md), it is the
foundational cubic map, equivalently the fixed rank-three map
`\mathcal K_3` in the
[universal relative theorem](UNIVERSAL_RELATIVE_KELLER_MAP.md).

Use the substitution convention `(g dot C)(U,V)=C(g(U,V))`.  This sends the
root divisor by the inverse projective transformation.  For
`g in GL_2`, put

\[
 \delta=\det g,\qquad D_g(C)=\ell(g\mathbin{\cdot}C).  \tag{4.3}
\]

On the target open `D_g!=0`, define

\[
\begin{aligned}
 C'&=\frac{g\mathbin{\cdot}C}{D_g(C)},\\
 L'&=\frac{D_g(C)}{\delta^2}(g\mathbin{\cdot}L),\\
 Q'&=\frac{\delta^2}{D_g(C)^2}(g\mathbin{\cdot}Q).
                                                               \tag{4.4}
\end{aligned}
\]

The binary resultant covariance

\[
 \operatorname{Res}(g\mathbin{\cdot}L,g\mathbin{\cdot}Q)
 =\delta^2\operatorname{Res}(L,Q)                    \tag{4.5}
\]

gives

\[
 \operatorname{Res}(L',Q')=1,\qquad
 L'Q'=C',\qquad\ell(C')=1.                            \tag{4.6}
\]

Changing the scalar lift of `g` cancels from every formula in (4.4), so this
is an intrinsic `PGL_2` transport.  For a second matrix `h`,

\[
 D_h\left(\frac{g\mathbin{\cdot}C}{D_g(C)}\right)
 =\frac{D_{gh}(C)}{D_g(C)},                           \tag{4.7}
\]

and the source scaling factors in (4.4) multiply to the factor for `gh`.
Here the order follows the substitution convention:
`h dot (g dot C)=(gh) dot C`.
Thus (4.4) satisfies the same exact cocycle as (2.5).

Given two normalized primitive presentations of one rank-three finite-etale
algebra, take the inverse of the projective matrix from Section 2 to match
the substitution convention in (4.3).  Equations (4.4)--(4.7) give an
isomorphism between the corresponding restrictions of the foundational map
to target neighborhoods containing the selected finite-etale fibers.

Therefore:

\[
\boxed{\text{the rank-three Keller incidence descends after target
localization at }D_g.}                                 \tag{4.8}
\]

In particular, the abstract cubic algebra, its full `S_3` collision frame,
and its marked-root Keller **germ** are compatible under arbitrary primitive
coordinate change.  The collision sheet supplies exactly the finite descent
data required for this statement.

## 5. The remaining global affine-space boundary

Formula (4.4) divides by the affine-linear target function `D_g(C)`.  It
therefore need not define a polynomial self-equivalence of the whole
`A^3 -> A^3` map.

In fact, within (4.4), polynomiality of the normalized target map forces
`D_g` to be constant.  If `D_g` were a nonconstant linear polynomial on
`H_ell`, each other coefficient of `g dot C` is also affine-linear.
Divisibility of all three normalized coordinates by `D_g` would make their
quotients constant.  The resulting target map would be constant, contradicting
invertibility of `Sym^3(g)`.  Thus no cancellation hidden in the target
coordinate ring enlarges the constant-denominator locus.

For

\[
 g=\begin{pmatrix}a&b\\c&d\end{pmatrix},\qquad
 C=\pi U^3+\beta U^2V+UV^2+hV^3,
\]

direct coefficient extraction gives

\[
\begin{aligned}
 D_g(C)={}&3\pi ab^2+\beta(2abd+b^2c)\\
           &+(ad^2+2bcd)+3hcd^2.                     \tag{5.1}
\end{aligned}
\]

This denominator is constant on `H_\ell` exactly when

\[
 ab^2=0,\qquad b(2ad+bc)=0,\qquad cd^2=0.             \tag{5.2}
\]

After saturating by `ad-bc`, the ideal in (5.2) is exactly

\[
\boxed{(b,c).}                                        \tag{5.3}
\]

Hence the projective stabilizer of the normalized tangent hyperplane is the
diagonal torus.  Writing `alpha=d/a`, its target action is

\[
\boxed{
 (\pi,b,c)\longmapsto
 (\alpha^{-2}\pi,\alpha^{-1}b,\alpha c).}              \tag{5.4}
\]

On the explicit fixed cubic Keller map this is realized globally by

\[
 (x,y,z)\longmapsto
 (\alpha x,\alpha^{-1}y,\alpha^{-2}z).                 \tag{5.5}
\]

Thus the affine scaling part of primitive-coordinate descent is an honest
global polynomial source--target equivariance.  General translation and
quadratic Tschirnhaus transitions require the divisor `D_g=0` in the
canonical transport.

This proves a boundary statement, not a global impossibility theorem:

- the finite `S_3` inertia and its cocycle are completely resolved;
- the full transition acts on the target-localized Keller map;
- within the canonical projective-root transport, only the torus is
  denominator-free on the entire affine chart;
- a different nonlinear polynomial automorphism, stabilization, or affine
  modification that cancels `D_g` has not been excluded.

Accordingly, the global rank-three arrow
\(BS_3\to\mathscr I_3\) from the universal-relative note remains open.

## 6. What this buys and the next focus

The result has three immediate uses.

1. Any construction depending only on the finite-etale cubic fiber, its
   ordered collisions, or the target-local marked-root germ is now
   presentation-independent.  Arithmetic decomposition, Frobenius data,
   normal closure, and local deformation calculations can be transported
   without choosing a preferred primitive generator.
2. A future global rank-three descent proof no longer needs to analyze
   `S_3` inertia or a mysterious Tschirnhaus cocycle.  Its exact remaining
   task is to cancel, absorb, or prove unavoidable the divisor `D_g=0` while
   retaining an affine-space polynomial Keller map.
3. Rank four has a genuinely new obstruction.  Its ordered-pair sheet has
   rank twelve and is not the full frame torsor; `Conf_3` has rank twenty-four
   and supplies the `S_4` frame.  The
   [rank-four continuation](RANK_FOUR_COLLISION_CROSS_RATIO.md) now computes
   the surviving cross-ratio defect exactly.  Canonical projective transport
   covers its zero hypersurface, while the complement requires a genuinely
   nonprojective Keller lift.

The
[all-rank continuation](ALL_RANK_COLLISION_PROJECTIVE_DESCENT.md) explains
this jump uniformly: the projective locus has codimension `N-3`, so the
cubic condition is automatic while the quartic has its first single
equation.

The recommended order is therefore:

1. test whether an affine modification or stable polynomial suspension
   cancels `D_g` in rank three;
2. attempt the analogous target-local projective lift on the rank-four
   cross-ratio hypersurface;
3. test the explicit primitive quadratic transition off that hypersurface
   against stable boundary and marked-normalization invariants.

## 7. Exact regression

Run

```bash
.venv/bin/python scripts/verify_rank_three_collision_descent.py
```

The checker verifies symbolically:

- all six cubic frames and the simply transitive `S_3` action;
- the signed-minor projective interpolation formula, Vandermonde
  determinant, simultaneous relabeling invariance, and three-presentation
  cocycle;
- the full quadratic Tschirnhaus symmetric formulas, primitive and
  normalization boundaries, discriminant transform, and coefficient
  composition law;
- resultant covariance, scalar-lift independence, and the normalized
  factorization transport;
- the saturated constant-denominator stabilizer ideal `(b,c)`; and
- the exact source--target torus equivariance of the fixed cubic Keller map.

These are exact polynomial or rational identities.  The checker does not
claim to classify nonlinear polynomial self-equivalences outside the
canonical factorization transport.
