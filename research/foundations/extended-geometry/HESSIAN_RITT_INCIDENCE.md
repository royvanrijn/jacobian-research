# Hessian incidence for vertical Ritt decompositions

Work over an algebraically closed field of characteristic zero.  Let

\[
 H\in\mathcal A_N,
 \qquad f_s(W)=H(W)-sW,
 \qquad K_H(W)=H''(W).
\]

This note replaces coefficient elimination in all variables of a vertical
decomposition by an incidence problem for the projective Hessian polynomial.
It also gives the complete first composite case, degree six.

## 1. The Hessian reduction

For a factor pair `ab=N`, define the projective Hessian-composition incidence

\[
 \mathcal C_{a,b}=
 \overline{\left\{
  [\operatorname{div}((A\circ B)'')]:
  \deg A=a,\ \deg B=b
 \right\}}
 \subset \mathcal Q_{N-2},                         \tag{1.1}
\]

where `Q_(N-2)` is the marked-Hessian quotient by `W -> lambda W`.  The
divisor retains multiplicities.

> **Hessian reduction lemma.**  On the exact-degree locus,
> \[
>  H-sW=A\circ B\quad\hbox{for some }s,A,B
> \]
> if and only if the marked projective Hessian of `H` belongs to
> `C_(a,b)`.

Indeed, one implication is obtained by differentiating twice.  Conversely,
after undoing the source dilation and projective scalar in (1.1), equality of
Hessians gives

\[
 H(W)=\lambda(A\circ B)(\mu W)+uW+v.                \tag{1.2}
\]

The constant `v` is absorbed into the outer polynomial and `s=u`.  Thus the
two integrations introduce exactly the vertical parameter and an irrelevant
target translation; they introduce no further condition.

This also proves the clean finite-over-incidence statement suggested by the
marked normalization.  Let

\[
 \mathfrak h_N:\mathcal A_N^{\rm hc}\longrightarrow\mathcal Q_{N-2}
\]

be the Hessian map on the D1 Hessian-clean open, shrunk as in D1/F2.  Then

\[
 \boxed{\mathcal D_{a,b}^{\rm hc}
 =\mathfrak h_N^{-1}(\mathcal C_{a,b})}.             \tag{1.3}
\]

The restriction of `h_N` to (1.3) is finite etale of degree `N-2` over its
image: vertical decomposability is preserved by every D1 rerooting.  After
adjoining the distinguished affine root sheet, F2 identifies the source with
the normalization of the marked image.  Hence reconstruction of `s,A,B` is
needed only after the divisor incidence test has passed.

The parameter count is now transparent.  Pairs `(A,B)` have dimension
`a+b+2`.  Quotienting the two-dimensional affine change of the intermediate
coordinate and the two-dimensional affine postcomposition leaves
`a+b-2` projective Hessian parameters.  Quotienting the marked source dilation
gives the unconditional bounds

\[
 \dim\mathcal C_{a,b}\le a+b-3,
 \qquad
 \operatorname{codim}_{\mathcal Q_{N-2}}\mathcal C_{a,b}
 \ge N-a-b=(a-1)(b-1)-1.                            \tag{1.4}
\]

In particular, the right side is positive for every factor pair of every
composite `N>=6`.  Since there are finitely many factor pairs, a generic seed
of every such degree has no imprimitive vertical specialization.  The quartic
pair `(2,2)` is the unique zero-codimension exception, agreeing with the
ubiquity described in the automorphism note.

Equality in (1.4) is the expected-codimension assertion.  The canonical
coefficient reconstruction in Section 6 proves it for every factor pair in
every degree.  The low-degree calculations below retain the explicit
incidence equations and intersection geometry.

## 2. The two degree-six Hessian hypersurfaces

Write

\[
 K(W)=k_4W^4+k_3W^3+k_2W^2+k_1W+k_0,
 \qquad k_4\ne0.                                    \tag{2.1}
\]

### The `3 o 2` locus

An inner quadratic is affine-equivalent to `(W-c)^2`.  Therefore a
`3 o 2` composition has Hessian even about `c`.  Conversely, integrating an
even quartic twice, with zero first derivative at its center, gives a cubic
in `(W-c)^2`.  The center is `c=-k_3/(4k_4)`, and the remaining odd
coefficient is zero exactly when

\[
 \boxed{\Phi_{3,2}
 =k_3^3-4k_4k_3k_2+8k_4^2k_1=0.}                   \tag{2.2}
\]

Thus `C_(3,2)` is a cubic hypersurface in the projective quartic-Hessian
space.

### The `2 o 3` locus

After completing the square in the outer quadratic, a `2 o 3` composition is
affine-equivalent on the target to `C(W)^2`, with `C` cubic.  Direct
differentiation and elimination of the four coefficients of `C` give

\[
\boxed{\begin{aligned}
 \Phi_{2,3}={}&3072k_0k_4^3-768k_1k_3k_4^2
 -320k_2^2k_4^2\\
 &+432k_2k_3^2k_4-81k_3^4=0.
\end{aligned}}                                      \tag{2.3}
\]

There are no hidden equations on `k_4!=0`.  After setting `k_4=1` and writing
`C=W^3+pW^2+qW+r`, the first three nonleading Hessian coefficients recover
`p,q,r` successively; the last coefficient is precisely (2.3).  Hence
`C_(2,3)` is a quartic hypersurface.

Both loci have dimension two in the three-dimensional marked-Hessian moduli
space, proving equality in (1.4) for `(a,b)=(2,3),(3,2)`.

For a normalized sextic

\[
 H(W)=\sum_{j=3}^6h_j(W^j-W^2),
 \qquad h_3+2h_4+3h_5+4h_6=-1,                     \tag{2.4}
\]

equations (2.2)--(2.3) become, after deleting nonzero scalar factors,

\[
 27h_3h_6^2-18h_4h_5h_6+5h_5^3=0                  \tag{2.5}
\]

and

\[
\begin{aligned}
0={}&32h_3h_5h_6^2+64h_3h_6^3+16h_4^2h_6^2
-24h_4h_5^2h_6\\
&+64h_4h_6^3+5h_5^4+64h_5h_6^3+64h_6^4.
\end{aligned}                                      \tag{2.6}
\]

These are the promised low-cost equations in the normalized seed chart.

## 3. The Ritt intersection in degree six

Translate to the center supplied by (2.2).  Then `k_3=k_1=0`, and (2.3)
reduces to

\[
 48k_0k_4-5k_2^2=0.                                 \tag{3.1}
\]

Consequently every point of the exact-degree intersection has, up to a
nonzero scalar,

\[
 K(x)=30x^4+24qx^2+2q^2
     =\bigl((x^3+qx)^2\bigr)''.                     \tag{3.2}
\]

The two integrations can be chosen even about the same center, so the two
vertical parameters are equal.  The common specialization is

\[
\begin{aligned}
 (x^3+qx)^2
 &= (Z^2)\circ(x^3+qx)\\
 &= \bigl(Z(Z+q)^2\bigr)\circ x^2.                 \tag{3.3}
\end{aligned}

This is the full Ritt-intersection diagram.  If `q=0`, it is the monomial
collision `x^6`.  If `q!=0`, a source dilation sends `q` to `-3/4`, after
which `x^3+qx` is a scalar multiple of `T_3`; (3.3) is the Chebyshev collision
`T_2 o T_3 = T_3 o T_2`, up to affine changes.  Therefore the unmarked
degree-six intersection contains no third, asymmetric Ritt type.

The marked intersection is nevertheless one-dimensional.  Translation of
the common center relative to the marked point cannot be removed by the
`G_m` action in `Q_4`.

## 4. Normalized intersection curve and its clean open

Put

\[
 g_{c,q}(W)=\bigl((W-c)^3+q(W-c)\bigr)^2.           \tag{4.1}
\]

The endpoint tangent condition `H(1)=0` is

\[
\boxed{\begin{aligned}
 E(c,q)={}&15c^4-20c^3+12c^2q+15c^2-8cq-6c\\
          &+q^2+2q+1=0.
\end{aligned}}                                      \tag{4.2}
\]

Let

\[
 D(c,q)=g_{c,q}'(1)-g_{c,q}'(0).
\]

On `E=0`, `D!=0`, the normalized seed and its common vertical parameter are

\[
 H_{c,q}(W)=-\frac{g_{c,q}(W)-g_{c,q}(0)-g_{c,q}'(0)W}{D(c,q)},
 \qquad
 s_0=\frac{g_{c,q}'(0)}{D(c,q)}.                    \tag{4.3}
\]

Equation (4.2) is quadratic in `q`, with discriminant

\[
 4c(3c-1)(7c^2-7c+2).                               \tag{4.4}
\]

Thus the normalized marked-seed intersection is a genuine curve rather than
a finite list of seeds.

The intersection with the D1/F2 clean open is explicit.  Write

\[
 Q_{c,q}(W)=\frac{g_{c,q}(W)-g_{c,q}(0)-g_{c,q}'(0)W}{W^2}.
\]

Besides `E=0`, delete the zeros of

\[
 q,\quad D,\quad Q_{c,q}(0),\quad
 \operatorname{disc}_W(Q_{c,q}),\quad
 g_{c,q}''(1)-2D.                                   \tag{4.5}
\]

Here

\[
 Q_{c,q}(0)=15c^4+12c^2q+q^2,                      \tag{4.6}
\]

\[
 \operatorname{disc}_W(Q_{c,q})=
 16c^2(c^2+q)^2(3c^2+q)
 (225c^4+285c^2q+64q^2),                            \tag{4.7}
\]

and `disc(g''_(c,q))=108380160 q^6`.  Conditions (4.5) say respectively:
the Hessian is squarefree, the endpoint normalization exists, zero is an
exact double root, the remaining root sheets are simple, and
`H''(1)!=-2`.  Their complement is the explicit boundary incidence to which
the H3 coarse marked-specialization theorem applies.  None of these factors
vanishes identically on (4.2), so the generic degree-six Ritt intersection is
inside the marked clean locus.

A rational clean witness is obtained from

\[
 c=\frac25,\qquad q=-\frac15.
\]

It gives

\[
 H(W)=-\frac{25}{4}W^6+15W^5-\frac{25}{2}W^4
      +4W^3-\frac14W^2,
 \qquad s_0=\frac7{125},                            \tag{4.8}
\]

with `H''(1)=-14`.  If

\[
 C(W)=(W-2/5)^3-(W-2/5)/5,
\]

then

\[
 H(W)-\frac7{125}W
 =-\frac{25}{4}\left(C(W)^2-\frac4{15625}\right),  \tag{4.9}
\]

and the same right-hand side is a cubic in `(W-2/5)^2`.  This witnesses the
nonempty clean intersection without numerical elimination.

## 5. Consequences for OP-RITT

The degree-six case changes the all-degree workflow:

1. Parametrize `(A o B)''` modulo the middle affine group, target affine
   postcomposition, and marked source dilation.
2. Implicitize only this Hessian map and compute its generic degree.
3. Pull the resulting incidence locus back through the finite D1 map.
4. Reconstruct the finite choices of `s,A,B` only over that pullback.
5. Restrict the explicit root/Hessian collision equations to the incidence
   locus and use H3 for their coarse marked limits.

For `N=6`, steps 1--4 are complete in (2.2)--(4.3), and the clean-boundary
intersection is (4.5).  Sections 6--9 complete the generic-finiteness problem,
degree eight, and the full degree-twelve Ritt-intersection diagram.

The [Gaussian invariant-transport theorem](GAUSSIAN_EXCEPTIONAL_MOMENT_GEOMETRY.md)
applies the triangular optimal-moment inverse to (2.5)--(2.6).  It gives two
explicit sextic hypersurfaces among normalized four-real-Gaussian `GMC`
witnesses.  Notably, the transported `2 o 3` hypersurface is exactly the
degree-six all-double nonsurjective component; the transported `3 o 2`
hypersurface is distinct, and their intersection contains the clean rational
moment point `(-1/4,4,-99/8)`.

The companion [degree-six Ritt atlas](DEGREE_SIX_RITT_ATLAS.md) pulls these
equations back to explicit normalized decomposition charts.  It identifies
the `2^3` and `3^2` omitted-value incidences, their four doubled type-`(6)`
collision points, and the factored coefficient equations where the
exact-double boundary-clean distinction of the affine root sheet degenerates.

## 6. Canonical Hessian atlas and exact codimension in every degree

There is a triangular replacement for elimination in every factor pair.
Work on the leading-coefficient chart of the projective Hessian and normalize
`K=H''` to be monic.  Twice integrate and rescale to the monic primitive

\[
 F_K(W)=N(N-1)\int_0^W\int_0^u K(v)\,dv\,du
       =W^N+\sum_{j=2}^{N-1}c_jW^j.                \tag{6.1}
\]

The missing constant and linear coefficients are exactly the irrelevant
target translation and vertical parameter.  Every exact `(a,b)` composition
has a unique middle-affine normalization

\[
 B(W)=W^b+\sum_{j=1}^{b-1}\beta_jW^j,
 \qquad
 A(Z)=Z^a+\sum_{m=1}^{a-1}\alpha_mZ^m.             \tag{6.2}
\]

Here `B(0)=A(0)=0`; the two conditions use the translation of the middle
coordinate and the target constant.

The coefficients of degrees

\[
 N-1,N-2,\ldots,N-b+1                             \tag{6.3}
\]

recover `beta_(b-1),...,beta_1` successively.  Once these are known, the
coefficients of degrees

\[
 (a-1)b,(a-2)b,\ldots,b                            \tag{6.4}
\]

recover `alpha_(a-1),...,alpha_1`.  Each step is linear in the new unknown.
All remaining coefficients of degrees `2,...,N-1` give exactly

\[
 (N-2)-(a+b-2)=N-a-b                              \tag{6.5}
\]

residual equations.

> **Canonical-atlas theorem.**  On the leading-coefficient chart, the
> Hessian-composition incidence `C_(a,b)` is rational and prime, the map
> (6.2) is birational onto its image, and
> \[
>  \dim\mathcal C_{a,b}=a+b-3,
>  \qquad
>  \operatorname{codim}_{\mathcal Q_{N-2}}\mathcal C_{a,b}=N-a-b. \tag{6.6}
> \]

The reconstruction (6.3)--(6.4) is the rational inverse, so it proves generic
finiteness and algebraic independence of the parameters without a Jacobian
guess.  The residuals (6.5) generate the incidence ideal on this chart.

## 7. Degree eight: the excess triple-quadratic intersection

For compact equations use the monic primitive

\[
 F=W^8+c_7W^7+c_6W^6+c_5W^5+c_4W^4+c_3W^3+c_2W^2. \tag{7.1}
\]

The `2 o 4` incidence is the codimension-two complete intersection

\[
\begin{aligned}
0={}&-128c_3+64c_4c_7+64c_5c_6-48c_5c_7^2
 -48c_6^2c_7+40c_6c_7^3-7c_7^5,\\
0={}&-512c_2+256c_4c_6-64c_4c_7^2+128c_5^2
 -256c_5c_6c_7+64c_5c_7^3\\
&\hspace{2.7em}-64c_6^3+144c_6^2c_7^2-60c_6c_7^4+7c_7^6.
\end{aligned}                                      \tag{7.2}
\]

The `4 o 2` incidence is

\[
\begin{aligned}
0={}&-32c_5+24c_6c_7-7c_7^3,\\
0={}&256c_3-128c_4c_7+20c_6c_7^3-7c_7^5.
\end{aligned}                                      \tag{7.3}
\]

Equivalently, if `K=F''` and `d=-k_5/(6k_6)`, equations (7.3) say

\[
 K'(d)=K'''(d)=0,                                   \tag{7.4}
\]

so the sextic Hessian is even about `d`.

Write the canonical `2 o 4` chart as

\[
 B=W^4+u_3W^3+u_2W^2+u_1W,
 \qquad A=Z^2+vZ.                                  \tag{7.5}
\]

Pulling (7.3) back through (7.5), both equations generate the single prime

\[
 J_8=8u_1-4u_2u_3+u_3^3.                           \tag{7.6}
\]

Indeed `C=B+v/2`, after translation by `-u_3/4`, has odd linear coefficient
`J_8/8`; its cubic coefficient already vanishes.  Thus `J_8=0` says that `C`
is a quadratic in a centered square.  Consequently

\[
 \boxed{\mathcal C_{2,4}\cap\mathcal C_{4,2}
        =\mathcal C_{2,2,2}.}                       \tag{7.7}
\]

Each primary locus has marked dimension three and codimension two in `Q_6`.
Their intersection has marked dimension two, one more than the transverse
expectation.  It is irreducible and has no primitive `2`-versus-`4` Ritt
component.

A rational clean normalized witness is obtained from

\[
 Q=W^2-3W,qquad R=Q^2-3Q,qquad f=R^2-100R.
\]

It gives

\[
 H=-\frac1{340}W^2(W-1)
 (W^5-11W^4+37W^3-17W^2-189W+519).                \tag{7.8}
\]

## 8. Degree twelve: all factor orders and pair intersections

There are four ordered factor loci.  The canonical atlas gives:

| locus | residual equations | marked dimension | codimension in `Q_10` |
|---|---:|---:|---:|
| `C_(2,6)` | 4 | 5 | 4 |
| `C_(3,4)` | 5 | 4 | 5 |
| `C_(4,3)` | 5 | 4 | 5 |
| `C_(6,2)` | 4 | 5 | 4 |

For `C_(6,2)`, the invariant form of the four equations is especially short.
If `K` is the degree-ten Hessian and

\[
 d=-\frac{k_9}{10k_{10}},
\]

then

\[
 K'(d)=K'''(d)=K^{(5)}(d)=K^{(7)}(d)=0.            \tag{8.1}
\]

For `C_(2,6)`, write the specialization as a square of a sextic; the four
residual equations are the degrees `5,4,3,2` left after the six coefficients
of the sextic have been reconstructed.  The two middle factor orders are
likewise the five residuals of (6.5).  This recursive form is substantially
smaller than their expanded Hessian equations.

Exact minimal-prime computation gives the complete pairwise diagram:

| intersection | marked dimension | components | generic description |
|---|---:|---:|---|
| `C_(2,6) cap C_(6,2)` | 3 | 1 | common refinement `2 o 3 o 2` |
| `C_(2,6) cap C_(4,3)` | 3 | 1 | common refinement `2 o 2 o 3` |
| `C_(3,4) cap C_(6,2)` | 3 | 1 | common refinement `3 o 2 o 2` |
| `C_(2,6) cap C_(3,4)` | 2 | 1 | degree-six Ritt collision composed on the right by `2` |
| `C_(4,3) cap C_(6,2)` | 2 | 1 | degree-six Ritt collision composed on the left by `2` |
| `C_(3,4) cap C_(4,3)` | 1 | 2 | power and Chebyshev components |

Every pairwise intersection is therefore nontransverse.  The first three are
ordinary associativity refinements.  The next two transport the unique
degree-six Ritt curve.  The last row is the coprime Ritt-second-theorem
collision and is the only reducible pair intersection.

To make the last row explicit, use the canonical `3 o 4` chart

\[
 B=W^4+b_3W^3+b_2W^2+b_1W,
 \qquad A=Z^3+a_2Z^2+a_1Z.                         \tag{8.2}
\]

The power component is the prime locus

\[
 a_2^2=3a_1,qquad 3b_3^2=8b_2,qquad
 b_2^2-3b_1b_3+4a_2=0.                             \tag{8.3}
\]

Its normal form is

\[
 f_{c,q}(W)=\bigl((W-c)^4+q(W-c)\bigr)^3.          \tag{8.4}
\]

Indeed (8.4) is a cube of a quartic and also a quartic in `(W-c)^3`.

The other prime is the Chebyshev component

\[
 f_{u,v}(W)\sim T_{12}(uW+v),qquad u\ne0,         \tag{8.5}
\]

with the canonical factors obtained from `T_3 o T_4=T_4 o T_3`.  Both
components have affine monic-Hessian dimension two and marked dimension one.
They meet along `q=0`, the source-affine orbit of the monomial `W^12`, which
is one point after the marked source dilation quotient.

## 9. Higher intersections and the clean marked pullback

The two transported degree-six collisions already carry a third factor
order:

\[
\begin{aligned}
 \mathcal C_{2,6}\cap\mathcal C_{3,4}
 &\subset\mathcal C_{6,2},\\
 \mathcal C_{4,3}\cap\mathcal C_{6,2}
 &\subset\mathcal C_{2,6}.                         \tag{9.1}
\end{aligned}
\]

More strikingly,

\[
 \boxed{\mathcal C_{2,6}\cap\mathcal C_{3,4}
 \cap\mathcal C_{4,3}\cap\mathcal C_{6,2}
 =\mathcal T_{12},}                                \tag{9.2}
\]

where `T_12` is the Chebyshev component (8.5).  Any triple containing both
`C_(3,4)` and `C_(4,3)` is already this all-four intersection.  The power
component meets either even factor-order locus only at its monomial point.

All components above meet the D1/F2 clean open.  The exact checker supplies
rational clean specializations for the three common refinements, the two
transported degree-six components, and the power component.  For the
Chebyshev component, set `v=1` in (8.5).  The endpoint tangent equation for
`u` is the squarefree degree-ten polynomial

\[
\begin{aligned}
P_{12}(u)={}&256u^{10}+3072u^9+16128u^8+48640u^7
+93024u^6+117504u^5\\
&+99008u^4+54912u^3+19305u^2+4004u+429.
\end{aligned}                                      \tag{9.3}
\]

It is coprime to the endpoint derivative, the marked Hessian value, both
Hessian and primitive-root discriminants, and the `H''(1)+2` condition.
Thus each of its ten roots gives a clean algebraic normalized seed.

By (1.3), every marked dimension in Sections 7--9 transfers unchanged to the
vertical decomposition loci in normalized seed space.  F2 reconstructs the
seed after adjoining the affine sheet, and H3 gives the unique coarse marked
specialization at the deleted boundary points.

The [degree-six boundary atlas](DEGREE_SIX_RITT_BOUNDARY_ATLAS.md) carries out
that deleted-boundary calculation completely in the first composite degree.
The two Ritt surfaces have respectively two and three irreducible
affine-boundary curves.  On their common marked curve the zero-cluster and
extra-root cuts are reduced of lengths six and six, while the four Hessian
collisions are exactly the type-`(6)` omitted-component support and do not lie
on the affine boundary.
