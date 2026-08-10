# Projective polar geometry of \(\mathrm{HC}_4\)

## Status and outcome

This is a first exact intersection-theoretic experiment, not a proof of
\(\mathrm{HC}_4\).  It establishes the degree formula with fixed
conventions, separates two projective maps that must not be conflated, and
gives a finite numerical atlas for quadratic, cubic, and quartic affine
gradients of geometric degree two or three.  It also gives an exact
leading-Hessian coverage matrix for the first open case, namely quintic
potentials.

The main correction is:

> The full polar map of the homogenized potential is not, in general, the
> projective compactification of the affine gradient map.

The compactification whose top projective degree equals
\(\deg(\nabla\Psi)\) is

\[
 \Gamma_\Psi=
 [X_0^m:G_1:\cdots:G_4]\colon
 \mathbf P^4\dashrightarrow\mathbf P^4,
 \qquad m=\deg(\nabla\Psi)=d-1,                     \tag{0.1}
\]

where \(G_i=X_0^m\Psi_i(X/X_0)\).  Its base scheme is automatically
contained in \(H_\infty=(X_0=0)\).  The notation
\(\deg_{\rm aff}(\nabla\Psi)\) below means the generic fiber cardinality,
not the polynomial degree \(m\).

For this map the exact formula is

\[
\boxed{
\deg_{\rm aff}(\nabla\Psi)
=m^4-6m^2\sigma_2-4m\sigma_3-\sigma_4,}             \tag{0.2}
\]

where

\[
i_*s(B_\infty,\mathbf P^4)
=\sigma_1H+\sigma_2H^2+\sigma_3H^3+\sigma_4H^4
\]

and \(\sigma_1=0\).  The later \(\sigma_i\) are signed Segre degrees;
only the first nonzero component is necessarily effective.

The low-degree numerical sieve gives:

| gradient degree \(m\) | affine degree \(2\) | affine degree \(3\) |
|---:|---:|---:|
| \(2\) | 9 signatures | 7 signatures |
| \(3\) | 72 signatures | 67 signatures |
| \(4\) | 319 signatures | 307 signatures |

If the infinity base scheme is zero-dimensional, its length is forced to
be \(14,13,79,78,254,253\), respectively.  These are necessary
configurations, not existence results.

There is an immediate exact exclusion on the first line.  Constant
nonzero Hessian determinant makes \(\nabla\Psi\) a Keller map, and Wang's
degree-two theorem says that every characteristic-zero quadratic Keller
map is a polynomial automorphism.  Therefore every quadratic-gradient
row of affine degree two or three is empty.  In particular, the proposed
length-\(13\) base scheme cannot occur.

The second line is also empty by the repository's theorem `HC4CQ1`.
After translating the midpoint of a hypothetical collision, every
degree-four potential has the form \(q_2+h_3+h_4\) with a nonzero
antipodal collision, exactly the case excluded there.  Ax--Grothendieck
then upgrades injectivity to polynomial invertibility.  Consequently all
155 listed affine-degree-two/three signatures are unrealizable, and any
four-variable constant-Hessian counterexample must have potential degree
at least five.

At gradient degree four, positivity and log-concavity alone leave 626
rows, so direct Hilbert-scheme realization is not a useful next
calculation.  If \(h_5\) is the leading quintic, its Hessian has rank at
most three.  Rank zero reduces to `HC4CQ1`.  In rank three, the aligned
branch is closed by `HC4CD5`, while the nonaligned branch reduces exactly
to a ternary quintic--cubic Schur divisibility problem.  A squarefree
Hessian determinant forces that cubic to vanish, so only the
nonsquarefree Hessian-discriminant locus remains.  Ranks one and two stop
at constant-kernel synchronization problems.  The four projective degrees
do not determine which Hessian-rank branch occurs, so this structural
matrix does not by itself exclude any of the 626 numerical rows.

The universal top-gradient/Rees sieve now refines that statement.  Rank
does not follow from the degree vector, but the rank and singularity type
of \(h_5\) determine the support codimension: essential rank one feeds only
the codimension-two columns; a squarefree binary rank-two top feeds only
codimension three; and a smooth ternary rank-three top feeds only the
single codimension-four row for each affine degree.  The generic smooth
top ideals are complete intersections of linear type, so further nonlinear
Rees restrictions must come from singular top strata or lower layers.
The first lower-layer calculation now closes the smooth rank-three packet:
the active vertex algebra has length \(256\), while a nonaligned missing
gradient component generates an ideal of length at least six.  Thus affine
degrees two and three cannot occur in codimension four.  The two isolated
rows are empty, leaving 318 and 306 numerical signatures.

The codimension-three calculation then closes a dense open part of the
smooth rank-two packet.  If \(h_4|_K\ne0\), the quartic determinant face
synchronizes a constant kernel direction.  A squarefree binary-quintic
Hessian determinant reaches `HC4CD5`; on the nonsquarefree remainder the
generic transverse length forces \(\sigma_3=16\).  In the rank-three
isolated-singularity packet, the Schur cubic must vanish at every ordinary
singular point.  These are conditional packet restrictions, not additional
unconditional atlas-row deletions.

The exact atlas is generated by
[`scripts/verify_hc4_projective_polar_atlas.py`](scripts/verify_hc4_projective_polar_atlas.py)
and stored in
[`artifacts/generated-results/hc4_projective_polar_atlas.json`](artifacts/generated-results/hc4_projective_polar_atlas.json).
Its degree/Segre transform now comes from the dimension-free implementation
in
[`jcsearch/projective_gradient_segre.py`](jcsearch/projective_gradient_segre.py).
The repository-wide conventions, integrability sieve, cotangent and stable
calibrations, and family attachment ledger are in
[`PROJECTIVE_GRADIENT_SEGRE_MACHINERY.md`](PROJECTIVE_GRADIENT_SEGRE_MACHINERY.md).

## 1. Two projective maps

Let \(\Psi\in k[x_1,\ldots,x_4]\) have degree \(d\), and put \(m=d-1\).
Its degree-\(d\) homogenization is

\[
 \widetilde\Psi(X_0,\ldots,X_4)
 =X_0^d\Psi(X_1/X_0,\ldots,X_4/X_0).
\]

The full polar map is

\[
 \mathcal P_\Psi
 =[\widetilde\Psi_{X_0}:\widetilde\Psi_{X_1}:
                  \cdots:\widetilde\Psi_{X_4}].
\tag{1.1}
\]

On \(X_0=1\), Euler's identity gives

\[
 \mathcal P_\Psi(x)
 =[d\Psi(x)-x\mathbin{\cdot}\nabla\Psi(x):
   \nabla\Psi(x)].                                  \tag{1.2}
\]

This is a projective Gauss map.  It is not the affine gradient
\([1:\nabla\Psi]\), whose closure is (0.1).  Consequently:

1. the top degree of \(\mathcal P_\Psi\) need not equal the generic degree
   of \(\nabla\Psi\);
2. the polar base scheme can have affine points satisfying
   \(\Psi=\nabla\Psi=0\); and
3. absence of affine ramification does not identify the two base schemes.

The distinction already occurs in the constant-Hessian family

\[
 \Psi_r=x_1x_2+x_3x_4+\frac{x_2^{r+1}}{r+1},
 \qquad r\ge2.                                      \tag{1.3}
\]

Its affine gradient is a polynomial automorphism.  Exact Macaulay2
calculation gives

\[
\begin{array}{c|c|c}
r& (g_0,\ldots,g_4)(\Gamma_{\Psi_r})
 &(g_0,\ldots,g_4)(\mathcal P_{\Psi_r})\\ \hline
2&(1,2,2,2,1)&(1,2,4,4,2)\\
3&(1,3,3,3,1)&(1,3,6,6,3).
\end{array}                                        \tag{1.4}
\]

Thus the quadratic full polar map is generically two-to-one although the
affine gradient is an automorphism; in the cubic calibration the two top
degrees are \(3\) and \(1\).  The calculation is certified by
[`scripts/verify_projective_polar_calibrations.m2`](scripts/verify_projective_polar_calibrations.m2).

The full polar map remains useful, but its multidegrees must be recorded as
a second invariant rather than substituted into (0.2).

## 2. Blow-up and Segre formula

Let

\[
 \phi=[F_0:\cdots:F_n]\colon\mathbf P^n\dashrightarrow\mathbf P^n
\]

be defined by forms of common degree \(m\), with base scheme \(B\).  Let
\(\pi:\widehat{\mathbf P^n}\to\mathbf P^n\) be the blow-up of the base
ideal, \(E\) its exceptional divisor, and \(H=\pi^*c_1(O(1))\).  The graph
linear system is \(mH-E\).  Its projective degrees are

\[
 g_i=\int_{\widehat{\mathbf P^n}}
 H^{n-i}(mH-E)^i.                                  \tag{2.1}
\]

Write

\[
 i_*s(B,\mathbf P^n)=\sum_{k=1}^n\sigma_kH^k.
\tag{2.2}
\]

The convention

\[
 i_*s(B,\mathbf P^n)
 =\pi_*\!\left(\frac{E}{1+E}\right)
 =\sum_{k\ge1}(-1)^{k-1}\pi_*(E^k)                 \tag{2.3}
\]

turns the binomial expansion of (2.1) into

\[
\boxed{
g_i=m^i-\sum_{k=1}^i
 \binom{i}{k}m^{i-k}\sigma_k.}                     \tag{2.4}
\]

For \(\Gamma_\Psi\), \(g_4=\deg_{\rm aff}(\nabla\Psi)\), because its
restriction to \(X_0\ne0\) is exactly \([1:\nabla\Psi]\) and a general
target point lies in that chart.  The base ideal contains \(X_0^m\).
It has no fixed divisor: any common divisor would have to divide \(X_0\),
but at least one derivative of the nonzero top homogeneous part of
\(\Psi\) is nonzero in characteristic zero.  Hence \(\sigma_1=0\), and
(2.4) gives (0.2).

For reference, the full \(\mathbf P^4\) degree list is

\[
\begin{aligned}
g_0&=1,\\
g_1&=m,\\
g_2&=m^2-\sigma_2,\\
g_3&=m^3-3m\sigma_2-\sigma_3,\\
g_4&=m^4-6m^2\sigma_2-4m\sigma_3-\sigma_4.
\end{aligned}                                      \tag{2.5}
\]

## 3. What codimension and Hilbert data determine

The leading nonzero Segre component is the top-dimensional normal-cone
cycle.  Therefore:

- \(\sigma_2>0\) detects a codimension-two component and records its
  cycle degree with generic multiplicities;
- if \(\sigma_2=0\) and \(\sigma_3>0\), the leading base has codimension
  three and degree \(\sigma_3\);
- if \(\sigma_2=\sigma_3=0\), the base is zero-dimensional and
  \(\sigma_4\) is its length.

The ordinary Hilbert polynomial does **not** determine the full Segre
class of an arbitrary singular or nonreduced scheme.  The normal cone, or
equivalently Rees-algebra data, is needed.  Thus codimension,
multiplicity, and Hilbert polynomial are a first filter, not complete
input for (0.2).

There are useful regular-embedding subcases.  If the base is a smooth
curve \(C\subset\mathbf P^4\) of degree \(e\) and genus \(g\), with no
other components, then

\[
 \sigma_2=0,\qquad
 \sigma_3=e,\qquad
 \sigma_4=2-2g-5e.                                 \tag{3.1}
\]

This follows from
\(s(C,\mathbf P^4)=c(N_{C/\mathbf P^4})^{-1}[C]\).
Its polar degree is therefore

\[
 g_4=m^4-(4m-5)e-2+2g.                             \tag{3.2}
\]

In this smooth-curve-only class, the numerical atlas and Castelnuovo's
bound leave:

\[
\begin{array}{c|c|c}
m&g_4=2&(e,g)\\ \hline
2&& (4,0)\\
3&& (11,0),(13,7),(15,14),(17,21),(19,28),(21,35)
\end{array}
\]

and

\[
\begin{array}{c|c|c}
m&g_4=3&(e,g)\\ \hline
2&&\text{none}\\
3&& (12,4),(14,11),(16,18),(18,25),(20,32).
\end{array}                                        \tag{3.3}
\]

These pairs pass only the displayed numerical tests.  They are not
asserted to occur as infinity base schemes of gradient Keller maps.

## 4. Quadratic-gradient atlas and exclusion

Projective degree sequences are positive, bounded by \(g_i\le m^i\), and
log-concave:

\[
 g_i^2\ge g_{i-1}g_{i+1}.                           \tag{4.1}
\]

For \(m=2\), these elementary constraints make the full list short.
The entries below are \((g_2,g_3;\sigma_2,\sigma_3,\sigma_4)\).

\[
\begin{array}{c|l}
g_4&\text{all numerical signatures}\\ \hline
2&
(2,2;2,-6,14),\
(3,3;1,-1,-2),\
(3,4;1,-2,6),\\
& (4,3;0,5,-26),\
(4,4;0,4,-18),\
(4,5;0,3,-10),\\
& (4,6;0,2,-2),\
(4,7;0,1,6),\
(4,8;0,0,14);\\[2mm]
3&
(3,3;1,-1,-3),\
(3,4;1,-2,5),\\
& (4,4;0,4,-19),\
(4,5;0,3,-11),\
(4,6;0,2,-3),\\
& (4,7;0,1,5),\
(4,8;0,0,13).
\end{array}                                        \tag{4.2}
\]

This is the first finite target list before using the Keller condition.  A
quadratic-gradient map of affine degree two or three would have to realize
one of these Segre signatures.  The next paragraph applies the
constant-Hessian hypothesis and removes all of them.

In fact the constant-Hessian condition removes every row at once.  The
map

\[
 F=\nabla\Psi:\mathbb A^4\longrightarrow\mathbb A^4
\]

has polynomial degree two and

\[
 \det DF=\det\operatorname{Hess}\Psi\in k^\times.
\]

Wang's quadratic Keller theorem therefore makes \(F\) a polynomial
automorphism.  Hence

\[
 g_4=\deg_{\rm aff}F=1                              \tag{4.3}
\]

and (0.2) forces

\[
\boxed{24\sigma_2+8\sigma_3+\sigma_4=15.}           \tag{4.4}
\]

The nine affine-degree-two rows and seven affine-degree-three rows in
(4.2) are consequently unrealizable by any characteristic-zero
constant-Hessian potential.  This includes both isolated candidates:

\[
\boxed{\text{no length-14 or length-13 zero-dimensional infinity base}.}
\tag{4.5}
\]

This is not a new proof of Wang's theorem.  It is the exact
intersection-theoretic consequence of that theorem for the infinity base
scheme.  The numerical degree-one atlas still has eleven log-concave
Segre signatures; determining which occur among quadratic gradient
automorphisms is an automorphism-geometry question, not an HC4
counterexample question.  The external input is S. S. Wang,
[*A Jacobian criterion for
separability*](https://doi.org/10.1016/0021-8693(80)90233-1),
Journal of Algebra **65** (1980), especially the quadratic-relations
case.

For \(m=3\), the corresponding elementary atlas has 72 and 67 rows.
They are retained in the generated JSON rather than printed here.  The
zero-dimensional rows have lengths \(79\) and \(78\).

These rows are likewise empty.  Suppose \(\deg\Psi\le4\) and
\(\nabla\Psi(a)=\nabla\Psi(b)\) with \(a\ne b\).  Translate the midpoint
\((a+b)/2\) to the origin, subtract the common affine-linear gradient
term, and decompose

\[
 \Psi=q_2+h_3+h_4.
\]

The quadratic part is nondegenerate because its Hessian is the Hessian of
\(\Psi\) at the translated origin and the determinant is a nonzero
constant.  The collision becomes a nonzero antipodal critical pair.
Theorem `HC4CQ1` in
[`HC4_MENG_DENSE_CUBIC_QUARTIC.md`](HC4_MENG_DENSE_CUBIC_QUARTIC.md)
excludes precisely this configuration.  Thus \(\nabla\Psi\) is injective
over the algebraic closure.  Ax--Grothendieck makes it a polynomial
automorphism, so

\[
 g_4=1,\qquad
\boxed{54\sigma_2+12\sigma_3+\sigma_4=80.}          \tag{4.6}
\]

Therefore all 72 affine-degree-two and 67 affine-degree-three
cubic-gradient signatures are unrealizable, including the isolated
lengths \(79\) and \(78\).  Combining (4.5) and (4.6) gives the concrete
degree bound

\[
\boxed{\text{an \(\mathrm{HC}_4\) counterexample must satisfy }
       \deg\Psi\ge5\quad(\deg\nabla\Psi\ge4).}       \tag{4.7}
\]

This lower bound is a projective restatement of the existing collision
theorem, not a new proof of `HC4CQ1`.

## 5. Cotangent lifts of quartic plane packets

On the isotropic cotangent branch,

\[
 \Psi_F(x,y,t,m)=tP(x,y)+mQ(x,y)+H(x,y),
\qquad F=(P,Q),
\tag{5.1}
\]

and

\[
 \nabla\Psi_F=
 ((DF)^T(t,m)^T+\nabla H,\ P,\ Q).                 \tag{5.2}
\]

If \(F\) is Keller, then over every point of a generic \(F\)-fiber the
first two equations determine \((t,m)\) uniquely.  Hence

\[
 \boxed{\deg_{\rm aff}(\nabla\Psi_F)
       =\deg_{\rm aff}(F).}                         \tag{5.3}
\]

A surviving plane quartic packet means geometric degree four.  Therefore
its cotangent lift has

\[
 6m^2\sigma_2+4m\sigma_3+\sigma_4=m^4-4.           \tag{5.4}
\]

In particular, the aggregate corrections would be

\[
\begin{array}{c|c}
m&\text{total correction}\\ \hline
2&12\\
3&77.
\end{array}                                        \tag{5.5}
\]

This does not produce an individual Segre signature.  The existing
quartic packet records finite-normalization and boundary data, not a
specific pair \(P,Q\) and common homogenization degree.  Different
polynomial representatives of the same geometric-degree packet can have
different infinity base schemes.  Moreover the repository's cited plane
degree bounds exclude actual noninvertible Keller maps in the low
coordinate-degree strata where the reductions use them.  Thus (5.5) is a
conditional target for a hypothetical packet, not a constructed row.

## 6. The \(R(P)=0\) bordered branch

For

\[
 \Psi=tP(x,y)+\Phi(x,y,m),
\]

the Schur coefficient gives

\[
 [t]\det\operatorname{Hess}\Psi=-\Phi_{mm}R(P),
\qquad
 R(P)=(\nabla P)^T\operatorname{adj}(\operatorname{Hess}P)\nabla P.
\tag{6.1}
\]

The equation \(R(P)=0\) alone does not determine an infinity base scheme
or a Segre class.  The constant-Hessian family

\[
 \Theta_r=tx+my+\frac{x^{r+1}}{r+1}                \tag{6.2}
\]

has \(P=x\), hence \(R(P)=0\), and

\[
 \det\operatorname{Hess}\Theta_r=1,\qquad
 \nabla\Theta_r=(t+x^r,m,x,y),
\]

which is a polynomial automorphism for every \(r\).  Nevertheless its
degree-\(r\) graph compactification has nontrivial and varying infinity
base scheme.  For \(r=2,3\), (1.4) gives Segre signatures

\[
 (0,2,-6,15),\qquad (0,6,-30,116),                 \tag{6.3}
\]

and aggregate corrections \(15\) and \(80\).

Thus “compute the contribution of the \(R(P)=0\) branch” is
underdetermined without the full \(\Phi\) and its lower homogeneous
layers.  What is exact in the current low-degree HC4 reductions is instead
a structural exclusion: the bordered equation forces a further constant
direction on the covered quadratic/cubic-leading strata, after which the
problem enters the cotangent or terminal \(\mathrm{HC}_2\) block.  The
relevant exact calculations remain in
[`scripts/verify_hc4_meng_dense_cubic_quartic_reduction.py`](scripts/verify_hc4_meng_dense_cubic_quartic_reduction.py)
and
[`scripts/verify_hc4_meng_triple_rank_three_reduction.py`](scripts/verify_hc4_meng_triple_rank_three_reduction.py).

## 7. Meng--Yang control before and after Schur descent

The six-variable doubling potential has degree \(8\), so its gradient
degree is \(m=7\).  Its affine gradient is generically three-to-one.  If

\[
 i_*s(B_6,\mathbf P^6)=\sum_{k=1}^6\sigma_kH^k,
\]

then the exact aggregate correction is

\[
 \sum_{k=1}^6\binom{6}{k}7^{6-k}\sigma_k
 =7^6-3
 =117646.                                          \tag{7.1}
\]

The Schur-descended five-variable potential has degree \(14\), hence
\(m=13\).  Its generic degree is also three.  Indeed, write

\[
 F(x,t)=F_0(x)+tF_1(x),\quad
 A=y^TF_1(x),
\]

and

\[
 \psi(x,y)=y^TF_0(x)+\frac{\lambda}{2}A^2+\mu A.
\]

Put \(t^*=\mu+\lambda A\).  Then

\[
 \nabla_y\psi=F(x,t^*),\qquad
 \nabla_x\psi=(D_xF(x,t^*))^Ty.                    \tag{7.2}
\]

For each point \((x,t^*)\) in a generic \(F\)-fiber, the remaining
equations together with
\(F_1(x)^Ty=(t^*-\mu)/\lambda\) form

\[
 DF(x,t^*)^Ty=
 \left(p,\frac{t^*-\mu}{\lambda}\right),
\]

which has a unique solution because \(F\) is Keller.  Schur descent
therefore preserves the generic fiber cardinality.  Consequently

\[
 \sum_{k=1}^5\binom{5}{k}13^{5-k}\sigma'_k
 =13^5-3
 =371290.                                          \tag{7.3}
\]

Equations (7.1)--(7.3) are the requested control totals.  Computing every
individual \(\sigma_k\) of these high-degree base ideals is a separate
Rees/Segre calculation and is not claimed here.  The source construction
and determinant identities are canonically recorded in
[`MENG_YANG_SCHUR_DESCENT_BRIDGE.md`](MENG_YANG_SCHUR_DESCENT_BRIDGE.md);
the external source is Meng--Yang,
[*A five-variable counterexample to the Hessian conjecture, and the
low-dimensional status of the Jacobian and Hessian
conjectures*](https://arxiv.org/abs/2607.22198).

## 8. Quintic coverage and gap matrix

Let a hypothetical degree-five collision be normalized at its midpoint:

\[
 \psi=q_2+h_3+h_4+h_5,
 \qquad
 H(\lambda)=H_0+\lambda A+\lambda^2B+\lambda^3D,     \tag{8.1}
\]

where \(A,B,D\) are the Hessians of \(h_3,h_4,h_5\).  The constant
Hessian identity says

\[
 \det H(\lambda)=\det H_0.                            \tag{8.2}
\]

Its top face is \(\det D=0\).  By the low-dimensional
Gordan--Noether theorem and iteration after removing inessential
variables, \(h_5\) has a constant kernel \(K\), with
\(\dim K=4-r\) when \(r=\operatorname{rank}D\).  In coordinates adapted
to \(K\), let \(C_r\) be the nonsingular Hessian block of \(h_5\).
The highest remaining determinant face is

\[
\begin{array}{c|c|c}
r&\lambda\text{-degree}&\text{face}\\ \hline
3&11&\det(C_3)\det(B|_K)\\
2&10&\det(C_2)\det(B|_K)\\
1&9 &\det(C_1)\det(B|_K).
\end{array}                                           \tag{8.3}
\]

For rank three write \(K=\langle\partial_t\rangle\) and let
\(u=(x,y,m)\).  Equation (8.3) first gives

\[
 D_t^2h_4=0,\qquad
 h_4=t\,s_3(u)+r_4(u).                                \tag{8.4}
\]

The next face, of degree ten, is exactly

\[
 \det(C_3)D_t^2h_3
 -
 (\nabla s_3)^{\mathsf T}
 \operatorname{adj}(C_3)\nabla s_3=0.                 \tag{8.5}
\]

If \(s_3=0\), equation (8.5) forces \(D_t^2h_3=0\).
Together with \(D_th_5=0\), these are precisely the hypotheses of
`HC4CD5` with \(h_6=0\), so this aligned branch is closed.  If
\(s_3\ne0\), every survivor must solve the new polynomial divisibility

\[
 \boxed{
 \det\operatorname{Hess}_u(h_5)
 \mid
 (\nabla s_3)^{\mathsf T}
 \operatorname{adj}(\operatorname{Hess}_u(h_5))
 \nabla s_3.}                                        \tag{8.6}
\]

The determinant has degree nine, the numerator degree ten, and the
quotient is the linear form \(D_t^2h_3\).  This is the quintic analogue
of the repository's sextic--quintic Schur face, but it is not covered by
those theorems: the present pair is a ternary quintic and a ternary cubic,
not a ternary sextic and a ternary quartic.

The generic branch nevertheless closes without coefficient elimination.
Put

\[
 \Delta=\det(C_3),\qquad
 w=\operatorname{adj}(C_3)\nabla s_3.                \tag{8.7}
\]

Suppose \(\Delta\) is squarefree.  At the generic point of each
irreducible component \(p\mid\Delta\), the matrix \(C_3\) has corank one:
over the discrete valuation ring at \(p\), squarefreeness says that
\(\det(C_3)\) has valuation one.  Hence
\(\operatorname{adj}(C_3)\) is a nonzero symmetric rank-one matrix modulo
\(p\).  Reducing (8.6) modulo \(p\) gives a square of one linear pairing,
so characteristic zero implies

\[
 w=0\pmod p.
\]

Every component \(p\) therefore divides every entry of \(w\).  Since
\(\Delta\) is squarefree, \(\Delta\mid w\).  But \(C_3\) has cubic
entries, so \(\deg\Delta=9\) and

\[
 \deg w=2\deg(C_3)+\deg(\nabla s_3)=6+2=8.
\]

Thus \(w=0\).  The adjugate identity \(C_3w=\Delta\nabla s_3\) then gives
\(\nabla s_3=0\), hence \(s_3=0\).  The branch reaches `HC4CD5`.
Consequently:

> **Theorem `HC4PPG5` — Squarefree rank-three obstruction.**  A
> nonaligned rank-three
> quintic candidate can exist only if the degree-nine ternary Hessian
> determinant \(\det\operatorname{Hess}_u(h_5)\) is nonsquarefree.

This is a genuine generic exclusion: nonsquarefreeness is a proper
discriminant condition in the coefficient space of ternary quintics.  An
exact witness outside it is

\[
 h_5=x^5+y^5+m^5+x^4y+y^4m+m^4x,                    \tag{8.8}
\]

whose degree-nine Hessian determinant has gcd one with its three partial
derivatives.  The checker verifies this squarefreeness over
\(\mathbb Q\).

Together with the numerical atlas, equations (8.3)--(8.6) give theorem
`HC4PPG4`.  Its exact coverage matrix is:

| rank of \(\operatorname{Hess}(h_5)\) | first structural consequence | status |
|---:|---|---|
| \(0\) | \(h_5=0\) by homogeneity | closed by `HC4CQ1` |
| \(1\) | \(\det(B|_K)=0\) on a constant three-plane | open constant-kernel synchronization |
| \(2\) | \(\det(B|_K)=0\) on a constant two-plane | open constant-kernel synchronization |
| \(3\), \(D_th_4=0\) | (8.5) gives \(D_t^2h_3=0\) | closed by `HC4CD5` |
| \(3\), \(D_th_4\ne0\), squarefree \(\det C_3\) | degree \(9>8\) adjugate obstruction | closed |
| \(3\), \(D_th_4\ne0\), nonsquarefree \(\det C_3\) | cubic Schur divisibility (8.6) | open exceptional locus |

For ranks one and two, (8.3) supplies a singular Hessian restriction over
the quotient function field.  Its kernel need not yet be a constant
direction over the ground field; proving that synchronization, or
classifying its moving exceptions, is the precise missing step.

The 626 numerical signatures split by leading base codimension as follows:

| affine degree | codimension \(2\) | codimension \(3\) | codimension \(4\) |
|---:|---:|---:|---:|
| \(2\) | 260 | 58 | 1 |
| \(3\) | 249 | 57 | 1 |

The codimension-four rows have lengths \(254\) and \(253\).  No entry in
this table is excluded by (8.3)--(8.6) alone, because those determinant
faces do not recover the rank of \(D\), or the alignment of its kernel,
from the four projective degrees.  The support theorem (9.6) supplies the
missing assignment, and the vertex-colength theorem in Section 10 then
excludes both rows.

There is nevertheless a small Hilbert-polynomial filter inside the
codimension-three column.  If the entire base is assumed to be one smooth
integral curve with no additional components, the Segre formula determines
its genus; Castelnuovo's bound leaves only 18 of the 58 affine-degree-two
rows and 18 of the 57 affine-degree-three rows.  This excludes the other
rows only for that smooth-curve-only base type, not for reducible,
nonreduced, or embedded infinity schemes.

The next bounded calculation should classify the repeated components of
degree-nine Hessian determinants of ternary quintics.  Local cone-vertex
Segre contributions should then be computed only on the surviving
nonsquarefree strata.  A universal realization search over all 626 base
schemes should not be attempted first.

## 9. Universal quintic infinity-gradient and Rees sieve

The numerical atlas does not yet use the presentation of the base ideal by
one gradient.  The first presentation-level calculation can nevertheless
be made without constructing the lower layers \(h_2,h_3,h_4\).

Write the universal top quintic as

\[
 h_5=\sum_{|\alpha|=5}c_\alpha x^\alpha,            \tag{9.1}
\]

with \(56\) coefficients.  Its four quartic derivatives satisfy the six
curl identities, Euler's identity, and its differentiated form

\[
 \sum_i x_i\partial_i h_5=5h_5,\qquad
 \operatorname{Hess}(h_5)x=4\nabla h_5.             \tag{9.2}
\]

The constant-Hessian equation first imposes

\[
 \det\operatorname{Hess}(h_5)=0.                    \tag{9.3}
\]

This genuinely removes the generic rank-four stratum: for
\(h_5=\sum_i x_i^5\), the determinant is
\(20^4x_1^3x_2^3x_3^3x_4^3\).  In the remaining low-dimensional
Gordan--Noether strata, choose coordinates so that \(h_5\) depends on
exactly \(r=1,2,3\) essential variables.  The universal essential forms
have respectively \(1,6,21\) coefficients.

On the open locus where the essential hypersurface
\((h_5=0)\subset\mathbf P^{r-1}\) is smooth, its \(r\) active partial
derivatives form a regular sequence.  The infinity gradient ideal is of
linear type.  Its Rees ideal has only:

1. the \(\binom r2\) Koszul equations
   \[
   (\partial_i h_5)U_j-(\partial_jh_5)U_i=0;         \tag{9.4}
   \]
2. the \(4-r\) linear equations for the inactive target coordinates.

Adding the compactifying generator \(X_0^4\) adds the corresponding
Koszul relations and no nonlinear Rees equation.  This is verified
independently with the Macaulay2
[`ReesAlgebra`](https://macaulay2.com/doc/Macaulay2/share/doc/Macaulay2/ReesAlgebra/html/index.html)
package.  Thus the generic
smooth top strata do **not** hide an additional nonlinear Rees constraint:
such constraints must occur on singular top quintics or after the lower
layers complete the base ideal.

The pure-top ideals

\[
 (X_0^4,\partial_1h_5,\ldots,\partial_rh_5)          \tag{9.5}
\]

are equal-degree complete intersections of codimension \(r+1\).  Their
calibration vectors are:

| \(r\) | pure-top \(g\)-vector | pure-top \(\sigma\)-vector |
|---:|---|---|
| \(1\) | \((1,4,0,0,0)\) | \((0,16,-128,768)\) |
| \(2\) | \((1,4,16,0,0)\) | \((0,0,64,-768)\) |
| \(3\) | \((1,4,16,64,0)\) | \((0,0,0,256)\) |

These are degeneration calibrations, not the Segre vectors of a completed
constant-Hessian potential.  Lower layers can change normal-cone
multiplicities.  What they cannot change is the reduced support on
\(X_0=0\):

\[
 \operatorname{Supp}B_\infty
 =V(\partial_1h_5,\ldots,\partial_4h_5).             \tag{9.6}
\]

For an essential rank-\(r\) cone this support is the join of the
constant-kernel vertex \(\mathbf P^{3-r}\) with
\(\operatorname{Sing}(h_5=0)\subset\mathbf P^{r-1}\).  Therefore the
leading base codimension, hence the first possible nonzero Segre degree,
is determined exactly by the singularity type:

| rank and essential singularity | base codimension | atlas rows for affine degree \(2,3\) |
|---|---:|---:|
| \(r=1\) | \(2\) | \(260,249\) |
| \(r=2\), squarefree binary quintic | \(3\) | \(58,57\) |
| \(r=2\), repeated binary root | \(2\) | \(260,249\) |
| \(r=3\), smooth ternary quintic | \(4\) | \(1,1\) |
| \(r=3\), isolated singularities | \(3\) | \(58,57\) |
| \(r=3\), positive-dimensional singular locus | \(2\) | \(260,249\) |

The support and codimension columns are the \(n=4,m=4\) instances of the
all-dimensional singular-stratum theorem `PGS3`.  That theorem also
explains why this table cannot assign a Segre multiplicity: along a
singular component, the lower active gradients determine a finite
\(k[[X_0]]\)-module whose generic rank and \(X_0\)-torsion orders enter the
transverse length.

The corresponding symbolic restrictions, with
\(\delta=g_4=\deg_{\rm aff}\nabla\Psi\), are

\[
\begin{array}{c|c|c}
\operatorname{codim}B_\infty&
\text{Segre vanishing}&\text{top equation}\\ \hline
2&\sigma_2>0&
96\sigma_2+16\sigma_3+\sigma_4=256-\delta\\
3&\sigma_2=0,\ \sigma_3>0&
16\sigma_3+\sigma_4=256-\delta\\
4&\sigma_2=\sigma_3=0&
\sigma_4=256-\delta.
\end{array}                                         \tag{9.7}
\]

In particular, a smooth essential rank-three \(h_5\) leaves only the
single signatures

\[
 (0,0,0,254),\qquad(0,0,0,253)                     \tag{9.8}
\]

for affine degrees two and three.  A squarefree essential rank-two binary
quintic leaves only the codimension-three columns.  Rank one lies only in
the codimension-two columns.  By itself this does not exclude a row
unconditionally, because the singular top strata together still meet all
three columns; it does replace the undifferentiated 626-row search by
explicit rank/singularity packets.  The smooth rank-three packet is
sharpened in the next section.

The midpoint collision supplies a second independent presentation
constraint.  At a nonzero antipodal point \(a\), parity gives

\[
 \nabla(h_2+h_4)(a)=0,\qquad
 \nabla(h_3+h_5)(a)=0.                              \tag{9.9}
\]

Equation (9.9) does not put \([a]\) in (9.6), but it couples the lower
normal-cone corrections to the top quintic before any Rees elimination.
The determinant faces (8.3)--(8.6) then impose the constant-determinant
conditions on those corrections.

> **Theorem `HC4PPG6` -- Top-gradient/Rees support sieve.**  For a
> collision-normalized quintic `HC4` candidate, the universal top gradient
> satisfies (9.2)--(9.4), rank four is excluded by (9.3), and its reduced
> infinity base support is the kernel-vertex/singular-locus join (9.6).
> Consequently the numerical atlas splits exactly by the six rows of the
> preceding table and obeys (9.7).  On each generic smooth essential
> rank stratum, the top ideal is a complete intersection of linear type,
> so no additional nonlinear top Rees equation exists.  The pure-top
> vectors above are calibrations only; no lower-layer Segre multiplicity or
> new `HC4` exclusion is asserted.

The singular-support part is precisely the `PGS3` specialization.  Its DVR
profile identifies the additional calculation needed on each remaining
singular packet.

The exact artifact is generated by
[`scripts/analyze_hc4_quintic_infinity_rees.py`](scripts/analyze_hc4_quintic_infinity_rees.py)
and independently checks the complete-intersection/Rees claims with
[`scripts/verify_hc4_quintic_infinity_rees_strata.m2`](scripts/verify_hc4_quintic_infinity_rees_strata.m2).

## 10. Smooth rank-three vertex colength

The isolated codimension-four packet admits a stronger local calculation.
Let \(t\) span the constant kernel of
\(\operatorname{Hess}(h_5)\), and write \(u=(u_1,u_2,u_3)\).  Thus
\(h_5=h_5(u)\).  Assume that the ternary quintic
\((h_5=0)\subset\mathbf P^2\) is smooth.  Its three quartic partial
derivatives form a regular sequence, and

\[
 B=k[u_1,u_2,u_3]/(\partial_{u_1}h_5,
                    \partial_{u_2}h_5,
                    \partial_{u_3}h_5)              \tag{10.1}
\]

is an Artinian complete intersection with

\[
 \operatorname{Hilb}_B(z)=(1+z+z^2+z^3)^3,\qquad
 \dim_k B=64,\qquad \operatorname{socdeg}B=9.       \tag{10.2}
\]

The first rank-three determinant face gives \(D_t^2h_4=0\), so

\[
 h_4=t\,s_3(u)+r_4(u),\qquad s_3=D_th_4.            \tag{10.3}
\]

Work on the vertex chart \(t=1\) and put \(\epsilon=X_0\).  The three
active homogenized gradient components are

\[
 G_i=\partial_i h_5+\epsilon\partial_i h_4
       +\epsilon^2\partial_i h_3+\epsilon^3\partial_i h_2,
                                                               \tag{10.4}
\]

while the missing component is

\[
 G_t=\epsilon s_3+\epsilon^2D_th_3+\epsilon^3D_th_2.            \tag{10.5}
\]

There is a useful flatness lemma behind the length calculation.  In
\(\widehat R=k[[\epsilon,u_1,u_2,u_3]]\), put

\[
 R'=\widehat R/(G_1,G_2,G_3),\qquad A=R'/(\epsilon^4).           \tag{10.6}
\]

The ideal \((\epsilon,G_1,G_2,G_3)\) is primary to the maximal ideal,
because its special fiber is (10.1).  Hence \(G_1,G_2,G_3\) are a regular
sequence, \(R'\) is one-dimensional Cohen--Macaulay, and the parameter
\(\epsilon\) is a nonzerodivisor.  Completeness and finite special fiber
make \(R'\) finite over \(k[[\epsilon]]\); torsion-freeness over this
discrete valuation ring makes it free of rank \(64\).  Consequently

\[
 \operatorname{gr}_\epsilon A
   \simeq B\otimes_k k[\epsilon]/(\epsilon^4),
 \qquad \dim_kA=4\cdot64=256.                       \tag{10.7}
\]

Now suppose \(s_3\ne0\), the nonaligned branch.  Its class in \(B\) is
nonzero because the Jacobian ideal in (10.1) is generated in degree four.
It cannot be a socle element: \(s_3\) has degree three, while the socle of
\(B\) is concentrated in degree nine.  Therefore

\[
 \dim_k(Bs_3)\ge2.                                  \tag{10.8}
\]

The initial form of \(G_t\) is \(\epsilon s_3\).  The associated graded
ideal of \(AG_t\) contains the ideal generated by this initial form, so

\[
\begin{aligned}
 \dim_k(AG_t)
 &\ge \dim_k\bigl((\epsilon s_3)
       (B\otimes k[\epsilon]/(\epsilon^4))\bigr)\\
 &=3\dim_k(Bs_3)\ge6.                               \tag{10.9}
\end{aligned}
\]

Because the top support is the single vertex, the completed infinity base
has coordinate algebra \(A/(G_t)\).  In codimension four,
\(\sigma_4=\dim_kA/(G_t)\), while the projective-degree formula gives
\(\sigma_4=256-\delta\), where
\(\delta=\deg_{\rm aff}\nabla\Psi\).  The exact sequence for the principal
ideal \(AG_t\) therefore yields the key identity

\[
 \delta=256-\sigma_4=\dim_k(AG_t)\ge6.              \tag{10.10}
\]

If \(s_3=0\), the aligned branch is already excluded by `HC4CD5`.  We have
therefore proved:

> **Theorem `HC4PPG7` -- Smooth rank-three vertex-colength
> obstruction.**  A collision-normalized quintic `HC4` candidate whose
> top quintic has essential Hessian rank three and defines a smooth
> ternary quintic has affine gradient degree at least six.  Hence
> it cannot have affine gradient degree two or three, and
> the two codimension-four atlas signatures
> \[
> (g_0,\ldots,g_4)=(1,4,16,64,2),(1,4,16,64,3)
> \]
> and
> \[
> (\sigma_1,\ldots,\sigma_4)=(0,0,0,254),(0,0,0,253)
> \]
> are empty.

This is the \((n,m,r)=(4,4,3)\) specialization of the all-dimensional
smooth-essential normal-slice theorem `PGS2` in
[`PROJECTIVE_GRADIENT_SEGRE_MACHINERY.md`](PROJECTIVE_GRADIENT_SEGRE_MACHINERY.md).

The dependency-free exact ledger is
[`scripts/verify_hc4_rank3_vertex_colength.py`](scripts/verify_hc4_rank3_vertex_colength.py).
The independent Macaulay2 replay
[`scripts/verify_hc4_rank3_vertex_colength.m2`](scripts/verify_hc4_rank3_vertex_colength.m2)
checks (10.2), the length \(256\), and exact Fermat and deformed filtered
calibrations.  Those representatives calibrate the length mechanism; the
universal result is the flatness/socle argument (10.6)--(10.10).

## 11. Codimension-three gradient strata

The codimension-three column has two geometrically different sources:

1. a smooth essential binary rank-two quintic, whose support is the kernel
   line \(\mathbf P(K)\);
2. an essential ternary rank-three quintic with isolated singularities,
   whose support contains the lines joining those singularities to the
   kernel vertex.

Both admit presentation-level restrictions, but they must be treated
separately.

### 11.1 Smooth rank two

Let \(K=\langle t,w\rangle\) be the constant two-dimensional kernel of
\(\operatorname{Hess}(h_5)\), and let \(U=(u_1,u_2)\) be the active
variables.  The first determinant face is

\[
 \det\operatorname{Hess}_K(h_4)=0.                  \tag{11.1}
\]

Suppose first that the pure-kernel restriction \(h_4|_K\) is nonzero.  A
binary quartic with zero Hessian is a fourth power of a linear form.  One
quick proof dehomogenizes it as \(y^4p(x/y)\):

\[
 \det\operatorname{Hess}(y^4p(x/y))|_{y=1}
 =3\bigl(4pp''-3(p')^2\bigr).                       \tag{11.2}
\]

At a root of multiplicity \(m\), the leading coefficient of the right-hand
side is proportional to \(m(m-4)\); hence every root has multiplicity four.
After a linear change and rescaling, take \(h_4|_K=t^4\).

Equation (11.1) synchronizes more than this restriction.  Over the quotient
field \(k(U)\), write

\[
 h_4=t^4+\sum_{i+j\le3}a_{ij}t^iw^j.                \tag{11.3}
\]

The eight coefficients of the Hessian determinant have radical containing

\[
 a_{02},a_{11},a_{03},a_{12},a_{21}.                \tag{11.4}
\]

The exact certificate is

\[
 a_{02}^2,\ a_{11}^3,\ a_{03},\ a_{12},\ a_{21}^3
 \in I_{\det\operatorname{Hess}_K(h_4)}.             \tag{11.5}
\]

Thus \(h_4=P(t,U)+w\,s_3(U)\), and the constant direction \(w\) satisfies

\[
 D_w^2h_4=D_tD_wh_4=0,\qquad s_3=D_wh_4.            \tag{11.6}
\]

Put \(C=\operatorname{Hess}_U(h_5)\),
\(\Delta=\det C\), and \(d=\nabla_U s_3\).  The next determinant
coefficient, in \(\lambda\)-degree nine, is

\[
 D_t^2h_4\,
 \bigl(\Delta D_w^2h_3-d^{\mathsf T}\operatorname{adj}(C)d\bigr)=0.
                                                               \tag{11.7}
\]

The first factor is nonzero over \(k(U,t,w)\).  If the degree-six binary
Hessian discriminant \(\Delta\) is squarefree, the same rank-one adjugate
argument as in Section 8 gives

\[
 \Delta\mid\operatorname{adj}(C)d.                  \tag{11.8}
\]

But \(\operatorname{adj}(C)d\) has degree \(3+2=5<6\), so it vanishes.
The adjugate identity then gives \(d=0\), and (11.7) gives
\(D_w^2h_3=0\).  Since also \(D_w^2h_5=0\), theorem `HC4CD5` closes this
branch.  The squarefree locus is nonempty: the binary quintic

\[
 x^5+x^4y+xy^4+y^5                                  \tag{11.9}
\]

is squarefree and has a squarefree degree-six Hessian determinant.

There is also an exact normal-slice consequence on the nonsquarefree
remainder.  At the generic point of \(\mathbf P(K)\), the two active
quartics form a \((4,4)\) complete intersection

\[
 B_K=k(\mathbf P(K))[u_1,u_2]/
       (\partial_{u_1}h_5,\partial_{u_2}h_5),
 \qquad \dim B_K=16.                                \tag{11.10}
\]

Before adjoining the two kernel-gradient components, the
\(\epsilon^4\)-truncated transverse algebra has length \(64\).
Because \(h_4|_K=t^4\), one kernel-gradient component is
\(\epsilon\) times a unit at the generic point of the line.  It kills
\(\epsilon\), leaving exactly (11.10).  Since the support line has degree
one,

\[
 \boxed{\sigma_3=16.}                               \tag{11.11}
\]

Equation (11.11) is the unit-penultimate law in `PGS2`, specialized to
\((n,m,r)=(4,4,2)\).

Only one codimension-three numerical row has this leading multiplicity for
each affine degree:

\[
\begin{array}{c|c|c}
\delta&g&\sigma\\ \hline
2&(1,4,16,48,2)&(0,0,16,-2)\\
3&(1,4,16,48,3)&(0,0,16,-3).
\end{array}                                         \tag{11.12}
\]

Thus the smooth rank-two packet splits into three precise pieces:

| lower-layer condition | status |
|---|---|
| \(h_4|_K\ne0\), \(\Delta\) squarefree | excluded by (11.7)--(11.8) and `HC4CD5` |
| \(h_4|_K\ne0\), \(\Delta\) nonsquarefree | open, but only the \(\sigma_3=16\) row |
| \(h_4|_K=0\) | open synchronization stratum |

### 11.2 Rank three with isolated singularities

Return to the rank-three notation of Section 8.  Let
\([p]\in\mathbf P^2\) be a singular point of the ternary quintic and assume

\[
 \operatorname{rank}\operatorname{Hess}(h_5)(p)=2. \tag{11.13}
\]

This includes an ordinary node.  Homogeneity gives
\(\operatorname{Hess}(h_5)(p)p=0\), so the kernel is \(kp\), and symmetry
gives

\[
 \operatorname{adj}(\operatorname{Hess}(h_5)(p))
 =c\,pp^{\mathsf T},\qquad c\ne0.                   \tag{11.14}
\]

Evaluating the Schur face (8.5) at \(p\) kills its determinant term.
Euler's identity for the cubic \(s_3\) then gives

\[
 0=c\bigl(p^{\mathsf T}\nabla s_3(p)\bigr)^2
   =9c\,s_3(p)^2,
 \qquad\boxed{s_3(p)=0.}                            \tag{11.15}
\]

Therefore the Schur cubic must pass through every ordinary isolated
singular point of the top quintic.  If the Hessian rank at \(p\) is at most
one, the adjugate vanishes and this value-level argument gives no
restriction.

> **Theorem `HC4PPG8` -- Codimension-three gradient-stratum sieve.**
> In the smooth essential rank-two packet, \(h_4|_K\ne0\) synchronizes a
> constant kernel direction.  If the binary quintic Hessian determinant is
> squarefree, that branch is empty.  On its nonsquarefree remainder the
> generic transverse multiplicity is forced to \(\sigma_3=16\), hence only
> the two rows (11.12) can occur.  The branch \(h_4|_K=0\) remains open.
> In the essential rank-three packet, the Schur cubic vanishes at every
> isolated singular point where the top Hessian has rank two.  No
> codimension-three row is excluded unconditionally, because the remaining
> exceptional rank-two and rank-three packets can still feed that column.

The exact coefficient, determinant-face, degree, normal-slice, and nodal
checks are in
[`scripts/verify_hc4_codim3_gradient_strata.py`](scripts/verify_hc4_codim3_gradient_strata.py).
The independent Macaulay2 replay
[`scripts/verify_hc4_codim3_gradient_strata.m2`](scripts/verify_hc4_codim3_gradient_strata.m2)
checks the radical powers (11.5), the transverse lengths \(64\to16\), and
the nodal Hessian-rank calibration.

### 11.3 Singular binary root partitions

The codimension-two essential-rank-two packet now has a direct `PGS3`
stratification.  Let \(p\) be a root of multiplicity \(e\ge2\) of the
binary quintic \(h_5\).  On a chart centered at \(p\),

\[
 h_5=x^e u(x),\qquad u(0)\ne0.                      \tag{11.16}
\]

The two top partial derivatives have respective orders \(e-1\) and \(e\),
the latter also following from Euler.  Hence the transverse top Jacobian
algebra is

\[
 B_p\simeq K[x]/(x^{e-1}),\qquad \mu_p=e-1.         \tag{11.17}
\]

Impose the following open lower-layer condition at every repeated root:
after eliminating the leading transverse active equation, one redundant
active gradient component has initial form \(\epsilon\) times a unit.
Then the active DVR profile consists of \(e-1\) order-one torsion
summands.  The order-one law in `PGS3` gives

\[
 \lambda_p=e-1.                                    \tag{11.18}
\]

If the binary quintic has \(q\) distinct roots, summing over its repeated
roots yields

\[
 \boxed{\sigma_2=\sum_p(e_p-1)=5-q.}                \tag{11.19}
\]

An essential binary quintic has at least two distinct roots, so the
singular partitions and their atlas intersections are:

| root partition | forced \(\sigma_2\) | affine degree \(2,3\) rows |
|---|---:|---:|
| \(2+1+1+1\) | \(1\) | \(51,50\) |
| \(3+1+1\), \(2+2+1\) | \(2\) | \(44,43\) |
| \(4+1\), \(3+2\) | \(3\) | \(37,36\) |

In particular, the generic point of the binary-quintic discriminant has
partition \(2+1+1+1\).  On the active-unit stratum its numerical packet
shrinks from \(260,249\) codimension-two rows to \(51,50\).

> **Theorem `HC4PPG9` -- Binary root-partition Segre sieve.**
> For an essential-rank-two singular binary quintic top satisfying the
> active-unit condition at every repeated root, the first Segre
> multiplicity is the total repeated-root excess \(5-q\), and the only
> compatible atlas rows are those in the preceding table.  No row is
> excluded unconditionally: failure of the active-unit condition is a
> proper lower-layer torsion stratum whose constant-Hessian determinant
> faces remain to be computed.

The exact root-partition and atlas ledger is
[`scripts/verify_hc4_binary_root_partition_segre.py`](scripts/verify_hc4_binary_root_partition_segre.py).
The Macaulay2 replay
[`scripts/verify_hc4_binary_root_partition_segre.m2`](scripts/verify_hc4_binary_root_partition_segre.m2)
checks multiplicities \(e=2,3,4\) with nontrivial active lower-layer
deformations.

## 12. Consequences and next finite tests

The first experiment changes the search order.

1. Use \(\Gamma_\Psi\), not the full polar map, for affine sheet counts.
   Record the full polar multidegrees separately.
2. Discard every \(m=2\), affine-degree-two or affine-degree-three row by
   Wang's theorem.  No Hilbert-scheme search is needed for lengths \(14\)
   or \(13\).
3. Discard every \(m=3\), affine-degree-two or affine-degree-three row by
   `HC4CQ1` and Ax--Grothendieck.  No searches for lengths \(79\) or \(78\)
   are needed.
4. For \(m=4\), use the leading-Hessian matrix above.  Ranks one and two
   begin with kernel synchronization; rank three reduces to the cubic Schur
   pair (8.6), and only its nonsquarefree Hessian-discriminant locus remains.
5. Before any coefficient search, apply the support sieve (9.6): rank one
   uses only codimension two; squarefree essential rank two uses only
   codimension three; and smooth essential rank three uses only the isolated
   codimension-four row.  Theorem `HC4PPG7` now removes both
   affine-degree-two/three rows in that last packet.
6. Compute nonlinear Rees corrections only on singular top strata or after
   adjoining the lower layers.  For smooth rank two, first apply
   `HC4PPG8`: the nonzero \(h_4|_K\), squarefree-Hessian branch is closed,
   and its nonsquarefree remainder has \(\sigma_3=16\).  At ordinary
   rank-three singular points impose \(s_3(p)=0\).  For a singular binary
   rank-two top, apply `HC4PPG9`: the active-unit stratum has
   \(\sigma_2=5-q\), and only its higher-torsion failure locus requires a
   new determinant-face calculation.
7. For a concrete cotangent representative, compute the actual base ideal
   \((X_0^m,G_1,\ldots,G_4)\); the abstract quartic normalization packet
   alone cannot supply its Segre class.
8. For the Meng--Yang control, compute projective degrees by general
   source/target slices first.  Recovering \(\sigma_k\) from (2.4) is
   triangular and avoids a full symbolic normal-cone presentation.

The first proposed exclusion target is now closed:

> There is no quadratic-gradient \(\mathrm{HC}_4\) counterexample of any
> affine degree.  In the projective atlas this excludes all sixteen
> affine-degree-two/three signatures, including the length-\(13\) row.

The stronger existing mixed cubic--quartic theorem also excludes every
\(m=3\) row.  For \(m=4\), 318 and 306 numerical signatures remain for
affine degrees two and three.  The next projective frontier is the
exceptional part of the codimension-three column—\(h_4|_K=0\) in rank two,
the nonsquarefree binary-Hessian \(\sigma_3=16\) row, and lower normal
cones at isolated ternary singularities satisfying (11.15)—together with
the higher-\(X_0\)-torsion failure locus of `HC4PPG9` in codimension two.
Rank-one synchronization remains behind it.

The affine reverse-Schur intersection is now finite as well. By `HC4RSD6`,
an affine singular pivot on the open rank-three packet requires the three
degree-eight components of
\(w=\operatorname{adj}(C_3)\nabla s_3\) to have a nontrivial constant
linear relation. Thus the affine-pivot-covered part of the nonsquarefree
Schur packet is cut out by the 3-by-3 minors of the 3-by-45 coefficient
matrix of \(w\). This coverage locus should be intersected with (8.6)
before any lower-face or collision elimination. See
[`HC4_AFFINE_PIVOT_COVERAGE_GATE.md`](HC4_AFFINE_PIVOT_COVERAGE_GATE.md).
On the diagonal nonsquarefree top, the intersection is explicit: the Schur
pair has three diagonal cubic channels and affine coverage forces at least
one channel to vanish.

The zero-metric branch does not survive this diagonal calibration.  Applying
the rank-one recognition scheme `HC4MR3` through the next immutable metric
face, `HC4MR4` shows that every nonaligned diagonal Schur prolongation has
empty projective constant-null-covector scheme, including the one-channel
and two-channel boundaries.  This does not exclude a direct quintic
candidate; it shows that the diagonal packet can enter the restricted
`JC2 <=> PHC4` theorem only through a higher-rank or nonlinear pencil
direction, or after a collision-preserving rechart.

For a marked collision, `HC4RSD7` makes the affine coverage intersection
empty before lower faces: every affine pivot fiber has a ternary
constant-Hessian restriction, and `HC3` separates its tangential gradient.
The coefficient-rank locus remains relevant to representation theory, but
not to inherited affine collision transfer. Direct degree-five exclusion
and nonlinear or mixed/coisotropic pivots remain live.

## 13. Reproduction

Run:

```bash
.venv/bin/python scripts/verify_hc4_projective_polar_atlas.py
M2 --script scripts/verify_projective_polar_calibrations.m2
.venv/bin/python scripts/analyze_hc4_quintic_infinity_rees.py
M2 --script scripts/verify_hc4_quintic_infinity_rees_strata.m2
.venv/bin/python scripts/verify_hc4_rank3_vertex_colength.py
M2 --script scripts/verify_hc4_rank3_vertex_colength.m2
.venv/bin/python scripts/verify_hc4_codim3_gradient_strata.py
M2 --script scripts/verify_hc4_codim3_gradient_strata.m2
.venv/bin/python scripts/verify_hc4_binary_root_partition_segre.py
M2 --script scripts/verify_hc4_binary_root_partition_segre.m2
```

The first command checks the triangular Segre inversion, enumerates the
log-concave degree lists through \(m=4\), tests the smooth-curve numerical
rows, verifies the finite consequences of Wang's theorem and `HC4CQ1`,
checks the rank-one, rank-two, and rank-three determinant faces
(8.3)--(8.5), checks the degree and adjugate identities in the squarefree
obstruction, verifies the cotangent and Meng--Yang aggregate corrections,
and writes the JSON artifact.  The second independently computes the graph
and full-polar projective degrees in (1.4) over \(\mathbb Q\).  The third
constructs the universal top quintic, checks the Euler, Hessian, curl,
Koszul, collision-parity, support, and atlas-intersection formulas, and
writes the rank-stratum artifact.  The fourth independently certifies the
generic complete-intersection, linear-type Rees, and pure-top projective
degree calculations.  The fifth records the universal flatness/socle length
bound and its exact atlas intersection.  The sixth independently checks the
complete-intersection Hilbert function and exact Fermat and deformed local
calibrations.  The seventh checks the rank-two radical synchronization,
Schur face, squarefree witness, \(\sigma_3=16\) normal slice, and
rank-three nodal incidence.  The eighth independently replays the radical
powers and local lengths in Macaulay2.  The ninth checks the binary
root-partition formula and exact atlas intersections; the tenth
independently calibrates its active-unit lengths.

The generated artifact currently has SHA-256

```text
350bc81b4ba7ac21289d7548f6d46de6526887c4e98a0813596cf20a454b240b
```

with the command and software assumptions above.

The generated quintic infinity/Rees artifact has SHA-256

```text
51ddbf2b7c0c2b9b3f2cd7c1a8dcb4bf0fe97a3e3bce306eeee031cd5c92b99d
```

under the repository Python lock; its independent Rees replay uses
Macaulay2 1.22 with `Cremona` and `ReesAlgebra` over \(\mathbb Q\).

The generated rank-three vertex-colength artifact has SHA-256

```text
c610f57af67061d0b4eb9523cb018569a7e8220a51dbd2350b71eb7007bfe473
```

and its independent calibration uses Macaulay2 over \(\mathbb Q\).

The generated binary root-partition artifact has SHA-256

```text
09f6a57c735b2751d0f890b8cd216822001bae875fd9a5156e2a27550f8e71ad
```

and its independent active-unit calibration uses Macaulay2 over
\(\mathbb Q\).

The generated codimension-three gradient-strata artifact has SHA-256

```text
8759875cf431d18f35321631984d9120c72a2335dcae31d107fa191ae539e5a3
```

and its independent calibration uses Macaulay2 over \(\mathbb Q\).
