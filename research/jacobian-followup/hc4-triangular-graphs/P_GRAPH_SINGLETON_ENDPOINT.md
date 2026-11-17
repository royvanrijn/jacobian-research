# One fixed p-graph realizes the entire singleton-zero row orbit

14 November 2026. **HC4MYGP1** identifies the exact remaining problem in
the p-chart after arbitrary nonlinear dual automorphisms. It is the
same singleton-zero plane Keller problem already retained in
[HC4MYOR1/2](../hc4-relative-row-orbit/PROOF.md), with an explicit genuine
graph realization. This does not solve that properness problem or JC2.

## 1. Every remaining p-chart coefficient pair has a singleton zero

The [rational-curve theorem](CURVE_TRACE.md) already proves that a
nonautomorphic unit-Hessian p-chart potential, after any base and dual
automorphisms preserving the fibration, must be globally affine-dual:

\[
 f=b(s,t)+P G_1(s,t)+Q G_2(s,t),\qquad
 \operatorname{Jac}(G_1,G_2)\in K^\times.           \tag{1}
\]

Its zero coefficient fibre is exactly the base point corresponding
to old (x,y)=(0,0). To prove this, fix a base point and invert the
specialized dual automorphism, so the retained ambient q,r are the
two fibre coordinates and p is a polynomial graph in them.

On u=1+xy=0 the graph disappears and the potential has nonzero q
coefficient 2N/x, so it cannot be constant on that fibre. On u nonzero,
the defining formulas give the complete identity

\[
 \Psi=L A^2+
 \left(M+\frac{Ny^2(3u+1)}{u^2}\right)A
 +Nyq+\frac{Nx(u+1)}{u^2}r.                       \tag{2}
\]

If this is constant in q,r, polynomial degree forces A to be constant:
its square otherwise has twice its positive dual degree, larger than
the other terms. Both displayed linear coefficients must then vanish.
They give y=0 and x=0. Conversely, at the origin use the parameter
dual coordinates of (1): the potential is Lp(0,P,Q)^2+Mp(0,P,Q).
Affineness in those same coordinates forces p(0,P,Q) to be constant
by polynomial degree, so the entire fibre potential is constant.
Constancy, unlike affineness, is preserved by the nonlinear dual
inverse. Thus G has exactly one zero over the algebraic closure.

Translate that base point to zero by a constant affine change. The
unit Jacobian makes its local zero ideal the maximal ideal with length
one. The quotient by (G1,G2) is a zero-dimensional algebra supported
only at that point; consequently

\[
 (G_1,G_2)=(s,t).                                \tag{3}
\]

This also descends to the original field by faithful flatness of field
extension. Gradient invertibility in (1) is exactly invertibility of G.
After a constant target normalization, its Jacobian can be set to
-2N^2, so it lies in the oriented row orbit constructed below.

## 2. An explicit fixed polynomial graph and polynomial frame

Return to the original base x,y and put v=xy, u=1+v. Fix the single graph

\[
 p=p_*=-2x(v^2+1)q-2x^3r.                        \tag{4}
\]

This graph has total degree six. Define the polynomial matrix

\[
 D=\begin{pmatrix}
 u^2(1-2v-2v^2-2v^3)&-x^2(2v^3+6v^2+6v+3)\\
 -y^2(2v^4+4v^3+3v^2-2)&1-v^2-4v^3-2v^4
 \end{pmatrix}.                                  \tag{5}
\]

It has determinant one. The original Meng--Yang linear forms on (4)
obey the literal polynomial identities

\[
 (t,w)^T=D(q,r)^T,\qquad A=xt,\qquad B=yt+2xw.    \tag{6}
\]

No localization at x, y or u is used: D and D^(-1) are polynomial.
Replace t=T and w=W-L*x*T^2/(2N). This is a polynomial dual shear
with a polynomial inverse. Equations (6) now give

\[
 \Psi\circ i_0=(Mx+Ny)T+2NxW.                    \tag{7}
\]

Thus a nonlinear dual parameterization of one fixed embedded graph
produces an affine-dual unit-Hessian potential. Its determinant is
4N^4 and its gradient is a linear automorphism.

## 3. Every oriented row is realized by the same embedded graph

For any C(x,y) in SL2(K[x,y]), set (T,W)^T=C(P,Q)^T before the shear.
The explicit parameterization is

\[
 (q,r)^T=D^{-1}
 \begin{pmatrix}T\\W-LxT^2/(2N)\end{pmatrix},
 \qquad p=p_*(x,y,q,r).                          \tag{8}
\]

Its retained-dual map is a polynomial automorphism over K[x,y]: invert
D, the shear and C successively. The image is exactly the fixed graph
(4), and its parameter inverse is polynomial on that graph. This is
a genuine polynomial closed embedding, not merely a Keller
parameterization or a map defined away from a divisor.

The actual pullback is

\[
 f_C=G_C\begin{pmatrix}P\\Q\end{pmatrix},\qquad
 G_C=(Mx+Ny,2Nx)C.                               \tag{9}
\]

The row theorem HC4MYOR1 identifies this orbit exactly as

\[
 \{G:(G_1,G_2)=(x,y),\ J(G)(0)=-2N^2\}.           \tag{10}
\]

Indeed let ell=(y,2x) and
H0=[[1,0],[M/(2N),1]]. Then (Mx+Ny,2Nx)=N*ell*H0 and det(H0)=1.
The retained oriented-row theorem for ell applies after division by N.
No inverse of G is assumed in this construction.

For every row in (10),

\[
 \det\operatorname{Hess}f_C=J(G_C)^2,
 \quad J(G_C)(0)=-2N^2.                          \tag{11}
\]

Thus unit Hessian is equivalent to G_C being Keller, and then its
determinant is exactly 4N^4. Its gradient is polynomially invertible
if and only if G_C is. The cotangent factorization proves both
directions, including collision lifting when the plane pair collides.

## 4. Exact equivalence and the remaining gate

Together, Sections 1 and 3 prove **HC4MYGP1**: universal gradient
invertibility for p-graphs after arbitrary base and retained-dual
automorphisms is equivalent, over C, to invertibility of every plane
Keller pair generating the origin ideal, up to the harmless orientation
normalization in (10). One fixed degree-six graph already realizes
the full possible coefficient orbit.

This is exactly the existing singleton-zero properness problem,
not an all-plane-JC2 equivalence claim. A nonautomorphic such pair
would be nonproper at its zero target; its finite normalization must
supply the omitted fibre length described in HC4MYOR2. Neither the
singleton affine fibre nor ideal equality proves that completeness.
Conversely a solution of that properness problem would close the entire
remaining nonlinear-dual p-chart endpoint.

Increasing the graph degree or dual word length does not evade this
identification. The graphs and reparameterizations in (8) are already
polynomial with explicit inverses. The earlier standard/affine-dual
p-chart obstruction remains valid: the shear used here is nonlinear,
and a nonlinear parameterization need not preserve a Hessian determinant.

Every parameterization (8) contains the known marked v2 points: zero
dual parameters give q=r=p=0. Their ambient collision transfers exactly
when G_C has equal values at the two marked base points, not for every C.
The general plane-map equivalence concerns arbitrary gradient collisions;
it does not assert that every singleton-zero Keller pair has a collision
at that fixed marked pair.

The [checker](verify_square_axis.py) recomputes (2), the frame (5)-(6),
the shear cancellation, inverse formulas, row orientation and Hessian
determinants from the original five-variable formulas. Polynomial
degree, singleton-ideal descent and the complete orbit classification
are written arguments, with HC4MYOR1 explicitly imported. The
reduction of non-affine-dual p potentials imports HC4MYGC3.
[Actual receipt](receipt-square-axis.json).

```sh
timeout 60s .venv/bin/python research/jacobian-followup/hc4-triangular-graphs/verify_square_axis.py
```

The q-chart extension is proved separately by the
[critical-coordinate first-jet theorem](CRITICAL_COORDINATE_JET.md),
which covers all polynomial base and retained-dual automorphisms.
Both native base-eliminated charts are subsequently closed by
[HC4MYGB2](BASE_ELIMINATED_GRAPHS.md). Nonlinear reparameterizations
of those charts and genuine base/dual mixing remain outside these
results. No general HC4 or JC2 solution, formal
verification, independent external review or novelty priority is claimed.
