# All five native coordinate-graph charts are closed in every degree

> The later [relative-coordinate theorem](RELATIVE_COORDINATE_ENDPOINT.md),
> `HC4MYGU0/1/2/3`, forbids unit Hessian after arbitrary retained-dual
> substitutions for a zero boundary slice. For every origin-avoiding
> base-only graph, all retained-dual automorphisms reduce exactly to an
> omitted-point plane Keller problem, with a fixed x=1 embedded realization.
> This does not solve that plane existence problem.

14 November 2026. **HC4MYGB2** closes the two remaining native charts
x=R(y,p,q,r) and y=R(x,p,q,r). Together with the earlier p,q,r results,
every polynomial graph eliminating one of the five original Meng--Yang
coordinates has an invertible gradient whenever its actual four-variable
Hessian determinant is a nonzero constant. Graph degrees are unrestricted.

The proof uses the [reproduced Tang alternative](../tang-hc4-audit/STRUCTURAL_REVIEW.md).
The y-chart argument also uses the retained
[Gwozdziewicz line theorem HC4MYGL0](CRITICAL_COORDINATE_JET.md).
The x-chart plane exclusion and the y-chart trace and rank-one arguments
are elementary. These are native graph results, stable under constant
affine parameter changes; arbitrary nonlinear reparameterizations of
the two base-eliminated charts are not asserted to preserve them.

Put H(z)=Lz^2+Mz, with LN nonzero, and use the original linear forms

\[
 A=(1+xy)^3p+3x(1+xy)^2q-x^3r,
\]
\[
 B=y^2(1+xy)(4+3xy)p+(y+3xy^2(4+3xy))q+(2x-3x^2y)r,
 \qquad \Psi=H(A)+NB.
\]

## 1. Classify every affine direction of the two boundary traces

An affine direction means a nonzero constant vector v in the three
retained dual variables with partial_v^2 h=0 identically. It does not
mean that the full four-variable potential has an affine two-plane.

**HC4MYGB0.** The following are exact characterizations.

For a graph y=R(x,p,q,r), write R0=R(0,p,q,r). Its x=0 trace is

\[
 h_y=H(p)+4NpR_0^2+NqR_0.                        \tag{1}
\]

Then v is an affine direction if and only if

\[
 v(p)=0,\qquad v(R_0)=0.                        \tag{2}
\]

To prove necessity choose a linear variable V along v. If p has
nonzero V coefficient and deg_V R0=m>=1, the term 4NpR0^2 has degree
2m+1, strictly larger than every other term of (1). If m=0, H(p)
has the uncancellable quadratic coefficient L*v(p)^2. Hence v(p)=0.
Now p is a nonconstant transverse linear polynomial. If m>=2,
4NpR0^2 again has strictly largest V degree. If m=1 and k is R0's
nonzero V coefficient, the quadratic coefficient is

\[
 Nk(4pk+v(q)).                                  \tag{3}
\]

It cannot vanish: 4pk=-v(q), a field constant, would make the nonunit
polynomial p divide a nonzero constant, or would force k=0.
Thus m=0. Conversely (2) makes (1) affine along v immediately.

For a graph x=R(y,p,q,r), write R0=R(0,p,q,r). Its y=0 trace is

\[
 h_x=H(p+3R_0q-R_0^3r)+2NR_0r.                  \tag{4}
\]

It has any nonzero affine direction only if R0 is a field constant c.
In that case its affine directions are exactly

\[
 v(p)+3c\,v(q)-c^3v(r)=0.                       \tag{5}
\]

Indeed, if deg_V R0=m>=1, the cubic term -R0^3r has degree 3m or
3m+1, according as v(r)=0 or not. Its leading coefficient is nonzero:
in the first case r is a nonzero transverse linear polynomial. The
other terms in its argument have degree at most m+1. The square in
H therefore has degree at least 6m and cannot cancel against its
linear term or 2NR0r. Hence R0 is V-independent. The quadratic
coefficient of (4) is then

\[
 L\bigl(v(p)+3R_0v(q)-R_0^3v(r)\bigr)^2.
\]

If R0 is nonconstant, it is transcendental over the coefficient field;
the displayed constant-coefficient cubic relation forces all three
components of v to vanish. This is impossible. For constant R0=c,
(5) is also sufficient. This proves HC4MYGB0 in every degree.

## 2. The x-eliminated chart has no surviving unit-Hessian affine plane

Let f=Psi(R(y,p,q,r),y,p,q,r) have a jointly affine constant two-plane W
and unit Hessian. Project W to the one-dimensional free base y.
The rank is zero or one, and W always contains a nonzero pure-dual
direction. Its restriction to y=0 is an affine direction of (4).
HC4MYGB0 therefore forces R0=c constant.

Set a=p+3cq-c^3r, a constant linear coordinate. The trace is

\[
 f|_{y=0}=H(a)+2Ncr.                            \tag{6}
\]

If the projection rank is one, constant affine coordinates put W in
the span of partial_y and a pure-dual partial_V, with

\[
 f=D(U,Z)+yF(U,Z)+VC(U,Z),\qquad
 \det\operatorname{Hess}f=\operatorname{Jac}(F,C)^2. \tag{7}
\]

The base-linear dual translation used here agrees with the original
duals at y=0. Equation (5) gives partial_V a=0, so (6) makes
C=partial_V f|y=0=2Nc*partial_V r a field constant. This contradicts
the unit Jacobian in (7).

If the rank is zero, W is a pure-dual two-plane. All its directions
annihilate a, so W=ker(da). Complete a to constant linear dual
coordinates a,U,V. Joint affineness gives

\[
 f=D(y,a)+U F(y,a)+V G(y,a),\qquad
 \det\operatorname{Hess}f=\operatorname{Jac}_{y,a}(F,G)^2. \tag{8}
\]

Equation (6) makes both F(0,a) and G(0,a) field constants. Their
a derivatives vanish there, so the Jacobian is zero on y=0,
contradicting unit Hessian. This excludes every W in the x chart.

## 3. The y-eliminated chart: pure-dual planes

Now f=Psi(x,R(x,p,q,r),p,q,r), and project W to the free base x.
If its rank is zero, HC4MYGB0 makes every direction of W annihilate
p and R0. Since dim W=2, W=ker(dp), so R0 depends only on p and

\[
 f=D(x,p)+qF(x,p)+rG(x,p),\qquad
 F(0,p)=NR_0(p),\quad G(0,p)=0.                 \tag{9}
\]

Unit Hessian is equivalent to a nonzero constant Jacobian for (F,G).
At x=0 it gives

\[
 -N R_0'(p)G_x(0,p)\in K^\times.                \tag{10}
\]

Both polynomial factors are nonzero constants, so R0=kappa*p+delta,
kappa nonzero. The plane map (F,G) is injective on x=0, where it is
(N*(kappa*p+delta),0). HC4MYGL0 makes it an automorphism.
The cotangent reconstruction then gives a polynomial inverse of grad f.

## 4. The y-eliminated chart: a moving base direction

If the projection rank is one, choose its pure-dual generator v.
HC4MYGB0 gives v(p)=v(R0)=0. If v(q)=0 as well, the coefficient of
this affine direction in f is zero on x=0 and hence identically zero
in the jointly affine normal form; its Hessian would be singular.
Normalize v(q)=1 and write v(r)=c.

The other generator has x component one; write its dual components
as (alpha,beta,gamma). The constant linear coordinates

\[
 p=U+\alpha x,\qquad q=V+\beta x,\qquad
 r=Z+cV+\gamma x
\]

make W=span(partial_x,partial_V). Put Y for the graph function in
these adapted variables. Then R=Y(0,U,V,Z) is independent of V and

\[
 f=D(U,Z)+xF(U,Z)+VC(U,Z),\quad
 D=H(U)+4NUR^2,\quad C=NR.                      \tag{11}
\]

Let E=Y_x|x=0, initially an arbitrary polynomial in U,V,Z. The full
first x jet from the original potential is

\[
 F=V\{3H'(U)+12NR^2+2Nc\}+N(V+8UR)E
       +2NZ+3UH'(U)R+7NUR^3
       +\alpha\{H'(U)+4NR^2\}+\beta NR.          \tag{12}
\]

The left side is V-independent. If E had positive V degree, the
monic factor V+8UR would give an uncancellable degree at least two.
Thus E is V-independent and its remaining linear coefficient forces

\[
 E=-3H'(U)/N-12R^2-2c.                          \tag{13}
\]

Substituting into (12) gives the exact normal coefficient

\[
 F=2NZ-21UH'(U)R-89NUR^3-16NcUR
          +\alpha\{H'(U)+4NR^2\}+\beta NR.       \tag{14}
\]

Consequently

\[
 \operatorname{Jac}_{U,Z}(F,C)
 =N\left[\{2L\alpha-(84LU+21M+16Nc)R-89NR^3\}R_Z
                          -2NR_U\right].       \tag{15}
\]

If m=deg_Z R>=1 with leading coefficient k(U), the uniquely highest
Z degree in (15) is 4m-1 from the cubic term, with coefficient
-89*N^2*m*k(U)^4. It cannot vanish in characteristic zero. All other
degrees are at most 2m-1 or m. A constant nonzero Jacobian therefore
forces R_Z=0, and then (15) forces

\[
 R=\kappa U+\delta,\qquad \kappa\ne0.           \tag{16}
\]

Now C=N*(kappa U+delta) and F=2NZ+F0(U). This is a triangular
automorphism, so grad f is explicitly invertible. In coordinates
(x,V,U,Z), recover U from f_V, Z from f_x, x from f_Z=2Nx, and V
from f_U after subtracting D'(U)+x F0'(U).

Sections 2-4 prove **HC4MYGB1**: under the jointly affine two-plane
hypothesis every actual unit-Hessian native base graph has an
invertible gradient; the x chart in fact has no such plane. No degree
or coefficient bounds occur in the trace and leading-degree arguments.

## 5. Global consequence and remaining scope

For a hypothetical nonautomorphic unit-Hessian gradient, HCTANG2
supplies the plane just classified. HC4MYGB1 contradicts nonautomorphy.
Thus **HC4MYGB2** closes both native base-elimination charts in every
degree and over every characteristic-zero coefficient field, by
scalar extension and polynomial-inverse descent.

Together with [HC4MYGL2](CRITICAL_COORDINATE_JET.md) for q,r and
[HC4MYGC3](CURVE_TRACE.md) for p, all five native coordinate-graph
charts are now closed as gradient-counterexample routes. The theorem
does not assert that every such unit-Hessian graph is nonexistent.
It does not require the marked collision to establish invertibility.

Arbitrary nonlinear reparameterizations need not preserve either native
affine directions or a Hessian determinant. HC4MYGU0/1/2/3 separately
handles zero boundary slices and reduces all base-only graph
reparameterizations to the exact omitted-point Keller endpoint.
General dual-dependent base-graph reparameterizations remain outside
those results. The full q,r
base/dual automorphism classes are separately closed. The full nonlinear
p-dual class is still exactly the unresolved singleton-zero Keller
problem in [HC4MYGP1](P_GRAPH_SINGLETON_ENDPOINT.md). Genuine nonlinear
base/dual mixing and uncovered embeddings also remain.

The [checker](verify_base_graphs.py) recomputes both trace identities,
the full first x jet, the forced coefficients -21,-89 and the residual
Jacobian, both projection normal-form determinants, singular x-chart
controls, and triangular gradient inverse compositions. The unbounded
degree and plane-exhaustion proofs are written arguments, not finite
coefficient searches. The line theorem and Tang remain explicit imports.
[Actual receipt](receipt-base-graphs.json).

```sh
timeout 60s .venv/bin/python research/hessian/hc4-triangular-graphs/verify_base_graphs.py
```

No unrestricted HC4 or JC2 solution, singleton-zero properness theorem,
formal verification, independent external review or novelty priority
is claimed.
