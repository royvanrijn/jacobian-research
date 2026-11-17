# An isolated trace point closes a larger mixing class

14 November 2026. The previous exceptional-fibre and cylindrical-trace
arguments needed information about generic base fibres. The isolated-point
criterion below removes that requirement. It closes the standard r graph
after **every** polynomial automorphism of `(y,p,q)` over K[x], allowing
old y to mix with p,q. It also closes arbitrary polynomial substitutions
whose restriction to one source hyperplane has the specified coordinate
triple, without a global inverse or embedding assumption.

All results use the retained Tang structural alternative and the retained
affine-critical-component theorem. There is no degree or word-length bound.
The unrestricted HC4 problem and both plane Keller endpoints remain open.

## 1. An isolated critical point gives a global trace bound (HC4MYHB0)

Let f(s,z1,z2,z3) be a complex four-variable polynomial with nonzero
constant Hessian determinant. Suppose its gradient is not an automorphism.
Let h=f(s0,-) have an isolated critical point a with Milnor number

\[
 \mu=\dim_C\mathcal O_{A^3,a}/(h_{z_1},h_{z_2},h_{z_3})\ge1.
\]

Here the dimension means the finite local length; analytic completion
gives the same number. Then the generic value fibre of this **same trace**
satisfies

\[
                         \chi(h=\lambda)\ge\mu+2.             \tag{1}
\]

Euler characteristic is compactly supported complex Euler characteristic.
Thus, if one such trace point has generic trace Euler at most mu+1,
unit Hessian forces a polynomial gradient inverse. No assumption about
generic base fibres, finite critical projection or genus is imposed.

### The whole critical curve and dominance

Write C_f=V(f_z1,f_z2,f_z3). Under unit Hessian it is smooth of pure
dimension one whenever nonempty: its three Jacobian rows are rows of
the invertible full Hessian. The trace critical scheme is exactly
`C_f intersect {s=s0}`. At the isolated point a, the component of C_f
through a cannot be contained in that hyperplane. It therefore projects
nonconstantly to s. In particular a dominating critical component exists.

The [HC4MYGR1 affine-component theorem](CRITICAL_COMPONENTS.md), using
its dominating-component alternative, now excludes **every** actual A1
component of C_f under hypothetical nonautomorphy. Each smooth affine
curve of positive Euler is A1: its completion has genus zero and exactly
one omitted point. Thus

\[
                              \chi(C_f)\le0.                  \tag{2}
\]

Vertical components are allowed in this conclusion. They cannot be
discarded from the whole-curve count.

### Both Tang ranks and the local intersection length

The retained [HCTANG2 alternative](../tang-hc4-audit/STRUCTURAL_REVIEW.md)
supplies a constant jointly affine source two-plane. If its projection
to s has rank zero, the fibration-preserving constant affine normal form is

\[
 f=D(s,U)+P F(s,U)+Q G(s,U),\qquad J(F,G)\in C^\times.
\]

All relative critical points lie over the finite common zero set of F,G.
Consequently every critical component has constant base value, contrary
to the dominating component just proved. This excludes rank zero without
calculating generic-base Euler.

In rank one the corresponding form is

\[
 f=E(U,Z)+sA(U,Z)+VB(U,Z),\qquad J(A,B)\in C^\times.        \tag{3}
\]

The **whole** C_f is isomorphic to C=B=0: solve the two remaining
critical equations for s,V using the unit determinant J(A,B). Denote
the resulting base function on C by sigma. On C,

\[
 dE=-\sigma\,dA,\qquad d(E+s_0A)=(s_0-\sigma)\,dA.          \tag{4}
\]

The Keller condition makes dA nowhere zero. The local ring of the smooth
C_f at a is a DVR. Since the trace Jacobian algebra is its quotient by
s-s0, isolation and the definition of mu give

\[
                   \operatorname{ord}_a(\sigma-s_0)=\mu.     \tag{5}
\]

In a local parameter on C, equation (4) therefore says that
`H0=E+s0*A` has local mapping degree mu+1 at a. Characteristic zero
ensures that integration increases the vanishing order by exactly one.
The function H0 is nonconstant on this component, and its generic degree
on the entire C, summed over nonconstant components, is at least mu+1.
Constant H0 components contribute zero to this generic count; their
finitely many values are excluded when choosing lambda.

### Compare the same entire trace fibre

At exactly s=s0 in (3), off B=0 the value equation solves V and gives
A2 minus C. Over B=0 each of the d0 generic points H0=lambda contributes
an affine V line. Therefore

\[
 \chi(h=\lambda)=1-\chi(C)+d_0
                \ge1+0+(\mu+1)=\mu+2.                       \tag{6}
\]

This proves (1). A compactification of C only interprets a nonconstant
function's degree; no boundary point is silently counted as affine.
For generic lambda all finite boundary values are avoided, so the
generic affine count equals that degree. No generic-rational-fibre theorem
or assumption that a special fibre is generic is used.

Actual polynomial retained-dual automorphisms over C[s] preserve the
whole critical scheme, its intersection length with s=s0 and the entire
value fibre by [HC4MYGK0](RELATIVE_CRITICAL_LOCUS.md). A polynomial
equivalence of the three-variable trace also preserves its isolated
critical point, Milnor number and Euler function. The bound therefore
applies throughout those orbits, without nonlinear invariance of the
four-variable Hessian.

## 2. Every power-plus-saddle trace is sufficient (HC4MYHB1)

Suppose at one scalar s0 the trace of an arbitrary polynomial f is
polynomially equivalent, with an actual inverse, to

\[
                   a(U-u_0)^d+bVW+c,
              \qquad ab\ne0,\quad d\ge2.                     \tag{7}
\]

It has exactly one critical point, of Milnor number d-1. Indeed its
Jacobian algebra has basis `1,(U-u0),...,(U-u0)^(d-2)`. Its generic
value fibre has Euler d: on V nonzero solve W with U free, giving
C* times A1 of Euler zero; on V=0, a generic value gives d distinct
U values, each with W free. The critical-value fibre has Euler one.

The generic value d is one less than the lower bound d+1 from HC4MYHB0.
Hence **unit Hessian forces gradient invertibility whenever (7) holds**.
This includes every nondegenerate three-variable quadratic trace, by
constant affine equivalence. All polynomial off-hyperplane terms are
unrestricted. The source hyperplane must be affine in the actual Hessian
variables; no nonlinear straightening of a source divisor is assumed.

The conclusion is invertibility, not absence of unit-Hessian members.
For example `f=sU+aU^d+bVW` has Hessian determinant b^2 and a direct
triangular gradient inverse for every d>=2.

## 3. The standard r graph after full three-variable mixing (HC4MYHB2)

For the scaled Meng--Yang potential Psi=L*A^2+M*A+N*B, LN nonzero,
take **any** polynomial graph r=R(x,y,p,q). Its full x=0 trace is
independent of the graph:

\[
             f(0,y,p,q)=Lp^2+Mp+4Npy^2+Nyq.                 \tag{8}
\]

Use the actual polynomial coordinate change of this three-variable trace

\[
      U=p+M/(2L),\quad V=y,\quad W=q+4py.
\]

It has polynomial inverse and carries (8) to

\[
                         LU^2+NVW-M^2/(4L).                  \tag{9}
\]

Its critical point is isolated, with mu=1 and relative Hessian determinant
-2*L*N^2. Its generic value fibre has Euler two. Now allow **every actual
polynomial automorphism of (y,p,q) over K[x]** before pulling back the
graph potential. Its x=0 specialization is still an actual three-variable
automorphism, so its trace remains polynomially equivalent to (9).
HC4MYHB1 proves polynomial gradient invertibility whenever the actual
four-variable Hessian is a nonzero constant.

This permits old y to depend on retained dual variables, extending the
earlier two-base-preserving r-chart closure. Only x is fixed here.
Degrees, mixed terms and word lengths are unrestricted. A graph followed
by such an automorphism is an actual polynomial embedding: its parameter
inverse is the projection to `(x,y,p,q)` followed by the retained
automorphism inverse. A marked collision would still need actual
pulled-back gradient-covector equality; none is asserted.

## 4. Only a hyperplane coordinate triple is needed (HC4MYHB3)

More generally, let `(s,z1,z2,z3)` be actual affine source coordinates,
and substitute arbitrary polynomials `(X,Y,P,Q,R)` into Psi. Require only

\[
 X(0,z)=0,\qquad
       (Y(0,z),P(0,z),Q(0,z))\text{ is an automorphism of }A^3. \tag{10}
\]

No global Jacobian or inverse condition is imposed on the substitutions;
R and all off-hyperplane terms are arbitrary. The actual pullback trace
is (8) composed with that coordinate triple, and hence is equivalent to
(9). Therefore **every unit-Hessian pullback satisfying (10) has a
polynomial gradient inverse**.

In particular one may take `X=s*E(s,z)` for arbitrary E, and arbitrary
s-divisible corrections to any polynomial trace coordinate triple. This
includes nonlinear changes of the old free base off the hyperplane.
The statement concerns the actual pulled-back potential; it does not
assert that these broader substitutions are embeddings.

As a concrete extension to the remaining y chart, take a native graph

\[
 y=R(x,p,q,r),\qquad R(0,p,q,r)=a r+b(p,q),\quad a\in K^\times. \tag{11}
\]

Then `(R(0,p,q,r),p,q)` is a trace coordinate triple, so (10) applies.
It continues to apply after **every** actual retained-three-dual
automorphism over K[x]. Thus this whole genuinely old-r-dependent
boundary class is closed, with arbitrary r dependence away from x=0.
The exact native trace-coordinate condition is (11): an automorphism of
K[p,q][r] fixing p,q must have r degree one and unit constant coefficient.

The same trace test applies to p or q graphs when their omitted-dual
trace is a coordinate complement to the other two entries of (Y,P,Q).
No merely Keller triple is substituted for an actual coordinate triple.
Consequently the existing unit-ideal and maximal-ideal plane endpoints
are not resolved by this criterion.

## Scope, verification and remaining work

First prove these statements over C. In characteristic zero, a candidate
and any required trace/retained inverse use finitely many coefficients;
embed that coefficient field into C. A polynomial gradient inverse then
descends by uniqueness, as in the preceding all-degree chart theorems.

The remaining old-r-dependent y bases must avoid (11), and other opposite
q/r-dependent bases, mixing outside (10), uncovered embeddings and both
distinct plane Keller endpoints remain open. This does not establish
unrestricted HC4 or JC2, formal verification, external review or novelty.

[verify_isolated_trace.py](verify_isolated_trace.py) and its
[receipt](receipt-isolated-trace.json) check the exact native trace,
coordinate inverses, Milnor/ramification identities, complete rank-one
reconstruction, rank-zero critical lines, and nonlinear mixing controls.
The whole-curve, local-length and unrestricted-degree arguments are
written proofs. Their external inputs are retained Tang and HC4MYGR1's
explicit classical plane-line inputs; no new classification is imported.
