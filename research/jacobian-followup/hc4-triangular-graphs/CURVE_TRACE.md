# Rational base-curve traces and two all-degree graph charts

> The later [square-axis theorem](Q_GRAPH_SQUARE_AXIS.md), `HC4MYGQ1/2`,
> closes native q graphs and their stated affine-dual/base extension.
> [HC4MYGP1](P_GRAPH_SINGLETON_ENDPOINT.md) now identifies the remaining
> nonlinear-dual p route exactly with the existing singleton-zero plane
> Keller problem, including a fixed degree-six embedded realization.
> The graph-independent q-curve limitation below remains valid, but is
> no longer an unresolved native-q determinant gate.
> The later [critical-coordinate first jet](CRITICAL_COORDINATE_JET.md),
> `HC4MYGL1/2`, closes q graphs under all polynomial base and retained-dual
> automorphisms, with Gwozdziewicz's line theorem as an explicit import.

14 November 2026. **HC4MYGC2** permits every polynomial automorphism of
the two base variables, followed by every polynomial dual automorphism,
in the standard r-graph chart. **HC4MYGC3** closes the standard chart
eliminating p, including arbitrary base automorphisms and relative affine
dual changes. Both follow from the rational-curve trace criterion
**HC4MYGC1** below. Degrees and word lengths are unrestricted.

The external inputs remain the reproduced [Tang alternative](../tang-hc4-audit/STRUCTURAL_REVIEW.md)
and the classical Rentschler endpoint retained in
[HC4MYGV0](AUTOMORPHIC_TRACE.md#2-the-hamiltonian-pencil-lemma).
No plane Jacobian conjecture or bounded degree theorem is assumed.

## 1. A rational-curve trace criterion

Work first over C. Write the four parameters as (s,t,P,Q), separating
the two base parameters from the two dual parameters. Let

\[
 \gamma(z)=(s(z),t(z))\in C(z)^2
\]

be a nonconstant rational base curve. Assume the trace

\[
 h(z,P,Q)=f(s(z),t(z),P,Q)
\]

is a coordinate of C(z)[P,Q]. Then the following conclusions hold.

**HC4MYGC1(a).** If gamma is regular at z=0 and

\[
 h(0,P,Q)=L a(P,Q)^2+M a(P,Q),\qquad L\ne0,
 \quad a\text{ nonconstant},                       \tag{1}
\]

every nonzero-constant-Hessian f has a polynomially invertible gradient.

**HC4MYGC1(b).** If the image of gamma is not contained in an affine
line, every nonzero-constant-Hessian f either has an invertible gradient
or is globally affine in P,Q:

\[
 f=b(s,t)+P A(s,t)+Q B(s,t).                       \tag{2}
\]

In particular its gradient is invertible if it has any nonlinear dual
dependence. In the remaining affine-dual case, unit Hessian is equivalent
to a constant nonzero plane Jacobian for (A,B), and gradient invertibility
is exactly invertibility of that plane map. Part (b) is a reduction on
the affine-dual branch, not a proof that every such branch is invertible.

All conclusions descend to characteristic-zero coefficient fields by
the usual finite-coefficient-field embedding and unique formal inverse.
For the dichotomy, if f is not affine-dual over the original field it
remains so after scalar extension, giving the invertibility conclusion;
otherwise the cotangent equivalence already holds over that field.

## 2. Replace a leading polynomial coefficient by a rational pole

Suppose A(U,V), B(U,V) have constant nonzero Jacobian and

\[
 b(U,V)+s(z)A(U,V)+t(z)B(U,V)
\]

is a coordinate over C(z). Every nonconstant rational base curve has
a pole somewhere on the projective parameter line. Choose a local
parameter epsilon there and let d>0 be the largest pole order of
s,t. Multiplying the trace by epsilon^d preserves its coordinatehood
over C(z). Its coefficients are regular at epsilon=0 and its value
there is

\[
 g=cA+eB,\qquad (c,e)\ne(0,0).                    \tag{3}
\]

The Hamiltonian of the scaled trace is locally nilpotent. Choose an
annihilating exponent for each generator U,V. Reducing those polynomial
derivation identities at epsilon=0 makes D_g locally nilpotent. The
derivation acts only on U,V, so reduction commutes with every iterate;
no commutativity of its coefficient derivations is assumed.

The constant-Jacobian pair makes every nonzero linear combination g
have a constant-Jacobian mate, and prevents (3) from being constant.
The Rentschler argument from HC4MYGV0 therefore makes that mate pair,
and hence (A,B), a polynomial automorphism. This is the rational-pole
version of the earlier polynomial Hamiltonian pencil argument. If the
pole is at infinity, epsilon=1/z gives the usual highest-degree term.

## 3. Exhaust planes by their base projection

Tang supplies a jointly affine constant two-plane W for any hypothetical
nonautomorphic unit-Hessian gradient. Let r be the rank of the projection
of W to the base (s,t) plane. This rank is 0, 1 or 2.

### Rank two

Constant linear changes of the dual variables by base-linear terms make
W the base plane. The potential then has the exact form

\[
 f=b(U,V)+sA(U,V)+tB(U,V).
\]

The trace's dual change is a translation over C(z), so it remains a
coordinate. Unit Hessian gives Jac(A,B) in C^times. Section 2 makes
the plane map invertible, and the usual cotangent inverse recovers
U,V from A,B and s,t using D(A,B)^T. Thus this rank cannot support a
gradient counterexample.

### Rank one

Constant affine base coordinates S,T and constant linear dual coordinates
U,V, with a base-linear dual translation if needed, put W in the span
of partial_S and partial_V. Hence

\[
 f=b(T,U)+S A(T,U)+V B(T,U).                       \tag{4}
\]

The trace is still a coordinate in U,V over C(z). Since it is affine
in V, its V coefficient B(T(z),U) is independent of U. This follows
from the elementary coordinate-Hamiltonian argument in HC4MYGV1:
on C(z)[U], its locally nilpotent Hamiltonian acts as
-B(T(z),U) partial_U, forcing that coefficient to be constant in U.

If T(z) is nonconstant, substitution T->T(z) is injective, so B_U=0
identically. Unit Hessian gives

\[
 \det\operatorname{Hess}f=(A_U B_T)^2\in C^\times.
\]

Both factors are nonzero constants; (A,B) is triangular and the gradient
is invertible. If gamma is not contained in a line, this is the only
rank-one case, proving that part of (b).

If T(z)=T0 is constant, the curve lies in the affine line T=T0.
The coefficient B(T0,U) is then a field constant independent of z.
Under hypothesis (1), affine V dependence of the special trace forces
a_V=0 and B(T0,U)=0 by polynomial degree. The entire trace therefore
depends only on U and must be affine in U over C(z). Its coefficients
are regular at z=0, so its U degree cannot increase there. But (1),
now a nonconstant polynomial in U, has degree at least two. This is a
contradiction. Thus (a) also excludes the flat rank-one case. No
embedding or linear parameterization of the rational curve is assumed.

### Rank zero

Here W is exactly the pure-dual plane, so f has form (2). Hypothesis
(1) excludes it because a nonconstant square cannot be affine in P,Q.
Without (1), this is precisely the possible residual branch of (b).
The determinant is Jac(A,B)^2. A polynomial source change in the dual
variables factors its gradient as (A,B) times identity, giving the
claimed equivalence without assuming JC2. This proves HC4MYGC1.

## 4. Every base and dual automorphism in the standard r chart

**HC4MYGC2.** Let (x,y)=beta(s,t) be an arbitrary polynomial base
automorphism. Let (p,q)=phi(s,t,P,Q) be an arbitrary dual automorphism
over K[s,t], and let r=R(s,t,P,Q) be any polynomial. For the scaled
Meng--Yang family, a nonzero constant Hessian determinant of the actual
pullback forces its gradient to be polynomially invertible.

Take the polynomial curve

\[
 \gamma(z)=\beta^{-1}(0,z).
\]

On this curve the graph disappears and

\[
 h=Lp_0^2+Mp_0+4Nz^2p_0+Nzq_0.
\]

The specialized dual pair is an automorphism over K[z]; hence h is a
coordinate over K(z), by the explicit triangular shear used in HC4MYGV2.
At z=0 it has form (1), with a a specialized dual coordinate and thus
nonconstant. HC4MYGC1(a) applies. This includes arbitrary compositions
of base automorphisms and base-dependent dual automorphisms that preserve
the base/dual fibration. Their degrees, coefficients and word lengths
are unrestricted. The full nonlinear Hessian correction is retained:
constant determinant is always a hypothesis on the actual potential.

## 5. The standard chart eliminating p

Write the ambient coordinates as (x,y,p,q,r), u=1+xy. The universal identity

\[
 x^2 B+(2-3xy)A=u(u+1)p+2x(2u+1)q
\]

shows that on the hyperbola (x,y)=(z,-1/z), independently of the graph p,

\[
 A=-z^3r,\qquad
 \Psi=Lz^6r^2+(5Nz-Mz^3)r+\frac{2N}{z}q.          \tag{5}
\]

This is a coordinate in q,r over K(z), since q has nonzero field-constant
coefficient. The graph p=R(x,y,q,r) does not occur in (5).

**HC4MYGC3.** Every standard polynomial p-graph has an invertible gradient
whenever its Hessian determinant is a nonzero constant. The same is true
after an arbitrary base automorphism beta and an arbitrary relative affine
dual automorphism of the retained q,r variables.

Indeed use gamma(z)=beta^(-1)(z,-1/z). Its image is not contained in a
line: beta is an isomorphism and the closed hyperbola is isomorphic to
G_m, whose coordinate ring has a nonconstant unit, whereas a closed
line has coordinate ring K[T] with only constant units. An irreducible
closed curve contained in a line would be that whole line.

Relative affine dual changes preserve the nonzero quadratic dual Hessian
of (5). Thus f is not globally affine-dual. HC4MYGC1(b) proves the claim,
in every graph degree and every degree of the base and affine dual changes.

For arbitrary nonlinear dual automorphisms of q,r, (5) is still a
coordinate trace, but may become affine in the dual parameters. The same
argument proves the exact dichotomy: every such unit-Hessian potential
has invertible gradient unless it is globally affine-dual as in (2).
In that branch it is precisely a plane Keller cotangent lift. The theorem
does not exclude that entire remaining endpoint.

For example q=Q-L*x^7*V^2/(2N), r=V makes (5) affine in Q,V. With p=0,
however, the full potential has nonzero V^4 coefficient
9*L^3*x^16*u^4/(4*N^2). This control checks that an affine curve trace
alone does not establish a globally affine-dual potential.

## 6. Why the same graph-independent curve does not close the q chart

For the omitted coordinate q, the two linear-form coefficients are

\[
 A_q=3x(1+v)^2,\qquad B_q=y(9v^2+12v+1),\quad v=xy.
\]

Their polynomial ideal is exactly (x,y). Indeed write
C(v)=9v^2+12v+1, F(v)=-27v/2 and G(v)=(3v+2)/2. Then

\[
 F(v)(1+v)^2+G(v)C(v)=1,
 \qquad v=F(v)yA_q/3+G(v)xB_q.
\]

This puts v in (A_q,B_q), and
x=A_q/3-x(2+v)v and y=B_q-y(12+9v)v put x,y in that same ideal.
The reverse containment is immediate. The common-zero locus is only
the origin; it contains no nonconstant base curve.

Because L is nonzero, independence of Psi from an omitted q first
requires A_q=0 from its quadratic q coefficient, and then B_q=0.
Thus a graph-independent base-curve restriction, as used for r and p,
cannot directly close the q chart. This diagnoses the limit of that
method; it is not an exclusion of q graphs with constant determinant.
A further argument must retain actual graph dependence or change the
geometric mechanism.

## Remaining scope and reproduction

Base automorphisms alone no longer reopen the standard r-graph route.
The standard p-graph route is also closed, with the stated affine-change
extension and nonlinear-dual dichotomy. HC4MYGL2 now also closes the
q chart under all polynomial base and retained-dual automorphisms.
The later [HC4MYGB2](BASE_ELIMINATED_GRAPHS.md) closes both native
base-eliminated charts. Nonlinear reparameterizations of those charts, the
remaining p-chart affine-dual endpoint under nonlinear dual changes,
genuinely nonlinear transformations mixing base and dual variables, and uncovered
embeddings or mixed/coisotropic reductions remain open. General HC4
and JC2 remain unresolved.

All conclusions are stable under constant affine parameter changes:
the Hessian transforms by constant congruence and the gradient by
invertible linear source and output maps.

Point containment still does not automatically transfer the marked
v2 collision: the pullback covector is -grad(p)/2 and must agree at
the two parameter points, including when p itself is the graph function.
The theorems prove invertibility of every actual unit-Hessian member
in their covered classes, regardless of a particular marked collision.

```sh
timeout 60s .venv/bin/python research/jacobian-followup/hc4-triangular-graphs/verify_curve_trace.py
```

The [checker](verify_curve_trace.py) verifies curve substitutions,
normal forms, block determinants, representative pole identities and
inverse constructions. The rank exhaustion and unbounded local-pole
argument are written proofs, not finite searches. [Actual receipt](receipt-curve-trace.json).
The imported Tang and Rentschler results have the assurance recorded in
their canonical sources; no formal verification, independent external
review or novelty priority is claimed.
