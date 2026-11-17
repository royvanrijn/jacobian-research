# Cylindrical traces close every y=R(x,p,q) in all degrees

14 November 2026. In the scaled Meng--Yang potential `Psi=L*A^2+M*A+N*B`,
with `L*N != 0`, take the genuine graph

\[
                   y=R(x,p,q),\qquad \deg_{p,q}R>0.
\]

Then apply any actual polynomial automorphism of the three retained duals
over K[x], with polynomial inverse. **HC4MYHA2: a nonzero constant
four-variable Hessian determinant forces a polynomial gradient inverse.**
Both graph degrees and automorphism word lengths are unrestricted. This
extends the closed p-only y chart to every q-dependent graph independent
of the old r, including arbitrary mixed p,q dependence.

The proof combines two general facts: a cylindrical exceptional trace
with a noncoordinate generic value fibre forces gradient invertibility;
and the generic plane fibre of `P(F^2-1)+QF`, for nonconstant F, has
nonpositive Euler characteristic. No enumeration of critical branches,
finite-base escapes, Newton faces or parameter words is needed.

## 1. A nonsingular plane polynomial has generic Euler at most one

Let h in C[u,v] have nowhere-vanishing gradient. Its generic fibre is
geometrically irreducible. Here is the normalization argument, to specify
the precise fact used rather than assume rationality of a generic fibre.

Let L be the relative algebraic closure of C(h) in C(u,v), a finite
extension. The integral closure D of C[h] in L is finite over C[h].
It lies in C[u,v]: each of its elements is integral over that normal
polynomial ring and lies in its fraction field. Thus B=Spec D receives
a dominant morphism from A2. Restricting to a general affine line gives
a dominant rational parameterization of B; Luroth makes B rational.
Moreover D has no nonconstant unit, since D is a subring of C[u,v].
The smooth rational affine curve B therefore omits exactly one point
of its P1 completion. It follows that D=C[g] for a polynomial g in u,v,
and h=H(g) for a one-variable polynomial H. If deg H>1, H' has a zero
alpha. The nonconstant polynomial g attains alpha (restrict to a line
on which g is nonconstant), contradicting dh nowhere zero. Hence
deg H=1 and L=C(h).

In characteristic zero this relative algebraic closure criterion makes
the generic fibre geometrically integral, and all fibres in a nonempty
open set are smooth and irreducible. Such a fibre is an affine smooth
curve with at least one missing point on its projective completion, so

\[
                 \chi(h=\lambda)=2-2g-r\le1                 \tag{1}
\]

for generic lambda. This argument proves irreducibility, not rationality;
arbitrary genus is permitted. Only standard finite normalization, Luroth
and the geometric-integrality criterion are used here.

## 2. A general cylindrical-trace theorem (HC4MYHA0)

Let f(s,z1,z2,z3) be a four-variable polynomial. Suppose that for generic
s,lambda the three-dual fibre `f(s,-)=lambda` has Euler characteristic
different from one. Suppose also that at one scalar s0, a polynomial
automorphism of the three-dual fibre carries the trace f(s0,-) to a
plane polynomial h(u,v), independent of its third coordinate w.
The automorphism must have an actual polynomial inverse. Then

**nonzero constant Hessian determinant implies that grad f is a
polynomial automorphism.** This conclusion holds throughout every actual
retained-dual automorphism orbit over C[s].

Assume instead a unit-Hessian member with nonautomorphic gradient exists.
Its whole relative critical scheme C_f is smooth of pure dimension one
when nonempty. The trace critical scheme at s0 is isomorphic to
`Crit(h) times A1`. A positive-dimensional Crit(h) would give a critical
subset of dimension at least two, impossible. If Crit(h) has any point,
its product with A1 gives a literal closed affine line in C_f. Under unit
Hessian it must be a component. The retained affine-component theorem
[HC4MYGR1](CRITICAL_COMPONENTS.md), using the generic Euler hypothesis,
then forces gradient invertibility, again a contradiction. Consequently
Crit(h) is empty, C_f misses s0, and (1) gives exceptional generic-value
Euler at most one.

The same HC4MYGR1 also shows that **no** component of C_f anywhere can
be A1 in this hypothetical nonautomorphic case. Each smooth affine curve
of positive Euler is A1: its completion has genus zero and exactly one
omitted point. Therefore chi(C_f)<=0, without a dominance or finiteness
assumption on the base projection. This strengthens the geometric step
of [HC4MYGZ0](EXCEPTIONAL_FIBRE_CLOSURE.md).

Apply retained [HCTANG2](../tang-hc4-audit/STRUCTURAL_REVIEW.md). Rank-zero
projection of its constant jointly affine two-plane to s gives coordinate
generic dual fibres, contradicting the generic Euler hypothesis. In rank
one, constant affine changes preserving the fibration give

\[
 f=E(U,Z)+sA(U,Z)+VB(U,Z),\qquad J(A,B)\in C^\times.        \tag{2}
\]

As in HC4MYGZ0, the **whole** C_f is isomorphic to the smooth nonempty
plane curve C=B=0. The remaining equations solve for s,V as polynomials
on C. Let sigma denote this base function. On C, `dE=-sigma*dA`.
Since C_f misses s0, sigma-s0 is a unit. Since dA never vanishes on C,
the exceptional restricted value function E+s0*A is nonconstant on
every component: its differential is `(s0-sigma)*dA`. Its generic degree
d0 is positive. At exactly s=s0, stratifying (2) by B=0 gives

\[
 \chi(f(s_0,-)=\lambda)=1-\chi(C)+d_0\ge2,                 \tag{3}
\]

contradicting (1). These arguments handle both Tang ranks and all
possible whole critical components, including vertical ones and an
initially empty scheme. All actual retained-dual automorphisms preserve
the trace equivalence, the whole value fibres and the critical schemes
by [HC4MYGK0](RELATIVE_CRITICAL_LOCUS.md). No change of Hessian away
from its critical scheme is asserted.

## 3. Plane value fibres with no generic affine-line component (HC4MYHA1)

Let F in C[P,Q] be nonconstant and define

\[
 Z=PF+Q,\qquad G=P(F^2-1)+QF=FZ-P.                         \tag{4}
\]

Only finitely many fibres G=lambda have an actual component isomorphic
to A1. To prove this, parameterize such a closed component by t, and put
a(t)=F on it and b(t)=Z on it. The identities in (4) imply

\[
                 P=ab-\lambda,\qquad
                 Q=(1-a^2)b+\lambda a.                     \tag{5}
\]

If both a,b are nonconstant with degrees A,B>0, the coordinate degrees
of this **closed embedding** are A+B and2A+B. The larger lies strictly
between the smaller and twice the smaller, so neither divides the other.
This contradicts the Abhyankar--Moh degree theorem, explicitly imported
from [Garcia Barroso--Ploski, Theorem 6.9](https://arxiv.org/html/1208.0913v1#S6).
The [new import record](sources/degree-embedding-import.json) reuses the
already retained primary v1 PDF and specifies this additional theorem.

Thus at least one of a,b is constant. Both constant would make the
image a point. The remaining possibilities are exhausted as follows.

- If a=c is constant, the component is a straight line
  `(c^2-1)*P+c*Q=lambda`. Its coordinate ring C[a,b]=C[b]=C[t]
  forces b to have degree one.
- If b=c is nonzero constant, it is the parabola
  `Q=c-P*(P+lambda)/c`. Again a has degree one.
- If b=0, lambda is nonzero and the component is the vertical line
  `P=-lambda`, on which Z=0. The case lambda=0 would be a point.

Any such parabola intersects every straight line of these types. Two
parabolas from different values lambda also intersect: their graph
difference is a nonconstant polynomial in P, quadratic when their c's
differ and linear when their c's agree. Distinct value fibres are disjoint,
so an infinite collection of A1 components over distinct values must
consist of pairwise parallel straight lines.

For a fixed direction, the normal vector `(c^2-1,c)` permits at most
two constants c. Hence infinitely many of those lines would make one
of the fixed nonzero polynomials F-c or Z vanish on infinitely many
distinct lines. This is impossible. Here F-c is nonzero since F is
nonconstant, and Z is nonzero since deg Z=deg F+1. The b=0 vertical
lines are charged to Z; they are not accidentally folded into the
constant-a case. This proves the asserted finiteness.

A generic plane fibre is smooth. Every one of its components then has
Euler at most zero, since none is A1. Consequently

\[
                        \chi(G=\lambda)\le0                 \tag{6}
\]

for generic lambda. Generic irreducibility is not assumed or needed;
all components are counted. A special A1 component is permitted and
is not substituted for the generic fibre.

## 4. The whole two-dual-dependent y chart (HC4MYHA2)

Fix generic x=s nonzero. Since R(s,p,q) is independent of r, the
coordinate `T=A` has r coefficient -s^3 and is an actual polynomial
coordinate of this fixed three-dual fibre. Set

\[
 F=3+2sR,\qquad c_0=M-13N/(2s^2),
\]
\[
 P=\frac{Np}{4s^2}-\frac{9N^2}{16Ls^4},\quad
 Q=\frac{2Nq}{s}-\frac{3Nc_0}{4Ls^2},\quad
 T'=T+\frac{c_0+3NF/(2s^2)}{2L}.
\]

The entire fixed-base transformation `(p,q,r)->(P,Q,T')` is polynomial
with polynomial inverse: first recover p,q by affine formulas and then
r by its nonzero constant coefficient. The exact native potential is

\[
 f_s=L(T')^2+P(F^2-1)+QF-C_s,
 \qquad C_s=\frac{9N^2}{16Ls^4}+\frac{c_0^2}{4L}.             \tag{7}
\]

Here F is a nonconstant polynomial in P,Q, since deg_{p,q} R>0 for
generic s. Project a value fibre of (7) to the P,Q plane. Off G=lambda+C_s
it is a finite unramified double cover; over that plane curve T'=0
gives one sheet. Additivity of compactly supported Euler gives

\[
 \chi(f_s^{-1}(\lambda))=2-\chi(G=\lambda+C_s)\ge2.          \tag{8}
\]

The nontrivial double cover need not be a product. Its Euler characteristic
is twice that of the complement by the degree-two covering formula.

At x=0 the native trace is exactly the plane polynomial

\[
 h(p,q)=Lp^2+Mp+4NpR(0,p,q)^2+NqR(0,p,q),                  \tag{9}
\]

independent of r. Apply HC4MYHA0 using (8) and this cylindrical trace.
This proves gradient invertibility under unit Hessian for every such
graph and every actual retained-three-dual automorphism over K[x]. No
condition on R(0,p,q), its critical points or its degree is imposed.

## 5. Endpoint and remaining scope (HC4MYHA3)

Allow deg_{p,q} R=0 as well. The entire family y=R(x,p,q), after all
actual retained-three-dual automorphisms, has a nonautomorphic unit-Hessian
member exactly when the retained base-only unit-ideal plane Keller endpoint
exists. Positive (p,q) degree contributes none by HC4MYHA2, and the
base-only case is [HC4MYGU1/2/3](RELATIVE_COORDINATE_ENDPOINT.md).
No unit-ideal pair or plane-surjectivity proof is supplied.

Over any characteristic-zero coefficient field, embed the finitely many
coefficients of a hypothetical graph and its actual automorphism/inverse
into C. A polynomial gradient inverse descends by uniqueness. Thus the
invertibility statement has the same unrestricted characteristic-zero
scope as the earlier chart closures.

This closes every q-dependent y graph independent of the old r. The
y graphs with genuine r dependence, opposite x graphs with q/r dependence,
actual free-base/dual mixing and uncovered embeddings remain outside
the proved class. The unit-ideal plane-surjectivity endpoint and distinct
maximal-ideal nonlinear-p singleton-zero properness endpoint remain open,
as do unrestricted HC4 and JC2. Invertibility is not relabeled universal
unit-Hessian nonexistence.

[verify_cylindrical_trace.py](verify_cylindrical_trace.py) and its
[receipt](receipt-cylindrical-trace.json) check the exact transformations,
inverse maps, degree alternatives, line/parabola intersection equations,
Euler identities and retained source pins. The normalization, whole-curve
and unrestricted finiteness arguments are written proofs. Tang, the
retained affine-component theorem and the additional explicitly pinned
degree-embedding statement are inputs, not theorems proved by finite
controls. No formal verification, independent external review or novelty
priority is claimed.
