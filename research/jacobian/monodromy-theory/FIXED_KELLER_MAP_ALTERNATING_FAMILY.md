# One fixed Keller map for a polynomial pencil

This note realizes all separable fibres of a ramified reparametrization of a
polynomial pencil inside **one** polynomial Jacobian-one map of affine
three-space. For the rational degree-ten alternating pencil this gives
infinitely many full A10 arithmetic fibres of one fixed map. The map's generic
arithmetic and geometric groups are S10. No polynomial realization of a
generic A10 extension is asserted.

The construction uses the retained [quadratic-gauge compiler](../../../unpublished/common-arithmetic-fibers/sections/02-quadratic-gauge.tex)
and its [whole-fibre reconstruction](../../../unpublished/common-arithmetic-fibers/sections/03-realization.tex).
The independent quadratic target coefficient is retained, so the
[universal quadratic-tilt theorem](../universal-keller-structure/PROOF.md)
continues to determine the generic symmetric group.

## Theorem 1. Uniform fixed-map realization

Let k have characteristic zero, H(T) in k[T] have degree N at least three,
and a belong to k. Choose tau in k such that c3 is nonzero, where

\[
 H(\tau+T)=H(\tau)+AT+\sum_{j=2}^N c_jT^j,
 \qquad \kappa=H(\tau)+a-\tau A.
\]

There is a single polynomial map F:A3_k -> A3_k with determinant one,
generic degree N, and a rational target curve Gamma:Gm -> A3 such that,
for every p in k* with H(T)-(A-p^5)T+a separable,

\[
 F^{-1}(\Gamma(p))\simeq
 \operatorname{Spec}k[T]/(H(T)-(A-p^5)T+a).
\]

This is an isomorphism of the entire scheme-theoretic fibre. It commutes
with field extension. If N is at least four, the coordinate degrees of F
are at most 20N-54. Its generic arithmetic and geometric groups are S_N.

Proof. Define fixed coefficient profiles

\[
 g_1=1,\quad g_2=0,\quad g_3=c_3,
 \qquad g_j(V)=c_j V^{2j-8}\quad(4\le j\le N).
\]

All exponents are nonnegative. In source coordinates (x,y,z), set

\[
 t=1+xy,\qquad q=t^2z+c_3^{-1}y^2(1+3t),\qquad V=tq,
\]
\[
 B=y+3c_3xq+
       \sum_{j=4}^N j c_j V^{2j-8}t^2x^{j-2}q^j,
\]
\[
 C=x(5-3t)-c_3x^3z-
       \sum_{j=4}^N(j-2)c_jV^{2j-8}(xq)^j,
 \qquad F=(V,-B/2,C).
\]

These are fixed polynomials. In the source chart t nonzero, put
S=x/t and Q=y+xq. Then

\[
 D=1-SQ+VS^2=t^{-1},\qquad
 G_V(S)=S+c_3VS^3+\sum_{j=4}^N c_jV^{3j-8}S^j.
\]

The formulas become

\[
 B=Q+\frac{G_V'(S)-1-VS^2}{S},\qquad
 C=2G_V(S)-BS^2.
\]

The displayed quotient is polynomial in S. At fixed V, the plane
Jacobian of (B,C) in (S,Q) is -2D. The source chart Jacobian of
(V,S,Q) in (x,y,z) is t. Consequently det JF=1 on this dense chart,
and therefore everywhere as a polynomial identity. Derivatives of the
profiles in V do not affect this determinant: the first target coordinate
is V itself.

For a target (v,Y,Z), the inverse polynomial in S is

\[
 E_{v,Y,Z}(S)=G_v(S)+YS^2-Z/2.
\]

Take

\[
 \Gamma(p)=\left(p,\frac{c_2}{p^2},
                  -\frac{2(\kappa+\tau p^5)}{p^8}\right).
\]

The decisive identity, valid as a polynomial identity after clearing p,
is

\[
 \boxed{p^8 E_{\Gamma(p)}(S)
       =H(\tau+p^3S)-(A-p^5)(\tau+p^3S)+a.}
\]

Its coefficients in degrees 0,1,2,3 are respectively
kappa+tau p^5, p^8, c2 p^6, c3 p^9; in every degree j at least four
both sides have coefficient cj p^(3j). Thus the identity is uniform
in N, rather than a finite interpolation check.

To verify fullness, on V=p the coefficient profiles equal the constants
gj(p): their differences belong to the ideal (V-p). Thus this is exactly
the constant-profile compiler fibre. Since p is nonzero, both t and q
are units on this fibre; no points lie outside the chart. There
Q=B-(G_p'-1-pS^2)/S, B=-2Y, and

\[
 D=G_p'(S)-BS=E_{\Gamma(p)}'(S).
\]

Separability makes D a unit in k[S]/E. The inverse formulas are

\[
 t=D^{-1},\quad x=S/D,\quad y=Q-pS,\quad q=pD,
 \quad z=D^2\{pD-c_3^{-1}y^2(1+3/D)\}.
\]

They reconstruct all source equations. Hence the quotient by the three
target equations is the whole algebra k[S]/E, with T=tau+p^3S.
This proves the asserted scheme isomorphism.

Over k(V), G_V has degree N and linear coefficient one. The independent
target coordinates Y,Z give the extension defined by
G_V(S)+YS^2-Z/2. The retained quadratic-tilt theorem gives arithmetic
and geometric S_N, and in particular degree N. Finally deg(t)=2,
deg(q)<=5 and deg(V)<=7. The j-th high terms in B and C have degrees
at most 20j-54 and 20j-56 respectively; the low terms have degree at
most seven. For N>=4 this proves the stated bound. The choice of tau
always exists because H''' is a nonzero polynomial in characteristic
zero and k is infinite.

## Proposition 2. Solvable reparametrization preserves perfect regular groups

Suppose the splitting field M/k(s) of H(T)-sT+a is regular with perfect
Galois group G. Then the pullback s=A-p^5 has splitting-field group G
over k(p), and is regular over k.

Proof. The normal closure of k(p)/k(s), with p^5=A-s, has solvable
Galois group: it is a subgroup of C5 semidirect Aut(C5), after adjoining
fifth roots of unity if necessary. Its intersection with M is a Galois
extension of k(s) whose group is both a quotient of G and solvable.
A perfect group has no nontrivial solvable quotient, so the intersection
is trivial. This proves linear disjointness and the asserted group.
Apply the same argument over an algebraic closure of k; it gives the
same geometric group G. Equality of arithmetic and geometric groups
proves regularity. This argument uses perfectness, not the unsupported
assertion that every radical base change preserves a Galois group.

## Corollary 3. One rational Keller map with infinitely many A10 fibres

Use the [degree-ten pencil](../../arithmetic/RATIONAL_DEGREE_TEN_ALTERNATING_PENCIL.md)

\[
 H(T)=9T^{10}-864T^8+28161T^6-350150T^4+1127500T^2,
 \qquad a=250000.
\]

With tau=1 the coefficients needed in Theorem 1 are

\[
 A=1016544,\quad\kappa=38112,\quad c_2=-574772,
 \quad c_3=-884684,
\]
\[
 (c_4,\ldots,c_{10})=(13675,122850,5859,-5832,-459,90,9).
\]

Thus the single fixed map in Theorem 1, of coordinate degrees at most
146, has complete fibres along

\[
 \Gamma(p)=(p,-574772/p^2,-2(38112+p^5)/p^8)
\]

isomorphic to the degree-ten root algebras at s=1016544-p^5.
There are infinitely many rational p for which those algebras are fields
with splitting-field group A10. The fixed map itself has generic group
S10.

Proof. The retained degree-ten theorem proves regular A10 for the
pencil. A10 is perfect. Proposition 2 proves regular A10 after the
displayed reparametrization. Hilbert irreducibility, in the exact
formulation already imported in the
[specialization note](../../arithmetic/ALTERNATING_PENCIL_SPECIALIZATIONS_AND_KELLER_FIBRES.md),
gives infinitely many full-group rational specializations after deleting
the finite branch set and p=0. Transitivity of the natural A10 action
gives fields of degree ten. Theorem 1 identifies their entire fibres.

## Evidence and remaining scope

The [exact checker](verify_fixed_keller_alternating_family.py) checks the
full two-variable polynomial identity, the symbolic chart Jacobians and
polynomial pullbacks, and whole-fibre reconstruction in a degree-ten
quotient. Its [receipt](fixed_keller_alternating_family_verification.json)
records internal replay. The quantified identity, scheme proof and
solvable-base-change proof are written above. No independent replay,
formal verification, or external review is claimed.

This settles the requested fixed-map realization for the degree-ten
pencil. It does not identify the generic S10 map with a generic A10
cover, and it does not supply one fixed map for arbitrary rational
coefficient families such as the higher-rank Hasse resolvents.
