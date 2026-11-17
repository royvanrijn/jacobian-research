# Whole value fibres and arbitrary source reparameterizations

14 November 2026. Affine trace criteria require actual affine source
hyperplanes. This note instead uses the full critical scheme and whole
three-dimensional value fibres, which survive **every** polynomial
source automorphism. It extends the zero-boundary graph obstruction
and both existing endpoint identifications to arbitrary mixing of the
four source variables.

In particular, full source mixing on a base-only graph does not create
a new endpoint: a graph through the base origin has no unit-Hessian
reparameterization, while an origin-avoiding graph gives exactly the
existing unit-ideal plane Keller problem. The fixed degree-six p graph
continues to give exactly the singleton-zero problem under all source
automorphisms. Neither endpoint is solved.

## 1. Global invariants and the whole value profile (HC4MYHG0)

Let h be a polynomial in four variables and Phi an actual polynomial
automorphism, with constant Jacobian determinant delta nonzero. Put
g=h composed with Phi. With I_h the full gradient ideal,

\[
 I_g=\Phi^*I_h,\qquad
 (\operatorname{Hess}g)|_{\operatorname{Crit}(g)}
  =D\Phi^T\,\Phi^*(\operatorname{Hess}h)\,D\Phi.             \tag{1}
\]

The second identity is on the full critical scheme: every extra chain
rule term has a coefficient h_i and vanishes there. Consequently

\[
 \det\operatorname{Hess}g
       \equiv\delta^2\Phi^*(\det\operatorname{Hess}h)
                                          \pmod{I_g}.       \tag{2}
\]

Also Phi identifies h=lambda and g=lambda for every lambda, including
their critical points and values. These are full source invariants;
there is no assertion of nonlinear Hessian invariance off Crit.

If g has unit Hessian, Crit(h) must be finite and reduced, possibly
empty, and every point must be Morse. Moreover the native Hessian
determinants at all those points must have one common nonzero value.
A positive-dimensional critical scheme, a degenerate critical point,
or two different critical Hessian values therefore forbids **every**
unit-Hessian member of this source-automorphism orbit. Empty critical
schemes give no critical-determinant constraint.

Now assume g has unit Hessian and a nonautomorphic gradient, over C.
Retained [HCTANG2](../tang-hc4-audit/STRUCTURAL_REVIEW.md) puts it, by
a constant linear source change, into

\[
 g=D(z_1,z_2)+P F(z_1,z_2)+Q G(z_1,z_2),
                              \qquad J(F,G)\in C^\times.   \tag{3}
\]

Let Z=V(F,G), n=#Z, and n_lambda=#\{z in Z:D(z)=lambda\}.
The full critical scheme is isomorphic to Z: its last equations solve

\[
                (P,Q)^T=-J(F,G)^{-T}\nabla D.               \tag{4}
\]

In particular all its points are reduced and their critical values
are D(z). Project the entire value fibre of (3) to the z plane.
On F nonzero solve P, with Q free. On F=0,G nonzero solve Q, with
P free. These are disjoint actual product strata, together giving
E=(t^2-n)t. At each z in Z with D(z)=lambda, the entire (P,Q) plane
is free. Thus for **every** scalar lambda,

\[
                  E(g=\lambda)=t^3-n t+n_\lambda t^2,
                                     \qquad t=uv.           \tag{5}
\]

The generic profile is t^3-n*t; the finitely many value jumps count
exactly the critical points at that value, and sum n_lambda=n. This
uses only the retained basic Hodge invariant, additivity and literal
product strata. It uses no arbitrary finite-cover multiplicativity.

Therefore, if the critical scheme and whole value-fibre polynomials
of h violate (5), a unit-Hessian reparameterization can only have an
automorphic gradient. If Crit(h) does not consist of one point, that
remaining alternative is impossible too, because a gradient
automorphism has exactly one zero. Formula (5) is asserted under
hypothetical nonautomorphy; it is not imposed here on every polynomial
with an automorphic gradient.

### The full linear-perturbation family measures the missing length

For (3), replacing g by g-alpha*P-beta*Q leaves its full Hessian
unchanged and replaces (F,G) by (F-alpha,G-beta). Write n(alpha,beta)
for that plane fibre's cardinality, and d for the generic degree of
the plane Keller map. Its generic value polynomial is

\[
                         t^3-n(\alpha,\beta)t.              \tag{6}
\]

By the retained [HC4MYOR2 finite-normalization argument](../hc4-relative-row-orbit/PROOF.md#6-properness-and-the-exact-missing-length),
the normalization N of the target in the source function field is
finite flat of degree d and contains the source as an open subset.
Its target fibre has total length d; each retained source point has
length one. Hence

\[
          d-n(\alpha,\beta)=
             \text{total omitted normalization length}.   \tag{7}
\]

The nonproperness locus is exactly where this quantity is positive:
the finite image of the closed omitted subset is precisely where
source properness fails. Equations (6),(7) identify the Hodge-coefficient
variation with that defect. They do not show that it vanishes.
Checking only generic perturbations misses the exceptional target set.
No topological local-triviality or absence of bifurcations at infinity
is inferred from the E-polynomial alone.

## 2. Two all-degree nonexistence consequences (HC4MYHG1)

### Zero-boundary base graphs after all source automorphisms

Use Psi=L A^2+M A+N B with L*N nonzero and H(p)=Lp^2+Mp. Consider

\[
 y=R(x,p,q,r),\quad R(0,p,q,r)=0,                           \tag{8}
\]

or the opposite x=R(y,p,q,r) with R(0,p,q,r)=0. Then **no polynomial
automorphism of all four source variables** can give unit full Hessian.
It need not preserve either physical base or any affine slice.

For (8) write R=xT and p_star=-M/(2L). At x=0 the trace is H(p).
Its three relative derivatives vanish on p=p_star. The remaining
full derivative, on that plane, is exactly

\[
              f_x=N\big(q T(0,p_\star,q,r)+2r\big).          \tag{9}
\]

It vanishes at q=r=0. Its r derivative there is 2N, so the full
critical scheme contains a smooth curve through that point, lying
in x=0,p=p_star. In particular this is a full critical point with
singular full Hessian. Equation (1), or simply critical-scheme
dimension invariance, forbids unit Hessian after any source automorphism.

For the opposite graph write x=yT. The corresponding remaining
derivative on y=0,p=p_star is

\[
              f_y=N\big(q+2r T(0,p_\star,q,r)\big).          \tag{10}
\]

Its q derivative at q=r=0 is N. The same conclusion follows.
All boundary jets, graph degrees and source automorphism degrees
are unrestricted. This extends HC4MYGU0's graph application beyond
automorphisms preserving the free-base fibration. The more general
composite-trace lemma alone does not extend that way: a normal
derivative might be nowhere zero on its relative critical set.

### A split one-variable potential cannot acquire unit Hessian in high degree

Let h=Q_3(z_1,z_2,z_3)+b(w), with Q_3 a nondegenerate quadratic form.
If deg b>=3, **no** polynomial source automorphism makes its Hessian
a nonzero constant. This needs neither Tang nor Hodge theory.

At a critical root rho of b, the full Hessian determinant is
det(Hess Q_3)*b''(rho). A multiple root of b' gives a degenerate critical
point and is forbidden by (1). Otherwise b' has deg b-1 distinct roots.
Unit Hessian after a source automorphism would force b'' to have one
common nonzero value c at all of them. But b''-c has degree at most
deg b-2, smaller than their number, so it vanishes identically.
This forces deg b<=2, a contradiction in every degree.

For constant b there is a critical line, also forbidding unit Hessian.
For nonconstant linear b, h is a polynomial coordinate. For quadratic
b of degree two, completing the square makes h a nondegenerate
four-variable quadratic plus a constant. Their two remaining full
source orbits are exactly the endpoints below, not new searches in
the degree of a reparameterization.

## 3. The two known endpoints are full source orbits (HC4MYHG2)

### The coordinate orbit is exactly the unit-ideal endpoint

Existence of a unit-Hessian polynomial **coordinate** on C^4 is
equivalent to existence of a unit-ideal plane Keller pair. A polynomial
coordinate means one component of an actual polynomial automorphism.
This statement imposes no preserved free base.

A coordinate has no full critical point, so its gradient omits zero
and is not an automorphism. Tang gives (3), and (4) makes Z empty.
Thus (F,G)=(1). Conversely every unit-ideal plane Keller pair gives
the explicit relative coordinate (3) and inverse completion in
[HC4MYGU2](RELATIVE_COORDINATE_ENDPOINT.md#3-classify-unit-hessian-relative-coordinates),
for arbitrary D. It is also an absolute coordinate on C^4.

By [HC4MYGU1/3](RELATIVE_COORDINATE_ENDPOINT.md), a base-only graph
x=h(y) or y=h(x) with h(0) nonzero has native potential a relative,
hence absolute, coordinate. Its orbit under **all** four-variable
source automorphisms is exactly the orbit of one coordinate, namely
all polynomial coordinates on C^4. It therefore gives exactly the
unit-ideal endpoint even after arbitrary base/dual mixing. If h(0)=0,
HC4MYHG1 instead forbids unit Hessian throughout the full orbit.

As before, absence of unit-ideal plane Keller pairs is equivalent to
plane Keller surjectivity. It is not identified here with full JC2.

### The nondegenerate quadratic orbit is exactly the singleton-zero endpoint

Let Q_4=z_1 P+z_2 Q. Universal gradient invertibility among unit-Hessian
members of its full polynomial source-automorphism orbit is equivalent
to invertibility of every singleton-zero plane Keller pair.

Indeed Q_4 has one reduced full critical point. If Q_4 composed with
Phi has unit Hessian and nonautomorphic gradient, Tang gives (3),
and (4) makes Z a singleton. Its ideal is the maximal ideal of that
point, giving the existing singleton-zero problem.

For the converse, translate a singleton fibre point to z=0. Given a
Keller row k=(F,G) with (F,G)=(z_1,z_2), the retained constructive
[HC4MYOR1 row theorem](../hc4-relative-row-orbit/PROOF.md#2-exact-classification-of-the-whole-orbit)
provides an invertible polynomial matrix C(z) such that

\[
                      (z_1,z_2)C=(F,G).                    \tag{11}
\]

One may see the harmless determinant normalization directly. Start
with any A satisfying z*A=k and set c=J(k)(0) nonzero. Since det A-c
is in (F,G), write it mu*F+nu*G. With k_z=(-z_2,z_1)^T, the rank-one
correction C=A+k_z*(nu,-mu) fixes z*A and has det C=c. Thus
(z,w) to (z,C(z)w) is an actual polynomial source automorphism and

\[
                Q_4(z,Cw)=P F(z)+Q G(z).                  \tag{12}
\]

This has unit Hessian J(F,G)^2, and its gradient is invertible exactly
when (F,G) is, by the cotangent factorization. No plane inverse is
assumed in constructing C.

More generally any hypothetical nonautomorphic unit-Hessian polynomial
with exactly one critical point is, after subtracting its critical
value, right-equivalent to Q_4: in (3), D-D(z_0) is in (F,G), so a
polynomial translation of P,Q absorbs it before (11). The zero-critical-
point case is a coordinate by the preceding paragraph. This classifies
these two possible critical-count strata under nonautomorphy, not all
HC4 polynomials.

The fixed degree-six p graph of [HC4MYGP1](P_GRAPH_SINGLETON_ENDPOINT.md)
has native potential carried by its actual polynomial frame and shear
to (Mx+Ny)T+2NxW. A constant invertible base change makes this Q_4.
Its **entire** four-variable source-automorphism orbit is therefore
the quadratic orbit just analyzed. Arbitrary base/dual mixing of that
fixed graph retains exactly the same singleton-zero endpoint; it does
not open a further graph-degree or automorphism-word frontier.

For reference, coordinate fibres have E=t^3 for every value. Q_4 has
E=t^3-t at nonzero values and t^3+t^2-t at zero. These are exactly
the n=0 and n=1 global profiles in (5). Neither profile supplies the
missing properness theorem.

## Scope and controls

The critical-scheme and determinant obstructions, including HC4MYHG1,
work over any characteristic-zero field by extension to its algebraic
closure. Hodge and endpoint equivalences are stated over C; inverse
conclusions for a stipulated characteristic-zero candidate descend by
uniqueness after embedding its finitely generated coefficient field.

[verify_global_value.py](verify_global_value.py) and its
[receipt](receipt-global-value.json) check the full chain rule and
critical determinant transport, global product-stratum formulas,
zero-boundary normal jets, split-potential obstructions and endpoint
row/inverse controls. Written arbitrary-degree arguments remain
separate from finite identities. No arbitrary finite set is claimed
to be realized as a plane Keller fibre.

Other dual-dependent base graphs and mixed embeddings remain open
where these global invariants and earlier affine-trace criteria do
not exclude them. The two distinct plane endpoints, full HC4 and
JC2 remain unresolved. No formal verification, independent replay,
external review or novelty priority is claimed.
