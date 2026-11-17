# Global sign twists and the exact two-point Keller gate

Theorems HC4MYJX0/1/2. Work over C. A unit Hessian means a nonzero
constant determinant in actual polynomial source coordinates. These are
all-degree statements; no source-degree or fibre-cardinality search is used.

## 1. Independent sign twists by an actual source automorphism (HC4MYJX0)

Let C(z) be nonconstant and squarefree, D(z) arbitrary, and

    h = 2xy + C(z)w + D(z).

For any chosen signs at the roots of C, choose p(z) by interpolation with
p=0 for sign +1 and p=1 for sign -1. Put a=1-2p. Then C divides
p(p-1). The polynomial matrix

    M(p) = [[1-p^2, 2p-p^3], [-p, 1-p^2]]
         = E12(p) E21(-p) E12(p)

has determinant one, M(0)=I, and M(1)=[[0,1],[-1,0]]. Write
(x',y')^t=M(p)(X,Y)^t. The polynomial

    q = 2x'y' - 2aXY

is divisible by p(p-1), hence by C. More explicitly,

    q = 2p(p-1)[(p+1)X^2 + (2p^2+2p-2)XY
                         + (p^3+p^2-2p-2)Y^2].

Set A=q/C. The map

    Psi(X,Y,z,W) = (x',y',z,W-A(X,Y,z))

is a polynomial automorphism with Jacobian one. Its actual inverse is:
first recover (X,Y)^t=M(p)^{-1}(x,y)^t, retain z, and then set
W=w+A(X,Y,z). Consequently

    h composed with Psi = 2a(z)XY + C(z)W + D(z).             (1)

In particular these are independent sign changes, not just simultaneous
scalings. They retain z but replace the embedded base C[x,z]. The unsquared
point Jacobians of the coefficient pair (2a(z)X,C(z)) are
2a(r)C'(r) at (X,z)=(0,r). Their squares agree with those for (2x,C).
Thus the unsquared orientation obstruction in
[JV2](TWO_SLICE_RECOVERY_AND_ORIENTED_COEFFICIENT_MODELS.md) is a
fixed-base obstruction, not an invariant of arbitrary source maps.

For C=z^2-1 choose p=(1-z)/2, so a=z. For every D this gives an explicit
z-preserving equivalence

    h_D = 2xy + (z^2-1)w + D(z)
       ~ C_D = 2zXY + (z^2-1)W + D(z).                      (2)

For D=z, C_D is exactly the earlier comparison model C2. This identifies
h1 with C2, not K2 with either of them. No nonlinear Hessian invariance is
being assumed: (1) is an equality of the actual pulled-back functions.

## 2. Every h_D encodes exactly a two-point Keller fibre (HC4MYJX1)

Let E2 mean: a polynomial map K:A2 -> A2 with nonzero constant Jacobian
has a target fibre consisting of exactly two points. This is the cardinality
of one actual fibre, not a claim about generic degree two.

For every polynomial D(z), the following existence statements are equivalent:

1. h_D has a unit-Hessian representative under a polynomial source automorphism.
2. Such a representative exists under a z-preserving source automorphism.
3. C_D has such a representative under a C[X,z]-affine automorphism of (Y,W).
4. E2 holds.

First suppose (1). The full critical equations of h_D give exactly

    x=y=0, z=r in {1,-1}, w=-D'(r)/(2r).

Both points are reduced Morse points, with Hessian determinant 16. Critical
scheme transport under source automorphisms preserves these two points.
A polynomial gradient automorphism cannot have two zeros. The retained
[Tang theorem](../tang-hc4-audit/STRUCTURAL_REVIEW.md) and
[HG0](GLOBAL_VALUE_ORBITS.md) therefore give actual coordinates (b1,b2,P,Q)
in which a hypothetical unit member is

    d(b) + P F(b) + Q G(b),   J(F,G) in C*.

Its full critical scheme is exactly the zero fibre of (F,G): for each
F=G=0 the invertible base Jacobian uniquely solves the two remaining
critical equations for P,Q. Thus this Keller fibre has exactly two points.
This proves (1) implies (4), for arbitrary source mixing.

Conversely suppose E2. Translate the target so that the chosen fibre is
zero. By an actual polynomial plane automorphism put its two source points
at (X,z)=(0,1),(0,-1). The finite-point interpolation construction and its
inverse are retained in JV2; its Jacobian is a nonzero constant. The resulting
Keller pair k=(k1,k2) has radical zero ideal

    I=(X,z^2-1).

Indeed the Keller condition makes every zero simple; Nullstellensatz then
identifies the ideal, not just its zero set. The pair

    f=(2zX,z^2-1)

also generates I: X=(z/2)f1-X f2. Its Jacobian is 4z^2, equal to 4 at both
zeros. If j=J(k1,k2), the two ratios J(k)/J(f) are the same unit j/4.
JV2 supplies M in GL2(C[X,z]) with k=Mf and det M=j/4.
For clarity, its elementary determinant repair starts with any k=M0 f.
Since u-det M0 lies in I=(k1,k2), write it as s k1-r k2. Adding
column(r,s) row(-f2,f1) to M0 fixes its product with f and changes its
determinant to u=j/4. This is a global polynomial matrix with a polynomial
inverse; no global Keller generator is inferred just from point data.
The globally Keller pair k is supplied by the E2 hypothesis itself.

Apply the fibre change (Y,W)^t=M^t(P,Q)^t. It fixes X,z, has inverse
(M^{-1})^t, and sends C_D to

    D(z) + P k1(X,z) + Q k2(X,z).

Its Hessian determinant is literally j^2, independently of D, by the
cotangent block determinant identity. This proves (4) implies (3).
Composing with the explicit map (2) proves (3) implies (2), and (2)
implies (1) immediately. Every inverse used is a polynomial inverse.

Hence the full-source unit-Hessian existence decisions for h0 and h1 are
identical, and are exactly E2. This does not identify h0 and h1 as functions:
their two critical values coincide for h0 and are distinct for h1. More
generally the proof does not change the critical values D(1),D(-1), and
does not construct E2 or resolve either existence decision.

There is a sharp source boundary. The older C[x,z]-preserving exclusion
for h_D remains valid; the sign twist changes x itself. The positive
conditional route fixes z and changes the hyperbolic pullbacks nonlinearly,
so it does not contradict the earlier x/y-preserving exclusions.

## 3. The universal affine test list contracts (HC4MYJX2)

Use the exact all-affine reduction in
[JQ0](TWO_NORMAL_MODELS_AND_RELATIVE_RANK_ONE_ALIGNMENT.md). Its universal
gradient-inversion property under actual unit Hessian was equivalent to
the four old tests (the coordinate endpoint, the nondegenerate quadratic
endpoint, h0, h1), the fixed model

    K2 = w(2xy+z^2-1)+z,

and the complete one-parameter three-point family K3_rho, equivalently its
actual suspension representative S_rho. Quantifiers remain over every
polynomial source automorphism, and over every complex rho.

By section 2 the h0 and h1 tests are the same requirement: absence of E2.
Any unit member of K2 has exactly two critical points (native critical
values +1,-1 and determinants 4). The same Tang argument gives E2. Thus
absence of E2 already excludes unit members of K2. The universal affine
property is therefore exactly equivalent to these four requirements:

- the coordinate endpoint has the universal gradient-inversion property;
- the nondegenerate quadratic endpoint has that property;
- E2 does not hold (equivalently, h0 has no unit source representative);
- S_rho has no unit source representative for every rho in C.

Necessity follows because each retained test is present in the old exact
list. Sufficiency follows from section 2 for h1 and from the two-critical-
point implication for K2, then JQ0. This is redundancy of K2 in a universal
test list, not an equivalence of its individual source orbit with h0, h1
or C2. If E2 holds, no K2 source representative is constructed here.
The K2 full-source orbit remains open and its previous one-coordinate
exclusion remains proved.

The two endpoint requirements are precisely the retained unit-ideal and
singleton-zero plane Keller endpoint questions in
[HG2](GLOBAL_VALUE_ORBITS.md). Do not replace the quadratic endpoint by
absence of unit representatives: the native quadratic has unit Hessian
and invertible gradient. A unit S_rho implies a three-point Keller fibre,
but the converse is not proved. The four requirements concern the admitted
affine graph/source family, not every HC4 embedding. Unrestricted HC4/JC2
and the original research goal remain open.

## Reproducibility and boundary

[The exact checker](verify_global_sign_twists.py) checks the universal
matrix factorization, divisibility and inverse identities; the explicit
four-coordinate twist and its inverse; critical schemes and determinants;
the generating-ideal identity, determinant repair, and symbolic cotangent
block formula. It binds the retained proof/checker/receipt bytes in
[the source context](sources/global-sign-twists-context.json).
[Its receipt](receipt-global-sign-twists.json) is a local symbolic replay,
not independent replay, formal verification or external review. The all-
degree arguments and the logical implications are proved above, not inferred
from symbolic samples. No new external theorem is imported in this step.
