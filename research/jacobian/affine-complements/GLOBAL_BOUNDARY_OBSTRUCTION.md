# Global boundary obstructions on the marked Wright surface

14 November 2026. PF2RULED32 applies at every pole height and every
finite field degree on the actual marked Wright model. It imports two
classical theorems of Nguyen Van Chau, with their hypotheses checked
below. It gives a global gate beyond the finite coefficient strata:

* Every polynomial coordinate on the target restricts to a nonconstant
  boundary polynomial with at least two distinct finite critical values.
* The boundary normalization must identify distinct parameter points.
* The finite-cover monodromy is S_d. Any meridional generating set for
  the image complement has at least d-1 elements.
* In degree four, the missing target values are exactly the boundary
  collision values. This set is finite and nonempty. Its cardinality
  equals the first Betti number of the boundary image.

These are necessary conditions. They do not exclude the remaining
colliding-image covers or supply a compatible completion of an arbitrary
normalization. The full goal remains active.

## The actual finite-map hypotheses

Keep the [marked ring and source open](WRIGHT_MARKED_MODEL.md)

    A=C[Y,XY,X^2Y,X^3Y+X] inside C[X,Y],
    N=Spec(A), U=N\R=A2_(X,Y), R=A1_v.

Suppose P,Q lie in A, J_XY(P,Q) is a nonzero constant, and
pi=(P,Q):N->A2 is finite of field degree d. The restriction F=pi|U
is a polynomial Keller endomorphism of the actual source plane.
These hypotheses include finiteness; membership of P,Q in A alone
does not establish it.

The [all-height image dictionary](../jc2-omitted-divisors/WRIGHT_IMAGE_RESTRICTION.md)
gives the finite birational normalization

    gamma:R=A1->C=pi(R), gamma(v)=(p0(v),q0(v)).

R is the sole ramification prime, with generic index two and residue
degree one. Outside C, pi is finite etale and its whole preimage lies
in U. It gives a d-sheeted covering there. At a point of C at least
one omitted local factor has length at least two, so the retained
source fibre has fewer than d points. Finite flatness follows from
the smooth source and target, or from the established fibre dictionary.
Consequently the entire exceptional/nonproper value set of F is
exactly C, not merely a curve contained in it. It is nonempty and
irreducible. Also d=3+deg gcd(p0',q0')>=3.

## The two primary literature inputs

Nguyen Van Chau, *A remark on Vitushkin's covering*, Acta Math.
Vietnam.24(1999),109-115, Theorem1.1: a polynomial Keller map which
is a covering outside a curve homeomorphic to C is bijective.
[Publisher PDF](https://math.ac.vn/uploads/files/9901109.pdf).

Nguyen Van Chau, *Two remarks on non-zero constant Jacobian polynomial
maps of C2*, Ann. Polon. Math.82.1(2003),39-44, Theorem2: a Keller
exceptional curve cannot be composed of polynomial images
(t^k,q(t)), with k positive. Corollary3 excludes a simply connected
exceptional curve.
[Publisher PDF](https://www.impan.pl/shop/en/publication/transaction/download/product/84523).

Both are imported results. The [source manifest](sources/global-boundary-topology/manifest.json)
pins the primary PDFs and exact theorem identifiers. Their proofs are
not independently re-proved or formally verified here. The deductions
below use these results only with their stated polynomial and global
covering hypotheses. No novelty priority is claimed.

## Every coordinate projection needs at least two critical values

Let L be a polynomial coordinate of the target plane, completed to
an automorphism (L,M). Its composition with F is again Keller and its
exceptional curve is (L,M)(C). If L(gamma(v)) were a shifted pure power
alpha*(v-beta)^k+c, rescale and translate the target coordinate and
translate the normalization parameter. The exceptional curve would
then be parametrized by (t^k,q(t)), contrary to Chau's Theorem2.
This works for every k, with no restriction on q, pole height, normal
corrections or field degree.

The restriction cannot be constant either. In that case C would be
the entire coordinate line L=c. The complementary coordinate M on
its finite birational normalization is affine linear; this is the
already excluded k=1 case.

For a nonconstant complex polynomial f of degree n>=2, having only
one finite critical value c forces f-c to have one distinct root.
Indeed, if its roots have multiplicities m1,...,mr, the derivative
has n-r zeros there counted with multiplicity. If these are all its
zeros, n-r=n-1, hence r=1. Thus f is a shifted pure power. A degree-one
polynomial is the excluded k=1 case as well. Therefore

    for every target coordinate L,
    L composed with gamma has at least two finite critical values.

In particular its degree is at least three and its derivative has at
least two distinct roots. This statement concerns target coordinates;
it does not assert the same condition for arbitrary elements of A or
arbitrary noncoordinate target polynomials.

## Image collisions are compulsory

If gamma were injective, its finiteness would make it a proper
continuous bijection from the complex line onto its image, and hence a homeomorphism.
The covering outside C supplied above would then satisfy Chau's
1999 Theorem1.1. F would be bijective, contradicting its nonempty
exceptional curve and d>=3. Thus gamma identifies distinct points.
This excludes every injective boundary image at every pole height,
without assuming a global pole-three root or solving a jet system.

For the quartic case the image dictionary gives

    #F^(-1)(y)=4-2*b(y)-1_(y=the cusp value),

where b(y) is the number of normalization parameters over y. The cusp
has b=1. Every collision has b=2 and both branches immersed; higher
branch counts or a cusp collision would give negative retained length.
Hence

    A2\F(A2) = { y : b(y)=2 }, a finite nonempty set K.

No nodality or transversality assumption is needed. In particular an
actual finite quartic Wright realization cannot be surjective.

Topologically C is the quotient of its normalization C by identifying
these finitely many pairs of points. One-point normalization fibres,
including the cusp, make no extra identifications. Identifying each
pair adds one independent circle, so

    b1(C)=|K|>=1, chi(C)=1-|K|<=0.

This does not prove the missing values cannot occur. Ruling out those
collisions from the actual finite algebra and source open is the
remaining global obstruction, rather than another local cusp test.

## Symmetric monodromy and a degree-free complexity bound

The complement of pi^(-1)(C) in N is connected, so the monodromy on
d sheets is transitive. A small meridian around a smooth point of C
acts as one transposition: R has e=2,f=1 and all other sheets there
are unramified. The complement group is normally generated by such
meridians. To see this, fill a loop by a disk in the simply connected
target plane and perturb the disk off the finite singular set and
transversely to the smooth part of C. Its finitely many intersections
express the loop as a product of conjugate meridians.

Thus the monodromy group is generated by transpositions. Form the
graph whose vertices are sheets and whose edges are these transpositions.
Its connected components are exactly the permutation orbits. Transitivity
makes the graph connected. Edge transpositions of a connected graph
generate all transpositions, using conjugation along a path, so the
monodromy is S_d. Moreover r meridians generating the complement group
would give a connected graph with at most r edges; hence r>=d-1.
This is a written proof for every d, not a subgroup census at small d.

## Why an unrestricted finite-jet campaign cannot replace the global gate

There is also a closed all-order control on this same ring. Write

    a=XY, B=X^2Y, D=X^3Y+X, z=B+1.

For every integer N>=0 and nonzero kappa define actual polynomials

    P=D,
    Q_N=-(kappa/2)*a^2*sum_(j=0)^N (j+1)*z^j.

In the omitted chart a=-u*(1-uv), z=uv, D=v. The finite geometric
sum gives the exact identity

    Q_N=-(kappa/2)*u^2*
        [1-(N+2)*(uv)^(N+1)+(N+1)*(uv)^(N+2)].

It follows for every N that

    J_uv(P,Q_N)=kappa*u modulo u^(N+2),
    J_XY(P,Q_N)=-kappa
       +(kappa/2)*(N+2)*(N+3)*z^(N+1)
       -(kappa/2)*(N+1)*(N+4)*z^(N+2).

The source Jacobian is nonconstant and vanishes on the retained line
X=0, where z=1. Both P and Q_N vanish on that entire line, so the map
even contracts a source curve. Nevertheless its boundary normal
identity passes arbitrarily many consecutive coefficients as N grows.
These are not Keller or finite-map examples.

There is no polynomial mate to D in A in any degree. The
[grading theorem](WRIGHT_GRADING_OBSTRUCTION.md) reduces any such mate
to Y*g(B). Its Jacobian with D is

    (3B+1)*g(B)+B*(B+1)*g'(B).

For nonzero g of degree e its leading coefficient is multiplied by
e+3, so it cannot be a nonzero constant. The formal limit of Q_N
near R is the rational mate -kappa/(2X^2); it has a pole on the actual
source line X=0. Formal boundary regularity does not remove that pole.

The generic field degree of (P,Q_N) is N+4: over C(v), Q_N is a
polynomial in u of that degree with nonzero leading coefficient. These
examples vary both field degree and pole height. They do not refute a
criterion that also uses fixed field degree, quasifiniteness or actual
finiteness.

Thus no uniform finite depth of boundary Jacobian tests, without a
pole bound and the global source conditions, suffices. The existing
fixed-pole completion cutoff remains valid: here the pole spaces of
Q_N grow, with the bound Q_N in V_(2N+6). A fixed-space proof and this
unrestricted all-order countercontrol have different hypotheses.

## Current research direction and verification boundary

Prioritize an obstruction to colliding exceptional images on the
actual finite Wright algebra, or a meridional-generator bound for the
forced image families. The coordinate-projection obstruction and the
identity between missing values and collisions provide global tests.
The stopped quadratic-M calculation is a
[historical checkpoint](sources/quadratic-remainder-checkpoint/STATUS.json),
with no quadratic-stratum exclusion or live process claimed.

Run the [scoped checker](verify_global_boundary.py):

```sh
.venv/bin/python research/jacobian/affine-complements/verify_global_boundary.py
.venv/bin/python research/jacobian/affine-complements/audit_pole_integration.py
```

The [receipt](global_boundary_receipt.json) verifies primary source
hashes, actual-ring identities, the symbolic all-N induction and
Jacobian formulas, local controls and finite graph controls. It does
not independently verify Chau's theorems or the written topological
proofs. No formal verification, independent whole-proof implementation
or external review is claimed. The colliding-image finite-map problem,
actual normalization-to-completion comparison and both quartic packets
remain open.
