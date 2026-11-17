# The full source boundary lattice and the marked cable-cover test

15 November 2026. PF2RULED61. Programme archived; full goal active.

The valuation matrix retained in PF2RULED50 is the actual integral
intersection lattice of the whole smooth affine normalization. It
therefore determines the homology and linking form of its end, together
with the canonical class marking. This is an all-degree bridge from
source valuations to the target's branched cover at infinity. Passing
an abstract S4 knot representation is no longer the end comparison.

For the infinity knot of PF2RULED60's six-node control, all 72 peripheral
S4 representations fall into three conjugacy classes. The four-sheeted
branched-cover homologies are 0, 0 and Z/3, and the ramification knot is
zero in homology in all three. This excludes that target infinity knot
whenever the source ramified dual square is divisible by 3. It also
excludes every odd-order complementary source ruling, uniformly in the
unramified pole order and every principal-part coefficient. The common
pencil and certain even complementary orders remain open.

## 1. An integral lattice for every smooth affine plane extension

Let N be a smooth affine complex surface containing U=A2 with complement
the disjoint affine lines R_1,...,R_b. No field-degree bound is imposed.
Take a smooth projective completion X exposing these valuations and
resolving the whole boundary of U into an SNC divisor Gamma. Its
remaining infinity divisor is T=X minus N. Each closure of R_i is a
leaf after separating its one point at infinity. The total boundary
Gamma is a rational tree and its components are a Z-basis of Pic(X)
and H2(X,Z). This follows by dominating the standard completion of A2
by boundary point blowups; every blowup adds one exceptional basis
class. The resulting intersection lattice is unimodular.

Write checkR_i for the integral Gamma-supported dual divisor:

    checkR_i.Gamma_j = delta_(R_i,Gamma_j),
    Q_ij = checkR_i.checkR_j.                           (1)

The integer normalization of these valuations is essential. Further
boundary blowups pull back the old duals orthogonally to the new
exceptional class, so Q does not depend on the chosen refinement.

Thom excision for (N,U) gives H2(N,U)=Z^b and zero other positive
relative homology. Since U is contractible, H2(N)=Z^b and H1(N)=H3(N)=0.
The generators intersect the R_i with Kronecker delta. Equivalently,
in X their images are exactly the classes orthogonal to every component
of T, which have basis checkR_i. Thus Q is the WHOLE intersection form
of N, not just an auxiliary necessary pair matrix.

Affineness makes T connected and gives an ample divisor supported on
T (use an affine projective closure and a sufficiently small exceptional
perturbation after resolution). Hodge index therefore makes Q negative
definite on its orthogonal complement. This also agrees with Xie's
negative-definiteness criterion imported in PF2RULED50, but here the
integral topology identifies what that matrix measures.

A sufficiently large compact core W of N has a product end and boundary
M. Poincare--Lefschetz duality and H1(W)=H3(W)=0 give

    0 -> Z^b --Q--> Z^b -> H1(M,Z) -> 0,
    H1(M,Z)=coker Q,  |H1(M)|=|det Q|.                (2)

In particular the end is a rational homology sphere. Its linking form
is represented by -Q^(-1) modulo Z, with the outward boundary convention
(reversing all orientations reverses the sign). A proposed boundary
cover with positive b1 cannot be the end of N. No rationality or
definiteness claim about the separate SIGN surface is used here.

## 2. The canonical marking survives at infinity

Suppose the source volume extends with divisor sum_i k_i R_i, as in
an actual smooth finite Keller normalization where k_i=e_i-1. In the
dual homology basis the canonical cohomology class is the covector

    k=(k_1,...,k_b),  Q_ii = k_i modulo 2.             (3)

The parity is the characteristic property of the canonical class on
the smooth rational projective completion. Modulo Q, k restricts to
the Poincare dual of sum_i k_i K_i on M, where K_i is the oriented end
knot of R_i. This is the relative Thom class of the proper divisor,
not a guess about its normal framing.

For the full quartic, b=2 and k=(1,0), hence write

    Q = [ -n_R   c  ],  n_R odd, n_D even,
        [   c  -n_D ],  D=n_R*n_D-c^2>0.              (4)

The ramification knot class is the covector e_R modulo Q. Consequently

    [K_R]=0  iff  D divides n_D AND D divides c.       (5)

When (5) holds, division of the determinant identity gives

    n_R*(n_D/D) - D*(c/D)^2 = 1,
    gcd(n_R,D)=1.                                    (6)

These are unbounded arithmetic constraints on the actual source.
For example, if an end cover has order D=l^a, a>=1, and zero marked
ramification class, l cannot divide n_R. If its order is 1 and l is
3 modulo 4, l|n_R is also impossible: (4) would give -c^2=1 modulo l.
Here l is prime. These are necessary marked-end conditions, not
affineness or finite-map existence certificates.

## 3. Computing the branched cover without deleting residual sheets

For any knot presentation, sheet action and meridian word, take its
connected d-sheeted presentation complex. Collapse a spanning tree of
the lifted one-skeleton and lift every relator from every sheet. For
each cycle of length m of the meridian permutation, add the lifted
meridian^m filling relation once. Abelianizing gives a finite integer
matrix whose cokernel is H1 of the closed branched cover. This is just
cellular covering theory and Dehn filling; it works at arbitrary finite
degree and presentation length. It does not require asphericity of
the chosen presentation complex. For a simple ramified cycle, lift
the preferred longitude from either sheet in that cycle to retain K_R.
Adding that row computes the quotient by its homology class.

For transparency fix positive winding w and cable parameter p=w+1
on the positive trefoil. Its companion presentation and peripherals are

    <a,b | a^2=b^3>,  M=a*b^-1, L=a^2*M^-6.

The annular cabling decomposition adjoins a w-th root z of the cable
slope h=M^(w+1)*L^w. A meridian and longitude are

    z^w=h, beta=M^-1*L^-1, mu=z*beta,
    lambda=h*mu^(-w*(w+1)).                            (7)

The Bezout identity (w+1)*1+w*(-1)=1 fixes this meridian. To check the
positive convention directly, when the companion is an unknot, L=1
and (7) becomes z^w=M^(w+1), mu=z*M^-1, the usual positive torus group.
Weidmann's Corollary1 proof gives the annular root amalgam; the slope
and peripheral normalization are derived here. Ito's complete PDF5
fixes the regular-fibre and cabling conventions. Boundary section
generators from differently oriented tori must not be substituted
without converting them to the common M,L basis.

An equivalent five-generator presentation uses a,b,beta,h,mu with

    a^2=b^3, [beta,h]=[mu,h]=1,
    (mu^-1*beta)^w*h=1,
    M=beta^w*h, L=beta^(-w-1)*h^-1.                    (8)

The checker verifies the branched homology and marked class with BOTH
presentations. These are two implementations of the same derived
presentation, not an independent topological review.

## 4. Exhaustive peripheral data for the six-node control's knot

Set w=3, so p=4. PF2RULED60 forces rho(mu) to be a transposition and
rho(lambda)=1. Equation (7) then forces rho(h)=1. Enumerate every
a,b,beta in S4 and every transposition mu satisfying (8), retaining only
image S4. This is exhaustive for this fixed knot, not a degree cutoff
on arbitrary target curves. There are 72 representations, with three
simultaneous-conjugacy orbits of size 24. In zero-based permutation
tuples, representatives (a,b,beta,mu) can be chosen as follows:

    (3120,3102,2103,0213) : H1=0,    [K_R]=0;
    (3120,3102,2103,1023) : H1=0,    [K_R]=0;
    (3201,1032,3201,3120) : H1=Z/3,  [K_R]=0.          (9)

Each tuple lists the images of 0,1,2,3. Smith forms of the complete
lift-and-fill matrices certify these groups; adjoining the longitude
row leaves each group unchanged, certifying the zero marked class.
The third orbit must be retained: considering only the explicit S3
companion representation in PF60 would miss it.

Thus EVERY actual full quartic with target infinity knot C_(3,4)(T2,3)
has D in {1,3} and (5). If 3|n_R, D=1 contradicts -c^2=1 modulo3,
while D=3 contradicts (6). Therefore

    3 divides n_R => this entire target infinity class is EXCLUDED. (10)

This applies also to common source pencils when their computed n_R
is divisible by3, and is independent of all target coefficients.

For complementary source pencils, PF2RULED50 gives c=q*p_D and either
r=1,q=2,n_R=5, or r>=2,q=r-1,n_R=r^2-r+1. The unramified diagonal
n_D is even at every order and with every characteristic support.
If r=1 or r is odd, q is even, so D is even, contradicting (9).
If r=2 modulo3, 3|n_R, so (10) applies. In particular the ENTIRE
minimal complementary r=2,q=1 row is excluded for this target knot,
with arbitrary unramified zero/pole orders and all Laurent coefficients.
Any surviving complementary row for this target knot must have

    r>=4 and r=0 or4 modulo6.                         (11)

No polynomial or coefficient search proves (10)-(11); they follow
from the actual full source lattice and exhaustive finite monodromy.

As a strictness control, at w=8,p=9 the same peripheral conditions have
168 S4 representations in seven conjugacy classes. Three classes have
H1=Z + Z/26. They pass the peripheral equations and determinant-nine
test at infinity but fail (2) for EVERY smooth affine N of Section1.
The other four classes, with finite homology, are not excluded by this
argument. No polynomial realization of this numerical infinity type
is asserted here. Its eta computed from the trefoil recurrence is zero.

## What remains

The six-node curve remains an appropriate warning against excluding
eta0 from numerical normal degree and a knot representation alone.
The new source comparison excludes its knot for substantial unbounded
source sectors; it does not exclude every source valuation matrix.
The common-pencil cases with 3 not dividing n_R, and the even orders
in (11), still require a global obstruction. The source lattice must
also match the full linking form and canonical marking, beyond its
determinant. General non-binomial keys and higher finite degrees remain
in the original goal.

Exploratory identifications of individual boundary groups with named
three-manifolds, and potential symplectic filling classifications, are
NOT used: contact structures, orientation and all representation
classes would have to be checked before importing a filling theorem.
No full quartic exclusion, affine S4 extension for the control, Keller
map, independent full reproof, formal verification, external review,
novelty priority or goal completion is claimed.

[Exact checker](verify_marked_end_lattice.py) and
[source manifest](sources/marked-end-lattice/manifest.json).
