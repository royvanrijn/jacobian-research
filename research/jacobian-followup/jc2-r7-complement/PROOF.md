# Proper nodal complement comparison and the closest quintic pencil

14 November 2026. PF2R7TOP1 supplies the missing proper comparison for
the known first image C on each of its two twelve-node opens. It gives
at most92 meridian-marked affine complement types across the complete
normalized gap-21 image family. PF2R7B5CONTACT1 gives an explicit contact
constraint if a second image has normalization degrees(3,5).

These are statements about the image curves and their complements,
not a classification of unions with unspecified additional images.
There is no finite-cover construction or JC2 conclusion. The geometry
is an authored written proof with named published inputs and exact
algebraic checks; no formalization or independent external review.

## 1. Uniform marked infinity data, including the exceptional parameters

Use the complete normalized coefficient charts L and K from PF2R7FAM1.
For L use its parameter d. For K one can use x nonzero, or the scaled
polynomial chart with z=x^3 nonzero. The first-image normalization is
finite birational, with monic degrees(6,10).

At infinity take T=p^(1/6), tau=1/T and write

    p=tau^-6, q=tau^-10*H(tau), H(0)=1.

The gap normalization gives no odd term of H below order21, and its
coefficient A at order21 is -rho11. In the projective chart about
[0:1:0] put X=p/q,Y=1/q. Then

    X=tau^4/H, Y=tau^10/H.

Choose the convergent uniformizing coordinate r=tau*H^(-1/4), so X=r^4.
A formal residue change of variables gives, for m>=10,

    [r^m]Y(r) = (6/(m-4))*[tau^(m-10)]H^((m-4)/4).       (1)

Indeed substitute r=tau*H^(-1/4) in the residue extracting [r^m]Y;
the derivative term is the derivative of H^((m-4)/4), and integration
by parts gives the displayed coefficient. This also works analytically
on a small parameter neighbourhood, since every leading unit is1.

The right side vanishes for odd m<31. At m31 it is (3/2)A. Consequently
the Puiseux characteristic is exactly(4;10,31) throughout both families.
In L the first odd projective coefficient is -1/13122. In the scaled K
chart it is z^6/19683. Both are units on their parameter domains.

The infinity semigroup is<4,10,41>, with conductor48 and delta24.
The infinity line Y0 has intersection multiplicity10. The union of
the curve branch and that marked line has delta34 and Milnor number67.
These marked local data are constant even at the affine exceptional
parameters. The correspondence between branch Puiseux data, intersections
and embedded type is recalled in
[Campillo--Greuel--Lossen, Section1, PDF page2](https://arxiv.org/pdf/math/0511407v2).
The numerical computations here are reconstructed by the checker.

The same expansion gives a uniform noncollision estimate. If two small
normalization parameters have the same X, they are r and zeta*r with
zeta^4=1. For zeta=+/-i the Y difference starts with2r^10. For zeta=-1
it starts with twice the nonzero order31 coefficient. After shrinking a
parameter neighbourhood, neither difference vanishes on its punctured
disk. All preimages of a sufficiently small infinity neighbourhood lie
in this disk, by compactness of the projective normalization family.
Thus no affine double point can escape into infinity in that neighbourhood.

## 2. The family is proper in the required ambient space

Homogenize the normalization map to degree10:

    [s:u] maps to [u^10*p(s/u):u^10*q(s/u):u^10].

There is no base point. At u nonzero the last coordinate is nonzero;
at u0 the middle coordinate is s^10. The unique source infinity maps
to [0:1:0], and there is no other intersection with the infinity line.
These facts are uniform over either parameter domain.

The projective normalization family is proper over the ambient product
P2 times the parameter base, so its image is a closed family of curves.
Work on the nodal opens

    UL = C minus the34 exceptional d values,
    UK = C minus {0 and the79 exceptional z values}.

Both bases are connected and path connected. Every fibre has precisely
twelve ordinary affine nodes by PF2R7COLL1 and PF2R7KCOLL1.

Near a fixed base point, each node continues as an analytic section:
the two distinct immersed normalization parameters solve the collision
equations transversely, so the implicit-function theorem applies.
All12 central nodes persist in disjoint neighbourhoods. Since every
nearby fibre also has exactly12 nodes, there are no additional node
sections. The uniform infinity estimate separately prevents escape.

Each node is locally a product of a transverse crossing and the base.
One may use local defining equations of its two smooth transverse
branches as the two ambient coordinates. Away from the nodes and infinity,
the normalization is a local immersion and the curve-family projection
is a submersion.

## 3. Direct Whitney check at infinity

The preceding construction gives locally finitely many smooth strata:
the ambient complement, the smooth curve part, the node sections, the
infinity line away from [0:1:0], and the infinity section. Projection to
the base is a submersion on every stratum. The node strata are products,
and the open ambient stratum satisfies the Whitney condition automatically.

For the remaining curve/section pair, write its parametrization as

    (r,b) maps to (r^4,Y(r,b),b),
    Y(r,b)=r^10+O(r^11).

The tangent planes on the punctured curve are spanned by

    (1,Y_r/(4r^3),0), (0,Y_b,1).

They tend to span(e_X,e_b), because Y_r/(4r^3)=O(r^6) and Y_b tends to0.
For a secant from (r^4,Y(r,b),b) to an arbitrary section point(0,0,b'),
the norm of its Y component divided by the total secant norm is at most
a constant times |r|^6. Hence every limiting secant lies in that same
limiting tangent plane. This proves Whitney(b), including varying b'.
The infinity-line/section pair is smooth and satisfies it directly.
There are no other frontier incidences at infinity.

The ambient projection P2 times U -> U is proper, since P2 is compact.
Thus all hypotheses of the first isotopy lemma have been checked for
this marked stratification. We import
[Mather, Notes on Topological Stability, Proposition11.1](https://webhomes.maths.ed.ac.uk/~v1ranick/surgery/mather.pdf),
PDF page27: a proper Whitney-stratified submersion admits local
trivializations preserving its strata. Its controlled-flow proof gives
the marked ambient comparison here.

It follows that (P2,C_b,L_infinity) is locally topologically trivial
on each nodal open. The affine complements are homeomorphic along paths,
with the curve-meridian conjugacy class preserved. One does not obtain
a global product over the punctured base or a fixed braid factorization
for a chosen projection; loops in the base may induce automorphisms.
No comparison of the joint union C with an unknown B or H is implied.

## 4. At most92 first-image complement types

Let omega^3=1. Exact polynomial identities give for L

    p_(omega*d)(t)=p_d(omega*t),
    q_(omega*d)(t)=omega^2*q_d(omega*t).

The same target transformation identifies K parameters x and omega*x.
It extends to a projective linear transformation preserving infinity.
Thus these coefficient points have meridian-marked homeomorphic complements.

Among L's34 exceptional points d0 is fixed and all other orbits have
size three, giving12 exceptional orbits. K has79 exceptional orbits.
Their three intersection points form one common orbit. There are
therefore90 exceptional coefficient orbits in the union. Section3 gives
at most two additional complement types from the two connected nodal opens.

This proves an upper bound of92 complement types, not92 distinct types.
It is a finite reduction for the first image alone. Replacing C by a
topologically equivalent representative need not preserve the algebraic
degree or position of an unspecified B under that homeomorphism.
Together Sections1--4 prove PF2R7TOP1.

## 5. A concrete closest-quintic construction

Let a possible additional irreducible target curve B have a birational
polynomial normalization of degrees(3,5) in the same affine coordinates
as C. This is an extra hypothesis: it is the smallest chart of PF2ATC1,
not a conclusion about every possible additional image.

Its implicit equation has ordinary degree5 and degree3 in V. Every
monomial U^i*V^j lies in the Newton triangle3i+5j<=15. To see the last
claim without an eliminant assumption, a larger weight among i+j<=5,
j<=3 has a unique exponent pair. Its pullback by a degree(3,5)
normalization could not cancel, a contradiction. The equation has a
nonzero V^3 coefficient and can be normalized to make it1.

There are13 monomials in this Newton triangle. For the known C, their
pole degrees are6i+10j. The only duplicated leading degree is30, from
U^5 and V^3. Starting with V^3-U^5, subtract the successive leading
monomials after restriction to(p(t),q(t)). The seven actual cancellations
have degrees28,24,22,18,16,12,10. The remaining polynomial has degree9.
Use the two lower-pole freedoms U and1 to set its t^6 and constant
coefficients to zero. This defines an explicit normalized quintic G_C.

The complete coefficients and the literal polynomial identity

    G_C(p(t),q(t))=R_C(t), deg R_C=9

are produced in [geometry.json](geometry.json). Its leading coefficient
is -1/6561 on L and2z^6/19683 on scaled K, both nonzero on their domains.
No coefficient parameter was divided out during the construction.

Replace the V^3 monomial in the13-dimensional Newton space by G_C.
The restrictions of this basis have distinct leading pole degrees

    30,28,26,24,22,20,18,16,12,10,9,6,0.

Consequently they are linearly independent at every parameter. An
implicit B equation with its V^3 coefficient1 has G_C coefficient1 in
this basis. Its restriction to C therefore has degree at least9, and
its exact possible pole degrees are contained in

    {9,10,12,16,18,20,22,24,26,28,30}.

Here the statement is necessary for actual degree(3,5) images; it does
not assert that every basis combination has that normalization.

## 6. Infinity contact and the next second-image gate

Bezout gives total intersection degree10*5=50. Both polynomial images
have their only infinity point at[0:1:0]. Normalization of C computes
the total finite intersection multiplicity as

    sum_(y affine) I_y(C,B)=deg_t B(p(t),q(t)).

Thus the finite sum is at least9, and

    I_infinity(C,B) in {20,22,24,26,28,30,32,34,38,40,41}.

The maximum41 is attained precisely when the implicit equation belongs,
up to a nonzero scalar, to the explicit pencil

    G_C(U,V)+alpha*U+beta.

Conversely every polynomial in this pencil has restriction degree9 and
contact41. The identity alone does not provide a polynomial normalization
of that quintic or its realization as a ramified Keller boundary. An
actual rational quintic would in addition need total affine delta4,
since its infinity cusp has type(2;5) and delta2. These are the next
concrete image tests. This is PF2R7B5CONTACT1.

The bound is on intersection multiplicity, not on distinct intersection
values. No transversality is inferred. Nor is degree(3,5) or maximal
contact derived for the actual unknown B. The other target charts,
the joint C/B/H complement relations, boundary attachment and source
realization remain part of the original research goal.

## 7. What the bounded braid probes establish

The exact witnesses d1 and z1 lie in their respective nodal opens.
Two Sage10.9/SIROCCO calculations reached their separate60-second caps
after constructing the degree10, six-strand implicit equations. Neither
returned a completed braid list. [bounds.json](bounds.json) retains the
log hashes and this terminal status; no live probe process remained.

They supply no generic fundamental-group presentation or S6 enumeration.
The old d0 braid sample remains a special triple-point member and cannot
be transported through the nodal opens by the theorem above. The proof
of proper comparison and the quintic identities are independent of the
capped braid calculations. The checker verifies exact algebra and source
custody; it does not formally reprove Whitney topology or Mather's theorem.
