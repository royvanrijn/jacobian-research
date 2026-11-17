# Matched meridian slopes close the entire common two-stage quartic class

15 November 2026. PF2RULED76 closes the residual r=3 family in PF75.
Consequently EVERY COMMON source in the actual primitive-binomial
TWO-STAGE QUARTIC class is excluded, at all first exponents, windings,
source orders and poles. With PF50's pencil dichotomy, any remaining
realization of that class must be in the COMPLEMENTARY source pencil.
That pencil and the full quartic/higher-degree goal remain OPEN.

The new mechanism has no degree bound: matched slopes outside the
INTERIORS of both actual pieces' L-space intervals obstruct gluing
to PF74's common-coordinate source end. A uniform Seifert endpoint
lemma below avoids a parameter sweep or full interval computation.
Published topology is imported with complete primary pages audited.
No independent reproof of those deep results is asserted.

## 1. A general Seifert endpoint lemma

Let W be the negative-definite Seifert plumbing with orientable base S2,
normalized data

    e0; (alpha_i,omega_i), i=1,...,nu,
    nu>=2, 0<omega_i<alpha_i, gcd(alpha_i,omega_i)=1,
    epsilon=e0+sum omega_i/alpha_i<0.

Let Y be the complement of an ordinary REGULAR fibre in W and suppose
Y is not a solid torus. Give its boundary Rasmussen's reverse-oriented
basis (fibre,-meridian); the ACTUAL meridian filling recovering W is
slope0. This is not the fibre slope infinity. Assume

    e0+nu>=1.                                              (1)

Then 0 is NOT in the INTERIOR L^o(Y) of its L-space interval.
It can be an L-space endpoint; failure to distinguish the interval
from its interior would invalidate the argument.

Rasmussen Proposition3.9, equation(27), gives

    y_+=min_(l>=1) (1-l*e0-sum ceil(l*omega_i/alpha_i))/l,
    L(Y)=[[y_-,y_+]].                                      (2)

First consider nu=2. Since epsilon<0 and e0+2>=1, e0=-1 and the
two normalized fractions x,y have x+y<1. For every positive integer l,

    ceil(lx)+ceil(ly)<=ceil(l(x+y))+1<=l+1.

Thus every numerator in (2) is nonnegative, while l=1 gives0.
Hence y_+=0, which is an endpoint, or the omitted point if y_-=y_+.
In either case 0 is outside the INTERIOR. No extension of Nemethi's
nu>=3 Pinkham criterion to a two-leg graph is used.

Now let nu>=3. If W's graph is non-rational, Nemethi Theorem1.3 says
W is not an L-space. Therefore 0 is not even in L(Y). If W is rational,
the negation of Nemethi equation(1.13) gives, for EVERY l>=1,

    sum floor(-l*omega_i/alpha_i)>l*e0-2,
    l*e0+sum ceil(l*omega_i/alpha_i)<=1.                    (3)

The integer inequality shows y_+>=0. At l=1 its numerator is
1-e0-nu<=0 by (1), so equality is forced: y_+=0. Again slope0
cannot be in L^o(Y). This proves the endpoint lemma without classifying
the rational cases or searching growing denominator intervals.

## 2. The degree-free matched-slope obstruction

Let Y_1,Y_2 be non-solid-torus graph manifolds with one torus boundary,
glued by an orientation-reversing actual map phi. If a slope gamma is
outside L^o(Y_1) and phi(gamma) is outside L^o(Y_2), Rasmussen
Corollary6.1 gives

    phi(L^o(Y_1)) union L^o(Y_2) != all slopes,
    M=Y_1 union_phi Y_2 is NOT an L-space.                  (4)

For a closed orientable IRREDUCIBLE graph manifold, the published
graph-manifold equivalence (Nemethi Theorem1.5) then makes pi1(M)
nontrivial and LEFT-ORDERABLE. This cannot be the end group of a source
satisfying PF74's common-coordinate hypotheses: every homomorphism
from that source end group to any left-orderable group is trivial.

This obstruction has no field-degree, cable-length or parameter bound.
It requires actual pieces, actual gluing, the two non-solid-torus
hypotheses, a genuinely matched slope, and the graph/irreducibility
hypotheses for the final orderability implication. It is not a theorem
about all smooth higher-degree normalizations or complementary pencils.
In particular matching meridians of two pieces satisfying Section1
is a sufficient obstruction.

## 3. The actual residual r=3 pieces and the meridian slope

PF75 leaves coprime e=2b,o=3t, with t odd, (b,t)!=(1,1), r=3.
Put A=bt,c=6bt-2b-3t,k=3c+1,D=k/2-A. Its ACTUAL main cable is
T2 times I, so the actual end glues

    E': base D, cones(b,2b,t), boundary(m,L),
    E_1: the actual positive T(2b,3t) knot exterior, boundary(n,L_1),
    m -> -n, L -> (k/2)*n+L_1.                             (5)

Order1 cones are omitted. The ramification mark(k,3) remains as in
PF75; the UNMARKED end already gives the contradiction here.
The crucial common slope is the genuine meridian m, mapped to -n,
which is the same projective slope as n.

Neither piece is a solid torus. The base orbifold of E_1 has Euler
1/(2b)+1/(3t)-1<0. For E', if b=1,t>=3 its genuine cones are2,t;
if t=1,b>=2 they are b,2b; otherwise they are b,2b,t. In every case
its disk-orbifold Euler number is negative. Both pieces are orientable
Seifert manifolds with incompressible boundary and are irreducible.
One standard product-cover argument verifies these facts: a finite
torsion-free cover of the disk orbifold is a surface Sigma of negative
Euler number with boundary; the pulled-back oriented circle bundle is
Sigma times S1. Boundary injection and irreducibility descend. Equivalently
the finite surface cover can be obtained from the kernel of the map
from the free product of cone cyclic groups to their finite product.

Gluing two irreducible manifolds along incompressible boundary tori
is irreducible, by the usual innermost-circle sphere argument. Thus
the actual closed end in (5) is an orientable irreducible graph manifold.
PF73 also gives finite H1=Z/D plus Z/b. Here D>=8 for b=1,t>=3,
and D>=5b-4>=6 for b>=2,t>=1, so that group is nontrivial.
No virtual companion filling is substituted into (5).

## 4. The ordinary torus-knot side

The actual positive torus-knot exterior E_1 is a regular-fibre complement
in the weighted Hopf fibration of S3. Its normalized negative Seifert
data have two genuine cones p=2b,q=3t,

    omega_p=-q^(-1) mod p, omega_q=-p^(-1) mod q,
    e0=-1, epsilon=-1/(pq).

Choose least positive residues. Since p,q are coprime,
q*omega_p+p*omega_q=pq-1, as follows from the congruences and
0<q*omega_p+p*omega_q<2pq-p-q. Thus the two fractions sum
to1-1/(pq)<1. Section1 applies: n, the regular-fibre meridian,
is NOT in L^o(E_1). This does not need a knot-surgery classification.

## 5. All cubic-companion meridians are outside the interval interior

PF73 identifies W=E'(m) with the oriented link of the explicit normal
cubic model

    w^3+u^t*w+v^b=0,

with primitive positive weights (2b,3t,bt). It has base S2, cone orders
b,2b,t and epsilon=-1/(2bt). The actual regular companion fibre is
H=3A*m+L and in H1(E') equals2A*m, with no torsion part. Since m
intersects H once, the filling core in W is an ORDINARY regular fibre.
Thus E' is a regular-fibre complement in this ACTUAL W.

The normalized negative-slice convention from PF66/73 gives the
following complete coefficients:

    u0=-(3t)^(-1) mod 2b, 0<u0<2b,
    omega_b=u0 mod b (0 if b=1),
    omega_(2b)=u0,
    omega_t=-(2b)^(-1) mod t (0 if t=1).                    (6)

For clarity, these are derived from the actual weighted action, not
from H1. At w=v=0 the isotropy is2b and transverse coordinate v has
weight3t. At v=0,w^2+u^t=0 the isotropy is b with the same transverse
weight. At u=0 the isotropy is t and transverse coordinate u has weight2b.
The oriented normal-slice Seifert coefficient solves transverse_weight
times omega=-1 modulo the isotropy, exactly the inverse convention
already calibrated against the complex A-chain in PF66/73. Coprimality
gcd(2b,3t)=1 makes every nontrivial pair primitive.

It remains to prove the simple inequality e0+nu>=1 in Section1 for
EVERY non-trefoil b,t, without a singularity classification. Define

    delta=floor(u0/b) in{0,1},
    a=(3t*u0+1)/(2b), N=ceil(a/t).

a is a positive integer less than3t. Formula(6) says
omega_t=N*t-a (including t=1), and omega_b=u0-delta*b. Therefore

    -e0=sum omega_i/alpha_i+1/(2bt)=N-delta.               (7)

If delta=0, then u0<b, so a<2t and N<=2. If delta=1, a<3t
gives N<=3. Thus -e0<=2 in every case, and -e0>=1 since it is
a positive integer.

If b,t>1, there are three genuine cones and e0>=-2, proving
e0+nu>=1. If b=1 then u0=1 and N=2,delta=1, so e0=-1.
If t=1,b>1, a is1 or2. For a=1 the formula u0=(2b-1)/3 gives
delta=0; for a=2, u0=(4b-1)/3 gives delta=1. Again e0=-1.
These cases have exactly TWO genuine cones, so e0+nu=1.
The excluded trefoil b=t=1 has only one genuine cone and E' is a
solid torus; it is precisely outside the gluing argument's hypotheses.

The slope coordinate is checked using the ACTUAL boundary maps:
in basis(H,-m), the rational longitude eta=H-2bt*m has slope1/(2bt),
which equals -epsilon as in Rasmussen equation(8). Hence m is exactly
slope0, with no hidden integer section translation. The endpoint
lemma now shows m is outside L^o(E') at all non-trefoil parameters,
whether W is rational or not.

## 6. Entire common-source closure, with the full goal still open

Both actual meridians in (5) are outside their respective L-space
INTERIORS, and the gluing matches them. The general obstruction (4)
therefore shows the actual r=3 end is not an L-space and its nontrivial
group is LEFT-ORDERABLE. It cannot be a PF74 common-coordinate source
end group. This excludes EVERY residual non-trefoil r=3 COMMON source,
at all b,t and every source order and pole, not only the bounded binary
controls or any selected target key.

PF75 excludes every other common winding. PF74 handles all trefoil
common sources, and PF68/71/72 handle the other representation sectors.
Thus the ENTIRE actual primitive-binomial TWO-STAGE QUARTIC COMMON
source class is CLOSED. PF50's common/complementary dichotomy makes
COMPLEMENTARY the only possible source pencil for any remaining
realization in that class. This does not exclude those remaining cases.

[The checker](verify_matched_meridian_common_exclusion.py),
[receipt](matched_meridian_common_exclusion_receipt.json) and
[manifest](sources/r3-seifert-gluing/manifest.json) retain exact modular,
Seifert, rounding, rational/non-rational and actual-boundary controls.
The proofs, not a bounded test or numerical interval sweep, give the
uniform conclusion. Complete primary PDF pages were visually inspected.
Rasmussen's gluing theorem and interval formula, Nemethi's rational-link
and graph-orderability equivalences and the retained PF74 source inputs
are IMPORTED deep topology, not independently reproved here.

The degree-free matched-slope and endpoint lemmas keep their actual
Seifert/graph/source hypotheses. Higher-degree normalizations are not
automatically smooth or common-pencilled. Passing complementary source
sectors, other/non-binomial keys, longer sequences, higher degree and
the full research goal remain OPEN. PF63 polynomial immersions still
exist; there is no finite Keller map, full quartic exclusion, formal
verification, independent full reproof, external review or novelty/
priority claim.
