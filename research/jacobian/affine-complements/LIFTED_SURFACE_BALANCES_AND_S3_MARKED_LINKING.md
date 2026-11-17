# Lifted surface balances and the full S3/S3 marked linking form

15 November 2026. PF2RULED73 supplies general all-degree peripheral
surface identities and a linking-form replacement theorem with
torsion retained. They resolve the two-boundary S3/S3 calculation
left open in PF72, at all actual primitive-binomial two-stage
quartic parameters and for every representation in that sector.

Write e=2b,o=3t,r=3s, with s,t odd positive, gcd(2b,3t)=1, and put

    c=6bt-2b-3t, k=3sc+1,
    A=bt, D=k/2-As=(s*(16bt-6b-9t)+1)/2.

Then D>=1, gcd(D,s)=1, and the full ORIENTED marked end is

    H1(M)=Z/D*g direct-sum Z/b*z,
    lambda = <-inverse(s)/D> orthogonal-sum <-t/b>,
    [K_R]=(sk*g,0).                                        (1)

Order1 factors are omitted. The quartic source requires gcd(D,b)=1.
Common-source self-linking and parity give further gates below.
For the trefoil first pair(2,3), common S3/S3 requires r=15 mod48
and r>=63. An UNBOUNDED family of whole-leaf source Gram controls
matches its full marked linking form. A stronger fundamental-group
check distinguishes that family from the target ends. Thus the
lattice/linking conditions alone leave compatible data, while the
explicit family fails the full boundary-group gate. Other source
configurations, and the whole common sector, are not closed.

## 1. Peripheral balances from oriented surfaces, at any degree

Let P->V be an orientation-preserving finite cover of a solid torus,
branched along a smooth interior link with the ordinary oriented
normal models. A meridian disk transverse to the branch link has
an oriented preimage surface. Its boundary is the sum of the
closed lifts of the oriented target meridian, one for each meridian
cycle, and is zero in H1(P). On a boundary torus T_i of covering
degree d_i, let m_i be the least positive power of the target
meridian that closes, and write the same symbol for its lifted
primitive curve. More explicitly, with that power denoted a_i,

    sum_i (d_i/a_i)*[m_i]=0.                               (2)

The preimage may have several connected components; summing their
boundaries proves (2). In particular a single meridian cycle on
all sheets of a connected P implies that its lift M^degree is
ZERO integrally. No Seifert, cyclic-cover or unit hypothesis is
needed for this conclusion.

Similarly, in any finite UNBRANCHED cover E of a knot exterior,
lift an oriented Seifert surface. If the target preferred longitude
has least closed power l_i on T_i and corresponding primitive
lift ell_i, then

    sum_i (d_i/l_i)*[ell_i]=0 in H1(E).                     (3)

The coefficients count the parallel boundary lifts. This is not
an individual-longitude vanishing theorem. Dividing its coefficients
requires an actual torsion argument. Both statements are proved
by the oriented preimage surfaces, independent of any upper bound
on degree, pole orders or cable length. Relations obtained before
filling persist afterwards. A map to the original exterior need
not persist through a meridian filling, as PF72's countercontrol warns.

## 2. Integral replacement WITH an orthogonal torsion linking summand

Let an oriented compact piece E have one torus boundary with basis
(m,l), and actual maps

    H1(E)=Z*[l] direct-sum T, [m]=0,

where T is finite. Glue E to a piece Q along that torus. Let Y_Q
be Q filled along m and let Y_E be E filled along l. If the resulting
closed manifolds have finite H1, Mayer--Vietoris gives

    H1(Q union E)=H1(Y_Q) direct-sum T,
    lambda_(Q union E)=lambda_Y_Q orthogonal-sum lambda_Y_E. (4)

The same assertion on torsion homology works if free homology is
present. Here T=H1(Y_E). The boundary map has NO T coordinates;
this is a necessary hypothesis, not an assumption about all pieces.

For the pairing on classes represented inside Q, use a multiple
bounding chain and cap its extra m-boundaries by integral surfaces
in E. These caps do not intersect the second cycle, kept inside Q;
the identical computation with meridian disks gives lambda_Y_Q.
A torsion cycle in T bounds a multiple entirely inside E, so its
pairing with a Q cycle is zero. Its pairings with other T cycles
are unchanged by filling l, giving lambda_Y_E. This proves the
orthogonal decomposition, not just a group-order statement.

Marks lying in Q retain their classes with zero T coordinate. An
arbitrary interior mark in E must retain its ACTUAL torsion part.
This all-degree theorem strengthens PF65/66's cyclic-companion
special cases when the displayed peripheral hypotheses can be
proved for a noncyclic piece. It asserts no homeomorphism or
contact equivalence.

A separate GENERAL group observation is useful. For any knot
exterior E, abelianization is a peripheral-compatible epimorphism
onto the solid-torus group: the knot meridian goes to its core,
and the preferred longitude to the identity. For any glued Q,
the universal property of the van Kampen pushout therefore gives

    pi1(Q union E) ->> pi1(Q(longitude of E)).              (4a)

More generally any epimorphism of piece groups respecting the
actual boundary maps induces such a quotient of the glued group.
No injectivity of the boundary maps is needed. Thus an infinite
group of the filling model obstructs a finite actual group even
when the replacement only computes homology and linking. Equality
of these groups is neither needed nor asserted.

As a useful separate consequence of (2), a branched cover P of V
with a single meridian cycle, free H1 rank1, and actual second
boundary generator l mapping to the core of V has precisely these
peripheral maps. Its lifted meridian is zero and H1(P)=Z*l+Tor H1(P).
The rank and the boundary generator must be checked. This removes
the Seifert/unit assumptions needed by the sufficient lemma in PF72.

## 3. The complete S3/S3 incidence and its two reduced pieces

In PF63's S3/S3 sector the companion and cable each act on three
sheets with one fixed sheet, and their fixed sheets DIFFER, since
the whole image is S4. Their intersection has two sheets, swapped
by the companion meridian M and cable beta. The whole orbit graph
is a three-edge PATH: a degree1 companion leaf, the main three-sheet
cable, the main three-sheet companion, and a degree1 cable leaf.
Preferred longitude L has identity monodromy, and the central
companion fibre H=a^e=b^o does too. The root is a 3-cycle, so 3|r;
M=beta^r is a transposition, so r is odd. The even torus generator
is a transposition and the odd one a 3-cycle, hence the parameters
in (1). This exhausts all64 periodic conjugacy classes, including
interchanging first exponents.

On the central degree2 edge use m=M^2,L as a basis. On each fixed-
sheet degree1 edge use M,L. Regular fibres are

    H=M^(eo)*L, h=M^k*L^r.

All branch meridian cycles are filled. As in PF72 the ramification
mark is h: the preferred longitude h*mu^(-rk) differs by a filled
ramified meridian power, since rk is even.

### Reduced companion E'

Fill the main companion's fixed-M boundary using the degree1 cable
solid torus. Its remaining boundary is (m,L). Its base is D with
cone orders b,2b,t. The Seifert presentation gives

    H1(E')=Z direct-sum Z/b.                               (5)

Indeed the torsion order of a base-D Seifert piece is product(alpha)
divided by lcm(alpha); here it is b, and only two cone orders can
share a prime, so that torsion is cyclic. Filling the remaining
meridian m gives the closed three-sheet companion branched cover W.
It has base S2, the same cone orders, and Euler number -3/(eo).
Its determinant is b. Comparing with the torsion order in (5)
shows that m is a PRIMITIVE free generator. This uses the actual
POST-filling Seifert manifold; the old exterior pushforward is not used.

The lifted Seifert surface BEFORE the fixed-M filling has relation
2L+L_1=0. After filling, L_1=H and H=(eo/2)m+L. Consequently
3(L+A*m)=0. Since 3 is coprime to b, (5) has no 3-torsion, so

    H1(E')=Z*m direct-sum Z/b,
    L=-A*m, H=2A*m.                                       (6)

The peripheral maps have zero torsion coordinates. The kernel
slope is eta=L+A*m. This is the required saturation of (3).

### Reduced cable P'

Glue the degree1 companion knot exterior to the fixed-M boundary
of the main cable. Its H1 is generated by M and its preferred
longitude is zero, so homologically and for linking of retained
marks this is longitude filling L_1=0, by PF66 (or (4) with T=0).
Denote this longitude-filled MODEL by P'; no fundamental-group
equality with the glued knot-exterior piece is asserted. The main
cable after branch fillings had two companion boundaries
and one cone r/3=s. Longitude filling introduces a cone of order k,
since h=k*M_1+r*L_1 on that degree1 boundary. Thus P' has base D
with COPRIME cones s,k and H1(P')=Z. Its regular fibre has class
sk times a primitive generator; killing it gives Z/s+Z/k of order sk.
Choose g so that h=sk*g.

The lifted meridian disk gives m+n=0, where n=M on the fixed
boundary. Also h=k*n after longitude filling, and on the degree2
boundary h=(k/2)m+rL. Since H1(P') is torsion free,

    n=s*g, m=-s*g, L=(k/2)*g, h=sk*g.                     (7)

The two-boundary relation (2), rather than an unjustified claim
that each meridian is zero, is essential here.

## 4. Full homology and the exact ramification mark

Glue (6) and (7) in the actual common boundary basis. The kernel
eta=L+A*m maps to D*g, so Mayer--Vietoris gives

    H1(M)=Z/D*g direct-sum Z/b, [K_R]=(sk*g,0),
    H1(M)/<K_R>=Z/gcd(D,sk) direct-sum Z/b.                 (8)

The quantity C=16bt-6b-9t is odd and >=7t-6>=1, since b>=1.
Thus D=(Cs+1)/2 is positive. Also 2D=1 mod s and s is odd, so
gcd(D,s)=1. From k=2D+2As and b|A, one gets

    gcd(gcd(D,sk),b)=gcd(D,b).

PF65's actual quartic marked-generator bound now gives the NEW
necessary gate

    gcd(D,b)=1.                                            (9)

The same condition makes the full group in (8) cyclic. R is zero
exactly when D divides k, equivalently D divides2A. It can be
nonzero or nonzero isotropic in other cases; no blanket vanishing
statement survives the two-boundary calculation.

## 5. The oriented linking forms and the sign

Apply (4) to E', with meridian eta=L+A*m and free longitude m.
The second closed summand is exactly W=E'(m). For its oriented
form, use the explicit normal cubic germ

    w^3+u^t*w+v^b=0 -> (u,v).

Its branch is the positive T(3t,2b) knot and its three-sheet simple
monodromy is the unique S3 torus class. The isolated singularity
model is normal (or smooth in the small cases), as in PF66's
quartic model. Primitive weights are (2b,3t,bt), with three orbit
orders b,2b,t, base S2, and Euler -1/(2bt). The order-b leg in
negative plumbing has beta=-(3t)^(-1) mod b. The same inverse-star
calculation and signed convention as PF66 give its self-linking
beta*(1+2t*beta)/b. Rescaling the generator by3t gives

    9t^2*beta*(1+2t*beta)=-t mod b.

It is a unit numerator and the group has order b. Therefore
lambda_W=<-t/b>. At t=1 the complex coordinate U=u+w^2 gives
UW+v^b, the A_(b-1) model with outward <-1/b>, independently
fixing this sign. A factor of order1 has the trivial form.

For the first summand, fill P' along eta. On its actual boundary,
(m,L) has intersection +1: this is the boundary orientation inherited
from the target companion SOLID TORUS. Its primitive null longitude
is lambda_0=(k/2)m+sL. Choose a curve gamma=u*m+v*L with
-su+(k/2)v=1, so [gamma]=g. Then (gamma,lambda_0) has intersection
-1, the knot-EXTERIOR convention. Write

    eta=D*gamma+Q*lambda_0, Q*s=1-D*v.

The integer Q exists since gcd(s,k/2)=1; in particular
Q=inverse(s) mod D. With the PF61 outward linking convention,

    lambda(g,g)=-Q/D=-inverse(s)/D mod Z.                   (10)

Here is a direct sign check. A surface in P' bounding lambda_0
meets the inward push of gamma positively once, since the boundary
basis is negatively oriented. Subtract Q copies of that surface
from the attached meridian disk, whose boundary is D gamma+Q lambda_0.
The resulting chain bounds D gamma and has intersection -Q with
the second gamma. This is also the -inverse-intersection convention
of the positive integral surgery trace when the slope is p gamma+lambda_0.
Thus the sign is not inferred from group order or a lens-space label.

Equations(4) and (10) prove the full orthogonal formula (1), even
before (9); no coprime-primary shortcut is required. The exact mark
therefore has self-linking

    lambda(K_R,K_R)=-s*k^2/D mod Z.                         (11)

## 6. Common-source gates, and the trefoil refinement

Let N be the least nonnegative residue of -s*k^2 mod D, and put
g0=gcd(D,sk). In the actual common source PF71 gives self-linking
strictly between0 and1/2. Its rank2 Gram G has diagonals n_R odd,
n_D even and determinant bD. Comparing the actual representative
of (11) with n_D/(bD) gives

    0<N<D/2, n_D=b*N, b*N even.                            (12)

The order of the R class also forces c=b*C0 and

    n_R*N-b*C0^2=D, gcd(D,N,C0)=g0, n_R odd.                (13)

These are necessary marked lattice equations, not sufficient source
or affine-cover conditions. The complete oriented form (1) must still
match. Complementary sources are not subject to the strict-half bound.

For the trefoil e=2,o=3, b=t=A=1, s=(r/3) odd, one has

    D=(s+1)/2, lambda=<1/D>, K_R=2g.

For D<=8, zero/half/large self-linking or odd n_D excludes the
common source. For D>8, (12) forces n_D=4. If D is even, any
marked isometry must send e_R to2g and e_D to xg. With source
off-diagonal c it would satisfy 2x=-c mod D and x^2=n_R mod D.
Thus c=2v, 4|D, and x^2=v^2 mod D. But determinant D=4n_R-c^2
gives n_R=v^2+D/4, a contradiction modulo D. Hence D is odd.
Characteristic parity n_R odd then gives D=3 mod8. Therefore

    trefoil common S3/S3 requires r=15 mod48 and r>=63.      (14)

The form matching in this argument is essential; self-linking alone
would leave even D controls that fail the full marked form.

## 7. An unbounded family that matches all these lattice/linking gates

For every positive ODD integer j set

    n_R=3j^2+7j+5, c=3j+4, delta=3j^2+4j+4,
    G=[[n_R,c],[c,4]], r_target=6delta-3.

Then det G=delta, delta=3 mod8, and G is positive definite with
the required canonical parity. On the target trefoil end, D=delta,
lambda=<1/delta>, K_R=2g. Define an isomorphism on the source
discriminant by

    e_R -> 2g, e_D -> ((delta-c)/2)*g.

Both G relations vanish, and2 is a unit modulo delta. The source
pairing G^(-1) differs entrywise by integers from the target pairing
on these images: its R diagonal is4/delta, mixed entry -c/delta,
and D diagonal n_R/delta. This proves a FULL ORIENTED MARKED
linking isometry, not merely matching group orders or self-linking.

The matrices also come from explicit whole-leaf cluster controls
in PF67's P1xP1 model:

    R: f=z^(j+1), x=2z^(-2(j+1))+z^(-j)+z^(1-j)*T;
    D: f=w,       x=2w^-2+W.

At the first three shared centres the R multiplicities are j+1
and the D multiplicities are1. After those centres use local
coordinates f and U=f*(x-2f^-2), up to a unit change. Along R,
U=z+Tz^2 and f=z^(j+1). The D creation has R multiplicity1;
j successive crossings with the old b=0 component have multiplicity1;
the final R creation is free, again with multiplicity1. D retains
its leaf throughout. These explicit Euclidean steps give the
displayed Gram for UNRESTRICTED j, with volume orders1,0 and both
retained components leaves. No finite cluster search proves the family.

At j=1, G=[[15,7],[7,4]], target r=63, H1=Z/11 and R=2g.
The target generator has linking1/11. Both source basis classes
map to2g; G^(-1) is [[4,-7],[-7,15]]/11, whose entries are all
4/11 modZ. This control passes the full marked comparison.

These are admissible boundary CLUSTER and lattice controls. Affineness
of the resulting open surface and a genuine finite target algebra
are NOT proved. No Keller map is constructed. The entire end group
is a stronger invariant, and in this family it can be checked exactly.

## 8. Full end groups distinguish this compatible family

The source infinity divisor T, after removing the R,D leaves, has
the following weighted tree. There is a negative chain

    [-2,-2,-(j+2), -2 repeated j times],

and at its SECOND vertex an extra leg [-2,0]. This follows from
the exact blowup sequence in Section7: the first two exceptional
vertices have weight-2; the shared b=0 vertex loses one intersection
to D and j to the subsequent crossings, giving -(j+2); each crossing
vertex loses one to its successor or the final free R, giving -2.
The original x-infinity section has weight-2 and the f-infinity
section weight0, giving the extra leg. No further compactification
components or fixed sheets are discarded.

For a rational-curve plumbing TREE with weights w_v, van Kampen
on the circle bundles over punctured P1s gives fibre generators h_v,
edge commutators [h_v,h_w]=1 and vertex relations

    h_v^(w_v)*product_(w adjacent v) h_w=1.                 (15)

The order of each product is fixed by the punctured base; there are
no stable letters since the graph is a tree. This derivation does
not require a definite matrix or incompressible gluing tori. In the
extra leg, the weight0 END relation gives h_x=1, and its weight-2
neighbour then expresses h_f as the inverse of the second chain
generator. Eliminating these two generators leaves exactly the
negative-chain presentation. Starting at one end of a chain, its
relations express every generator as a power of the first, so the
group is cyclic. The chain determinant is

    (3(j+2)-2)*(j+1)-3j = delta.

Thus the source end group is Z/delta. The boundary orientation is
opposite to the tube of T; this does not change its fundamental group.
This is a group calculation from (15), not a homeomorphism inferred
from its intersection matrix or linking form.

For the corresponding target, b=t=1 makes E' a solid torus: it
has base D with just one nontrivial cone of order2. Its meridian
is eta=L+m. The OTHER degree1 companion leaf is still an actual
trefoil exterior. Replacing it by longitude filling was justified
for homology and linking, not for homeomorphism or fundamental-group
equality. For the present argument its abelianization M_1->1,L_1->0
gives a SURJECTION onto the replacing solid-torus group, compatible
with both peripheral generators. Van Kampen therefore gives a
surjection from the ACTUAL target end group onto the group of the
longitude-filled model P'(eta). Filling P' along eta extends its
Seifert fibration with third order |det(h,eta)|=3delta-2. This
QUOTIENT model has base orbifold

    S2(2delta-1, 6delta-2, 3delta-2).

All three orders exceed2 for delta>=11, and their reciprocal sum
is less than1. Its triangle group is infinite: a hyperbolic triangle
with those angles gives an infinite reflection group and its
orientation-preserving subgroup. The actual target end group thus
surjects onto this infinite triangle group, and is itself INFINITE.

The source and target end groups cannot be isomorphic. This excludes
EVERY source configuration in the explicit odd-j family as a source
for its paired trefoil target, despite its full marked linking match.
It does NOT exclude every source configuration for those target keys,
nor prove that all common-source ends have cyclic fundamental group.
The useful next gate is the full end group/geometry for the other
source configurations, or the actual finite affine algebra.

## Evidence and remaining work

[verify_lifted_surface_s3.py](verify_lifted_surface_s3.py) audits
the complete64 S3/S3 periodic cycle/incidence types, exact peripheral
relations, the unrestricted gluing and marked-group formulas,
orientation identities and cubic linking congruence, ordinary full
satellite controls, the trefoil gate, symbolic matching family and
its full boundary-group mismatch.
[lifted_surface_s3_receipt.json](lifted_surface_s3_receipt.json)
pins the proof, checker, predecessor evidence and
[source manifest](sources/lifted-surface-s3/manifest.json).

The general surface and torsion-linking replacement theorems apply
at all degrees under their stated hypotheses. The explicit formula
(1) is for the actual-binomial two-stage QUARTIC S3/S3 sector only.
It does not assume general lower-weight keys are binomial, or
smoothness and common coordinate hypotheses in arbitrary higher
degree. PF63 polynomial immersions still exist. The whole research
goal remains open, including common S3/S3 beyond the displayed
necessary gates, complementary passing sectors, other keys, longer
sequences and higher degree. No independent full reproof, formal
verification, external review or novelty/priority claim is asserted.
