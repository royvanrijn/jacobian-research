# Trefoil caps and complementary boundary-fibre fusion

15 November 2026. PF2RULED78 excludes EVERY trefoil S3/S3 complementary
source in the actual primitive-binomial TWO-STAGE QUARTIC class.
Together with PF76 and PF77 this closes the ENTIRE S3/S3 representation
sector in BOTH source pencils, at every first pair and winding.
Other representations of the same target keys can remain open.

The new source reduction has no degree, order, support or characteristic
depth bound. It requires PF77's actual complementary coordinate pair,
smooth full source, exactly TWO entire omitted lines, one fold R and one
unramified D, and algebraically exact volume divisor R. It does not
assert these hypotheses for arbitrary higher-degree packets. No full
quartic or higher-degree exclusion is claimed.

## 1. The first D fibre and the only possible corner fusion

Keep PF77's actual corner basis (h_x,h_f), its R-side models and the
WHOLE unramified chart

    x-lambda=w^sigma,
    f=h(w)+w^(1-sigma)*V,
    h=sum c_j*w^-j,  sigma+1<=j<=p,
    c_p!=0, gcd(sigma,support h)=1.                       (1)

Let g=gcd(p,sigma). At the first contact the branch (u=1/f,x-lambda)
has leading orders (p,sigma). The Euclidean extraction with the ORIGINAL
u=0 axis retained gives, at the corner boundary, the primitive first
D-node fibre and its initial side-cone order

    F_D=-(sigma/g)*h_x+(p/g)*h_f,   alpha_0=p/g.           (2)

These are actual torus slopes. One way to keep the framing is to start
with the adjacent normal fibres h_x,h_f and transport across each
degree-two vertex of weight -a by

    h_next=a*h_current-h_previous.                       (3)

The Euclidean algorithm for (p,sigma) ends at the primitive vector
(-sigma/g,p/g); the off-continuation chain has determinant p/g.
The negative coefficient is forced by (3) and the positive plumbing
exchange. This also treats integral contact p/sigma: the initial node
then has side cone p/sigma, and later characteristic nodes follow it.
If sigma=1 there is only the one disk piece with cones (p,p).

For sigma>1, a nonintegral first contact with g>1 leaves one further
characteristic continuation. The initial D piece is then an annulus
with its one cone p/g. If g=1 there is no later characteristic block:
the disk piece has cones

    (p,kappa),  kappa=p+1-sigma>=2.                       (4)

The R-boundary side is PF77's annulus with cone q=rho-1 for rho>=3,
and its corner fibre is

    F_R=q*h_x-rho*h_f.                                   (5)

For rho=2 it is a product collar, to be removed. For rho=1 it is instead
an annulus with cone2 and fibre 2*h_x-h_f. Its original R filling has
distance3 from the regular fibre, since R was the order3 leg core.

For rho=1, equality of the two corner fibre slopes would require
2p=sigma, impossible by (1). For rho>=3, the corner pieces can merge
in their given orientable-base fibrations EXACTLY when

    q*p=rho*sigma.

Since gcd(q,rho)=1, this means

    p=rho*a, sigma=q*a, g=a.                             (6)

In particular p/sigma=rho/q is NONINTEGRAL. There is no earlier integral
contact node. If a>1 the merged boundary piece is an annulus with TWO
CONSECUTIVE genuine cones

    (q,rho)=(q,q+1).                                     (7)

If a=1, D has a single characteristic block and the merger gives a
disk with THREE genuine cones (q,rho,kappa). This distinction retains
the number of boundary circles as well as the cone orders.

## 2. Internal D edges and the Klein exception

Every exceptional boundary vertex in the minimal D extraction is
created before D and is subsequently blown up at least once: the next
centre lies on the most recently created exceptional curve. Its final
weight is therefore at most -2. In particular all degree-two vertices
BETWEEN D nodes have weight at most -2. The first original section can
be handled in the initial collar of (2); it is not between later nodes.

Across a bridge with internal weights -a_1,...,-a_l, equation (3)
makes the coefficient of the second initial basis vector in the final
node fibre the determinant of the chain [a_1,...,a_l]. For l=0 it is1.
For a_i>=2 that determinant is positive: its recurrence is
D_l=a_l*D_(l-1)-D_(l-2), with D_0=1, D_1=a_1, and
D_l>D_(l-1)>0. Thus adjacent D fibres are NEVER parallel.
This is an unrestricted negative-chain argument, including integral
contact and every later characteristic block, rather than a sampled
strict-cabling assertion. The final filling leaves its node fibre
unchanged and adds the cone kappa_D=j_*+1-sigma>=2.

All annulus-cone pieces and every terminal disk piece other than (2,2)
have unique boundary fibre slopes by Hatcher Lemma1.15. Consequently
their internal tori cannot disappear. The remaining terminal disk(2,2)
is the twisted I-bundle K over the Klein bottle. It MUST be treated.

Here is an elementary way to dispose of this exception for the present
target, without assuming that its alternative fibre cannot align.
Write

    pi1(K)=<a,b | a*b*a^-1=b^-1>,
    pi1(boundary K)=<a^2,b>.

A primitive boundary cyclic subgroup normal in pi1(K) is generated
either by a^2 or by b: conjugation by a changes (m,n) to (m,-n)
in a^(2m)*b^n, so normality of a primitive rank-one subgroup requires
m=0 or n=0. These are the disk(2,2) and Mobius-band fibre slopes.
The first is the original plumbing fibre. If its edge is retained,
the terminal JSJ piece is K. If that edge can be removed by matching
the OTHER slope, the merged piece has a NONORIENTABLE Mobius-band
base with at least one genuine cone, contributed by its annulus-cone
neighbour. Further mergers keep that nonorientable base and its cones.
There is no product or solid-torus absorption here: both sides have
incompressible boundary and the neighbour has a genuine cone.

Neither outcome can be a target terminal trefoil piece. K has
H1=Z plus Z/2, whereas a trefoil exterior has H1=Z. A Mobius-base
piece with a genuine cone has TRIVIAL group centre: its base group
is Z*C_alpha*..., with trivial centre; any central element must be
an infinite-fibre power, and a Mobius loop inverts that power.
The trefoil group, and the target annulus two-cone group, instead have
infinite cyclic centre. These facts also exclude the case in which
all source pieces merge into one nonorientable-base piece.
Hence a source homeomorphic to one of the target exteriors below
CANNOT have a terminal Klein piece. All remaining D edges survive,
and (6) is the only possible source merger after collar removal.

For the orientable pieces used here, their fibre slope is unique and
their orbifold cone orders are homeomorphism invariants: quotienting
the group by the infinite central fibre gives a free product of the
finite cyclic cone groups and the base free group. Maximal finite
subgroups recover the genuine orders with multiplicities. This also
distinguishes disks with three cones from disks with two cones.

## 3. The actual trefoil target after its solid companion is capped

Put e=2,o=3,r=3s, with s positive odd. PF73/75 give

    k=3s+1, D=(s+1)/2, H1(closed target end)=Z/D.          (8)

The main companion E' is an ACTUAL solid torus, with meridian
eta=m+L. Before attaching it, drill ONLY R from the main cable:
its base has THREE boundary circles and one cone s. On the central
torus the regular fibre is (k/2)*m+r*L. Therefore the actual E' cap
has meridian-fibre distance

    J=r-k/2=(3s-1)/2.                                   (9)

The target R exterior is now an annulus piece with cone orders(s,J),
attached to the actual trefoil exterior E_1. For s>=3 both cones are
genuine. The fibres at the remaining edge are

    h=k*n+r*L_1, H=6*n+L_1,
    Delta(h,H)=6r-k=15s-1>0.                            (10)

Both pieces have negative orbifold Euler characteristic and unique
boundary fibres. This is exactly ONE canonical JSJ edge; the marked
boundary piece is the annulus with TWO cones(s,J), and the terminal
piece is the trefoil disk(2,3). The R meridian meets h once.

For s=1 both orders in (9) are1. The annulus piece is a product, so
the actual R exterior is E_1 itself. PF75's actual r3 collar gluing
m=-n,L=2n+L_1 sends the solid companion meridian to n+L_1.
Equivalently the regular R fibre becomes the winding-one core of
that filling; after drilling, the product transports its meridian
to n+L_1. The marked meridian-fibre distance is consequently

    Delta(n+L_1,6n+L_1)=5.                              (11)

No companion exterior has been replaced on the basis of H1 alone.
Only the ACTUAL trefoil E' solid torus is capped in this argument.

## 4. Peripheral recognition applies also at s=1

Use PF77's stable full exterior peripheral comparison for the FULL
two-prime N. All source pieces, and the target pieces above, are
irreducible and boundary-incompressible; injective tree gluing gives
the same for each exterior. The target is aspherical. For s>=3 its
internal torus is incompressible. For s=1 the trefoil's vertical annulus
separating its two cones is incompressible and boundary-incompressible
by Hatcher Proposition1.13. It is a properly embedded two-sided
incompressible surface, so this target is sufficiently large as well
(Waldhausen Definition1.1.7, not a requirement for a closed torus).

Realize the peripheral group isomorphism by a map of pairs with a
homeomorphism on the single boundary torus. Waldhausen Theorem6.1
gives a boundary-preserving homeomorphism: the exceptional product
over a closed orientable surface would have TWO boundaries, and the
covering has degree1 because the map on pi1 is an isomorphism.
The normal R meridian is preserved. Thus the canonical rooted pieces,
their genuine cone orders and meridian-fibre distances must match.
Filling that meridian also preserves the closed end homology.

## 5. The smallest winding s=1 is impossible

The target exterior is the single trefoil piece. Section2 excludes
terminal Klein alternatives. For rho=1 the nonparallel corner leaves
an internal JSJ torus. For rho>=3 an unmatched corner also leaves one;
a matched corner either leaves a later D edge, or produces a disk with
three genuine cones. None is a trefoil exterior.

The only remaining source possibility is the rho=2 collar with a
SINGLE D node. Its disk cones must be (2,3). Smooth sigma=1 gives
(p,p), impossible. Equation (4), with p>sigma and kappa<=p, forces

    sigma=2, p=3, kappa=2.                              (12)

At its boundary the actual R meridian is h_x-3h_f (PF77 with q=1),
while (2) gives fibre -2h_x+3h_f. Their distance is3, contradicting
(11)'s distance5. Independently the full source matrix is

    Q=[[-3,3],[3,-8]], det Q=15,

whereas (8) gives target end homology order1. These are checks on
the forced residual model, not a bounded search for all models.

## 6. All larger windings reduce to s=3, and it too is impossible

For s>=3 the boundary piece must be an annulus with TWO genuine cones.
After excluding terminal Klein alternatives, rho=1 gives only one
boundary cone and cannot merge. The rho=2 collar gives a D annulus
with ONE cone if there is continuation, or a DISK if there is none.
For rho>=3 an unmatched corner again gives only one cone. Thus (6)
must hold with a>1, and (7) must match the target pair (s,J).

Since J-s=(s-1)/2, these target orders are consecutive ONLY at s=3.
Then J=4, q=3,rho=4, and

    sigma=3a, p=4a, a>=2.                              (13)

The target has only one remaining JSJ piece, its terminal trefoil.
All internal D edges survive, so there must be EXACTLY two D nodes:
the first gcd drops from 3a to a and the second from a to1.
The terminal D disk has cone orders (a,kappa_D), which must be (2,3).
There are just two possibilities, forced without bounding a:

* a=2,kappa_D=3. Then sigma=6 and the last drop would be
  j_*=sigma-1+kappa_D=8. But p=8 already gave the first drop to2;
  index8 cannot break that remaining gcd2. Impossible.
* a=3,kappa_D=2. Then sigma=9,p=12,j_*=10. Support must be
  EXACTLY {12,10}: all indices are at least10, index11 must vanish
  or the last drop would occur there instead. Both coefficients are
  arbitrary nonzero; their sizes do not affect the contact orders.

In the second case the conjugate contact sum is

    S_D=(9-3)*12+(3-1)*(24-10)=100,
    n_D=S_D+2p-sigma+1=116.

PF50/61's FULL complementary source matrix has n_R=13 and cross
q*p=36. Hence

    det Q=13*116-36^2=212.                              (14)

But the target at s=3 has homology order D=2. This contradiction
finishes every trefoil S3/S3 complementary winding. PF76 supplies
the common closure; PF77 supplies all non-trefoil complementary cases.
Thus EVERY S3/S3 representation in the stated actual-binomial two-stage
quartic class is excluded in BOTH pencils.

## Evidence and limits

The [checker](verify_trefoil_complementary_fusion.py) seals PF77 and
its dependencies, checks actual corner slope transport and internal
node bridges on integral/multiple-characteristic/Klein controls, exact
fusion and target-cap identities, both forced residual lattices, and
ordinary full-satellite trefoil groups/marks. The
[receipt](trefoil_complementary_fusion_receipt.json) and
[manifest](sources/trefoil-complementary-fusion/manifest.json) retain
the primary page audits. Hatcher PDF20/21 and Waldhausen PDF5/23 were
visually rechecked completely; PF77's JSJ uniqueness audit is retained.
The written Euclidean, negative-chain, Klein and cone-count arguments
prove the unbounded reduction. The finite controls do not prove it.
Deep topology is imported, not independently reproved.

The earlier trefoil next-gate checkpoint in tmp/wright77 remains an
EXPLORATORY historical input, not the adopted certificate. Other
complementary representation sectors, other/non-binomial keys, longer
sequences, higher degrees and the full research goal remain OPEN.
No finite Keller map, formal verification, independent full reproof,
external review or novelty/priority claim is made. PF63 polynomial
immersions still exist and are not excluded by this source obstruction.
