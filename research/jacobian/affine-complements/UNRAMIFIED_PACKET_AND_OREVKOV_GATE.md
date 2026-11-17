# The remaining quartic packet and an all-degree Orevkov gate

14 November 2026. PF2RULED46 retains the original finite normalization.
It sharpens the next gate after PF2RULED45; it does not exclude it.

Every hypothetical geometric-degree-four plane Keller map now has exactly
two omitted primes, R of index two and D unramified. Both have residue
degree one and zero excess. Their polynomial image normalizations are
immersive, their image supports C and H are distinct, and C has a double
image value. In particular every such quartic map is non-surjective.

In every degree, if the sole ramified prime has index two and immersive
image C, the Orevkov invariant O(C) surjects onto S_d. Thus every primary
criterion making O(C) abelian is an obstruction to the FULL finite map,
even when other, unramified omitted primes are present. Positive braids
at infinity and the precise Neumann--Norbury sufficient class below are
excluded. A separate explicit Nori calculation retains the cost of
tangential double points in the remaining quartic case.

The local argument is more general: with any number of ramified primes
and any generic indices, if all their image normalization maps are
birational and immersive, the Orevkov invariant of the FULL branch union
surjects onto the nonabelian finite monodromy group. Consequently the
same abelianity obstructions apply to that entire union, without an S_d
assumption or deleting any prime.

We also classify all smooth affine one-line plane extensions with
nowhere-vanishing source volume. This describes every individual
unramified omitted-prime open in every field degree. Its exactness test
does not restore finiteness after other primes are removed.

## 1. The quartic ramification alternatives are already exhausted

Use the original normal affine finite normalization pi:N->A2 of a Keller
map F on its marked U=A2. Omitted primes are disjoint. The retained
Orevkov prime-indexed Euler identity is

    d-1=sum_i(e_i+E_i), E_i=sum_(a in R_i)(ell_a-e_i)>=0.       (1)

At least one prime ramifies. Orevkov's final Remark in the three-sheeted
paper excludes a component of index d-1. At d=4, an index at least four
already exceeds (1), and index three is excluded by that Remark.
Consequently there is one e=2 prime R, and the remaining budget is one.
The alternatives are precisely

    sole R, E_R=1; or R and D, e_D=1, E_R=E_D=0.               (2)

This is QPK4, not a new high-index classification. PF2RULED45 excludes
the first alternative, with PF2RULED44. Only the second remains.
In particular higher-index QUARTIC packets are not an additional open
case. Higher indices in higher degrees remain open. The earlier current
navigation's broader quartic wording is corrected here without altering
sealed historical scopes or receipts.

PF2OMIT6/7 give smooth N, two disjoint affine lines R and D, residue
degrees one, and immersive image normalizations in the second row. Locally pi
at R is a fold (s,t)->(s^2,t); at D it is etale. All special fibre
factors at R and D have lengths two and one respectively. Tangencies
between distinct image branches do not increase either local length.

## 2. A finite-cover quotient that keeps all unramified primes

First retain arbitrary field degree d>1 and any number of ramified
primes, with arbitrary generic indices. Assume every ramified prime has
residue degree one and immersive image normalization. Equivalently in
this Keller-boundary setting, every ramified prime has zero excess
(PF2OMIT6). Let B be the reduced FULL branch-image union. At each point
of N on a ramified prime the local cover is cyclic over its one smooth
image branch, by PF2OMIT2's local cyclic lemma. This remains true when
several image branches meet: their normalization points upstairs are
distinct, and each local finite factor ramifies over only its own branch.
Every local branch meridian consequently acts by cycles on disjoint
finite factors; different local branch meridians commute. Thus the full
monodromy representation factors as

    O(B) ->> G,                                            (3a)

where G is the actual transitive finite monodromy group. The indices
and the number of ramified primes are unrestricted in this statement.
We do not assert G=S_d for arbitrary cycle profiles or shared supports.

Here G must be nonabelian. Indeed, a transitive abelian permutation
group is regular: an element fixing one sheet fixes all sheets, hence
is the identity. Abelian G would therefore make the original function
field extension Galois. If C0 were any branch-image component, every
prime over its generic point in the finite Galois normalization would
have the same ramification index greater than one. But let b0 define C0.
The polynomial b0(P,Q) on the marked source A2 is nonconstant, hence
has a curve of zeros. The Keller map is etale and quasi-finite, so each
such curve dominates C0 and supplies an UNRAMIFIED generic prime over
it in the same normalization. This is a contradiction. Thus the Galois
cover has no branch component. Purity makes it finite etale over A2,
and connected finite etale triviality forces degree one, contrary to d>1.
This is the elementary Galois Keller obstruction, proved here to avoid
assuming symmetric monodromy in (3a).

Therefore O(B) must have a nonabelian finite quotient for EVERY
zero-excess ramified-boundary packet. In particular it cannot be abelian.
The hypothesis concerns the full branch union, not the nonproperness
union, and not the separate Orevkov invariants of its components.

Specialize now to d>=3, with full ramification divisor of pi a
sole index-two prime R with residue degree one, and its normalization
maps immersively onto its image C. Other omitted primes are allowed,
but must be unramified. Their images need not equal C. Purity implies
that pi is finite etale over A2 minus C, INCLUDING its points on those
other omitted primes. We do not replace N by U union R here.

The connected cover over A2 minus C has a transitive monodromy group.
A generic meridian of C is a single transposition. The target plane is
simply connected, so these meridians normally generate its complement
group. Their conjugate transpositions generate a connected sheet graph,
and hence the monodromy is S_d. This argument does not require degree
four or a ruling on N.

At a singular target value a of C, each preimage on R is a separate
smooth fold point. The local finite cover splits into these degree-two
factors and unramified degree-one factors. Thus every local branch
meridian acts on a DISJOINT sheet pair. The entire local monodromy image
is abelian, even if the image branches are tangent. Its local commutator
subgroup is killed. Therefore monodromy factors as a surjection

    O(C)=pi1(A2 minus C)/<<[G_a,G_a]:a in Sing(C)>> -> S_d,    (3)

where G_a is the small-ball local complement group. This is exactly
Neumann--Norbury Definition1, PDF page1. In particular O(C) cannot be
abelian. In this setting the Euler identity also says

    number of unramified omitted primes = d-3,                 (4)

since the immersive simple prime has E_R=0 and unramified primes have
zero excess. Formula (4) is a count of omitted primes, not image supports.

This is an application of classical topology, with no novelty assertion.
It uses the full finite map to justify both global and local monodromy;
an abstract colliding polynomial curve or a deleted affine piece does
not supply that map.

## 3. Primary abelianity tests and their exact scope

Use Neumann--Norbury's splice convention for the TARGET knot at infinity:
(p_i,q_i), with q_i the successive winding, and define their canonical
number (called d_i in that paper, not the field degree here) by

    a_1=p_1+q_1,
    Delta_i=p_i-p_(i-1)*q_(i-1)*q_i,
    a_i=q_i*a_(i-1)+Delta_i.                               (5)

Section5.2, PDF12-13 identifies a_n=-K_X.Cbar-1 on its specified
efficient infinity compactification X. One must keep this infinity
convention; the SOURCE characteristic knots in PF2RULED43 are different.

Theorem2, with Section5.2 Theorem5.2 and the proof on PDF15, makes O(C)
abelian in the stated sufficiently negative class. We use the explicit
sufficient conditions

    a_n>=p_1*q_2*...*q_n-1 and a_n>1.                       (6)

Both inequalities are retained, including small/degenerate conventions.
Thus (6) is impossible in (3). Their positive-braid theorem also makes
O(C) abelian for every positive braid at infinity (Introduction and
Section5.3). Such links are excluded in ALL finite field degrees under
Section2, independently of any number of unramified primes or collisions.
Polynomially parametrized knots at infinity need not be positive braids.
We do not infer the omitted positivity hypothesis from polynomiality.

For (3a), apply the same primary theorems to the ENTIRE branch union B:
a positive braid for that whole link, or conditions (6) on every path in
its simultaneous infinity diagram, make O(B) abelian and are impossible.
Separate component diagrams do not retain their mutual infinity data.
In particular a multi-prime zero-excess packet cannot escape this
obstruction merely by adding unramified omitted primes or by changing
its ramification indices.

For a nodal C there is a stronger direct test. On X its normalization
is P1 and its only finite singularities are k nodes. Adjunction gives

    Cbar^2=a_n-1+2k.                                       (7)

The infinity boundary is transverse to Cbar. If a_n>1, Nori
Proposition3.27 gives an abelian kernel for filling Cbar. The filled
complement is A2 and is simply connected, so the whole complement group
is abelian, contrary to S_d. Therefore a nodal C necessarily has

    a_n<=1.                                               (8)

For example the numerical recurrence for a two-stage splice with first
pair (p,q) and second (s,r) is a_2=s+r*(p+q-p*q). Thus nodal candidates
require s<=1+r*(p*q-p-q). These are necessary topological inequalities;
no arbitrary integer splice is asserted to have a polynomial realization.

We use the nodal form of Nori directly. We do not use its bare nodal
self-intersection inequality for a curve with tangential double points.
For those points the next section retains every finite blowup loss.
Orevkov's separate 1988-paper Added-in-2003 withdrawal of a claimed
general Jacobian reduction is not a premise of this argument.

## 4. An explicit all-contact Nori budget for the remaining quartic C

In (2)'s remaining row at most two parameters of R lie over any value,
because each costs two of the four sheets. Every singularity of C is
therefore two smooth branches with some intersection multiplicity m>=1,
analytically A_(2m-1). No cusp, triple point or nonimmersive branch occurs.
Let k be its number of ordinary nodes, and put

    T=sum m over tangential double values (m>=2).             (9)

This sum counts intersection multiplicity, not the number of tangencies.
Its affine delta is k+T. Thus on the infinity compactification

    Cbar^2=a_n-1+2k+2T.

At a tangency use branches v=0, v=u^m. The first m-1 blowups leave
two transverse strict branches at a node of the removed boundary. One
additional blowup separates those branches from that boundary node.
There are exactly m centres of TOTAL C multiplicity two. Include every
exceptional curve in E, the removed divisor. The self-intersection cost
is 4m, the tangency disappears from the strict curve, and its two
branches meet the new exceptional at different smooth points. Ordinary
nodes away from these centres are left unchanged. Hence

    Ctilde^2=a_n-1+2k-2T, r(Ctilde)=k.                       (10)

Xtilde minus E is A2 minus finitely many tangency points; it is simply
connected. Xtilde minus(E union Ctilde) is the original A2 minus C.
Nori Proposition3.27 would make that whole group abelian if
Ctilde^2>2k. Its S4 quotient forbids this. Consequently

    a_n<=1+2T.                                              (11)

This is an unbounded contact budget. In particular a_n>1 forces actual
tangencies with T>=ceil((a_n-1)/2), not merely additional nodes. Passing
(11) or failing the sufficient abelianity condition (6) does not construct
a cover or prove nonabelianness. No equality case is excluded.

## 5. Distinct images and actual missing values in degree four

Write C=pi(R), H=pi(D). If C=H, a value with b normalization branches
would consume 3b omitted sheets: two at R and one at D per branch.
Thus b<=1 everywhere. The common immersive finite birational
normalization A1->C is then an isomorphism onto a closed embedded
affine line. The Abhyankar--Moh line theorem makes C a coordinate line,
so pi1(A2 minus C)=Z. This cannot surject onto S4. Therefore C!=H.

The same argument, using only C and the FULL finite cover pi, shows
that C cannot itself have an injective normalization even when H is
different. Thus C has at least one double value. A line theorem is
applicable because the image is immersive; we are not applying a
Keller covering theorem to F outside C while ignoring H.

Let c(a),h(a) count normalization parameters of C,H over a, zero if
absent. The exact whole-source fibre rule is

    #F^(-1)(a)=4-2c(a)-h(a).                                (12)

It gives the complete finite branch-count table:

* c=0: h=0,1,2,3,4, retained lengths 4,3,2,1,0;
* c=1: h=0,1,2, retained lengths 2,1,0;
* c=2: h=0 only, retained length zero.

Thus every double C value misses the source and is disjoint from H.
Every hypothetical degree-four Keller map is non-surjective. Its whole
missing set is finite and equals the union of: double C values; fourfold
H values off C; and C/H intersections with h=2. Singleton values are
exactly c=0,h=3 or c=1,h=1. Images are immersive but need not be nodal;
all pairwise contact multiplicities survive. No claim that C and H must
intersect, or that they admit arbitrary independent prescribed contacts,
is made.

## 6. Complete all-degree geometry of an unramified affine piece

Let M be smooth affine, with marked U=A2, sole complement D=A1 and
div_M(dx wedge dy)=0. PF2RULED45 Section1 gives a regular source
coordinate without assuming exactness. It is constant on D, since a
source coordinate line cannot acquire an endpoint in an affine surface.
Shift its value to zero and call it y; let r=v_D(y)>=1. Its complementary
coordinate x has a pole. PF2RULED38's arbitrary-automorphism lemma says
that y is unique up to nonzero affine change among regular coordinates.

The full marked normal form is

    y=z^r, x=z^(1-r)*t_i+h_i(z),
    h(z)=sum_(j=r)^q c_j*z^(-j), c_q!=0,
    h_i(z)=h(zeta^(-i)*z), gcd(r,{j:c_j!=0})=1.              (13)

For r=1 the gcd condition is automatic and q>=1. For r>=2 it forces
q>=r+1: the only possible support at q=r would have gcd r.

Necessity. Normalize M in y=z^r, retaining the entire source plane
A2_(x,z). The root is branched along the retained source line y=0,
etale near D, and has r disjoint boundary copies of A1. As in
PF2RULED40, remove that retained line and all but one boundary copy.
This is affine and its reduced z fibres are A1. Miyanishi's affine-line
fibration theorem gives the WHOLE plane A2_(z,t_i), not a formal chart.
The overlap automorphism is x=z^k*t_i+h_i; the volume order zero gives
k=1-r. Gauge removes every exponent of h_i at least 1-r, leaving
exactly j>=r. A nonzero part remains, or both source coordinates would
be regular. Deck equivariance gives (13); equal conjugate Laurent parts
violate separatedness, giving the gcd condition.

Sufficiency. Every datum (13) is actual smooth affine geometry. Glue
the source plane and r boundary planes along z!=0 by (13). Dubouloz
Section2.8 applies with source weight zero and boundary weights r-1.
Both source/boundary directional translations h_i and -z^(r-1)h_i
have a pole, of orders q and q-r+1. Between boundary charts the
translations z^(r-1)(h_i-h_j) have a pole too: their highest differing
index is at least r+1, since r cannot distinguish characters. Thus
EVERY pair satisfies the affineness criterion. Take the finite cyclic
quotient, with z->zeta*z on the source, permuting the boundary planes
and t_(i+1)=zeta^(r-1)*t_i. It is free along those boundaries and is a
reflection along the source z=0. The quotient is smooth affine with
the specified U,D, and its volume has no zero or pole, because

    dx wedge dy=r*dt_i wedge dz.                            (14)

This proves the full marked classification. No target algebra or finite
Keller map is supplied by its affine converse.

## 7. Exactness, equivalence and controls in every boundary order

The two polynomial primitives on the root planes are

    theta_src=r*z^(r-1)*x dz, theta_i=r*t_i dz.

Their difference is r*z^(r-1)*h_i dz, with common residue r*c_r.
Thus exactness implies c_r=0. Conversely zero residue gives Laurent
primitives of all these differences; the same affine Cech argument as
PF2RULED40 glues a regular global primitive on the root surface.
Average over the deck group and descend, including across the retained
branch line by the invariant-form calculation in PF2RULED40. Therefore

    dx wedge dy is globally algebraically exact iff c_r=0.   (15)

In particular an exact r=1 model has q>=2. Every r>=2 has exact models
already with h=z^(-(r+1)). These models are not claimed to admit Keller
maps. Marked equivalence is exactly equality of r and

    c'_j=b*a^j*c_j, a,b nonzero,                             (16)

induced by y'=a^r*y, x'=b*x+p(y). ALL polynomial p are allowed, including
constants when r=1: in the boundary chart

    t'=b*a^(r-1)*t+a^(r-1)*z^(r-1)*p(z^r)

is regular, with regular inverse. Coordinate uniqueness proves the
converse. The actual locally nilpotent field partial_x has kernel C[y]
and boundary expression z^(r-1)*partial_t; it is free at r=1 and has
zero divisor (r-1)D otherwise. No target-algebra invariance is inferred.

For r=1 there is also an explicit ring. Put H=y^q*h(y), H(0)!=0. Then

    M=Spec C[y,W,V]/(y^q*V-W*(W-H(y))),
    W=y^q*x, V=x*(y^q*x-H(y)).                             (17)

Deleting D={y=0,W=H(0)} gives the entire source plane; the other zero
line is retained. The source chart is W=y^q*x, V=x*(y^q*x-H), while the
boundary chart is W=H+y^q*t, V=H*t+y^q*t^2. The two charts cover M.
At y=0 the W derivative is a unit; for y!=0 solve for V. Thus M is
smooth affine. For example y^2*V=W*(W-1) is an actual exact-volume
control with U=A2 and one extra unramified line. Global exactness and
existence of the source additive action cannot by themselves exclude D.

For an actual Keller normalization with several primes, BGV
Corollary11.4(2) makes M=U union D affine. It is smooth and has precisely
these hypotheses; the primitive P*dQ/c restricts to it. So (13)-(16)
apply with c_r=0. Coordinates for different primes may differ; their
individual affine pieces are not finite covers of the target.

## 8. Next gate and assurance

The quartic task is now the actual TWO-prime R,D normalization, with
distinct immersed polynomial images, missing double C values, an S4
quotient of O(C), the abelianity exclusions (6) and the contact budget
(11). Its index-two and unramified individual affine pieces have complete
geometric normal forms. Their two valuations and the common actual
finite target algebra must be compared inside the same source C[x,y].
Deleting D cannot transfer the one-prime finite-cover exclusion.

Higher-degree single-prime cases with a nonimmersive ramified image are
outside Section2's cyclic-factor hypothesis: their cusp local monodromy
can be nonabelian. Multiple ramified primes ARE covered by (3a) when
all ramified images have residue degree one and are immersive, using
the full branch union and its actual monodromy, which need not be S_d.
Positive excess, ramified residue maps or other noncyclic local factors
require new local data. Higher-index and arbitrary higher-degree
exclusions and the full objective remain open.

Primary statements, locations and byte seals are in
[sources/unramified-orevkov/manifest.json](sources/unramified-orevkov/manifest.json).
The relevant new complete Orevkov, Neumann--Norbury and Nori pages were
rendered and visually inspected. The whole-chart and affine-converse
imports reuse the precisely retained PF2RULED40/43 primary statements.
The checker verifies algebraic identities, exhaustive small-degree
packet/branch tables and permutation controls, not unbounded topology.
Formal verification, independent full reproof, external review and
novelty priority remain false. The goal remains active.
