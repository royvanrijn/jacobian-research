# The borderline regular-coordinate sector: marked ends and quartic exclusion

14 November 2026. PF2RULED41 closes the entire q=r-1 sector left by
PF2RULED40, for finite field degree four. It also gives its marked
end in every field degree and a stronger general quartic obstruction.

For every existing smooth affine one-line normalization with
div(dx wedge dy)=R and regular-coordinate boundary order r>=2,
PF2RULED40 gives q>=r-1 and excludes q=r for any finite Keller target.
Here the whole borderline q=r-1 has a unique marked model up to
rescaling, namely the explicit affine quotient N_r of PF2RULED40.
Put q=r-1 and n=qr+1=r(r-1)+1. Its end pair has:

    closed end: a lens space with H1=Z/n;
    omitted knot: primitive, with exterior X_(T(q,r));
    filling slope: n*mu+lambda on that torus-knot exterior.     (1)

Descriptions are up to a simultaneous orientation/mirror convention.
For r>=3 the exterior is boundary-incompressible and atoroidal,
and is not a solid torus. Thus the omitted knot is not a Heegaard
core. The actual punctured graph and filling, not H1 alone, give (1).

A new topological theorem excludes any simple full-S4 cover of a
finite iterated torus knot with identity preferred longitude, odd
nontrivial lens source, primitive ramification knot and atoroidal
ramification exterior. It extends the PF2RULED39 lens/core theorem.
The non-core case needs its own final three-sheeted companion argument.

Consequently no N_r admits a finite quartic Keller extension.
In all finite degrees the torus-infinity lemma still excludes torus
or unknot boundary infinity, coprime boundary degrees and prime
boundary-coordinate degrees on these models.

For every remaining actual quartic normalization with exactly one
index-two omitted prime and a regular source coordinate, one now has

    r>=2,   q>=r+1.                                       (2)

The cases q<r-1 were geometrically impossible, q=r was excluded by
volume, and q=r-1 is excluded here. Higher q, the absence of a regular
coordinate, extra omitted primes and higher field degrees remain open.
These are written proofs with primary imports and exact supporting
controls, not formal verification or independent external review.

## 1. The entire borderline model has a two-open toric description

The normalized Laurent part in PF2RULED40 has support j>=r-1.
If q=r-1 it is c*z^(1-r), c!=0. Rescale x to make c=1.
The full r+1-plane gluing and its deck action in PF2RULED40 are
then uniquely fixed. Thus every existing affine marked model with
these invariants is the explicit quotient N_r constructed there.

Use its root ring

    f(v)=v^r-1,
    B_r=C[z,v,v*f(v)^j/z^j : 1<=j<=q], q=r-1,
    y=z^r, x=v/z^q.

On D(v), let

    s=z*v=x*y,  t=f(v)/(z*v)=w_1/v^2,  a=v^r=1+s*t.

All higher generators are v^(j+1)*t^j. Therefore

    (B_r)_v=C[s,t,v,v^-1]/(v^r-1-s*t).

The cyclic action fixes s,t and sends v to zeta^-1*v. Its quotient is

    W=Spec C[s,t,(1+s*t)^-1].

Together with the entire source U=A2_(x,y), this covers N_r, since
D(v) contains all omitted root lines. On the overlap,

    a=1+s*t,   x=a/s^q,   y=s^r/a.                       (3)

The omitted line is s=0, with every finite t present. This is a
global open chart, not only an equality of generic valuations.

## 2. Compactify and keep the punctured graph

Start with P1_x times P1_y. Its deleted original boundary is
D_x={x=infinity}, D_y={y=infinity}; both have square zero.
At the corner x=infinity,y=0, refine the toric fan by the rays

    (1,1), (1,2), (2,3), ..., (q,q+1).

After the first two insertions, each new ray is the sum of
(1,1) and the previous last ray. The target divisor E has ray
(q,r). On its open torus chart the coordinates are

    s=x*y,   a=x^r*y^q.

The exponent matrix is unimodular since r-q=1. Blow up the free
point s=0,a=1. In the chart a=1+s*t, the last exceptional is
Rbar={s=0}. The complement of the strict transform of E is exactly
the whole W in Section1 near Rbar, including every finite t.

The deleted boundary B consists of D_x,D_y and all the toric
exceptionals; it excludes Rbar. Hence S minus B is the actual N_r.
The deleted graph, in spatial order, is the linear graph

    D_y(0)--D_x(-1)--E_(1,1)(-r)--E(-2)--[-2]_(q-1),     (4)

and Rbar(-1) meets E once transversely. There are q+2 blowups in all.
For q=1 the last tail is empty. The retained strict transform of
y=0 has square -2 and is not a deleted boundary component.

Apply the retained, visually audited Neumann Proposition2.1 R6
move at D_y--D_x. There is one remaining branch, one attachment,
genus zero and no puncture on the removed pair, so there are no
extra connected-sum summands. Both the closed graph and the graph
with a fiber deleted at E reduce to

    [-r]--[-2]_(q entries), with the puncture at the first -2. (5)

The continued fraction is

    [r,2,...,2]_(q twos)=r-q/(q+1)=(qr+1)/r.             (6)

Thus the closed end is L(n,r) up to orientation. To check the
actual marked class, orient the linear-chain meridians mu_j,
0<=j<=q, so the positive matrix has negative off-diagonal ones.
Its relations give mu_j=(1+j*q)*mu_0 and n*mu_0=0.
The marked knot is the circle fiber at E, namely mu_1=r*mu_0.
Since gcd(r,n)=1, it is primitive.

For the knot exterior, the punctured vertex has two filled arms.
Their orders are r and q: the first arm is [-r], and the second
is a length q-1 negative-two chain. Thus its Seifert base is a disk
with these two cones (order-one entries are discarded).

Here are also the filling slopes. In one consistent orientation,
let h be the regular fiber and x_1,x_2 the arm generators. The
exterior homology has relations

    r*x_1+h=0,   q*x_2+(q-1)*h=0.

It is Z with generator m and

    h=qr*m, x_1=-q*m, x_2=-r*(q-1)*m.
    mu_R=2*h+x_1+x_2=(qr+1)*m.                        (7)

The fiber and the meridian mu_R of the omitted knot intersect once.
The slope mu=mu_R-h has homology m and intersects h once. Filling
it gives a Seifert lens space with H1 zero, hence S3, with regular
fiber T(q,r). Therefore the punctured graph is the torus-knot
exterior, h=qr*mu+lambda, and mu_R=n*mu+lambda.
This proves (1) with the actual filling tracked.
For r>=3 it cannot be a Heegaard core; its exterior is noncyclic.

Primary inputs: Neumann's closed and punctured plumbing calculus,
and the retained Seifert presentations. Their exact retained bytes
and the new rigidity import are listed in
[sources/borderline-regular-end/manifest.json](sources/borderline-regular-end/manifest.json).

## 3. Atoroidality and comparison with a target-radial end

A disk-base Seifert space with two genuine coprime cones is atoroidal.
By Hatcher Proposition1.11 an incompressible torus is vertical or
horizontal. A closed horizontal surface cannot cover a base with
nonempty boundary. A vertical torus enclosing at most one cone is
compressible; one enclosing both cones is boundary-parallel.
It is boundary-incompressible and aspherical, and a properly
embedded essential annulus makes it sufficiently large.

For an actual finite map, the source target-radial end and the
computed plumbing end give cofinal product neighborhoods of (N_r,R).
Their stable systems identify pi1 of the closed sections, pi1 of
the knot exteriors, the peripheral Z2 subgroups and the meridian
of the omitted knot. This is stronger than matching H1.

Retained elliptization recognizes the actual closed radial section
as a lens space from its cyclic group Z/n. It is irreducible.
Since the knot is nonzero in H1, its exterior is irreducible too:
a sphere bounds a ball in the lens space, and that ball cannot
contain the entire knot without making its class zero.

For r>=3 the exterior group is the noncyclic torus-knot group.
A boundary compression would make the exterior a solid torus,
contradicting its group. Thus it is boundary-irreducible.
The stable peripheral group isomorphism is realized by a map of
pairs to the model exterior, with a homeomorphism on the boundary:
the target exterior is aspherical, so the compatible group and
boundary data extend to a map.

Apply Waldhausen1968 Theorem6.1, printed page77. Its target is the
irreducible, boundary-irreducible, sufficiently large model exterior;
the source is irreducible and boundary-irreducible and has nonzero pi1.
The exceptional product-line-bundle alternative has two boundary
components and cannot occur here. The covering in the conclusion
has degree one because pi1 is an isomorphism. Thus the actual
radial exterior is homeomorphic to the model, with its peripheral
marking. This justifies transferring atoroidality, rather than
assuming a product-end group already identifies the three-manifold.
For r=2 use the retained solid-torus/core recognition instead.

## 4. A stronger odd-lens quartic theorem

Consider a connected simple degree-four cover M->S3 of a finite
iterated torus knot K, with:

* full S4 monodromy and identity preferred longitude;
* M a lens space with H1=Z/n, n>1 odd;
* its unique ramification knot k primitive in H1(M);
* E(k) atoroidal, allowing the solid-torus case.

No such cover exists. The exact PF2RULED35 census of all thirty
S4 subgroups and its peripheral residues is reused unchanged.
Write the outer cable as winding p>=2 and signed slope a coprime to p.
Its fiber is p*lambda_J+a*mu_J. We use the retained pattern notation:
B is the branched pattern piece, Y a companion piece after the
unbranched caps, and alpha its exceptional monodromy.
Only the marked-exterior and final companion arguments change.

Two elementary facts will be used explicitly.

First, if a regular fiber h is primitive in a cyclic H1, the
orbifold abelianization after quotienting by h is zero. A lens
Seifert filling with two cone orders A,B has this quotient C_gcd(A,B).
Thus a primitive regular fiber rules out gcd(A,B)>1.
Fiber filling separately kills h and contradicts its nonzero class.

Second, gluing an incompressible companion exterior to an annulus
Seifert cable space with genuine cone order A>=2 gives an essential
non-boundary-parallel torus in E(k). It is not boundary-parallel
because the side containing the sole exterior boundary is that
cable space, not T2 times I. Atoroidality excludes this gluing.
Boundary incompressibility and compression-ball alternatives are
handled as in PF2RULED35: a ball containing the essential marked
knot would make its homology class zero.

### 4.1 Four-letter and three-letter transposition patterns

If B is boundary-incompressible, uncapped companion covers inject
Z2 into the closed lens group. The only viable four-letter solid
pattern is p=3, alpha cycles (1,3), with one cone of order three.
Its marked regular-fiber exterior is a genuine order-three cable
space attached to an incompressible companion cover. The second
fact excludes it. The p=2,(2,2) alternative is already impossible
by the literal peripheral subgroup conditions. Four-cycle patterns
have annular incompressible B and were excluded before.

In the three-letter transposition pattern, B has disk cones p/2,p.
For p>=4, compression forces Y solid; nonfiber lens filling leaves
those two cones, whose gcd=p/2>1. The first fact contradicts
primitivity of the marked regular fiber. Fiber filling kills it.

For p=2, B is a solid torus with one cone of order two. If Y were
boundary-incompressible, its gluing to that genuine cable exterior
would violate atoroidality. Thus Y is solid. The retained exact
meridian lattice m0=3mu_J,l0=lambda_J replaces B by the standard
index-three branch solid torus, giving a lens/core A_n cover of J.
PF2RULED39 excludes every such odd-order alternating type.

### 4.2 Two-letter patterns

Here p is odd, the companion has full S4 and identity longitude.
For even a, B has disk cones p,p and the knot is regular.
Compression makes Y solid. The first fact excludes its primitive
regular fiber (gcd=p>1), and fiber filling kills the marked class.

For odd a, B has disk cones p,2, and the knot is the order-two fiber.
Y is solid by compression. Fiber filling kills 2[k], nonzero since
n is odd and k is primitive. Nonfiber lens filling adds no genuine
third cone. The actual boundary homology map on m0=2mu_J,l0=lambda_J
agrees with that of the standard simple-branch solid torus because
H1(B)=Z for coprime cone orders p,2. Replacement preserves H1=Z/n
and gives a lens/core S_n cover of J. PF2RULED39 excludes it.
No equality of closed manifolds is inferred from equal homology.

### 4.3 The remaining three-letter three-cycle pattern

If the companion image is S4, its dual gluing graph has a cycle,
giving infinite H1. The remaining companion image is S3 on the
branched pair and the outside fourth letter. The fixed companion
sheet is an actual copy of X_J. After the one unbranched cap, the
other three sheets give a one-boundary piece Y. The branched B has
annulus base with one possible cone p/3.

If Y were incompressible, all sides would inject into M. Hence Y
is solid. Its meridian cannot be a fiber slope, since that kills
the marked regular knot. After filling Y, the B side becomes a
disk-base Seifert space. Its gluing to X_J can be a lens space
only if this filled B is solid. If it has a genuine cone, deleting
the marked regular fiber leaves a genuine cable space, contradicting
atoroidality. Therefore both cone orders must be one:

    p=3,  distance(meridian(Y),h_B)=1,
    E(k)=X_J.                                          (8)

Since J is itself a finite iterated torus knot, atoroidality of X_J
makes J a torus knot: every genuine cabling stage gives a
non-boundary-parallel incompressible cabling torus.

### 4.4 A solid partially filled three-sheeted torus cover forces trefoil

The representation on the three companion sheets is full S3 with
transposition meridian and identity longitude. Fill its missing
degree-two boundary in the usual branched way, leaving the existing
unbranched cap in place; call the closed cover M_3. It is Y union
the standard ramification solid torus.

For a full-S3 simple cover of a torus knot, the central fiber dies.
The two exceptional monodromies must be a three-cycle and a
transposition. Their cycle counts are one and two, so the lifted
base is S2 even without first assuming finite H1. Writing the
target exceptional orders as 3A,2b gives source orders A,2b,b.
The lifted Euler number is nonzero, so H1 is finite of order b.
These are the d=3 computations of PF2RULED34, now with base genus
zero checked directly.

Since Y is solid, M_3 is a lens space and its ramification knot is
a Heegaard core. The torus lemma makes that knot zero in H1;
primitivity of a Heegaard core forces b=1. Thus M_3=S3, and the
ramification exterior is a disk Seifert space with cones A,2.
It is the solid Y only when A=1. Hence J is the trefoil or its mirror.

### 4.5 The actual boundary slopes exclude the trefoil case

Let epsilon=1 for a positive trefoil and -1 for its mirror. On its
degree-two boundary use the actual basis

    m0=2mu_J, l0=lambda_J.

The three-sheeted trefoil cover has S3 source and regular-fiber
ramification an unknot. The target fiber 6epsilon*mu_J+lambda_J
lifts to 3epsilon*m0+l0. The source Seifert fibration has exceptional
orders 2,1, so its regular-unknot fiber framing is 2epsilon.
Subtracting that framing times m0 gives the actual meridian of Y:

    meridian(Y)=epsilon*m0+l0.                          (9)

This uses the covering fiber framing, not a guessed preferred
longitude, and follows directly from the oriented Seifert cover.

For p=3 the branched pattern B is T2 times I and a is even.
Let mu,lambda be the basis on the degree-one boundary going to X_J.
The map H1(B)->H1(V) sends both m0 and mu to zero and both l0
and lambda to a generator. Thus across B,

    m0=eta*mu, l0=lambda+c*mu, eta in {1,-1}.

The common fiber has slopes (a/2)*m0+3*l0 and a*mu+3*lambda.
Equality and gcd(a,3)=1 force eta=-1 and c=a/2. Consequently

    meridian(Y) -> (a/2-epsilon)*mu+lambda,
    distance(meridian(Y),h_B)=|a/2-3epsilon|.             (10)

Condition (8) makes a/2=3epsilon+delta with delta in {1,-1}.
The induced filling slope on the trefoil exterior is therefore

    (2epsilon+delta)*mu+lambda.

Its distance from the trefoil fiber 6epsilon*mu+lambda is
|delta-4epsilon|, which is three or five. Filling a trefoil exterior
along this slope has three genuine Seifert cones 2,3,3 or 2,3,5,
with a noncyclic triangle-group quotient. It cannot be a lens space.
This excludes the last pattern, with both mirrors retained.

The terminal target torus knot is already excluded by PF2RULED34
because its ramification class would be zero. An unknot cannot have
full S4 monodromy. All genuine cable cases were exhausted above.
The stated topological theorem follows for every odd n>1.

## 5. A genuine lens-cover control that fails the required marking

The atoroidal marking cannot be dropped. There is an actual simple
S4 topological cover of the (3,16)-cable of the positive trefoil with
lens source of H1 order seven, identity longitude and primitive
ramification class. It does not have the marked exterior of N_3.

One literal monodromy, with letters 0,1,2,3, is

    tau=(01), alpha=(012), delta=alpha*tau=(02),
    companion b=(023), companion A=(23).

Here A^2=b^3=1 and A^-1*b^2=delta is the trefoil meridian.
The companion image is S3 on {0,2,3}, while the pattern uses {0,1,2};
together they generate S4. The outer slopes have delta^3=delta
and delta^-16=1, and the inner preferred longitude is identity.
Thus the pieces define the claimed connected simple branched cover.

Formula (10) gives surgery slope 7mu+lambda, a lens filling of the
trefoil, but distance |8-3|=5 on the pattern fiber. The ramification
knot winds five in that filled solid torus. Its exterior contains
the essential torus separating the order-five cable space from X_J.
Its class is 16mu=2mu modulo seven, primitive; atoroidality fails.
Matching lens homology alone would have produced a false counterexample
to the borderline exclusion.

This target knot is also compatible with polynomial infinity type.
For example p(t)=t^6+t^2, q(t)=t^9+t^7+t^3 has that link at infinity.
After weighted radial scaling the leading terms trace the trefoil
three times, and q(t)^2-p(t)^3 has leading term 2t^16, giving the
normal-meridian winding sixteen. The perturbations have smaller
order, so the tubular braid is the (3,16) cable.
Its derivative gcd is t. The parametrization is birational: its
normalization degree divides gcd(6,9)=3, while degree three would
factor through a cubic polynomial whose degree-two derivative divides
both derivatives, impossible. Its point t=0 is an ordinary cusp, and
that cusp has no other normalization preimage. These checks do not
claim all finite collision profiles or a finite Keller realization.
They retain a control against dismissing this end solely from its
algebraic knot type or the one-cusp count.

## 6. Apply to actual finite normalizations

For a finite Keller map on N_r, its target-radial cover has the
full S_d monodromy and the computed primitive marked class.
The all-degree torus lemma excludes every torus/unknot at infinity.
The retained polynomial-infinity and target-shear arguments then
exclude coprime boundary degrees and every prime boundary-coordinate
degree; each such degree is composite and at least four.

In degree four, PF2RULED36 gives identity preferred longitude
intrinsically, and the boundary curve has one place at infinity,
hence finite iterated cabling. The marked end comparison in Section3
gives a lens source with primitive knot and atoroidal exterior.
The new theorem excludes it. This applies uniformly to every r,
not a bounded list of pole pairs.

Together with PF2RULED39, PF2RULED40 and the coordinate-pole invariant,
every remaining actual quartic normalization with exactly one omitted
prime of index two and a regular source coordinate has r>=2 and q>=r+1.
All Laurent coefficients, support-gcd conditions, c_r=0, norm and
conductor constraints must still be retained. This does not establish
an affine converse for arbitrary data or solve the surviving q>=r+1
finite-algebra problem. The no-regular-coordinate case, extra primes,
higher finite Wright/model degrees and full goal remain open.

The new rigidity theorem is imported from pinned primary bytes and
was visually checked on its complete relevant page. All graph,
homology, slope and literal permutation controls are exact, but the
three-manifold and all-depth proofs are written, not formalized.
