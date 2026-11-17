# The regular-coordinate quartic sector is excluded

14 November 2026. PF2RULED44 proves:

> No finite quartic plane Keller map has an actual normalization whose
> sole omitted prime has index two and for which a source coordinate
> extends regularly to that normalization.

PF2RULED39 handles regular-coordinate order one. PF2RULED42 handles
higher order with one source characteristic pair, including every prime
boundary order. This note excludes every remaining multiple-characteristic
case, uniformly in all Laurent coefficients and cabling depths.

The first part is independent of field degree: the marked surgery model
of PF2RULED43 is now recognized as the actual target-radial pair. The
quartic proof uses the retained exhaustive S4/A4 peripheral census,
not a larger coefficient search. The absence of a regular coordinate,
additional omitted primes and higher finite Keller degrees remain open.
No general Jacobian-conjecture or full-goal completion is asserted.

## 1. Setup and the source data that must be preserved

Let U=A2_(x,y) be the marked source in smooth affine N, with sole
complement R=A1 and div(dx wedge dy)=R. For r=v_R(y)>=2, write the
whole cyclic root charts in the now sufficient form PF2RULED43:

    y=z^r, x=z^(2-r)*t+h(z),
    h=sum_(j=r-1)^q c_j*z^(-j), c_q!=0,
    gcd(r,support h)=1, q>=r-1.                         (1)

Finite Keller exactness forces c_r=0. The source reference branch is
Gamma=(w=1/h(z),y=z^r). Read its strict descending support gcd drops
j_i, g_0=r, g_i=gcd(g_(i-1),j_i), g_s=1, and put

    beta_i=2q-j_i,
    2delta=sum_i(g_(i-1)-g_i)*beta_i-r+1,
    F=sum_i(g_(i-1)-g_i)*beta_i+beta_s,
    kappa=j_s-r+2, n=F+kappa=2delta+2q+1.               (2)

The compactification end pair is n-surgery on Gamma with the omitted
knot k its surgery dual. H1=Z/n, n odd>=3, and k is primitive.
In the knot exterior E_Gamma, use its S3 meridian/longitude mu,lambda:

    meridian(k)=n*mu+lambda,
    outer regular fiber h_Gamma=F*mu+lambda,
    Delta(meridian(k),h_Gamma)=kappa!=2.                 (3)

These slopes are part of the marked data, not recoverable from H1 alone.
Assume henceforth Gamma has s>=2 characteristic pairs. Its exterior is
a nontrivial chain of cable spaces and a torus-knot exterior. Write
p_s for the last winding, alpha_s for its cable slope and F_prev for
the previous cabling fiber framing. The positive algebraic inequality is

    alpha_s>p_s*F_prev, F=p_s*alpha_s.                 (4)

The source Gamma must never be confused with the target nonproper curve
or its knot at infinity.

## 2. Actual radial recognition in every finite field degree

### 2.1 A free-product obstruction for the closed model group

We use an elementary group lemma, with its proof included. A connected
tree of injectively amalgamated noncyclic groups, each with a nontrivial
center, and with nontrivial edge groups, has no nontrivial free-product
splitting.

Indeed construct the usual free-product tree for a hypothetical A*B:
vertices are factor cosets, edges are group elements; edge stabilizers
are trivial. A nonidentity elliptic element fixes a unique vertex.
A hyperbolic element has a unique axis and cyclic centralizer: a
centralizer preserves and translates its axis, and the kernel would
fix an edge and is trivial. Thus a noncyclic subgroup with nontrivial
center fixes a vertex. Neighboring vertex groups in our amalgam fix
the same vertex, since their intersection contains the nontrivial
edge group. The whole group then fixes a vertex, contradicting a
nontrivial free-product action. This is a normal-form/tree argument;
no inference from finite group controls is used.

For the n-surgery model M_0, if kappa>=3, Ito's cable formula retains
the terminal Seifert D2(p_s,kappa) piece attached to the companion
exterior. Every piece has noncyclic group and central infinite fiber,
and every torus edge injects. The gluing graph is a chain. The lemma
therefore makes pi1(M_0) freely indecomposable.

If kappa=1, the outer cable collapses to companion filling n/p_s^2.
PF2RULED43 gives

    D=n-p_s^2*F_prev>=p_s*(p_s+1)+1>=7.                 (5)

For a multiple-characteristic companion its terminal piece is
D2(p_(s-1),D), with injective gluing to the rest of the chain; the same
lemma applies. For a torus-knot companion the filled model is closed
Seifert. Its central fiber is nonzero in H1: the filling numerator
is n, the meridian generates Z/n, and its fiber class is F_prev times
that generator, with 0<F_prev<n. A nontrivial free product has trivial
center, so this group too is freely indecomposable.

All Seifert boundary pieces here have negative orbifold Euler
characteristic; the distance-two and product exceptions are absent.
These facts also prove model irreducibility: glue irreducible pieces
along incompressible tori, or use the closed Seifert criterion. The
exceptional reducible Seifert spaces have infinite or noncyclic even
first homology, incompatible with Z/n of odd order.

### 2.2 From stable peripheral groups to the marked pair

For an actual finite Keller map, the retained cofinal end-pair comparison
gives isomorphisms of the closed group, exterior group and full peripheral
system, with the omitted meridian preserved. Apply the preceding lemma
to this closed group. Kneser--Milnor prime decomposition allows at most
one nontrivial fundamental-group factor. Perelman/Poincare identifies
any simply connected prime summand with S3. Finite H1 excludes S1xS2.
Thus the actual closed radial section M is irreducible.

The marked knot is nonzero in H1(M), so its exterior E is irreducible:
a ball bounded by a sphere cannot contain that knot. Its group is the
noncyclic knot group pi1(E_Gamma); boundary compression of an irreducible
one-torus-boundary manifold would make it a solid torus, impossible.
Hence E is also boundary-irreducible.

The target E_Gamma is irreducible, boundary-irreducible, aspherical and
sufficiently large: it contains a genuine incompressible cabling torus.
The peripheral isomorphism is realized by a map of pairs to this
aspherical target, with a homeomorphism on its boundary. Apply
Waldhausen1968 Theorem6.1. Its exceptional source product over a closed
orientable surface has two boundary components; E has one. The resulting
cover has degree one since pi1 is an isomorphism. The theorem preserves
the boundary map throughout the homotopy. We obtain a homeomorphism

    E ~= E_Gamma preserving meridian(k),
    (M,k) ~= (M_0,surgery dual).                        (6)

It extends over that marked filling. This closes the general radial
recognition gap left expressly open in PF2RULED43. It does not classify
a three-manifold from H1 or from an arbitrary abstract group. The
irreducibility and peripheral hypotheses have each been proved.
The argument has no quartic field-degree assumption.

## 3. Essential annuli determine the source boundary fiber

Record three elementary facts about an iterated torus-knot exterior.
All its Seifert pieces have negative orbifold Euler characteristic:
a torus leaf is a disk with two coprime genuine cones; a cable piece
is an annulus with one genuine cone. The two fiber slopes on any
cabling torus are distinct, since equality would require a=p*F_prev,
contrary to gcd(a,p)=1. They form its canonical minimal JSJ chain.

Put an essential annulus or incompressible torus in minimal position
with this chain. Innermost-circle arguments remove inessential
intersections; boundary-parallel annular components are removed by
minimality. Remaining annuli in each Seifert piece are vertical by
Hatcher Lemma1.14 (none of its product/Klein exceptions occurs).
Across a chain torus they would require the two distinct fiber slopes
to agree. Consequently an essential annulus with both boundary curves
on the exterior boundary lies in the outer piece, with boundary slope
its regular fiber. Every incompressible torus is isotopic to a chain
torus. Horizontal annuli/tori cannot appear because their Euler
characteristic zero contradicts the negative base-orbifold characteristic.

Now let Z be a Seifert submanifold containing the exterior boundary,
whose other boundary tori are incompressible, and whose base is an
annulus with genuine cones or a three-holed sphere with possible
cones. Unless Z=T2xI, a vertical annulus based at that exterior boundary
is incompressible and boundary-incompressible by Hatcher Proposition1.13.
It stays essential after gluing the other pieces: a compressing or
boundary-compressing disk can be made disjoint from the incompressible
gluing tori by innermost-circle reduction. Hence its fiber slope is
h_Gamma.

When Gamma has s>=2, such a Z must be exactly the outer cable piece,
up to collars. All its other tori are nested chain tori. It cannot
contain two consecutive chain pieces, since their nonmatching fibers
preclude a single Seifert fibration (Hatcher Lemma1.15). Multiple
nested inner boundary tori cannot all bound the connected side
containing the exterior boundary. Thus Z has one inner boundary and
one genuine cone of order p_s. In particular, on the marked boundary,

    Delta(meridian(k),fiber(Z))=kappa,
    number of genuine cones in Z=1.                   (7)

If a one-boundary companion side compresses, it is a solid torus.
To see this inside the irreducible E, its exterior-boundary component
lies outside that side. Spheres on the side bound balls on the same
side, since the component containing the exterior boundary cannot
lie in such a ball. Thus the side is irreducible; compression gives
a solid torus. This excludes the compression-ball alternative without
a closed lens-space or closed atoroidality assumption.

A fiber filling is treated separately. In a simple S4 pattern its
fiber represents k or 2k; in an A4 pattern it represents k or 3k.
For primitive k in H1=Z/n, n>3 odd, these classes are nonzero except
that no assertion about 3k is needed at n=3. Fiber filling cannot
kill them in the cases used below.

### All-degree local ramification-fiber bound

The multiplicity bound used later is itself independent of degree and
of preferred-longitude identity. In any connected branched component
above a target cable piece, let tau be its simple branch transposition
and let h be the central regular-fiber permutation. Commutation with
tau preserves its two-letter support, so h^2 fixes both marked letters.
Since h^2 is central on the transitive component, it fixes every letter.
The fiber degree is therefore one or two. The marked branch core has
residue degree one, so its Seifert multiplicity is that fiber degree:

    multiplicity of the simple ramification fiber is in {1,2}.

If this lifted piece survives as an essential nonproduct piece at the
source exterior boundary, Section3 forces kappa in {1,2}; PF2RULED43
then forces kappa=1. This is an all-degree conditional obstruction.
Do not conclude that the branched orbit has only two letters in the
nonidentity-fiber case without an additional longitude condition:
h may also permute residual letters. For example h=(12)(34) commutes
with tau=(12) and with alpha=(13)(24), on a transitive four-letter
pattern orbit. Its marked fiber still has multiplicity two. The
quartic census has the stronger identity-longitude hypothesis and
therefore legitimately uses h=tau^(p*a).

## 4. No unbranched degree-four cover between iterated knot exteriors

There is no connected unbranched degree-four cover E_L -> E_J when
both L,J are nontrivial finite iterated torus knots (arbitrary signs).
This is a cover obstruction, not a statement about arbitrary graph
manifolds or links.

A finite cover pulls back the canonical JSJ chain. Negative orbifold
Euler characteristic excludes new product/Klein pieces, and distinct
fiber slopes remain distinct under torus covering. Its pullback is
therefore still the minimal JSJ decomposition. A source cable piece
has two boundary tori, one above each target cable boundary; a source
torus leaf cannot cover a cable base, which has at least two boundary
preimages. Since the source graph is a chain with one exterior boundary
and one leaf, it cannot branch. A two-boundary lift of the target torus
leaf would reflect the chain back toward the target exterior boundary;
each intervening lifted cable has one boundary above each target
side, so the reflected path must reach a second physical source
boundary. That is impossible. Thus it cannot fold either.
There is exactly one component
above each target piece, each of covering degree four. In particular
a torus-knot exterior covers a torus-knot exterior in degree four.

Let its regular-fiber degree be f and its base-orbifold degree b,
so f*b=4. The target cone orders A,B>=2 are coprime. The source is a
disk with two coprime genuine cone orders. The two base monodromies
alpha,beta on b letters satisfy

    c(alpha)+c(beta)=b+1;
    source cones are A/len(alpha-cycle), B/len(beta-cycle).

For b=2, one monodromy is identity, producing two equal genuine cones;
this cannot be the required two coprime cones. For b=4 the cycle counts
are (1,4) or (2,3), up to exchange. The first produces four equal genuine
cones. In the second, a (2,2) and a (2,1,1) force both A,B even. A (3,1)
and a (2,1,1) instead produce at least three genuine cones, including
one of order A and two of order B. All are impossible.

For b=1,f=4 the source has the same orders A,B. In meridian/zero-longitude
bases its boundary map is

    mu_L -> u*mu_J+v*lambda_J,
    lambda_L -> w*lambda_J, |u*w|=4.

The fiber +/-AB*mu_L+lambda_L maps to +/-4 times the target fiber.
Thus |u|=4, |w|=1, and AB divides 4-1 or 4+1. But AB>=6. This excludes
the remaining case and proves the unbounded-depth cover obstruction.
Its small permutation control checks only the finite base-degree layer.

## 5. Solid partial three-sheeted covers force the trefoil

Say an iterated torus knot has absolute algebraic growth if at every
genuine cabling stage its signed slope a and winding p satisfy

    |a|>p*|F_prev|.                                    (8)

Positive algebraic knots and their mirrors have this property.
Every iterated torus knot with exterior homeomorphic to E_Gamma has
it as well: the JSJ chain, each cable cone order p and the absolute
homology coefficient |F| of its boundary fiber are intrinsic. These
recover |a|=|F|/p and all previous |F_prev| recursively. This argument
uses exterior topology and does not assume knot-complement recognition.

For a full S3 representation of its knot group with transposition
meridian, fill the unramified degree-one boundary of the three-sheeted
exterior cover. The result Y has one remaining degree-two boundary.
If Y is a solid torus, the knot must be the trefoil or its mirror.

Here is an induction on arbitrary finite cabling depth. On a torus
knot, central fiber image is identity. Exceptional monodromies are a
three-cycle and a transposition, giving target orders 3A,2b and lifted
closed orders A,b,2b. The ramification fiber is regular, so Y has disk
base with these cone orders. A solid torus has at most one genuine
cone. Thus b=1 and A=1, giving the trefoil (with either orientation).

At a genuine outer cable let B contain the ramification knot. If its
branched orbit has three letters, removing that knot leaves a Seifert
piece with incompressible boundaries, glued only to unbranched companion
exterior covers, also boundary-incompressible. Their groups inject Z2
into pi1(Y), impossible for a solid torus. This includes both the
transposition and three-cycle exceptional monodromies; in the latter
case the companion image may be S3 or C2, and both are retained.

For a two-letter orbit, p is odd. The remaining sheet caps the companion
cover to its own partial three-sheeted cover Y_J. Incompressibility
would again inject Z2 into the solid Y; by the compression argument
of Section3, Y_J is solid. Induction makes J a trefoil or its mirror.
For even a, the drilled pattern has annular base with two genuine
p-cones; neither a fiber filling nor a nonfiber filling makes it solid.
For odd a, it has annular base with one p-cone, with ramification the
order-two fiber. Fiber filling has noncyclic C_p*Z quotient. A nonfiber
solid filling therefore has distance one from the regular fiber.

Use the actual degree-two boundary basis m0=2mu_J,l0=lambda_J. For a
signed trefoil epsilon=+/-1, the retained exact meridian is

    meridian(Y_J)=epsilon*m0+l0,
    pattern fiber = a*m0+2p*l0.

Distance one forces |a-2p*epsilon|=1, hence |a|<=2p+1<6p. This contradicts
(8) over a trefoil. No genuine cable is possible. Notice that the lemma
did not assume finite H1 for the closed three-sheeted cover.

The growth condition is essential. The positive (3,5)-cable of the
positive trefoil has an actual full-S3 simple representation: use the
two-letter pattern and the trefoil companion, whose meridian maps to
the same transposition and longitude to identity. They agree on the
outer pattern boundary since p=3 is odd. Its partial exterior is solid
because |5-2*3|=1. But 5<3*6, so it fails (8). This is a topological
control, with no finite Keller interpretation.

## 6. An auxiliary A4 fiber-distance bound

Consider an actual transitive A4 four-sheeted cover of a finite iterated
torus knot, with three-cycle meridian. Suppose its marked exterior is
homeomorphic to a nontrivial algebraic-knot exterior E_L, with meridian
m, primitive marked knot, and H1(M)=Z/n for n>3 odd. Then

    Delta(m,outer fiber of E_L) belongs to {1,3}.        (9)

The longitude is automatically identity: C_A4(three-cycle)=C3 and
A4'=V4. For a torus target the lifted ramification fiber is regular,
and uniqueness of the boundary fiber gives distance one. An unknot
cannot have full A4 image.

For a cable the retained exhaustive A4 census gives just these types.
A nonidentity fiber has a three-letter branched orbit: B is a disk
with cones p,3 and the marked knot has multiplicity three. With trivial
fiber, a three-letter orbit has disk cones p,p,p and regular marked
knot. A four-letter orbit has an annulus base and two outer boundaries,
with no unbranched caps. Drill the marked knot and compress/fill any
compressible one-boundary companion side as in Section3. Fiber filling
would kill k or 3k, nonzero at n>3. Any remaining negative-base Seifert
piece supplies an essential annulus at the marked boundary, forcing its
fiber to be the intrinsic boundary fiber of E_L. Thus the distances
are three or one, respectively. If every other side was solid, the
whole exterior is Seifert; it is either a nontrivial torus-knot exterior
with the same unique fiber, or is incompatible with the assumed E_L.
There is no product exception: at least a genuine p-cone remains in
the three-letter cases, and the four-letter drilled base has three
boundaries whose companion covers are all incompressible.

This is not a claim that every A4 cover is impossible. It is a bound
for its actual marked meridian against the intrinsic source fiber.

## 7. First quartic reduction: only terminal distance one can survive

A hypothetical finite quartic Keller map in our sector gives full S4
monodromy, transposition meridian and identity preferred longitude by
PF2RULED36. Its target knot is neither torus nor unknot by PF2RULED43,
so choose a genuine outer cable of winding p>=2 and signed coprime slope
a over nontrivial companion J. Use the exact PF2RULED35 pattern census.
The local marked fiber multiplicity is one when the fiber monodromy
is identity, and two in the nonidentity two-letter case.

Drill the marked knot. Every incompressible nonproduct branched piece,
after solid companion fillings, has its boundary fiber equal to
h_Gamma by Section3. Thus kappa is one or two; (2) excludes two.
If a one-boundary companion is solid and the entire exterior becomes
Seifert, it cannot equal our multiple-characteristic E_Gamma.
The S4 companion in the three-letter/three-cycle pattern has a dual
graph cycle and infinite closed H1, so it is impossible.

The only remaining product possibility in the census is the
three-letter/three-cycle pattern with S3 companion, p=3 and distance-one
filling of its partial three-sheeted companion Y. Its drilled piece is
T2xI and E_Gamma=E_J, the fourth companion sheet being an actual identity
copy of E_J. By Section5 J has absolute algebraic growth, and solid Y
forces it to be a trefoil. But E_Gamma has multiple characteristic pairs
and a nonempty JSJ chain. This product case is impossible too.
Consequently

    kappa=1.                                          (10)

## 8. The distance-one source peels to another regular-coordinate model

By PF2RULED43, (10) is equivalent to c_(r-1)!=0 and

    g=gcd(r,{j>r-1:c_j!=0})>1.

Since Gamma has at least two characteristic pairs, g<r. Thus

    h(z)=H(z^g)+c*z^(-(r-1)), c!=0, 1<g<r,
    r'=r/g>=2, q'=q/g,
    gcd(r',support H)=1, coefficient of index r' in H is zero.

All higher indices are multiples of g. Exactness removes the index r,
so q>=r+g and q'>=r'+1. The shortened datum H gives an ACTUAL smooth
affine regular-coordinate model by PF2RULED43. Its source reference
knot Gamma' is precisely the companion of the last cable of Gamma:
its characteristic exponents are beta_i/g for i<s.

Let F',kappa',n'=F'+kappa' be its invariants. The preceding support
drop is at j_(s-1)>=r+g, so kappa'>=3. The final cabling and surgery
identities are

    alpha_s=g*(n'-2)+1,
    Gamma = (g,g*(n'-2)+1)-cable of Gamma',
    n=g^2*(n'-2)+g+1,
    meridian on E_Gamma' after collapse = n*mu'+g^2*lambda',
    D=n-g^2*F'=g^2*(kappa'-2)+g+1>=7.                 (11)

These are all-degree geometric identities, not a descent of a finite
Keller map to the shortened affine model. In particular n>=31.
We now exclude each quartic pattern without assuming such a descent.

## 9. Exhaust every distance-one quartic cabling pattern

Write C_g for the actual outer cable piece of E_Gamma, an annulus
with one genuine cone of order g. Section3 identifies each essential
nonproduct marked pattern with C_g, not just its fiber slope.

### 9.1 Four-letter branched orbit

There are no unbranched caps. A four-cycle exceptional monodromy gives
a drilled three-boundary piece with all companion boundaries
incompressible; it cannot be C_g. Two two-cycles give two genuine
cones unless p=2, whose literal peripheral subgroup obstruction was
already proved in PF2RULED35. A (1,3) exceptional monodromy has cone
orders p,p/3. Exactly one genuine cone requires p=3, g=3. The other
side is then an unbranched degree-four companion exterior cover

    E_Gamma' -> E_J.

Both are nontrivial iterated torus-knot exteriors. Section4 excludes it.
Thus no four-letter orbit remains.

### 9.2 Three letters, transposition exceptional monodromy

The drilled annulus has cones p/2,p. If the partially capped companion
Y were solid, the whole marked exterior would be Seifert, impossible.
It is therefore incompressible and the drilled piece is C_g. Exactly
one genuine cone requires p=2 and g=2. The companion image is A4.

The original B is a solid torus. On its actual boundary lattice
m0=3mu_J,l0=lambda_J, its meridian is m0; this is the exact integral
identification of PF2RULED35, not a guessed surgery slope. Replace B
by the standard index-three branch solid torus W_3. The boundary
meridians agree, so the closed manifold is unchanged. The new marked
knot k' is its core and the new cover of J has full A4 monodromy.
The old regular marked fiber winds twice around that core. Since n
is odd and the old class is primitive, k' is primitive too.

Its exterior Y is E_Gamma'. Its filling meridian is, by (11),

    m'=n*mu'+4*lambda',
    Delta(m',h_Gamma')=n-4*F'=4*(kappa'-2)+3>=7.

But the actual A4 fiber-distance bound (9) permits only one or three.
Contradiction. No equality of finite Keller algebras is asserted by
this topological replacement.

### 9.3 Two-letter branched orbit

The companion side has one boundary after its two unbranched caps.
If solid, the full marked exterior is Seifert, impossible. Otherwise
it is incompressible. For even a the drilled pattern has two genuine
p-cones, and cannot be C_g. For odd a it has one p-cone but marked fiber
multiplicity two, forcing kappa=2, already impossible. Fiber fillings
kill k or 2k and were excluded separately. Thus no two-letter orbit
remains.

### 9.4 Three letters, three-cycle exceptional monodromy

The full-S4 companion has the already excluded dual graph cycle.
For the S3 companion, one sheet is exactly E_J and the other three,
after their unbranched cap, form a one-boundary Y. Before filling Y,
the drilled branched piece has three boundaries and a possible cone
c=p/3. If Y were incompressible it could not be C_g, so Y is solid.
A fiber filling kills k. A nonfiber filling gives an annulus with
cone orders c and b=Delta(meridian(Y),fiber), ignoring order one.
The product c=b=1 was excluded in Section7. Exactly one genuine cone
leaves only the following two possibilities.

* c>1,b=1: g=c and p=3g.
* c=1,b>1: p=3 and g=b.

The identity companion sheet is homeomorphic to E_Gamma'. Its S3
meridian and zero longitude agree up to signs with those of Gamma'.
Here is a direct justification, avoiding a knot-complement theorem.
An exterior homeomorphism sends longitude to +/-longitude and
meridian to +/-meridian+v*longitude. Boundary fiber uniqueness and
its homology coefficient give |F_J|=F'>=6. Comparing the two fiber
coefficients forces F_J*v to be zero or +/-2; hence v=0. On this
boundary the target pattern fiber is a*mu_J+p*lambda_J; the source
outer cable fiber is alpha_s*mu'+g*lambda'. Thus their longitude
coefficients force

    p=g.                                              (12)

The first possibility p=3g contradicts (12). In the second, g=3.
Since E_J is homeomorphic to E_Gamma', it has absolute algebraic
growth. Section5 makes J a signed trefoil, and its exact partial-cover
meridian gives

    b=|a/2-3*epsilon|.

But b=g=3 makes a divisible by three, contrary to gcd(p,a)=1 at p=3.
The evenness of a and the displayed meridian lattice are precisely the
retained S3 pattern conditions. This excludes the last case.

The real (3,16)-trefoil S4 cover retained in PF2RULED41 is consistent
with the new obstruction: its drilled filled pattern has cone order
b=5, whereas the identity-boundary fiber has longitude coefficient
p=3. It fails (12), so its marked exterior is not an iterated S3
knot exterior. Its lens H1 of order seven and primitive ramification
therefore supply no counterexample to this theorem.

Every S4 peripheral state and every nonidentity-fiber pattern has now
been covered. Together with Sections7-8, no multiple-characteristic
regular-coordinate quartic model exists. PF2RULED39 and PF2RULED42
supply the other regular-coordinate sectors. The theorem follows.

## 10. Remaining objective, sources and assurance

A remaining actual quartic normalization with a sole index-two omitted
prime has NO regular source coordinate. The present argument does not
show that such a coordinate must exist. Additional omitted primes and
higher-index packets must retain their full finite algebra: deleting
one generally loses finiteness. The radial-recognition theorem,
low-degree exterior-cover lemma, solid three-sheeted-cover lemma and
source peeling identities and the local ramification-fiber bound are
all-depth structure; only the final
S4/A4 census-based exclusion is quartic. Higher finite Wright degrees
and the full objective remain open.

Use the finite target algebra, the all-degree residual/sign-cover data,
or the geometry of valuations with no regular coordinate for the next
step. Do not resume quartic coefficient searches on any regular-coordinate
Laurent model, including the distance-one controls: their finite quartic
realizations are now globally excluded. The smooth affine controls and
their exact primitives remain valid geometric examples.

All primary statements, exact bytes and new visual page audits are in
[sources/regular-quartic-exclusion/manifest.json](sources/regular-quartic-exclusion/manifest.json).
Hatcher's essential-surface and fiber-uniqueness criteria were checked
on complete printed pages19-20; Waldhausen Theorem6.1 on printed77
was rendered and visually checked, including its boundary-relative
conclusion. Retained Scott fiber-cover and Ito cable-filling formulas,
Perelman only for simply connected summands, and Lokteva's algebraic
cabling inequality are imported with their stated limits.

[verify_regular_quartic_exclusion.py](verify_regular_quartic_exclusion.py)
writes [regular_quartic_exclusion_receipt.json](regular_quartic_exclusion_receipt.json).
It seals the predecessor proofs, checks the small cover/peripheral
permutation layers, the marked replacement lattice, source peeling,
covering boundary equations and signed final slopes. The universal
three-manifold arguments are written proofs, not inferred from those
controls. Formal verification, independent full reproof, external
review and novelty priority remain false. No finite Keller map is
constructed, and the research goal remains active.
