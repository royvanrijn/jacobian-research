# Proper finite maps identify the marked end; component-open audit

16 November 2026. PF2RULED98 closes two geometric comparison gates in
every degree under the stated hypotheses. First, the individual affine
component opens used by PF45/46/50 are exactly within BGV Corollary11.4.
Second, the actual target-radial section of a whole finite normalization
is the marked compactification plumbing end, by a direct geometric
construction. No recognition of a three-manifold from an abstract group
is needed for this construction.

The second result does not prove that an arbitrary abstract plumbing is
a finite normalization. It begins with the ORIGINAL proper finite map,
and uses its extension over the target projective plane. Source chart,
canonical-volume, monodromy and finite-completion hypotheses remain
separate obligations. The full research objective remains active.

## 1. The component opens satisfy the primary hypotheses

Let N be a connected normal affine complex surface with a nonsurjective
marked open embedding U=A2 into N. BGV, *On endomorphisms of affine
spaces and the Jacobian problem*, author manuscript BGV160,
Corollary11.4(1)--(2), printed/PDF43, states:

* every connected component of N minus U is an irreducible affine line;
* for each such prime D, U union D is an OPEN AFFINE subvariety of N.

In particular the primes are disjoint. These conclusions do not require
a finite target map, a chosen source ruling, smoothness of N, an etale
map on the chosen component, or a restriction on field degree. They
apply to the original finite normalization because finite over A2 is
affine and normality and the marked open are part of its construction.
The corollary's proof removes the other primes by Nagata's theorem,
then uses the one-place statement of Theorem11.2(4) to rule out
intersections. Its use here does not delete a prime from the FINITE MAP.

The complete primary PDF pages39--43 were rendered and visually read,
including Lemma11.1, the full Theorem11.2 and proof, and the corollary.
Theorem11.2 permits a multiple fiber; smoothness alone does not make
that multiplicity one. We do not import a trivial-bundle assertion in
place of the retained PF45 regular-coordinate existence theorem.

For a smooth original Keller normalization with target coordinate
pullbacks P,Q and Jacobian constant c, one has, as rational forms and
therefore everywhere they are regular,

    omega=dx wedge dy=c^(-1)*dP wedge dQ
         =d(c^(-1)*P*dQ).                            (1)

Restriction to U union D preserves this ACTUAL primitive. If D is
unramified its volume divisor on that open is0; if D is index2 it isD.
This verifies the affine/componentwise/exact-volume input in PF45/46.
The smoothness and index assertions still come from their own local
boundary proofs; Corollary11.4 alone does not supply them. Nor does
affineness of U union D make its restriction to the target finite.

## 2. A compactification of the original finite map

Let pi:N->A2 be a finite dominant morphism of connected smooth complex
affine surfaces, of ANY finite degree. Let A_1,...,A_b be specified
proper algebraic curves in N, smooth and mutually disjoint outside a
compact set. For the Wright application these are all omitted primes;
the formulation does not assume any particular ramification indices.

Normalize P2 in the function field C(N). This gives a normal projective
surface Nbar with a finite morphism pibar:Nbar->P2. Over the target
affine chart the normalization is precisely N, by uniqueness of finite
normalization in the same field. Thus

    pibar^(-1)(A2)=N.                                (2)

Resolve the boundary and the closures of the A_i, leaving N unchanged,
to get a smooth projective X and a morphism g:X->P2. Put T=X minus N.
It can be arranged that T together with those closures is SNC near T and each
curve end meets a smooth point of T transversely, away from crossings.
For the omitted affine lines there is one such point per line. This
uses resolution at infinity, not a modification of the marked source.
One can further dominate any preassigned smooth source completion by
resolving the birational comparison entirely at infinity.

Equation(2) implies that the effective divisor

    H=g^*(target line at infinity)

has support EXACTLY T, with a STRICTLY POSITIVE integer multiplicity
on EVERY component of T. A divisor at infinity mapping into the finite
target chart would contradict(2). This positivity is the crucial
proper finite-map input.

Use homogeneous target coordinates [Z0:Z1:Z2], and define

    alpha=g^*( |Z0|^2/(|Z0|^2+|Z1|^2+|Z2|^2) ).       (3)

This is a globally defined smooth nonnegative function on X, zero
exactly on T. On N it is

    alpha=1/(1+|P|^2+|Q|^2).                         (4)

Locally at a smooth T point, or at an SNC crossing, it has the form

    alpha=c(z,w)*|z|^(2a)*|w|^(2b), c>0,
    a>=1, b=0 at a smooth point; a,b>=1 at a crossing. (5)

The positive factor is smooth; after choosing local coordinates all
its derivatives are bounded on a sufficiently small compact chart.
Every actual large target-radial section is therefore a small positive
level of a function with exactly this positive normal-crossing form.

The divisor T is connected. Indeed H is nef, H^2=deg(pi)>0. If its
support split into disjoint nonempty parts, write H=H_1+H_2 accordingly.
Then H_1.H_2=0 and H_i^2=H.H_i>=0, with at least one strictly positive.
Hodge index makes the orthogonal complement of that positive class
negative definite. The other effective nonzero H_j is not numerically
zero (it meets any ample class positively), a contradiction. Hence its
plumbing boundary, and the section in(8), are connected. No surjectivity
theorem for the infinity knot group is needed for this connectedness.

## 3. Local rounding gives the ACTUAL marked plumbing boundary

Here is the elementary neighborhood argument, including what must be
preserved. Take the oriented real blowup of X along T: replace each
normal coordinate z by polar coordinates (r,theta), and at a crossing
replace both normal coordinates. This is a compact manifold with
corners with interior N. Its boundary faces are normal circle bundles
over the punctured components of T; at a crossing the two circle
coordinates are exchanged. These are exactly the disk-bundle plumbing
pieces. Neumann's Section1, complete PDF4--5/printed302--303, fixes
this base/fiber exchange, orientations and punctured boundary tori.

In a corner chart (5) becomes

    alpha=c*r^(2a)*s^(2b).                          (6)

An inward vector field with positive r component near r=0 and positive
s component near s=0 satisfies

    V(log alpha)=2a*V(r)/r+2b*V(s)/s+V(log c)>0       (7)

in a sufficiently thin boundary neighborhood. At a smooth face only
the first term is needed. Indeed derivatives of log c stay bounded,
whereas the positive normal term tends to infinity at each approached
face. Choose the local fields only on sufficiently thin face/corner
charts, then patch them by a nonnegative partition of unity. Positivity
is preserved. This constructs a field inward at every boundary face
and transverse to all sufficiently small positive alpha levels.

Round the corners by any fixed sufficiently small smooth inward
boundary transverse to this field. Each field orbit intersects this
rounded boundary and each smaller sufficiently near-boundary regular
alpha level exactly once. For a fixed small alpha level, move the
rounded boundary farther inward if necessary; compactness makes the
required choices uniform. The hitting-time map is a diffeomorphism
between smooth boundaries. Existence and uniqueness follow from
strict increase(7), the inward collar and the fact that alpha tends
to0 at its boundary endpoint. Smooth dependence follows from
transversality. Hence the positive multiplicities and the smooth
positive factor alter the shape of the rounding, not its boundary
three-manifold or the plumbing attaching maps.

The curve closures meet only smooth boundary faces. Near each mark
choose coordinates T={z=0}, A_i={w=0}; choose the field in (7) tangent
to w=0. The disjoint marked circles in the faces persist under both
rounding and hitting-time maps. The construction extends to their
small tubular neighborhoods. It therefore identifies the PAIR with
all marked curves, the positive normal meridians, the full boundary
tori and any specified geometric peripheral maps by transport.
It does not equate a target preferred longitude with an independently
chosen source framing. That equality, when needed, must be computed
from the actual boundary covering map, as in PF60/80/92.

Orient a section M as the OUTWARD boundary of the compact core in N.
This agrees with its target branched-cover orientation. It is the
opposite of the boundary orientation of a tubular neighborhood of T
inside X. Use this same outward-core orientation for the source
plumbing end. The comparison then preserves orientation and the
sign of every linking form. Do not silently compare the two opposite
boundary conventions.

This proves, for all sufficiently large R,

    (pi^(-1)(S3_R), pi^(-1)(S3_R) intersect union A_i)
       ~= (outward-core plumbing end of X minus T,
                                    all curve-end marks). (8)

It is an actual geometric diffeomorphism of marked pairs, not merely
an isomorphism of stable fundamental groups, an equality of determinants,
or an appeal to irreducibility of an arbitrary three-manifold. The
construction works for any chosen compatible completion after common
refinement; it does not infer geometric realizability of that completion.

## 4. The target branched cover and its filling slopes

For a sufficiently large target sphere, the reduced branch divisor B
is smooth near that sphere and transverse to it. Its finite singular
set and the isolated exceptional points of a generic local branch
model are contained in a compact set. Algebraic curve ends have this
eventual transversality; alternatively pull back the sphere radius to
their normalization parameters and use the leading pole terms at their
finitely many points at infinity.
There is no transversality-at-projective-infinity assumption here.

Above a small smooth branch disk the normal finite cover has local
cyclic factors (z,w)->(z^e,w), including e=1 unramified factors.
Therefore pi^(-1)(S3_R)->S3_R is precisely the connected degree-d
branched covering with the ACTUAL restricted monodromy. Removing the
branch-link tubes and filling the lifted positive meridians recovers
the actual closed section in(8). A meridian cycle of length e lifts
the normal filling slope mu^e. All residual sheets remain present.

For arbitrary degrees, preferred-longitude monodromy may permute
meridian cycles; boundary tori are indexed by the full peripheral
orbits. Their longitudinal slopes and residue degrees must be retained.
No identity-longitude theorem is claimed in that generality. When
rho(lambda)=1, each meridian cycle supplies one boundary torus,
with actual rows lift(mu^e), lift(lambda), exactly the PF92 routine.
For the quartic PF60 packet, mu=(12), lambda=1; there are three
boundary tori and the index2 ramification-core mark is the lifted
longitude AFTER the indicated meridian filling. Adding a meridian
multiple does not change that filled core class, but would change
an unfilled framing matrix; the full peripheral data stay explicit.

There is a further degree-free consequence. If B has one component
and its generic meridian is a SINGLE transposition (the sole index2,
residue1 ramified-prime setting), the infinity action is the FULL S_d.
The section is connected, so the action is transitive. A knot group
is normally generated by its meridian: meridian filling gives S3.
Its image is thus generated by conjugate transpositions. Put an edge
between the two letters of each such transposition. Transitivity makes
this graph connected; transpositions along a connected graph generate
S_d, by expressing any transposition using the edges of a joining path.
This proves full infinity monodromy directly from the original finite
map. It does not identify S_d from mere transitivity in a general
branched cover, nor require preferred-longitude identity. With several
branch-link components the same argument works when EACH generic
component meridian is a single transposition. More general ramification
profiles need not have symmetric monodromy.

For general peripheral monodromy there is also an exact integral
normal form, useful beyond the quartic longitude1 case. Let mu,lambda
act by commuting permutations on a sheet set. For each orbit O of
the FULL peripheral group, choose a sheet s. Let e be the length of
its mu cycle and f the length of the induced lambda cycle on mu
cycles. There is a unique a in {0,...,e-1} such that

    mu^a lambda^f fixes s.

The stabilizer lattice in Z*mu plus Z*lambda has basis

    (e,0), (a,f),       |O|=e*f.                    (9)

Indeed any stabilizing word has lambda exponent divisible by f;
subtract the corresponding multiple of(a,f), then its remaining
mu exponent is divisible by e. Commutation makes e constant on O
and the condition independent of the sheet chosen in that orbit.
The normal filling slope is(e,0), and the filled core has longitudinal
degree f, retaining the actual shear a. This is an integral boundary
cover calculation, not a rational homology substitute for periphery.

For example pi(z,w)=(z^e,w^f) has Hopf branch link. Above its first
component there are f meridian cycles of length e, but lambda rotates
them transitively: there is ONE boundary torus with lattice(e,0),(0,f).
Counting only meridian cycles would give the wrong answer. For e=f=2
the whole action is transitive V4, not S4: each branch meridian is a
DOUBLE transposition, so the single-transposition hypothesis above is
essential. These are actual finite map controls, not Keller examples.

Thus the actual finite-map radial end can be compared directly with
the GEOMETRIC whole-source plumbings of PF43/61/77, preserving the
R mark and every needed filling map. The retained JSJ and slope
recognition arguments within those explicit geometric ends still
have their own hypotheses; (8) does not classify graph manifolds
from unmarked homology or dispense with those internal arguments.

## 5. Audit outcome and the remaining scope

Two formerly listed premise checks now have explicit evidence:

| Gate | Evidence and exact outcome |
| --- | --- |
| Individual U union R_i is affine; omitted primes disjoint entire A1 | BGV Corollary11.4 has exactly the original connected normal affine marked-plane hypotheses, without a degree restriction. |
| Exact source volume restricts to each component open | The actual global primitive in(1), retaining both target coordinates. |
| Actual finite-map radial pair equals geometric source plumbing pair | Construction(2)--(8), positive polar divisor on a common boundary resolution, local corner rounding with all marks. |
| Full infinity S_d for single simple branch meridians | Connectedness from Hodge index and positive polar divisor; knot/link meridian normal generation and connected transposition graph. No separate infinity-group epimorphism is needed. |
| Individual affine opens are themselves finite target covers | NOT inferred; the whole original normalization is essential. |
| Every algebraic source chart/plumbing has a finite Keller realization | NOT inferred; the finite map is an input to(2). |

The earlier PF97 logical implication remains visible: PF46's passport,
PF45/46 whole charts, PF50 pencils, PF59 normal bound, PF60 actual
periphery, PF96 genuine-class closure and PF97 target coverage would
give an all-quartic exclusion if every scope implication is valid.
This note verifies component-open use and supplies a direct end
comparison. It does not independently reprove the all-polar key-form
coordinate bridge of PF45, every local passport lemma, the normal
deformation theorem, or the entire later JSJ exclusion chain.

NEXT audit the exact regular-coordinate/key-form bridge and whole-chart
necessity, then consolidate the remaining global passport and exclusion
implications; retain the original forced-F1 finite-map-to-plane-Keller
and marked-normalization requirement. A premise-by-premise audit must
state any actual failed implication precisely. An unexamined inherited
proof is an audit obligation, not evidence of a new family of knots.
Do not reopen closed unit-stage or A4 reset coefficient searches.

The [general peripheral lattice helper](peripheral_boundary_lattices.py)
implements(9) without a longitude-identity assumption. It supplies
boundary components and integral slopes, not the lifted edge-flow rows
needed for general homology or linking. The sealed PF92 helper retains
its identity-longitude condition; no unsupported extension of that
linking routine is inferred.

The [checker](verify_proper_finite_marked_end.py),
[receipt](proper_finite_marked_end_receipt.json), and
[manifest](sources/proper-finite-end/manifest.json) retain the primary
page audits, signed local monomial/blowup identities, finite cyclic
map controls and the properness countercontrol. Finite controls do
not prove the global topology or an algebraic realization. No full
quartic/Keller theorem, arbitrary higher-degree exclusion, independent
whole-chain reproof, external review, formal verification or full-goal
completion is registered here.
