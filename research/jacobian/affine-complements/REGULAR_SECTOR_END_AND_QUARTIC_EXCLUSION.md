# Lens ends and quartic exclusion of the simple regular-coordinate sector

14 November 2026. PF2RULED39 closes the finite quartic realization
problem in the whole simple regular-coordinate sector of PF2RULED38.
All principal parts, source pole orders and target degrees are retained.

For the entire marked surface

    y=u, x=u*t+h(u), h=sum_(j=1)^q c_j*u^(-j), c_q!=0,

the end is a lens space of homology order 2q+1, and the omitted-line
knot has solid-torus exterior. This is proved by an explicit compactification
and the punctured plumbing graph, not just its intersection determinant.
It also makes q an abstract topological invariant in this sector.

A strengthened topological theorem excludes simple S4 covers of any
iterated torus knot with identity preferred longitude, odd nontrivial
lens-space source and a ramification core. Two changes extend the earlier
proof from homology orders divisible by three to every odd order.
The four-sheeted pattern is excluded using the marked knot exterior;
a surviving alternating pattern now descends to its companion instead
of requiring a divisibility-by-three contradiction.

The intrinsic quartic longitude of PF2RULED36 therefore applies.
There is no finite degree-four Keller map on any of these surfaces.
In every finite field degree, their boundary knot at infinity cannot
be a torus knot, and its target-coordinate degrees must be composite.

This does not exclude finite maps of every higher field degree, regular
source coordinates of higher boundary order, the case with no regular
source coordinate, or arbitrary multi-prime normalization packets.
The geometric and topological proofs are written, with explicit primary
imports and exact supporting controls, not formally verified or externally
reviewed.

## 1. A compactification uniform in every principal part

Keep the full actual plane charts from
[PF2RULED38](COORDINATE_POLE_AND_REGULAR_SECTOR.md).
Put H(u)=u^q*h(u), so H(0)=c_q!=0.
Start with P1_x times P1_y, with boundary

    D_x={x=infinity}, D_y={y=infinity}, D_x^2=D_y^2=0,
    D_x.D_y=1.

At the point (x=infinity,y=0), take w=1/x and u=y. The actual
omitted chart gives

    w=u^q/(H(u)+u^(q+1)*t),
    w_0(u)=u^q/H(u),
    w-w_0=-u^(2q+1)*t/[H(u)*(H(u)+u^(q+1)*t)].          (1)

Thus the omitted valuation is produced by 2q+1 successive blowups
along the smooth reference germ w=w_0(u). The first q centers
follow its contact of order q with D_x. After those q blowups
the reference point lies on E_q away from D_x and from the earlier
exceptional curve. The next q+1 centers are free points on the
successive last exceptional curves. Denote all new curves by
E_1,...,E_(2q+1), in creation order.

After the first q blowups the graph is

    D_y -- D_x -- E_q -- E_(q-1) -- ... -- E_1.

The free blowups attach E_(q+1),...,E_(2q+1) starting from E_q.
The last curve Rbar=E_(2q+1) has self-intersection minus one.
The deleted boundary is

    B=D_x+D_y+E_1+...+E_(2q).

Its weights are D_y^2=0, D_x^2=-q, and E_i^2=-2 for 1<=i<=2q.
Its edges, besides D_x--D_y and D_x--E_q, are E_i--E_(i+1)
for 1<=i<2q. Thus E_q is the sole possible branching vertex.

The open S minus B is the actual N_h, not merely a surface with
a matching generic valuation. Both contain the original U; the only
new points are Rbar minus its point on E_(2q). The last blowup
coordinate z=(w-w_0)/u^(2q+1) satisfies

    z=-t/[H(u)*(H(u)+u^(q+1)*t)],
    t=-z*H(u)^2/[1+z*H(u)*u^(q+1)].                    (2)

Near u=0 these are inverse regular coordinates, and on Rbar they
give z=-t/c_q^2. Every finite t is included; the one missing point
is exactly the intersection with B. Away from R the identification
is the identity on the entire source plane. This verifies the
whole marked compactification.

The normal-crossing tubular neighborhoods of B have a plumbed
three-manifold cross-section M_B. Their shrinking collars give
the cofinal end of N_h, up to its irrelevant orientation convention.
The omitted line meets B transversely at E_(2q); its end knot K_R
is the circle fiber there.

## 2. Compute both the closed end and the marked knot exterior

Use Neumann, *A calculus for plumbing*, Proposition2.1 move R6
(splitting), printed pp305--306. Apply it at the pair D_y--D_x:
D_y has weight zero, D_x has genus zero and has just one remaining
branch, attached by a single edge. The R6 count of extra summands
is zero, and D_x has no deleted fiber boundary. The result is the
linear chain

    E_1(-2) -- E_2(-2) -- ... -- E_(2q)(-2).             (3)

By Neumann Theorem6.1(III), its continued fraction is

    [2,...,2]_(2q entries)=(2q+1)/(2q).

Therefore, up to orientation,

    M_B=L(2q+1,2q),  H1(M_B)=Z/(2q+1).                (4)

A matrix determinant alone would not identify the omitted knot.
For its exterior, puncture the graph at E_(2q): in Neumann's
notation put r_(2q)=1. This deletes the circle fiber corresponding
to the actual transverse R. The same R6 move still applies, with
the puncture on the surviving branch; it removes neither that
puncture nor an unrecorded summand.

The surviving punctured chain has one disk-base circle-bundle
piece at E_1 and annulus-base pieces at every other vertex.
The former is a solid torus and the latter are T2 times I.
Plumbing them in a line leaves a solid torus. Consequently

    M_B minus int nu(K_R)=S1 times D2.                 (5)

Thus K_R is a Heegaard core and is primitive in H1(M_B).
The core argument already proves this intrinsic primitivity. In the
linear model the meridian relations also give mu_i=i*mu_1 and
(2q+1)*mu_1=0; its endpoint mu_(2q) is a generator.

The primary source is pinned and the full relevant diagrams were
rendered and visually inspected:
[sources/regular-sector-end/manifest.json](sources/regular-sector-end/manifest.json).
The input is Proposition2.1 R6, not R3 or an inference that equal
determinants identify knots. The punctured version explicitly uses
the boundary convention on printed p303.

For q>=2 the surface cannot even be abstractly isomorphic to the
Wright surface, whose end homology has order three. Moreover

    q=(|H1(end N_h)|-1)/2                              (6)

recovers the minimum source-coordinate pole from the abstract end
in this classified sector. Lower Laurent coefficients change the
positions of blowup centers but not this graph or end pair.

## 3. Strengthen the quartic topological obstruction to every odd order

For an odd integer n>1, let S_n mean a four-sheeted cover of a knot
with transitive S4 monodromy, transposition meridian, identity preferred
longitude, source a lens space with H1=Z/n, and ramification knot
having solid-torus exterior. Let A_n have transitive A4 monodromy,
three-cycle meridian, identity longitude, and the same lens/core
conditions. A core is essential and primitive because n>1.

No iterated torus knot supports either type.

Use the complete literal S4 subgroup and cable-pattern census of
[PF2RULED35](QUARTIC_CABLING_EXCLUSION.md). Its nine S states,
four A states, companion normal-generation conditions and unbounded
peripheral congruences are independent of n. The compression and
incompressibility arguments there also apply to any lens space with
an essential marked knot.

Here p>=2 is cable winding and a is the signed coprime second
parameter. The letter q remains the source pole in Sections 1--2.

### 3.1 The alternating type has a new removal step

The torus case A_n has H1 order 2uv, even, by PF2RULED35's
Seifert computation, so it is excluded for odd n. The unknot
cannot have noncyclic A4 monodromy.

For a genuine cable with 3 dividing p*a, the old fiber-trivial
cases need no condition on n. The three-letter pattern is impossible
by trivial companion monodromy, or has disk base with three genuine
cones p,p,p. Compression makes the partially filled companion a
solid torus; fiber filling kills the essential marked knot and
nonfiber filling cannot be a lens space. The four-letter pattern
has an annular base and all companion gluing tori incompressible.

It remains to treat 3 not dividing p*a. The branched pattern orbit
has three letters; its complement is the unbranched fourth sheet.
The companion image is A4, with identity longitude by the
three-cycle centralizer and A4'=V4. The branched piece B has
Seifert disk base with cones p,3. In particular gcd(p,3)=1.
After capping the fourth sheet the companion piece Y has one
boundary. Boundary-incompressibility of B and the compression
argument in the lens space force Y to be a solid torus. The
compression-ball alternative is excluded by the essential knot
on the B side, exactly as in PF2RULED35.

On the branched boundary the actual covering lattice is

    m0=3*mu_J,  l0=lambda_J.                            (7)

H1(B)=Z without torsion: the two cone orders are coprime.
Its map to H1(V)=Z sends l0 to a generator, so is an isomorphism
and sends m0 to zero. This is precisely the boundary homology map
of the usual index-three branch solid torus W_J.

Replace B by W_J, retaining Y and its actual gluing. Mayer--Vietoris
therefore preserves H1=Z/n. Both Y and W_J are solid tori, so the
new source is a lens space; the companion ramification knot is its
core. The companion cover has A4 monodromy and identity longitude.
It is type A_n over the strictly simpler knot J.

No identification of closed manifolds from homology is made.
Nor is the homological kernel m0 called a compressing meridian
of B: B generally is not a solid torus. Its replacement is allowed
to change the lens-space parameter. One can track the old knot
as [K_R]=p[l0], forcing gcd(p,n)=1, but this does not contradict
odd n on its own.

Every surviving A cable step thus removes one cabling stage while
preserving A_n. Finite descent to its excluded torus or unknot
base proves that no A_n exists, for every odd n>1.

### 3.2 The four-letter symmetric pattern uses the marked exterior

There are no unbranched caps in this case. If the pattern piece B
is boundary-incompressible, its gluing to the unbranched companion
cover injects Z2 into the lens-space group, impossible. The only
remaining patterns have B a solid torus and the marked knot a
regular Seifert fiber. The old order-three argument only rejected
one of them by nonprimitivity when 3 divided n.

Instead remove the marked knot. The exterior in B of a regular
fiber is an annulus-base Seifert piece with at most one cone.
Its two boundaries are incompressible, including the cone-order-one
case T2 times I. The companion cover Y is connected with one
boundary and has incompressible boundary, since it is an unbranched
finite cover of a nontrivial knot exterior. There are no caps
to change that statement. Gluing Y to B minus the marked knot
injects Z2 into the whole marked-knot exterior. But that exterior
is a solid torus. This excludes every four-letter symmetric
pattern, independently of any divisibility of n by three.

### 3.3 All other symmetric cases retain or decrease depth

The three-letter transposition patterns with p>=4 and the
three-letter three-cycle patterns are excluded by the same
marked-exterior incompressibility or companion dual-graph cycle
as in PF2RULED35. Those arguments use an essential marked knot
and a lens source, not n=3.

The remaining three-letter p=2 case uses the exact lattice
m0=3mu_J,l0=lambda_J. B is a solid torus and this m0 is its actual
meridian. Its partially filled companion is a solid torus, so
the usual replacement gives A_n over J with a ramification core.
Section 3.1 excludes that type.

In the two-letter pattern p must be odd. Even a is excluded by
the regular-fiber marked-exterior argument. For odd a, B has disk
base cones p,2, Y is a solid torus, and its boundary lattice is

    m0=2*mu_J,  l0=lambda_J.                            (8)

A fiber filling would kill twice the primitive marked class;
that class is nonzero for every odd n>1. The surviving filling
and exact boundary map in PF2RULED35 replace B by the standard
simple-branch solid torus. This gives S_n over J, preserving
H1=Z/n, a lens source, its ramification core, full S4 and identity
longitude, while decreasing cabling depth. The old knot represents
p[l0], so gcd(p,n)=1; here one must not impose gcd(p,6)=1 when
3 does not divide n. For example p=3 is allowed locally at n=5.

The terminal simple torus case is excluded by PF2RULED34:
its ramification class must be zero, whereas the present core is
primitive and n>1. The unknot has cyclic complement group.
Finite descent proves nonexistence of S_n for every odd n>1.

Two locally consistent controls prevent replacing these descents
by a false pointwise contradiction. At n=5, cable (p,a)=(2,1)
in the new A step and (3,1) in the S step both have fiber slope
m0+6l0 and solid-torus filling slope m0+5l0. The distance is one
and the marked classes are 2 and 3 modulo five. These are local
pieces, not actual global covers; only depth descent excludes
the hypothetical global type.

## 4. Transfer the computed end to the actual finite map

Suppose N_h has an actual finite Keller map pi:N_h->A2, with
source plane U and degree d. Its sole omitted prime is R, and
the volume calculation gives generic index two. PF2OMIT7 and
PF2RULED36 give residue degree one and the entire branch curve
C=pi(R), with finite birational normalization R=A1->C.

The actual large target-sphere inverse image M_pi is a simple
d-sheeted branched cover of S3 along the knot C at infinity,
with full S_d monodromy by the retained infinity-group input.
Both the compactification collars and the radial pair collars
are cofinal neighborhoods of the same marked end. Their stable
groups and marked class are consequently

    pi1(M_pi)=C_(2q+1),
    pi1(M_pi minus K_R)=Z,
    [K_R] primitive in H1(M_pi).                        (9)

For the quartic application one needs a lens/core source, not
merely these homology groups. The end-recognition argument of
PF2RULED35 applies: the imported elliptization theorem recognizes
a closed oriented manifold with finite cyclic fundamental group
as a lens space. The primitive knot is essential, so its exterior
is irreducible; its torus boundary cannot inject into the cyclic
exterior group. The loop theorem and compression make that exterior
a solid torus. This identifies the required conditions for the
actual radial cover, without assuming an arbitrary collar
homeomorphism or using homology as a lens-space recognition theorem.

## 5. Consequences in all degrees and in field degree four

In every possible finite degree d, a torus knot at infinity is
impossible: the simple symmetric finite-H1 torus lemma of
PF2RULED34 forces the ramification class to vanish, contradicting
(9). The unknot is included. Thus no target coordinate presentation
has coprime boundary-normalization degrees. The prime-degree shear
argument of PF2RULED34 also applies: any prime boundary-coordinate
degree could be reduced to a coprime pair by shearing the other
coordinate, or to an excluded coordinate line. Hence every
target-coordinate restriction to R has composite degree at least
four. These are boundary degrees, not source polynomial degrees
or poles. No bound on d or q was used.

Now set d=4. PF2RULED36 derives identity preferred longitude from
the intrinsic residual quadratic cover, without a Wright marking.
The boundary curve has one place at infinity, so its knot is a
finite iterated torus knot by the retained primary input. The
radial cover is S_(2q+1). Section 3 excludes it for every q>=1.
Therefore

    no N_h in PF2RULED38 admits a finite quartic Keller map. (10)

All principal-part coefficients, q, tangency orders and source
degrees are covered. No Groebner or coefficient elimination is
used. In particular the q>=2, c_1=0 locus, which genuinely has
no reduced source-coordinate ruling, is now closed in field
degree four. Its underlying geometric models and exact volume
primitives remain valid controls.

The Euler/collision budget already forces d>=4 in the one-simple-
prime case. Thus any still possible finite Keller map on these
surfaces would have d>=5; existence in those degrees is not claimed.
The pole and normal-form classification stays valid in all degrees.

## 6. Remaining objective and assurance

The quartic finite-map question in the whole simple regular-coordinate
sector is closed. Any actual quartic normalization with exactly one omitted prime of
index two and a regular source coordinate must therefore have boundary
vanishing order at least two after subtracting its boundary value.
PF2RULED36 also excludes a pole-one source coordinate in this quartic
setting. By PF2RULED38 the complementary pole is the intrinsic minimum
when a regular coordinate exists. Thus the remaining regular-coordinate
profiles have both boundary vanishing order and minimum pole at least two.
The possibility of no regular source coordinate remains separate,
as do additional omitted primes and other packets. Higher finite
degrees on these models and on Wright also remain open.

The full original goal is active. The source topology, odd-order
cabling and comparison proofs are written and rely on explicit primary
imports; exact graph, arithmetic and peripheral controls do not
constitute a formal or independent implementation of the three-manifold
proof. No external review or novelty priority is claimed.
