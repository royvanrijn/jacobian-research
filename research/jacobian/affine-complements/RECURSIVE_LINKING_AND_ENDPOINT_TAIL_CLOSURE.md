# Recursive peripheral linking and endpoint tail exclusion

15 November 2026. PF2RULED92 proves an all-degree peripheral
reconstruction lemma, the all-length actual-D3 linking recursion,
and an endpoint version of the tail gluing obstruction. It closes
the exceptional common S3 L-space seed tails left by PF91 and
further restricts the remaining solid histories. It does not close
all later S3/A4 tails or the full goal.

## 1. A degree-free unfilled peripheral reconstruction

Let E be a connected finite cover of an oriented knot exterior,
with preferred longitude having IDENTITY monodromy. Suppose the
meridian has cycles of lengths a_1,...,a_c. Each gives a boundary
torus with actual positive normal m_i=lift(mu^a_i) and preferred
lift L_i=lift(lambda). Fill ALL m_i to obtain W, with cores K_i
oriented along their images in the original knot. REQUIRE every
[K_i]=0 in INTEGRAL H1(W), not merely rational nullity.

The homology sequence for drilling these null-homologous curves,
and intersection with their integral Seifert surfaces, give

    H1(E)=direct-sum_i Z*m_i plus H1(W).

The meridians are primitive and independent. Each L_i maps to
[K_i]=0 after filling, so it is an INTEGER combination of m_j,
with NO added-module coordinate. Write L_i=sum F_ij*m_j.
For i!=j, F_ij is the mutual linking lk(K_i,K_j); hence F is
symmetric, including in W with positive first Betti number.
The diagonal entries keep the ACTUAL lifted preferred framing.

Lift a Seifert surface of the original knot. Its boundary contains
a_i parallel copies of L_i on boundary i. Consequently

    sum_i a_i*L_i=0,   F*(a_1,...,a_c)^T=0.             (1)

For TWO components, set g=gcd(a,b), a'=a/g,b'=b/g.
Symmetry, integrality and(1) force, for a unique integer t,

    F=t*[[-b'^2, a'b'],[a'b',-a'^2]].                  (2)

In particular mutual linking is divisible by a'b'. To prove this,
write the off-diagonal entry ell. The two equations imply
a' divides ell and b' divides ell; these are coprime, so
ell=a'b't, giving both diagonal entries. This is valid at EVERY
degree under the stated peripheral and integral-null-core conditions.
The parameter t may be0; then F=0. Thus the resulting matrix has
rank at most one, with the displayed primitive rank-one generator.
It asserts neither those hypotheses nor formula(2) for an arbitrary
S4/A4 cover. With more components, the symmetric matrix/kernel
condition(1) is retained; one coefficient need not determine it.

## 2. Every retained actual D3 history has the square recursion

Take PF90's ACTUAL connected three-sheet S3 companion cover,
with meridian cycles(2,1). Its partial fixed-meridian cap Esharp
has primitive m2 and lifted L2=-A*m2 with zero added-module
coordinates; PF90 proves its index2 filling core is integrally zero
at every retained length. Filling both meridians in the lifted
surface balance2L2+L1=0 also makes K1 integrally zero.
Apply Section1 with a=2,b=1. Capping m1 fixes t=A, so

    L2=-A*m2+2A*m1,
    L1=2A*m2-4A*m1,
    lk(K2,K1)=2A.                                    (3)

Initially T(p,q) with S3 image has A=pq/6 by PF73. If the first
S3 image is instead born from cyclic C2 or C3 at(r,k), PF65's
actual peripheral cyclic-companion replacement works BEFORE either
new boundary filling. It replaces the older pieces by the unknotted
torus pattern, adding unmarked modules. Thus A_birth=rk/6; the
coprime exponents are divisible by2 and3 in opposite order.
Every subsequent state2 cyclic2 loop of odd winding w uses PF90's
UNFILLED transfer and gives

    A_new=w^2*A,   ell_new=w^2*ell.                    (4)

Equations(3)--(4) hold for every finite retained history, including
NONsolid companions, even cable parameters, all allowed added
homology modules and all unbounded parameters. They concern the
specified COLORING/monodromy history. They are not a single invariant
value attached to an uncolored knot, and the full four-sheet end
has not been called dihedral. No novelty/priority claim is made.

## 3. Independent computation and comparison with the authors' code

The new [composed peripheral helper](composed_peripheral_cover.py)
returns ALL unfilled normal-meridian and preferred-longitude rows
for connected finite covers with longitude monodromy1. It tracks
permutations and lifted edge flows and uses binary powers. It never
expands recursive words and never substitutes a meridian filling
for an exterior. Both ordinary and reversed maximal trees are checked.
An exact Hermite basis for the SAME row lattice, followed by unit
relation cancellation, precedes Smith reduction. Small independent
Smith controls check the integral group, including free summands.
Two generic-Smith diagnostics were stopped on rapid coefficient/
memory growth; their traces and terminal records remain in
`research/.scratch/wright92/smith_stop.json`. They are not mathematical evidence.
The completed replay, including a sixteen-loop cover and a
ten-level large-coefficient solid chain, passes under a512MiB limit.
The exact checker verifies both equations(3) integrally, the whole
closed homology module and the weighted balance at every boundary.

The [external comparison](compare_dihedral_linking.py) executes
UNMODIFIED dihedrallinking.py from the authors' repository at the
commit pinned in the manifest. It constructs actual closed-braid
diagrams, solves Fox colorings and records the full input arrays.
For a trefoil(r,k) cable the braid is the r-parallel of sigma1^3,
followed by residual twist k-3r; this corrects its blackboard writhe.
One positive Reidemeister-I curl makes the crossing count even
when needed by the external code's configuration convention.

For cyclic2 loops, top colors are constant within each r-strand
bundle and different between bundles. This selects the old S3
companion and cyclic cable restriction BEFORE computing linking.
The authors' algorithm then gives18 for EACH of the actual
trefoil cables(3,7),(3,5),(3,2), agreeing with(4). Five independent
torus examples and the C3 birth(2,3) also agree.
Additional degree5/7 torus-cover comparisons check the multi-component
weighted-kernel calculation; only off-diagonal entries are compared.

The SAME knot C_(3,2)(trefoil), with a DIFFERENT coloring selected
by equal old bundle-meridian reflections, has first S3 BIRTH(3,2)
and linking2. This independently agrees with A_birth=1 rather than
the loop A=9. The initially unmatched coloring is retained as a
countercontrol; the coloring is never selected by its linking output.
Only off-diagonal linking is compared. The authors' chosen diagonal
push-off framings are not identified with our preferred lifts.

This uses the concrete framework of Cahn et al.2112.14790v2
Sections3--6 and the actual-D3 scope of Cahn--Kjuchukova1611.10330v2.
No ribbon, slice, Xi or cobordism hypothesis is inferred.

## 4. A missing-INTERIOR seed suffices for all cyclic2 tails

Let E_R be an irreducible, nonsolid, boundary-incompressible graph
manifold, b1(E_R)=1, with actual basis(m,L). Suppose

    m is OUTSIDE the INTERIOR of its L-space interval. (5)

Its m filling may be an L-space ENDPOINT. Attach an ACTUAL cyclic2
tail with odd winding w>=3, coprime cable parameter k, and all
genuine fixed-sheet caps, as in PF82. Its filled pattern P has
disk cones(w,w) if k even, or(w,2) if k odd. Both are nonsolid,
irreducible, boundary-incompressible graph pieces. The integral
peripheral calculation makes m the rational longitude of P, hence
m is outside the INTERIOR of P's L-space interval as well.

The retained published graph gluing theorem applies to this MATCHED
slope and makes E_R union P an irreducible non-L-space graph manifold.
If m was the rational longitude of E_R, its filling already has
positive b1 and PF82's homological transfer preserves that obstruction.
Otherwise the actual marked Mayer--Vietoris transfer keeps the new
closed H1 finite (adding only Z/w or0), so the new exterior has b1=1.
Drilling its new R leaves an annulus with cones(w,w), or cone w,
attached to E_R. This is again irreducible/incompressible graph,
and its normal meridian filling is the non-L-space just obtained.
Thus(5) holds at the next step. Induction excludes EVERY nonempty
finite tail of such a seed as a common source, using PF74.

This strengthens PF82 by replacing its sufficient non-L-space
seed filling with(5). It is independent of the field degree whenever
the actual tail has that cyclic2 active orbit and genuine fixed caps.
It is NOT a claim about all higher-degree tail representations.

## 5. The exceptional L-space seed is now closed through ALL tails

PF91's exception is the solid S3 all-plus history, A=W^2,
c=2g(J)-1=2A-1, arrival r=3,k=6A-2. Its R exterior is the ACTUAL
ordinary E(J), and its normal R meridian is c*n+L1. This is the
INCLUDED genus endpoint of the positive L-space interval[c,infinity],
so it satisfies(5). Its boundary and graph hypotheses hold.

Section4 therefore excludes EVERY nonempty later simple tail.
PF91 already excludes the final arrival by its marked linking.
Its H1Z17/linking8/17 tail countercontrol is retained: it passes
that particular linking test, while the endpoint/gluing argument
excludes it. No false invariance of linking under odd squares is used.

Consequently EVERY incoming S3 COMMON history with arrival eta<=0
is now excluded at arbitrary companion depth and later tail length.
The prefix is never assumed to have its own affine source.

## 6. Further all-depth reduction of the positive-arrival frontier

PF91's two-denominator argument for r>=9 actually requires only

    l0=(2A-1)*r-k>0,                                 (6)

not eta<0. It implies Jcap=Ar-k/2>(r-1)/2>1; the same exact
floors/ceilings put the ordinary meridian outside Y's L-space
INTERIOR, so the closed seed is non-L-space and all later tails
are excluded. No arrival eta sign is needed for this extension.

After a LAST MINUS solid loop of odd winding w, PF91 gives
B=2A-w, D=2A-c=2w-1. Every genuine next arrival satisfies
k<r*B<r*(2A-1), so(6) ALWAYS holds at r>=9. In the complementary
pencil the same inequality gives Jcap>rw/2>1; with any later tail
the genuine two-cone annulus becomes internal and is forbidden.
Final arrivals were already excluded by PF91. Thus EVERY surviving
LAST-MINUS S3 history in EITHER pencil must arrive with r=3.

For r=3 the seed filling is S3_t(J), t=k/2-A. Whenever t<=c,
that slope lies outside the INTERIOR of the ORDINARY companion's
positive L-space interval. If t<c, PF91's non-L-space/positive-b1
argument and actual R-exterior checks apply. If t=c and D=1,
Section5 applies. At t=c,D>1 the R exterior still has an extra
annulus cone D; the ordinary endpoint is NOT silently transported
to its R boundary, and that case remains open.
For a last-minus history, the strict infinity bound
k<3B, with k even and B odd, gives

    t-c <= (w-3)/2.

Therefore a surviving COMMON last-minus history has odd w>=3 and

    r=3, 0<=t-c<=(w-3)/2,                            (7)

and at least one later tail. All-plus common survivors instead
have POSITIVE arrival eta, outside all prior seed exclusions.

For complementary histories at r>=9, PF91 already forces absolute
cap order1 if any later tail follows. The minus case was just excluded,
so only ALL-PLUS histories remain there, with EXACTLY

    k=2Ar +/-2,
    eta_arrival=r+1 or r-3.                           (8)

These are necessary possibilities, not realized affine data.
Complementary r=3 tails and PF89's later A4 ordered terminal pairs
(2,4),(2,8),(4,10) retain their prior restrictions and remain OPEN.

## Evidence and limits

The [checker](verify_recursive_linking_endpoints.py),
[receipt](recursive_linking_endpoints_receipt.json), and
[manifest](sources/recursive-linking-endpoints/manifest.json)
retain the all-length written proofs, all unfilled integer boundary
maps, independent external code and coloring inputs, and exact
endpoint/frontier identities. Deep graph gluing/LO and peripheral
rigidity are retained imports. The external code is a comparison,
not a proof of our universal statements or an independent reproof
of its complete algorithm. No full later-S3/A4, longer-stage,
quartic, higher-degree, degenerate or finite Keller-map exclusion,
formal verification, external review or novelty claim is made.
