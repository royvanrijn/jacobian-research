# Iterated peripheral states and absorbing cabling tails

15 November 2026. PF2RULED82 gives a finite-state reduction for ALL
lengths of actual iterated torus-knot exteriors, and an obstruction
transfer through an arbitrary number of outer cablings. It does not
identify a longer target with its two-stage prefix.

The permutation construction works in every finite degree. In degree4
the necessary tree/planar test has exactly TEN states and29 transition
profiles. Its final simple full-S4 state is absorbing: every further
passing cable has ODD winding and cyclic order2 image. The ACTUAL
filled cable preserves a zero R mark, an order4 mark with self-linking1/2,
positive first Betti number, or an irreducible non-L-space graph end.
Each obstructs the common source. No complementary-source tail theorem
or complete longer-stage exclusion is asserted.

## 1. Exact peripheral transition in any finite degree

Let K have its ACTUAL preferred pair(M,L), and G be the image of its
whole exterior group in S_d. For a genuine positive(r,k)-cable,
r,k>=2, gcd(r,k)=1, choose ru-kv=1. PF80's ordinary cable presentation
and preferred pair apply recursively to a companion knot:

    F=M^k L^r=z^r, beta=M^u L^v,
    mu=beta z^(-v), lambda=mu^(-rk) F.                 (1)

The new image is <G,z>. Include EVERY z with z^r=F. Nothing beyond
(G,M,L) is needed for the next permutation state. Every representation
restricts this way; conversely the root-amalgam presentation extends
a chosen companion representation for every such z. This converse
concerns knot groups, not affine maps.

For E=exp(S_d), the permutation data depend on r,k modulo E.
Modular coprimality gcd(r,k,E)=1 is necessary. Modular Bezout solutions
lift to exact ones for ACTUAL coprime r,k. Two solutions differ by
(kt,rt); beta changes by F^t and z^(-v) by F^(-t), leaving mu unchanged.
Thus this is an exact period reduction, not a bounded-winding inference.

A state is the simultaneous conjugacy class of(G,M,L). There are
finitely many subgroups and commuting pairs at every fixed degree.
Starting states come from torus pairs using PF80's torus presentation.
Genuine splice inequalities and polynomial semigroup conditions remain
EXTRA conditions: accepting a residue path does not realize a curve.

## 2. The whole graph and base genera propagate

Assume the existing lifted decomposition is a forest of planar-base
pieces. Its connected components are exactly G-orbits on sheets.
For the next cable set C=<beta,z>, B=<M,L>. The new incidence graph
has old-component vertices indexed by G-orbits, cable vertices by
C-orbits, and edges by B-orbits, on the SAME based sheet set.
Its new cycle count is

    #B-orbits - #G-orbits - #C-orbits + #<G,z>-orbits.   (2)

Contracting an existing TREE component loses no cycles. Add(2) to
the previous count; checking individual cable edges would be invalid.

On each cable component first quotient sheets by the ACTUAL F-orbits.
The induced cone/boundary permutations come from z,beta,mu; PF69's
Riemann--Hurwitz formula computes the genus. The geometric cone loop
z^(-v) has the same quotient cycles as z, since gcd(v,r)=1 and its
quotient order divides r. F need not be identity. Keep all fixed sheets.

PF62/69 give an integral quotient of final closed-end homology onto
a free abelian group of rank b1(graph)+2 sum(base genera). Earlier cycles or
positive genera cannot be repaired by later attachments or branch
fillings. The incremental test is exact for these two necessary gates.
It does not track integral cone orders, slopes or incompressibility.

## 3. Complete degree4 saturation

For E=12 there are96 admissible ordered residue pairs at each stage.
The starting table has3648 labelled planar torus rows and TEN conjugacy
states. Closing them under ALL transitions(1) passing(2) and the new
genus test produces NO new state. There are1596 accepted transitions
from canonical states,184 root choices rejected by cycles and60
further choices by genus. Grouping accepted transitions by source,
destination and cable-image order gives29 profiles.

| Label | Image G | Meridian M | Longitude L |
| --- | --- | --- | --- |
| 0 | trivial | identity | identity |
| 1 | C2, two fixed sheets | transposition | identity |
| 2 | S3, one fixed sheet | transposition | identity |
| 3 | S4 | transposition | identity |
| 4 | S4 | 4-cycle | identity |
| 5 | S4 | 4-cycle | square of M |
| 6 | C3, one fixed sheet | 3-cycle | identity |
| 7 | A4 | 3-cycle | identity |
| 8 | C2, no fixed sheet | double transposition | identity |
| 9 | C4, transitive | 4-cycle | identity |

Only state3 meets the FINAL quartic conditions G=S4, mu transposition,
lambda=1. All ten can reach it in this necessary finite model. Retain
state5's nontrivial longitude. This saturation is not a depth cutoff.
Counts are canonical-state root choices, not conjugacy classes of
all generator words. The certificate retains all residue transitions.
Every one of PF80's2752 two-stage classes agrees independently:
448 graph cycles,128 further positive-genus failures and2176 accepted
classes. This checks the contraction against the original full graph.

## 4. The final state forces a cyclic odd-winding tail

Suppose state3: M=tau,L=1,G=S4. For r odd, F=tau^k.
If k odd, z^r=tau forces z=tau in S4. For k even, z^r=1
allows z=1 or a3-cycle when3|r. In the latter case every nontrivial
3-cycle gives a forbidden incidence cycle or positive cable genus.
The complete root table checks this for every residue. Therefore

    z=tau^k, mu=tau, lambda=1, C=<tau>.                (3)

For r even, k must be odd, but no even power in S4 is a transposition.
No root exists. EVERY passing continuation thus has r odd, at least3,
and remains state3. This holds however late state3 was first reached.

There are64 admissible odd-r residue pairs and exactly one passing root
for each. Repeating(r,k)=(3,2) any number of times passes these gates:
permutations, tree incidence and genus do not bound characteristic depth.

## 5. Actual tail piece and integral marked homology

Let W be the actual simple degree4 branched end for the current knot,
and R its ramification knot. Delete R to obtain E_R. This fills the
two unramified meridian boundaries of the exterior cover. Under(3),
the degree1 cable pieces perform exactly those same fillings.
The remaining degree2 piece P glues to E_R in the basis m=M^2,L.
Its map is to the actual companion SOLID TORUS, where m dies and
L generates H1; it is not a map to the companion exterior after filling.

For k even, F=1 and P is disk Seifert with cones(r,r).
PF72's integral unit-transfer lemma applies to F=(k/2)m+rL,
since gcd(k/2,r)=1. Thus

    H1(P)=Z*L plus Z/r, m=0, F=rL, [R_new]=rL.        (4)

R is an ordinary fibre. For k odd the generic fibre degree is2;
true branch filling adds cone2, and P has cones(r,2). Since r odd,
its H1 is torsion-free and the map to the solid torus is onto via L.
PF80's local odd-fibre argument, independent of its torus companion,
then gives

    H1(P)=Z*L, m=0, [F_lift]=2rL, [R_new]=rL.         (5)

The last equality uses torsion-freeness, not division in finite H1.
The actual R-meridian distances are1 in(4),2 in(5). Both P are
irreducible and boundary-incompressible, having two genuine cone
orders(r,r) or(r,2), with r>=3.

Apply PF73's integral replacement theorem with torsion retained:

    H1(W_new)=H1(W) plus T,
    [R_new]=(r[R_old],0),
    lambda_new restricted to H1(W)=lambda_old,         (6)

where T=Z/r for k even,0 for k odd. The boundary image in P's
splitting has NO T component. For finite H1 the linking forms split
orthogonally; no oriented formula for the unused T factor is asserted.
In particular R self-linking is multiplied by r^2. This homological
replacement does not require E_R incompressible.

Across any finite tail, [R_final] is the odd winding product times
[R_seed], with zero added-torsion coordinates; self-linking scales by
the product squared. Positive b1 and zero R persist. An order4 mark
with self-linking1/2 retains BOTH its order and its self-linking.
PF67/71 exclude these marks for actual common sources, at all source
orders/depths under their smooth/exact/whole-chart hypotheses. Positive
b1 contradicts finite source-end H1. The shorter prefix need not itself
have any source packet: these are intrinsic end obstructions.

## 6. Irreducible non-L-space seeds persist too

Assume W closed orientable IRREDUCIBLE graph and non-L-space,
and REQUIRE its actual ramification exterior E_R irreducible,
nonsolid, boundary-incompressible and graph. Irreducibility after
drilling a knot is not inferred from irreducibility of W alone.
Its m filling W is non-L-space, so m is outside its L-space INTERIOR.

By(4)-(5), m is P's rational longitude and outside its L-space set.
Rasmussen Corollary6.1 makes the ACTUAL gluing E_R union P non-L-space.
It is irreducible along incompressible tori and remains graph.
Drilling new R replaces P by an annulus piece with cone orders(r,r)
in the even case, or cone r in the odd case. Both have negative
orbifold Euler characteristic and incompressible boundaries. Gluing
to E_R proves the new R exterior irreducible and boundary-incompressible.
Thus all seed hypotheses persist and induction applies.
Némethi Theorem1.5 gives a nontrivial LO group; PF74 excludes the common
source. These topological hypotheses do not follow from the ten states.

For example, for ANY genuine two-stage seed in PF81's strictly
negative EVEN range, the cyclic/A4/C2 sectors have zero R, S4/S4
has order4 and linking1/2, and S3/S3 has an irreducible non-L-space end
except the r3 trefoil b1=1 case. In the non-trefoil S3 cases, drilling
R leaves a pair-of-pants base, with cone s omitted when s=1, attached
to the same two nonsolid companions. All boundaries are incompressible.
The trefoil cap instead gives an annulus with genuine cones(s,J).
Thus the seed exterior hypotheses hold.
EVERY such seed with ANY number
of later genuine cables satisfying the necessary gates is excluded
for the COMMON source. This is an unbounded length consequence.
Complementary sources are not excluded by these common-source gates.

## 7. A many-boundary balance lemma without a degree bound

Let P have n torus boundaries, primitive slopes m_i and integral
relation sum a_i[m_i]=0, all a_i nonzero. Attach actual nonsolid
irreducible boundary-incompressible one-boundary graph pieces A_i,
with matched m_i outside L-space INTERIORS.
For a specified ordering define X_j by attaching A_i for i<j and
filling m_i for i>j, leaving only boundary j. REQUIRE each X_j
irreducible, nonsolid, boundary-incompressible graph with b1=1,
and all used gluing tori incompressible.

In X_1, a_1[m_1]=0, so m_1 is its rational longitude and its filling
is not L-space. Gluing A_1 gives a non-L-space closed C_1.
For j>1, X_j(m_j)=C_(j-1); hence m_j is outside L(X_j).
Gluing A_j gives non-L-space C_j by the same theorem. Induction
closes all n boundaries. No degree, n or source-depth bound is used,
and coefficients need not be units.

PF73's lifted meridian surfaces supply balances, with coefficients
d_i/a_i counting parallel lifts. They do NOT prove the X_j hypotheses.
A solid intermediate cap invalidates this argument and requires a
separate treatment. The tail proof establishes its own hypotheses.

## 8. Frontier and evidence

First arrival at state3 can occur after many earlier states. Retain
loops at states1,2,4,5,6,7,8,9 as actual geometric stages. They cannot
be deleted or replaced by torus companions. Next propagate intrinsic
obstructions through those states, including state5's nontrivial
longitude, and compare longer COMPLEMENTARY target rooted JSJ
with PF79's source endpoint.

Not every prefix has one of the obstructions in(6) or is irreducible
non-L-space. Complete longer-stage quartic exclusion, higher-degree
packet classification, affine extensions and the full goal stay OPEN.
The construction is degree-free; only the ten-state table is d4.

The [checker](verify_iterated_peripheral_tail.py),
[receipt](iterated_peripheral_tail_receipt.json) and
[manifest](sources/iterated-peripheral-tail/manifest.json) seal the
complete finite saturation after the written period proof, parity
transfer and ordinary full-cover controls at three/four stages.
Finite controls support identities; the written topology proves
unbounded transfer. Complete primary audits from PF80/81 are retained.
Rasmussen's related [algebraic-link satellite paper](https://arxiv.org/abs/1703.06874)
was consulted for its recursive perspective; no additional theorem
from it is imported. Deep gluing, orderability and retained topology
are imported, with no formal verification, external review or
novelty/priority claim.

