# Isolation and connectivity with one 3-cycle

Let \(N=2r+1\), and fix an \(N\)-cycle \(P\). The mixed polynomial
domain consists of tuples with product \(P\), one 3-cycle, and \(r-1\)
double transpositions (DTs). The 3-cycle can occupy any position. Hurwitz
moves use ordinary composition, with the rightmost permutation applied first:
\[
 Q_i(A,B)=(B,B^{-1}AB),\qquad
 Q_i^{-1}(A,B)=(ABA^{-1},A).                         \tag{1}
\]
The absolute quotient identifies simultaneous sheet relabelings commuting
with \(P\). The exact subgroup for a first-position marked 3-cycle is
\(H=\langle Q_2,\ldots,Q_{r-1},Q_1^2\rangle\) for \(r\ge2\);
at \(r=1\) it is the trivial group.

The main argument isolates the 3-cycle with two exclusive leaves. Deleting
it gives a pure-DT tuple. The existing pure-DT classification then gives
connectivity, with a separate degree-nine base to handle the exceptional
pure degree-seven cores. Double-transposition extraction and arbitrary
two-leaf restoration follow without a circular mixed induction.

## Domain and primary-source comparison

Expanding the 3-cycle as two transpositions, and each DT as its two disjoint
transpositions, gives a tree: its product is transitive and its \(N-1\)
edge occurrences must therefore be distinct and connected. Deleting the
3-cycle leaves three complementary trees, each containing exactly one of
its support labels. These trees do not depend on the chosen two-edge
expansion of that 3-cycle.

Liu--Osserman, [Proposition 1.1, Theorem 1.2 and Section 2](https://arxiv.org/html/math/0609118v1),
require every branch permutation to be a single cycle. A DT fails that
hypothesis. Ripoll, [Definition 0.3 and Theorem 0.4](https://arxiv.org/html/0903.3604v3),
allow only one nonreflection block; here every finite branch factor has
reflection length two. The actual statements were inspected. Neither
theorem supplies the paired DT action used below. Splitting the DTs into
independent transposition strands would change that action. The argument
below supplies its own tree and braid proofs and makes no priority claim.

## Lemma M1. A single contact enlarges either other component

Write the unique 3-cycle first as \(T=(a\,b\,c)\), and let its
complementary components be \(V_a,V_b,V_c\). Suppose a DT \(B\) meets
its support only at \(a\). Then a braid can strictly enlarge either
\(V_b\) or \(V_c\) as a component at the new 3-cycle.
Consequently a representative maximizing the largest complementary-component
size over its entire orbit cannot have such a single contact unless
\(V_a\) is the unique largest component.

**Proof.** A selected complementary factor can be brought to the second
position by inverse moves, conjugating that selected factor by each crossed
DT and leaving crossed factors unchanged. Each crossed DT preserves the
three component sets. Thus this version of gathering retains their sizes,
although it need not retain the selected literal DT.

For a single-contact DT that is to remain literal, instead bring it left
by positive moves. The crossed factors are conjugated by it. They still
preserve each component, so the component sets and both selected literal
permutations remain unchanged. Use this latter gathering here.

Now \(B=(a\,e)(f\,g)\) with six distinct labels. Direct substitution gives
\[
 \begin{split}
 Q_1^2(T,B)&=((e\,b\,c),(a\,c)(f\,g)),\\
 Q_1^{-2}(T,B)&=((a\,e\,c),(b\,e)(f\,g)).
 \end{split}                                                   \tag{2}
\]
Delete \(ae\) from the tree on \(V_a\), obtaining nonempty pieces
\(A_a,A_e\). The first identity replaces the complementary edge \(ae\)
by \(ac\); the new three components are
\(A_e,V_b,V_c\cup A_a\). The second replaces it by \(be\), giving
\(A_a,V_b\cup A_e,V_c\). All other edges, including \(fg\), are unchanged.
Both enlargements are strict. If an original largest component is different
from \(V_a\), choose that enlargement. This also covers ties. \(\square\)

## Lemma M2. Complementary twists synchronize equal component trees

Choose an orbit representative maximizing
\(L=\max(|V_a|,|V_b|,|V_c|)\), and put
\(R=B_1\cdots B_{r-1}\) after front placement of \(T\).
The three restrictions \(R_i\) are full cycles. A full twist on the
complementary strands conjugates their factors by \(R^{-1}\), leaving
\(T\), its component sets and \(L\) unchanged.

At any root whose component is not uniquely largest, a DT wholly within
that component is impossible. If a DT pairs edges in two such components,
their sizes are equal and their edges correspond under the single rooted
cyclic bijection
\[
 f_{ij}(R_i^k a_i)=R_j^k a_j.                                  \tag{3}
\]
If one of the paired components is uniquely largest and the other smaller,
such a factor is also impossible.

**Proof.** Every complementary edge is internal to one of the three
components, and its component edge product is a full cycle by the tree
leaf-insertion argument. The full-twist identity follows by cycling the
complementary entries, as in
[Lemma 9 of the terminal-mass note](TERMINAL_MASS_HURWITZ_DESCENT.md#lemma-9-complementary-twists-develop-the-three-cuts-simultaneously).
Front placement conjugates the 3-cycle within the three component sets and
retains their sizes. Finiteness of the tuple space makes the maximum attained.

For an edge-pair factor define \(E_i\) to be the set of phases at which
the root in component \(i\) hits its support. An internal DT has a
nonempty phase set in its sole component and no other hit. Lemma M1 excludes
it whenever that component is not uniquely largest. For a factor pairing
components \(i,j\), Lemma M1 gives \(E_i\subseteq E_j\) when
\(i\) is not uniquely largest. Their densities are \(2/|V_i|\)
and \(2/|V_j|\), respectively. If \(j\) is larger, this containment
is impossible. Otherwise both containments hold, their sizes are equal,
and their endpoint residue sets coincide. This is precisely (3), with
the same roots and cycles for every factor. All conclusions hold after
every complementary braid: it retains \(T\), the component sets and
the product \(R\). \(\square\)

## Lemma M3. A unique largest component forces isolation

At the maximum in Lemma M2, if one component is uniquely largest, the other
two are singletons.

**Proof.** Denote the smaller sizes by \(s,u\) and the largest by \(\ell\).
No edge in a smaller component can be paired internally or with the largest
component. If either smaller tree is nontrivial, all its edges must therefore
be paired with the other smaller tree. Lemma M2 forces \(s=u>1\) and
identifies the two rooted edge-colored trees by (3).

Choose a nonroot leaf in the first tree and its corresponding leaf in the
second. They are exclusive leaves of one DT. Bring that DT to the second
position by inverse moves among the complementary factors. Crossed factors
fix its two exclusive leaves, so those leaves stay exclusive. The three
component sets and \(T\) stay unchanged. By Lemma M2 the two trees are still
identified by (3). A complementary full twist can put the two nonleaf
endpoints of this DT at the two smaller attachment roots simultaneously.
The two leaf endpoints remain distinct from those roots. Thus the adjacent
pair has the form
\[
                     ((a\,b\,c),(a\,v)(b\,w)),                \tag{4}
\]
where \(v,w\) are exclusive leaves and \(c\) lies in the largest tree.
If necessary interchange the two smaller trees to use this cyclic notation.
The full twist returns factor positions, so adjacency is retained.

Two positive moves give
\[
 ((a\,b\,c),(a\,v)(b\,w))
       \longmapsto ((c\,v\,w),(a\,c)(b\,v)).                 \tag{5}
\]
After deleting the new 3-cycle, one complementary component contains the
entire old largest tree and all \(s-1\) vertices of \(V_a\setminus\{v\}\).
Its size is \(\ell+s-1>\ell\), contrary to maximality. Hence both smaller
trees were singletons. \(\square\)

## Lemma M4. Three equal nontrivial components are impossible

At the maximum in Lemma M2 the three components cannot all have equal size
\(q>1\).

**Proof.** Every DT must pair two distinct component trees, and (3)
identifies all three vertex sets with the same cyclic coordinate set
\(\mathbb Z/q\). Let \(x_{ij}\) count factors of type \(ij\).
Their edge counts give
\[
              x_{12}=x_{13}=x_{23}=(q-1)/2.                    \tag{6}
\]
Thus \(q\) is odd and every type is nonempty. Each component tree projects
bijectively on vertices and edges to a tree on these coordinates; its edges
are the union of the two incident types.

Two projected edges of different types cannot share a vertex. Their types
share one component; identical projected edges would repeat its literal
edge, contrary to the tree property. If they meet at exactly one projected
vertex, their DTs have a single literal common label in that component.
Gather those two factors unchanged while retaining \(T,R\), and conjugate
one by the other. Only one of its two component edges changes its shared
endpoint. Its two projected edges consequently become different, violating
(3). This contradiction is legitimate because maximality and the same
rooted cyclic identifications persist under every complementary braid.

The two nonempty types in any one component tree therefore have disjoint
vertex sets. Their union cannot be a connected tree, contradicting that
component's connectedness. \(\square\)

## Lemma M5. The two-largest tie has a uniform seven-sheet escape

If two largest components have size \(q>1\), and the third is smaller,
then the third is a singleton. For \(q\ge3\) a braid isolates the 3-cycle
with two exclusive leaves.

**Proof.** No component is uniquely largest, so Lemma M2 forbids all internal
DTs and all pairings of unequal component sizes. A nontrivial smaller tree
would have no possible edge partner. It is therefore a singleton \(\{a\}\),
and the other two trees are paired by the same rooted bijection (3).

Hurwitz moves on these paired factors act diagonally as ordinary
transposition Hurwitz moves on either tree. These moves can make that tree
a star at any chosen root. Here is a direct induction proving the needed
statement. Bring a transposition \((b\,v)\) incident to the chosen root
\(b\) first, unchanged, by positive moves. Deleting it splits the remaining
tree into rooted parts at \(b,v\); transpositions in the two parts commute,
so sort them into two blocks. Normalize both blocks recursively to stars
at their roots. Put the \(v\)-block just after \((b\,v)\), then move
\((b\,v)\) right across it by inverse moves:
\[
              ((b\,v),(v\,z))\longmapsto((b\,z),(b\,v)).      \tag{7}
\]
All resulting transpositions contain \(b\). The induction includes singleton
parts. Applying exactly this word diagonally makes the two paired trees
stars at their attachment roots \(b,c\), while \(T=(a\,b\,c)\) stays fixed.

For \(q\ge3\) choose two star satellites
\(D=(b\,x)(c\,y)\), \(E=(b\,z)(c\,w)\). All four satellite
labels are exclusive leaves. Commuting swaps are unnecessary: take the first
two complementary factors after the star normalization. The literal word
\((-1,-1,-1,2,2)\) on the first three entries gives
\[
 \bigl((a\,b\,c),(b\,x)(c\,y),(b\,z)(c\,w)\bigr)
 \longmapsto
 \bigl((a\,b)(c\,y),(x\,y\,z),(b\,y)(c\,w)\bigr).        \tag{8}
\]
Every other factor is unchanged and fixes \(x,z\). They are now exclusive
leaves of the 3-cycle. No condition on the other star satellites is omitted.
\(\square\)

## Theorem M6. Mixed-factor isolation and its exact boundary

Every mixed polynomial DT/3-cycle braid orbit on odd \(N\ge7\) sheets
contains a representative whose 3-cycle has two exclusive leaves. Deleting
those leaves and that factor gives a pure-DT full-cycle tuple on \(N-2\)
sheets. The statement also holds trivially at \(N=3\), where deletion leaves
the empty tuple on one sheet, with identity product. It fails at \(N=5\).

**Proof.** Maximize the largest complementary-component size, equivalently
minimize
\[
                    \rho=N-\max_i|V_i|-2.                     \tag{9}
\]
If \(\rho=0\), two component trees are singleton exclusive leaves of the
3-cycle, as required. A unique largest component forces this by Lemma M3.
Lemma M4 excludes a three-way nontrivial tie. In a two-way largest tie,
Lemma M5 gives sizes \((1,q,q)\); since \(N=2q+1\ge7\), it also gives
isolation, contradicting a positive minimum. These exhaust the cases.
Deleting the isolated factor leaves a tree, hence a full-cycle product.

At degree five the only possible positive maximum pattern is \((1,2,2)\).
More directly, isolation would leave one DT on three sheets, which is
impossible. There are five literal first-marked tuples for a fixed 5-cycle;
the full two-strand orbit has ten tuples, all of pattern \((1,2,2)\).
The exact finite verification is described below. \(\square\)

## Theorem M7. All-degree mixed connectivity, extraction and restoration

For every odd \(N=2r+1\ge3\), the mixed full-cycle domain has one literal
unmarked braid orbit. Its first-marked domain has one orbit under the exact
subgroup \(H\). The same assertions hold in the absolute quotient.

For every \(N\ge5\), some representative has a removable DT with one
exclusive leaf on each transposition. Deletion leaves the same mixed domain
on \(N-2\) sheets. Conversely, adjoining a DT with any two distinct old
anchors and two fresh leaf sheets always reaches the canonical mixed
component in the enlarged degree. This is valid starting at the degree-three
3-cycle. Adjoining an isolated 3-cycle to a pure-DT core, with one old anchor
and two fresh sheets, likewise has no extra attachment or orientation orbit.

**Proof.** At \(N=3\) there is just the tuple \((P)\). Degree five is the
single ten-tuple literal orbit above. At degree seven, Theorem M6 reduces
to the single pure degree-five orbit. Degree nine is the separate complete
mixed base: the retained census, replayed by the checker below, has 78,732
literal unmarked tuples in one orbit, of which 19,683 have their 3-cycle
first. This base is indispensable for this proof because the pure degree-seven
domain has four orbits, including two proper-group orbits.

For \(N\ge11\), Theorem M6 reduces to pure degree \(N-2\ge9\). That domain
has one orbit: the retained complete bases at degrees nine, eleven and
thirteen and
[Theorem 20](TERMINAL_MASS_HURWITZ_DESCENT.md#theorem-20-every-residual-pure-dt-orbit-admits-extraction),
with the existing two-leaf absorption induction, establish all higher
degrees. Move the isolated 3-cycle to the end by positive moves. Crossed
DTs fix its two leaves, so only its old anchor is conjugated. Normalize the
pure prefix while leaving the 3-cycle untouched.

Only its old anchor remains to be normalized. A full twist of the pure
prefix conjugates that prefix by its cycle product \(C^{-1}\). Relabel the
old sheets by \(C\), fixing the two fresh leaves. This restores the prefix
and moves the anchor by \(C\). Its powers act transitively on old sheets.
The two orientations of the 3-cycle are interchanged by swapping its two
fresh leaves, which fixes the entire prefix. Thus all isolated extensions
of the unique pure component form one absolute orbit. Every product-preserving
relabeling of the original tuple is a power of \(P\), and the full twist
realizes it, so absolute connectivity is also literal fixed-product
connectivity. The degree-seven argument uses exactly the same one-anchor
restoration over its unique degree-five pure core.

The unique 3-cycle is a distinguished branch type. A braid taking it from
first position back to first position fixes that strand under the strand
permutation. The
[exact subgroup identification and marking theorem](UNMARKED_POLYNOMIAL_HURWITZ_REDUCTION.md#theorem-4-unmarked-connectivity-suffices-for-both-marked-spaces)
therefore give the marked conclusion in its stated range. The small marked
bases are also checked directly with these exact generators; at \(r=1\)
the singleton domain needs no generator.

For removable DTs, degree five already has its two non-3-cycle labels as
exclusive leaves on different DT edges; placing both on one edge would
disconnect the tree. In degrees seven and nine the exact single-orbit bases
contain a two-leaf DT extension of a lower mixed tuple. For \(N\ge11\),
normalize the pure core to its canonical tuple, which has a removable DT,
then use the anchor transport just proved to attach the isolated 3-cycle
away from that DT's two leaves. It remains removable in the full mixed tuple.
Explicitly, in the tuple \(W_s\) of the
[absorption theorem](UNIFORM_TWO_LEAF_DT_ABSORPTION.md#the-statement-and-metric),
the last-index satellite has exclusive leaves \(a_{s-1},b_{s-1}\) for
\(s\ge3\); the distinguished factor only uses the first satellite's label.
Moving it to the end by positive moves retains its two leaves, and deleting
it leaves a lower mixed full-cycle tree.

Finally, every asserted fresh-leaf extension is itself a mixed full-cycle
tuple: its tree acquires exactly two leaves. The already proved mixed
connectivity normalizes it. This restoration conclusion is a consequence
of isolation and pure connectivity, not a premise in their proof. \(\square\)

## Monodromy, inner components and lifting invariants

For odd \(N\ge7\), every mixed tuple has full \(A_N\). The full-cycle
product makes it transitive. A nontrivial block size in odd degree is at
least three; a DT or 3-cycle moving two such blocks would need to move at
least six sheets. All generators would therefore fix each block, contrary
to their full-cycle product. The group is primitive and contains a 3-cycle
fixing at least three points. Jordan's precise prime-cycle theorem, retained
as [Jones, Theorem 1.1](sources/jones-primitive-cycle-1209.5169v1.pdf), gives
\(A_N\); all generators are even. The degree-five finite orbit has generated
group \(A_5\), and the degree-three tuple has group \(A_3\).

The fixed-product inner and absolute quotients agree: the centralizer
\(C_{S_N}(P)=\langle P\rangle\) consists entirely of even permutations.
If the product is instead allowed to range over both \(A_N\)-conjugacy
classes of long odd cycles and relabeling is restricted to even permutations,
there are exactly two components, distinguished by the product class.
An odd relabeling exchanges them. This distinction is removed by the stated
absolute convention, and is not an additional fixed-product orbit.

Connectivity is established by actual braids and complete small bases, so
any well-defined lifting invariant on the stated tuples must be constant
on the single orbit. It is not assumed away. In the standard Clifford
double cover, a DT lift \(d\) has \(d^2=-1\), and its two signs are conjugate:
for \(D=(12)(34)\), the even centralizer element \((13)(24)\) interchanges
the two orthogonal transposition vectors and sends their product to its
negative. Thus the sign of an arbitrarily chosen DT lift is not a canonical
invariant of these tuples. The checker verifies this calculation exactly.
No assertion about a refined Nielsen class with independently specified
lifted branch elements is made. In particular the exceptional central
extensions in degree seven do not justify omitting its explicit orbit check.

For comparison, the genuine small pure-DT exceptions remain:

| Pure degree | Absolute orbit sizes | Generated groups |
|---:|---|---|
| 3 | empty | no DT exists |
| 5 | 1 | dihedral of order 10 |
| 7 | 21, 21, 7, 7 | two \(A_7\) orbits and two order-168 orbits |
| odd \(N\ge9\) | one orbit | \(A_N\) |

The pure degree-seven orbits with the same group are not silently identified.
The separate mixed degree-nine base is what bridges their isolated 3-cycle
extensions. There are no extra mixed fixed-product exceptional orbits.

## Exact evidence and review boundary

The [checker](verify_mixed_hurwitz_isolation.py) and its
[receipt](mixed_hurwitz_isolation_verification.json) retain deterministic
small-domain enumerations, literal local identities, all 3,102 labeled-tree
orders at paired-tree sizes three through five, larger paired paths through
size twenty in both 3-cycle orientations, 103 full-tree leaf-alignment and
strict descent controls for a unique largest component, and exact lifting controls. All
56 oriented, anchored degree-seven pure-core extensions have literal words
to the common mixed degree-nine seed. Wrong branch classes, repeated edges
and an outside satellite contact are detected by negative controls.
The universal results rest
on the written component, density and braid arguments above, together with
the explicitly named pure bases and absorption theorem. The finite controls
are not an extrapolation to unrestricted degrees.

The [source ledger](mixed_hurwitz_isolation_source_ledger.json) pins the proof,
inputs and current replay. Independent mathematical review, independent replay,
formal verification and historical priority are not asserted. Theorems M6--M7
concern precisely one 3-cycle, all other finite factors DTs, and a full-cycle
product. The [arbitrary-mixture theorem](INDEX_TWO_HURWITZ_MIXTURES.md)
extends connectivity to any number of 3-cycles by a collective maximum,
and also proves fully marked connectivity. Refined lifted domains remain
separate questions.
