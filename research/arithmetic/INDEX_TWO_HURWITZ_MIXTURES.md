# Arbitrary index-two Hurwitz mixtures

Fix nonnegative integers \(k,m\), put \(r=k+m\ge1\) and \(N=2r+1\),
and fix an \(N\)-cycle \(P\). Let \(X_{k,m}\) consist of tuples with
product \(P\), exactly \(k\) 3-cycles and \(m\) double transpositions
(DTs), with all placements allowed. An **ordered passport** fixes the set of
positions occupied by each type. Fully marked connectivity means transitivity
on that fixed ordered passport under the pure braid group
\(PB_r=\ker(B_r\to S_r)\). All permutations are composed right to left,
with \(Q_i(A,B)=(B,B^{-1}AB)\).

There is one literal orbit whenever \(k\ge1\), including under \(PB_r\)
for any fixed ordered passport. The only multiple orbits in the entire
index-two tree domain are the four pure-DT degree-seven orbits. The proof
below uses a collective maximum over all 3-cycles, not repeated use of a
uniquely distinguished 3-cycle.

## Literature and exact scope

The all-3-cycle corner is already known. Mühle--Nadeau,
[Proposition 3.1, Theorem 1.2 and Proposition 6.2](sources/muehle-nadeau-alternating-3-cycles-2019.pdf),
give \(\ell_3(g)=(N-\operatorname{ocyc}(g))/2\) and transitivity on
reduced 3-cycle factorizations of an odd full cycle. A DT has 3-cycle length
two, so replacing the \(m\) DTs by 3-cycle pairs gives length \(r+m\),
not the reduced length \(r\). Their theorem does not supply the paired
mixed action. Liu--Osserman,
[Proposition 1.1 and Theorem 1.2](https://arxiv.org/html/math/0609118v1),
already cover genus-zero pure-cycle passports, including this all-3-cycle
corner and higher single-cycle factors. A DT fails their single-cycle
hypothesis. These actual statements and their hypotheses were inspected;
no priority claim for the present elementary argument is made.

The [one-3-cycle theorem](MIXED_HURWITZ_ISOLATION_AND_CONNECTIVITY.md#theorem-m7-all-degree-mixed-connectivity-extraction-and-restoration)
is the induction base for \(m>0\). It includes its genuine degree-five
isolation boundary and its separate degree-nine bridge across the exceptional
pure degree-seven cores. The pure-DT classification is retained separately.

## Lemma A1. Comparing two connected factors

Expanding every factor into its two transpositions gives a tree. Delete a
3-cycle \(T\), with complementary components \(V_1,V_2,V_3\), and set
\(L_T=\max_i|V_i|\). Every other 3-cycle \(C\) lies wholly in one
component, say \(V_i\). Deleting \(C\) leaves a component containing
\(V_j\cup V_l\) and the attachment root of \(T\) in \(V_i\).
Consequently
\[
 L_C\ge N-|V_i|+1;
 \qquad |V_i|\text{ not uniquely largest}\Longrightarrow L_C\ge L_T+2.
                                                               \tag{1}
\]
If \(k\ge2\), a representative and a 3-cycle maximizing \(L_T\)
over **all** representatives and **all** 3-cycle factors therefore have a
unique largest component, containing every other 3-cycle.

**Proof.** The full-cycle product makes the expanded graph connected; its
\(N-1\) edge occurrences must be distinct and form a tree. Each 3-cycle's
two-edge support is connected, so the support of \(C\) is in one component
after deleting \(T\). All edges outside \(V_i\), including the two edges
of \(T\), remain after deleting \(C\). They connect the two other component
sets to the root in \(V_i\). These disjoint sets have total size
\(N-|V_i|+1\). If another component attains \(L_T\), its size, the
third component and that root give at least \(L_T+2\). The tuple space
is finite, so the collective maximum is attained. Equation (1) excludes
every placement of \(C\) outside a unique largest component. \(\square\)

## Lemma A2. DT pruning with other 3-cycles in the large component

At the collective maximum of Lemma A1, the two smaller complementary trees
are singletons.

**Proof.** Front-place the chosen \(T\) by inverse moves. Each crossed
factor is unchanged and preserves its three component sets, so front placement
retains the cut sets. Write \(T=(a\,b\,c)\), with smaller components
\(S,U\) and unique largest component \(M\). By Lemma A1 the smaller trees
contain no 3-cycle factors. All their edges therefore belong to DTs.

The DT single-contact identities remain valid with arbitrary unchanged
outside factors:
\[
 \begin{split}
 Q^2((a\,b\,c),(a\,e)(f\,g))&=((e\,b\,c),(a\,c)(f\,g)),\\
 Q^{-2}((a\,b\,c),(a\,e)(f\,g))&=((a\,e\,c),(b\,e)(f\,g)).
 \end{split}                                                   \tag{2}
\]
Deleting \(ae\) splits its component into nonempty root-side and e-side
pieces. The first identity joins the root side to the component at \(c\);
the second joins the e-side to that at \(b\). A single contact at either
smaller root can thus enlarge the largest component. Gather the selected DT
unchanged by positive moves first; every crossed factor remains within the
same cut components. Uniqueness of \(T\) is not needed.

Let \(R\) be the complementary product. Each restriction to a cut component
is a full cycle, including when that component contains other 3-cycles:
expand those factors into their consecutive two transpositions and use the
tree-product argument. A complementary full twist conjugates by \(R^{-1}\)
and preserves \(T,R\), the component sets and the collective maximum.
For a DT pairing edges in components of sizes \(s,\ell\), absence of a
single contact at the smaller root gives containment of phase sets with
densities \(2/s\le2/\ell\). Thus an edge in a smaller tree cannot be
paired with the larger tree. A DT internal to a smaller tree likewise
develops a forbidden single contact. All smaller-tree edges must therefore
be paired across \(S,U\).

If either is nontrivial, both sizes equal some \(s>1\); their DT edges
correspond under the same rooted cyclic bijection. This follows by applying
the single-contact exclusion at both roots, giving equality of the two phase
sets. The matching identifies the two entire rooted trees. Choose a nonroot
leaf and its matched leaf, which are exclusive to one DT. Bring that DT
second by inverse complementary moves, preserving its two leaves, then use
a complementary twist to align its two anchors with \(a,b\). Their phases
agree under the rooted bijection. The largest tree can contain arbitrarily
many other factors, all of which fix these two leaves.

The resulting pair is \(((a\,b\,c),(a\,v)(b\,w))\), with \(v,w\)
exclusive leaves. Two positive moves give
\[
 ((a\,b\,c),(a\,v)(b\,w))
       \longmapsto ((c\,v\,w),(a\,c)(b\,v)).                 \tag{3}
\]
The new 3-cycle has a complementary component of size
\(|M|+s-1>|M|\). This contradicts the collective maximum. Hence
\(|S|=|U|=1\). \(\square\)

## Theorem A3. Collective 3-cycle isolation

For \(k\ge2\), every orbit in \(X_{k,m}\) contains a 3-cycle with two
exclusive leaves. Its deletion leaves \(X_{k-1,m}\) on \(N-2\) sheets.
Combined with the one-3-cycle result, a tuple with \(k\ge1\) admits such
an isolation if and only if \((k,m)\ne(1,1)\). The \((1,0)\), degree-three
case leaves the empty tuple on one sheet.

**Proof.** Lemma A1 makes the collective maximizer uniquely largest when
there are at least two 3-cycles. Lemma A2 makes its two smaller components
singletons, which are precisely exclusive leaves of the chosen factor.
Deleting that factor and those leaves leaves a tree and full-cycle product.
The remaining 3-cycle count is \(k-1\ge1\). The cases with \(k=1\)
are exactly Theorem M6 of the one-cycle note; its degree-five exception
cannot leave a DT on three sheets. \(\square\)

## Theorem A4. Connectivity for arbitrary mixtures containing a 3-cycle

For every \(k\ge1,m\ge0\), \(X_{k,m}\) has one literal full braid orbit,
and one absolute orbit. An isolated 3-cycle can be restored at any old anchor
with either orientation; no extra attachment orbit is introduced.

**Proof.** First work up to arbitrary simultaneous relabeling, allowing the
product to range over the conjugacy class of \(P\). This is equivalent to
the absolute fixed-product quotient. Induct on \(k\), with the already
proved one-3-cycle theorem as base. For \(k\ge2\), isolate a factor by Theorem A3 and move it last
by positive moves. Every crossed factor fixes its two exclusive leaves, so
only its old anchor changes. Delete it and normalize the lower tuple by
induction, leaving the last factor untouched.

For any normalized lower tuple with product \(C\), a lower full twist
conjugates it by \(C^{-1}\). Relabel old sheets by \(C\), fixing the two
fresh leaves. This restores the entire lower tuple and moves the new anchor
by \(C\), whose powers are transitive. Swapping the two fresh leaves reverses
the new 3-cycle while fixing the lower tuple. Hence all isolated extensions
are absolutely connected. Product-preserving sheet relabelings are powers
of \(P\), realized by the full twist of the entire tuple. Thus the conclusion
is also literal fixed-product connectivity. The argument includes the
degree-three lower base and introduces no mixed two-leaf-DT restoration
assumption. \(\square\)

## Lemma A5. A five-move exchange of two 3-cycle markings

Two 3-cycle supports in this tree domain intersect in at most one sheet.
For five distinct labels,
\[
      Q^5\bigl((a\,b\,c),(a\,d\,e)\bigr)
              =\bigl((a\,b\,c),(a\,d\,e)\bigr).             \tag{4}
\]
This literal stabilizer exchanges the two strand markings. If the factors
are adjacent, every other permutation, and hence the forest obtained by
deleting the two selected factors, is fixed throughout the local loop.
For separated factors, gather them unchanged, apply (4), then undo gathering.
The resulting loop fixes the entire literal tuple and exchanges just those
two strand markings. Gathering need not fix every outside edge at its
intermediate steps; its endpoint loop fixes them all.

**Proof.** Two connected three-vertex subtrees sharing two vertices would
repeat an edge or form a cycle, both impossible. Equation (4) is a literal
permutation identity on five labels; its endpoint substitutions are retained
in the checker. Each \(Q\) exchanges its two strands, so the fifth power
has odd strand permutation. Literal gathering moves the right selected
factor left by positive moves without crossing the left selected factor.
Conjugating the local stabilizer by this gathering word gives exactly the
transposition of the original selected strand labels. For disjoint 3-cycles
the local period is two, so there is no odd local literal stabilizer; global
normalization below puts the 3-cycles at a common anchor instead. \(\square\)

## Theorem A6. Fully marked connectivity

For \(k\ge1\), the pure braid group is transitive on every fixed ordered
passport in \(X_{k,m}\). In particular any two marked 3-cycle factors can
have their relative strand markings normalized without changing the
underlying literal tuple at the end of a stabilizer loop.

**Proof.** A convenient tuple has every 3-cycle at the same anchor 0. For
\(m=0\) use \((0\,1\,2),(0\,3\,4),\ldots\). For \(m=1\), use the
degree-five core \(((0\,1\,2),(0\,3)(1\,4))\) and append the remaining
3-cycle buds at 0. For \(m\ge2\), use the pure tuple \(W_m\) of the
[absorption theorem](UNIFORM_TWO_LEAF_DT_ABSORPTION.md#the-statement-and-metric)
and append all \(k\) fresh 3-cycle buds at 0. These are tree tuples in the
same component by Theorem A4.

At this tuple, Lemma A5 supplies a stabilizer transposition between any two
3-cycle strands, independently of the DT strands. For the DTs of \(W_m\),
the distinguished factor \(A\) and first-index satellite \(D_1\) have
product order five; each pair \(D_1,D_j\), \(j\ge2\), has product order
three. For DT involutions \(D,E\), the identity
\[
 Q^s(D,E)=(D(DE)^s,D(DE)^{s+1})                               \tag{5}
\]
shows that each odd order gives a literal stabilizer transposing those two
strands. Gather and undo as in Lemma A5. These transpositions form a connected
graph on all DT strands, generating \(S_m\). The cases \(m\le1\) need no
DT transposition. The literal stabilizer thus maps onto
\(S_k\times S_m\), the full group of type-preserving strand permutations.

Conjugating these loops by a connecting braid transports this full
type-preserving stabilizer image to any tuple. Given two tuples with the
same ordered passport, Theorem A4 supplies a braid between them. Its strand
permutation preserves the two type sets. Multiply by a stabilizer at the
starting tuple with the inverse type permutations;
the adjusted braid has identity strand permutation and the same endpoint.
This is an element of \(PB_r\). \(\square\)

## Theorem A7. Removable factors and complete index-two exceptions

If \(k\ge1,m\ge1\), every orbit has a removable DT; deletion leaves
\(X_{k,m-1}\). Every extension of a tuple with \(k\ge1\) by a DT at any
two distinct old anchors and two fresh leaves restores to the unique enlarged
component. The corresponding one-anchor 3-cycle restoration is Theorem A4.

**Proof.** Normalize a one-3-cycle core with the same DT count to a removable
DT representative, using Theorem M7. Add the remaining 3-cycle buds one at
a time at an old anchor outside that DT's two leaves. Anchor transport in
Theorem A4 permits this while restoring the lower literal tuple; its two
leaves stay exclusive. Connectivity then gives extraction in every orbit.
Positive moves bring the removable DT last without moving its leaves in
the crossed factors. Deletion leaves a lower full-cycle tree. Every asserted
fresh-leaf extension is a valid tree tuple in the already connected enlarged
domain, proving restoration without circularity. \(\square\)

For completeness the pure-DT results and genuine exceptions are:

| Domain | Literal or absolute full-braid orbits | Removable factor |
|---|---|---|
| \(k=0,m=1,N=3\) | empty | none |
| \(k=0,m=2,N=5\) | one, group \(D_5\) of order 10 | no removable DT |
| \(k=0,m=3,N=7\) | two \(A_7\) and two order-168 orbits; absolute sizes \(21,21,7,7\) | DT in both \(A_7\) orbits; none in the order-168 orbits |
| \(k=0,m\ge4\) | one, group \(A_N\) | removable DT |
| \(k\ge1\) | one, group \(A_N\) (including \(A_3\)) | T as in A3; DT whenever \(m\ge1\) |

There is no extra pure-braid splitting of the pure-DT rows. In degrees five
and seven this is checked on their complete finite domains; above degree
seven the same odd DT stabilizers in \(W_m\) give the assertion. Pure-DT
restoration from degree five to seven retains its two genuine full-\(A_7\)
attachment-distance orbits from the absorption theorem. A 3-cycle changes
this situation: the separate mixed degree-nine base already joins extensions
of all four pure seven-sheet cores. With no factors, the one-sheet empty
tuple is a further trivial case.

## Proposition A8. Monodromy, inner classes and spin

If \(k\ge1\), the group is \(A_3\) at degree three and \(A_N\) for
all odd \(N\ge5\). At \(N\ge7\), odd-degree transitivity and support
at most four imply primitivity: a generator moving two blocks of size at
least three would move at least six sheets. Jordan's 3-cycle theorem then
applies. Degree five has just the \((k,m)=(1,1),(2,0)\) cases, whose exact
generated-group closures have order 60. This accounts for monodromy without
assuming it in the isolation argument.

Fixed-product inner and absolute quotients agree because
\(C_{S_N}(P)=\langle P\rangle\subset A_N\). If \(P\) is allowed to
range over both alternating conjugacy classes and only even sheet relabelings
are allowed, there are two product-class components. They are exchanged by
odd relabeling; they are not extra fixed-product braid orbits.

When \(m>0\), a DT's lifts in the standard double cover have order four,
and their two signs are conjugate by an even centralizer element, as checked
in the one-cycle note. Thus an arbitrary DT lift sign is not canonical.
When \(m=0\), canonical odd-order lifts do exist. With the inverse cycle
at infinity included, their spin product is the single forced value
\[
 s=\widehat T_1\cdots\widehat T_r\widehat{P^{-1}}
       =(-1)^{\binom r2}.                                    \tag{6}
\]
Indeed the tuple has genus zero: its finite indices sum to \(N-1\) and
the infinity index is \(N-1\). Fried's
[Invariance Corollary 2.3](sources/fried-lifting-invariants-math-0611591v4.pdf)
gives exponent \(r+(N^2-1)/8\), which differs from \(\binom r2\) by
\(2r\). The actual statement and its proof were inspected. Exact Clifford
controls check the lift convention for small \(r\). Wrong prescribed spin
products therefore have no lift of this original pure-cycle tuple; they do
not define a second orbit of the original domain. Other independently
specified lifted branch elements are a refined domain, not silently identified
here. Connectivity proves that every invariant well defined on the original
fixed-passport tuples is constant on the stated single orbit.

## Lemma A9. A higher-index dihedral family

For every odd \(n\ge5\), on \(\mathbb Z/n\) put
\[
                    A(x)=-x,\qquad B(x)=-x-1.                 \tag{7}
\]
Each has type \(2^{(n-1)/2}1\), their product is the full cycle
\(x\mapsto x+1\), and their indices sum to \(n-1\). Their expanded
graph is a tree, but their group is the dihedral group of order \(2n\),
proper in \(S_n\). It is contained in \(A_n\) when \(n\equiv1\pmod4\).
For composite \(n\) it is imprimitive. Nevertheless this fixed passport
has one literal \(B_2\)-orbit of size \(n\), hence one absolute orbit.

**Proof.** Both permutations are involutions with one fixed point. They
invert their product \(P\), giving the dihedral group. All factorizations
of this product by two involutions satisfy \(APA=P^{-1}\), so
\(A(x)=c-x\); the \(n\) choices of \(c\) are exactly these reflections.
Hurwitz moves cycle them. For any proper divisor \(d\) of \(n\), residue
classes modulo \(d\) are nontrivial blocks. Parity is
\((-1)^{(n-1)/2}\). Connectedness of the edge graph follows from the product,
and its \(n-1\) distinct edge occurrences make a tree. \(\square\)

This gives proper monodromy at arbitrarily large higher indices, including
proper primitive dihedral groups in prime degree: a transitive action of
prime degree is primitive. It prevents extending the index-two small-support
argument by degree alone.

## Lemma A10. A higher-index imprimitive family containing a 3-cycle

For odd \(b\ge3\), set \(n=3b\) and
\[
 T=(0\ b\ 2b),\qquad
 B=(0\ 1\ \cdots\ b-1)(b\ b+1\ \cdots\ 2b-1)
                   (2b\ 2b+1\ \cdots\ 3b-1).                 \tag{8}
\]
Then \(TB=(0\ 1\ \cdots\ 3b-1)\), their indices sum to \(n-1\),
and their expanded graph is a tree. Both factors are even, yet the group
is \(C_3\wr C_b=C_3^b\rtimes C_b\) of order \(3^b b\), with the
nontrivial blocks \(\{i,i+b,i+2b\}\), \(0\le i<b\).

**Proof.** Track the three ends of the displayed b-cycles to check the
full-cycle product. The conjugates \(B^iTB^{-i}\), \(0\le i<b\), are
disjoint 3-cycles on exactly the displayed blocks. They generate \(C_3^b\),
and \(B\) cyclically permutes those blocks with order \(b\). This proves
the claimed semidirect product. The two edges of \(T\) connect the three
disjoint b-vertex paths in the expansion of \(B\), giving a tree. \(\square\)

Thus even a distinguished 3-cycle does not force alternating monodromy
after disconnected higher-index factors are admitted.

## Lemma A11. A genuine lifting split outside the minimal tree domain

On five sheets, factorizations of a fixed 5-cycle into **three** 3-cycles
have 125 literal tuples in exactly two full braid orbits, of sizes 75 and
50. Both groups are \(A_5\), but their standard spin products are \(+1\)
and \(-1\), respectively. Their absolute orbit sizes are 15 and 10.

**Proof.** The retained checker directly enumerates all ordered pairs of
the 20 3-cycles, completes the third from the product equation, and checks
full signed braid closure. Exact group closure and rational Clifford lift
products give the stated two values. A same-order odd lift is conjugation
equivariant, so its product is preserved by Hurwitz moves. The fixed-product
centralizer of order five acts freely, since a nonidentity 5-cycle power
cannot preserve a three-element support. This proves the finite assertion.
The exact inputs and the two orbit digests are retained in the receipt.
\(\square\)

This example has finite index sum six, exceeding \(N-1=4\), and genus
one after the full cycle at infinity is included. It is deliberately outside
the present hypothesis \(N=2r+1\). A full-cycle product alone is insufficient
to assert the tree representation or the connectivity theorem.

Single-cycle higher-index
genus-zero passports are already covered by Liu--Osserman. Ripoll's primitive
block theorem covers one nonreflection block with all remaining factors
reflections, not arbitrary collections of disconnected higher-index factors.
For those collections, the remaining research question is to identify a
component comparison and restoration theorem appropriate to their number
of cut components, together with actual block systems and reduced lifting
invariants. No unrestricted higher-index connectivity theorem is asserted.

## Evidence and review scope

The [checker](verify_index_two_hurwitz_mixtures.py),
[receipt](index_two_hurwitz_mixtures_verification.json) and
[source ledger](index_two_hurwitz_mixtures_source_ledger.json) separate the
universal written proofs from finite controls. They retain collective-cut
comparisons, exact literal marking loops, complete small fixed-product
mixture and pure-braid domains, spin normalization checks and the dihedral
and imprimitive-family controls, and the complete nonminimal lifting split
of Lemma A11. The unrestricted theorem is not an extrapolation of those
finite domains. Independent mathematical review, independent replay, formal
verification and priority are not asserted.
