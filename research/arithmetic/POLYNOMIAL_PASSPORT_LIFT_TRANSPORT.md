# Return tuples and component transport for polynomial passports

The [component problem](POLYNOMIAL_PASSPORT_COMPONENTS.md) fixes a full
cycle \(P_N=(0\,1\,\ldots\,N-1)\), finite branch factors
\(g=(g_1,\ldots,g_r)\), and
\[
                 \prod_i g_i=P_N,\qquad
                 \sum_i\operatorname{ind}(g_i)=N-1.           \tag{1}
\]
Fix a preserved congruence partition with \(d\) blocks of size \(b\),
so \(N=bd\) and block \(a\) is \(\{a+dt:0\le t<b\}\).
Its quotient tuple \(q=(q_1,\ldots,q_r)\) has product \(P_d\).
Identity quotient factors are retained as branch positions.

The return-budget theorem H3 gives quotient index \(d-1\) and total
return index \(b-1\). The results below put those returns in a single
literal degree-\(b\) tuple, determine the induced braid transport, and
prove sufficiency of the decorated quotient component when \(b=2,3\).
An explicit degree-eight example shows that this sufficiency fails at
\(b=4\), even when the separate inner and outer components are unique.
The arbitrary-passport classification remains open.

## Primary comparison and conventions

Jones--Zvonkin [Section 5.4, retained pp155--159](sources/jones-zvonkin-braid-cacti-2002.pdf)
build imprimitive cacti from copies of a quotient graph and attach
remaining edges in specified blocks. Their exact example makes the
choice of attachment data essential. The present argument gives a
general tree normalization and transport rule for arbitrary cycle types.
No historical novelty is asserted.

Müller's [Appendix, decomposition setup](sources/mueller-primitive-polynomial-monodromy-author.pdf)
relates polynomial decompositions to intermediate groups. Separate
decomposition factors do not specify a braid action on attachments.
König--Neftin--Rosenberg,
[Theorem 1.1](https://neftin.net.technion.ac.il/files/2024/02/version-2.pdf),
prove a monodromy largeness statement for compositions over \(\mathbb Q\)
with degrees at least five and alternating/symmetric constituent groups.
That theorem is not a component theorem and does not supply the
small-degree or cyclic-fibre transport here. The actual primary statements
and these hypothesis differences were inspected.

Composition is right to left and
\(Q_i(A,B)=(B,B^{-1}AB)\). Full braids allow all placements of the
prescribed original cycle types. A fixed ordered passport can instead
be included as a tuple of tags in the quotient action; no unrestricted
marked transitivity is presumed.

## Lemma L1. The canonical quotient tree and its fibre coordinates

Let \(q_1\cdots q_r=P_d\) and
\(\sum\operatorname{ind}(q_i)=d-1\). There are
\[
             M=\sum_i\#\operatorname{cycles}(q_i)=(r-1)d+1    \tag{2}
\]
quotient cycle slots, including fixed points. For each nontrivial cycle,
its cycle-incidence tree has a unique root-side block, with root block
zero. Start its cyclic order at that block and delete its final directed
edge returning to it. The remaining directed edges form a tree on the
\(d\) blocks. A fixed-point cycle has just its deleted loop.

Write a lifted permutation as
\[
       (a,t)\longmapsto(q_i(a),A_{i,a}(t)),\qquad A_{i,a}\in S_b.
                                                                  \tag{3}
\]
There is a unique choice of fibre charts \(F_a\), with \(F_0=1\), that
makes every retained tree-edge transport the identity. The residual map
at a deleted edge of a slot \(O\) is
\[
             h_{i,O}=F_{q_i(a)}^{-1}A_{i,a}F_a               \tag{4}
\]
where \(a\) is its source. Its cycle type is the return type of H3.
It is a permutation of one common labeled \(b\)-element set.

**Proof.** The quotient cycle-incidence graph is a tree by H1, after
identity factors and fixed-point leaves are omitted. Rooting it fixes
the indicated parent block of every nontrivial cycle. Replacing its
cycle vertex by the directed path obtained after deleting the last
edge preserves connectedness and gives exactly \(d-1\) edges.
On a retained edge \(a\to z=q_i(a)\) define recursively
\(F_z=A_{i,a}F_a\). The tree makes this well defined and unique.
Changing all old fibre labels changes the normalized maps by the same
root conjugation. Formula (4) is the normalized map at the only
remaining edge of its quotient cycle. Going once around that cycle
proves the return-type assertion. \(\square\)

## Theorem L2. A literal degree-reducing return-tuple bijection

Starting at block zero, apply the quotient word \(q_r,\ldots,q_1\)
\(d\) times, following the current block. This boundary walk uses every
directed edge \((i,a)\) exactly once. List its deleted edges in reverse
chronological order and place the corresponding \(h_{i,O}\) in that order.
The resulting tuple \(R_q(g)\) has length \(M\) and satisfies
\[
              \prod R_q(g)=P_b,\qquad
              \sum_{h\in R_q(g)}\operatorname{ind}(h)=b-1.    \tag{5}
\]

For a fixed literal quotient \(q\), this is a bijection between literal
lifts with fixed product \(P_N\) and degree-\(b\) tuples (5), subject only
to the following local original-passport constraint: a slot of quotient
length \(\ell\), with return partition \(\nu\), contributes the parts
\(\ell\nu\) to its original factor's cycle partition. Identity entries
in \(R\) are retained. Original identity factors are excluded when
required by the passport.

**Proof.** In successive copies of \(q_r,\ldots,q_1\), each intermediate
block is a fixed permutation of the blocks \(P_d^k(0)\). Thus each
edge \((i,a)\) occurs exactly once. Every tree-edge map in the normalized
coordinates is the identity. The remaining chronological maps compose
in the reverse order, so their product is the first return of the full
lifted product after \(d\) turns, restricted to root block zero.
Because the original product is \(P_N\) and \(F_0=1\), this is exactly
\(P_N^d|_0=P_b\). The return-budget theorem proves the index sum.

For the inverse construction, assign the supplied returns to the deleted
edges and identity maps to the tree edges. Let the product of these
normal-form factors have transports \(C_a\) over \(P_d\).
Set \(E_0=1\) and \(E_{a+1}=C_aE_a\) for \(0\le a<d-1\).
The boundary product in (5) says \(C_{d-1}E_{d-1}=P_b\).
Relabel fibre \(a\) by \(E_a^{-1}\). The full product then has identity
transport from block \(a\) to \(a+1\), except at the last block where
it has \(P_b\). In labels \(a+dt\), this is \(P_N\).

The local cycle rule follows from the return-cycle calculation in H3.
The index sum of the reconstructed factors is
\(b(d-1)+(b-1)=bd-1\). Finally, the tree normalization of the
reconstructed tuple recovers the same returns. Conversely, if one starts
with fixed product \(P_N\), the product charts \(E_a\) are precisely the
inverse normalization charts \(F_a^{-1}\). Reconstruction therefore
recovers the entire literal tuple, not only its relabeling class.
\(\square\)

## Theorem L3. Exact transport and the quotient-stabilizer action

In coordinates (3), multiplication and inverse are
\[
 \begin{split}
 (q,A)(s,B)&=(qs,(A_{s(a)}B_a)_a),\\
 (q,A)^{-1}&=(q^{-1},(A_{q^{-1}(a)}^{-1})_a).                \tag{6}
 \end{split}
\]
Lift a return tuple by L2, apply a signed original Hurwitz generator
using (6), and extract its new return tuple using L1--L2. This gives an
explicit invertible transport
\[
             (q,R)\longmapsto(Q_i^{\pm1}q,\mathcal T_i^{\pm1}R).
                                                                  \tag{7}
\]
It preserves exactly the original passport, full product and minimality.
For a quotient braid orbit \(\Omega\) and a chosen \(q_0\in\Omega\),
original components above \(\Omega\) are in bijection with the orbits of
\(\operatorname{Stab}_{B_r}(q_0)\) on the admissible return tuples over
\(q_0\), under (7).

A finite set of generators of this action is explicit. Choose a path
\(w_q\) from \(q_0\) to every quotient vertex. For every signed edge
\(q\xrightarrow{\beta}q'\), use the loop
\[
                         w_q\,\beta\,w_{q'}^{-1}.            \tag{8}
\]
Its return action is obtained by (7). This is a component reduction to
permutations of \(b\) letters and a quotient of degree \(d\); it is not
a claim that all inner braids are available.

**Proof.** Formula (6) is direct substitution in (3). L2 is a literal
bijection, so (7) is precisely the original braid action in those
coordinates and has the stated inverse. If two lifts over \(q_0\)
are connected, the projecting braid fixes \(q_0\). Conversely every
such stabilizer word connects its endpoints. Every lift above a different
quotient vertex can first be transported back along \(w_q^{-1}\).
Factoring a closed word edge by edge and inserting these paths writes
it as a product of the loops (8). Thus they generate the complete
stabilizer action, including all its coupling constraints. \(\square\)

A global original full twist conjugates by \(P_N^{-1}\). Its \(d\)-th
power fixes the literal quotient and conjugates every return by
\(P_b^{-1}\). Indeed one left-to-right sweep
\(Q_1\cdots Q_{r-1}\) moves the first factor last as
\(P_N^{-1}g_1P_N\); repeating the sweep \(r\) times conjugates all
factors. The return assertion then follows from the root charts and
\(P_N^d|_0=P_b\). These root phases are therefore available braid
moves, rather than free identifications silently added to the problem.

## Theorem L4. Decorated quotient components are sufficient for \(b=2,3\)

Decorate every quotient cycle slot by its return cycle partition.
Under a positive move on adjacent factors \(q_i=A,q_{i+1}=B\), cycles
of \(B\) move unchanged to position \(i\), while a cycle \(O\) of \(A\)
moves to position \(i+1\) as \(B^{-1}O\). The return partitions follow
those cycles. The inverse rule is obtained by reversing this move.

For block size two or three, the decorated quotient braid component is
a complete component invariant of the original lifts with those blocks.
There are no further gluing components above one decorated component.

**Proof.** The decoration rule is conjugation of a return map between
the corresponding blocks, hence preserves its cycle type.
For \(b=2\), the return budget is one: there is one transposition
return and every other return is the identity. Its value is forced by
the fixed product \(P_2\), so each decorated quotient has one literal
lift.

For \(b=3\), either there is one 3-cycle return, forced to be \(P_3\),
or there are two transposition returns in distinct slots. In the latter
case their ordered product is \(P_3\), giving exactly three possibilities.
Simultaneous conjugation by \(P_3\) is transitive on these three pairs.
The \(d\)-th full twist supplies that conjugation. Transport a connecting
decorated quotient braid first, then correct this remaining root phase.
This proves sufficiency; necessity follows from the decoration rule.
\(\square\)

More precisely, the literal lift multiplicity over a decorated quotient
is \(f=1\) for \(b=2\), \(f=1\) for the cubic 3-cycle case, and
\(f=3\) for the cubic two-transposition case. If its cyclic quotient
stabilizer has order \(s\), then the lifted tuple stabilizer has order
\(bs/f\). Indeed every quotient rotation fixing the decorations can
be corrected by a root phase, whose kernel has order \(b/f\).
Consequently corresponding absolute braid orbits have the same size,
although their literal sizes and deck groups can differ.

The threshold is sharp: L8 gives a failure at \(b=4\).

## Theorem L5. Complete lifting over a free cyclic branch

Fix \(b\ge2\) and a genus-zero polynomial quotient passport of degree
\(d\ge2\), with \(s\) nonidentity factors. In degree \(N=bd\), prescribe
one additional single \(b\)-cycle \(\tau\); every other factor's cycle
partition is \(b\) copies of its quotient partition. Restrict to the
union of components preserving blocks of size \(b\).

There is exactly one such original component above each quotient
component. If a quotient component has \(L\) literal tuples and group
\(K\), its lifted component has
\[
 \begin{array}{c}
 \text{group }C_b\wr K,\qquad
 \text{literal size }d(s+1)L,\qquad
 \text{deck-group order }b,\\[2mm]
 \text{absolute orbit size }(s+1)L .
 \end{array}                                                   \tag{9}
\]
The stabilizer action on the \(d\) choices over a fixed quotient tuple
with the identity position for \(\tau\) is exactly \(K\), in its
original degree-\(d\) action.

**Proof.** A factor supported on only \(b\) sheets cannot move a block
of size \(b\) to a different block, since this would move at least
\(2b\) sheets. Thus \(\tau\) fixes the quotient and lies in one block.
Its return index \(b-1\) uses the entire fibre budget. Every other
return is therefore the identity. Their cycle partitions must be
the prescribed repetitions, and their quotient types are uniquely
recovered by dividing the multiplicities by \(b\).

For a fixed literal quotient with its identity branch position, the
only choice is the marked block supporting the sole return \(P_b\).
L2 gives exactly \(d\) literal lifts. Normalize the identity position
to be last. Gather it next to any quotient factor by positive moves;
the other quotient permutations and the marked block stay literal.
Two moves on the adjacent pair return the quotient tuple and move the
marked block by that factor's inverse. Undo gathering. These loops
generate the action of \(K\), which is transitive. Conversely every
transport of this mark is a word in the quotient permutations, so
the image is exactly \(K\), not just a group containing it.
Thus each quotient component has one lifted component.

The quotient orbit with its identity inserted has \((s+1)L\) vertices.
Each has \(d\) lifts, proving the literal count. In normal tree
coordinates all unmarked factors act diagonally as \(K\), and the
remaining factor is a local \(b\)-cycle. Its \(K\)-conjugates generate
the \(d\) independent cyclic base factors, giving \(C_b\wr K\).
The relabeling that straightens the full product consists of cyclic
fibre powers, so this also gives the standard embedded wreath group
in the fixed labels \(a+dt\).

A commuting power of \(P_N\) must preserve the marked block, so it is
a power of \(P_N^d\). All \(b\) such powers commute with the independent
cyclic base factors and diagonal \(K\). This proves the deck order and
the absolute count in (9). \(\square\)

Other components of the same original passport can lack this block
system. They are not identified with this classified subdomain.
In particular H7 has an additional primitive \(A_{15}\) component.

## Lemma L6. Monodromy containing a transposition

Any finite transitive permutation group \(G\) containing a transposition
is a full imprimitive wreath product \(S_a\wr K\), including the case
\(a=N,K=1\) giving \(S_N\). Its blocks are the connected components
of the graph of supports of the \(G\)-conjugates of that transposition.
Their common size is \(a\ge2\), and the subgroup generated by those
conjugates is \(S_a^e\), where \(N=ae\).

This construction is independent of the chosen transposition in \(G\).
The natural centralizer is trivial when \(a\ge3\), and is the diagonal
\(C_2\) when \(a=2\).

**Proof.** Conjugation makes the support graph \(G\)-invariant.
Transpositions along a connected graph generate the full symmetric
group on its vertices. Thus its edge generators give \(S_a^e\);
transitivity of \(G\) makes all component sizes equal. The block action
is some transitive \(K\), and
\[
             S_a^e\le G\le S_a\wr K.
\]
Surjectivity onto \(K\) and containment of the full base imply equality.
Every transposition in the wreath product acts in one block, and
all such transpositions are conjugate using the base and \(K\).
This proves independence. A centralizer must preserve each block
and centralize its \(S_a\); its blockwise values must then be equal
because \(K\) is transitive. The stated centralizers follow.
\(\square\)

The graph, group and centralizer statements need no genus hypothesis.
For a polynomial tuple they combine with H3 and L2, but they do not
assert that its wreath-group label determines a component.

## Corollary L7. An infinite higher-index component catalogue

Take a degree-\(d=2(k+m)+1\) index-two quotient passport, \(k+m\ge1\), with \(k\)
3-cycles and \(m\) DTs. Replace every quotient factor by \(b\) disjoint
copies, add a single \(b\)-cycle, and retain the components preserving
the \(d\) blocks of size \(b\). Theorem L5 and the
[index-two classification](INDEX_TWO_HURWITZ_MIXTURES.md) give the complete
component list in this subdomain for every \(b\ge2\).

For \(k\ge1\) there is one component, with group \(C_b\wr A_d\)
(including \(A_3\)). For pure DT quotients: degree three has no tuples,
degree five gives one \(C_b\wr D_5\) component, degree seven gives two
\(C_b\wr A_7\) and two \(C_b\wr\operatorname{PSL}_3(2)\) components,
and every odd degree at least nine gives one \(C_b\wr A_d\) component.
The counts are (9), with the actual quotient orbit size for each case.
Writing \(s=k+m\), a connected nonempty quotient has
\[
                 L=\binom{s}{k}d^{s-1}(s-1)^m,              \tag{10}
\]
where the last factor is one when \(m=0\). This is the
[exact class-character count](HURWITZ_EXACT_CENSUS_AND_REMOVABLE_FACTORS.md#the-class-product-as-a-finite-difference)
for one ordered pattern, summed over its \(\binom{s}{k}\) placements.
For the four pure degree-seven components, use \(L=147,147,49,49\)
instead of their total.

In particular the passport on \(7b\) sheets
\[
         (2^{2b}1^{3b})^3,\qquad (b\,1^{6b})                 \tag{11}
\]
has four components in the specified block-size-\(b\) subdomain:
two of literal size \(4116\), absolute size \(588\), and group
\(C_b\wr A_7\); and two of literal size \(1372\), absolute size \(196\),
and group \(C_b\wr\operatorname{PSL}_3(2)\).
Each deck group has order \(b\). These formulas hold for all \(b\ge2\);
the literal sizes are independent of \(b\). The two full-\(A_7\)
quotient components yield the same standard embedded wreath group
but remain distinct components. \(\square\)

## Theorem L8. The first nontrivial return-coupling obstruction

On eight sheets, the two-factor passport
\[
                      (4,2,2),\qquad (2,2,1^4)              \tag{12}
\]
has three full components, of literal sizes \(16,16,8\).
The two size-sixteen components generate the identical subgroup
\[
 H=\{(u,v;\sigma)\in S_4\wr S_2:
                     \operatorname{sgn}(u)\operatorname{sgn}(v)
                           =\operatorname{sgn}(\sigma)\},    \tag{13}
\]
of order 576, preserving the two congruence blocks of size four.
Their deck groups are trivial. The third component has group order
64 and deck-group order two.

For both size-sixteen components the literal quotient is in the same
orbit of \(((0\,1),1)\). Every quotient cycle has transposition return
type. The inner degree-four tuple is in the single full braid component
of minimal factorizations of \(P_4\) into three transpositions.
Nevertheless the coupled action has two components.

**Proof.** Directly enumerate the 210 DTs on eight sheets and complete
the other factor using the product equation. There are exactly twenty
ordered pairs of the indicated types. Their complete full closures
have the three stated sizes and exhaust both type placements.

For \(q=((0\,1),1)\), there are three cycle slots. The all-transposition
return sector consists of the sixteen literal triples with product
\(P_4\). Direct completion enumerates all sixteen, and the full inner
\(B_3\) acts transitively on them. The outer stabilizer is generated by
\(Q_1^2\). Its return action, computed by (6)--(7), is a permutation
with two cycles of length eight. Thus this sector gives two original
full orbits of size sixteen. Original full tuples in these orbits
satisfy the parity condition (13), and exact generated-group closure
has order 576 in each case. Since the parity kernel has that order,
both groups are literally \(H\). Their block and deck data are also
checked exactly.

One may see their distinction without a return computation: in one
component the unique DT's two transpositions have cyclic distances
\(\{2,2\}\) along \(P_8\), and in the other they have
\(\{2,4\}\). In the indicated ordering, two braid moves only rotate
the pair, and an odd move puts the same literal DT first. These
distance multisets cannot change.

The remaining eight-tuple component has a different return decoration:
one transposition, one DT and an identity. Its exact group closure
and deck check give the stated order 64 and order two. The receipt
retains every representative and orbit digest. \(\square\)

This example shows why one quotient component and one inner full-braid
component, even together with the return partitions and actual monodromy
group, do not generally determine a glued component. The restricted
stabilizer action in L3 is indispensable.

## Proposition L9. A canonical divisor prevents decomposition overcounting

For an imprimitive tuple in (1), let
\[
 d_*=\max\{d:1<d<N,\ \mathcal B_d\text{ is preserved}\}.
                                                                  \tag{14}
\]
This chooses one proper quotient, of degree \(d_*\), for every original
braid component. Both \(d_*\) and the entire block-divisor lattice are
braid invariants. Classify the return action of L3 at each candidate
\(d\), keeping only reconstructed tuples with \(d=d_*\). The resulting
component lists are disjoint and exhaust the imprimitive domain.
Thus multiple polynomial decomposition chains require no extra
identification between these canonically selected lists.

**Proof.** H2 identifies every invariant partition and shows that their
divisor set is unchanged by a braid. The finite nonempty set of proper
divisors in (14) therefore has a unique largest numerical member.
The predicate \(d=d_*\) is invariant under the transport in L3.
Every imprimitive component occurs in that selected list and in no
other selected list. Its quotient and fibre degrees are both smaller
than \(N\).

In L8, both order-576 groups have only the proper block-count divisor
two, so their selected fibre size really is four. Choosing a different
decomposition cannot remove that coupling obstruction. The order-64
component has additional blocks and is selected at fibre size two.

For the full classification, identity quotient entries can be deleted
before classifying the quotient: they cross any factor unchanged, so
their insertion does not create new full braid components. Consequently
complete classifications of primitive \(S_N/A_N\) passports with at
least three finite values, the known rigid and exceptional cases, and
all the remaining restricted return actions would give the full
polynomial component classification by induction on \(N\). This is a
conditional reduction; those remaining classifications are not
established here. \(\square\)

## Evidence and remaining scope

The [checker](verify_polynomial_lift_transport.py),
[receipt](polynomial_lift_transport_verification.json) and
[source ledger](polynomial_lift_transport_source_ledger.json) retain the
literal inverse constructions, signed transport controls, direct
small-degree comparisons, complete small-block return domains,
quotient-loop permutations, the degree-eight obstruction and lifts
of all four pure degree-seven components. The unbounded results use
the written proofs, not extrapolation from these domains.

The full component-classification goal remains open. L3 replaces the
unspecified gluing action by an exact action on smaller-degree data.
L4 resolves its extra-fibre obstruction for \(b=2,3\), and L5--L7
classify a further all-degree cyclic-inner subdomain. Arbitrary
primitive \(S_N/A_N\) passports with at least three finite values,
and the remaining restricted return actions for larger fibres, still
require a complete classification. Independent review, independent
replay and formal verification are not asserted.
