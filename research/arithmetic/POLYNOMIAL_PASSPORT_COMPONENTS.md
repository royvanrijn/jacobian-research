# Components of polynomial passports with disconnected factors

The problem is to classify components for every finite passport, not to
assert connectivity after discarding proper monodromy. Fix an \(N\)-cycle
\(P=(0\,1\,\ldots\,N-1)\) and nonidentity permutations
\(g_1,\ldots,g_r\in S_N\) with
\[
 g_1\cdots g_r=P,\qquad \sum_i\operatorname{ind}(g_i)=N-1,
 \qquad \operatorname{ind}(g)=N-\#\operatorname{cycles}(g).       \tag{1}
\]
Fixed points count among the cycles. This is the genus-zero polynomial
domain: the cycle at infinity contributes another \(N-1\) to
Riemann--Hurwitz. Full braids use
\(Q_i(A,B)=(B,B^{-1}AB)\), with all placements of the prescribed finite
cycle types allowed. The absolute quotient is simultaneous conjugation
by \(\langle P\rangle\). Ordered or individually marked branch values
require a smaller braid group and are distinguished below.

The arbitrary-index-two subdomain is classified in the
[mixture theorem](INDEX_TWO_HURWITZ_MIXTURES.md). The full problem (1)
remains open in this repository. The results here supply a complete
two-factor classification, general block invariants and exact component
classifications of selected higher-index passports. They do not constitute
the complete classification for arbitrary \(r,N\).

## Primary literature and the scope of the obstructions

Jones--Zvonkin, [Theorem 1.1 and its discussion, Section 1.3 and Theorem 3.1](sources/jones-zvonkin-braid-cacti-2002.pdf),
identify polynomial topological classes with unrooted-cactus braid orbits,
give exact cactus counts and the central full-twist action, and explicitly
warn that a monodromy group need not determine an orbit. Their treatment of
two finite critical values already explains the rigid tree case below.
Their historical statement about the general classification is not used
as a claim about the current literature. No priority is asserted here.

Bétréma--Zvonkin,
[Sections 2--3, Main Theorem and the degree-partition lemma](sources/betrema-zvonkin-plane-trees-shabat-1996.pdf),
give the Shabat-polynomial correspondence with bicolored plane trees,
including the cyclic orders and distinction between the two colorings.
The algebraic two-strand calculation below keeps the precise quotient
and marking conventions explicit.

Müller's [Section 2.2 theorem, with the primitive hypothesis supplied by
Section 2.3](sources/mueller-primitive-polynomial-monodromy-author.pdf),
classifies primitive polynomial monodromy groups and possible finite
branch orders. The primitive/indecomposable qualifier, missing from the
isolated printed theorem sentence, is necessary and explicit in the
surrounding argument. Jones--Zvonkin Theorem 1.3 states its consequence
for three or more finite critical values. These classify possible
**groups**, rather than all components with those groups.

The earlier dihedral family has one braid orbit; it obstructs automatic
alternating monodromy, not connectivity of that passport. The five-sheet
three-3-cycle spin split has genus one and is outside (1). Genuine
genus-zero splitting is established by the examples below, including
splitting with the same primitive group and the same standard spin.

## Lemma H1. The cycle-incidence tree

Form a bipartite graph with \(N\) sheet vertices and one further vertex
for every nontrivial disjoint cycle of every \(g_i\). Join a cycle vertex
to each sheet in that cycle, and retain its cyclic order and factor label.
This incidence graph is a tree. Two cycle vertices with the same factor
label have disjoint sheet neighbors.

**Proof.** Its connected components are exactly the orbits of
\(\langle g_1,\ldots,g_r\rangle\), so it is connected. If there are
\(c\) cycle vertices, its edge count is
\(\sum (\ell-1)+c=N-1+c\), one less than its \(N+c\) vertices.
It is therefore a tree. The last assertion is disjoint-cycle
decomposition. Conversely, a tree of this form, with a cyclic order at
each cycle vertex and disjoint cycles within each factor, expands to
\(N-1\) transpositions on a tree. Their product in any order is a full
cycle: remove a leaf transposition, move it to the end while conjugating
its old endpoint, and insert its leaf into the lower cycle. Thus these
data represent precisely the minimal transitive cycle-support structures.
Conjugating one whole factor acts on all its cycles simultaneously;
independent braiding of cycle vertices changes the problem. \(\square\)

## Lemma H2. The invariant block-divisor lattice

Every block system for the generated group in (1) is, for a divisor
\(d\mid N\), the partition
\[
 \mathcal B_d=\{\{a,a+d,\ldots,a+(b-1)d\}:0\le a<d\},
 \qquad b=N/d.                                               \tag{2}
\]
It is invariant exactly when
\[
       g_i(x+d)\equiv g_i(x)\pmod d\quad\text{for every }i,x. \tag{3}
\]
The admissible divisors include \(1,N\), are closed under gcd and lcm,
and are unchanged by every Hurwitz move. Thus this divisor lattice is
an exact test for imprimitivity, not merely a test for selected blocks.

**Proof.** Every block system for the generated group is invariant under
its regular cyclic subgroup \(\langle P\rangle\). In a regular cyclic
action the block containing zero is a subgroup: the translations carrying
that block to itself are exactly its members. Its cosets give (2).
Condition (3) is precisely preservation of these congruence classes.
The intersection and join of congruence partitions are those for lcm and
gcd. Hurwitz moves preserve the generated subgroup as a literal subgroup
of \(S_N\), so they preserve all its block systems. \(\square\)

## Theorem H3. The quotient and fibre index budget

Fix an admissible divisor \(d\), and let \(q_i\in S_d\) be the induced
permutation on blocks. For each cycle \(O\) of \(q_i\), of length
\(\ell_O\), choose one of its blocks and let \(h_{i,O}\in S_b\) be
the return permutation \(g_i^{\ell_O}\) on that block. Its cycle type is
independent of the chosen block and of its labeling. Then
\[
 \begin{split}
 \operatorname{ind}(g_i)
     &=b\operatorname{ind}(q_i)+\sum_O\operatorname{ind}(h_{i,O}),\\
 \sum_i\operatorname{ind}(q_i)&=d-1,\qquad
 \sum_{i,O}\operatorname{ind}(h_{i,O})=b-1.                    \tag{4}
 \end{split}
\]
The quotient tuple has product \(P_d\) and is itself a genus-zero
polynomial tuple after identity entries are omitted. Projection commutes
with Hurwitz moves. The quotient braid component, and the multiset of
decorated branch types
\[
        \bigl\{(\ell_O,\operatorname{type}(h_{i,O})):O\bigr\}, \tag{5}
\]
are necessary component invariants. The statements do not assert that
these invariants are sufficient or that the return permutations can be
independently braided as a single inner tuple.

**Proof.** A return cycle of length \(s\) on one block corresponds to a
cycle of length \(\ell_Os\) on the union of its \(\ell_O\) blocks.
The contribution to the index is therefore
\(b(\ell_O-1)+\operatorname{ind}(h_{i,O})\). Sum over \(O\).
Projection is a homomorphism and sends \(P\) to \(P_d\). Factoring a
full \(d\)-cycle into the projected factors requires total index at least
\(d-1\): expand each into a minimal transposition word, whose edge graph
must be connected. The first identity and (1) give
\[
 d-1\le\sum_i\operatorname{ind}(q_i)
       \le\left\lfloor\frac{bd-1}{b}\right\rfloor=d-1.
\]
The two sums in (4) follow. Projection commutes with conjugation, hence
with both signed Hurwitz moves. Conjugation by a block-preserving
permutation only permutes the block cycles and conjugates their return
maps; this proves (5). \(\square\)

## Theorem H4. Complete classification with two finite branch values

Let \(r=2\), with arbitrary cycle types satisfying (1). Write
\(Q=Q_1\). Every literal full braid orbit is exactly
\[
 \bigl\{P^{-a}Q^\epsilon(g_1,g_2)P^a:
                   0\le a<N,\ \epsilon\in\{0,1\}\bigr\}.    \tag{6}
\]
The lexicographic minimum of this explicit set is a complete component
label. No braid search is needed. The pure braid orbit is the conjugacy
orbit under \(\langle P\rangle\).

Equivalently, for a fixed ordered passport the absolute pure-braid
components are the bicolored plane trees with the two prescribed degree
partitions. If the types are distinct, every full component meets the
chosen type ordering in exactly one such tree. If the types coincide,
full components additionally identify the two colorings of a tree.
Cyclic orders are retained; plane reflection is not an allowed
identification.

**Proof.** Substitution using \(g_1g_2=P\) gives
\[
           Q^2(g_1,g_2)=(P^{-1}g_1P,P^{-1}g_2P).             \tag{7}
\]
Since \(B_2=\langle Q\rangle\) and \(PB_2=\langle Q^2\rangle\),
(6) and the pure statement follow. For the tree description, make one
black vertex per cycle of \(g_1\), one white vertex per cycle of \(g_2\),
and one edge labeled \(x\) between the cycles containing \(x\). The
number of vertices is \(N+1\), by (1), and transitivity makes the graph
connected, hence a tree. The two permutations give its cyclic orders.
Conversely these rotations reconstruct the pair and its full-cycle face
permutation. Relabeling that cycle to \(P\) leaves exactly its cyclic
centralizer as residual ambiguity. A half twist interchanges the two
vertex colors, and (7) accounts for all remaining action. \(\square\)

If \(s\) is the order of the pair's cyclic centralizer, its pure literal
orbit has size \(N/s\). Its full literal orbit has size \(N/s\) when
\(Q(g_1,g_2)\) is already a cyclic conjugate, and \(2N/s\) otherwise.
This includes symmetric trees; division by \(N\) without checking \(s\)
would give incorrect absolute counts.

For one finite branch value the only fixed-product tuple is \((P)\).
With no finite branch values, (1) forces \(N=1\) and the empty tuple.

## Proposition H5. Standard spin is forced for every odd-order passport

Suppose all finite factors in (1) have odd order. Then \(N\) is odd,
and the standard spin invariant including infinity is the single value
\[
 s=(-1)^{\sum_i w(g_i)+w(P^{-1})},\qquad
 w(g)=\sum_{\ell\text{ a cycle length of }g}\frac{\ell^2-1}{8}.
                                                                  \tag{8}
\]
In particular this standard spin cannot distinguish components of any
one such genus-zero passport, even when its factors are disconnected.

**Proof.** All finite indices are even, so \(N-1\) is even. The tuple
\((g_1,\ldots,g_r,P^{-1})\) is transitive, product one, genus zero,
and consists of odd-order permutations. These are exactly the hypotheses
of Fried's [Invariance Corollary 2.3 and its proof, Section 2.3](sources/fried-lifting-invariants-math-0611591v4.pdf).
Its formula gives (8), whose exponent depends only on the cycle types.
Other central covers and separately prescribed lifted branch elements
are distinct questions. \(\square\)

## Theorem H6. Six primitive components with the same standard spin

For \(N=11\), the two finite types
\[
                    (3^2 1^5),\qquad (3^3 1^2)              \tag{9}
\]
give exactly six absolute full braid components. Each has 22 literal
tuples when both type placements are included, and 11 in the indicated
ordering. Every component has group \(A_{11}\), no nontrivial
block divisors and standard spin \(+1\). Thus the monodromy group,
imprimitivity data and standard spin together do not classify components.

**Proof.** Directly enumerate all 18,480 permutations of type
\(3^2 1^5\), and complete the other factor as \(g_1^{-1}P\). Exactly
66 completions have type \(3^3 1^2\). Formula (6) partitions these
and their reverse placements into six disjoint 22-tuple orbits. Their
complete literal representatives and orbit digests are retained. For
each representative an explicit group word in the two factors is a
single 3-cycle. The degree is prime, so the transitive group is primitive;
Jordan's theorem implies \(A_{11}\), since both generators are even.
There are no nontrivial block divisors. The two finite indices are
\(4,6\), so this really is genus zero. Formula (8) gives
\(s=(-1)^{2+3+(11^2-1)/8}=+1\). Exact Clifford controls check the
same lift convention. Each finite cycle class is also a single
\(A_{11}\)-class: interchanging two of its fixed points is an odd
centralizer element. The long cycle is fixed throughout, so these
components do not arise from its two alternating conjugacy classes.
No assertion about an independent review is made.
\(\square\)

In degree nine, the equal-type passport \((3^2 1^3),(3^2 1^3)\)
similarly has three full components of literal sizes \(18,9,9\), and
four color-ordered absolute pure components. All three generated groups
have order 1512 and all standard spins are \(+1\). This illustrates the
color-interchange clause of Theorem H4; it is not needed for H6.

## Theorem H7. A genus-zero primitive/imprimitive split on fifteen sheets

The unordered three-factor passport
\[
                    (3^3 1^6),\ (3^3 1^6),\ (3 1^{12})      \tag{10}
\]
has exactly two full braid components:

| Group | Literal orbit, all three placements | Absolute orbit size | Nontrivial block divisors |
|---|---:|---:|---|
| \(A_{15}\) | 58,725 | 3,915 | none |
| \(C_3\wr A_5=C_3^5\rtimes A_5\), order 14,580 | 75 | 15 | \(d=5\) |

Both have standard spin \(-1\). On any fixed ordered passport their
pure braid components have 19,575 and 25 literal tuples, respectively.
There is no further splitting from individually marking the two
equal-type factors.

**Proof.** The independent character calculation below gives 19,600
literal tuples in one type ordering, and 58,800 over its three placements.
The checker closes the two retained seeds under every signed generator;
their disjoint orbit sizes add to 58,800. This proves completeness.
For the first seed, Theorem H2 excludes both possible nontrivial divisors
\(3,5\). Its literal 3-cycle and Jordan's theorem therefore give
\(A_{15}\). For the second seed, first work on five blocks of three
points. Lift the two permutations \((0\,1\,2),(2\,3\,4)\) diagonally
to the blocks, and add one 3-cycle inside block zero. The product is a
15-cycle, with the types in (10). The diagonal generators give \(A_5\)
on blocks; conjugating the local 3-cycle gives the five independent base
generators of \(C_3^5\). This proves the claimed wreath product, also
checked by exact finite group closure. The block and quotient controls
give only \(d=5\); the quotient types are \(3,3,1\) and its index sum
is four, while the fibre return index sum is two, as in (4).

The first group's centralizer is trivial. For the wreath product,
commuting with the five independent block 3-cycles forces a centralizer
element into \(C_3^5\); commuting with the transitive diagonal \(A_5\)
then makes its five powers equal. The centralizer is the diagonal
\(C_3\), of order three. The checker also tests every centralizer power
on every tuple of both complete finite orbits. Thus the cyclic
conjugacy classes have sizes fifteen and five, giving the absolute
sizes in the table. This deck symmetry must not be lost by simply
dividing both literal counts by fifteen. Complete pure braid closure
has one third of each full orbit, exhausting its fixed-pattern slice.
Finally (8) gives \((-1)^{3+3+1+(15^2-1)/8}=-1\).
\(\square\)

## Lemma H8. Refinement can bridge rigid components in a different passport

Split any one of the three disjoint 3-cycles off the second factor in
each ordered representative of (9). The resulting 18 triples all have
the fixed product \(P_{11}\) and ordered types
\((3^2 1^5),(3^2 1^5),(3 1^8)\). They lie in one full braid orbit of
size 3,267, which is the entire three-placement domain of that passport.
The receipt supplies an ordinary braid word from each triple to a common
seed. Coalescing its last two factors recovers its original rigid pair.

**Proof.** The disjoint cycles commute, so splitting preserves the product
and minimal total index. Complete signed braid closure, compared with
the exact independent total, proves the assertion and retains the bridge
words. This is a finite statement about the refined passport. The
intermediate braid does not preserve the coalescing constraint, so it
does not connect the six original two-factor components. \(\square\)

## Proposition H9. The primitive group reduction

For \(r\ge3\), a primitive group in (1), other than \(S_N,A_N\),
can occur only for \(r=3\) in these cases, up to permutation of the
finite factors:

| Degree | Finite passport | Proper primitive group |
|---:|---|---|
| 7 | \((2^2 1^3)^3\) | \(\operatorname{PSL}_3(2)\) |
| 13 | \((2^4 1^5)^3\) | \(\operatorname{PSL}_3(3)\) |
| 15 | \(2^6 1^3,\ 2^4 1^7,\ 2^4 1^7\) | \(\operatorname{PSL}_4(2)\) |

**Proof.** This is Jones--Zvonkin Theorem 1.3, with the finite-order
list in Müller's Section 2.2 theorem and its primitive hypothesis in
Section 2.3. The actual statements and the reduction context were
inspected. No component uniqueness is inferred from a group on this
list. Jones--Zvonkin Sections 5.2--5.4 give the separate component
classifications of these exceptional passports. \(\square\)

## Theorem H10. Imported exceptional-passport component classification

The following are complete absolute full-braid orbit lists in
Jones--Zvonkin Sections 5.2--5.4, with all distinct type placements
included:

| Degree and finite passport | Monodromy group | Absolute orbit sizes |
|---|---|---|
| \(7:\ (2^2 1^3)^3\) | \(A_7\) | \(21,21\) |
| same | \(\operatorname{PSL}_3(2)\) | \(7,7\) |
| \(13:\ (2^4 1^5)^3\) | \(A_{13}\) | \(35\,620\) |
| same | \(\operatorname{PSL}_3(3)\) | \(13,13,13,13\) |
| \(15:\ 2^6 1^3,\ 2^4 1^7,\ 2^4 1^7\) | \(A_{15}\) | \(377\,835\) |
| same | \((S_3\wr S_5)\cap A_{15}\) | \(135\) |
| same | \(\operatorname{PSL}_4(2)\) | \(15,15\) |

The totals are \(56\), \(35\,672\) and \(378\,000\), respectively.
All tuple stabilizers under cyclic sheet relabeling are trivial,
so literal sizes are the displayed sizes
times the degree. In degree fifteen, this full domain includes the
three placements of the \(2^6 1^3\) factor.

**Source and assurance.** These are imported component classifications,
not a new internal replay of the degree-thirteen or degree-fifteen
computations. The actual source statements, their cactus-count coverage
arguments and their scope were inspected. The degree-seven table also
agrees with the repository's retained independent code path for that
complete small domain. The groups and orbit sizes in this theorem are
different from the degree-fifteen odd-order passport in H7.

## Exact counting and retained evidence

For arbitrary finite cycle classes \(C_i\subset S_N\), the derivation of
the [class-product formula](HURWITZ_EXACT_CENSUS_AND_REMOVABLE_FACTORS.md#the-class-product-as-a-finite-difference)
gives, for a fixed long cycle,
\[
 \#\{g_1\cdots g_r=P\}
 =\frac{\prod_i|C_i|}{N!}
   \sum_{j=0}^{N-1}(-1)^j
       \frac{\prod_i\chi_j(C_i)}{\binom{N-1}{j}^{r-1}},       \tag{11}
\]
where \(chi_j\) is the \(j\)-th exterior power of the standard
representation. For a cycle partition (lambda),
\[
 \sum_j\chi_j(\lambda)x^j
            =\frac{\prod_{\ell\in\lambda}(1-(-x)^\ell)}{1+x}.
                                                                  \tag{12}
\]
The retained checker uses exact rational arithmetic in (11)--(12).
For comparison, the Goulden--Jackson expression stated in
Jones--Zvonkin Section 1.3 is
\(N^{r-1}\prod_i (c_i-1)!/\prod_\ell m_{i,\ell}!\).
The individual factors can be fractional: in (10) each higher factor
contributes (28/3). Rounding these factors would corrupt the total.
Completeness of H7--H8 uses (11), not that comparison formula.

The [checker](verify_polynomial_passport_components.py),
[receipt](polynomial_passport_components_verification.json) and
[source ledger](polynomial_passport_components_source_ledger.json)
retain exact representatives, group 3-cycle words, orbit digests,
two-strand normal forms, centralizer checks, block quotients and return
budgets, spin controls and the 18 refinement words. These are complete
finite classifications of the stated passports and controls for the
written unbounded reductions. They do not infer a classification of
other passports by extrapolation. Replay and proof audit are internal;
independent mathematical review and formal verification remain open.

## The remaining complete-classification obligation

The target still ranges over every \(N,r\) and every passport in (1).
Theorem H4 resolves the rigid \(r\le2\) portion, including arbitrary
plane-tree types. For \(r\ge3\), the primitive group list has only
\(S_N,A_N\) and the three exceptional primitive cases in
Jones--Zvonkin Theorem 1.3. H10 imports their complete exceptional
component tables; this is distinct from an independent replay. Those
tables do not resolve arbitrary \(S_N,A_N\) passports.

The unresolved structural tasks are:

1. Determine all components for primitive \(S_N,A_N\) passports with
   at least three finite branch values. A group or standard-spin label
   is not a sufficient theorem; H6 already disproves that proposed
   criterion when the rigid boundary is included.
2. Classify the joint braid action on lifts above each quotient component
   in H3, keeping the complete block-divisor lattice and decorated return
   data. Separate outer and inner component lists do not supply this
   gluing action, and polynomial decomposition chains need not be unique.
3. Prove extraction and restoration within these actual components, or
   replace extraction by a complete normal form for the cycle-incidence
   tree under whole-factor braids. H8 demonstrates a possible refinement
   mechanism but does not prove a general descent or a coalescing rule.
4. Determine any remaining reduced lifting invariants in their stated
   domains, without importing a positive-genus spin split into (1).

A terminating orbit search for one supplied passport is useful evidence.
It is not by itself the requested structural classification for all
passports. No completeness claim for this remaining scope is made here.
