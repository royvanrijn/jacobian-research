# Cycle folding for alternating Laurent Hurwitz tuples

Every transitive tuple in the even-degree Laurent passport can be braided to
a graph consisting of a doubled edge at the fixed sheet and an attached tree.
The doubled edge can be placed in two adjacent branch permutations, with no
other edge incident to the fixed sheet. Coalescing these two permutations
then gives an odd-degree alternating polynomial tuple.

This is a uniform normal-form reduction. It does not prove one complex
Hurwitz component in degree twelve or higher: connectivity of the remaining
colored-tree forms, and compatibility of their expansions, remain open here.

## Theorem 1. The unicyclic graph and fixed-sheet trajectory

Let n=2b>=8 and let sigma_1,...,sigma_b be double transpositions on
{0,...,n-1}, generating a transitive group, whose product P is an
(n-1)-cycle fixing 0. Use composition pq=p composed with q.
Expand each sigma_i into its two disjoint transpositions and make a
multigraph on the n symbols, with one edge occurrence for each
transposition.

This graph is connected and has exactly one cycle, counting a parallel
pair as a cycle of length two. That cycle contains 0, has length at
most b, and has at most one edge in each sigma_i. Its edges appear in
cyclic order when the fixed-sheet trajectory is followed through the
product from right to left.

**Proof.** The orbits of the generated group are exactly the connected
components of this graph: each edge is a swap in one generator, so its
endpoints lie in the same group orbit. Thus it is connected, with n
vertices and 2b=n edge occurrences. Its cycle rank is one.

Follow the image of 0 while applying the expanded product right to left.
It starts and ends at 0. Every edge occurrence can be traversed at most
once. The resulting nonempty closed trail cannot use a bridge, since
returning across the bridge would require a second occurrence of the
same edge. Consequently it traverses the unique cycle once and uses
no tree edge. It is nonempty because the graph has an edge incident to
0. Thus 0 lies on the cycle and the trajectory gives its cyclic edge
ordering.

If two edges of one double transposition lay on this cycle, they would
be adjacent in the expanded word and hence consecutive in its trajectory.
They would share a vertex, contradicting their disjointness. There is
therefore at most one cycle edge per block and at most b cycle edges. ∎

The graph uses edge occurrences. A doubled edge contributes two of the
n occurrences; its underlying simple graph is a tree.

## Theorem 2. A uniform doubled-edge normal form

Under the hypotheses of Theorem 1, there is a braid word of length at
most (b+1)(b-2) taking the tuple to one whose graph has:

* a unique doubled edge (0,a);
* its two occurrences in adjacent double-transposition blocks;
* no further edge incident to 0;
* a tree on the other n-1 symbols after deleting 0.

The braid preserves both the product P and the generated group.

**Proof.** The positive and negative Hurwitz moves are

\[
 (x,y)\longmapsto(y,yxy),\qquad
 (x,y)\longmapsto(xyx,x).
\]

They are inverse and preserve product and generated subgroup.

First observe a transport operation. Suppose x carries one cycle edge
and the next block y carries none. The negative move keeps x unchanged
and moves it one position to the right. All old cycle edges remain in
the graph: no edge removed from y was on that cycle. The new graph
is still connected with n edges and n vertices. Its unique cycle
therefore remains exactly the old one. In particular the new first
block xyx has no cycle edge. A cycle block can thus be transported
past any intervening tree-only blocks without changing the cycle.

Suppose now its length k is at least three. Choose two consecutive
cycle edges in the product's linear edge order. Their shared cycle
vertex is not 0, which is the start and end of the fixed-sheet
trajectory. Transport the earlier block until the two blocks are
adjacent. Write their cycle edges as (a,v),(v,c), with v!=0, and the
blocks as

\[
                         x=(a\,v)p,\qquad y=(v\,c)q.
\]

Here p,q are tree edges, disjoint from their own displayed cycle edges.

If q does not meet a, a positive move replaces the first cycle edge
by (a,c). The unchanged outside path from a to c closes a cycle of
length k-1. Connectedness and cycle rank one imply this is the new
unique cycle. If p does not meet c, the same shortening follows
from a negative move instead.

The remaining case is p=(c,d), q=(a,e). The tree property implies
that d,e are distinct and lie outside the cycle: a path of tree edges
joining two old cycle vertices would create a second cycle.
Two positive moves give

\[
 \bigl((a\,v)(c\,d),(v\,c)(a\,e)\bigr)
 \longmapsto
 \bigl((c\,e)(v\,d),(a\,c)(d\,e)\bigr).
\]

Again the chord (a,c) shortens the old cycle by one. All factors
remain double transpositions because these are actual Hurwitz moves.
Repeating gives a doubled edge incident to 0.

There are at most b-2 shortening steps. Each uses at most b-2
transport moves and two local moves, hence at most b moves.
Finally transport one doubled-edge block until the two occurrences
are adjacent, using at most b-2 additional moves. This proves the
stated length bound.

Apply Theorem 1's trajectory to this final tuple. An edge incident
to 0 outside the doubled edge could only occur strictly between
the two doubled-edge occurrences: before the first and after the
second, the trajectory is at 0 and would traverse that tree edge.
Adjacency leaves no intervening block, and the secondary transpositions
in the two blocks are disjoint from 0. Hence no such tree edge exists.
Deleting 0 and the doubled pair leaves a connected graph with n-1
vertices and n-2 edges, namely a tree. ∎

This normalization applies equally to proper transitive groups in the
exceptional degree eight. It preserves the given group rather than
changing it to A8.

## Corollary 3. Reduction to odd-degree alternating polynomial tuples

In the normal form, write the adjacent doubled-edge blocks as

\[
                      (0\,a)p,\qquad(0\,a)q.
\]

Coalesce them to pq and delete the fixed symbol 0. The resulting
(n-1)-point tuple has b-1 entries, product an (n-1)-cycle, and every
entry has permutation index two. All entries are double transpositions
except possibly the coalesced entry, which can be a 3-cycle.

For n>=10 this polynomial tuple has full group A_(n-1). In degree
eight its group is left to the corresponding finite polynomial analysis.

**Proof.** The secondary edges p,q are distinct, since equality would
give another parallel cycle. They are either disjoint or meet once.
Their product is accordingly a double transposition or a 3-cycle.
Every other block fixes 0, by Theorem 2. The remaining graph is the
tree just constructed, and the product is P restricted to the other
symbols. The tuple is transitive because that product is a full cycle.

The index sum is 2(b-1)=n-2. Together with the (n-1)-cycle at infinity,
it gives genus zero by Riemann--Hurwitz. Riemann existence therefore
realizes it as a polynomial cover of degree n-1.

Put N=n-1. Its action is primitive. Any nontrivial block size in odd
degree N is odd and at least three. A permutation moving at most four
points cannot permute two such blocks, which would move at least six
points. Thus every finite branch generator would fix each block
individually. Their product is a full N-cycle, which cannot fix a
nontrivial block system pointwise on its set of blocks.

For n>=10 we have N>=9 and at least one double transposition in the
tuple. The retained primitive double-transposition theorem, Jones,
Proposition 6.1(3), gives A_N. All entries are even, so equality
follows. If the coalesced entry is a 3-cycle, Jordan's prime-cycle
theorem gives the same conclusion already for N>=6. The exact
imported statements and primary inspections are in the
[priority ledger](alternating_hurwitz_priority_ledger.json) and
[earlier Hurwitz source ledger](alternating_hurwitz_source_ledger.json). ∎

The coalescence is a reduction of tuples. It does not identify the
braid orbits of the lower-degree polynomial tuples, or show that all
expansions of a coalesced tuple lie in one orbit.

## Corollary 4. A uniform link between the two coalescence types

For every even n=2b>=8, there is a full-A_n tuple in this passport
with an adjacent doubled-edge pair coalescing to a double transposition,
and the same two-letter braid word converts it to an adjacent pair
coalescing to a 3-cycle.

On {0,...,7} start with

\[
 ((0\,1)(2\,3),\ (0\,1)(4\,5),\
                  (3\,4)(6\,7),\ (1\,2)(5\,6)).
\]

For j=4,...,b-1 append

\[
                              (1\,2j)(6\,2j+1).                \tag{1}
\]

The required braid word is (-2,1), denoting Q_2^-1 followed by Q_1.
It changes the first three entries to

\[
              ((3\,5)(6\,7),\
                 (0\,1)(2\,5),\ (0\,1)(4\,5)).                 \tag{2}
\]

Thus the first two entries of the initial tuple coalesce to a double
transposition, while the second and third entries in (2) coalesce
to a 3-cycle. The common doubled edge remains (0,1).

**Proof.** After deleting 0 and one copy of the doubled edge, the initial
graph on {1,...,n-1} consists of the path 1-2-3-4-5-6-7 and the two
leaf edges supplied by each appended block (1). It is a tree.
The coalesced product is a long (n-1)-cycle.

Here is the elementary tree-product fact used in that assertion.
A product using every edge transposition of a tree exactly once is a
cycle on its vertices, regardless of their order. Remove a leaf l,
whose unique factor is (l v), and write the product as L(l v)R.
Both L and R fix l, and

\[
                       L(l\,v)R=(l\,L(v))(LR).
\]

Induction makes LR a cycle on the remaining vertices. The displayed
transposition inserts l into it. This proves the fact. It applies
to the expanded coalesced tuple, since each remaining tree edge
occurs once.

Consequently the original product is an (n-1)-cycle fixing 0, and
the generated group is transitive. The two Hurwitz moves give (2)
directly and leave every appended block unchanged.

For completeness, the group is A_n even in degree eight. Coalesce
the pair in (2) and delete 0. This is an odd-degree transitive
polynomial tuple with a 3-cycle and all other entries double
transpositions. The small-support block argument of Corollary 3
makes it primitive, and Jordan's 3-cycle theorem makes its group
A_(n-1), since n-1>=7. The original group contains this A_(n-1)
fixing 0 and is transitive on n points. Its order is therefore at
least n*(n-1)!/2=n!/2. Since all original generators are even, it
is exactly A_n. ∎

This links at least one expansion of each coalescence type uniformly.
It does not by itself connect different marked polynomial components
or different attachment choices.

## Corollary 5. A precise uniform connectivity criterion

Let n=2b>=10, N=n-1 and r=b-1. On N symbols fix a full cycle P,
write D for the double-transposition class and T for the 3-cycle class,
and mark the first entry in the two fixed-product spaces

\[
 \mathcal P_D=\{(\delta,g_2,\ldots,g_r)\in D^r:
                                      \delta g_2\cdots g_r=P\},
\]
\[
 \mathcal P_T=\{(\delta,g_2,\ldots,g_r)\in T\times D^{r-1}:
                                      \delta g_2\cdots g_r=P\}.
\]

Suppose the operations Q_2,...,Q_(r-1),Q_1^2 and their inverses
are transitive on each space. Suppose also that, at one base tuple
in each space, cabled lifts of polynomial stabilizer loops are
transitive on the choices a outside the support of the marked delta,
after allowing the internal two-factor splitting moves.
Then the original degree-n absolute full-A_n Hurwitz space has
one braid orbit.

**Proof.** Theorem 2 reduces every original tuple to an adjacent doubled
edge. Bring those two blocks first, coalesce them and delete the fixed
sheet. The result is one of the marked polynomial tuples above.
Its expansions are exactly

\[
              (\star\,a)p,\quad(\star\,a)q,\qquad pq=\delta,
                              \quad a\notin\operatorname{supp}\delta.
\]

An ordinary marked-preserving Q_j with j>=2 lifts to Q_(j+1).
The marked operation Q_1^2 lifts to the four-letter word (2,1,1,2)
on the expanded first three blocks; inverse words lift as inverses.
Thus polynomial orbit paths carry every expansion to an expansion
of the chosen base tuple.

For a fixed a the internal B2 action is transitive on the two
transposition splittings of a double transposition or the three
splittings of a 3-cycle. The assumed stabilizer action connects
the remaining attachment choices. Each coalescence type therefore
gives one original orbit. Corollary 4 joins these two orbits. ∎

The stabilizer hypothesis is about its actual cabled action. It cannot
be replaced simply by transitivity of the abstract centralizer of delta
on its complement. The theorem in
[degree twelve](DEGREE_TWELVE_HURWITZ_ORBIT.md) verifies precisely
these hypotheses for the two marked degree-eleven polynomial spaces.

## Corollary 6. Uniform attachment transitivity

For every r>=3 and N=2r+1, each of the two marked polynomial
passports D^r and T D^(r-1) has an explicit base tuple with
full group A_N whose cabled stabilizer acts transitively on all
unramified attachment choices outside the marked permutation.
Consequently the attachment hypothesis of Corollary 5 requires no
further degree-dependent proof.

**Proof.** A positive marked move Q_1^2 sends an adjacent pair
(delta,gamma) to

\[
 (\gamma^{-1}\delta\gamma,\,
   \gamma^{-1}\delta^{-1}\gamma\delta\gamma)
\]

and sends the attachment a to gamma(a). These formulas follow either
from the doubled-edge cable (2,1,1,2) or by two Hurwitz moves. A
negative marked move transports a by delta gamma delta^-1.
Moves involving only unmarked factors leave the attachment unchanged.

In the following table code 1 denotes Q_1^2, code 2 denotes Q_2,
and a negative code denotes the inverse. Words act from left to right.
Each row is a stabilizer loop of the indicated literal base tuple;
the final column is its full transport permutation on seven symbols.

| Base tuple, with its first entry marked | Word | Attachment transport |
|---|---|---|
| ((01)(34),(12)(45),(23)(56)) | (2,1,1,1,-2) | (34)(56) |
| Same DT base | (1,2,2,-1,-1,-2,-2,-1) | (01)(25) |
| ((012),(23)(45),(34)(56)) | (2,1,-2) | (34)(56) |
| Same 3-cycle base | (1,1,2,-1,-2,1,-2,1,2) | (012)(36)(45) |

These identities are finite direct substitutions in the displayed
Hurwitz and transport formulas. The
[exact checker](verify_uniform_hurwitz_attachment.py) also replays
every literal cable, for every attachment and every two-transposition
splitting of the marked factor. No enumeration or presumed orbit
connectivity is used in these four identities.

The DT base has outside set {2,5,6}. Its two restrictions are
(56) and (25), so they act transitively. The 3-cycle base has
outside set {3,4,5,6}; its restrictions (34)(56) and (36)(45)
also act transitively.

Both bases have full-cycle product. Expand the 3-cycle into (01)(12)
in the second base: the six transposition edges in either base are
exactly the path 0-1-2-3-4-5-6. The tree-product argument from
Corollary 4 applies.

To increase r by one, append two fresh symbols x,y and the DT

\[
                         \eta=(a_0\,x)(a_1\,y),
\]

where a_0,a_1 are distinct old outside symbols. Use the fixed anchors
(2,5) for the DT family and (3,4) for the 3-cycle family. The edge
graph acquires two leaves, so the enlarged product remains a full
cycle. Each old stabilizer word still acts on the first three factors
and retains its old transport, fixing every new symbol.

The appended eta commutes with the marked delta. Move eta unchanged
to the second position by positive moves among unmarked factors,
apply Q_1^2, and undo those moves. This is a polynomial stabilizer
loop whose attachment action is exactly eta. If eta initially has
zero-based factor index j>=3, the compressed loop is

\[
                  (j,j-1,\ldots,2,\ 1,\ -2,\ldots,-j).
\]

Its two transposition edges connect x,y to the old transitive
attachment orbit. Repeating the construction gives attachment
transitivity for every r>=3.

For completeness these bases have full alternating groups.
The 3-cycle base already contains a 3-cycle. In the DT base,
the marked delta multiplied by the two transport permutations in
its table gives the 3-cycle (2 6 5). Those transports belong to the
generated polynomial group, since every transport step is a word
in that group. Both base groups are transitive of prime degree seven,
hence primitive; Jordan's 3-cycle theorem gives A_7.
The appended families remain transitive and contain this 3-cycle.
At their odd degrees the small-support block argument of Corollary 3
makes them primitive, so the same theorem gives A_N. ∎

The construction is a uniform statement about one explicit base
tuple in each marked space. It does not prove that every marked
polynomial tuple can be connected to that base. Combining
Corollaries 4-6 shows that, for n>=10, **transitivity of the two
marked polynomial spaces alone suffices** for full connectivity of
the original degree-n passport. The attachment condition is no
longer a separate frontier.

The [attachment receipt](uniform_hurwitz_attachment_verification.json)
records the four exact seed identities, all attachment/splitting cable
checks through polynomial degree thirteen, and deterministic
leaf-extension controls through polynomial degree fifty-one.
Uniformity follows from the written leaf-extension proof, not
from the finite range of controls.

## The remaining connectivity problem

The original fixed-product space can now be studied through colored
trees with one marked adjacent doubled-edge expansion. To obtain
uniform braid transitivity by degree induction one still needs, for
example, a valid operation producing a block with two leaf edges whose
leaf labels occur nowhere else, followed by control of the reduced
orbit and all expansion choices. No such operation is proved here.

The polynomial-cactus viewpoint is established prior work. Jones and
Zvonkin's [primary article](https://www.labri.fr/perso/zvonkin/Research/jones-zv.pdf),
Theorem 1.1 and its discussion, identifies polynomial topological types
with braid orbits and explicitly warns that equal monodromy groups need
not mean equal orbits. The reduction above therefore cannot be completed
merely by identifying A_(n-1) as the lower group's name. The
[source ledger](alternating_hurwitz_cycle_folding_source_ledger.json)
pins the inspected version; no universal connectivity statement for this
exact passport was verified from it.

An inspected earlier uniqueness theorem does not fill this gap.
Khovanskii and Zdravkovska's
[Theorem 3, printed page 66](https://www.math.utoronto.ca/askold/1996-J-Knot-theor-5.pdf)
requires the defect, summed over nonsimple finite critical values,
to be at most the polynomial degree plus one. A double-transposition
value contributes 4 and a 3-cycle value contributes 3. Our degree
N=2b-1 passports therefore have defects 4b-4 and 4b-5, respectively;
both exceed N+1=2b for b>=3. Their unmarked uniqueness theorem consequently
does not apply, and does not supply a cabled attachment calculation.

The marked polynomial conditions are verified in
[degree twelve](DEGREE_TWELVE_HURWITZ_ORBIT.md).
Corollary 6 makes the attachment condition uniform.
The [degree-fourteen theorem](DEGREE_FOURTEEN_HURWITZ_ORBIT.md)
uses the uniform reduction and complete finite censuses to prove one
absolute full-alternating orbit. Degree sixteen and uniform connectivity
remain open.

## Exact verification

The [checker](verify_alternating_hurwitz_cycle_folding.py) constructs and
replays explicit braid words. Its
[receipt](alternating_hurwitz_cycle_folding_verification.json) records:

* exhaustive normalization of all 49,392 transitive degree-eight DT
  quadruples with fixed product of type (7,1);
* initial cycle counts 28,224 of length two, 18,144 of length three,
  and 3,024 of length four;
* exact product, tree and adjacent-block checks for every resulting word;
* 1,000 deterministic exact braid-translate controls in each of degrees
  ten, twelve and twenty-six, and their coalesced polynomial tuples.

The degree-eight domain includes the 47,334 A8 tuples and 2,058 tuples
in proper transitive AGL3(2) subgroups. Its completeness is the same
two-pair product enumeration used by the
[degree-eight orbit theorem](DEGREE_EIGHT_HURWITZ_ORBIT.md).
A tuple is transitive exactly when some generator moves the unique
fixed sheet of P. The checker uses this exact criterion.

The higher-degree tuples are controls generated from the retained
[real matching constructors](verify_real_paired_laurent_matching.py).
They are a finite sample within their known orbit, not a census.
The all-degree result follows from the graph and explicit local-move
proofs above.

Run from the repository root:

~~~sh
.venv/bin/python research/arithmetic/verify_alternating_hurwitz_cycle_folding.py
~~~

Only the Python standard library is required. The computation and proof
checking are internal. No formal verification, independent external
replay, external human review or priority claim is made.

