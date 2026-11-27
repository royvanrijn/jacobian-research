# Terminal mass and Hurwitz descent with changing candidate leaves

There is a decreasing integer complexity for a useful part of the arbitrary
DT-tree problem. An overlap at the outer endpoint of a minimizing factor
gives strict descent, even when the two factors are initially far apart.
A separate two-contact move can create a removable factor which was not
minimizing. Arbitrary external trees in a two-factor calculation can be
compressed to five positive integer weights.

These statements do **not** yet prove extraction from every full-alternating
tuple. The remaining assertion is escape from a positive-complexity plateau.
The [overlap-one theorem](PAIRED_HURWITZ_LEAF_EXTRACTION.md#theorem-5-every-full-alternating-dt-tree-reaches-overlap-one)
does not locate its overlap at the required outer endpoint, or control this
complexity during its braid. The mixed three-cycle problem remains separate.

## Domain and literature comparison

Let \(t=(g_1,\ldots,g_r)\) consist of DTs on \(N=2r+1\) sheets,
with full-cycle product. Expand every DT into its two disjoint edges.
The resulting graph \(T\) is a tree: the product makes it connected,
and it has exactly \(N-1\) edge occurrences. A connected graph needs at
least \(N-1\) distinct edges, so no occurrence is repeated.

The descent lemmas below require only this domain, not full alternating
monodromy. The proposed global application imposes \(\langle t\rangle=A_N\)
and \(r\ge7\), as in the
[removable-factor obligation](UNIFORM_TWO_LEAF_DT_ABSORPTION.md#corollary-3-recursive-extensions-and-the-buffered-families).

Two primary results help delimit the literature interface.
[Ripoll, arXiv:0903.3604v3, Definition 0.3 and Theorem 0.4](https://arxiv.org/html/0903.3604v3)
classifies primitive block factorizations of a Coxeter element: all but one
block must be reflections, and conjugacy of the long factors characterizes
the orbits. In \(S_N\), the present shape is \(2^r\), with every DT a
nonreflection. Splitting the other DTs into reflections changes the allowed
braid action; re-pairing those reflections is precisely an additional
obligation. Ripoll's Theorem 5.5 identifies general block-factorization
orbits with connected components of a stratum, but does not prove that
the required stratum or its DT sublocus is connected.
[Liu--Osserman, arXiv:math/0609118v1, Proposition 1.1, Theorem 1.2 and Section 2](https://arxiv.org/html/math/0609118v1)
requires every branch permutation to be a single cycle. A DT has two
ramified points over its branch value and fails that hypothesis. These
are inspected comparisons, not imported extraction theorems or a claim
that no applicable theorem exists elsewhere.

## Lemma 1. Terminal mass detects exactly the removable locus

For a factor \(A\), delete its two edges from \(T\). There are three
components. The **middle component** meets both deleted edges; the other
two, \(S_A,U_A\), are terminal. Orient the edges as

\[
 A=(a\,b)(c\,d),\qquad
 a\in S_A,\quad b,c\in M_A,\quad d\in U_A .
\]

Define

\[
 \lambda_T(A)=|S_A|+|U_A|-2=N-|M_A|-2,\qquad
 \mu(t)=\min_{A\in t}\lambda_T(A).                    \tag{1}
\]

Then \(0\le\mu\le N-4\), and \(\mu=0\) if and only if some DT
has one exclusive leaf on each transposition.

**Proof.** The two edges of a DT are disjoint, so its middle component
has at least two vertices and each terminal component has at least one.
Equality \(\lambda_T(A)=0\) means both terminal components are singletons.
Their vertices are exactly the two exclusive leaves. The converse is
the same tree deletion. This also proves the bounds. \(\square\)

The minimum in (1) is recomputed after every move. Neither the minimizing
factor nor either prospective leaf is fixed. Even if another factor
already has two leaves, \(\lambda_T(A)\) for a chosen \(A\) can be large.

## Lemma 2. Bringing another factor preserves a chosen cut

Fix \(A\) and another factor \(B\). One can bring \(B\) next to \(A\),
keeping both literal permutations unchanged, in \(|i-j|-1\) moves when
their original positions are \(i,j\). Throughout, the three component
sets \(S_A,M_A,U_A\) remain unchanged.

**Proof.** Every edge other than the two edges of \(A\) lies inside one
of these components. Thus every other factor, including \(B\), preserves
all three sets. If \(B\) is on the right, move it left by positive moves
\(Q(D,B)=(B,BDB)\). If it is on the left, move it right by negative
moves \(Q^{-1}(B,D)=(BDB,B)\). Every changed edge still lies inside
its original component. The unchanged \(A\) supplies the only edges
between the three sets. The graph after each move is still a tree,
so each of these sets must still be connected after deleting \(A\).
They are therefore exactly its component sets. \(\square\)

Moving \(A\) unchanged through the intervening factors would not justify
this proof: the conjugating permutation would then cross the cut.

## Lemma 3. An outer one-contact gives strict descent

Suppose \(A\) minimizes (1), \(\mu(t)>0\), and another factor \(B\)
meets its support at exactly one label, which is an outer endpoint of
\(A\). Then a braid of length at most \(|i-j|\le r-1\) reaches a
strictly smaller value of \(\mu\). The word can be stopped at its first
decrease, so every preceding step has the original value of \(\mu\).

**Proof.** Relabel the common label as \(a\), and write

\[
 A=(a\,b)(c\,d),\qquad B=(a\,e)(f\,g),                \tag{2}
\]

with seven distinct labels, \(a,e\in S_A\), and \(d\in U_A\).
Bring \(B\) adjacent by Lemma 2. At every intermediate tuple,
\(\lambda(A)\) retains its original value, hence \(\mu\) cannot
exceed its original value. Stop immediately if it decreases.

Otherwise use the appropriate one-move replacement

\[
 (A,B)\xrightarrow{Q^{-1}}(ABA,A),\qquad
 (B,A)\xrightarrow{Q}(A,ABA),\quad
 ABA=(b\,e)(f\,g).                                    \tag{3}
\]

The only uncolored edge change is \(ae\mapsto be\).
In the tree after bringing, let \(S_a,S_e\) be the two pieces of
\(S_A\) after deleting \(ae\),
with \(a\in S_a,e\in S_e\). After (3), the terminal components of
the unchanged \(A\) are precisely \(S_a,U_A\); its middle component
is \(M_A\cup S_e\). Consequently

\[
             \lambda_{T'}(A)=\lambda_T(A)-|S_e|<\mu(t).
                                                               \tag{4}
\]

This proves strict descent and the word bound. \(\square\)

The statement controls arbitrary external edges: their incidences need
not fix the middle labels or be individually evacuated. The cut argument,
rather than a support assumption on those edges, controls the bringing.

## Lemma 4. Two degree-two contacts can create two new leaves

For six distinct labels suppose adjacent factors are

\[
             A=(a\,b)(c\,d),\qquad B=(a\,e)(c\,f),    \tag{5}
\]

and every other factor fixes \(a,c\). Then one inverse move gives

\[
           (A,B)\longmapsto ((b\,e)(d\,f),A).          \tag{6}
\]

Thus \(a,c\), both initially of tree degree two, become exclusive leaves
of \(A\). For the reversed adjacent order use a positive move. Deleting
\(A\) and those leaves leaves a full-cycle tuple on \(N-2\) sheets.

**Proof.** Formula (6) is conjugation by the involution \(A\). The other
entries do not change. Deleting two leaves leaves a tree. A product of
the edge transpositions of any tree, each used once in any order, is a
full cycle: delete a leaf inductively, and its transposition inserts it
into the lower cycle, regardless of its position in the word. This
proves the asserted lower product. \(\square\)

This is an adjacent-pair statement. Bringing such a pair together can
spoil the degree-two condition, so Lemma 2 is not an unconditional way
to arrange (5). A factor captured by (6) need not have minimized (1).

## Lemma 5. Five weighted components control every adjacent-pair rotation

Fix any adjacent \(A,B\). Remove all four of their edges from \(T\),
leaving five component sets \(V_1,\ldots,V_5\), of weights
\(w_i=|V_i|>0\). Remember the component containing each of the at most
eight labels in \(\operatorname{supp}(A)\cup\operatorname{supp}(B)\).

For every integer \(k\), replace the pair by \(Q^k(A,B)=(A_k,B_k)\).
The four new edges join these same five sets as a tree. The masses of
\(A_k,B_k\) are obtained exactly by deleting their respective two edges
in this five-vertex weighted tree, adding the weights in its two terminal
components, and subtracting two. In particular each is a linear form
in the five weights minus two.

It suffices to inspect at most six rotations: if \(C=AB\), then

\[
 Q^k(A,B)=(AC^k,AC^{k+1}),\qquad
 \operatorname{ord}(C)\in\{2,3,4,5,6\}.                \tag{7}
\]

**Proof.** Hurwitz moves retain the two DT conjugacy classes and the
full product. Every changed endpoint stays in the original pair support.
All other edges are unchanged and connect each \(V_i\). Since the whole
graph is still a tree, the four new edges cannot lie within one \(V_i\),
and their quotient is a tree on these five components. Deleting either
factor commutes with this contraction, including identification of the
middle component. This proves the weighted formula.

For (7), \(ACA=C^{-1}\) gives the displayed identity by induction.
The union of the original four edges is a forest of alternating paths.
Its path-length partitions of four are
\(4,3+1,2+2,2+1+1,1+1+1+1\). The corresponding product orders are
\(5,4,3,6,2\). Thus (7) returns to the original pair after at most
six steps. \(\square\)

This is a complete compression of the external interactions for the
two rotating factors. It does not compute the changed masses of every
outside factor, or prove that some pair must descend. If a predicted
pair mass is smaller than the current \(\mu\), it certifies strict
descent. For descent through a plateau one must also check the
intermediate rotations; an endpoint inequality alone does not give
a nonincreasing path.

## An input already past the overlap gate can need a different candidate

Consider this degree-fifteen tuple:

\[
\begin{split}
 (&(67)(10\,11),\ (34)(9\,10),\ (23)(89),\ (12)(56),\\
  &(11\,12)(13\,14),\ (45)(78),\ (01)(12\,13)).        \tag{8}
\end{split}
\]

Its tree is the path \(0,1,\ldots,14\), so its product is a 15-cycle.
Its first two supports already overlap in one label. Their product
squared is \((9\,10\,11)\). The group is primitive: an odd-degree
nontrivial block has size at least three, and a DT cannot move two such
blocks; transitivity therefore excludes a nontrivial block system.
Jordan's 3-cycle theorem gives \(A_{15}\), using the precise formulation
already retained in the
[Jones source ledger](alternating_hurwitz_source_ledger.json).

The factor masses in order are \((9,7,7,9,11,10,1)\). The only minimizer
is the last factor, whose outer endpoints are \(0,13\). No other factor
has a one-label intersection at either endpoint. Thus the hypothesis
of Lemma 3 does not follow immediately from an adjacent overlap-one pair,
even in the residual degree-fifteen domain.

Nevertheless \(Q_2^{-1}\) applies Lemma 4 to the second and third factors.
They become \((24)(8\,10),(34)(9\,10)\). Sheets \(3,9\) are now leaves
of the latter, so \(\mu\) drops from one to zero. Neither sheet was
initially a leaf, and this factor was not the original minimizer.
This is a positive extraction certificate, not an orbit obstruction.

There is also a genuine plateau for every power of a single initially
adjacent pair. Take

\[
\begin{split}
 (&(12)(13\,14),\ (67)(10\,11),\ (34)(11\,12),\\
  &(45)(9\,10),\ (89)(12\,13),\ (01)(23),\ (56)(78)).    \tag{9}
\end{split}
\]

Again the tree is the path \(0,1,\ldots,14\). The second and third
supports already overlap in one label, and their product squared is
\((10\,12\,11)\); the same primitivity and Jordan argument gives
\(A_{15}\). Its factor masses are \((1,9,5,8,9,11,11)\), and there
is no minimizing outer one-contact.

The periods of \(Q_i\) on (9), for \(i=1,\ldots,6\), are
\((2,6,6,6,2,2)\). Every state in each of these six periodic orbits
has mass exactly one. Thus no word supported on just one of these
adjacent pairs can lower \(\mu\), with no word-length restriction.
These are 24 pair/power cases in the retained certificate, including
the initial state for each pair. Lemma 5 reduces all integer powers to
this complete finite list.

In contrast, the word \((-3,-4)\) has mass trajectory \((1,1,0)\).
Its affected three-factor window is exactly

\[
\begin{split}
 &( (34)(11\,12),\ (45)(9\,10),\ (89)(12\,13))\\
 &\quad\longmapsto
 ((35)(9\,10),\ (34)(11\,12),\ (89)(12\,13))\\
 &\quad\longmapsto
 ((35)(9\,10),\ (89)(11\,13),\ (34)(11\,12)).
\end{split}                                                     \tag{10}
\]

The last displayed factor now has exclusive leaves \(4,12\), both
initially of degree two. Its original mass was five. The first move
keeps the global mass at one while reducing this candidate's mass to
two; the second prunes a terminal branch of size two and makes its
mass zero. This is the tree-pruning calculation (3)--(4) at a
nonminimizing factor, certified by the weighted criterion of Lemma 5.
It demonstrates a need for a flat step and a different candidate when
using these local moves, rather than an obstruction to general descent.

## The remaining plateau-escape obligation

A sufficient global statement is the following, currently **open**.
For every full-\(A_{2r+1}\) tuple, \(r\ge7\), of positive mass \(m\),
a braid on which \(\mu\) stays equal to \(m\) reaches either:

1. a minimizing factor with the outer one-contact of Lemma 3; or
2. an adjacent pair having a rotation whose intermediate masses are at
   most \(m\) and whose final rotating factor has mass less than \(m\).

Under this statement, Lemmas 3 and 5 give a strict decrease, stopping at
the first decrease if necessary. Repeat with a newly chosen minimizer.
There are at most the original \(\mu\le N-4\) strict decreases before
Lemma 1 supplies a removable DT. No fixed sheet is required to remain a
leaf. The lower-degree induction and arbitrary two-anchor restoration
then use the separately proved
[two-leaf extension theorem](UNIFORM_TWO_LEAF_DT_ABSORPTION.md).

The plateau assertion is a sufficient strategy, not a proved consequence
of the arbitrary-input overlap gate, and not asserted equivalent to
unrestricted extraction. A braid to extraction might increase \(\mu\).
A positive global orbit minimum would, however, forbid every Lemma 3
contact and every Lemma 5 endpoint decrease in every minimizing
representative. In particular each factor touching an outer endpoint
of a minimizer must also touch at least one other support label of that
minimizer. These coupled contacts are the precise next structural target.

## Exact evidence

The [checker](verify_terminal_mass_hurwitz_descent.py) reuses the existing
complete 392-tuple degree-seven domain only to test the new mechanisms.
It checks cut-preserving bringing in both orders, the exact terminal-size
drop, two-contact capture, all pair rotations by the independent weighted
quotient, the proper order-168 controls, cyclic higher-degree controls,
and the literal degree-fifteen certificate (8). Deliberately invalid
trees and corrupted capture hypotheses are rejected. It also replays
all 24 initial pair/power cases and the two-move escape (9)--(10).

The [receipt](terminal_mass_hurwitz_descent_verification.json) and
[source ledger](terminal_mass_hurwitz_descent_source_ledger.json) record
the replay and exact inputs. No new degree census, universal plateau
escape, arbitrary-input extraction, mixed connectivity, formal
verification, independent replay, external review or priority claim is
asserted.
