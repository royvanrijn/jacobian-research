# Terminal mass and Hurwitz descent with changing candidate leaves

There is a decreasing integer complexity for a useful part of the arbitrary
DT-tree problem. An overlap at the outer endpoint of a minimizing factor
gives strict descent, even when the two factors are initially far apart.
A separate two-contact move can create a removable factor which was not
minimizing. Arbitrary external trees compress to five weights for two
factors and seven weights for three factors, retaining the distinct
endpoint labels and their component memberships. A complete three-factor
calculation below gives eight inequalities valid at every weighting of
the plateau example's local pattern.

These statements do **not** yet prove extraction from every full-alternating
tuple. The principal remaining objective is to exclude a positive minimum
over an entire braid orbit. This needs no monotonicity along a braid.
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

## Lemma 6. Seven weighted components control the complete three-factor orbit

Fix three consecutive DTs \((A,B,C)\). Delete their six edges, leaving
seven components \(V_0,\ldots,V_6\), with weights \(w_i=|V_i|\).
Let \(L\) be their original support, \(s=|L|\le12\), and retain both
the literal permutations on \(L\) and the owner map
\(\pi:L\longrightarrow\{0,\ldots,6\}\). Then:

1. The entire local \(B_3\)-orbit can be computed on \(L\), without
   any external permutations. It has at most
   \([3\binom{12}{4}]^2=2{,}205{,}225\) states, independently of \(N\).
2. In every local state the six edges project to a tree on the same seven
   components. Each of its three factor masses is exactly
   \[
      \lambda=\sum_{i=0}^6 c_i w_i-2,\qquad c_i\in\{0,1\},             \tag{11}
   \]
   where \(c_i=1\) precisely on the two terminal pieces after deleting
   that factor's two projected edges.
3. Write \(\mathcal C\) for the distinct coefficient vectors in this
   complete orbit. Delete any vector strictly containing another vector
   coordinatewise, leaving \(\mathcal C_{\min}\). For all positive weights,
   \[
       \min_{\text{local states and their factors}}\lambda
       =\min_{c\in\mathcal C_{\min}}c\cdot w-2.                       \tag{12}
   \]

**Proof.** Every local Hurwitz move conjugates by a permutation supported
on \(L\). The support union, three DT classes and product \(ABC\) stay
fixed. There are \(D_s=3\binom{s}{4}\) DTs on \(L\); with product fixed,
the first two entries determine the third, giving the bound \(D_s^2\).
Breadth-first exploration under \(Q_1^{\pm1},Q_2^{\pm1}\), stopped only
at closure, therefore computes the complete orbit.

All external edges remain unchanged and connect each literal \(V_i\).
The full product and number of edge occurrences remain unchanged, so
the new full graph is again a tree, by the argument in the domain
definition. Its six changed edges project to a tree on the seven
components. Deleting a factor commutes with contracting the unchanged
components, proving (11). Positivity of \(w_i\) proves every strict
coordinatewise domination and (12). \(\square\)

There are finitely many combinatorial templates: an ordered DT triple
on at most twelve support labels, together with a surjective owner map
onto seven vertices whose six projected edges form a tree. Every actual
window supplies one. Its weights satisfy
\(w_i\ge |\pi^{-1}(i)|\), and \(\sum_iw_i=N\) is odd. Extra conditions
may be required to realize its external edges in DT pairs. A template
with positive integer weights is not by itself a full-alternating tuple.

The owner map does **not** permit replacing every label by its component
vertex before computing the braid. Several distinct labels may belong
to one component, and conjugation must distinguish them. The weights
and endpoint incidence data, rather than seven fictitious sheets, are
the compressed input. No external factor masses are computed by (11);
a local candidate of mass below the current global minimum is already
sufficient to certify a decrease.

## Proposition 7. The plateau window has eight uniform inequalities

Consider the following endpoint template, with nine distinct labels
\(x_0,\ldots,x_8\):

\[
 A=(x_0\,x_1)(x_6\,x_7),\quad
 B=(x_1\,x_2)(x_4\,x_5),\quad
 C=(x_3\,x_4)(x_7\,x_8),                              \tag{13}
\]

and owners

\[
 (\pi(x_0),\ldots,\pi(x_8))=(0,1,2,2,3,4,4,5,6).      \tag{14}
\]

The quotient tree is the path \(0,1,\ldots,6\), whose edge colors are
\(A,B,C,B,A,C\). In particular \(x_2,x_3\) and \(x_5,x_6\) are distinct
labels joined through external trees, not equal labels. The weight
lower bounds are \((1,1,2,1,2,1,1)\).

The complete local orbit has **162 states**, giving 27 distinct
coefficient vectors. Its minimal coefficient supports are exactly

\[
 \{0,3\},\ \{0,5\},\ \{0,6\},\ \{1,3\},\ \{1,5\},\ \{1,6\},
 \ \{3,5\},\ \{3,6\}.                                \tag{15}
\]

Consequently its smallest candidate mass, for every admissible weighting
with no upper bound, is

\[
 \min(\alpha+\beta,\alpha+\gamma,\beta+\gamma)-2,
 \quad \alpha=\min(w_0,w_1),\quad\beta=w_3,\quad
 \gamma=\min(w_5,w_6).                                \tag{16}
\]

The state count also has a direct proof. On each of the three disjoint
two-edge support paths, one factor is the identity and the other two
are transpositions with a fixed three-cycle product. There are three
possible ordered transposition pairs for that product. The three
identity positions must be distinct: each factor acts on exactly two
paths, and conjugation within the original support components preserves
this property. Thus there are at most \(3!\,3^3=162\) states.
Initially \(Q_1^2\) cycles the three choices on the \(AB\) path and
fixes the other two path restrictions; \(Q_2^2\) does the same on the
\(BC\) path. The word \(Q_2Q_1^2Q_2^{-1}\), read left to right,
cycles only the \(AC\) path. These three actions independently reach
all 27 choices with the identity positions fixed. Ordinary generators
realize every permutation of those positions. The upper bound is
therefore attained.

**Exact finite calculation of the coefficient vectors.** The
[three-factor checker](verify_terminal_mass_three_factor.py) starts
from the literal nine-label input (13)--(14), uses endpoint substitution
to explore all four signed generators to closure, and checks every
transition separately against permutation conjugation. It checks all
486 candidate coefficient vectors, retains their 27 distinct values,
and verifies that the eight sets (15) occur and that every other set
strictly contains one of them. It retains a braid word and factor index
for every set in (15), and a deterministic digest of the complete
sorted orbit, in the
[receipt](terminal_mass_three_factor_verification.json). This proves
the finite assertion and, by Lemma 6 and positivity, (16) for all the
stated weights. The calculation enumerates local permutations, not
values of \(N\) or possible weights.

There is a uniform exclusion for a minimizing factor. The original
three terminal supports and some strictly smaller attainable supports
are

| Original factor | Original support in (11) | Attainable proper subset | Strict drop |
|---|---|---|---|
| \(A\) | \(\{0,5,6\}\) | \(\{0,5\}\) | \(w_6\) |
| \(B\) | \(\{0,1,4,5,6\}\) | \(\{0,5\}\) | \(w_1+w_4+w_6\) |
| \(C\) | \(\{0,1,2,6\}\) | \(\{0,6\}\) | \(w_1+w_2\) |

Thus none of the three original factors of this template can attain the
global mass of an orbit-minimizing representative. The smaller candidate
can be any factor of a new local state; it need not be the same literal
factor as the original one. This exclusion is uniform in all external
trees with this owner data. It does not assert that every minimizing
factor participates in this template.

For (9), the window starting at the third factor is (13) with
\[
 (x_0,\ldots,x_8)=(3,4,5,8,9,10,11,12,13).
\]
The seven components, in order, are
\[
 \{0,1,2,3\},\ \{4\},\ \{5,6,7,8\},\ \{9\},\ \{10,11\},
 \ \{12\},\ \{13,14\},
\]
with weights \((4,1,4,1,2,1,2)\). Formula (16) gives zero.
The local word \((-1,-2)\), or global word \((-3,-4)\), realizes
\(w_1+w_5-2=0\), with terminal components \(\{4\},\{12\}\).
The original masses in this window were \((5,8,9)\); the global
minimizer of mass one was elsewhere. No assumption that the window
contains the initial minimizer was needed to certify this extraction.

All five initially consecutive windows of (9) have been computed to
closure. Their orbit sizes are \((18,144,162,18,18)\), with smallest
candidate masses \((1,0,0,1,1)\). Thus two windows rule out this tuple
as an orbit minimizer even though every initial adjacent-pair power
is flat. Identifying the distinct endpoint labels in (14) before
braiding instead gives 147 states, illustrating why that substitution
would compute a different problem.

The support forest in (13) consists of three disjoint two-edge paths,
with color pairs \(AB,BC,CA\). Keeping that forest fixed, the calculation
also exhausts **all 189 partitions of its nine labels into seven owner
components for which the projected edges form a tree**. Owner vertex
names are identified, while the nine support labels remain literal.
To enumerate them, identify two pairs of labels and retain exactly the
partitions with tree quotient. Every partition of nine labels into
seven blocks has either one triple block or two double blocks, so this
enumeration is complete. Equivalently, a tree quotient is obtained by
identifying one label from each of the three path components (27
choices), or by joining the three path components through two distinct
pair identifications (162 choices).

Every owner partition has 27 distinct coefficient vectors; 162 have
eight minimal vectors and 27 have twelve. The number of original
factors that can minimize the local candidate mass is distributed as
follows:

| Number of possible original minimizing factors | Owner partitions |
|---|---|
| 0 | 44 |
| 1 | 111 |
| 2 | 33 |
| 3 | 1 |

Here “can minimize” refers exactly to feasibility of the isolated
inequalities (18)--(19), with positive mass, the endpoint-label weight
lower bounds and odd total weight. This feasibility has an elementary
test. An original coefficient vector containing a smaller attainable
vector cannot minimize, by positivity. If it contains no smaller one,
give its selected components weights \(|\pi^{-1}(i)|+1\), and all other
components one common weight larger than their selected total. Every
other distinct coefficient vector then has larger value. Increase one
unselected weight by one if necessary to make the total odd. The
selected mass is positive. The checker verifies these witnesses.

Thus **387 of the 567 marked owner templates are uniformly excluded**;
the other **180** have positive isolated linear witnesses. All owner
maps, initial vectors, minimal vectors and possible original minimizing
factors are retained in the receipt. This exhausts one support forest,
not every possible three-factor support forest. The surviving linear
witnesses are not claimed to extend to global orbit minima or even to
full-alternating DT completions.

## Proposition 8. Orbit minima replace monotone descent

For a full-cycle DT tuple on fixed \(N\), let
\[
                  m(\mathcal O)=\min_{u\in\mathcal O}\mu(u),          \tag{17}
\]
where \(\mathcal O\) is its entire braid orbit. Then extraction of a
removable factor somewhere in \(\mathcal O\) is equivalent to
\(m(\mathcal O)=0\). If \(u\) attains a positive orbit minimum \(m\),
every three-factor template arising at any minimizing representative
must satisfy the finite system
\[
       c\cdot w\ge m+2\quad(c\in\mathcal C_{\min}).                  \tag{18}
\]
If a factor in the initial window minimizes \(\mu(u)\), its initial
coefficient vector \(c_A\) additionally satisfies
\[
                         c_A\cdot w=m+2.                           \tag{19}
\]
A violation of any endpoint inequality contradicts minimality,
irrespective of intermediate masses along its realizing braid.

**Proof.** The fixed-degree tuple space is finite and \(\mu\) is a
nonnegative integer, so (17) is attained. Lemma 1 proves the equivalence.
Every local braid endpoint belongs to \(\mathcal O\), hence each of its
factor masses is at least \(m\). Lemma 6 gives (18), and the assumed
initial minimizing factor gives (19). \(\square\)

One may also gather any three selected factors, including a minimizing
\(A\), into a consecutive block while keeping all three literal factors
and the cut of \(A\). Keep \(A\) stationary. Bring the selected factors
from the right leftward by positive moves, nearest selected factor
first; bring those from the left rightward by negative moves, again
nearest first. Stop beside the selected block and never cross a
selected factor. Every crossed factor preserves the cut of \(A\), as
in Lemma 2, so each changed external edge stays within its old component.
Consequently \(\lambda(A)=m\) is retained. At an orbit minimum the
gathering states also have \(\mu=m\). The seven-component owner data
must be computed **after** this gathering: its other external cuts
need not agree with those of the original separated triple.

The principal open objective is therefore to show that **every
hypothetical positive global orbit minimum in the full-alternating
residual domain has a minimizing representative and a three-factor
window violating (18)**, or to find another contradiction to a
positive minimum. That would prove extraction directly by (17),
without a globally nonincreasing extraction algorithm. It is a
sufficient three-factor strategy; no equivalence with this restricted
local strategy is asserted.

There are finitely many endpoint templates and finitely many linear
inequalities per template, but their mere finiteness does not prove
infeasibility of (18)--(19). The missing global step is to force a
violating template or incompatible systems across windows, using the
tree and full-alternating hypotheses. Feasible weights for an isolated
template need not extend to any actual global tuple, and full alternating
monodromy is not encoded by seven weights. Conversely, a local candidate
test may miss a decrease in an outside factor. The positive local
minimum in the proper order-168 degree-seven control shows that an
unqualified assertion about all DT trees would be false.

## The remaining plateau-escape obligation

A stronger monotone strategy is the following, still **open** and no
longer needed as the principal proof objective.
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
contact and every Lemma 6 endpoint decrease in every minimizing
representative, as made explicit in Proposition 8. In particular each
factor touching an outer endpoint
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

The [three-factor checker](verify_terminal_mass_three_factor.py) and its
[receipt](terminal_mass_three_factor_verification.json) give the
complete local orbits above, 1,245 comparisons of compressed candidate
masses against full-tree cuts, a twelve-distinct-label support control,
and a proper order-168 control with positive local minimum. Corrupted
DTs, insufficient component weights, quotient loops and an incomplete
orbit are detected. No degree census is extended. All 189 owner
partitions of the specified nine-label forest are computed, including
387 uniform exclusions and 180 positive isolated linear witnesses.
The inventory for other support forests and a global contradiction to
(18)--(19) remain open.

The [receipt](terminal_mass_hurwitz_descent_verification.json) and
[source ledger](terminal_mass_hurwitz_descent_source_ledger.json) record
the replays and exact inputs. No new degree census, exclusion of every
positive global orbit minimum, universal plateau escape,
arbitrary-input extraction, mixed connectivity, formal
verification, independent replay, external review or priority claim is
asserted.
