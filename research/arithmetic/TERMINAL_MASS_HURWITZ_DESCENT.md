# Terminal mass and Hurwitz descent with changing candidate leaves

There is a decreasing integer complexity for a useful part of the arbitrary
DT-tree problem. An overlap at the outer endpoint of a minimizing factor
gives strict descent, even when the two factors are initially far apart.
A separate two-contact move can create a removable factor which was not
minimizing. Arbitrary external trees compress to five weights for two
factors and seven weights for three factors, retaining the distinct
endpoint labels and their component memberships. A complete three-factor
calculation below gives eight inequalities valid at every weighting of
the plateau example's local pattern. Complementary full twists impose
simultaneous cyclic contact conditions on a minimizing factor's three
cut components. Exact phase sets and reapplication of the cut-size
restriction to new minimizing factors exclude the residual ratio three
for equal terminals, without controlling middle-middle attachments.
For a singleton terminal with equal middle and other-terminal sizes,
minimality forces two paths and their ordered edge pairing. Commuting
swaps then expose a two-contact leaf capture, excluding this cut for
every size at least three.

An elementary comparison with all the other factor masses excludes every
terminal-to-terminal pair at a minimizing cut. In particular the equal-terminal
family in Corollary 16 is impossible in the residual degree range. Theorem 17
leaves only singleton cuts of sizes \((1,t;2t)\), with even \(t\), so a
positive minimum requires \(N\equiv1\pmod6\) and equality in \(3m\ge N-4\).
Lemma 18 identifies the antipodal middle coordinates and the exact projected
terminal tree in this last family. Lemma 19 uses the fact that this projection
is fixed under every complementary braid: incident terminal edges must have
the same middle lift, and conjugating a middle-middle factor then contradicts
this common lift. Theorem 20 excludes the last family and proves removable
factor extraction in every odd degree at least fifteen.

This proof uses a minimum over the entire braid orbit and does not require
monotonicity along a braid. It also supplies the extraction input to the
existing pure-DT connectivity induction. The mixed three-cycle problem
remains separate.
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

One sufficient local strategy would show that every hypothetical positive
global orbit minimum in the full-alternating residual domain has a minimizing
representative and a three-factor window violating (18). Theorem 20 instead
excludes positive minima using complementary braid stability. A direct
three-factor-window exclusion remains a separate optional strategy; no
equivalence with that restricted local strategy is asserted.

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

## Lemma 9. Complementary twists develop the three cuts simultaneously

Suppose a tuple attains a positive global orbit minimum \(m\), and
\(A\) is a minimizing factor. There is a representative with \(A\)
in the first position, possibly conjugated, with its three component
sets and mass unchanged. Write it as
\[
 (A,B_1,\ldots,B_{r-1}),\qquad
 A=(a\,b)(c\,d),\quad
 a\in S,\ b,c\in M,\ d\in U,
 \qquad R=B_1\cdots B_{r-1}.                         \tag{20}
\]
The restrictions \(R_S,R_M,R_U\) are full cycles on their respective
sets, including the identity cycle on a singleton. For every integer
\(k\), the tuple
\[
 (A,R^{-k}B_1R^k,\ldots,R^{-k}B_{r-1}R^k)             \tag{21}
\]
is in the same braid orbit, has these same cuts for \(A\), and has
\(\mu=m\). None of its other factors has an outer one-contact with
\(A\).

**Proof.** Move the selected factor left by inverse moves
\((D,A)\mapsto(DAD,D)\). Each crossed \(D\) is unchanged and preserves
the three component sets of the current cut. The new edges of \(A\)
therefore still connect \(S\) to \(M\) and \(M\) to \(U\). The
unchanged complement connects each set, so these remain precisely the
three components after deleting the new \(A\). Its mass remains \(m\).
This placement conjugates \(A\); it differs from the literal-factor
bringing of Lemma 2.

The complement edges in each set form a tree. Their product uses each
edge once, so the leaf-insertion proof gives a full cycle there.
The ordinary full twist on the last \(r-1\) strands simultaneously
conjugates their entries by \(R^{-1}\). Explicitly its word is
\((Q_2Q_3\cdots Q_{r-1})^{r-1}\); for one complementary strand it
is the empty word. This follows by cycling the first complementary
entry to the end \(r-1\) times. It proves (21), also for negative
powers. Conjugation by \(R\) preserves each component set and its
connected complement tree, so the cut and \(\lambda(A)=m\) persist.
Global minimality forces \(\mu=m\), even if individual intermediate
moves in this full twist might have behaved differently without that
hypothesis. Lemma 3 excludes every asserted outer contact. \(\square\)

## Theorem 10. Cyclic compatibility excludes closed terminal trees

Use (20), and put \(s=|S|\), \(u=|U|\), \(\ell=|M|\). At a positive
global orbit minimum the following restrictions hold simultaneously.

1. No complementary DT has both its edges in \(S\), or both in \(U\).
2. A DT pairing an \(S\)-edge with a \(U\)-edge forces \(s=u\).
   Moreover all such pairs agree with the **same** bijection
   \[
       f:S\longrightarrow U,\qquad f(R_S^k a)=R_U^k d.               \tag{22}
   \]
   Its two edges are \(e\) and \(f(e)\).
3. At least one DT pairs a terminal edge with a middle edge.
   Every DT pairing an \(S\)-edge with an \(M\)-edge satisfies
   \[
          \ell\le2s,\qquad \ell/\gcd(s,\ell)\le4.                 \tag{23}
   \]
   The corresponding inequalities with \(u\) hold for every
   \(U\)-to-\(M\) pair.

These restrictions require only the polynomial DT tree and global
minimality, not full alternating monodromy.

**Proof.** For a fixed complementary factor \(B\), define the periodic
sets of integers
\[
 \begin{split}
 E_S(B)&=\{k:R_S^k a\in\operatorname{supp}(B)\},\\
 E_U(B)&=\{k:R_U^k d\in\operatorname{supp}(B)\},\\
 E_M(B)&=\{k:R_M^k b\text{ or }R_M^k c\in\operatorname{supp}(B)\}.
 \end{split}                                                       \tag{24}
\]
They record intersections with the **fixed** support of \(A\) in
(21). Their common phase \(k\), and the three cycles defining it,
come from one complementary tuple. They are not independently chosen
weights or witnesses for separate windows.

Every complementary edge lies wholly inside one component. If both
edges of \(B\) lie in \(S\), then \(E_U=E_M=\varnothing\), whereas
\(E_S\ne\varnothing\) because \(R_S\) is a full cycle. At such a
phase the intersection with \(A\) is precisely \(\{a\}\), contrary
to Lemma 9. This proves assertion 1, also for \(U\).

For an \(S\)-to-\(U\) pair, \(E_M=\varnothing\). Neither outer
endpoint may occur alone, so \(E_S=E_U\). Their densities in a common
period are \(2/s\) and \(2/u\), respectively. Thus \(s=u\).
Equality of these sets modulo \(s\) says exactly that (22) carries
the two endpoints of the \(S\)-edge to those of the \(U\)-edge.
The roots and cycles in (22) are fixed for all factors, proving the
simultaneous, rather than pairwise independently chosen, matching.

For an \(S\)-to-\(M\) pair, the absence of an outer one-contact says
\(E_S\subseteq E_M\). The former has density \(2/s\). The latter
is a union of at most four residue classes modulo \(\ell\): two
middle endpoints can meet either \(b\) or \(c\). Hence
\(2/s\le4/\ell\). In addition, even one residue class modulo \(s\)
in \(E_S\) runs through \(\ell/\gcd(s,\ell)\) distinct residues
modulo \(\ell\). They must all belong to \(E_M\), proving (23).

It remains to exclude a closed pair of terminal trees, in which no
terminal edge is paired with a middle edge. By assertion 1, all
\(s-1\) edges of \(T_S\) are then paired with all \(u-1\) edges
of \(T_U\). If both sizes are one then \(m=0\), already a
contradiction. Otherwise there is a paired edge, so assertion 2
identifies the two edge-colored trees by \(f\), with \(f(a)=d\).
A finite tree of size at least two has a leaf \(v\ne a\). Then
\(f(v)\ne d\) is also a leaf, and their two incident edges belong
to one DT. Only \(A\) connects these terminal trees to the middle,
at \(a,d\). Thus \(v,f(v)\) are exclusive leaves in the full tree
on different edges of that DT. Lemma 1 gives \(\mu=0\), the final
contradiction. \(\square\)

## Corollary 11. Positive minima have restricted mass and terminal sizes

Every positive global orbit minimum of a polynomial DT tuple satisfies
\[
                         3m\ge N-4.                               \tag{25}
\]
At every minimizing factor, unequal nonsingleton terminal sizes are
possible only with
\[
       (\{s,u\},\ell,N)=(\{2,3\},4,9)
          \quad\text{or}\quad(\{3,4\},6,13).                       \tag{26}
\]
Consequently, in the residual degree range \(N\ge15\), the terminal
sizes are equal or one is a singleton. If the sizes are \(1,t\),
\(t>1\), their middle size satisfies one of the four exact conditions
\[
 \begin{array}{c|c}
 \ell & \text{necessary divisibility}\\ \hline
 t & \text{none}\\
 2t & 2\mid t\\
 3t/2 & 4\mid t\\
 4t/3 & 6\mid t .
 \end{array}                                                       \tag{27}
\]
For equal terminal sizes \(s=u>1\), \(\ell\) is odd and
\(\ell/\gcd(s,\ell)\) is either 1 or 3.

**Proof.** Let \(x,y,z,q\) count complementary factors of types
\(SU,SM,UM,MM\), respectively. Theorem 10 exhausts the possible types,
and the tree edge counts give the coupled equalities
\[
 x+y=s-1,\qquad x+z=u-1,\qquad y+z+2q=\ell-1.                      \tag{28}
\]
At least one of \(y,z\) is positive. For its terminal size \(v\),
(23) gives \(\ell\le2v\le2(m+1)\), since \(s+u=m+2\).
Now \(N=\ell+m+2\le3m+4\), proving (25).

If \(s\ne u\), assertion 2 gives \(x=0\). Suppose both sizes exceed
one and put \(k=\min(s,u)\). Then \(y,z>0\), so (23) and (28) give
\[
   m+1\le\ell\le2k\le s+u-1=m+1.
\]
All are equal: the sizes are \(k,k+1\), \(\ell=2k\), and \(q=0\).
Apply the gcd inequality in (23) to the larger terminal component:
\[
              \frac{2k}{\gcd(k+1,2k)}
              =\frac{2k}{\gcd(k+1,2)}\le4.
\]
For \(k\ge2\), this leaves only \(k=2,3\), proving (26).

For sizes \(1,t\), (28) gives \(\ell=t+2q\), while (23) gives
\(\ell\le2t\) and \(\ell/\gcd(t,\ell)\le4\). Write
\(t=ag,\ell=bg\) with \(\gcd(a,b)=1\). The conditions
\(a\le b\le2a\), \(b\le4\) leave exactly
\((a,b)=(1,1),(1,2),(2,3),(3,4)\). Integrality of
\(q=(\ell-t)/2\) gives the divisibilities in (27).
For \(s=u\), (28) gives \(y=z>0\) and odd \(\ell\). Its positive
odd divisor ratio in (23) is therefore 1 or 3. \(\square\)

This is a global compatibility obstruction: a positive isolated
seven-weight witness must admit one common complement, its three
cycles, all covers (24), the shared matching (22), and the edge counts
(28). No assertion is made that these necessary conditions suffice
for a global orbit minimum. In particular no new count of surviving
marked templates is inferred.

The proper order-168 degree-seven orbit has \(m=1\), attaining
\(3m=N-4\); dropping the monodromy hypothesis from the desired complete
exclusion would remain false.

## Lemma 12. Exact three- and four-phase covers align the endpoints

Use cyclic coordinates on the middle component, with \(b=0\) and
\(c=h\ne0\) in \(\mathbb Z/\ell\). For a terminal-middle factor whose
middle edge has distinct endpoints \(p,q\), its middle contact set is
\[
                    E=\{p,q,p-h,q-h\}.                            \tag{29}
\]
Put \(g=\gcd(v,\ell)\), where \(v\) is the size of its terminal
component. Each terminal endpoint phase modulo \(v\) develops into
a whole coset of the order-\(\ell/g\) subgroup modulo \(\ell\).
If \(\ell/g=3\), covering even one such coset forces
\[
 \ell=3g,\qquad h\in\{g,2g\},\qquad q-p\in\{g,2g\},
 \qquad E=p+\{0,g,2g\}.                                           \tag{30}
\]
Both terminal endpoint coordinates are congruent to \(p\) modulo
\(g\). At some common twist phase, this factor has support intersection
three with the minimizer: one outer endpoint and both middle endpoints.
If \(\ell/g=4\), the same cover forces all four points of (29) to
be one coset modulo \(g\), and again both terminal endpoint coordinates
are congruent modulo \(g\).

**Proof.** An integer phase fixed modulo \(v\) runs through a full
coset modulo \(\ell\), by the elementary congruence criterion modulo
\(g\). For an order-three coset, if \(h\not\equiv0\pmod g\), each
of \(\{p,p-h\}\) and \(\{q,q-h\}\) contributes at most one point
to that coset. This cannot cover three points. Thus \(g\mid h\).
Neither two-point set alone covers three points, so both lie in the
required coset. Nonzero \(h\), distinct \(p,q\), and \(\ell=3g\)
now give (30). Every covered coset must equal this same \(E\), proving
the assertion for the two terminal endpoints. The two middle pairs
\(\{0,h\}\) and \(\{p,q\}\) are translates inside this order-three
coset. A phase making them coincide is congruent to either terminal
endpoint phase modulo \(g\). The congruence criterion therefore gives
a common integer phase with an outer contact and both middle contacts.
For order four, a four-element coset contained in the at most four
points of \(E\) equals \(E\); in particular \(p,q,p-h,q-h\) all
have its residue modulo \(g\). \(\square\)

## Corollary 13. Ratio three is excluded in the residual degree range

At a positive global orbit minimum, every minimizing factor with
\(N\ge15\) satisfies one of the following necessary conditions:
\[
 \begin{array}{c|c}
 \text{terminal sizes}&\text{middle size}\\ \hline
 (s,s),\ s>1&\ell\text{ odd and }\ell\mid s\\
 1,t\ (t>1)&\ell=t\text{, or }\ell=2t\text{ with }2\mid t.
 \end{array}                                                       \tag{31}
\]
These conclusions still require no full alternating monodromy assumption.

**Proof.** First suppose both terminals have size \(s>1\) and the
remaining ratio in Corollary 11 is \(\ell/\gcd(s,\ell)=3\).
A terminal-middle factor exists by Theorem 10. Lemma 12 supplies a
complementary twist phase with
\[
                  A=(a\,b)(c\,d),\qquad B=(e\,a)(b\,c).
\]
Bring \(B\) beside \(A\), retaining its literal permutation and the
cut of \(A\). Their four edges form the path \(e,a,b,c,d\).
Delete them and let the five component weights in this order be
\(w_0,\ldots,w_4>0\). Thus
\[
             w_0+w_1=s,\quad w_2+w_3=\ell,\quad w_4=s.
\]
The five pair rotations have the following middle components; the
entries are sets of component indices, not identified sheet labels:
\[
 \begin{array}{c|cc}
 k&\text{first factor of }Q^k(A,B)&\text{second factor}\\ \hline
 0&\{2,3\}&\{1,2\}\\
 1&\{0,3\}&\{2,3\}\\
 2&\{0,4\}&\{0,3\}\\
 3&\{1,4\}&\{0,4\}\\
 4&\{1,2\}&\{1,4\}.
 \end{array}                                                       \tag{32}
\]
This table follows directly by conjugating the two displayed DTs.
By Lemma 5 the factor mass is \(N\) minus its middle weight minus
two. Global minimality therefore bounds every middle weight in (32)
by the original middle weight \(\ell\). In particular
\[
          \ell\ge s+\max(w_0,w_1)\ge3s/2.                        \tag{33}
\]
Write \(s=ag,\ell=3g\), with \(\gcd(a,3)=1\). The old bound
\(\ell\le2s\) and (33) force \(a=2\). Since \(\ell\) is odd,
\(g\) is odd, \(N=7g\), and \(N\ge15\) gives \(g\ge3\).
Equality in (33) gives
\[
 w_0=w_1=g,\quad w_4=2g,\quad w_2+w_3=3g.
\]
The other two middle weights in (32) give \(g\le w_2,w_3\le2g\).
The first factors at rotations two and three both have middle weight
\(3g\), hence both attain the same global minimum. Their terminal
sizes are respectively
\[
                 (g+w_3,w_2),\qquad(g+w_2,w_3).                   \tag{34}
\]
All four sizes exceed one. Corollary 11 must apply to **both new
minimizing factors**. It forces \(g+w_3=w_2\) and
\(g+w_2=w_3\), an impossibility. This excludes ratio three.
Corollary 11 leaves exactly \(\ell\mid s\) for equal terminals.

For singleton terminals, no terminal-to-terminal factor is possible.
Thus every edge of the size-\(t\) terminal tree is paired with a
middle edge. In either of the old cases \(\ell=3t/2\) and
\(\ell=4t/3\), apply Lemma 12 to **each** of those edges, using
the one fixed cyclic coordinate system. Its terminal endpoints have
the same residue modulo \(g=\gcd(t,\ell)\). The terminal tree is
connected, so all its vertices have one residue modulo \(g\).
Its cyclic coordinates exhaust \(\mathbb Z/t\), forcing \(g=1\).
But Corollary 11 requires respectively \(t=2g\) with \(4\mid t\),
or \(t=3g\) with \(6\mid t\). Both force even \(g\), a
contradiction. The other two cases in (27) give (31). \(\square\)

The ratio-three argument uses only complementary twists, pair rotations
and reapplication of a global restriction to tied minimizers. No condition
on the individual middle-middle factors, nor a three-factor completion
of an isolated seven-weight witness, is needed for this exclusion.

## Lemma 14. The tree support code has an exact invariance test

For a polynomial DT tree let \(C\subseteq\mathbb F_2^N\) be the span
of the factor support indicators. Then \(\dim C=r=(N-1)/2\), and
\(C\) is invariant under the generated permutation group if and only
if every pair of factor supports has even intersection.

**Proof.** The vertex incidence map from the binary edge space of a
tree to the even-weight vertex space is an isomorphism: its kernel
would be a nonempty subforest of even degrees, contradicted by a leaf,
and both spaces have dimension \(N-1\). Each factor support is the
image of the sum of its two edges. These disjoint edge pairs give
\(r\) independent vectors. A single edge incidence cannot lie in
\(C\): its unique preimage has unequal coefficients in that pair.

For a DT \(D\) acting on another factor's support indicator \(v\),
\(D(v)+v\) is the sum of those edge incidences of \(D\) whose
edges have exactly one endpoint in that support. An odd intersection
gives exactly one such edge, whose incidence is outside \(C\).
Hence invariance implies even intersections. Conversely, an even
intersection gives zero or both edge incidences, so \(D(v)+v\)
is zero or the support of \(D\). This proves invariance. \(\square\)

Consequently this proper support code is excluded by full alternating
monodromy, as in Lemma 4 of the
[paired-tree note](PAIRED_HURWITZ_LEAF_EXTRACTION.md#lemma-4-the-even-overlap-code-is-a-hurwitz-invariant).
Exact ratio-three covers actually produce an odd three-label contact;
they do not by themselves establish code invariance. Corollary 13
resolves that case through new minimizers instead.

The support-code strategy for the singleton family in Theorem 17 would force
even intersections for **all** factor pairs or an endpoint decrease.
The one-cut covers (24) constrain contacts with the minimizer, but do
not constrain every middle-middle attachment or the overlaps between
two complementary factors. They alone do not imply an invariant code:
take two seven-vertex paths \(M=(0,\ldots,6)\), \(U=(7,\ldots,13)\),
the singleton \(S=\{14\}\), and
\[
 A=(14\,0)(1\,7),\qquad
 B_i=(i\,i+1)(7+i\,8+i)\quad(0\le i<6).
\]
The paired paths give the same rooted cyclic coordinates, so every
terminal hit is accompanied by a middle hit in every complementary
twist. Nevertheless \(A,B_0\) overlap in three labels. The last
factor is removable, so this is not a positive global minimum and
does not refute the proposed global bridge. Lemma 19 instead uses stability
of the antipodal projection under complementary braids to force an endpoint
decrease in the remaining saturated family. This weaker phase-compatible
example lacks that stability. Mixed connectivity remains separate.

## Lemma 15. Equal-middle singleton cuts force a removable factor

At a positive global orbit minimum, no minimizing factor has terminal
sizes \((1,t)\) and middle size \(t\) for any \(t\ge3\). This does
not require full alternating monodromy. The exception \(t=2\) is real.

**Proof.** Place the minimizing factor first by Lemma 9 and write
\[
 A=(a\,b)(c\,d),\qquad S=\{a\},\quad |M|=|U|=t,
 \quad b,c\in M,\ d\in U.
\]
Its mass is \(m=t-1\), and \(N=2t+1\). By Theorem 10 and the edge
counts (28), every complementary factor pairs one edge \(e\) of the
tree on \(M\) with one edge \(f\) of the tree on \(U\); there are
no middle-middle factors. For such a pair, let \(M_c(e)\) be the
component containing \(c\) after deleting \(e\) from the middle
tree, and let \(U_d(f)\) be the component containing \(d\) after
deleting \(f\) from the terminal tree. The middle component after
deleting this factor from the full tree is precisely these two sets,
joined by \(cd\), together with \(a\) if \(b\in M_c(e)\).
Its factor mass is at least \(m\), so
\[
 |M_c(e)|+|U_d(f)|+\boldsymbol1_{b\in M_c(e)}\le t.       \tag{35}
\]

If a leaf \(v\) of the middle tree is different from \(b,c\), take
its incident edge for \(e\). The left side of (35) is then at least
\((t-1)+1+1=t+1\), a contradiction. Thus the middle tree has no
leaves except \(b,c\), and is a path with these endpoints. For any
leaf \(w\ne d\) of the terminal tree, its incident edge \(f\) has
\(|U_d(f)|=t-1\). Inequality (35) forces
\(M_c(e)=\{c\}\). The only possible paired edge is therefore the
unique middle edge incident to \(c\). Distinct terminal leaf edges
have distinct partners, so there is at most one leaf different from
\(d\). Every nontrivial tree has at least two leaves; hence the
terminal tree is a path with endpoint \(d\).

Write the paths in their respective orders as
\[
 M=(v_0=b,v_1,\ldots,v_{t-1}=c),\qquad
 U=(w_0=d,w_1,\ldots,w_{t-1}).
\]
If the middle edge \(v_iv_{i+1}\) is paired with terminal edge
\(w_jw_{j+1}\), (35) becomes
\((t-i-1)+(j+1)\le t\), or \(j\le i\). The pairing is a
permutation of \(0,\ldots,t-2\). Summing these inequalities forces
equality at every index. Thus the complementary factors, in some
order, are exactly
\[
 B_i=(v_i\,v_{i+1})(w_i\,w_{i+1}),\qquad 0\le i\le t-2.   \tag{36}
\]

Since \(t\ge3\), there is a pair \(B_i,B_{i+1}\). Among all such
pairs choose one whose positions in the complementary block are
closest. Every intervening factor commutes with both: a factor that
failed to commute would have a neighboring index and give a closer
pair. Hurwitz moves across these commuting factors are literal swaps,
so bring the chosen pair adjacent without changing any factor or any
tree degree. Their two common labels \(v_{i+1},w_{i+1}\) are internal
vertices of their paths, are fixed by \(A\), and occur in no other
factor. If the adjacent order is \((C,D)\), one positive move gives
\[
                    (C,D)\longmapsto(D,DCD).                \tag{37}
\]
Conjugation removes both common labels from \(DCD\). They become
exclusive leaves on the two different edges of the unchanged \(D\),
as in Lemma 4. Hence the new tuple has mass zero, contradicting the
positive global minimum.

For \(t=2\), the pair \(((a\,v_0)(v_1\,w_0),
(v_0\,v_1)(w_0\,w_1))\) is a counterexample to extending the lemma.
Its Hurwitz orbit has five states, all of mass one; its generated
group is dihedral of order ten. There is only one complementary
factor, so the last step is unavailable. \(\square\)

The proof uses global minimality only for Theorem 10 and the ordinary
factor inequalities (35). Once the pairing (36) holds, arbitrary
complementary factor order is handled by commuting swaps; no general
claim about bringing a two-contact pair is used.

## Corollary 16. The remaining positive-minimum cuts

At a positive global orbit minimum, every minimizing factor with
\(N\ge15\) satisfies one of the following necessary conditions:
\[
 \begin{array}{c|c}
 \text{terminal sizes}&\text{middle size}\\ \hline
 (s,s),\ s>1&\ell\text{ odd and }\ell\mid s\\
 (1,t),\ t>1&\ell=2t\text{ with }2\mid t.
 \end{array}                                                   \tag{38}
\]

**Proof.** Corollary 13 leaves these cases and \((1,t;t)\). In that
last case \(N=2t+1\ge15\) implies \(t\ge7\), so Lemma 15 excludes
it. \(\square\)

## Theorem 17. Only the equality singleton family can remain

At a minimizing factor \(A\) of **any** polynomial DT tuple, no other
factor pairs an edge of \(S_A\) with an edge of \(U_A\). This assertion
requires neither a global orbit minimum nor any monodromy hypothesis.

Consequently, at a positive global orbit minimum with \(N\ge15\), every
minimizing factor has
\[
 (|S_A|,|U_A|;|M_A|)=(1,t;2t),\qquad 2\mid t,\qquad
 N=3t+1,\qquad m=t-1=\frac{N-4}{3}.                 \tag{39}
\]
In particular every polynomial DT orbit on \(N\ge15\) sheets with
\(N\not\equiv1\pmod6\) contains a removable factor. This is extraction
in those degrees, not an assertion of mixed or all-degree connectivity.

**Proof.** Suppose a factor \(B\) pairs an edge in \(S_A\) with an
edge in \(U_A\). Let \(S_0\) and \(U_0\) be the sides of these
edges containing the attachment roots \(a\) and \(d\) of \(A\).
Deleting \(B\) leaves middle component
\[
 M_B=M_A\cup S_0\cup U_0.
\]
Its terminal components are \(S_A\setminus S_0\) and
\(U_A\setminus U_0\). Both root sides are nonempty. Hence
\[
 \lambda(B)=\lambda(A)-|S_0|-|U_0|
             \le\lambda(A)-2,                       \tag{40}
\]
contradicting that \(A\) minimizes the tuple's factor masses.
Thus \(x=0\) in (28).

For the equal-terminal family of Corollary 16 this gives
\(y=z=s-1\), so \(\ell-1=2s-2+2q\ge2s-2\).
But \(\ell\mid s\) gives \(\ell\le s<2s-1\) for \(s>1\),
a contradiction. Corollary 16 leaves exactly (39). The degree congruence
and the value of \(m\) follow immediately. If the degree congruence
fails, Proposition 8 forces the attained orbit minimum to be zero,
and Lemma 1 supplies a removable factor. \(\square\)

The ordinary all-factor comparison (40) is stronger here than attempting
to construct an invariant from the synchronized terminal forests. The
earlier two-family statement remains a valid intermediate reduction.

## Lemma 18. Saturated singleton covers force antipodal roots

In the remaining family (39), front-place the minimizer as in (20).
There are exactly \(t-1\) terminal-middle factors and \(t/2\)
middle-middle factors. For even \(t\ge4\), use the fixed cyclic coordinates
\[
 R_M^k b\leftrightarrow k\pmod{2t},\qquad
 R_U^k d\leftrightarrow k\pmod t.
\]
Then \(c=R_M^t b\). The projection
\[
 \pi:M\longrightarrow U,\qquad
 \pi(R_M^k b)=R_U^k d                              \tag{41}
\]
has precisely the antipodal pairs as fibres. Every terminal-middle
factor pairs a middle edge \(pq\), with \(\pi(p)\ne\pi(q)\),
with the terminal edge \(\pi(p)\pi(q)\). These \(t-1\) projected
edges are exactly the terminal tree, each used once.
At every complementary twist phase such a factor has support
intersection zero or two with the minimizer. An odd intersection with
the minimizer can therefore only be a single middle contact from a
middle-middle factor.

**Proof.** The counts follow from (28), with \(s=1,u=t,\ell=2t\).
Put \(c=h\ne0\) and let \(p,q\) be the middle endpoints of a
terminal-middle factor. The two distinct terminal endpoint residues
develop to four distinct phases modulo \(2t\). Their union is contained
in
\[
 E=\{p,q,p-h,q-h\}.
\]
Thus \(E\) has four elements and equals that union. In particular it is
invariant under addition of \(t\). The two sets \(\{p,q\}\) and
\(\{p-h,q-h\}\) are disjoint. At a covered phase exactly one middle
root and one terminal root meet this factor; at every other phase
neither does. This proves the even-intersection assertion already,
without assuming \(h=t\).

If \(h\ne t\), translate to \(p=0\) and write \(q=e\).
The antipode of zero in \(E=\{0,e,-h,e-h\}\) is \(t\).
It cannot be \(-h\), so either \(e=t\) or \(e-h=t\).
In the second case \(e=h+t\), and the antipode \(e+t=h\)
must also belong to \(E\). It cannot be \(0,e,e-h=t\);
equality with \(-h\) forces \(h=t\), contrary to assumption.
Consequently \(e=t\): every terminal-middle edge would be a diameter
of this same \(2t\)-cycle.

There cannot be two such edges in the complementary middle tree.
Here is the needed elementary noncrossing fact, with its factor-order
justification. Expand the complementary middle permutations into their
individual transpositions. They form a reduced factorization of \(R_M\).
Two selected disjoint transpositions can be brought unchanged to the
first two positions by ordinary transposition Hurwitz moves. If their
four endpoints alternated on \(R_M\), left multiplication of \(R_M\)
by their product would still be a full cycle: the first transposition
splits the cycle and the second rejoins the two pieces. That full cycle
cannot be a product of the remaining \(2t-3\) transpositions, since a
\(2t\)-cycle requires at least \(2t-1\). Thus their endpoints do not
alternate. Distinct diameters have four alternating endpoints and are
excluded. There are \(t-1\ge2\) terminal-middle edges, so \(h\ne t\)
is impossible.

With \(h=t\), equality of the phase sets says exactly that the two
terminal endpoint coordinates are \(p,q\pmod t\), and they are
distinct. Every terminal edge has exactly one middle partner, proving
(41) and the tree assertion. A middle-middle factor can meet the support
of \(A\) only at \(b,c\), so an odd intersection is exactly one label.
\(\square\)

The antipodal pairing in (41) is **not** an invariant pairing for the
generators: a terminal-middle factor swaps only one endpoint from each
of two distinct antipodal pairs in \(M\), while their mates are fixed.
It therefore breaks that pairing. More generally no nonempty collection
of disjoint sheet pairs can be invariant under this transitive odd-degree
group. The paired vertices would be an invariant nonempty subset and
hence all sheets, which is impossible in odd degree. Allowing one
unpaired fixed sheet does not avoid this transitivity obstruction.

Lemma 18 by itself controls the minimizer against terminal-middle factors.
Its additional force is that the *same* projection must work after every
complementary braid. The next lemma exploits that requirement without
asserting an invariant matching or code. The \(t=2\), degree-seven order-168
controls show why the lower size restriction matters.

## Lemma 19. Complementary braids cannot preserve a saturated projection

Let a polynomial DT tuple have first factor
\(A=(a\,b)(c\,d)\), with cut components
\(\{a\},M,U\), where \(|M|=2t\), \(|U|=t\), and \(t\ge4\) is even.
Suppose its complementary factors consist of \(t-1\) terminal-middle
factors and \(t/2\) middle-middle factors. Fix a two-to-one map
\(\pi:M\to U\). It is impossible that, after **every** braid on the
complementary factors, each terminal-middle factor has the form
\[
                (x\,y)(\pi(x)\,\pi(y)),\qquad \pi(x)\ne\pi(y).
                                                               \tag{42}
\]
Thus some finite complementary braid violates (42).

**Proof.** Complementary moves keep all edges within their original
components and retain \(A\). The graph remains a tree, so its edges on
\(U\) form a tree \(T_U\), each supplied by exactly one terminal-middle
factor. The middle edges of those factors will be called the lifted edges.
Any two selected complementary factors can be brought adjacent while keeping
both literal permutations unchanged: move the right selected factor left
through intervening factors by positive moves. Only crossed factors are
conjugated. Either adjacent order permits conjugating either selected factor
by the other while retaining that other factor. All these operations leave
\(A\) and the complementary product fixed.

Assume (42) holds throughout this complementary braid orbit. If two edges
\(uv,uw\) of \(T_U\) use distinct middle lifts \(x,x'\) of \(u\),
their factors have the form
\[
 C=(x\,y)(u\,v),\qquad D=(x'\,z)(u\,w),
 \quad \pi(y)=v,\quad\pi(z)=w,\quad x\ne x'.
\]
Their middle supports are disjoint, since \(u,v,w\) are distinct.
Gather \(C,D\) unchanged and conjugate \(D\) by \(C\). The middle
edge \(x'z\) is unchanged, but the terminal edge becomes \(vw\).
This violates (42). Therefore all lifted edges incident at any terminal
vertex use the same lift. Since \(T_U\) is connected and has no isolated
vertex, the lifted edges form an isomorphic tree \(F\) on
\[
                  L=\{x_u:u\in U\}\subset M,
                  \qquad \pi(x_u)=u.                         \tag{43}
\]
This conclusion holds in every complementary state.

Let \(B\) be any middle-middle factor in one such state. For each lifted
edge \(x_ux_v\), gather its factor \(C\) and \(B\), then conjugate
\(C\) by \(B\). Its terminal edge is unchanged, so (42) requires
\[
 \{\alpha(u),\alpha(v)\}=\{u,v\}\quad(uv\in T_U),
 \qquad \alpha(u)=\pi(Bx_u).                                  \tag{44}
\]
The selected factors, \(\pi\), and the displayed labels are unchanged by
gathering, so (44) constrains the factors in the original state.
Every vertex with two distinct neighbors is fixed by \(\alpha\), by
intersecting the two endpoint sets in (44). Every leaf has such a neighbor
when a tree has at least three vertices; its own equation then fixes that
leaf. Consequently
\[
                         \pi(Bx)=\pi(x)\quad(x\in L).          \tag{45}
\]
We have proved (45) for every middle-middle factor in every state.

Suppose \(B\) moves some \(x\in L\), and choose a lifted edge
\(xy\) with factor \(C=(x\,y)(u\,v)\). Gather \(B,C\)
unchanged. The resulting lifted tree may differ elsewhere, but it still
contains \(xy\). By (45), \(Bx=x^*\), where \(x^*\ne x\)
is the other sheet above \(u\); it lies outside the current lifted
vertex set. Now conjugate \(B\) by \(C\), retaining \(C\).
Every terminal-middle factor remains unchanged by this single move,
so the lifted tree and its vertex set remain fixed. The new middle-middle
factor \(B'=CBC\) satisfies
\[
                   B'y=CBx=Cx^*=x^*,
                   \qquad \pi(B'y)=u\ne v=\pi(y),             \tag{46}
\]
contradicting (45) in the new state. Thus every middle-middle factor
fixes \(L\) pointwise. Every terminal-middle factor preserves \(L\).
No complementary edge can therefore join the nonempty sets \(L\) and
\(M\setminus L\), whereas the complementary edges on \(M\) form a
connected tree. This is the contradiction. \(\square\)

The edge-set argument (44) is where \(t\ge3\) is essential. On a
two-vertex terminal tree a map can exchange its endpoints and still
preserve its sole edge. No invariant antipodal matching is used in this
proof; the fibres of \(\pi\) serve only as fixed coordinates against
which to test the changed factors.

## Theorem 20. Every residual pure-DT orbit admits extraction

Every polynomial DT braid orbit on odd \(N\ge15\) sheets contains a
factor with one exclusive leaf on each of its two transpositions. No
monodromy hypothesis is needed for this extraction statement. In particular
the remaining equality singleton minimum is impossible, and the existing
degree-thirteen base and two-leaf absorption theorem imply pure-DT
connectivity in every residual degree, both unmarked and with the first
branch marked.

**Proof.** Suppose the orbit minimum is positive. By Theorem 17 every
minimizing representative has a minimizing factor with cuts
\((1,t;2t)\), even \(t\), and \(N=3t+1\). Thus \(t\ge6\)
in the stated degree range. Front-place this factor by Lemma 9. Every
braid on the complementary factors retains its three component sets,
its mass \(t-1\), and their product \(R\). Global minimality therefore
makes every such endpoint a minimizing representative. The cycles
\(R_M,R_U\) and the attachment roots \(b,d\) remain literal, so Lemma 18
gives (42) in every state for the **same** map
\(\pi(R_M^k b)=R_U^k d\). Its factor counts also remain unchanged.
Lemma 19 is a contradiction. Proposition 8 and Lemma 1 now give extraction.

For connectivity use the retained complete degree-thirteen pure-DT base
and Corollary 3 of the
[two-leaf absorption theorem](UNIFORM_TWO_LEAF_DT_ABSORPTION.md#corollary-3-recursive-extensions-and-the-buffered-families).
Move a removable factor to the end by positive moves. Every crossed factor
fixes its two leaves, so only its old anchors are conjugated. Delete those two
leaves and that factor. The remaining tree has full-cycle product and
degree \(N-2\); induction normalizes it to the canonical component.
Two-leaf absorption restores the deleted factor at its two distinct old
anchors. This proves absolute unmarked connectivity for every odd
\(N\ge15\). The full twist realizes every product-preserving sheet
relabeling, and Theorem 4 of the
[marking-removal note](UNMARKED_POLYNOMIAL_HURWITZ_REDUCTION.md#theorem-4-unmarked-connectivity-suffices-for-both-marked-spaces)
then gives connectivity under its exact marked subgroup. \(\square\)

There is also an explicit escape interpretation of the contradiction.
If a complementary braid first breaks (42), choose a terminal endpoint
residue absent from the two projected middle endpoints of that factor.
A complementary full twist placing this residue at the fixed root \(d\)
gives precisely one outer contact with \(A\): neither antipodal middle
root is hit. Lemma 3 lowers the mass. An earlier decrease already suffices.
The proof uses pair conjugations, cut-preserving gathering, and a
complementary twist; it need not establish invariance of the binary support
code. Mixed connectivity, independent mathematical review and formal
verification remain separate.

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

The [global compatibility checker](verify_terminal_mass_global_compatibility.py)
and its [receipt](terminal_mass_global_compatibility_verification.json)
replay front placement and complementary full twists in the retained
392-tuple degree-seven domain. Its two proper-group orbits have positive
minimum one, and all their minimizing marked inputs satisfy the cyclic
covers and terminal-size restrictions. Both full-\(A_7\) orbits have
minimum zero. The degree-fifteen plateau example (9) develops an outer
one-contact after one complementary full twist; either that endpoint
already lowers \(\mu\), or Lemma 3 then does. Four synchronized-tree
controls have removable factors, while an unequal edge matching between
equal-sized terminal trees and a DT wholly inside a terminal component
are detected as incompatible with the phase covers. These are internal
controls of the written global argument, not a higher-degree census or
a computational proof of complete extraction.

The same checker also exhausts the modular three- and four-coset
implications for \(\ell=3g,4g\), \(1\le g\le12\), and checks
the congruence lifts for the tested terminal sizes. It independently
recomputes all ten middle and terminal masks in (32), rejects a corrupted
mask, and checks the tied-minimizer contradiction for odd
\(3\le g\le31\). Binary support-span enumeration checks the exact
code criterion on the retained degree-seven domain and the displayed
degree-fifteen phase-compatible counterexample to the weaker implication.
These bounded controls supplement the universal written proofs of
Lemmas 12 and 14 and Corollary 13; they are not the reason those
statements hold for unbounded sizes.

The [singleton-cut checker](verify_terminal_mass_singleton_equal_cut.py)
and its [receipt](terminal_mass_singleton_equal_cut_verification.json)
provide dedicated exact controls for Lemma 15. They enumerate all
labeled middle and terminal trees and edge pairings for \(2\le t\le4\),
test the forced path pairing and arbitrary complementary order, and
replay the final braid with direct exclusive-leaf and lower-cycle
checks. The positive degree-five dihedral orbit checks the necessary
exception. The receipt gives the exact larger path-only ranges and
negative controls. These finite controls supplement the universal
proof; they do not prove it by extrapolation or settle the remaining
singleton family (39).

The global compatibility checker also tests (40) by direct tree cuts
in all 126 terminal-bridge instances in the retained degree-seven
domain. It checks 966 disjoint-edge noncrossing comparisons in the
complementary component cycles, and exhausts translation-normalized
saturated two-coset covers for even \(2\le t\le40\). The controls
retain valid single-diameter covers with nonantipodal roots; the written
noncrossing argument is what excludes multiple such edges in a tree.
They also verify that a terminal-middle generator breaks the natural
antipodal pairing. These bounded tests support Theorem 17 and Lemma 18;
the universal statements follow from their written arguments.

The [projection checker](verify_terminal_mass_projection_escape.py) and its
[receipt](terminal_mass_projection_escape_verification.json) check all 1,076
lift assignments over labeled terminal trees of sizes two through four,
detecting split lifts by the explicit pair conjugation. They also test all
82,245 labeled-tree/middle-middle pairs for terminal sizes three through five,
including every allowed active attachment and its conjugation obstruction.
A two-vertex edge exchange checks the necessary limitation of (44). Five
deterministic saturated full-tree inputs at \(t=4,6,8,12,16\) have exact
projection-breaking braid words, complementary twist phases, outer contacts
and strict pruning certificates. These are internal bounded controls; the
unrestricted conclusions follow from the written proof, not extrapolation.

The [receipt](terminal_mass_hurwitz_descent_verification.json) and
[source ledger](terminal_mass_hurwitz_descent_source_ledger.json) record
the replays and exact inputs. Extraction for every odd \(N\ge15\) and the
resulting pure-DT connectivity are written mathematical conclusions. No new
degree census, universal monotone plateau escape, mixed connectivity, formal
verification, independent replay, external review or priority claim is asserted.
