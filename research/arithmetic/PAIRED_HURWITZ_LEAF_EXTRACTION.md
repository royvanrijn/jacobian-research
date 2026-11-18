# Paired Hurwitz trees: an overlap gate and leaf extraction

Every full-alternating polynomial DT tuple can be braided to two
adjacent factors sharing exactly one support label. A binary-code
obstruction and a five-cycle rotation prove this for arbitrary input.
At that pair, three Hurwitz moves exchange the secondary edges.

A separate four-move identity extracts two leaves under explicit
external-support conditions. A block-evacuation braid repairs a
broader class of external interactions. The cyclic tuples from the
[unmarked reduction](UNMARKED_POLYNOMIAL_HURWITZ_REDUCTION.md)
then satisfy an exact lower-degree recursion.

Turning the arbitrary-input overlap gate into arbitrary-input
leaf extraction and full connectivity remains open.

## Lemma 1. Pure four-move extraction

Suppose a polynomial branch tuple has full-cycle product and contains
three adjacent DT entries

\[
 A=(02)(34),\qquad B=(24)(56),\qquad C=(46)(78),
\tag{1}
\]

where the eight displayed labels are distinct. Assume every other
entry fixes both 2 and 4. The entries outside this window may be
DTs or a single 3-cycle as in the two polynomial passports.

Apply, from left to right, the four moves

\[
                       W=(-1,-2,-2,-1)
\tag{2}
\]

relative to this window. The result in that window is

\[
              E=(23)(45),\qquad
              A'=(03)(56),\qquad B'=(36)(78).
\tag{3}
\]

The first new DT has exclusive leaf ends 2 and 4. Delete that
entry and those two symbols. The remaining tuple has a full cycle
on the remaining symbols as its product, with one fewer DT.

**Proof.** Use Q_i^-1:(a,b)->(aba^-1,a). The calculation is

\[
 (A,B,C)\longmapsto(A',A,C)
       \longmapsto(A',B',A)
       \longmapsto(A',(02)(46),B')
       \longmapsto(E,A',B').
\]

This proves the exact identity, leaves every entry outside the
window unchanged, and preserves the original tuple product.
The strand permutation of (2) is

\[
                         (12)(23)^2(12)=1.
\]

Thus W is a pure braid.

Afterwards only E moves either 2 or 4. In a transposition-edge
expansion of the polynomial tuple, these are therefore leaf
vertices on its two different edges. For a 3-cycle entry use either
of its two-edge transposition expansions.

For completeness such an expanded graph is a tree: it is connected
because the full-cycle product lies in the generated group, and its
number of edge occurrences is exactly the number of vertices minus
one. A connected graph needs at least that many distinct edges,
so these edges are all distinct and form a tree. Removing the two
leaves and their incident edges leaves a tree. A product using
every transposition edge of a tree exactly once is a full cycle,
regardless of the edge order, by the leaf-insertion proof in
Corollary4 of the
[cycle-folding note](ALTERNATING_HURWITZ_CYCLE_FOLDING.md).
The remaining tuple therefore has the claimed product. ∎

For a window beginning at factor i, the word is
(-i,-(i+1),-(i+1),-i). If the marked branch lies outside this
window, its position and label are retained. If the extracted
DT itself was marked, deleting it requires choosing a new marked
branch for the smaller problem. The present unmarked application
needs no such choice.

## Corollary 2. Exact cyclic recursion down to degree nine

Let G_r be the degree-2r+1 cyclic DT tuple from Theorem2 of the
[unmarked reduction](UNMARKED_POLYNOMIAL_HURWITZ_REDUCTION.md),
with r>=4. Its first three entries are exactly (1).
No other entry moves 2 or 4. Applying (2), deleting E and its
two leaf labels, and relabeling

\[
                 0\mapsto0,\quad1\mapsto1,\quad3\mapsto2,
                 \qquad j\mapsto j-2\quad(j\ge5)
\tag{4}
\]

gives precisely G_(r-1).

**Proof.** The next unchanged entry, when present, is (68)(9,10),
and the later entries continue the same pattern. The final entry
is (2r-2,2r)(01); these all fix 2 and 4.
Under (4), A' and B' become

\[
                           (02)(34),\qquad(24)(56).
\]

Every later entry becomes the next cyclic entry, including the
wraparound (01) in the final one. This proves equality of the
entire lower tuple, not only its passport or generated group. ∎

For r>=5 the lower group is still A_(2r-1), by the explicit
3-cycle and primitivity proof for G_(r-1). Repeated extraction
therefore reaches the full A9 core in exactly r-4 deletion
steps and four local Hurwitz moves per step.

This is a recursion on successive smaller tuple problems.
Reinserting earlier deleted factors can attach them to vertices
removed at a later stage; the statement does not assert that
all extracted factors simultaneously have exclusive leaves
in the original full tuple.

## Lemma 3. Evacuating a block of external pivot interactions

The strict pivot-fixing condition of Lemma1 can sometimes be repaired
by a further explicit braid. After its four-move identity, write the
tuple as

\[
              (E,U_1,U_2,\eta_1,\ldots,\eta_m,\text{rest}),
\]

where E=(23)(45), U_1=(03)(56), U_2=(36)(78), and put
T=eta_1...eta_m. Let F be the common fixed set of U_1,U_2.
Suppose

\[
 \eta_j(3)=3,\quad\eta_j(5)=5\quad\text{for every }j,
 \qquad T(2),T(4)\in F,
\tag{6}
\]

and every remaining entry fixes 2 and 4. The eta_j may move the
pivots. They may be DTs or the single 3-cycle permitted by the
polynomial passport.

There is a further braid of length 3m that makes 2 and 4
exclusive leaf ends of E. Thus the lower-degree deletion applies.

**Proof.** First bring the eta block unchanged to the left across
U_1,U_2 using positive Hurwitz moves. It becomes

\[
          (E,\eta_1,\ldots,\eta_m,T^{-1}U_1T,T^{-1}U_2T,
                    \text{rest}).
\]

Then move E to the right across that block using negative moves.
The result is exactly

\[
 (E\eta_1E,\ldots,E\eta_mE,\ E,\
                 T^{-1}U_1T,T^{-1}U_2T,\text{rest}).
\tag{7}
\]

The first stage uses, for j=1,...,m, the two positive moves
Q_(j+2),Q_(j+1). The second uses Q_1^-1,...,Q_m^-1.
This proves the formulas and the length assertion.

Because E exchanges 2 with 3, its conjugate E eta_j E fixes 2
exactly when eta_j fixes 3; the same assertion holds for 4 and 5.
Furthermore T^-1 U_i T fixes a pivot p exactly when T(p)
belongs to Fix(U_i). Conditions (6) therefore make every factor
other than E fix both pivots. The tree deletion argument of
Lemma1 applies. ∎

A convenient sufficient version of (6) requires every eta_j
to preserve F setwise and to fix 3,5. The actual product condition
in (6) is weaker; the individual factors need not preserve F.

The repair allows arbitrarily many external interactions. For example,
begin with the cyclic G_5, insert after its first three factors the
DT block

\[
                  \eta_j=(2\,x_j)(4\,y_j),\qquad1\le j\le m,
\]

with all x_j,y_j fresh and distinct, and then retain its last two
factors. The initial pivots have arbitrarily large degree.
All eta_j fix 3,5 and send the pivots within F, so (7) evacuates
every interaction and extracts E. This is a uniform family of
interactions satisfying the repair hypotheses, rather than a claim
that all external interactions do so.

## Exact limitations of stronger extraction claims

The condition that the other entries fix 2 and 4 is essential.
Take G_5 on eleven symbols, extend every entry to fix 11,12, and
append

\[
                             (2\,11)(0\,12).
\tag{5}
\]

This still has a full-cycle product and full group A13. The tree
has acquired two leaves, and its odd-degree small-support argument
makes the generated group primitive; it contains the old A11,
so Jordan's 3-cycle theorem applies.

The four-move identity remains valid, but after it the degrees at
2,3,4,5 in the full edge graph are respectively 2,3,1,2.
Consequently the edge (23) is not a leaf edge. The output E has
only one exclusive leaf end, and the asserted two-leaf deletion
does not apply after those four moves alone. This example satisfies
Lemma3 with m=3: take the block consisting of G_5's final two
factors followed by (5). All three fix anchors 3,5, and their
product sends 2,4 to 11,4 in F. The additional nine moves repair it. Thus the original local pattern
needs either its external-support condition or an explicit cleanup,
even for full alternating monodromy.

The anchor-fixing condition in the repair has its own genuine blocker.
Instead insert eta=(2 11)(3 12) after the first three factors of
G_5. Its tree extension still has full group A13. It satisfies
T(2),T(4) in F, but moves anchor 3. Formula (7) then contains

\[
                         E\eta E=(3\,11)(2\,12).
\]

Consequently another factor still moves pivot 2, and the edge (23)
has no leaf end. This disproves an unconditional version of the
block-evacuation claim. It does not prevent a different braid word
from extracting some removable factor.

There is also a low-degree full-group issue. For G_4 the extraction
returns G_3 in degree seven. Its group has order 168, rather than
|A7|=2520. The checker certifies this by complete finite group
closure: it constructs 168 distinct elements and verifies closure
under all generators. Thus successful leaf extraction need not
preserve full alternating monodromy in this boundary degree.
Stopping the uniform full-group recursion at degree nine avoids
that issue.

## Lemma 4. The even-overlap code is a Hurwitz invariant

For a tuple of DTs on N labels let s_i in F2^N be the support
indicator of its i-th factor. If all support intersections have even
size, their span C is self-orthogonal for the ordinary binary dot
product. Every Hurwitz move preserves the literal subspace C and
the even-overlap condition.

Moreover the generated permutation group preserves C. Consequently
a tuple generating A_N, with N>=5, must contain two supports with
odd intersection.

**Proof.** Let d be one DT and s another support indicator with
even intersection with supp(d). Each of d's two transposition edges
has zero, one or two endpoints in s. An even number of edges,
either zero or two, have exactly one endpoint in s. Therefore

\[
                            d(s)=s
                     \quad\hbox{or}\quad d(s)=s+s_d.
\tag{8}
\]

A Hurwitz move replaces a support pair by (s_j,s_i) or
(s_j,s_i+s_j). These are invertible row operations on the
generating vectors, so the span is unchanged. Its self-orthogonality
also makes all new intersections even. The same argument works
for every generator acting on every support vector, hence the
generated group preserves C.

If the group were A_N, its transitivity on four-element subsets
would put every weight-four vector into C. Two such vectors can
have intersection three when N>=5, contradicting self-orthogonality.
This proves the last assertion. ∎

The code is literal when the sheet labels are fixed. Simultaneous
relabeling carries both the product P and C to their corresponding
relabelings. Thus the pair (P,C), up to simultaneous relabeling,
is the corresponding absolute invariant on this even-overlap locus.

For the cyclic G_3 core, complete closure gives order168 and its code
has dimension three: its eight vectors have weights 0 once and 4
seven times. The full A7 seed

\[
                         ((01)(34),(12)(45),(23)(56))
\]

has an odd support intersection and lies outside this invariant locus.
These are exact controls of the scope. The code does not furnish
an invariant distinguishing full-alternating orbits.

## Theorem 5. Every full-alternating DT tree reaches overlap one

Let a polynomial tuple of DTs have full-cycle product on N>5 labels
and generate A_N. It is Hurwitz equivalent to a tuple with two
adjacent factors whose supports intersect in exactly one label.

This is an arbitrary-input statement. It uses neither the cyclic
normal form nor a higher-degree enumeration.

**Proof.** Lemma4 supplies an odd support intersection. Different
factors in the transposition tree share no edge, so an odd
intersection is one or three.

Any two chosen factors can be brought adjacent while leaving both
literal permutations unchanged: move the factor farther to the
right leftwards by positive Hurwitz moves, stopping just after the
other chosen factor. Each crossed factor is conjugated; the two
chosen ones are retained. The product and generated group remain
the same, and the expanded graph is still a tree.

An overlap-one pair already gives the conclusion. Otherwise bring
an overlap-three pair adjacent. Their union has five vertices and
four distinct edges, hence is a path; their product C is a 5-cycle.
Call its five labels the core.

Because N>5 and the whole graph is connected, some other DT gamma
touches the core. No edge of gamma can have both endpoints there:
it would close a cycle with the core path. Therefore gamma meets
the core in exactly one or two vertices.

If it meets the core in one vertex, at least one of the two DTs
contains that vertex, so its overlap with gamma is one. If the
intersection is {p,q}, use even powers Q^(2k) on the adjacent pair.
The move Q^2 simultaneously conjugates both factors by C^-1.
Each factor has a unique fixed vertex within the five-element
core, and conjugation makes the first factor's fixed vertex run
through all five vertices. Choose k in {0,1,2,3,4} so that its
fixed vertex is p. The first support then contains q but not p,
and its intersection with gamma is exactly one.

Finally bring that pair adjacent by the same unchanged-factor
operation. This proves the theorem. The only local rotation used
has at most eight moves, in addition to the two bringing
operations. ∎

## Lemma 6. An overlap-one pair exchanges its secondary edges

For seven distinct labels, an adjacent pair

\[
                   A=(uv)(pq),\qquad B=(vw)(rs)
\]

has the exact three-move replacement

\[
                  Q^3(A,B)=((uv)(rs),(vw)(pq)).
\tag{9}
\]

Thus the secondary edges are exchanged while the two primary edges
and the entire uncolored edge set are retained.

**Proof.** The successive pairs are

\[
 ((vw)(rs),(uw)(pq)),\qquad
 ((uw)(pq),(uv)(rs)),\qquad
 ((uv)(rs),(vw)(pq)).
\]

This is direct substitution in Q(a,b)=(b,b^-1ab). ∎

The move Q^3 swaps its two branch labels; it is allowed in the
unmarked problem. Bringing a selected pair adjacent can change
the surrounding colored tree. Formula (9) preserves the
uncolored tree at the stage where the pair is already adjacent.

## What remains to prove for uniform connectivity

Theorem5 reduces every full-alternating DT input to an adjacent
overlap-one pair, and Lemma6 gives a concrete edge exchange there.
A global leaf-pair extraction theorem must still use such moves
to produce a removable block, for example the complete pattern and
support conditions of Lemma1 or its repaired version in Lemma3.
No such leaf-extraction theorem is proved here. The recursion for G_r does not imply that every
tuple is in its orbit.

After a general deletion, one would also need control of all
ways of restoring the two leaves and their attachment vertices.
Identifying the smaller generated group alone does not supply
that extension connectivity. The two unmarked polynomial
connectivity questions, first unresolved in polynomial degree
thirteen, remain open in this work.

The [checker](verify_paired_hurwitz_leaf_extraction.py) verifies the
literal identity, its pure strand action, every recursive step in
retained cyclic controls through r=50, arbitrary allowed external
leaf attachments on the first relevant control, the counterexample
(5), the exact external-block conjugation formula and its paired
pivot-anchor blocker, the complete degree-seven group closure,
the binary-code update rule, all294 full-A7 degree-seven tuples in
the overlap-one gate, and the exact secondary-edge exchange.
The [receipt](paired_hurwitz_leaf_extraction_verification.json) and
[source ledger](paired_hurwitz_leaf_extraction_source_ledger.json)
record the exact evidence. Uniformity follows from the written
identity and relabeling proof. No higher-degree exhaustive census,
uniform connectivity theorem, formal verification or priority
claim is made.

