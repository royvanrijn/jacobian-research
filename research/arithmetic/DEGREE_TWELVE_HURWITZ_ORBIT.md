# The degree-twelve alternating Hurwitz orbit

There is one absolute full-A12 braid orbit for six double-transposition
finite branch values and a distinguished infinity value of profile (11,1).
The proof combines the uniform cycle-folding theorem with two complete
degree-eleven polynomial censuses and explicit braid words connecting
every expansion fibre. It does not enumerate the much larger original
space of six-tuples directly.

This is a finite computational theorem with a written reduction and short
word certificates. It extends the completed degree-eight and degree-ten
cases. Uniform connectivity in higher degrees remains open.

## Theorem 1. One absolute degree-twelve component

Let C be the double-transposition class in S12. All transitive tuples

\[
                  (x_1,\ldots,x_6)\in C^6,\qquad
                     x_1\cdots x_6\text{ of type }(11,1)
\]

belong to one braid orbit up to simultaneous sheet relabeling.
They all generate A12. Thus the absolute Hurwitz space with passport

\[
                   ((2,2,1^8)^6,(11,1))
\]

and distinguished infinity branch value has one geometric component.

**Proof.** We first give the algebraic reduction, then the exhaustive
finite steps and the literal words joining them.

### Reduction to marked polynomial tuples

The [uniform cycle-folding theorem](ALTERNATING_HURWITZ_CYCLE_FOLDING.md)
braids every tuple to a doubled edge at the fixed sheet, in two adjacent
blocks, with no other edge incident to that sheet. Rename the fixed sheet
star and the others 0,...,10. Move the adjacent pair to the first two
positions by cabling. It has form

\[
                         (\star\,a)p,\quad(\star\,a)q,
\]

where p,q are distinct transpositions and a lies outside their supports.
All other blocks fix star. Coalescing these two blocks gives
delta=pq, either a double transposition or a 3-cycle.

Delete star and relabel the long product to
P=(0 1 2 3 4 5 6 7 8 9 10). The result belongs to one of

\[
 \mathcal P_D=\{(\delta,g_2,g_3,g_4,g_5)\in D^5:
                                       \delta g_2g_3g_4g_5=P\},
\]
\[
 \mathcal P_T=\{(\delta,g_2,g_3,g_4,g_5)\in T\times D^4:
                                       \delta g_2g_3g_4g_5=P\},
\]

where D is the class of 990 double transpositions in S11 and T the
class of 330 3-cycles. The first entry is marked. These are genus-zero
polynomial branch tuples: every finite entry has index two and the
product is an eleven-cycle.

Conversely every tuple in either set expands by choosing a fixed sheet
star, an attachment a outside the support of delta, and an ordered
factorization delta=pq into two transpositions. Replace delta by
(star a)p,(star a)q. The expanded tuple is transitive: its group
contains the coalesced group, transitive on eleven symbols, and moves
star to a. Its product has type (11,1).

The 11-cycle forces primitivity in degree twelve, and a primitive group
of degree at least nine containing a double transposition contains the
alternating group. Since all entries are even, every such expansion
has group A12. The exact imported theorem is Jones, Proposition 6.1(3),
in the [retained source ledger](alternating_hurwitz_priority_ledger.json).

### Complete marked polynomial connectivity

Use Hurwitz moves

\[
 Q_i:(a,b)\longmapsto(b,b^{-1}ab)
\]

on entries i,i+1; inverses have (a,b)->(aba^{-1},a).
The subgroup of explicitly allowed operations

\[
                         Q_2,Q_3,Q_4,Q_1^2                   \tag{1}
\]

and their inverses keeps the marked entry first. It suffices to use
this subgroup; no assertion that it is the entire marked braid
stabilizer is needed.

The complete exact census proves:

| Domain | Number of tuples | Orbits under (1) |
|---|---:|---:|
| P_D | 14,992,384 | 1 |
| P_T | 3,748,096 | 1 |

For coverage, store every ordered pair of D-elements by its product.
Also store every ordered pair from D x D or T x D, according to the
first class. For a first-pair product m and a middle entry c in D,
the remaining ordered pair must have product c m^-1 P.
This enumerates every five-tuple exactly once. Breadth-first traversal
checks every move in (1), its inverse, and membership of every neighbor
in that exhaustive domain. It visits the complete respective set.

This marked transitivity is stronger than an unmarked polynomial
orbit calculation. In particular it carries the coalesced entry and
its expansion data throughout the path.

### Lifting the marked moves

For an expanded tuple with the adjacent marked pair first, an ordinary
Q_j with j>=2 lifts to the original six-strand Q_(j+1).
The marked operation Q_1^2 lifts to the literal word

\[
                              (2,1,1,2),                     \tag{2}
\]

where a positive integer i denotes Q_i applied in the written order.
Its inverse lifts to (-2,-1,-1,-2).
Indeed the first two moves carry the pair past its neighbor gamma;
the last two carry that neighbor past the pair. The coalesced tuple
undergoes Q_1^2, while the attachment transforms by gamma^-1=gamma.
For Q_1^-2 it transforms by delta gamma delta^-1.

These identities hold for every factorization delta=pq, in both
classes, so every path in either marked polynomial census lifts to
an actual B6 path between expansions.

It remains to connect different expansions over one base tuple.
For a fixed attachment, the internal B2 action connects all ordered
two-transposition factorizations of delta: there are two when delta
is a double transposition and three when it is a 3-cycle.

Attachment transitivity requires additional loops. It does not follow
merely from the coalesced group's name. The retained certificate gives
four explicit marked polynomial stabilizer loops for each census base.
Their induced permutations act transitively on all attachments, seven
in the D case and eight in the T case. The
[word verifier](verify_degree_twelve_hurwitz_orbit.py) independently
replays every loop and its literal cable (2), for every attachment.
Each returns the polynomial tuple exactly, fixes the remaining four
expanded blocks, and changes the common auxiliary edge from (star,a)
to (star,h(a)) as claimed.

Thus, within each of the two coalesced classes, all expansions lie
in one original braid orbit. Internal B2 moves remove any dependence
on which transposition splitting was chosen.

### A literal link between the two classes

On twelve sheets take

\[
 \begin{split}
 x_1&=(0\,1)(2\,3),&x_2&=(0\,1)(4\,5),\\
 x_3&=(3\,4)(6\,7),&x_4&=(1\,2)(7\,8),\\
 x_5&=(5\,6)(9\,10),&x_6&=(8\,9)(10\,11).
 \end{split}
\]

The product is

\[
                    (1\ 3\ 5\ 7\ 8\ 10\ 11\ 9\ 6\ 4\ 2)(0).
\]

The first two blocks share (0,1) and coalesce to a double transposition.
Apply the two-letter B6 word (-2,1). The first three entries become

\[
       ((3\,5)(6\,7),\ (0\,1)(2\,5),\ (0\,1)(4\,5)).
\]

The adjacent second and third entries now coalesce to a 3-cycle.
All other entries and the total product are preserved.
Both tuples are in the required transitive passport; alternatively
their connected doubled-edge graph and the displayed long product
verify this directly.

This links the two expansion orbits. The uniform normal form covers
every initial tuple, completing the proof. ∎

## The attachment certificates

The [retained JSON certificate](degree_twelve_hurwitz_certificate.json)
contains both five-permutation polynomial bases, the eight stabilizer
loops and their attachment permutations, and the twelve-sheet link.

For compactness its polynomial codes use +/-1 for Q_1^(+/-2), and
+/-j, 2<=j<=4, for Q_j^(+/-1). They are:

| Coalesced class | Four stabilizer-loop codes |
|---|---|
| Double transposition | (-2,-1,-1,-1,2); (2,-1,-1,-1,-2); (-3,-2,-1,2,3); (-4,-3,-2,-1,2,3,4) |
| 3-cycle | (-2,-1,2); (-3,-2,-1,2,3); (3,2,-1,-2,-3); (-4,-3,-2,-1,2,3,4) |

The receipt additionally records their expanded ordinary B6 words.
The codes are never treated as unexpanded B5 generator letters by
the verifier.

During the complete census, one transport permutation is stored along
each breadth-first tree path while attachment orbits are being computed.
A non-tree edge yields a literal stabilizer loop and its transport
permutation. Union-find on the seven or eight attachments verifies
transitivity. The separately retained words then verify this conclusion
without trusting those transport labels or the union-find computation.

The entire proof therefore separates two finite obligations: exhaustive
coverage of the marked polynomial domains, and independently replayable
short words for the expansion fibres.

## Arithmetic and remaining scope

The absolute passport is invariant under Galois, and its unique geometric
component is therefore defined over Q in the usual absolute Hurwitz
moduli problem. This does not supply a Q-rational point or a rational
degree-twelve coefficient model.

The elementary polynomial and tree reductions are uniform, but their
remaining connectivity step was exhaustively checked here only in
degree eleven. There is no all-degree transitivity theorem in this note.
The [degree-fourteen theorem](DEGREE_FOURTEEN_HURWITZ_ORBIT.md)
establishes one absolute full-alternating orbit in degree fourteen.
Degree sixteen and uniform connectivity remain open.

A local three-block normalization cannot be assumed without proof:
degree six has distinct full-A6 braid orbits. The proof above instead
uses the exact degree-eleven marked spaces and actual expansion loops.
No claim about the degree-six phenomenon is needed as a premise.

## Reproduction and assurance

Run the full proof replay from the repository root:

~~~sh
.venv/bin/python research/arithmetic/verify_degree_twelve_hurwitz_orbit.py
~~~

The [checker](verify_degree_twelve_hurwitz_orbit.py) embeds its C++17
enumeration and compiles it with assertions enabled. It needs g++
and about 1.6 GB RAM. It verifies class completeness, duplicate-free
enumeration, annular inverse and product identities, closure of every
breadth-first edge, both complete orbit sizes, attachment transport,
all retained words and the cross-class link. Wrong-product tuples
and corrupted word certificates are rejected.

The [receipt](degree_twelve_hurwitz_orbit_verification.json) pins the
checker, embedded source, certificate, compiler and
[uniform reduction source](ALTERNATING_HURWITZ_CYCLE_FOLDING.md).
Its coverage_replayed field is true for the retained full receipt.
An optional --certificates-only run verifies the short words but
explicitly does not establish exhaustive coverage.

This is internal exact computation and mathematical checking.
Formal verification, external human review, independent external
replay and priority are not claimed.

