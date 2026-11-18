# Natural D5 lobe trees lie in one canonical Hurwitz orbit

An ordered family built by attaching five-sheet D5 actions along one
existing sheet can be normalized uniformly. The shared sheet may be
any existing sheet, and it may be an interior vertex of the new
lobe's transposition path. A partial full twist moves the old
attachment to the canonical root; at most two internal moves put
the new pair in the rooted form. The established bounded-core
absorption then applies.

This gives reachability for the whole stated family in every
degree, with a coarse quartic ordinary-word bound. It does not
classify arbitrary double-transposition trees.

The full-twist conjugation used in Lemma 3 is a classical Hurwitz-action
identity; compare Bessis, *The dual braid monoid*,
[arXiv:math/0101158v2, Sections 1.5--1.6, printed pp.6--8](https://arxiv.org/pdf/math/0101158v2).
Our generators are the inverses of Bessis's convention, giving conjugation
by the inverse total product, as derived below. His Proposition 1.6.1
also gives transitivity for reduced reflection factorizations of a
parabolic Coxeter element. Expanding our double transpositions gives a
reduced transposition factorization of a full cycle, but that theorem
does not preserve the adjacent commuting pairs needed to return to
double-transposition moves. It therefore credits the full-twist mechanism
and the expanded-factorization background; it does not replace the
pair-preserving absorption argument, finite controllers, or word bound.

## The ordered one-vertex class

A **natural D5 pair** is an ordered pair \((a,b)\) of double
transpositions supported on a five-element set with \(ab\) a 5-cycle.
It generates the natural dihedral group of order10.

Start with one such pair on five sheets. At each later step append
one natural D5 pair whose support contains exactly one existing
sheet \(v\) and four fresh sheets. The shared sheet \(v\) can be
any existing sheet and any vertex of the new pair's support.
After \(t\) pairs there are \(2t\) factors on \(4t+1\) sheets.
The order of these adjacent pairs is part of the construction.
Call this an **ordered one-vertex D5 lobe tree**.

The term tree describes the expanded transposition graph: the two
matchings of each pair form a four-edge path, and each new path is
glued at exactly one vertex. It does not require the support
intersection graph of the lobes to be a simple tree; a bouquet
with several lobes at the same sheet is included.

Let \(W_r\) be the canonical wreath-leaf tuple
\[
 W_r=\bigl((2r\,\,2)(01),\
 (0\,\,2r-2)(1\,\,2r-1),\ldots,(02)(13)\bigr)
\]
from the [wreath normal-form theorem](IMPRIMITIVE_LEAF_DELETION_NORMAL_FORM.md).

**Theorem 1.** For every \(t\geq2\), every ordered one-vertex
D5 lobe tree is taken to \(W_{2t}\), up to simultaneous sheet
relabeling, by an ordinary signed Hurwitz word of length at most
\[
 6656+(t-2)M+
       \sum_{j=2}^{t-1}
       \left((2j)^2(2j-1)+14(2j)-38+2\right),
 \qquad M=446265984.                                      \tag{1}
\]
In particular this is \(O(t^4)\). For each fixed \(t\), the entire
stated class lies in one absolute Hurwitz orbit and has full
group \(A_{4t+1}\).

The bound imports the already proved finite degree-nine and
degree-thirteen controller existence. Those fixed words have
not been recovered or replayed. Full-twist and pair-normalization
words are explicit. The checker replays all noncontroller moves
and substitutes only the certified finite-controller endpoints.

Our convention is
\[
 Q_i(a,b)=(b,b^{-1}ab),\qquad
 Q_i^{-1}(a,b)=(aba^{-1},a),
\]
and words act from left to right.

## 1. The new pair can be rooted at any shared sheet

**Lemma 2.** Let \((a,b)\) be a natural D5 pair, and let \(v\) be
any of its five sheets. There is a unique \(j\in\{-2,-1,0,1,2\}\)
for which the second factor of \(Q^j(a,b)\) fixes \(v\).
After this move, labeling only the four other sheets gives
\[
                         (v\,\alpha)(\beta\,\gamma),\qquad
                         (\alpha\,\beta)(\gamma\,\delta).   \tag{2}
\]

**Proof.** Put \(C=ab\), so \(b=aC\). Since \(aCa=C^{-1}\),
direct substitution gives
\[
                         Q^j(a,b)=(aC^j,aC^{j+1}).         \tag{3}
\]
This holds for signed \(j\), and \(Q^5\) fixes the literal pair.
The five reflections \(aC^j\) are distinct double transpositions.
They fix distinct sheets: otherwise their product would be a
nonidentity power of the 5-cycle \(C\) fixing that same sheet.
Thus each sheet is fixed by the second reflection in exactly
one of the five signed representatives.

The two matchings have no common edge, since such an edge would
give fixed sheets of their 5-cycle product. Their union is connected,
because their group is transitive. It has five vertices, four edges
and maximum degree2, so it is a path. When the second matching
fixes \(v\), that sheet is an endpoint of the path. Following
alternately the first and second matchings from \(v\) gives
\[
 v,\ \alpha=a'(v),\ \beta=b'(\alpha),\
          \gamma=a'(\beta),\ \delta=b'(\gamma).
\]
These five sheets are distinct and give exactly (2).
At most two ordinary moves are used. \(\square\)

If \(v\) is the only old sheet, the relabeling in this lemma
acts only on the four fresh sheets and fixes every other old
sheet. Future attachments at these four sheets are relabeled
with them; their one-old/four-fresh structure is preserved.

## 2. An arbitrary old attachment can be moved to the root

**Lemma 3.** Suppose an old prefix is \(W_r\), and the appended
natural D5 pair shares an arbitrary old sheet \(v\).
Using at most \(r^2(r-1)\) ordinary moves on the old prefix,
followed by a common sheet relabeling fixing the four fresh sheets,
the old prefix returns to \(W_r\) and the shared sheet becomes0.
No appended or later factor is touched by the ordinary word.

**Proof.** Let \(P=\prod W_r\), a full \((2r+1)\)-cycle.
The cycling word \(Q_1\cdots Q_{r-1}\) sends a tuple of product
\(P\) to
\[
                  (g_2,\ldots,g_r,P^{-1}g_1P).
\]
Its \(r\)-th power \(\Delta_r^2\), of length \(r(r-1)\),
therefore conjugates every old factor by \(P^{-1}\).

Choose the signed exponent \(-r\leq k\leq r\) satisfying
\(P^k(0)=v\). Apply \((\Delta_r^2)^{-k}\) to the first \(r\)
factor positions. Its endpoint is
\[
                  P^kW_rP^{-k},
\]
and every later factor is unchanged literally.
Extend \(P^{-k}\) to fix all fresh sheets and relabel the
entire tuple by it. The old prefix returns to \(W_r\), while
the shared sheet becomes \(P^{-k}(v)=0\).
The cost is at most
\[
                         |k|r(r-1)\leq r^2(r-1).
\]
Future attachments are conjugated but remain attachments at
existing sheets. Every still unused sheet is fixed. \(\square\)

After Lemma3, apply Lemma2 to the appended pair and relabel only
its four fresh sheets to the next four labels. The prefix now
has exactly the root0 form required by the
[root-preserving absorption theorem](D5_LOBE_BOUNDED_CORE_ABSORPTION.md).
That theorem costs at most \(M+14r-38\) and returns the prefix
to \(W_{r+2}\). No common-root constraint on subsequent lobes
is needed in the present class.

## 3. The first two lobes and the induction

The expanded transposition graph of every ordered lobe tree is
a tree. A product of its edge transpositions, in any order, is a
full cycle: each new edge joins two components of the preceding
forest and merges their two permutation cycles. In particular
the first two lobes have full 9-cycle product.

They also generate \(A_9\). Choose one reflection from each lobe
moving their shared sheet. Their four-sheet supports intersect
only at that sheet, so the square of their product is a 3-cycle.
The group is transitive. In odd degree a nontrivial block has
size at least3; a DT cannot move two such blocks, so transitivity
forces primitivity. Jordan's 3-cycle consequence in
[Jones, Theorem1.1](sources/jones-primitive-cycle-1209.5169v1.pdf)
gives \(A_9\).

The explicit complete degree-nine premise in
[the absorption theorem, Corollary2](D5_LOBE_BOUNDED_CORE_ABSORPTION.md)
applies to every four-DT tuple of full 9-cycle product, not only
to the common-root bouquet. The independent literal total is59049,
the free cyclic quotient has6561 states, and the retained complete
BFS reaches all of them. An absolute path of length at most6560
and at most eight four-strand full twists of length12 give a
literal base controller of length at most6656.

Thus the first two lobes can be normalized to \(W_4\), with all
later factor positions untouched by the ordinary word. Extend
the endpoint relabeling to fix all unused sheets. This merely
updates the later attachment labels and retains the ordered
one-vertex construction.

Now suppose \(j\) lobes have been normalized. The old factor count
is \(r=2j\). Lemma3, Lemma2 and root0 absorption handle the next
lobe at cost at most
\[
                    (2j)^2(2j-1)+2+M+14(2j)-38.
\]
At every stage, relabelings act within the active sheet set and
fix unused sheets. Consequently each later pair still shares
exactly one active sheet when it is reached. Summing these
bounds proves (1).

Common relabeling commutes with every Hurwitz generator, so
all intermediate relabelings may be pushed to a single final
one without changing the ordinary word length.
The endpoint is \(W_{2t}\), whose group is \(A_{4t+1}\) by
the wreath normal-form theorem. This also proves full alternating
group for every original tuple in the class. \(\square\)

## Exact controls and scope

The [checker](verify_d5_lobe_tree_reachability.py) verifies all120
ordered natural D5 pairs on five sheets and every choice of
shared sheet:600 local cases. It checks the literal period5,
the unique signed exponent of absolute value at most2,
the full product, and the alternating path formula.

It checks old-root transport for every old sheet of \(W_r\)
for \(4\leq r\leq8\), including both signs of the full twist.
Whole-construction controls use four branching/chain patterns
at \(t=2,3,4,5,8,12\), with internal shared vertices in the new
lobes and all future labels transported. They verify the
canonical endpoint and the stated cost.
The [receipt](d5_lobe_tree_reachability_verification.json)
records actual noncontroller word replays and explicitly
declares the finite-controller endpoint substitutions.

Uniformity follows from Lemmas2-3 and the construction-preserving
induction, not the finite control range. The quartic bound is
coarse, and no optimal cost is asserted. The [common-root bouquet](QUADRATIC_HURWITZ_BOUQUET_DISTANCE.md)
already has a quadratic lower bound, so a universal linear
ordinary-word bound for this class would be impossible.

This is the stated ordered one-vertex class. It does not cover
arbitrary DT tuples, overlapping new lobes at several old sheets,
arbitrary branch interleavings of an untracked lobe decomposition,
or general extensions by one DT. It does not settle the next
unclassified degree16 original passport or uniform all-tuple
Hurwitz connectivity. No formal verification, external human
review, independent external replay or novelty claim is made.
