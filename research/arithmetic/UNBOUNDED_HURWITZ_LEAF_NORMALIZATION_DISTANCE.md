# Unbounded cost of wreath-leaf normalization

A bounded local buffer cannot normalize every polynomial DT tuple to a
wreath leaf. For every \(r\ge4\), there is an explicit full
\(A_{2r+1}\) tuple for which any such braid word has length at least
\[
                         \left\lfloor r/2\right\rfloor-1.
\]
The initial tuple has a symmetric deletion at every leaf.

The obstruction is an odd transitive subgroup that survives deletion.
Several such subgroups can be preserved independently. One elementary
Hurwitz move can disturb at most one of them. This proves a uniform
lower bound on a particular normalization strategy; it leaves algorithms
with linearly growing or larger cost entirely open.

## An odd orbit cannot preserve a pairing

**Lemma 1.** Suppose a permutation group \(J\) is transitive on an odd set
\(\Omega\), \(|\Omega|>1\), and fixes every other letter. Then \(J\)
cannot preserve a partition of the whole letter set into pairs.

**Proof.** If a pair is \(\{a,b\}\) with \(a\in\Omega\), \(b\notin\Omega\),
an element of \(J\) moving \(a\) gives a different pair intersecting it
at the fixed letter \(b\). This contradicts the disjointness of blocks.
Thus every pair meeting \(\Omega\) is contained in \(\Omega\). This would
partition an odd set into pairs. ∎

Consequently, if a leaf \(x\) is outside \(\Omega\), deletion at \(x\)
cannot have group \(C_2\wr S_r\) whenever the retained factors contain
such a \(J\). The odd subgroup is unchanged by that deletion.

## The explicit family

Fix \(t\ge2\). Use the common root \(0\) and disjoint four-element sets
\[
        (a_i,b_i,c_i,d_i)=(4i+1,4i+2,4i+3,4i+4),
                          \quad0\le i<t.
\]
Define the two factors of the \(i\)-th lobe by
\[
        A_i=(0\,a_i)(b_i\,c_i),\qquad
        B_i=(a_i\,b_i)(c_i\,d_i).
 \tag{1}
\]
Its lobe support is \(\Omega_i=\{0,a_i,b_i,c_i,d_i\}\).
The subgroup \(J_i=\langle A_i,B_i\rangle\) is dihedral of order ten,
transitive on \(\Omega_i\), and fixes its complement pointwise:
the two generators are involutions and their product is a 5-cycle.

There are two versions:
\[
 \mathcal T_t=(A_0,B_0,\ldots,A_{t-1},B_{t-1}),
       \quad N=4t+1,\quad r=2t,
 \tag{2}
\]
and
\[
 \mathcal T_t^+=(A_0,B_0,\ldots,A_{t-1},B_{t-1},E),
 \quad E=(0\,y)(a_0\,z),\quad
 y=4t+1,\ z=4t+2,
       \quad N=4t+3,\quad r=2t+1.
 \tag{3}
\]
The buffer uses two fresh letters and two distinct old anchors.

**Proposition 2.** Each tuple in (2)--(3) has a full \(N\)-cycle product
and generates \(A_N\). Every initial leaf deletion has group \(S_{N-1}\).
Tuple (2) has no removable DT. Tuple (3) has precisely one, the buffer.

**Proof.** The expanded transposition graph of each lobe is the path
\[
                          0,a_i,b_i,c_i,d_i.
\]
The lobes meet only at zero, so their union is a tree. The buffer adds
two leaf edges, retaining a tree. Multiplying its distinct edge
transpositions in any order gives a full cycle. Indeed, at each step
an edge joins two components of the preceding forest, and hence merges
their two permutation cycles.

The generated group is transitive. In odd degree, any nontrivial block
size is at least three. A DT moving any block would have to move at least
two whole blocks, hence at least six letters. All the DT generators would
therefore preserve every block, contradicting transitivity. Thus the
group is primitive. Factors from two different lobes satisfy
\[
                           (A_0A_1)^2=(0\,a_0\,a_1).
 \tag{4}
\]
This is a 3-cycle with at least three fixed letters. Jordan's theorem,
in the retained [Jones formulation, Theorem 1.1](sources/jones-primitive-cycle-1209.5169v1.pdf),
gives \(A_N\), and all generators are even.

The common root is moved by at least two different lobe factors, so is
not a leaf. Any other leaf is in at most one lobe support. Another
\(J_i\) therefore avoids it and survives deletion, so Lemma 1 excludes
a wreath deletion. The
[deletion dichotomy, Theorem 1](IMPRIMITIVE_LEAF_DELETION_NORMAL_FORM.md)
then forces \(S_{N-1}\).

In (2) the only leaves are the \(d_i\). Their owners \(B_i\) have the
other edge \((a_i\,b_i)\) with no leaf endpoint. In (3) the same is true
for the old factors, while \(y,z\) are leaves on the buffer's two
different edges. Hence exactly that buffer is removable. ∎

## Tracking the preserved odd groups

Keep labels on the original factor colors. Let \(Q_j\) be an ordinary
positive elementary Hurwitz move:
\[
               Q_j(a,b)=(b,b^{-1}ab),\qquad
               Q_j^{-1}(a,b)=(aba^{-1},a).
 \tag{5}
\]
Colors travel with their factors. A positive move conjugates only the
first color; a negative move conjugates only the second color.

Call a lobe **exposed** when one of its two colors has been conjugated
by a color outside that original two-color lobe. Internal moves between
its two colors do not expose it: they preserve the generated group
\(J_i\) literally, even if they alter its individual generators.
As long as a lobe is unexposed, its two current factors therefore
still generate its original odd transitive group \(J_i\).

One elementary move exposes at most one previously unexposed lobe.
Moving a buffer color by conjugation exposes no old lobe; conjugating
an old lobe color by the buffer can expose only that lobe.

**Theorem 3.** For either tuple (2) or (3), every braid word whose
endpoint admits any wreath leaf deletion has length at least \(t-1\)
in the generators \(Q_j^{\pm1}\). At least \(t-1\) different original
lobes must have received an external conjugation.

For (3), if the endpoint has wreath deletion at the protected fresh
letter \(y\), the word has length at least \(t\), and every old lobe
must have received an external conjugation. The statement applies to
the fresh letter \(z\) as well.

**Proof.** If two lobes remain unexposed, their groups are still odd
transitive groups on their original supports. Their common root cannot
be a leaf: each group must have at least one of its two factors moving
that root, giving two different moving factors. Every other leaf
belongs to at most one of the two supports. Thus a surviving odd
transitive group avoids the selected leaf and is unchanged by deletion.
Lemma 1 again forbids a pairing, so the deletion is not wreath.

An endpoint with a wreath leaf can therefore have at most one
unexposed lobe. At least \(t-1\) lobes must have been exposed, which
needs at least \(t-1\) elementary moves.

For a protected fresh leaf, even one unexposed old lobe is an
obstruction: all old supports avoid \(y,z\). All \(t\) lobes must
therefore have been exposed. The word need not keep the protected
letter a leaf at intermediate stages. ∎

The theorem allows arbitrary sheet relabelings when comparing absolute
tuples. Common conjugation commutes with all Hurwitz moves and preserves
the wreath-deletion property, so it can be moved to the end of a word
without changing its length. The proof can therefore be carried out
with fixed sheet labels. For the stronger protected-fresh-letter bound,
the selected letter is marked and is transported by every relabeling;
one may not rename a different old leaf as the protected letter.
The unmarked t-1 bound permits arbitrary common relabelings.
Free relabeling under these respective conventions cannot evade the bound.

**Corollary 4.** For every \(r\ge4\), a full \(A_{2r+1}\) polynomial DT
tuple exists with symmetric deletion at every initial leaf, for which
any braid word to a wreath-leaf representative has length at least
\(\lfloor r/2\rfloor-1\).

There is no repair of uniformly bounded word length, nor one confined
to uniformly bounded numbers of original lobe colors, that always
reaches that locus. Adding one fresh DT buffer does not remove this
obstruction. Arbitrarily many internal braids within unexposed lobes
are permitted in this assertion.

**Proof.** Take \(t=\lfloor r/2\rfloor\) and choose (2) for even \(r\),
(3) for odd \(r\). Apply Proposition 2 and Theorem 3. ∎

This is not a distinction between unrestricted Hurwitz orbits. It does
not exclude a repeated local rule that propagates through a growing
number of lobes, a procedure with \(O(r)\) or larger cost, a different
normal form, or a route through the marked 3-cycle polynomial type.
It does not assume that the wreath locus is reachable in a degree
where connectivity is still unclassified; unreachable targets may
equivalently be assigned infinite distance.

## Exact controls

The [checker](verify_hurwitz_leaf_distance.py) verifies the family
identities for \(t=2,3,4,5,12,25\), including the literal 3-cycle certificate
and every initial leaf deletion. It exhausts every signed ordinary
word through depth \(t-2\) for both variants and \(2\le t\le5\), checking
all endpoint leaves and the actual unexposed odd groups. Separate
complete word balls through depth \(t-1\), \(2\le t\le4\), test the
protected fresh-leaf bound.

The deletion detector uses the unique pairing preserved by a full
\(2r\)-cycle: its antipodal pairs. For a full cycle Q, the stabilizer of a two-element
block in the cyclic group is its unique subgroup of order two, generated
by Q^r; hence the partner of a is Q^r(a). This is independent of a sampled
group-generation test. A positive control accepts the retained wreath
normal form; the degree-five one-lobe control has group \(D_5\), so is
not incorrectly claimed to be full alternating. An internal-braid
control changes the lobe generators while retaining its unexposed
group. The [receipt](hurwitz_leaf_distance_verification.json) and
[source ledger](hurwitz_leaf_distance_source_ledger.json) record the
actual finite replay. The unbounded lower bound rests on the written
group argument. No novelty, external human review, independent external
replay or formal verification is asserted.
