# Normal form and extraction for imprimitive leaf deletions

A whole arbitrary-input case of the polynomial DT connectivity problem
has a uniform normal form. If deleting a chosen leaf gives an imprimitive
group, that group is exactly C2 wr S_r. An explicit braid and sheet
relabeling put every such tuple into the same star form. This controls
all attachment choices in this case and produces r-2 removable DTs.

The complementary deletion group is always S_(2r). That symmetric
case remains the frontier. The result does not assert that every
full-alternating orbit reaches the imprimitive-deletion case.

## Setting and the deletion dichotomy

Let a tuple of r DTs on N=2r+1 symbols have full-cycle product P.
Choose a leaf x: exactly one factor moves it. Write that factor

\[
                         A=(x\,u)e,\qquad e=(v\,w).
\]

Delete x by replacing A with e, as in the
[fixed-leaf obstruction](FIXED_LEAF_HURWITZ_OBSTRUCTION.md).
The lower tuple has degree 2r, one transposition and r-1 DTs,
and full-cycle product. Its group H is transitive.

## Theorem 1. There are exactly two possible deletion groups

For r>=2, the lower group is either S_(2r) or the natural full
wreath product C2 wr S_r on r blocks of size two.

In the imprimitive case the block system is unique and can be
found without enumerating the entire group: take the orbit of
the support pair of e under the lower generators, regard its
pairs as graph edges, and take connected components.

**Proof.** If H is primitive, conjugates of its transposition e
give a nonempty invariant edge graph. Components are blocks, so
the graph is connected. Its edge transpositions generate S_(2r),
which proves H=S_(2r).

Suppose H is imprimitive with block size d. A transposition cannot
move a nontrivial block. If d>=3, a DT cannot move blocks either,
because moving at least two such blocks would move at least six
symbols. Then every generator would fix every block setwise,
contradicting transitivity. Hence d=2.

The lower full cycle induces a full r-cycle on blocks. Each DT
induces either the identity or a transposition. There are r-1 DTs,
and an r-cycle needs at least r-1 transpositions; therefore all
project to transpositions. In particular each is a block swap,
with its four moved symbols exactly those of its two blocks.
The r-1 projected edges form a connected tree.

Each block swap is parallel or crossed after labeling the two
symbols of each block. Orient a root block arbitrarily and
propagate orientations along the projected tree so every swap
becomes parallel. The tree has no cycle, so these choices are
consistent. The parallel swaps generate a diagonal S_r because
their projected tree is connected. Conjugating the single block
flip e by it supplies independent flips in all blocks. Thus H
is the whole C2^r semidirect S_r.

In this wreath group the transpositions are precisely its r
individual block flips. Their conjugacy orbit graph has the r
blocks as its components. In the symmetric case the graph is
connected. This proves the algorithm and uniqueness. ∎

## Lemma 2. A constructive transposition star normalization

A reduced transposition factorization of a full cycle C on r symbols
can be braided to the unique star factorization centered at any
chosen symbol b.

If C=(b,c_1,...,c_(r-1)), that star is

\[
                   ((b\,c_{r-1}),\ldots,(b\,c_1)).
\tag{1}
\]

**Proof and algorithm.** The input edge graph is a tree. Choose a
leaf l different from b. The target star also has l as a leaf.
In each factorization move its unique l edge to the final position
by positive Hurwitz moves. Every crossed factor remains unchanged;
the chosen edge is conjugated but still contains l. The common
product C forces the final factor to be

\[
                              (l\,C^{-1}(l)).
\]

Deleting this last factor and the symbol l gives two reduced
factorizations of the same cycle with l omitted. The target one
is the star centered at b on the remaining symbols. Apply induction,
then undo the target's bringing word. This is an explicit finite
algorithm. The cases r=1,2 are immediate. ∎

This classical transposition argument is used only after passage to
the two-element blocks. Its lift preserves every paired DT: a conjugate
of a block swap is again a block swap, and all such lifts fix x.

## Theorem 3. A uniform normal form for every imprimitive deletion

Let r>=3 and suppose the chosen-leaf deletion is imprimitive.
Through elementary Hurwitz moves that keep x a leaf, followed by
a simultaneous relabeling carrying x to 2r, the entire tuple becomes

\[
 \boxed{
 A_0=(2r\,\,2)(01),\qquad
 D_j=(0\,\,2j)(1\,\,2j+1),\quad j=r-1,r-2,\ldots,1.}
\tag{2}
\]

Thus all tuples admitting such a chosen leaf are in one absolute
Hurwitz orbit. For a fixed product P and fixed chosen sheet x,
they are connected literally while x remains a leaf.

**Proof.** First move the unique x factor to the first position
by inverse moves

\[
                        (B,A)\longmapsto(BAB,B).
\]

The B factors fix x, so this always retains the leaf. Both the
lower flip block and its attachment block are conjugated by B;
they remain distinct. Let R be the current flip block and U
the block containing the attachment u.

Braid only the r-1 other factors. Their projections are a reduced
transposition factorization of a full r-cycle C. Apply Lemma2
with center R. Lift exactly the same word to the block swaps.
Every intermediate factor is still a DT fixing x, and A is
unchanged.

Put R'=C^-1(U), so R' differs from U. If R'=R, no adjustment
is needed. Otherwise the star has a unique swap B between
R and R'. Bring it next to A unchanged, using positive moves
among the block swaps, and apply Q_1^2 to (A,B).

The swap B fixes the attachment u, because neither of its blocks
is U. The exact result is

\[
                       (BAB,\,(BAB)B(BAB)).
\tag{3}
\]

At both intermediate stages x remains a leaf. The first factor's
flip block moves from R to R'; its attachment u is unchanged.
The second factor fixes x and is a block swap with the same
projected transposition as B. Thus U and the projected product C
have remained fixed.

Apply Lemma2 again to make the block tree a star centered at R'.
Now its flip block is the predecessor of U in C. Label the
blocks cyclically along C starting at R', so the root is block0
and U is block1. Orient the root arbitrarily and use the star
swaps to orient every other block in parallel. Reversing all
orientations simultaneously, if needed, makes u the even symbol
2 in block1. The flip is (01).

The projected star factorization of the labeled cycle
(0 1 ... r-1) is uniquely (1). Hence its parallel lifts and A
are exactly (2). Every step has been specified constructively.

For the literal fixed-product assertion, let P_0 be the product
of (2). A relabeling taking P to P_0 and x to 2r is unique:
the centralizer of a full cycle is its cyclic group, acting
regularly on the symbols. Consequently two inputs with the same
P,x have the same final relabeling. Their pre-relabeling endpoints
are identical, so their actual leaf-preserving braid paths meet.
If x is not fixed, conjugation by powers of P moves its label
transitively and is realized by the central full twist, as in the
[unmarked reduction](UNMARKED_POLYNOMIAL_HURWITZ_REDUCTION.md).
Thus the imprimitive-deletion locus lies in one literal fixed-product
braid orbit as well. ∎

The statement controls the complete imprimitive-deletion locus,
including its extension and attachment choices. It does not say
that this locus is invariant under unrestricted Hurwitz moves;
a tuple in its orbit can have only symmetric deletions in another
representative.

## Corollary 4. Full alternating monodromy and recursive extraction

Every tuple in Theorem3 generates A_(2r+1). Its normal form has
r-2 simultaneous removable DTs, corresponding to blocks
2,...,r-1. Deleting one of these blocks leaves the same normal
form with r replaced by r-1.

For r>=4 the remaining group is A_(2r-1). Thus this arbitrary-input
case admits repeated extraction with full alternating cores down
to degree seven. The next degree-five core has group D5.

**Proof.** The upper product is a full cycle, so its group is transitive.
The odd-degree small-support block argument makes it primitive.
For a literal 3-cycle, use from (2)

\[
 A=A_0,\quad B_1=(02)(13),\quad B_2=(04)(15),\quad C=AB_1.
\]

On the seven labels 0,1,2,3,4,5,x, with x=2r,

\[
 C=(0\,x\,2\,1\,3),\qquad
 CAC^{-1}=(12)(x3),\qquad
                       ((CAC^{-1})B_2)^2=(1\,2\,5).
\tag{4}
\]

Jordan's prime-cycle theorem gives A_(2r+1), including r=3.
The precise imported formulation is retained in the
[primary-source ledger](alternating_hurwitz_source_ledger.json).

For j>=2, only D_j moves either symbol 2j,2j+1. These are
exclusive leaf ends on its two edges. Delete D_(r-1) and its
two symbols, and relabel the old x=2r as the new x=2r-2.
The remaining tuple is literally (2) for r-1.
Apply the same group proof while r-1>=3. At r-1=2, two
involutions with 5-cycle product generate the dihedral group
of order10. ∎

## The remaining symmetric-deletion gate

Within the full-alternating polynomial DT space, all inputs
admitting an imprimitive leaf deletion lie in the single orbit
of (2). Therefore unrestricted connectivity is equivalent to
the following precise remaining assertion:

\[
 \text{Every full-alternating orbit has a representative with
 an imprimitive chosen-leaf deletion.}
\]

No proof of that assertion is supplied here. The dichotomy
makes its complement exact: in any representative outside this
locus, every leaf deletion has group S_(2r).

The degree-seven input

\[
                 ((01)(34),(12)(45),(23)(56))
\tag{5}
\]

already shows that one cannot simply choose an imprimitive
leaf in an arbitrary full-alternating tuple. It has full group
A7, no removable DT initially, and both leaf deletions have
group S6.

The complete degree-seven control is sharper. Of its two full-A7
orbits of147 tuples, one has126 imprimitive and294 symmetric
chosen-leaf inputs; the other has420 symmetric chosen-leaf
inputs and none imprimitive. After normalizing its full cycle by relabeling, (5) lies in the latter orbit.
The two proper order168 orbits have only symmetric deletions.

These are degree-seven boundary controls. They do not assert a
higher-degree orbit distinction. In the previously classified
higher cases, an input with only symmetric deletions can still
lie in the orbit of (2).

For the residual range \(N=2r+1\ge15\), the equivalent removable-factor
route is now organized by
[global terminal-mass minima](TERMINAL_MASS_HURWITZ_DESCENT.md#proposition-8-orbit-minima-replace-monotone-descent).
Theorem 10 and Corollary 11 there impose simultaneous complementary cyclic
covers, \(3m\ge N-4\), and equal terminal sizes or a singleton at every
minimizer of a hypothetical positive-minimum orbit. These necessary
conditions leave a full-alternating compatibility problem. In particular
a singleton terminal in such a cut is a leaf, but does not imply an
imprimitive deletion or membership in the normal-form locus.

## Exact verification

The [checker](verify_imprimitive_leaf_normal_form.py) implements
the normal-closure block finder, the recursive star word, its
actual DT lift, the root-adjustment word (3), and the complete
sheet relabeling. It checks all imprimitive leaf choices in
the complete degree-seven domain, along with broader ordered
and twisted block-tree inputs. Every intermediate upper move
is required to retain the chosen leaf.

The [receipt](imprimitive_leaf_normal_form_verification.json)
and [source ledger](imprimitive_leaf_normal_form_source_ledger.json)
pin the exact data. The all-degree result follows from the
written classification, lifting and normalization proofs.
No arbitrary symmetric-deletion normalization, complete
higher-degree census, full uniform connectivity theorem,
formal verification or priority claim is made.

