# Quadratic Hurwitz distance for the D5 bouquet family

The canonical wreath-leaf normalization of the unbuffered D5
bouquet requires quadratically many ordinary Hurwitz moves.
A historical support interval for each sheet gives a quadratic
lower bound, uniform over every possible endpoint leaf.
The root-preserving bounded-core absorption theorem supplies
a matching quadratic upper bound.

Thus a linear ordinary-word algorithm for this explicit family
is impossible. The family nevertheless lies in the canonical
wreath-leaf orbit in every degree under consideration.
This is a distance theorem within that orbit, not an orbit distinction.

## The family and metric

For \(t\geq2\), on \(4t+1\) sheets put
\[
 \begin{split}
 \mathcal T_t&=(A_0,B_0,\ldots,A_{t-1},B_{t-1}),\\
 A_i&=(0\,a_i)(b_i\,c_i),\qquad
 B_i=(a_i\,b_i)(c_i\,d_i),\\
 (a_i,b_i,c_i,d_i)&=(4i+1,4i+2,4i+3,4i+4).
 \end{split}                                                \tag{1}
\]
The [lobe theorem](UNBOUNDED_HURWITZ_LEAF_NORMALIZATION_DISTANCE.md)
proves that these polynomial tuples have full group \(A_{4t+1}\)
and a full-cycle product. Every initial leaf deletion is symmetric.

A **wreath leaf** of a tuple is a sheet moved by exactly one factor,
such that deleting it from that double transposition gives the
natural full \(C_2\wr S_{2t}\) action on the remaining \(4t\) sheets.
Let \(d_t\) be the minimum ordinary signed Hurwitz word length taking
\(\mathcal T_t\) to a tuple with some wreath leaf. Common sheet
relabelings are allowed for free. Arbitrary permutations of branch
positions are not free: they must be effected by adjacent generators.

**Theorem 1.** Put
\[
                     g(t)=\max\{0,\,2\lceil t/2\rceil-3\}.
\]
Then
\[
 \boxed{
 \left\lceil\frac{(2t-1)g(t)}8\right\rceil
 \ \leq d_t\ \leq\
 6656+(t-2)\,446265984+14t^2-52t+48 .}                        \tag{2}
\]
The lower bound holds for every endpoint wreath leaf, including
one at the original common root0. In particular
\[
                 d_t=\Theta(t^2).
\]
Since the number of factors is \(r=2t\), the lower bound also gives
\(d_t\geq r^2/16-O(r)\). It excludes a total \(O(r)\)
normalization procedure in this ordinary adjacent-generator metric.

## 1. Historical intervals move at bounded speed

Fix the sheet labels. Along a braid word, let \(I_a\) be the smallest
interval of branch positions at which sheet \(a\) has ever been moved
by a factor, including the initial tuple. Every sheet of (1) occurs
initially, so these intervals are nonempty. Their lengths can only
increase. Write
\[
       H=\sum_a\bigl(|I_a|-|I_a^{\,0}|\bigr),
 \qquad |[i,j]|=j-i,
\]
where \(I_a^{\,0}\) is the initial occurrence interval.

**Lemma 2.** Each ordinary signed Hurwitz move increases \(H\)
by at most8. Consequently a word of length \(\ell\) has \(H\leq8\ell\).

**Proof.** The move affects adjacent positions \(i,i+1\), occupied
by double transpositions \(a,b\). Both new factors fix every sheet
outside \(\operatorname{supp}(a)\cup\operatorname{supp}(b)\),
a set of size at most8.

Any affected sheet already occurs in at least one of positions
\(i,i+1\). A new occurrence at the other position can therefore
enlarge its historical interval by at most1. If it already occurred
at both positions, neither is outside its historical interval.
No other position changes. This proves the bound for each sheet
and hence for their sum. The argument applies equally to
\(Q_i(a,b)=(b,b^{-1}ab)\) and its inverse. \(\square\)

**Lemma 3.** Suppose two initial occurrence intervals \(I_a^0,I_b^0\)
are separated by a gap \(q\geq0\). If the two historical intervals
eventually intersect, their combined length increase is at least \(q\).

**Proof.** If \(I_a^0\) ends at \(u\) and \(I_b^0\) begins at
\(u+q\), intersection requires the rightward extension of the first
plus the leftward extension of the second to be at least \(q\).
Their total growth is at least that sum. \(\square\)

## 2. Every endpoint block has a shared factor

Let \(P\) be the fixed product of (1). Suppose an arbitrary endpoint
has wreath leaf \(x\). Delete the unique edge containing \(x\)
from its owner factor. The lower tuple has one transposition,
\(2t-1\) double transpositions, and product the \(4t\)-cycle
obtained by omitting \(x\) from \(P\).

The deletion formula holds at any factor position. Indeed if
the upper product is \(L(x\,u)eR\), all other factors fix \(x\).
Its lower product is \(Q=LeR\), and
\[
                          P=(x\,L(u))Q.
\]
The displayed transposition inserts \(x\) into the nontrivial
cycle of \(Q\); omitting \(x\) recovers it. This is also the
[fixed-leaf deletion identity](FIXED_LEAF_HURWITZ_OBSTRUCTION.md).

A full even cycle has only one invariant partition into pairs:
its antipodal pairs. To see this, write a pair as
\(\{a,P_x^j(a)\}\). Translating it by \(P_x^j\) gives another
block meeting the original one, so these blocks are equal.
Thus \(P_x^{2j}(a)=a\) and \(j=2t\).

Now project the lower tuple to its \(2t\) two-sheet blocks.
The single transposition fixes every block. A double transposition
can induce either the identity or one block transposition.
The projected product is a \(2t\)-cycle, which requires at least
\(2t-1\) transpositions. Since exactly that many double
transpositions are available, every one induces a block swap.
Their projected edges form a connected tree.

In particular each block is incident with some double transposition,
which moves both its sheets. These factors are unchanged by deletion.
Thus for every lower block \(\{a,b\}\), both sheets occur jointly
in some endpoint branch position, and
\[
                              I_a\cap I_b\ne\varnothing.    \tag{3}
\]
No assumption is made about which sheet \(x\) is chosen, or about
it remaining a leaf during the word.

## 3. Antipodal partners start far apart

The product of the two factors in each lobe is
\((0\,a_i\,c_i\,d_i\,b_i)\). Hence the full product is
\[
 P=(0,\ a_{t-1},c_{t-1},d_{t-1},b_{t-1},\
          \ldots,\ a_0,c_0,d_0,b_0).                         \tag{4}
\]
The nonroot sheets of each lobe form four consecutive letters of
this cycle. The root contributes one additional letter between
the last and first blocks.

Consider any antipodal pair of the lower cycle with \(x\) omitted.
Its two directed distances in the original \(P\) are \(2t\) and
\(2t+1\): exactly one of the two arcs contains the deleted sheet.
Suppose neither member of the pair is the original root0,
and their lobe indices are \(i,j\). Put
\[
                        q=\min(|i-j|,t-|i-j|).
\]
There is an arc between them crossing \(q\) lobe transitions.
Its length is at most \(4q+4\): each full transition accounts
for four letters, the endpoint offsets add at most3, and the
root can add at most1. Therefore
\[
                     2t\leq4q+4,\qquad
                         q\geq\lceil t/2\rceil-1.           \tag{5}
\]

Initially every nonroot sheet in lobe \(i\) appears only at
positions \(2i+1,2i+2\). Consequently the two initial occurrence
intervals of these partners have gap at least
\[
                 \max(0,\,2|i-j|-1)\geq g(t).              \tag{6}
\]
Of the \(2t\) lower pairs, at most one contains the original root.
If the deleted sheet itself is0, none does. Thus there are always
at least \(2t-1\) disjoint pairs to which (6) applies.

For each of them, (3), Lemma3 and (6) force combined historical
growth at least \(g(t)\). The pairs have disjoint sheets, so their
contributions add:
\[
                              H\geq(2t-1)g(t).              \tag{7}
\]
Lemma2 now gives the lower bound in (2).

This argument is unchanged by free simultaneous sheet relabelings.
They commute with every Hurwitz move and may be pushed to the end.
Pulling the endpoint labels back gives exactly the fixed-product
pairing argument above. This does not make branch-position
permutations free.

## 4. Matching upper bound and remaining scope

The [root-preserving bounded-core absorption theorem](D5_LOBE_BOUNDED_CORE_ABSORPTION.md),
Corollary2, takes the entire bouquet (1) to the canonical
\(W_{2t}\) tuple with the upper bound in (2). That construction
uses one fixed degree-nine and \(t-2\) fixed degree-thirteen
controllers whose existence follows from complete finite
connectivity theorems. Their literal words have not been recovered;
their explicit length bounds and all endpoint conversions are proved.

The upper bound is \(O(t^2)\), since its controller constants
are fixed independently of \(t\). The lower bound is
\(t^2/4-O(t)\), proving the claimed quadratic order.
This proves both reachability and failure of a linear ordinary-word
algorithm for this family.

The result does not classify arbitrary polynomial DT tuples,
prove or disprove universal braid connectivity, treat the fresh
two-sheet buffered odd-factor variant, or establish an optimal
leading constant for \(d_t\). It permits arbitrary endpoint
wreath leaves and free common sheet relabelings throughout.

## Exact evidence

The [checker](verify_quadratic_hurwitz_bouquet_distance.py)
checks the local historical-span statement for all210
double transpositions on eight letters, all ordered pairs,
and both move signs:88200 cases. Any pair of DTs in any degree
has support at most8, so this covers every local support pattern.
The written interval argument handles arbitrary previous history.

For every \(2\leq t\leq40\), it checks the product (4),
every possible omitted sheet, every lower antipodal pair,
the original cycle distances, lobe separation and the initial
interval gaps. Complete short word balls also update the actual
historical intervals and check every wreath-leaf endpoint found.
The [receipt](quadratic_hurwitz_bouquet_distance_verification.json)
pins these controls and the imported upper theorem.

The uniform result rests on the written argument. The bounded
controls do not substitute for the all-leaf quantifiers, the
historical-growth proof, or the imported finite-controller existence.
There is no formal verification, external human review,
independent external replay or novelty claim.
