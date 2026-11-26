# Minimum dictionaries and exact scaling in triangle transfers

A critical collection of fixed-fibre triangles from a minimum-dictionary
abc graph produces a new minimum-dictionary family. One can preserve
the full exponent of an arbitrary selected triangle relation while
fixing the exact scaling gcd and the original residue classes.
After pinning one original endpoint, the transfer is injective.
This retains an explicit affine link between the two endpoint families.

## The minimum-dictionary input and the tagged relation

Use the height-first convention with \(n=X^L\) and write \(c_0=1/6\).
Take actual positive primitive triples \(a+b=c\) in the shell
\(X/2<c\le X\) and with \(\operatorname{rad}(abc)\le c^L\).
Suppose the family has normalized exponent \(2c_0-o(1)\)
and its three complete-letter dictionaries have exponents at most
\(c_0+o(1)\). The input has balanced individual radical caps
\(n^{2c_0+o(1)}\). These are hypotheses on the complete values
and the family. They allow arbitrary nonatomic mass in the
limiting degree measures.

Fix distinct roles \(T,W\) and call the remaining role \(V\).
A tagged triangle comprises a complete \(T\)-value and three
increasing complete \(W\)-values \(x_1<x_2<x_3\) occurring with it.
Each pair \((T,x_i)\) determines its original row.
Put \(h_{ij}=x_j-x_i\) and
\[
        g=\gcd(h_{12},h_{23}),\qquad D_{ij}=h_{ij}/g .
\]
The triple \((D_{12},D_{23},D_{13})\) is primitive positive.
Let \(\mathcal Q\) be any selected relation of such tagged triangles
of normalized exponent \(4c_0-o(1)\), the full tagged size.

## Deleting rare physical differences

Write \(N=n^{c_0+o(1)}\) for an upper bound on each original
dictionary size, and write \(A,B\) for the complete dictionaries
on \(W,V\). For \(h>0\), let \(M_h\) be the number of increasing
\(W\)-pairs of gap \(h\) on common \(T\)-fibres, including their
\(T\)-tags. Let \(r_A(h)\) and \(r_B(h)\) count pairs of gap \(h\)
in the two complete dictionaries. The original signed equation
implies
\[
        M_h\le r_A(h)r_B(h),\qquad
        \sum_h M_h\le N^3,\qquad
        \sum_h r_A(h),\ \sum_h r_B(h)\le N^2 .
\]
For the first inequality, the two \(V\)-values have difference
of absolute value \(h\) with orientation fixed by the signs.
The first \(W\)-value and its first \(V\)-value determine \(T\).
For any fixed \(0<\xi<c_0\), call a gap rare if
\(M_h<n^{2c_0-\xi}\). Then
\[
 \sum_{h\ \mathrm{rare}} M_h
   \le n^{c_0-\xi/2}\sum_h\sqrt{r_A(h)r_B(h)}
   \le n^{c_0-\xi/2}N^2 .
\]
After choosing one of the three pair positions, a tagged pair
has at most \(N\) possible third values. Thus the rare sides
account for only \(O(n^{4c_0-\xi/2+o(1)})\) tagged triangles.
The popular gaps number at most \(n^{c_0+\xi+o(1)}\).
For every fixed \(\xi\) the entire deletion has a power gap.
Sending \(\xi\) to zero diagonally with the other allowances
gives a relation at full exponent using only a gap dictionary
of exponent \(c_0+o(1)\).
[Earlier popularity extraction](ABC_POPULAR_SHIFT_RADICAL_AMPLIFICATION.md#a-bound-for-rich-shifts-with-a-radical-cap)

## Fixing the exact scaling and the three residue classes

The following conclusion does not require a cap on the raw
triangle radicals. For every sufficiently small fixed \(\xi>0\)
the relation of full exponent has a subclass of exponent at least
\(4c_0-O(\xi)-o(1)\) with one exact gcd \(g=G\) and all its
original rows in fixed residue classes on all three roles
modulo \(G\). In particular the exponent loss tends to zero
as \(\xi\) is subsequently decreased. The proof keeps the
selected triangles throughout.

For one modulus \(g\le X\) and residue \(\rho\) modulo \(g\),
the cell is the subset of the original \(W\)-dictionary in
that residue. Classify each tagged triangle by its *exact*
\(g=\gcd(h_{12},h_{23})\) and the cell containing all its values.
Let \(H=n^{c_0-\xi}\). A triangle with cell size below \(H\)
is counted by fixing \(T,x_1,x_3\), then a divisor \(g\mid h_{13}\),
then one of at most \(H\) middle values. These triangles number
at most \(N^3 H X^{o_X(1)}\). The divisor bound is subpower
since the budget is fixed at the height limit. A cell of
size at least \(H\) contains \(\gg H^2\) distinct ordered pairs.
Every ordered pair of distinct values can belong to only
divisor-many cells, by again dividing the difference.
There are therefore at most \(N^2 X^{o_X(1)}/H^2\) thick cells.
Selecting a largest exact cell from the nonrare triangle class
loses at most \(n^{2\xi+o(1)}\) in exponent. The choice fixes
both \(G\) and the original \(W\) residue.

Now an original \(T\)-residue determines an original \(V\)-residue,
bijectively, from the signed equation and the fixed \(W\) residue.
If the \(T\)-class has fewer than \(H'\) original values, the
number of compatible original edges in the paired residue
classes is bounded by \(H'\) times the size of the \(V\)-class.
Summing over these \(T\)-classes gives at most \(H'N\) edges.
Apply the same argument if the \(V\)-class has fewer than
\(H'\) values. Each original row extends to at most \(N^2\)
tagged triangles with itself among the three legs, up to a bounded
orientation factor. Take \(H'=n^{c_0-4\xi}\) to make this
second deletion smaller than the retained relation. The
remaining pairs of \(T,V\) residue classes number at most
\(N/H'\). Select a largest one, costing at most
\(n^{4\xi+o(1)}\). This proves the assertion with total
normalized loss at most \(6\xi+o(1)\) after the rare deletion.
One can, for example, use independent tolerances so that
the popularity cutoff has a larger exponent loss than
the original shell deficit. All count comparisons are made
first at fixed tolerances and fixed budget.
The original rows used as the triangles' legs still have
exponent \(2c_0-O(\xi)-o(1)\), by division by \(N^2\).
Their original radical, height and primitivity constraints
are retained simultaneously.

## Pinned transfer at the critical radical cap

After the preceding selection, choose one \(T\) and one
first complete value \(x_1\) retaining a relation of exponent
\(2c_0-O(\xi)-o(1)\) by averaging over at most \(N^2\) pins.
On the retained relation,
\[
          x_2=x_1+G D_{12},\qquad x_3=x_1+G D_{13}
\]
recovers both other endpoints from the primitive output.
Thus the pinned transfer is injective. Each of its three
output coordinate dictionaries has exponent at most
\(c_0+o(1)\), since each normalized side is the quotient of a
popular physical gap by this single \(G\).

If in addition the normalized radical on \(\mathcal Q\) is at
most \(n^{6c_0+o(1)}\), the output is a balanced small-budget
endpoint with all its complete dictionaries at the minimum
exponent. In this application discard short normalized heights
and small side radicals before the exact-cell and pin selections. Fibres of size below \(n^{c_0-\xi_0}\) contribute at
most \(N n^{3(c_0-\xi_0)}\) tagged triangles for any fixed
small \(\xi_0>0\). On the others the fixed-fibre rank bound
gives a uniform normalized height floor \(X^\theta\)
outside a power-smaller class. The fixed-fibre radical
estimate also gives all three normalized side radicals
at least \(n^{2c_0-o(1)}\), by taking the threshold for the
fibre exponent to \(c_0\) along fixed tolerances.
The total cap gives the matching individual upper caps.
On a dyadic output height shell, the absolute raw radical
exponent remains \(6c_0 L+o(L)\) in base \(X\); the relative
budget is this exponent divided by \(\log_X Y\ge\theta\).
Choose its fixed-budget mesh before the height limit. Thus
the one-third coefficient gives the matching total output
exponent \(2c_0\), and a fixed output letter has fibre
exponent at most \(c_0\) by (FL1). The complete dictionaries
have exponent exactly \(c_0\). This argument gives an
output with minimum dictionaries even if the output
nonatomic valuation mass is larger than its minimum
allowed value.
[Normalized triangle height and side bounds](ABC_FIXED_FIBRE_TRIANGLE_RADICALS.md#primitive-triangles-of-three-values-on-the-fibre)
· [Fixed-letter caps (FL1)](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)

The pinned output also has a joint arithmetic constraint
which a count on the output alone would lose. Its first and
third coordinates are sent by the *same* affine map
\(z\mapsto x_1+Gz\) into the fixed original \(T\)-fibre.
For each of the two images the original complementary
letter is another affine function with the same fixed
\(T,x_1,G\). All four complete image values retain
their original caps. The output equation itself retains all
three of its primitive radical caps. The unresolved
forcing step is whether the original minimum-dictionary
families necessarily support a full-size tagged relation
at a bounded or subcritical primitive triangle radical.

## A single fibre at an arbitrary positive exponent

The exact-gcd argument works in a single fixed complete-letter
fibre with \(n^{r+o(1)}\) distinct complete \(W\)-values for any
fixed \(r>0\). On an arbitrary triangle relation of exponent
\(3r-o(1)\) there is a subclass of full limiting exponent
with one exact gcd \(G\) and one \(W\)-residue modulo \(G\).
In the cell argument take \(N=n^{r+o(1)}\) and
\(H=n^{r-\xi}\). There is no sum over \(T\).
Thin cells contribute at most \(N^2 H X^{o_X(1)}\)
triangles. Thick cells number at most
\(N^2 X^{o_X(1)}/H^2\). These two bounds give the
claimed selection, with an exponent loss at most
\(2\xi+o(1)\) and with \(\xi\) arbitrarily small.

Pin the first value \(x_1\) from that cell. The normalized
output has \(n^{2r-o(1)}\) distinct triples. Both sides
incident to the pin have dictionaries of size at most
\(n^{r+o(1)}\), since the corresponding original endpoints
vary in the fixed \(W\)-dictionary and \(G\) is fixed.
The target of the output is one of these two sides.
There is no assertion on the size of the remaining
dictionary from a single fibre alone.

If the primitive triangle radical is at most \(n^{6r+o(1)}\),
the output has the balanced one-third exponent at this
radical scale and reaches **maximum target multiplicity**
there. In detail, pre-discard the short normalized
triangles and the triangles with a side smaller in
radical exponent than \(2r\) by a chosen fixed margin,
using the fixed-fibre counts. Choose the cell-loss
tolerance below the gaps in these discards, or discard
before the cell selection. For the new relative budget
use a finite grid above \(6rL\log X/\log(Y/2)\) at
each fixed inner budget in the dyadic output shell
\(Y/2<D_{13}\le Y\). The total normalized exponent is
\(2r-o(1)\), each output radical exponent is \(2r+o(1)\),
and each output target fibre has size at most \(n^{r+o(1)}\)
by (FL1). The target dictionary has that same exponent
by pinned counting and the matching lower bound from
(FL1). To extract a full maximum-multiplicity cohort,
discard target fibres of size less than \(n^{r-\eta}\)
for fixed \(\eta>0\). They account for a power-smaller
count because there are at most \(n^{r+o(1)}\) targets.
The same holds for the target-radical fibres since
each fixed radical gives only subpower many complete
values, and the complete-target lower multiplicities
already attain \(n^{r-o(1)}\).

In particular the critical triangle relation in a complete
subfibre at any positive intermediate least-repair
multiplicity transfers to maximum target multiplicity
at the new radical scale. This does not use a cap on the
radical of the original physical differences or their
exact gcd. The powerful-output exclusion applies to
this pinned subfamily as in the unpinned triangle transfer.
[Fixed-fibre counts and high-powerful output saving](ABC_FIXED_FIBRE_TRIANGLE_RADICALS.md#transfer-of-the-high-powerful-saving)

## Consequences of the fixed original progression width

Write \(\gamma=\lim\log_X G\) along a height-first
subsequence of the exact tagged selection. Each original
role now lies in a fixed residue class modulo \(G\).
The common parameter in a star of actual rows can be
written \(x=x_0+G t\), with \(t\) in an interval of length
\(O(X/G)\). Applying the compact prime-packing bound with
this numerical width gives, on any original role,
\[
 \operatorname{supp}\nu_{W,\mathrm{na}}
                  \subseteq[0,6(1-\gamma)],\qquad
 \operatorname{supp}\nu_W\subseteq[0,\max\{6,12(1-\gamma)\}]
\]
for the finite scaled-degree measures in the original
full-exponent subfamily. The second statement does not
bound escaped numerical moment. For the nonatomic measure the cost overshoot can be
arbitrarily small. At scaled orders above six, an
individual token already wider than the query interval
can be used on its own at affordable cost. Removing
that case bounds the remaining overshoot in terms of
the same interval. Mixed packets from different roles
give the sharper cutoff and indivisible-atom alternative
in (MPK3)–(MPK4).
[Compact prime packing with a fixed progression](ABC_MINIMUM_DICTIONARY_COMPACT_PRIME_PACKING.md)
· [Mixed packing across roles](ABC_CROSS_ROLE_COMPACT_PACKING.md)
