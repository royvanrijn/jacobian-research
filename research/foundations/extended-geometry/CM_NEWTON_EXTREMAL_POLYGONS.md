# Extremal split CM Newton polygons and sharp numerical p-rank

Let a CM action have \(e\) conjugate pairs and relative rank \(h\).
Orient each pair by the smaller Hodge multiplicity \(m_i\), with
\(0\le m_i\le H=\lfloor h/2\rfloor\). A component is a lower convex
polygon of length \(h\), total slope sum \(m_i\), slopes in \([0,1]\),
and integral slope-change vertices. Its conjugate has the complementary
slopes in reverse order. The results below are exact for these numerical
polygons and the simultaneous relative norm inequalities from the
[arithmetic manuscript](../../../unpublished/cm-newton-polygon-constraints/main.tex).
A polygon passing these tests is a numerical possibility. Its occurrence as a
reduction of a fixed abelian variety is a separate arithmetic question.

## Least polygon with prescribed endpoint slopes

Fix one component of total \(m\), with \(z\) slope-zero terms and \(o\)
slope-one terms. The case \(q=h-z-o=0\) is ordinary and requires
\((z,o)=(h-m,m)\). If \(q>0\), deleting the zero and one runs leaves a
polygon on \([0,q]\) ending at \((q,w)\), where \(w=m-o\) is an integer
strictly between \(0\) and \(q\). Put \(a=q-w\) and \(b=w\).

At any internal vertex \((u,v)\), the positive interior slopes and
integrality give \(v\ge1\). The remaining slopes are strictly less than
one, so \(w-v<q-u\), and integrality gives \(v\ge u-a+1\). If \(b=1\),
there are no internal vertices because \(v<b\). If \(a=1\), there are
none because \(v<u\) conflicts with \(v\ge u\). Thus either case has only
the chord with slope \(w/q\).

If \(a,b\ge2\), the broken line through \((0,0),(a,1),(q,w)\) is convex,
with slopes \(1/a\) and \(1-1/b\). It lies below the height-one line
to the left of \(a\) and below the line of height \(u-a+1\) to its right.
Consequently it lies below every internal vertex of any admissible
polygon. Each intervening segment of the latter lies above the proposed
line by convexity of the proposed line. This also covers an isoclinic
competitor, which has only the common endpoints. The proposed polygon
is admissible, since all its displayed vertices are integral, and has
the prescribed endpoints and no extra zeros or ones. When \(a=b=2\)
the two proposed segments coincide.

If \(S_j\) is the sum of the first \(j\) slopes in the component, the
conjugate prefix is \(S_j^*=j-m+S_{h-j}\). Conjugation therefore preserves
pointwise order. Replacing a polygon by the least one with its endpoint
counts can only decrease both rounded prefixes in the relative
coefficient inequality.

The order convention rests on integral vertices as in
[Oort, *Newton polygons and p-divisible groups*, Astérisque 298, p. 257](https://www.numdam.org/item/AST_2005__298__255_0.pdf#page=4).
For a symmetric height-\(2g\) polygon with \(z=o=f\), the formula specializes
to the extremal symmetric polygon in
[Achter–Pries, *Generic Newton polygons for curves of given p-rank*, Definition 2.4](https://arxiv.org/pdf/1311.5846#page=6).
Their order compares the same ordinates in the reverse direction.
Here the components can have different numbers of zeros and ones, and
the rounded constraints couple many conjugate pairs.

## Endpoint profiles and the ordinarity classification

For a nonordinary endpoint choice the strict inequalities
\(0\le z<h-m\) and \(0\le o<m\) are necessary and sufficient. For the
least polygon just constructed, its rounded prefix sum exceeds the
ordinary prefix \(\max(0,j-(h-m))\) by
\(\mathbf 1_{z<j<h-o}\) at each integer \(j\). Indeed, between \(z\) and
\(h-m\) a positive displacement is at most one, and between \(h-m\) and
\(h-o\) the displacement decreases from one to zero. If \(a=1\) or
\(b=1\), the chord has this same rounding property by direct substitution.
The formula extends to the ordinary choice with an empty interval.
Conjugating and exchanging \(z,o\) gives the pair's exact minimum
penalty over its ordinary rounded sum

\[
 \Delta_j(z,o)=\mathbf 1_{z<j<h-o}+\mathbf 1_{o<j<h-z}.
\]

The ordinary slack before subtracting one is
\(D_j=\sum_i\min(j,h-j,m_i)\), by adding the ordinary prefixes of each
component and conjugate. Hence an endpoint vector is numerically
admissible under all relative norm constraints precisely when
\(\sum_i\Delta_j(z_i,o_i)\le D_j-1\) for all \(1\le j<h\), if
some pair is mixed. Its total \(p\)-rank is \(\sum_i(z_i+o_i)\).

The minimal nonordinary penalty comes from \(z=h-m-1\) and
\(o=m-1\). This change inserts two half-slopes and gives
\(\Delta_j=\mathbf 1_{j=m}+\mathbf 1_{j=h-m}\). Any nonordinary choice
has both endpoint counts at most these counts, so its two intervals
contain the intervals of this change. There is consequently only one
pointwise smallest penalty vector to test for each changed pair.

With no mixed pairs there is no change to test. With a single
multiplicity-one mixed pair, the slack \(D_1=1\) cannot absorb even
one penalty. At relative rank two with exactly two mixed pairs,
\(D_1=2\) cannot absorb the double midpoint penalty. At rank four
with one balanced pair, \(D_2=2\) also cannot absorb that penalty.
For all other signatures a single change is admissible. If it is
made on an unbalanced pair of positive \(m\), the two affected
indices are distinct, and \(D_m=D_{h-m}\ge m+B-1\), where
\(B\) is the number of mixed pairs. This is at least two when
\(B\ge2\) or \(m\ge2\). If all mixed pairs are balanced, the
midpoint charge is two and \(D_m=Bm\). Outside the stated cases
\(Bm\ge3\) for this balanced case. All unaffected indices have
positive slack. These are explicit nonordinary numerical witnesses
in every remaining signature.

For \(eh\le30\), the four disjoint classes have sizes \(111,81,14,7\)
in the order given above. The first is the number of pairs of positive
integers \((h,e)\) with \(he\le30\), and the second omits \(h=1\)
from that count. The third counts \(h=2\) with \(2\le e\le15\),
and the fourth counts \(h=4\) with \(1\le e\le7\). The full
signature count is

\[
 \sum_{h=1}^{30}\;\sum_{e=1}^{\lfloor30/h\rfloor}
       \binom{e+\lfloor h/2\rfloor}{\lfloor h/2\rfloor}=1188 .
\]

## Sharp numerical minimum p-rank in closed form

Let \(B>0\) be the number of mixed pairs. The definite pairs contribute
precisely \(h(e-B)\) to the minimum. Define for \(1\le j\le H\)

\[
 L_j=\max\left(0,\,2B+1-\sum_{i:m_i>0}\min(j,m_i)\right).
\]

These are nonincreasing and \(L_1=B+1\). The minimum contribution of
the mixed pairs is as follows.

| Signature condition | Minimum mixed contribution |
| --- | --- |
| \(h=2\) | \(2\lceil(B+1)/2\rceil\) |
| All mixed \(m_i=1\), odd \(h\ge3\) | \(H(B+1)+1\) |
| All mixed \(m_i=1\), even \(h\ge4\) | \((H-1)(B+1)+\lceil(B+1)/2\rceil+1\) |
| \(h=4\), all mixed \(m_i=2\) | \(B+3\) |
| Some mixed \(m_i\ge2\), odd \(h\) | \(\sum_{j=1}^H L_j\) |
| Some mixed \(m_i\ge2\), even \(h\), excluding the preceding rank-four case | \(\sum_{j=1}^{H-1} L_j+\lceil L_H/2\rceil\) |

To prove these formulas let \(P_j=\sum_{i:m_i>0}(2-\Delta_j(z_i,o_i))\),
so that the inequalities ask for \(P_j\ge L_j\). For each nonordinary
pair its contribution for \(j\le H\) is exactly

\[
  \mathbf1_{z\ge j}+\mathbf1_{o\ge j}+\mathbf1_{j\ge h-z}.
\]

Here \(o\le m-1 < H\) and \(o+z<h\), so the tail \(j\ge h-z\)
does not overlap \(j\le o\). An ordinary pair has constant contribution
two. If \(h\) is odd, summing the contributions for a nonordinary
pair gives exactly its cost \(z+o\). An ordinary pair contributes
\(2H\), one less than its cost \(h\). Thus the mixed contribution
is at least \(\sum_j L_j\). If \(h\) is even, let \(R\) be the number of
nonordinary pairs with \(z\ge H\), and \(C\) the number of ordinary mixed
pairs. The sum of a nonordinary pair's contributions is its cost plus
one exactly for a pair counted by \(R\), and an ordinary pair contributes
its cost. The midpoint sum is \(P_H=2(R+C)\). Therefore the mixed cost
is at least

\[
 \sum_{j<H}L_j+P_H-R
 =\sum_{j<H}L_j+R+2C
 \ \ge\ \sum_{j<H}L_j+\lceil L_H/2\rceil .
\]

When some \(m_i\ge2\), the last two rows attain these lower bounds
using nonordinary pairs only. Set \(T_j=L_j\) for \(j<H\) and also
at \(j=H\) for odd \(h\). For even \(h\), set \(T_H=\lceil L_H/2\rceil\).
Choose the two endpoint lengths \(z_i,o_i\le H\) so that exactly
\(T_j\) lengths are at least \(j\). Such lengths exist if at least \(T_j\)
of the allowed upper bounds are at least \(j\). The required lengths
have \(T_j-T_{j+1}\) copies of \(j\), and arranging lengths and upper
bounds in decreasing order gives a feasible assignment. The bounds are
\(h-m_i-1,m_i-1\), with each truncated at \(H\). At the first layer,
all \(B\) zero slots and at least one unit slot are available.
For \(j\ge2\) the existence of \(m_i\ge2\) implies \(L_j\le B\).
All \(B\) zero slots are available if \(j<H\), and also if \(h\) is odd.
When \(h\) is even and \(H\ge3\), the zero slots at \(H\) exclude
the \(B_{\mathrm{bal}}\) balanced pairs. The sum of the excesses
\(\sum_i(m_i-1)\) over mixed pairs is at least \(B_{\mathrm{bal}}+1\).
Indeed this follows from some \(m_i\ge2\) if \(B_{\mathrm{bal}}=0\),
and the balanced pairs each have excess at least two otherwise.
It follows that \(L_H\le B-B_{\mathrm{bal}}\).
For \(h=4\) with both mixed types present, writing \(u\) for the
number of multiplicity-one pairs gives \(L_H=u+1\) and
\(\lceil L_H/2\rceil\le u\), the number of available zero slots.
The lengths then give \(P_j=T_j\) at \(j<H\) and also at the odd
midpoint, while an even midpoint has \(P_H=2T_H\).
The cost is \(\sum_jT_j\), as required.

At \(h=4\) with all mixed \(m_i=2\), a nonordinary pair has
\(z,o\le1\). The midpoint therefore requires at least one ordinary
pair. If \(C\) pairs are ordinary the layer-one constraint forces
cost at least \(B+1+2C\). With one ordinary pair, put one endpoint
on each of the other \(B-1\) pairs to attain \(B+3\).

When all mixed \(m_i=1\), we have \(L_j=B+1\) throughout. At
\(h=2\), a nonordinary pair has no endpoints, and
\(\lceil(B+1)/2\rceil\) ordinary pairs are necessary and sufficient.
For \(h\ge3\), every nonordinary pair has \(o=0\) and \(z\le h-2\).
At \(j=1\) at least one ordinary pair is necessary. For a fixed
number \(C\ge1\) of ordinary pairs the mixed nonordinary part must
supply \(T=\max(0,B+1-2C)\) at each layer. There are at least
\(T\) nonordinary pairs available. If \(h\) is odd, choose \(z=H\) on
\(T\) of them. If \(h\) is even, put \(z=H\) on \(\lceil T/2\rceil\)
and \(z=H-1\) on the other selected pairs. Use \(z=0\) on all others.
The parity lower bounds above are attained, giving costs \(hC+HT\)
or \(hC+(H-1)T+\lceil T/2\rceil\), respectively.
Each expression strictly increases with \(C\), because increasing
\(C\) by one decreases \(T\) by at most two. Setting \(C=1\) gives
the two stated formulas.

## Exact replay and scope

The [endpoint checker](verify_cm_newton_endpoint_bounds.py)
recomputes every one of the 1,188 minima in the
[original census](cm_newton_polygon_census.json), as well as every
ordinary-forcing classification. Its independent exhaustive endpoint
comparison runs through dimension 20. To generate all exact
minimum rows through dimension 100 or the corresponding aggregate,
run the following commands from this directory.

    python3 -I verify_cm_newton_endpoint_bounds.py --all-minima
    python3 -I verify_cm_newton_endpoint_bounds.py

The new minimum computation has 98,729 signatures through dimension
100, with 938 ordinary-forcing cases. The extended aggregate replay
with \(\texttt{--bound 200}\) covers 5,475,208 signatures, with 2,145
ordinary-forcing cases. The SHA-256 digests of the newline-terminated
exact minimum rows are
\(\texttt{6b82f5ce9c1933b129c3ef78b4db086cc14cb59189a7d31647b4ac6efbb47c19}\)
through dimension 100 and
\(\texttt{fcaf6e7f71ae8108da41b5ddd439625112bb40f39f1031328fed514d381dee53}\)
through dimension 200. The original count of all
surviving polygons cannot be recovered by retaining just the least
polygon for each endpoint tuple. For example, with a single pair
of rank \(8\), \(m=4\), and \((z,o)=(2,1)\), the extremal component
has slopes \(0,0,1/2,1/2,2/3,2/3,2/3,1\) and satisfies every
inequality with its dual. The component
\(0,0,3/5,3/5,3/5,3/5,3/5,1\) has the same endpoints, but at the
midpoint its top rounded-down sum is two and its bottom rounded-up
sum is two, so it fails. The endpoint method determines feasibility
and sharp minimum rank. Multiplicities require the full profiles.
