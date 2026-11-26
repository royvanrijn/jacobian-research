# All-order closure for overlapping opposed factorial windows

Fix any ordinary two-monomial core. For integers \(D'\ge D>K\ge e\ge1\)
and \(1\le e'<D'\), append the four selectable row types
\[
 (M,M+D),\quad (M+e,M),\quad
 (M+K+D',M+K),\quad (M+K,M+K+e').
 \tag{1}
\]
Their counts \(x,y,z,w\) are independent nonnegative integers and some
moving row is present.

**Theorem.** The original factorial moments in adjacent positive
orders are coprime above a computable height \(M_0\), uniformly over
all such counts, in every degree and over every characteristic-zero
field. The lower and upper factorial ramps in the two height
windows overlap. In particular the fixed core may itself have
a consecutive collision.

In this family the first mark of the second window records a
contribution from the already active increasing \(D\)-row of the
first window. That contribution gives an exact constraint on the
counts. With all four counts positive the signed-gap bound then
bounds the counts outright. When the fourth type is absent, the
remaining infinite count rays have full divisor rank by a
strict-ratio argument.

## First marks and signed-gap bounds

Use an integer \(J\ge16\) at least the
[absolute eventual cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL.md).
Only pairs \(n,N=n+1\), \(15\le N<J\), remain by the
[terminal-fourteen theorem](MULTIVARIATE_TERMINAL_FOURTEEN_ALL_PROFILES.md).
At each such pair take the marks \(jM+1\) and \(j(M+K)+1\)
for \(j=n,N\). The normalized offsets of these four affine
forms are distinct. At each first mark, a companion-order
factor cannot be proportional, since \(\gcd(n,N)=1\) and
the mark offset within its order is nonzero modulo that order.
The [exact slope-cone theorem](MULTIVARIATE_EXACT_SLOPE_CONE_MULTIPLICITIES.md#all-but-one-mark-gives-an-exact-slope-constraint)
supplies isolated primes for at least three of the marks beyond
a computable height. In particular each of the two windows
supplies at least one isolated mark. At such a prime the
companion Newton polygon is flat, and any unique minimum
of the marked ordinate array excludes a shared slope.

At the first mark of the first window the ordinates are
\(y,x+y,\ldots,x+y,x\). No type in the second window has
reached that mark. If \(x\ne y\), the minimum is unique,
so a collision forces \(x=y\), also if one of these counts
is zero.

At the first mark of the second window, the short row of gap
\(e\) in the first window contributes nothing because \(e\le K\).
The increasing \(D\)-row contributes only at some positive indices
and does contribute at the last index because \(K<D\). The two
endpoints of the marked array are therefore \(z,x+w\).
If \(w>0\), all interior ordinates are at least \(z+w>z\).
If the endpoints differed, one endpoint would be the unique
minimum. A collision with \(w>0\) therefore forces \(z=x+w\).

Let \(\Delta_C\) be the core signed gap. The
[signed imbalance theorem](MULTIVARIATE_ALL_DEGREE_SIGNED_IMBALANCE.md)
bounds \(\Delta_C+(D-e)x-D'z+e'w\) at a collision. In the case
\(w>0\), this is \(\Delta_C-(D'-D+e)x-(D'-e')w\), so \(x\) and \(w\)
are bounded. The other two counts are then bounded as well.
For these finitely many count vectors the
[all-direction affine ray theorem](MULTIVARIATE_AFFINE_BASE_RAY_ALL_DIRECTIONS.md)
supplies computable height cutoffs, with no condition on the core.

If \(w=0\), possible count vectors have
\((D-e)x-D'z\) in a fixed finite set of integers. Put
\(g=\gcd(D',D-e)\). These count solutions lie on finitely many
exact shifted rays whose repeated direction for the first three
types is \(D'/g,D'/g,(D-e)/g\); the count \(y\) always equals \(x\).
The repeat direction has signed gap zero. Integer offsets can be
obtained by choosing the representative of \(x\bmod (D'/g)\).
The following ranks apply equally to each of those shifted rays.

## Full divisor rank on every remaining ray

Use the raw divisor coordinates of the
[overlapping-window rank theorem](MULTIVARIATE_OVERLAPPING_AND_SEPARATED_PAIR_DIVISOR_RANK.md#raw-counts-and-the-normalization-map).
In order \(j\), put \(S_l(k)=1_{\{k\ge l\}}\) on \(1\le k\le j\).
For a potential left kernel write \(A_l=\sum_{k\ge l}\alpha_k\).
The factor with normalized offset \(q=t/j\), \(j\nmid t\),
occurs only at its own order, so its linear constraint can be
read in this block of raw coordinates. The raw divisor of an
increasing row that starts sampling a factor at \(l\) is \(S_l\).
The raw divisor of its reflected row is \(-S_{j-l+1}\).

For \(D<q\le K+D'\), only the decreasing long row can contribute.
Its local offset is \(h=t-jK\). Each block
\((l-1)D'<h\le lD'\) with
\(l\ge \lfloor j(D-K)/D'\rfloor+1\) has a choice with \(q>D\)
and \(j\nmid h\). In the first such block take \(h=j(D-K)+1\);
in each later block use any offset not divisible by \(j\).
The annihilation equations give
\[
                 A_k=0\quad(1\le k\le j-\lfloor j(D-K)/D'\rfloor).
\]
Since \(D'\ge D\), this range includes \(\lceil jK/D\rceil\).
For \(l>j-\lfloor j(D-K)/D'\rfloor\), choose
\(t\in((l-1)D,lD]\) with \(j\nmid t\). This mark has
\(K<q\le D\), so both long rows contribute and the short
row contributes nothing since \(K\ge e\). Its constraint is
\[
            A_l=\frac{D-e}{D'}
                   A_{\,j-\lceil(t-jK)/D'\rceil+1}. 
\]
This defines one link for each index in the remaining range.
Iterating links reaches a previously zero index or a directed
cycle. In a cycle the multiplier around it has modulus less
than one. Thus all \(A_l\) vanish.

All chosen divisors are exclusive to their own order.
This proves full raw spanning in both orders, and the
[normalization map](MULTIVARIATE_OVERLAPPING_AND_SEPARATED_PAIR_DIVISOR_RANK.md#raw-counts-and-the-normalization-map)
is surjective onto the affine divisor coordinates. The
[effective Ax affine theorem](MULTIVARIATE_AX_AFFINE_COLLISION_BRANCHES.md)
therefore applies to every shifted count ray, including at
zero curvature strength. Its uniform height cutoffs and the
finite-count cutoffs prove the claimed theorem.

The [separated-window integer count theorem](MULTIVARIATE_EXACT_SLOPE_CONE_MULTIPLICITIES.md#all-degree-closure-for-two-opposed-windows)
uses disjoint windows to force two balance equations. Here
one balance equation is replaced by the contribution from the
overlapping long row.
