# Codimension two for wide ordinary factorial profiles

Fix \(r\) unequal exponent rows with positive integral gaps \(d_i\)
and fixed orientations. Let \(A_i\) be their independent smaller
exponents. Equal rows may be appended without affecting any common
zeros. Write \(J\) for any absolute terminal cutoff from the
[effective eventual theorem](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL.md).

**Theorem.** For every terminal \(N\ge4\) with \(r\ge 2N\), the
ordinary nonnegative integer base vectors at which orders \(N-1,N\)
have a common zero are contained in a complex algebraic set of
dimension at most \(r-2\). If \(r\ge2(J-1)\), the same conclusion
holds for the set of base vectors failing SFC(2) at any order.

Consequently, for these fixed signed gaps the all-order exceptional
fraction on \(\{0,\ldots,H\}^r\) is \(O((H+1)^{-2})\), with a
constant depending on the fixed signed gaps and \(J\). From any
ordinary starting vector there is a strictly positive integral
direction along which all positive integer translations are
regular at every order. The initial vector itself is unrestricted.
These strengthen the [degree and grid bounds](MULTIVARIATE_ALL_ORDER_EXCEPTIONAL_LOCUS.md)
in sufficiently many unequal rows.

## Differential of the paired coefficient map

Put \(n=N-1\), \(q=n+N=2N-1\), and use the normalized factorial
polynomials \(F_j(A,z)\) in the
[resultant degree proof](MULTIVARIATE_TRANSLATED_RESULTANT_DEGREE.md).
Let \(a_{j,k}=[z^k]F_j(A,z)/F_j(A,0)\). On the complex open
set where all sampled factors \(jA_i+h\) are nonzero, these
coefficient ratios are nonzero. The \(q\) consecutive ratios
\[
 y_{j,k}=\frac{a_{j,k}/\binom jk}{a_{j,k-1}/\binom j{k-1}},
       \qquad j=n,N,\quad 1\le k\le j,
\]
give equivalent coordinates for the coefficient pair. All ordinary
nonnegative bases lie in this open set. For an increasing row of
gap \(d\), its column in the logarithmic differential has entries
\[
                  g_{j,k,d}(A)
              =\sum_{h=(k-1)d+1}^{kd}\frac1{A+h/j}.
                                                               \tag{1}
\]
A decreasing row reverses the \(k\) indices separately for each
order and changes their signs. The rank of the logarithmic
differential equals that of the ordinary differential because all
\(y_{j,k}\) are nonzero.

For \(d\ge2\), the \(q\) rational functions in (1) are linearly
independent. Each block contains a noninteger \(h/j\), since
\(j\ge2\) and the block contains consecutive integers. This
gives a pole unique to that block. Blocks in different orders
can only share poles at integers, by \(\gcd(n,N)=1\).
For \(d=1\) there is exactly one repeated pole across the two
orders, at \(h/j=1\). Thus the common relation space for an
increasing unit-gap column is the line spanned by
\(e_{n,n}-e_{N,N}\), where \(e_{j,k}\) denotes the coordinate
vector at \((j,k)\). For a decreasing unit-gap column the
line is spanned by \(e_{n,1}-e_{N,1}\).

Let \(s_+,s_-\) count the increasing and decreasing unit-gap rows.
Rank failure gives a projective vector \([c]\in\mathbb P^{q-1}\)
annihilating every column at the chosen bases. Away from the two
relation points identified above, each column equation is a
nonzero rational function in its own base. Hence each base has
only finitely many choices for fixed \([c]\). The incidence there
has dimension at most \(q-1\). At either of the two relation
points, only the corresponding unit rows are free, so the
dimension is at most \(s_+\) or \(s_-\). The closure of the
rank-defective set therefore has dimension at most
\[
                         \max(q-1,s_+,s_-).
\]
If \(r\ge q+1\) and \(s_+,s_-\le r-2\), this is at most \(r-2\).

## Common conjugate pairs in the full-rank range

For ordinary real bases, the even moment is strictly positive on
the real line by its real integral of an even power. Any common
root is therefore nonreal and its conjugate is a distinct common
root. Both have nonzero modulus.

In the complex coefficient chart take normalized polynomials
with constant term one. The pairs with two distinct nonzero
common roots form a constructible set of dimension at most
\(q-2\). Indeed, for each lower polynomial there are finitely many
choices of its two distinct nonzero roots \(\zeta_1,\zeta_2\).
Their vanishing imposes two independent linear conditions on the
upper coefficients. In particular the evaluation matrix on the
variables at \(z,z^2\) has determinant
\(\zeta_1\zeta_2(\zeta_2-\zeta_1)\ne0\). Restricting to the open
coefficient chart leaves the dimension bound intact.

At points where the paired coefficient map has full rank, it is
locally a complex submersion. Its inverse image of the
constructible coefficient locus therefore has dimension at most
\(r-2\). This follows in local holomorphic coordinates from the
product of the target locus and the smooth fibre of dimension
\(r-q\). The constructible inverse image has the same dimension
as its Zariski closure. Together with the critical-locus bound
this proves the claim when \(s_\pm\le r-2\).

## The case of a single dominant unit orientation

Suppose at least \(r-1\) rows are unit gaps in one orientation.
By reversing both exponent columns if necessary, assume it is
the increasing orientation. Keep the bases of the other rows as
separate coordinates, denoted by \(B\), and write \(s=s_+\).
Since \(r\ge q+1\), we have \(s\ge q\).

For these \(s\) rows the \(q-1\) functions in (1) after deleting
the duplicate are reciprocal linear functions with distinct
poles. A nonzero linear combination has numerator of degree at
most \(q-2\) over their common denominator. Consequently the
differential in these bases has rank \(q-1\) whenever they take
at least \(q-1\) distinct values. The configurations with at most
\(q-2\) distinct values lie in a finite union of coordinate
diagonals and have dimension at most \((r-s)+(q-2)\le r-2\).

Each selected unit-gap row contributes \(N(A_i+1)\) to \(y_{N,N}\)
and \(n(A_i+1)\) to \(y_{n,n}\). Thus the image of the map
\((A,B)\mapsto(B,y)\) lies in the smooth chart
\[
                         y_{N,N}/y_{n,n}=C(B)
\]
for a nonzero rational function \(C(B)\) incorporating the other
rows and the constant \((N/n)^s\). On the rank-\((q-1)\) set the
map is a submersion to this chart, since the \(B\)-coordinates
also vary independently. For fixed \(B\) and fixed lower
coefficients the relation fixes only the ratio of the upper
coefficients at \(z^N,z^{N-1}\). Since \(N\ge4\), its upper
coefficients at \(z,z^2\) remain free. The two distinct common
root evaluations again have independent linear conditions in
those two coordinates. The same submersion argument gives the
dimension bound after adjoining the diagonal configurations.

## Integer points and regular directions

The [already closed low terminals](MULTIVARIATE_TERMINAL_TWELVE_ALL_PROFILES.md)
handle \(N=2,3\). The effective eventual theorem handles
\(N\ge J\). If \(r\ge2(J-1)\), the union of the finitely many
other exceptional algebraic sets just constructed has dimension
at most \(r-2\).

For an affine algebraic set of dimension \(d\), generic linear
Noether normalization gives a finite linear projection onto a
space of dimension \(d\), and the linear forms can be taken
with integer coefficients. One can choose them by requiring
that their projective center at infinity miss the projective
closure of the set at infinity. This is a nonempty open
condition on the forms, so rational forms suffice and can be
cleared of denominators. A finite projection has uniformly
bounded fibre cardinality, since the coordinate ring is a
finite module. The projection sends an integer box of side
\(H+1\) into an integer box of side \(O(H+1)\), so at most
\(O((H+1)^d)\) of its integer points lie in the algebraic set.
Applying this with \(d\le r-2\) proves the grid estimate.

Finally fix a starting vector \(A^0\) and let \(V\) denote the
exceptional algebraic set. The closure of the image of
\(V\times\mathbb A^1\) under \((x,t)\mapsto t(x-A^0)\) has
dimension at most \(r-1\). Strictly positive integer vectors
are Zariski dense, so choose one \(v\) outside this closure.
If \(A^0+kv\) were in \(V\) at an integer \(k\ge1\), its image
with \(t=1/k\) would be \(v\), a contradiction. Hence each
of these ordinary profiles is regular at every order.
