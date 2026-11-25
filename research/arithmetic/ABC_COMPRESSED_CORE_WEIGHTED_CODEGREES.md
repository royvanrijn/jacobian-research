# Weighted codegrees with a compressed core on one edge

A weighted difference at two fixed columns can use a compressed
high-degree core on just one of the edges. The other target
retains its entire original valuation list. This gives a
per-column cofactor count, a reduced-ratio count when both cores
are specified, and a bound on multiplicatively consistent
rectangles. The estimates allow polynomial-height cofactors
and arbitrary nonzero column determinants.
[Exact-power weighted codegrees](ABC_SHIFT_DICTIONARY_WEIGHTED_CODEGREES.md#cofactor-ratio-dictionaries-for-two-fixed-columns)
· [Compressed term in a fixed-shift fibre](ABC_FIXED_SHIFT_STRICT_HALF_SAVING.md#a-compressed-term-in-an-arbitrary-fixed-shift-fibre)

## Setting and finite-width hypothesis

Use a simple bipartite graph with row label \(x\) and column
label \(j=(y_j,k_j)\), where all labels are positive integers
and each column pair is primitive. On an edge let
\[
    T_{xj}=Fxy_j+\sigma Gk_j>0,\qquad
    F,G>0,\quad \sigma\in\{1,-1\}.
\]
The fixed factors, labels and targets are at most \(X^K\)
for fixed \(K\). In column \(j\) the target has a chosen
factorization from original singleton coordinates,
\[
    T_{xj}=A_j\prod_i u_{xji}^{e_{ji}},
\]
with \(A_j\) fixed on the column. The positive integral roots
within a term are pairwise coprime. At each fixed small budget
\(\lambda\), the degree lists and the choices of positions and
representatives below come from a finite family fixed before
the height limit. Assume \(u_{xji}\le X^{\alpha_{ji}+o_X(1)}\)
and \(\sum_i\alpha_{ji}\le C_0\lambda\) in each column. The
dictionary of root position \(i\) on column \(j\) has size
at most \(X^{W_{ji}+o_X(1)}\), with \(W_{ji}\le\alpha_{ji}\).
The height and dictionary allowances are uniform over columns.
Fix a scaled cutoff \(a>0\). Use \(B\) as a uniform bound
on the sum of costs of positions of degree less than
\(a/\lambda\) over any two columns.

On each edge split the target factorization into a cofactor
\(W_{xj}\) and the exact product \(R_{xj}\) over
a subset of original positions, so \(T_{xj}=W_{xj}R_{xj}\).
The cofactor includes \(A_j\). An edge is eligible if we
can choose at most \(J\) disjoint groups of positions in
\(R_{xj}\), each with a positive integer representative
\(q_g\), satisfying
\[
  \sum_g\ \sum_{i\in g}|e_{ji}-q_g|\alpha_{ji}
      +\sum_{\substack{i\ \mathrm{in}\ R_{xj}\\
                        i\ \mathrm{outside\ all\ groups}}}
                  e_{ji}\alpha_{ji}
       \le \omega.                                                   \tag{WC1}
\]
Here \(J\) is fixed. Each representative can be chosen from a
different valuation window. The core is a product of at most
\(J\) powered roots of height \(X^{O(\lambda)}\) and a
rational correction of height at most \(X^{\omega+o_X(1)}\).
This follows by multiplying the original roots in each group
and moving the exponent differences into the correction.

For every fixed \(\eta>0\) there is an \(\omega>0\) and a small
budget threshold depending on \(K,C_0,a,J,\eta\), such that the
following counts hold for every fixed budget below the threshold.
Choose the width so that the fixed-width estimate for (SHF2) also
holds at twice this width. Put \(H=X^{B+\eta\lambda+\varepsilon}\).
Constants may depend on the fixed budget, finite choice family
and \(\varepsilon>0\). They are uniform in the fixed coefficients,
labels and cofactors. In a specific fibre the low-position
costs on its remaining positions can replace the uniform \(B\).

## One-sided records and two-sided ratios

Fix distinct columns \(j,j'\). Write \(f_j\) for the number of
cofactor values on eligible edges of column \(j\).

1. For any fixed cofactor on column \(j\), there are \(O(H)\)
   common rows of an eligible edge at \(j\) and an arbitrary
   edge at \(j'\). Thus there are \(O(H f_j)\) such common rows.
2. For any fixed reduced cofactor ratio, there are \(O(H)\)
   common rows with at least one eligible edge. Let
   \(V_{jj'}\) count the distinct ratios on common rows
   with both edges eligible. The eligible codegree is
   \(O(H\min(f_j,f_{j'},V_{jj'}))\).
3. There are \(O(H)\) common rows with an eligible edge at
   \(j\) and reduced ratio numerator at most \(X^\omega\).
   The denominator can be arbitrary, and there is no charge
   for the number of ratios in this branch.

Eliminate \(x\) to obtain
\[
 y_{j'}T_{xj}-y_jT_{xj'}
       =\sigma G(y_{j'}k_j-y_jk_{j'})=:D_{jj'}\ne0.                  \tag{WC2}
\]
Distinct primitive positive pairs have distinct slopes.
In the first case, fix \(W_{xj}\) so the left term is
\(y_{j'}W_{xj}R_{xj}\). The other term is the entire product
\(y_jT_{xj'}\) with fixed coefficient \(y_jA_{j'}\).

In the second case write \(W_{xj}/W_{xj'}=r/s\) in lowest
terms. With \(W_{xj}=g r\) and \(W_{xj'}=g s\), the
integer \(g\) divides \(|D_{jj'}|\). Fix \(g\) at divisor cost
and write
\[
     y_{j'}rR_{xj}-y_jsR_{xj'}=D_{jj'}/g.
\]
In both cases the remaining equation is a fixed nonzero
signed shift with fixed polynomial-height coefficients.

For either shift \(t\), fix the parts of all remaining
original roots supported on the primes of \(t\). At each
fixed budget this costs \(X^{o_X(1)}\), uniformly in
polynomial-height \(t\), by the fixed-support reduction
for (NS1). Fix also the actual gcd of the weighted complete
terms at divisor cost, since it divides \(t\). On dividing
by it, the complete terms satisfy the primitivity condition
in (SHF2). The stripped root positions retain their original
degrees and their raw caps and dictionary costs can only
decrease. After dyadic refinement, (WC1) still holds with
\(\omega+o_X(1)\), since every summand is a nonnegative
multiple of a raw cap. The supported parts and the gcd
division enter the two fixed integer coefficients.
[Nonprimitive shift normalization](ABC_NONPRIMITIVE_FIXED_SHIFT_ENERGY.md#two-dictionary-bounds-for-any-nonzero-shift)

Apply (SHF2) to the resulting equation. Its vanishing-error
version implies the fixed-width formulation by contradiction.
If no positive width tolerance and positive budget threshold
worked for a fixed normalized margin, choose width and
budget tending to zero along violating families. At each
budget first take arbitrarily large heights with the degree
and grouping lists fixed and the inner exponent allowance.
The compression then has vanishing numerical correction,
contradicting (SHF2) at the fixed normalized margin.
Use a smaller normalized margin and inner allowance for
(SHF2) to absorb support, divisor and dyadic enumerations.
Sum over the finite choice family at each fixed budget.
A recovered target on either column specifies \(x\)
from its edge equation.

For the third bound apply the same support stripping directly
to (WC2), with shift \(D_{jj'}\), leaving the full original
root lists on both sides. Since \(\gcd(W_{xj},W_{xj'})\)
divides \(|D_{jj'}|\), the part of \(W_{xj}\) coprime to
\(D_{jj'}\) divides the reduced numerator \(r\). In particular
the product of stripped original roots with their powers
from the complement of the core at \(j\) is at most
\(X^\omega\). After dyadic refinement this gives the same
bound on the complement's numerical weight up to \(o_X(1)\).
Take all these original positions into the ungrouped
correction in (WC1). The complete variable part on \(j\)
then satisfies the compression criterion at width
\(2\omega+o_X(1)\). Apply (SHF2) to the two variable parts.

## Graph consequences

Let \(A\) be the number of rows and \(M\) the number of
eligible edges. Cauchy–Schwarz on the eligible row degrees
gives
\[
 \frac{M^2}{A}-M
     \ll H \sum_{j\ne j'}\min(f_j,f_{j'},V_{jj'}).
                                                                    \tag{WC3}
\]
With \(n=X^\lambda\), if the row and column caps have
normalized exponents \(p,b\) and every relevant \(f_j\) or
every relevant \(V_{jj'}\) has exponent at most \(v\), the
eligible edge exponent is at most
\[
                 \max\bigl\{p,\ b+(p+v+B/\lambda)/2\bigr\}
\]
up to the normalized allowances. In the complete-label
minimum-mass setting, \(p=b=c=1/6\). For \(B/\lambda\)
sufficiently small, a one-sided local cofactor deficit or a
two-sided ratio deficit below \(c\) gives a strict saving on
eligible edges. If the cofactors are multiplicatively
separable, \(W_{xj}=P_xQ_j\) with positive rational \(P_x,Q_j\),
one has \(V_{jj'}\le1\) even if the marginal cofactor sets
are large. The same inequality applies for arbitrary part
caps whenever they give a strict margin.
[Complete-label reduction](ABC_DIFFUSE_FACTOR_GRAPH.md#the-complete-label-reduction-close-to-minimum-mass)

Count ordered rectangles with distinct rows and columns,
at least one eligible edge, and
\[
                   W_{xj}W_{x'j'}=W_{xj'}W_{x'j}.                  \tag{WC4}
\]
Their number is
\[
                 \ll H\sum_x d(x)(d(x)-1),                         \tag{WC5}
\]
where \(d(x)\) is the degree of row \(x\) in the full graph.
Orient the rectangle in one of at most four ways so an eligible
edge is in row \(x'\). For fixed \(x,j,j'\) the identity fixes
that row's ratio, and the ratio bound needs just one eligible
edge there. Sum over the first-row stars. If there are at least
\(n^{2c-o(1)}\) eligible edges and both parts have caps
\(n^{c+o(1)}\), these rectangles number at most
\(n^{3c+B/\lambda+\eta+\varepsilon/\lambda+o(1)}\).
The total rectangle count has lower exponent \(4c\)
in the same limits. Indeed the sum of the ordered
column-pair codegrees is at least \(M^2/A-M\), and
Cauchy–Schwarz on those codegrees then counts at least
\(n^{4c-o(1)}\) ordered rectangles. For \(B/\lambda<c\)
with a fixed margin and sufficiently small allowances,
multiplicatively consistent rectangles
form an exponentially small fraction on the normalized
logarithmic scale.

A cofactor dictionary cardinality cap at the unsaved scale
\(n^c\) cannot imply a ratio deficit by itself. On the
complete graph with both parts indexed by \(\mathbb Z/m\mathbb Z\),
let \(W_{xj}\) be the distinct prime labelled by \(x+j\).
Each column has \(m\) cofactors, and every distinct column
pair has exactly \(m\) ratios by unique factorization.
This is a finite dictionary obstruction. A global use of
(WC3) still requires arithmetic forcing a cofactor or
ratio deficit, or a different use of (WC5). The width tolerance
is chosen after fixing the number of groups. Extending the
criterion by letting the number of groups grow would require
additional estimates on this dependence.
