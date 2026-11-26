# Separation rank in randomized tests for additive occurrence

An affine hash reduces a multiplicity test to an occurrence test
at one randomly sampled bucket. Further collision gates can lower
the false-pair rate while keeping each summand separate from the
other after the target is fixed. The required number of nonnegative
separated terms has a sharp universal lower bound. Factoring all
three letters separately has a different lower bound. A signed
version measures the false energy through the actual distribution
of target multiplicities.

## A sharp count for two separated legs

Let \(S\) be any finite set of \(J\ge2\) distinct labels, and let
\(K_\omega:S^2\to[0,1]\) be a random acceptance kernel. Each outcome
is a sum of at most \(d\ge1\) products \(f_i(x)g_i(y)\) of
nonnegative functions. Put
\[
 p=\frac1J\sum_x\mathbb E K_\omega(x,x),\qquad
 e=\frac1{J(J-1)}\sum_{x\ne y}\mathbb E K_\omega(x,y).
\]
Then
\[
              e\ \ge\ \frac{Jp^2/d-p}{J-1}.                 \tag{SR1}
\]
Indeed \(f_i(x)g_i(x)\le1\). For one outcome, Cauchy–Schwarz gives
\[
 \left(\sum_x f_i(x)\right)\left(\sum_y g_i(y)\right)
      \ge\left(\sum_x\sqrt{f_i(x)g_i(x)}\right)^2
      \ge\left(\sum_x f_i(x)g_i(x)\right)^2.
\]
Sum over the \(d\) terms, apply Cauchy–Schwarz, then take expectations
and apply Jensen to the trace. Subtract the expected trace \(Jp\).

The same bound holds if the matrix has ordinary rank at most \(d\)
and entries in \([0,1]\). More generally, if \(K_\omega\) is signed,
has ordinary rank at most \(d\), and its diagonal belongs to
\([0,1]\), then the mean off-diagonal square is at least the
right side of (SR1). The nuclear norm bounds the absolute trace
and is at most \(\sqrt d\,\|K_\omega\|_{\rm F}\), by the singular
value decomposition. The sum of diagonal squares is at most the
trace. Apply Jensen as above. For entries in \([0,1]\), the sum of
squares is at most the sum of the entries.

For integer labels and \(p\in[0,1]\), the upper construction
uses two independent affine maps \(h_j(x)=\alpha_j x+t_j\pmod1\)
with independent Haar variables \(\alpha_j,t_j\in\mathbb R/\mathbb Z\).
Retain \(x,y\) only when both first hashes land in a fixed
interval of length \(p\) and the second hashes land in the same
member of a partition into \(d\) intervals of length \(1/d\).
This is a sum of \(d\) nonnegative products. True diagonal entries
are accepted with probability \(p\), and each off-diagonal entry
with probability \(p^2/d\). In fact for \(x\ne y\) the map to
\((h_j(x),h_j(y))\) has nonzero integer determinant and is onto the
two-dimensional torus, so its two coordinates have independent
Haar laws. This gives the scale of (SR1) when \(d\) is small compared
with \(Jp\).
[Strongly universal hashing, Wegman–Carter](https://doi.org/10.1016/0022-0000(81)90033-7)
· [Nonnegative factorization and randomized communication in expectation](https://arxiv.org/abs/1105.4127)

When testing \(a+b=c\) at fixed \(c\), the labels on the two
legs are \(x=a\) and \(y=c-b\). To retain each true pair with
probability about \(1/H\) and to give a uniform false rate at most
\(1/B_0\), (SR1) requires on a sufficiently large full dictionary
an order of at least \(B_0/H^2\) terms. More precisely an upper
bound \(e_0\) on each false acceptance forces
\[
                      d\ge
       \frac{p^2}{(1-1/J)e_0+p/J}.
\]
Thus the free-bucket factors in the affine hash implement the
two-leg tradeoff. This is a statement for a full equality grid.
An arithmetic estimate may instead bound the actual eligible
false pairs without giving a uniform bound on every cell.

## Three separated letters

Let \(A,B,C\) be finite dictionaries, with \(M=|A||B||C|\).
Fix a subset \(\mathcal T\subseteq A\times B\times C\) with
\(T\ge1\) members such that each two coordinates determine at
most one member. In particular \(\mathcal T\) can be the true
additive incidence set. Suppose \(K_\omega\in[0,1]\) is a sum
of at most \(d\) products \(f_i(a)g_i(b)h_i(c)\) of nonnegative
functions. If the average expected acceptance on \(\mathcal T\)
is \(p\), then
\[
   \sum_{\mathbf x\notin\mathcal T}\mathbb E K_\omega(\mathbf x)
                 \ \ge\ \frac{(pT)^{3/2}}{\sqrt d}-pT.       \tag{SR2}
\]
To prove this for one term, take its three maximum values
\(v,w,z\) and its sums \(F,G,H\). If the term is nonzero
then \(vwz\le1\) since its maximum on the full product is
at most one. Its sum \(S_i\) on \(\mathcal T\) satisfies
\(S_i\le \min(zFG,wFH,vGH)\). This follows by forgetting
one role and using the other pair's unique completion.
Multiply the three upper bounds to get \(S_i^3\le(FGH)^2\).
The full-cell mass \(FGH\) is therefore at least \(S_i^{3/2}\).
Sum over the terms and use convexity twice, first on the
terms and then on the random outcome.

On the addition table of a finite group of size \(J\), with
three copies of the group, \(T=J^2\) and \(M=J^3\). The false
acceptance rate is at least
\[
               \frac{p^{3/2}/\sqrt d-p/J}{1-1/J}.
\]
There are sharp examples in finite vector spaces. Let
\(G=\mathbb F_q^s\), choose a uniform random linear form
\(\ell:G\to\mathbb F_q\), and take independent uniform
shifts \(t_A,t_B\). Use the three residues
\[
        \ell(a)+t_A,\quad \ell(b)+t_B,\quad
        \ell(c)+t_A+t_B .
\]
Choose \(d\le q^2\) fixed ordered residue pairs \((r,s')\);
accept when the three residues equal \((r,s',r+s')\) for
one selected pair. The resulting disjoint rectangles give
at most \(d\) separated products. For each true triple the
probability is \(p=d/q^2\). For each false triple it is
\(d/q^3=p^{3/2}/\sqrt d\). The shifts make the first two
residues uniform independently of the discrepancy residue
\(\ell(a+b-c)\), which is uniform for a nonzero discrepancy.
Take the dimension \(s\) large if needed to make the finite
\(J\) correction negligible.

## Signed energy for actual target multiplicities

The rank constraint has a version directly on a genuine
integer incidence set. Use complete target values in \(C\)
and an arbitrary subset of eligible true pairs \((a,b)\)
at each \(c\). The sets of pairs for different targets are
disjoint by exact addition. Let the multiplicities, in
decreasing order, be \(m_{(1)},m_{(2)},\ldots\), with zeros
inserted as needed. Flatten a signed, possibly complex-valued
test as a matrix \(K_{c,(a,b)}\) on the full product of
the target and pair dictionaries. For every matrix of
rank at most \(d\), with
\(S_K=\sum_{a+b=c,\ {\rm eligible}}K_{c,(a,b)}\), one has
\[
        |S_K|^2\ \le\
      \left(\sum_{i=1}^d m_{(i)}\right)
      \sum_{c,(a,b)} |K_{c,(a,b)}|^2.                        \tag{SR3}
\]
Let \(J_0\) be the eligible true-incidence matrix, of
zeros and ones. Its row supports are disjoint, so
\(J_0J_0^*\) is the diagonal matrix of multiplicities.
Project \(J_0\) orthogonally onto the column space of \(K\).
The Hilbert–Schmidt inner product with \(K\) is unchanged.
The squared Frobenius norm of the projection of \(J_0\)
is at most the sum of the \(d\) largest diagonal entries
of \(J_0J_0^*\). In fact its trace formula is a weighted
sum of these entries with weights in \([0,1]\) and total
weight at most \(d\). Cauchy–Schwarz proves (SR3).
For random \(K\) the same inequality with
\(\lvert\mathbb E S_K\rvert^2\) on the left and the expected
energy on the right follows by Jensen.

This count can retain the three radical caps, the exact
primitive quotient and its least repaired height in defining
\(J_0\). The energy on the right of (SR3) is over the full
matrix on which the rank is asserted. A real or arithmetic
eligibility mask can alter that rank. In particular these
bounds do not infer large noise inside a short real window.
A saving still requires an arithmetic estimate on eligible
incidences, or a signed estimate that avoids complete
full-grid energy.
[Filtered signed popularity](ABC_CAPPED_SIGNED_POPULARITY.md#capped-comparisons-at-the-filtered-prime-scale)
