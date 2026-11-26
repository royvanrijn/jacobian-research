# Two-sided cancellation in auxiliary-prime tree comparisons

The two diagrams in the finite-field comparison can
both be projected onto their cycles. This gives a
product of the two small projected norms. For sparse
masks, a further gain is available when a cycle crosses
the bottom sibling pairs, visits more than one bottom
quartet, or when there are two separate sibling cycles.
The refinements apply to the uniform subgroup law of
the original comparison.

## Hypotheses and the two-sided projection

Use the tree recursion and zero convention
`eq:src36` of the inverse Goldbach manuscript,
summarized in (AT1). The depth \(\ell\ge2\) is fixed,
\(r=2^\ell\) is the number of leaves in each diagram,
and the field is \(\mathbb F_q\) for an odd prime \(q\).
In both trees partition the leaves into \(b\) nonempty
blocks, pair the blocks across the trees, and use the
uniform law where products in paired blocks agree.
Write \(\mathcal T\) for the vector of their common
products. Conditioned on \(\mathcal T\), the two trees'
leaves are independent, and each tree separately
has independent uniform unit leaves before conditioning.
All unit parameters and the leaf conjugation patterns
are as in the source.
[Tree comparison and the subgroup law, Section 6](https://github.com/openai/math/blob/main/preprints/the-additive-indecomposability-of-the-primes-September-24-2026/build/sections/06-tree.tex)
· [Tree recursion and auxiliary law (AT1–AT2)](ABC_SIGNED_AUXILIARY_TREE_COMPARISON.md#the-tree-estimate-and-the-law-of-shared-slots)

Form the bipartite multigraph in each diagram whose
vertices are its \(b\) blocks and \(r/4\) bottom
quartets, with one edge for each leaf. Write
\(\epsilon_j\) for the mixed-character bound of the
masks used in diagram \(j\).
Each mask vanishes at zero and has probability
\(L^2\)-norm at most one. The formulas first hold
with a single mask per diagram and opposite
conjugations within each bottom pair. The extension
to multiple masks is given below.

**CT1.** If both graphs contain cycles, then
\[
 \left|\mathbf E W_1\overline{W_2}\right|
   \ll_\ell
       (\epsilon_1+q^{-1/4})(\epsilon_2+q^{-1/4}).
                                                               \tag{CT1}
\]
In particular \(b\le3r/4\) suffices. Indeed there
are \(r\) edges and \(b+r/4\le r\) nonisolated
vertices in each graph.

For each tree, the alternating cycle action preserves
both the block products and the quartet products.
Let \(\mathcal P_j\) be averaging over the action in
diagram \(j\). For either diagram, its projection-norm
proof in the source gives
\(\|\mathcal P_j W_j\|_2\ll_\ell\epsilon_j+q^{-1/4}\).
The norm is in the unconditional marginal of its
own uniform leaves. The action preserves the
conditional law given \(\mathcal T\), so
\(\mathbf E(W_j\mid\mathcal T)
=\mathbf E(\mathcal P_jW_j\mid\mathcal T)\).
Conditioned on \(\mathcal T\) the two values
are independent. Apply conditional-expectation
contraction and Cauchy--Schwarz to obtain
\[
 |\mathbf E W_1\overline{W_2}|
 \le \|\mathbf E(W_1\mid\mathcal T)\|_2
      \|\mathbf E(W_2\mid\mathcal T)\|_2
 \le \|\mathcal P_1 W_1\|_2 \|\mathcal P_2 W_2\|_2 .
\]
This uses the source's projection estimate in both
marginals. No additional independence of reconstructed
arguments is imposed.
[Source formulas `eq:tree-projection-reduction`, `eq:tree-cycle-relation`, and the completion of the cycle estimate](https://github.com/openai/math/blob/main/preprints/the-additive-indecomposability-of-the-primes-September-24-2026/build/sections/06-tree.tex)

## A stronger projection for specific cycle shapes

Separate each diagram's mixed-character bounds into
\[
 p=\max_a|\widehat g(a)|,\qquad
 \nu=\max_{\chi\ne1,\ a}|\widehat{g\chi}(a)|,\qquad
 \epsilon=\max(p,\nu).
\]
Take maxima over all its masks if there are several.
The Fourier norm is the probability normalization
from (AT1). Write \(\chi=1\) for the principal
character extended by zero at zero, which does not
change \(\widehat g\) since \(g(0)=0\).

**CT2.** In one diagram the projected norm may
be bounded by \(O_\ell(p^2+\nu+q^{-1/4})\) if
its graph has any of the following:

- A simple cycle with at least two quartets.
- A parallel-edge cycle whose two leaves lie in
  different bottom pairs of the same quartet.
- Two distinct parallel-edge cycles, each with its
  leaves in a single bottom pair.

Thus it suffices that its cycle space have dimension
at least two. In particular both diagrams satisfy
this condition if \(b<3r/4\). In (CT1) the improved
factor may be used independently on either diagram
with the requisite graph. The qualifying projection
may average over two commuting cycle actions.

### A quartet with moving leaves in different bottom pairs

In the source's bottom-pair autocorrelation notation,
let \(w_\chi(a)=\widehat{P_\chi}(a)\) and
\(r_\chi=\sum_a w_\chi(a)^2\). The formula for \(w_\chi\)
is \(q/(q-1)\) times the squared Fourier modulus of
\(g_*\overline\chi\) at a signed frequency, where
\(g_*\) is the pair's mask or its conjugate. Its mass
is \(O(1)\) and \(\sum_\chi r_\chi\ll1\) by
`eq:tree-total-r`. Therefore
\[
 r_1\ll p^2,\qquad
 \max_{\chi\ne1}r_\chi\ll\nu^2,\qquad
 \sum_\chi r_\chi^2\ll p^4+\nu^2.
\]
The last inequality uses the single principal square
and the bounded sum of the other \(r_\chi\).
[Bottom-pair autocorrelations and `eq:tree-total-r`](https://github.com/openai/math/blob/main/preprints/the-additive-indecomposability-of-the-primes-September-24-2026/build/sections/06-tree.tex)

In the two-bad case of
`eq:tree-twisted-convolution`, the
nonprincipal fourth moment before summing \(\chi\)
is at most \(2r_\chi^2+O(q^{-2})\). Thus replace
the source's bound \(\ll\epsilon^2+q^{-1}\) in
`eq:tree-all-nonprincipal-fourth` by
\(\ll p^4+\nu^2+q^{-1}\). The subsequent Cauchy--Schwarz
step allows the two bottom pairs to use different
arrays \(w_\chi\). If exactly one twist is principal
its contribution is already \(O(q^{-1})\), and if
both are principal it is \(O(\epsilon^4+q^{-2})\).
The friendly cases already cost
\(O(\epsilon^4+q^{-1/2})\) by `eq:src42`.
Since \(\epsilon\le1\), the source's majorant
for a cross-pair quartet consequently satisfies
\[
            \sup_\rho b_\rho
                     \ll p^4+\nu^2+q^{-1/2}.
\]
Here \(b_\rho\) is the \(y\)-average of the quartet
coefficient majorant, as in the completion of the
cycle estimate. All the bounded-sum estimates in
that proof persist. In the cycle constraint choose
this quartet as the one whose character is
determined by the others.
[Equations `eq:src42`, `eq:tree-twisted-convolution`, `eq:tree-fourth-moment`, and `eq:tree-principal-correction`](https://github.com/openai/math/blob/main/preprints/the-additive-indecomposability-of-the-primes-September-24-2026/build/sections/06-tree.tex)

### Cycles using only same-pair movements

For a same-pair quartet,
`eq:tree-same-pair-majorant` involves
\(\sum_{\chi^2=\rho^{\pm1}} |P_\chi(d)|^2\).
Its means satisfy \(b_1\ll p^2+\nu^2\) and
\(\sup_{\rho\ne1}b_\rho\ll\nu^2\): the principal
\(\chi\) can occur only when \(\rho=1\).
If a cycle visits at least two quartets and all
its movements are within bottom pairs, split the
final cycle sum into the all-principal contribution
and its complement. The first is at most
\(O_\ell((p^2+\nu^2)^2)\). In the complement some
quartet has a nonprincipal index. Sum over which
quartet this is and take it as the index determined
by the cycle relation. This gives \(O_\ell(\nu^2)\).
Taking the square root yields the stated norm bound.

The only cycle untreated by the preceding
arguments is a loop of two parallel edges
belonging to one bottom pair. Two distinct
such cycles use disjoint bottom pairs. Their
actions commute and each preserves the
conditioning by \(\mathcal T\). Project both,
requiring the principal index in each.
If the pairs are in different quartets,
use the two bounds \(b_1\ll p^2+\nu^2\) and
the horizontal density estimate of the source.
If they are in the same quartet, condition on
its product and expose the bottom pair products.
The relative argument \(d/e\), with \(d-e=y\),
is uniform on a square coset on valid points
up to exclusion of the invalid values, as in
`eq:tree-same-pair-majorant`.
Within the two pairs the splits are independent.
Projection in each pair gives a squared mean at
most \(2\sum_{\chi^2=1}|P_\chi(d)|^2\), or its
right-pair counterpart at \(e\), by
`eq:tree-square-pullback`. Average their
product over \(d-e=y\) and then over \(y\).
Dropping the restriction \(d\ne e\) bounds this
by the product of the two means, hence by
\(O((p^2+\nu^2)^2)\). The horizontal density
bound supplies the same estimate in the full
diagram. This proves (CT2) in the last case.
Any graph with at least two independent cycles
has either a cycle of one of the first two types
or two distinct parallel sibling cycles, since
its cycle space is generated by simple cycles.
[Same-pair majorant, square pullback, and the horizontal density estimate](https://github.com/openai/math/blob/main/preprints/the-additive-indecomposability-of-the-primes-September-24-2026/build/sections/06-tree.tex)

## Multiple masks and occurrence scales

These estimates allow a different mask for each
leaf. First use a common mask within each bottom
pair, allowed to vary by pair. All the preceding
bounds on \(P_\chi,w_\chi,r_\chi\) apply separately
to each pair. Conditional second moments and
the horizontal density bound multiply their
nonnegative factors. The friendly estimate
`eq:src42` uses only its own pair's
mask, and the two-bad Cauchy--Schwarz step uses
the separate fourth-moment estimates.
For two different masks \(f,g\) in one pair,
polarize with \(h_\omega=(f+\omega g)/2\),
\(\omega\in\{1,i,-1,-i\}\). The identity
\[
 f(x)\overline{g(y)}
       =\sum_{\omega}\omega h_\omega(x)\overline{h_\omega(y)}
\]
reduces to common-pair masks. They satisfy the
same separate bounds \(p,\nu\) by the triangle
inequality, have \(L^2\)-norm at most one,
and vanish at zero. Conjugation can absorb
any prescribed leaf sign. Expanding all
pairs increases the constant only in terms
of fixed \(\ell\).
[Different-mask polarization and the centered masks (AT3)](ABC_SIGNED_AUXILIARY_TREE_COMPARISON.md#different-masks-at-different-leaves)

For centered transforms of residue sets of sizes
\(s_i=n^{d+o(1)}\) at a prime \(q\asymp n^\beta\)
with \(\beta>2d\), (AT3) gives
\[
       p\le n^{-d/2+o(1)},\qquad
       \nu\le n^{-(\beta-d)/2+o(1)}.
\]
Thus (CT1) gives paired cancellation
\(n^{-d+o(1)}\) throughout its cycle range.
A diagram satisfying (CT2) instead contributes
the factor \(n^{-\min(d,\beta/4)+o(1)}\)
in the product. In particular two such
diagrams give \(n^{-2d+o(1)}\) for
\(\beta\ge4d\). These are fixed-depth statements,
and the test primes grow at fixed small budget
before the outer limit. They can use separate
masks for the three original roles.

The paired arithmetic histories would have to
produce this subgroup law or an approximation
whose signed error is below the required power
scale. They also need a positive statistic
forced by true primitive equations with the
three original radical caps and their repair
cutoffs in place. The finite-field bounds do
not supply that realization or lower bound.
[Arithmetic transfer obligations for the paired histories](ABC_SIGNED_AUXILIARY_TREE_COMPARISON.md#what-the-positive-input-would-need-to-retain)

## Direct row splits do not realize the uniform law

**CT3.** Fix a collection of \(M\) actual unit additive
rows at the test prime. Suppose a proposed history
realizes every internal split of a diagram by using
one of those rows, specifically by putting
\[
       y_+=y_j B/C,\qquad y_-=-y_j A/C.
\]
This is a valid additive split \(y_+-y_-=y_j\), but its
child ratio is in a set of at most \(M\) elements.
Under the uniform single-diagram leaf law, the event
that every one of the \(r-1\) internal splits meets
this ratio restriction has probability at most
\[
                   \left(\frac{2M}{q-1}\right)^{r-1}.
\]
The validity convention only lowers this probability.
In particular, when \(q\gg M^2\) by a fixed power,
any probability law supported on these histories has
total variation distance tending to one from the
uniform law in (CT1). The same conclusion holds
if every original row is also allowed arbitrary
common unit multiples.

To prove this, expose the product of all leaves
and then expose the two child products at each node
from the root down. With the parent product fixed,
the first child product is uniform on the unit group
and the second is determined. The ancestor data and
the node's current entries are already fixed. On
valid points the recursion gives \(y_+/y_-\) equal
to a fixed unit times the quotient of the two
child products. Since their product is fixed,
the ratio is a fixed unit times the square or inverse
square of the first child product. A specified ratio
therefore has conditional probability at most
\(2/(q-1)\). Multiply the successive conditional
bounds. The paired subgroup law has this uniform
single-diagram marginal. The prescribed-support
event has probability one for the proposed histories,
and common row multiplication does not change
\(-B/A\). This proves the claim.
[Arithmetic conditions on the paired history law](ABC_SIGNED_AUXILIARY_TREE_COMPARISON.md#what-the-positive-input-would-need-to-retain)
