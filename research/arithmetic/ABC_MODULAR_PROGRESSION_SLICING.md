# Fixed-rank modular progressions on regularized complete letters

Use the same degree-regularized complete letter \(\mathcal U\) and
small-budget limit as the
[all-scale progression estimates](ABC_DIFFUSE_FACTOR_GRAPH.md#an-all-scale-progression-bound-for-complete-letters).
Write \(n=X^L\), \(c=1/6\), and let \(s_*\) be the larger of the two
opposite nonatomic support suprema. Fix \(R\geq1\), \(\kappa\geq0\),
and an interval of width at most \(X^h\). In \(\mathbb Z/M\mathbb Z\),
with \(M\geq X^\beta\) of polynomial height, let \(P\) be given by a
coefficient-box presentation of rank at most \(R\) and volume at most
\(n^\kappa\). Then the normalized exponent of the distinct values of
\(\mathcal U\) in the interval reducing into \(P\) is at most
\[
 \min\!\left(c,(1-1/R)\kappa+\frac{(h-\beta)_+}{s_*}\right)
\]
up to arbitrary extraction loss.

Indeed, if the actual box volume is \(V\), its largest side has
length at least \(V^{1/R}\). Fix the other coefficients. This leaves at
most \(V^{1-1/R}\) slices, each an ordinary modular progression of at
most \(X^{\kappa L}\) terms. For any fixed \(\zeta>0\) these slices
have at most \(X^\zeta\) terms for all sufficiently small \(L\).
Apply (AC2) uniformly, sum over the slices, then send \(\zeta\) to
zero. The count permits noninjective presentations and nonunit steps.
In the complete ambient height range \(h=1\), \(\beta\geq1-o(1)\)
gives upper exponent \(\min(c,(1-1/R)\kappa)\), by taking fixed
lower bounds for \(\beta\) tending to one. For \(R\geq2\), capturing
exponent \(v>0\) in this range consequently requires
\(\kappa\geq Rv/(R-1)\). At rank one the intersection has exponent
zero.

