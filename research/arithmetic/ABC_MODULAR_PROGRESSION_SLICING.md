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

## A primitive additive graph with uniform progression avoidance

Fix \(c=1/6\) and a sufficiently small \(L>0\), and write \(n=X^L\).
There are three finite sets of distinct positive integers \(A,B,C\)
within \([X/5,X]\), each of size \(n^{c+o_X(1)}\), with
\(\gg n^{2c+o_X(1)}\) primitive edges \(a+b=c'\), \(c'\in C\).
All three doubling constants are \(n^{o_X(1)}\).
Nevertheless each part meets every interval-progression pair having at
most \(X^{1-\sqrt L}\) progression lengths in only \(n^{o_X(1)}\)
values. The bounds below also handle fixed-rank modular progressions.
The construction prescribes the complete integers and their additive
graph. Their prime-exponent spectra are free.

### The coefficient boxes

Take \(R=\lfloor\sqrt{\log n}\rfloor\),
\(D=\lfloor n^{c/R}\rfloor\), and \(H=\lfloor X/(100RD)\rfloor\).
Choose the integers \(g_1,\ldots,g_R\) independently and uniformly in
\([1,H]\). The doubled coefficient box
\(\mathcal B=(\{0,\ldots,2D\})^R\) has size \(n^{c+o_X(1)}\),
\(H=X^{1-o_X(1)}\), and every image under
\(u\mapsto g\cdot u\) is in \([0,X/50]\).
For each pair of distinct coefficient vectors the probability of a
collision is at most \(1/H\), by fixing all \(g_i\) except one with a
nonzero difference coefficient. A union bound gives injectivity
on \(\mathcal B\) with probability \(1-o_X(1)\) since \(2cL<1\).

An affine subspace of dimension \(j\) meets \(\mathcal B\) in at most
\((2D+1)^j\) coefficient vectors, since projection onto some \(j\)
pivot coordinates is injective on the subspace. Thus any fibre of
\(t>2(2D+1)^{r-1}\) coefficients contains at least \((t/2)^{r+1}\)
ordered affinely independent \((r+1)\)-tuples. For such a tuple
\(u_0,\ldots,u_r\), the difference matrix \(u_i-u_0\) has rank \(r\).
Condition on every generator except the \(r\) indexed by a nonzero
square minor. Any list of \(r\) outputs
\(g\cdot(u_i-u_0)\) determines those \(r\) generators uniquely.
In particular if each output lies in a fixed dictionary of at most
\(X^{a+o_X(1)}\) integers, the conditional probability is at most
\(X^{-r(1-a)+o_X(1)}\). Here \(r\) is fixed as \(X\to\infty\).

### Ordinary progressions

Fix \(0\le\tau\le1\). If all the tuple images belong to one interval
\(I\) and one residue class modulo \(Q\), with
\(|I|/Q\le X^\tau\), each difference lies among at most
\(2X^\tau+1\) known multiples of \(Q\). The original interval
position and residue do not enter this test. Moduli \(Q>2X\) give
singletons. Summing the preceding tuple expectation over the
\(O(X)\) possible smaller moduli gives at most
\[
 X^{1+(r+1)Lc-r(1-\tau)+o_X(1)}
\]
independent tuples satisfying those difference tests, counted with
the modulus. In particular for
\(\tau=1-\sqrt L\) choose a fixed \(r\) with
\(r(\sqrt L-cL)>1+cL\). With probability \(1-o_X(1)\) there is
no such independent tuple for any modulus, so every fibre has
at most \((2D+1)^{r-1}=n^{o_X(1)}\) values.

More generally, for any fixed \(\delta>0\), let
\(v=\max(0,c-(1-\tau)/L)+\delta\).
A violating fibre with \(t\ge n^v\) gives at least
\((t/2)^{r+1}\) tuples for its modulus. Markov's inequality and
the displayed expectation bound its probability by
\[
 X^{1+(r+1)L(c-v)-r(1-\tau)+o_X(1)}
 \le X^{2-(r+1)L\delta+o_X(1)}.
\]
Taking sufficiently large fixed \(r\) shows that all fibres have
exponent at most \(v\), capped at \(c+o_X(1)\) by box size.
This may be made simultaneous over any finite grid of
\(\tau,\delta\), and extended to prescribed fixed
\(\tau\) by approximation and a diagonal selection.
In particular, a fixed \(\xi>0\) gives at most
\(n^{(c-\xi)_++\delta}\) values in any given class modulo any
\(Q\ge n^\xi\). By the divisor bound it follows uniformly
that, for any fixed polynomial-height integer \(P\), the number
of values with \(\operatorname{rad}\gcd(P,u)\ge n^\xi\)
is at most \(n^{(c-\xi)_++\delta+o_X(1)}\).
Choose a squarefree divisor of \(\operatorname{rad}(P)\)
of size at least \(n^\xi\) for every such value and sum over
its \(n^{o_X(1)}\) possible choices.

### Modular progressions of fixed rank

Fix a rank bound \(J\), a polynomial-height exponent \(K\),
and parameters \(h,\beta,\zeta\ge0\). Consider any interval of
width at most \(X^h\) and any modular progression modulo
\(X^\beta\le M\le X^K\) with at most \(J\) coefficient directions and
coefficient-box volume \(\prod N_i\le X^\zeta\). For each difference
between two eligible values, choosing coefficient differences
uses at most \(\prod(2N_i+1)\le3^J X^\zeta\) residues.
Each of those residues has at most \(1+2X^h/M\) representatives
among differences of absolute value at most \(X^h\). Thus one may
take \(a=\zeta+(h-\beta)_+\) in the tuple estimate.
The modulus and directions, reduced modulo \(M\), have at most
\(O(X^{(J+1)K})\) possibilities by a crude count.
The coefficient lengths have at most \(O_J(X^{J\zeta})\)
possibilities. The modular origin and interval location cancel
in the differences, so no union bound over either is needed.

The same Markov argument shows that all these intersections
have exponent at most
\[
 \min\!\left(c,\max(0,c-(1-a)/L)+\delta\right)+o_X(1).
\]
For \(a<1-cL\), take \(r\) large enough that the expected number
of tested independent tuples already tends to zero. This gives
\(n^{o_X(1)}\) values in every such intersection. In particular
it applies to fixed-rank modular progressions of
\(n^{O(1)}\) coefficient-box volume in the complete height
range for any fixed modulus-weight lower bound sufficiently
close to one when \(L\) is small.

### Making the graph primitive

Choose generators meeting the preceding tests for the desired
finite parameter lists. Let \(S\) be the image of the small
box \((\{0,\ldots,D-1\})^R\). Take independent translations
\(\lambda_A,\lambda_B\) uniformly from the integers in
\([X/4,X/3]\), put \(A=\lambda_A+S,\ B=\lambda_B+S\), and
\(C=\lambda_A+\lambda_B+(S+S)\). The doubled-box injectivity
gives the claimed size for each part and all pairs of
\(A\times B\) yield an edge. It also gives doubling
\(n^{o_X(1)}\) for \(A,B\). For \(C\) the same upper bound
uses its coefficient-box side length \(2D-1\), its known
size, and the trivial at-most \((4D-3)^R\) possible
coefficient sums, without needing injectivity on that
larger box. Translation preserves every progression and
overlap bound established uniformly in the location and the
residue, since the proof used only differences.

For any fixed pair of coefficient representatives and a
prime \(p\le X\), the probability that \(p\) divides both
shifted values is at most \((1/p+1/J_X)^2\), where
\(J_X\asymp X\) is the translation interval cardinality.
This uses independence of the shifts. The union bound over
primes gives probability at most
\(\sum_{j\ge2}j^{-2}+O(\log X/X)<1\) for a nonprimitive
pair, uniformly in the fixed coefficients. Averaging the
number of primitive edges therefore yields some shifts
with a fixed positive fraction of all edges primitive.
A corresponding triple is primitive since \(c'=a+b\).
All three parts remain within \([X/5,X]\), and pair
uniqueness is automatic for numerical additive triples.

### Quantitative inverse theorems at the endpoint scale

The additive-energy argument can lead to a set with doubling
\(K=n^{o(1)}\). Raghavan's Theorem 1.4 bounds the dimension of a
covering convex coset progression by
\(O(\log(16K)\log\log(16K)\log\log\log(16K))\) and the logarithmic
overhead in its cardinality and the number of translates by the same
quantity. The condition \(\log K=o(\log n)\) alone allows this overhead
to exceed \(\log n\), and the dimension can grow in the fixed-budget
height limit. The coefficient-box construction above meets the
endpoint edge exponent even with much stronger ordinary progression
bounds. [Raghavan, Theorem 1.4](https://arxiv.org/html/2512.11217v1)

Bloom's Theorem 7 gives a quantitative lower bound on the larger of
the sumset and product set when each value has few distinct prime
factors. An application forcing additive expansion would also need
a product-set restriction on the same selected complete values.
The complete-label graph supplies an additive relation and prime
factor records, and it permits a large product set.
[Bloom, Theorem 7](https://arxiv.org/html/2512.04931v1)
