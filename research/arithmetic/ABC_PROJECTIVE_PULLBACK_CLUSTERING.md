# Projective clustering in fixed binary parameter families

A congruence on a projective slope pulls back to boundedly many
parameter classes at each prime power, at a depth controlled by
the degree of the two forms. The pullback remains valid at primes
of bad reduction and at ramified slopes. For the false Pythagorean
carriers that pass the first two support-prime layers, it gives
a bound on every joint real and projective cluster. An exponent-full
family of those carriers has no common projective modulus of
fixed positive numerical weight.

[Divisor lattice and the gcd of the minors (DL4–DL5)](ABC_DIVISOR_LATTICE_LINE_OCCUPANCY.md#dl4-the-gcd-of-the-actual-projective-minors)
· [Pythagorean carriers and auxiliary unit tests (AU1)](ABC_OCCURRENCE_AUXILIARY_UNIT_TESTS.md#au1-uniform-supply-of-false-candidates)

## A projective pullback bound at every prime power

**PP1.** Fix two linearly independent integer binary forms
\(F,G\) of degree \(D\ge1\), and choose a nonzero \(2\)-by-\(2\)
minor \(R\) of their two coefficient rows.
For \(Q\ge1\), let \((a,b)\) be a representative with
\(\gcd(a,b,Q)=1\). The primitive integer parameter pairs
\((u,v)\) satisfying
\[
                  bF(u,v)-aG(u,v)\equiv0\pmod Q
\]
lie in at most \((2D)^{\omega(Q)}\) projective classes modulo
\[
 Q_*=\prod_{p^j\parallel Q}
             p^{\lceil (j-v_p(R))_+/D\rceil}
              \ \ge\ (Q/|R|)^{1/D}.
 \tag{PP1}
\]
In particular, for \(0<v\le C H\) and \(u/v\) in an interval
of length \(d\), the number of such pairs is at most
\[
                 (2D)^{\omega(Q)}
                       \bigl(1+C^2 dH^2/Q_*\bigr).
 \tag{PP2}
\]
The constant \(C\) here is any given positive height multiplier.
If \(Q\le H^K\) for fixed \(K\), the factor \((2D)^{\omega(Q)}\)
is \(H^{o(1)}\), uniformly in \(Q\).

To prove (PP1), take \(p^j\parallel Q\). The two coefficient
positions defining \(R\) show, by Cramer's rule, that the
coefficient content of \(bF-aG\) at \(p\) is at most
\(p^{v_p(R)}\), since at least one of \(a,b\) is a \(p\)-unit.
If \(j\le v_p(R)\) there is no restriction to show. Otherwise
divide by the actual content. Primitivity of the parameter
pair ensures that it lies in either the chart \(v\) a unit
or the chart \(u\) a unit. In each chart we have one univariate
polynomial over \(\mathbb Z_p\), of degree at most \(D\)
and with at least one unit coefficient. Its value on the chart
coordinate is zero modulo \(p^{j-v_p(R)}\).

In general a polynomial of this kind with values divisible by
\(p^a\), \(a\ge1\), occupies at most \(D\) solution classes
modulo \(p^{\lceil a/D\rceil}\). If \(D+1\) distinct such
classes had representatives \(t_i\) satisfying the congruence,
the denominator of each Lagrange interpolation term would
have valuation at most \(D(\lceil a/D\rceil-1)<a\).
Its numerator polynomial has integral coefficients and its
scalar value \(f(t_i)\) is divisible by \(p^a\). Every term
would have coefficients in \(p\mathbb Z_p\), contradicting
the unit coefficient of \(f\). The two charts and CRT give
(PP1). For (PP2), distinct primitive pairs with positive
denominator represent distinct rationals. In one projective
class modulo \(Q_*\), the nonzero cross determinant of two
pairs is a multiple of \(Q_*\). Their rational slopes are
therefore at least \(Q_*/(C^2H^2)\) apart.

## Two support-prime layers without projective clustering

Fix \(e\ge5\) coprime to \(6\) and \(M\ge2e\), and take
the AU1 pairs for auxiliary modulus one. Write
\[
 \ell=u^2-v^2,\qquad m=2uv,\qquad w=u^2+v^2,\qquad
 (A,B,C)=(\ell^{2e},m^{2e},w^{2e}),
\]
where \(3H\le u\le4H\) and \(H\le v\le2H\).
All three carrier values are comparable to \(X=H^{4e}\).
By AU1, at least a fixed positive multiple of \(H^2\) of
these parameter pairs give primitive complete representatives
with
\[
 \lambda=\frac{3}{2e},\qquad
 \operatorname{rad}(ABC)<C^\lambda,\qquad
 H_M(C)=w^M\le20^M H^{2M},\qquad
 \operatorname{rad}(ABC)^2\mid A+B-C.
\]
The individual radicals are \(\ell,m/2,w\), all comparable
to \(C^{\lambda/3}\). The actual defect is nonzero. In fact
\(\ell^2+m^2=w^2\) and strict convexity give \(A+B<C\).
Each complete value has at most \(H^{o(1)}\) parameter
preimages in its role. For \(\ell\) use the divisor count for
\((u-v)(u+v)\), for \(m\) use that for \(2uv\), and for \(w\)
use the two-squares bound \(r_2(w)\le4\tau(w)\).
Greedy deletion leaves \(H^{2-o(1)}\) rows sharing no
complete labels within a role. One may also make all mixed
triples across the three dictionaries pairwise coprime by
global role colouring, at \(H^{o(1)}\) loss, since each
original row is primitive. The target radicals on the
retained rows are all distinct, and the true repair cutoff
still counts these prescribed false representatives at the
full exponent \(1/M\).
[Radicals, exact repair and second-layer divisibility in AU1](ABC_OCCURRENCE_AUXILIARY_UNIT_TESTS.md#radicals-repair-and-target-multiplicity)

**PP3.** For these carriers, choose any ordered ratio of
different roles, a real interval \(I\) of length \(d\), and
one projective class modulo \(Q\le H^K\), with \(K\) fixed.
Even before the thinning, their joint occupancy is bounded by
\[
                  H^{o(1)}
                      \left(1+\frac{dH^2}{Q^{1/(4e)}}\right).
 \tag{PP3}
\]
The bound is uniform in the class, the interval, and the
modulus. In particular \(Q\) may contain support primes.

Indeed \(t=u/v\in[3/2,4]\). The three ratios of the
unpowered forms, in one orientation each, are
\[
     \frac{t^2-1}{2t},\qquad
     \frac{t^2-1}{t^2+1},\qquad
     \frac{2t}{t^2+1}.
\]
Their derivatives are, respectively,
\[
     \frac{t^2+1}{2t^2},\qquad
     \frac{4t}{(t^2+1)^2},\qquad
     \frac{2(1-t^2)}{(t^2+1)^2}.
\]
Each derivative has constant sign and absolute value bounded
below on the specified interval. All three ratios stay in
positive compact intervals. Raising a ratio to \(2e\), or
taking its reciprocal, therefore preserves these properties,
with constants depending on \(e\). This gives a real
parameter interval of length \(O_e(d)\) for all six ordered
ratios. The relevant carrier forms in \((u,v)\) have degree
\(4e\) and are linearly independent. Apply (PP1)–(PP2).
If one restricts to \(\gcd(Q,2eABC)=1\), the above
derivatives are also units at every applicable local root.
The powered derivative is a unit because the underlying
ratio is a unit and \(p\nmid 2e\). Hensel lifting then
allows \(Q\) in place of \(Q^{1/(4e)}\) in (PP3).

For a subset of at least two rows let \(Q_0\) be the gcd
of their \(A,B\) projective minors. These carriers have
pairwise coprime entries, so their pairs \((A,B)\) are
primitive. All are consequently in the projective class
of any anchor modulo \(Q_0\). The minors are nonzero on
distinct parameter pairs by real monotonicity, and
\(Q_0\ll_e H^{8e}\), in the polynomial range of (PP3).
The same argument applies to either other choice of
two roles, each with its own gcd. If the subset has
at least \(H^{2-\eta+o(1)}\) rows with \(0\le\eta<2\),
its ratio span \(d>0\) and its respective gcd satisfy
\[
            \frac{\log_H Q_0}{4e}+\log_H(1/d)
                             \le \eta+o(1).
 \tag{PP4}
\]
At full exponent and fixed \(e\), the gcd is \(H^{o(1)}\)
and the span is at least \(H^{-o(1)}\).
With \(e\) varying only after the height limit, an exponent
loss tending to zero still implies
\(\log_X Q_0\le\eta+o(1)\) on these bounded real arcs.

## Scope alongside the support-depth tests

The construction uses a single additive residual and two
support-prime layers. The original full orders on these
carriers are generally \(2e\). The higher-order test (AS1)
has a different conclusion for a bounded-resultant list
at aggregate depth three. The counterexample (AS2) reaches
full depth with a growing list of residuals and different
binary forms. The ratio pullback depends on the actual
carrier forms and does not rely on either support-depth
test. In particular (PP1) is an upper bound for every
projective cluster here and gives no count of parameters
that pass the actual additive equation.
[Aggregate depth and growing-list false carriers (AS1–AS2)](ABC_POLYNOMIAL_AGGREGATE_SUPPORT_TESTS.md)

The global-field inverse sieve of Menconi, Paredes and
Sasyk can force a bounded-degree equation on a positive
proportion of a height-bounded rational set under
sufficiently sparse reductions at a common prime family.
The cap and support tests here do not force such a
common prime family. Moreover the actual additive
projective slope is in \(\mathbb P^1\), so the codimension
hypothesis in their projective Theorem 1.6 gives no
positive-dimensional sparsity exponent on that line.
A separate bridge from true occurrence to arithmetic
clustering would still be needed.
[Menconi–Paredes–Sasyk, Theorem 1.6](https://arxiv.org/html/1907.02049#S1)
