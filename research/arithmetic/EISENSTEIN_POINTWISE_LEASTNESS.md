# Pointwise rough Eisenstein pairs and the cost of leastness

At roughness \(N^{0.116}\) the existing quaternionic distribution
already gives a pointwise lower bound on the distinct-value shell.
An elementary progression sieve bounds the composite mass missed by
a least-prime cutoff just below \(N^{1/2}\). The remaining marked counts
need a different treatment of leastness. Even a square-root saving over
all direct raw vector-sieve moduli would leave too few useful cross marks.

## The rough shell

Let \(N\equiv2\bmod6\), \(z=N^\eta\), and let \(\mathcal R_N\) consist of
ordered pairs \(a+b=N\) with \(a,b>1\), \(ab\) squarefree, and all primes
of \(ab\) split in \(\mathbb Q(\sqrt{-3})\) and greater than \(z\). Write
\(R_N=|\mathcal R_N|\), and put
\[
 \mathfrak S(N)=\prod_{\substack{p\mid N\\p>2}}\frac{p-1}{p-2},
 \qquad X_N=\frac{N\mathfrak S(N)}{(\log N)^2}.
\]

**Theorem (pointwise rough shell).** For every fixed
\(0<\eta<1/(2\cdot4.267)\), uniformly over all sufficiently large
\(N\equiv2\bmod6\), one has \(R_N\asymp_\eta X_N\).
The lower bound already holds on the fixed central triangle
\(|a-b|<N/2\).

For the lower bound take the odd coordinate shell and the triangular
weight \(b_{1/2}\) from
[EisJointCosetDistribution](EISENSTEIN_BOUNDED_PRIME_FACTORS.md#all-target-valuations-and-joint-cofactor-levels).
The main mass is \(R_o(N)/4\), where
\(R_o(N)=12(\sigma(N)-\sigma(N/4))\), with the second term zero if
\(4\nmid N\). Uniformly in the two-adic valuation, this is comparable
to \(\sigma(N)\). Also \(\sigma(N)/N\asymp\mathfrak S(N)\), since
each odd local Euler-factor ratio is \(1+O(p^{-2})\), uniformly in its
positive valuation. In the notation of that theorem the sifting
density is \(4/p+O(p^{-2})\) at split primes and \(O(p^{-2})\) at
inert primes, including primes dividing \(N\). Relative to the
fixed-progression Mertens product the errors converge absolutely,
and the local factors at small primes stay bounded away from zero.
This supplies the dimension-two product condition uniformly and
\(V_N(z)\asymp_\eta(\log N)^{-2}\).

Choose \(4.267\eta<\xi<1/2\) and use the dimension-two DHR lower
bound at \(D=N^\xi\). Its main term is \(\gg_\eta X_N\). With the
cofactor budget set to one, (J4) bounds the weighted remainder by
\(O(N^{1/2+\xi+o(1)})=o(X_N)\). The factors exponential in
\(\omega(d)\) cost \(N^{o(1)}\) uniformly for \(d\le D\).
After this sift, tuples with a square dividing a component, or
with a common prime divisor, have a prime \(p>z\) with \(p^2\mid a\)
or \(p^2\mid b\), or a divisor \(p>z\) of \(N\) dividing both.
There are \(O_\eta(N/z+\sqrt N)\) such value pairs by elementary
divisibility, and each has \(N^{o(1)}\) coordinate witnesses.
Inert primes in an Eisenstein norm occur to even exponent, so the
same deletion removes all remaining inert factors. On the surviving
squarefree split pairs the coordinate multiplicity is
\(36\,2^{\Omega(ab)}\le36\,2^{2/\eta}\) and \(b_{1/2}\le1\).
This gives the lower bound with the full target factor.

## Uniform least-prime boundary estimates

Here is the progression bound used for the upper bound and the
least-prime marks. Set \(M=1,p,q,\) or \(pq\) for zero, one, or two
prescribed split primes \(p,q>z\) dividing the designated summands.
In the two-mark case take different primes, prime to \(N\). For one
mark on the coprime shell the marked prime is also prime to \(N\).
With \(M\le N^c\) for any fixed \(c<1\), the conditions put \(a\)
in one residue class modulo \(M\). Over this progression the number
of pairs with both summands \(z\)-rough is
\[
               \ll_{\eta,c} \frac{X_N}{M}.                 \tag{1}
\]
The progression length is \(Y=N/M+O(1)\). For odd primes \(r\le z\)
the two linear forms have two forbidden classes, or one when
\(r\mid N\). Sift at \(y=Y^{\min(\eta/2,\,1/100)}\) with
upper-sieve level \(Y^{1/4}\). The dimension-two upper bound
at this fixed ratio applies with the uniform product condition
from the two root counts. The remainder uses \(O(2^{\omega(d)})\)
from the progression and is \(O(Y^{1/4}(\log Y)^{O(1)})\).
The sifting product is at most
\(O(\mathfrak S(N)/(\log y)^2)\), by the ordinary Mertens product
and the extra factor \((r-1)/(r-2)\) when \(r\mid N\).
This proves (1). The \(M=1\) case proves the upper bound for \(R_N\).

**Proposition (pointwise boundary control).** Fix \(\eta>0\) in the
range above. Uniformly for \(\eta\le u\le v\le1/2\), the fraction
of \(\mathcal R_N\) with \(N^u<P^-(a)\le N^v\) is
\[
              O_\eta\!\left(\log(v/u)+(\log N)^{-1}\right).
                                                               \tag{2}
\]
For any fixed \(c<1\), uniformly for \(2\eta\le u\le v\le c\),
the fraction with \(N^u<P^-(a)P^-(b)\le N^v\) and both least
primes \(>z\) is
\[
              O_{\eta,c}\!\left(v-u+(\log N)^{-1}\right).       \tag{3}
\]
In (2), sum (1) over \(p\) and apply the prime reciprocal Mertens
estimate. In (3), sum (1) over \(p,q\) restricted by the product
bound. For each \(p\le N^{c-\eta}\), the interval for \(q\) has
lower log-exponent at least \(\eta\) and width at most \(v-u\).
The accompanying sum of \(1/p\) is \(O_\eta(1)\).
Inclusion or exclusion of the interval endpoints costs
\(O_\eta(N^{-\eta})\) in these harmonic bounds.

This also allows the factor bins for cross marks to be stated
using \(\log p/\log N\) pointwise. At most \(2N/(\log N)^4+O(1)\)
pairs have either summand below \(N/(\log N)^4\). Elsewhere the
two log ratios for a least prime differ by
\(O(\log\log N/\log N)\). For a fixed list of bin endpoints below
\(1/2\), apply (2) to all the resulting thin bands. The overall
change of the cross mark is \(o(R_N)\), uniformly in \(N\).

Let \(C_N(B)\) count pairs whose first summand is composite with
least prime at most \(B=N^{1/2-\delta}\), where \(0<\delta<1/10\).
Every composite has least prime at most \(\sqrt N\).
By (2) its omitted band \(H_N=C_N-C_N(B)\) satisfies
\[
                       H_N\ll_\eta (\delta+1/\log N)R_N.      \tag{4}
\]
Every pair in this band has a semiprime first summand. This
estimate also follows by writing \(a=pr\), with the smaller prime
\(p\) in the thin band and applying (1) to \(r,N-pr\).
The constants in (1)–(4) are independent of \(\delta\).

## Signed reduction and the leastness budget

Use \(\eta=0.116\) and the bins of the
[budget comparison](EISENSTEIN_BUDGET_INTERSECTION.md#the-one-sided-and-two-sided-shape-limits),
now stated in \(\log p/\log N\). Let \(Q_N\) be their cross mark.
Choose a fixed sufficiently small \(\delta>0\) so that
\(B>N^{0.495}\) and \(H_N\le0.0002R_N\) eventually, using (4).
A sufficient signed pointwise estimate is
\[
                        R_N-2C_N(B)+Q_N\ \ge\ 0.001R_N.      \tag{5}
\]
The exact inclusion–exclusion inequality then gives
\(G_N\ge0.0006R_N\) on every sufficiently large admissible target.
The almost-all limiting margin of the left side divided by \(R_N\)
is at least \(0.004498\) for the stated bins, because truncating
the negative marginal can only increase it and the bin change is
\(o(R_N)\). The lower bound on \(R_N\) and the cutoff loss in (5)
are supplied pointwise above.

Leastness in (5) has an exact switched expression. Write
\(a=pr\) for the first mark and \(b=qs\) for the second mark.
In the single count sum over \(z<p\le B\) with \(r>1\),
\(P^-(r)>p\). In the cross count sum over the ordered bins with
\(pr+qs=N\), \(r,s>1\), \(P^-(r)>p\) and \(P^-(s)>q\).
The cofactor variables retain split support. The inequality \(r>1\)
in the single count is necessary to omit prime summands, and
\(r,s>1\) are automatic in the cross count outside the negligible
endpoint range. Each accepted pair has exactly one least prime
per marked summand. The coordinate multiplicity removed by
prescribing \(p,q\) accounts for the two orientations at the
marked primes only. The residual normalization divides by
\(2^{\Omega(r)+\Omega(s)}\) for the two-mark count.
Thus cancellation for (J4) at the original cutoff \(z\) alone
does not evaluate (5).

For example, let \(x=\log p/\log N\), \(y=\log q/\log N\).
A direct raw coordinate vector sieve for the two leastness
conditions has dimension one on each side. Away from the target
and the marked primes, the density of a split prime dividing a
specified component is \(2/p+O(p^{-2})\), with a summable inert
contribution. Using linear lower sieves on both sides in the
usual vector-sieve combination asks for levels above \(N^{2x}\)
and \(N^{2y}\). The product of these two levels and the mark
product has exponent \(3(x+y)\). Even hypothesizing a signed
remainder \(\sqrt N\) times the square root of that product,
the direct product budget only covers \(x+y<1/3\).
The [missed-mass bound](EISENSTEIN_BUDGET_INTERSECTION.md#does-normalization-lower-the-sieve-dimension)
at \((1/8,0.34)\) gives missed mass \(>1.37250\), so the
smaller cap at \(\eta=0.116\) has a negative shape margin even
if the complete marginal were known.

**Proposition (a common leastness cutoff also loses the margin).**
Suppose instead we grant level estimates for the distinct-value
local model, and impose both leastness conditions by one
dimension-one sieve at \(\max(p,q)\). In the direct product
budget, even a hypothetical square-root saving in the signed
aggregate leaves a negative limiting margin for every
\(\eta\le1/(2\cdot4.267)\), granting the complete marginal.

Indeed a normalized dimension-one lower sieve at that cutoff
uses exponent exceeding \(2\max(x,y)\). For the hypothesized
saving, the allowed cross shape has \(x+y+2\max(x,y)\le1\).
This already implies \(x+y\le1/2\). Enclosing this entire
shape at the largest allowed \(\eta=500/4267\) gives a
normalized upper margin below \(-0.0261\). In unnormalized
units its sign is the sign of one minus the missed
two-composite mass. The latter is nondecreasing as \(\eta\)
decreases, so the upper sign holds throughout the stated range.

One can also give each side its own cutoff in the normalized
local model, with dimension \(1/2\) per side. Take a common
upper/lower level on each side, with ratios \(s_a,s_b\)
relative to \(p,q\). The standard vector-sieve lower main
term can be positive only when
\(f_{1/2}(s_a)/F_{1/2}(s_a)+f_{1/2}(s_b)/F_{1/2}(s_b)>1\).
The direct product ceiling with the same hypothetical saving
requires \(x s_a+y s_b\le1-x-y\). For \(\eta\ge0.08\) and
\(t=x+y\ge0.47\), this leaves at most \(1-2t\le0.06\)
beyond the two baseline ratios one, so both ratios are below
two if both lower terms are positive. If either lower term
is zero the vector coefficient is nonpositive.
[Iwaniec's formulas](sources/iwaniec-half-dimensional-sieve-1976.pdf)
for \(1\le s\le2\) give
\[
 \frac{f_{1/2}(s)}{F_{1/2}(s)}
    =\operatorname{arsinh}\sqrt{s-1}\le\sqrt{s-1}.
\]
Cauchy–Schwarz now makes positivity require
\(x(s_a-1)+y(s_b-1)>xy/t\).
But \(xy/t\ge0.08(1-0.08/t)>0.06\), since both variables
are at least \(0.08\). Thus this recipe permits only \(t<0.47\).
The full \(0.47\) cap at \(\eta=500/4267\) already has
normalized upper margin below \(-0.0089\), as certified below.
For \(\eta<0.08\) the cheaper obstruction is the full
\(0.5\) cap from the earlier missed-mass bound. This argument
uses the same upper and lower level for a given side.

**Proposition (unequal upper and lower vector levels).**
In the normalized half-dimensional model, allow different
upper and lower levels on each side. The usual vector-sieve
lower bound with the same hypothetical signed square-root
saving and direct product accounting still has a negative
limiting margin for square-root cross marks, throughout
\(0<\eta\le500/4267\), with the full marginal granted.

To see this, let the upper ratios to the two prime cutoffs
be \(u_a,u_b\), and the lower ratios \(s_a,s_b\). With
\(t=x+y\) and \(L=1-t\), the two mixed products require
\[
               xs_a+yu_b\le L,\qquad xu_a+ys_b\le L.
\]
Dropping the upper-upper product only improves this feasibility
test. Positivity still requires
\(f_{1/2}(s_a)/F_{1/2}(u_a)+f_{1/2}(s_b)/F_{1/2}(u_b)>1\).
If either \(s\le1\) its term vanishes and the other term is
at most one. For \(x,y\ge0.08\), it suffices to exclude the
boundary \(t=0.475\), since any larger sum can be reduced to
that value keeping \(x,y\ge0.08\) and weakening both level
constraints. On that boundary, \(u_a,u_b<1.625\) by \(s_a,s_b>1\).
For \(0<u\le2\), the sieve ratio is
\(\sqrt{u}\,g(s)\), where \(g(s)=f_{1/2}(s)/(2\sqrt{e^\gamma/\pi})\).
At \(1<s<3\),
\[
           g(s)=\frac{\operatorname{arsinh}\sqrt{s-1}}{\sqrt s},
\]
which is increasing there by differentiation. For \(s\ge3\)
it is at most \(2/3\), since \(f_{1/2}(s)\le1\) and
\(e^\gamma>9\pi/16\).

The checker bounds the maximum sum of the two ratios by \(0.983\).
It uses \(x\le y\) by symmetry and covers
\(0.08\le x\le0.2375\) in bins of width \(1/400\).
With \(x_0,y_0\) the lower bounds of \(x,y\) on a bin, a
pair of upper-ratio cells with lower endpoints \(u_{a,0},u_{b,0}\)
has lower ratios bounded above by
\((L-y_0u_{b,0})/x_0\) and \((L-x_0u_{a,0})/y_0\).
The enclosing calculation uses those bounds on \(g\) and
the upper endpoints for the square roots of \(u_a,u_b\).
It rounds square roots upward and uses the positive atanh
series with its geometric tail for \(g\).
Thus the surviving product cap is at most \(0.475\).
Its full profile has normalized upper margin below
\(-0.0054\) at \(\eta=500/4267\), hence has negative sign
throughout \(0.08\le\eta\le500/4267\) by monotonicity of
the missed mass. When \(\eta<0.08\), the full \(0.5\) cap
itself has negative sign by the earlier missed-mass bound.

These product-budget arguments concern direct lower sieves
for leastness. A switched estimate on the cofactor sequence,
or an estimate that combines unequal sieve cutoffs without
charging these products, would access different shapes.

## Exact replay and sources

The [exact checker](verify_eisenstein_pointwise_leastness.py)
uses the rational CDF enclosure from the budget comparison
to bound the polygon \(3x+y\le1,\ x+3y\le1\) and the full \(0.47\) cap. As a further test,
it restricts the credit to bins with both marks below \(N^{1/3}\).
At \(\eta=0.117\), with product cap \(N^{0.499}\) and bin width
\(1/2000\), those smaller marks still give a limiting margin
in \([0.001654,0.001868]\). Thus the cross term can omit the
branches whose leastness requirement forces the cofactor to
be prime just from its size. The visible marginal still needs
control. At \(\eta=0.116\) the semiprime first-summand band with
least exponent \(>1/3\) has limiting fraction between
\(0.143759\) and \(0.143775\), so omitting all of it would
use more than the margin.

The progression upper sieve uses the standard dimension-two
upper-bound setup in
[Johnston, v3, (1.12)–(1.15)](sources/johnston-2510.10853v3.pdf),
with residue-counting remainders in place of prime distribution.
The lower input is [the same reference, (6.11)–(6.12)](sources/johnston-2510.10853v3.pdf).
[Rosser–Schoenfeld, Theorems 5, 7 and 8](sources/rosser_schoenfeld_1962.pdf)
supplies the ordinary prime reciprocal and product estimates.
The local coordinate densities and all-valuation summation used
above are specified in (3)–(8) and (J1)–(J5) of
[EisJointCosetDistribution](EISENSTEIN_BOUNDED_PRIME_FACTORS.md).