# Factor counts and divisor marks in rough Eisenstein sums

The [one-third comparison](EISENSTEIN_ROUGH_ONE_MARK.md) has a
positive lattice margin for almost all targets. At the starting
cutoff \(N^{1/5}\) the answer depends on what is counted.
Even an exact composite test on the raw quaternion shell has
the wrong sign for almost all targets. After dividing out
lattice multiplicity, a half-weighted proper-divisor mark has
the right sign, while the first prime-factor moment still has
the wrong sign. The four-factor summands decide the last
comparison.

## The factor-count distribution

Fix \(0<\eta<1/2\) and let \(J=\lceil 1/\eta\rceil-1\).
On the pure Eisenstein branch \(N\equiv2\bmod6\), take ordered pairs
\(a+b=N\), \(a,b>1\), with \(ab\) squarefree and all its prime factors
\(1\bmod6\) and \(>N^\eta\). Let \(A_{ij}\) be the number with
\(\Omega(a)=i,\Omega(b)=j\). Each index is at most \(J\).
For \(j\ge2\) put
\[
 c_j(u)=\frac1{(j-1)!}
   \int_{\substack{x_1,\ldots,x_{j-1}\ge u\\
                    x_1+\cdots+x_{j-1}\le1-u}}
                    \prod_{k=1}^{j-1}\frac{dx_k}{x_k},
 \qquad c_1(u)=1.                                             \tag{1}
\]
In particular \(c_2(u)=\log((1-u)/u)\) for \(u<1/2\).
For \(ju\ge1\), \(c_j(u)=0\) when \(j\ge2\).

**Theorem (fixed-factor cells).** For every \(A>0\), except for
\(O_{\eta,A}(X/\log^A X)\) targets \(N\in(X,2X]\) with
\(N\equiv2\bmod6\), one has \(A_{11}\gg N/(\log N)^2\) and
simultaneously for \(1\le i,j\le J\)
\[
 \frac{A_{ij}}{A_{11}}
  =2^{\,2-i-j}c_i(\eta)c_j(\eta)
        +O_{\eta,A}\!\left(\frac{\log\log N}{\log N}\right).
 \tag{2}
\]
These constants extend the one-third prime/semiprime comparison
to the full bounded factor count. Almost-all Goldbach in the
fixed progression supplies the prime-pair baseline. The new
comparison concerns the factor cells on those same targets.

### Prime-product exponential sums

For \(N\in(X,2X]\), freeze the cutoff at either \(Z=X^\eta\)
or \(Z=(2X)^\eta\). Write \(S_{j,Z}(\alpha)\) for the sum of
\(e(\alpha n)\) over \(n\le2X\) which are products of \(j\)
distinct primes \(>Z\) and \(1\bmod6\), where \(e(v)=e^{2\pi iv}\).
Use \(P(\alpha)\) from the [one-third comparison](EISENSTEIN_ROUGH_ONE_MARK.md#circle-comparison-for-the-three-factor-cells)
for degree one, allowing all primes \(1\bmod6\) up to \(2X\).
The resulting small-prime terms at any target cost \(O_\eta(X^\eta)\)
because each small prime fixes the other summand. For arbitrary
degrees, shared primes must divide \(N\). There are \(O_\eta(1)\)
of them above \(Z\), and the cost to delete them, including the
diagonal where applicable, is \(O_\eta(X^{1-\eta})\).
Before deletion the two frozen cutoffs bracket every moving cell.

The minor-arc estimate (3) in the one-third comparison extends
to each \(S_{j,Z}\). In its bilinear argument write \(n=pr\)
with \(p\) the least prime and \(r\) the product of the other
\(j-1\) distinct split primes. The coefficient on \(r\) is an
indicator by unique factorization. The condition \(p<P^-(r)\)
adds an upper endpoint to the interval for \(p\). After
Cauchy–Schwarz and expansion in \(r_1,r_2\), the intersection
is still one interval, so the same geometric-series bound
applies. On each nonempty dyadic block \(p\asymp U,r\asymp V\)
above the discarded products \(X/\log^{K+2}X\), both variables
are \(\gg_\eta X^\eta\) and \(UV\gg X/\log^{K+2}X\).
The divided square bound is again
\[
 \ll V^{-1}+q^{-1}+\frac{\log(2q)}U
                       +\frac{q\log(2q)}{UV}
\]
for the Dirichlet approximation with
\(\log^B X<q\le\lceil X/\log^B X\rceil\).
Taking \(B\) large in terms of \(K\) and summing the dyadic
blocks gives \(\sup_{\mathfrak m}|S_{j,Z}|\ll_{\eta,K}X/\log^K X\).
The degree-one estimate is unchanged.
By the same Bessel argument (4), all products of two such sums,
including \(P\), at both cutoffs have minor-arc error
\(O(X/\log^3 X)\) outside one \(O_{\eta,A}(X/\log^A X)\) set.

On a major arc \(a/q+\xi\), first replace the distinct products
by products of \(j\) ordered prime variables and divide by \(j!\).
Repeated primes cost \(O_j(X\sum_{p>Z}p^{-2})=O_j(X/Z)\):
their product is divisible by \(p^2\), and a product has at
most \(j!\) orders. Apply Siegel–Walfisz in each variable with
modulus \([6,q]\). For a replacement with the other factors
fixed and product \(d\), the remaining variable lies in
\((Z,2X/d]\). Partial summation with phase \(e(\xi d v)\)
costs \(O(1+|\xi|X)\), so the per-class error is
\(O_C((X/d)(1+|\xi|X)(\log X)^{-C})\) for any \(C\).
Summing or integrating \(1/d\) over the other variables costs
at most \(O_\eta((\log X)^{j-1})\), even with integer harmonic
sums in the unreplaced variables. The number of residue choices
and \(1+|\xi|X\) cost fixed powers of \(\log X\).
The product of \(j\) independent uniform elements in the
subgroup of units \(1\bmod6\) is uniform in that subgroup.
Thus the arithmetic multiplier is precisely the same
\(\kappa_q(a)\) as for \(P\) and \(S_{2,Z}\).
[Siegel–Walfisz in Green, Proposition 4.2](sources/green-roth-primes-2005.pdf)

The associated continuous transform integrates \(e(\xi\prod v_i)\)
against \((2^jj!)^{-1}\prod dv_i/\log v_i\) on
\(v_i>Z,\prod v_i\le2X\). With \(t=\prod v_i\) and
\(x_i=\log v_i/\log t\), its density at \(t\) is
\[
 v_{j,Z}(t)
  =\frac{1}{2^j\log t}\,\frac1{j!}
     \int_{\substack{x_1,\ldots,x_{j-1}\ge u\\
                    x_j=1-\sum_{k<j}x_k\ge u}}
       \frac{dx_1\cdots dx_{j-1}}{x_1\cdots x_j},
 \qquad u=\frac{\log Z}{\log t}.
\]
The displayed simplex integral divided by \(j!\) equals (1).
Indeed \(\sum x_i=1\), and each of the \(j\) integrals obtained
by omitting a denominator factor is the same by symmetry of
simplex coordinates. The coordinate permutations, including
exchanging one coordinate with \(x_j\), have Jacobian of
absolute value one. Set the density to zero below \(Z^j\).

The functions \(c_j(u)\) decrease in \(u\) and are Lipschitz
on any compact positive range: the integrand is bounded on
the largest simplex and the volume removed under a change
of \(u\) is \(O(|\Delta u|)\). In particular each \(v_{j,Z}\)
has bounded variation uniformly in \(X\), since
\(c_j(\log Z/\log t)\) increases in \(t\) while \(1/\log t\)
decreases. On \([X/(\log X)^4,2X]\) the ratio
\(v_{j,Z}(t)/v_P(t)\) is
\[
 2^{\,1-j}\left(c_j(\eta)
                  +O_\eta\!\left(\frac{\log\log X}{\log X}\right)\right).
\]
In the convolution at \(X<N\le2X\), either argument below
\(X/(\log X)^4\) contributes \(O_\eta(X/(\log X)^5)\) since
the other is of order \(X\). The major-arc contribution for any
pair of degrees is then the real convolution times the same
\(\Sigma_{Q_0}(N)\) from equation (6) of the one-third proof.
There \(\Sigma_{Q_0}(N)\gg1\) and \(v_P*v_P(N)\gg X/\log^2 X\).
Taking the ratios, using the minor-arc exceptional set above,
and removing the shared and small primes proves (2).

## What the normalization changes

Write \(W=1\) for the value weight and
\(W(a,b)=6^2 2^{\Omega(ab)}\) for the lattice weight. Let
\[
 F_{N,W}(t)=\sum_{a+b=N}W(a,b)t^{\Omega(ab)},\qquad
 P_{\eta,w}(t)=\sum_{j=1}^J c_j(\eta)(wt)^{j-1},
\]
where the sum is over the restricted pairs, and \(w=1/2\) for
values, \(w=1\) for lattice witnesses. Let \(G_W=[t^2]F_{N,W}\).
Equation (2) gives, on the same good targets, uniformly for
\(0<t\le1\),
\[
 \frac{F_{N,W}(t)}{G_W t^2}
 =P_{\eta,w}(t)^2+
   O_{\eta,A}\!\left(\frac{\log\log N}{\log N}\right).            \tag{3}
\]
To see the uniformity, expand the bounded number of cells
after division by \(t^2\). Every surviving error is multiplied
by \(t^{i+j-2}\le1\). The lattice multiplicity
\(2^{\Omega(ab)}\) exactly cancels the split-prime density
penalty \(2^{2-i-j}\) after normalization by \(G_W\).
In particular \(F_{N,\mathrm{lat}}(1/2)=36 F_{N,\mathrm{val}}(1)\)
even before taking an asymptotic. This is a fixed multiplier
that realizes the value weight on the shell.

For a fixed \(w\), let \(r_j=w^{j-1}c_j(\eta)\) and
\(T=\sum_j r_j\). On the good targets the share of first
summands of degree \(j\) tends to \(r_j/T\).
The composite share therefore tends to \(1-1/T\).
The sign of the first-factor-moment certificate
\(3F_{N,W}(1)-F'_{N,W}(1)\) is governed by
\[
 H(\eta,w)=1-\sum_{j\ge2}(2j-3)r_j.                            \tag{4}
\]
The ratio of that certificate to total mass tends to \(H/T\).
Lower bounds from per-divisor tests can use a different weight.
For a retained first summand \(a\) define
\[
 B(a)=\frac12\sum_{N^\eta<p\le a/N^\eta}\mathbf1_{p\mid a}.
\]
If \(a\) is prime the sum is zero. If \(a\) is composite
each prime factor lies in the indicated range, since its
complementary factor is \(>N^\eta\). Hence \(B(a)=\Omega(a)/2\)
on composites and \(B(a)\ge\mathbf1_{a\mathrm{\ composite}}\)
on the whole family. Symmetry and inclusion–exclusion give
\[
 G_W\ge R_W-2\sum W(a,b)B(a),\qquad R_W=F_{N,W}(1).
 \tag{5}
\]
The ratio of this half-weighted proper-divisor certificate
to \(R_W\) tends to
\[
 \frac{1-\sum_{j\ge2}(j-1)r_j}{T}.                             \tag{6}
\]
The derivative assigns the same charge to semiprimes but
charges more for three or more factors. The mark (5)
uses only prime divisibility and the bound \(p\le a/N^\eta\).
Estimating it in two ranges gives an explicit switching
form: below \(\sqrt a\) it is a single norm divisibility
test, whereas above \(\sqrt a\) the substitution \(a=dp\)
leaves \(d<\sqrt N\). At \(\eta=1/5\), that cofactor has
at most two prime factors.

**Certified sign comparison.** At \(\eta=1/5\), the raw
lattice composite share lies asymptotically in
\([0.643727,0.643837]\). Even the exact composite indicator
fails the half-mass criterion there. The value-normalized
composite share lies in \([0.443339,0.443404]\).
For values, the limit in (6) is positive, between \(0.0545\)
and \(0.055\). Yet \(H(\eta,1/2)\) is in
\([-0.007339,-0.006701]\). Omitting the four-factor summands
would put it in \([0.001982,0.002586]\).
By monotonicity of the \(c_j\), the unique change in sign
of (4) as \(\eta\) increases in \([1/5,1/2)\) is in
\((0.200,0.202)\) for values and in \((0.290,0.292)\)
for raw lattice witnesses. The zero at the transition is
understood as a zero of the limiting certificate; this
comparison gives no strict sign at that exact cutoff.

Thus a pointwise one-sided estimate for the proper-divisor
mark at \(1/5\) has a positive target main term only after
controlling the multiplicity normalization. The result
isolates a viable choice of mark, rather than obtaining the
required pointwise estimate. The existing rough-shell and
marked cofactor distribution estimates do not control its
remainders. The [pair-modulus comparison](EISENSTEIN_DIVISOR_MODULUS.md)
tests whether this mark can be replaced by one- and two-prime
tests at a shorter modulus level.

## Why smaller thinning does not fix the exceptional targets

Keep a frozen cutoff as in the circle proof and let
\[
 T_t(\alpha)=P(\alpha)+
       \sum_{j=2}^J(\lambda t)^{j-1}S_{j,Z}(\alpha),
 \qquad \lambda=1\ \text{for values},\quad
        \lambda=2\ \text{for lattice witnesses}.
\]
The prime-pair coefficient of the square is then unchanged
by \(t\). Parseval gives
\(\|P\|_2,\|S_{j,Z}\|_2\ll_\eta (X/\log X)^{1/2}\), using the
one-variable estimates at frequency zero. Therefore
uniformly in \(N\) and \(0<t\le1\),
\[
 \left|\int_{\mathfrak m}
      \big(T_t(\alpha)^2-P(\alpha)^2\big)e(-N\alpha)\,d\alpha\right|
 \ \le\ \|T_t-P\|_2\,\|T_t+P\|_2
 \ \ll_\eta \frac{tX}{\log X}.                                 \tag{7}
\]
At \(t=o(1/\log X)\), all the composite terms in the
normalized minor arcs are uniformly negligible at the
prime-pair scale. What remains is precisely the ordinary
binary prime minor-arc integral in this progression. The
uniform small-\(t\) discriminant comparison transports that
quantity but supplies no estimate for its sign. The
almost-all statement (3) therefore calibrates the weighting
and the mark, while leaving the same pointwise gap.

## Numerical certificates and finite examples

Run `python3 research/arithmetic/verify_eisenstein_factor_moments.py`
from the repository root. The standard-library checker partitions
the simplex in (1) into boxes with step \((1-j\eta)/2048\).
Products of the right-endpoint reciprocal bounds on wholly
contained boxes give a lower bound. Products of the left-endpoint
bounds on boxes meeting the simplex give an upper bound.
It rounds each scaled factor down or up by exact integer division
before summation. The \(j=2\) value also uses the exact positive
series for \(2\operatorname{atanh}((r-1)/(r+1))=\log r\), with
a geometric bound for the tail. The asserted signs use rational
enclosures. The display rounds outwards to six decimals.

At \(\eta=1/5\), the run encloses
\(c_2\) in \([1.386294,1.386295]\),
\(c_3\) in \([0.405689,0.406495]\), and
\(c_4\) in \([0.014859,0.014914]\).
The numerator in (6) for values is in
\([0.098013,0.098436]\). The check also enumerates distinct norm
pairs for four sample targets using a least-prime-factor sieve.
At \(N=999998\), the unweighted first-factor-moment certificate
divided by total mass is \(+0.29722\), while its raw lattice
version is \(-0.10928\). The positive unweighted value at this
finite target contrasts with its asymptotic density-one sign.
The finite examples illustrate the slow convergence and play no
role in the sign proof.
