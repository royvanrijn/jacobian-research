# Limits of character-divisor truncation for simultaneous norms

The simultaneous Gaussian--Eisenstein weight has an exact character-divisor
expansion. Truncating its divisors at any fixed power of the component size,
however, leaves a full-order error in both the first and second absolute
moments. The truncation is also not a pointwise minorant, even inside the
required residue class. These facts identify a cancellation problem in a
natural route to the pointwise binary correlation; they do not rule out a
refined divisor method.

The unrestricted pointwise lower bound and fixed-prime marginals in the
[pointwise sieve note](SIMULTANEOUS_NORM_POINTWISE_SIEVE.md) remain open.

## The exact expansion

Let \(h(n)=r_G(n)r_E(n)/24\), where \(r_G,r_E\) count all ordered integer
representations by \(x^2+y^2\) and \(u^2-uv+v^2\), respectively. Set

\[
 h_6(n)=h(n)\mathbf1_{(n,6)=1},\qquad k=h_6*\mu.
\]

Here \(*\) is Dirichlet convolution. Thus \(h_6=1*k\). A supported value
of \(h_6\) is \(1\pmod{12}\). The Euler factors established in the
[case study](SIMULTANEOUS_GAUSSIAN_EISENSTEIN_SUMS.md#what-the-global-correlation-actually-is)
give, for \(\Re s>1\),

\[
 \sum_{n\ge1}\frac{k(n)}{n^s}
 =(1-2^{-2s})(1-3^{-2s})
 \frac{L(s,\chi_{-4})L(s,\chi_{-3})L(s,\chi_{12})}
      {L(2s,\chi_{12})}.
 \tag{1}
\]

The characters are the primitive quadratic characters of the indicated
discriminants. For \(p\ge5\), the local coefficients are

\[
 k(p^e)=
 \begin{cases}
 2e+1,&p\equiv1\ (12),\\
 e+1,&p\equiv5,7\ (12),\quad e\text{ even},\\
 -e,&p\equiv5,7\ (12),\quad e\text{ odd},\\
 (-1)^e,&p\equiv11\ (12).
 \end{cases}
 \tag{2}
\]

Indeed \(k(p^e)=h_6(p^e)-h_6(p^{e-1})\). At \(p=2,3\),
\(k(p)=-1\) and \(k(p^e)=0\) for \(e\ge2\).
In particular \(k(p)=3\) on \(p\equiv1\pmod{12}\) and \(-1\) on
the other prime classes. Define

\[
 h_D(n)=\sum_{\substack{d\mid n\\d\le D}}k(d).
 \tag{3}
\]

## Full-order remainder at every fixed power cutoff

**Theorem.** For each fixed \(0<\gamma<1\), put \(D=X^\gamma\).
As \(X\to\infty\),

\[
 \sum_{\substack{X/2<n\le X\\n\equiv1\ (12)}}
       |h_6(n)-h_D(n)|\gg_\gamma X,
 \qquad
 \sum_{\substack{X/2<n\le X\\n\equiv1\ (12)}}
       |h_6(n)-h_D(n)|^2\gg_\gamma X(\log X)^3.
 \tag{4}
\]

These are the natural first- and second-moment scales of \(h_6\).
Projecting \(h_D\) onto the residue class \(1\pmod{12}\) does not change
either lower bound.

To prove this, use the one-variable Frobenian mean theorem,
[Loughran--Matthiesen, arXiv:1904.12845v2, Lemma 2.8](https://arxiv.org/html/1904.12845v2#S2).
The functions \(h_6\) and \(h_6^2\) have prime means \(1\) and \(4\):
their prime values are \(4\) and \(16\), respectively, on a prime set of
density \(1/4\), and zero elsewhere. Both are nonnegative multiplicative
functions with polynomial prime-power growth and divisor bounds. Thus,
for positive constants \(c_1,c_2\), partial summation gives

\[
 \sum_{m\le y}\frac{h_6(m)}m\sim c_1\log y,
 \qquad
 \sum_{m\le y}\frac{h_6(m)^2}m\sim\frac{c_2}{4}(\log y)^4.
 \tag{5}
\]

Choose \(0<\varepsilon<\min(\gamma,1-\gamma,1/2)\). Consider
\(m\le X^\varepsilon\) and primes \(p\equiv1\pmod{12}\) in

\[
 X/(2m)<p\le X/m.
\]

For sufficiently large \(X\), these satisfy \(p>D\), \(p>\sqrt X\),
\(m<D\) and \(p>m\). Every divisor of \(pm\) at most \(D\) is
therefore a divisor of \(m\), and all divisors of \(m\) are included.
Consequently

\[
 h_D(pm)=h_6(m),\qquad h_6(pm)=4h_6(m).
 \tag{6}
\]

When \(h_6(m)>0\), the integer \(pm\) is \(1\pmod{12}\).
Distinct pairs \((m,p)\) give distinct integers \(pm\), since \(p>\sqrt X\)
is their unique prime factor above \(\sqrt X\).
The fixed-progression prime number theorem, equivalently
[Loughran--Matthiesen, Lemma 2.4(1)](https://arxiv.org/html/1904.12845v2#S2),
gives \(\gg X/(m\log X)\) such primes uniformly for \(m\le X^\varepsilon\).
This uniformity uses only that all interval endpoints tend to infinity at
least as fast as \(X^{1-\varepsilon}/2\); the modulus is the fixed integer 12.
Each selected integer contributes \(3h_6(m)\), or \(9h_6(m)^2\), to the
corresponding absolute-error moment. Summing and using (5) proves (4).

The same argument gives, for each fixed \(s\ge1\),

\[
 \sum_{\substack{X/2<n\le X\\n\equiv1\ (12)}}
 |h_6(n)-h_D(n)|^s
 \gg_{\gamma,s}X(\log X)^{4^{s-1}-1},
 \tag{7}
\]

because \(h_6^s\) has prime mean \(4^{s-1}\). Thus the missing large-prime
cofactor terms carry full-order mean mass, not just rare large values.

## The truncation is not a pointwise minorant

A concrete square-root example is

\[
 n=5\cdot7\cdot11\cdot17\cdot29=189805\equiv1\pmod{12},
 \qquad\lfloor\sqrt n\rfloor=435.
\]

All five prime factors are forbidden, so \(h_6(n)=0\). Among its divisors
at most 435 are the unit, all five single primes, nine prime pairs
(only \(17\cdot29=493\) is excluded), and one triple
\(5\cdot7\cdot11=385\). No divisor with more factors occurs. By (2),

\[
 h_{\sqrt n}(n)=1-5+9-1=4>0.
 \tag{8}
\]

In fact, for every fixed \(0<\gamma<1\), infinitely many squarefree
\(n\equiv1\pmod{12}\) satisfy \(h_6(n)=0<h_{n^\gamma}(n)\).
Choose all primes below in a common interval \([Y,2Y]\); the fixed-modulus
prime number theorem supplies the required distinct primes as \(Y\to\infty\).
For \(\gamma<1/2\), take two primes \(5\pmod{12}\); no nonunit divisor
is included, giving sum 1. For \(\gamma=1/2\), take five primes in classes
\(5,5,5,7,11\); exactly the subsets of size at most two are included,
giving \(1-5+10=6\). For \(\gamma>1/2\), choose an even integer \(\ell\)
with \(1<(1-\gamma)\ell<2\), and take \(\ell\) primes \(5\pmod{12}\).
Their product is \(1\pmod{12}\), and precisely subsets of size at most
\(\ell-2\) are included. Their sum is
\(\sum_{j=0}^{\ell-2}(-1)^j\binom\ell j=\ell-1>0\).
Such an even \(\ell\) exists because its permitted interval has length
greater than two. Strict gaps between powers of \(Y\) justify the stated
subset cutoffs despite the fixed factors \(2^\ell\).

## The absolute CRT rounding budget

Expanding a truncated binary convolution gives

\[
 \sum_{a+b=N}h_D(a)h_D(b)=\sum_{d,e\le D}k(d)k(e)
 \#\{a,b>0:a+b=N,\ d\mid a,\ e\mid b\}.
\]

The inner count vanishes unless \((d,e)\mid N\). Otherwise it equals
\(N(d,e)/(de)+O(1)\). Replacing each rounding error by its absolute
constant bound incurs the budget

\[
 \sum_{\substack{d,e\le D\\(d,e)\mid N}}|k(d)k(e)|
 \asymp D^2\log D,
 \tag{9}
\]

uniformly in the positive integer \(N\). This is a statement about that
absolute-error majorant, not a lower bound for the actual signed error.

Here is a justification of both bounds. Put \(G=|k|\). By (2), it is
nonnegative Frobenian with prime mean \(3/2\), \(G(n)\le\tau(n)^2\),
and \(G(pn)\le3G(n)\) at every prime. Lemma 2.8 gives
\(\sum_{d\le D}G(d)\asymp D(\log D)^{1/2}\), proving the upper bound.
For the lower bound, restrict to coprime pairs, which are included for every
\(N\). After excluding all primes at most a fixed \(Y\), the remaining
\(G\)-mass is asymptotic to \(c(Y)D(\log D)^{1/2}\), where
\(c(Y)>0\) and \(c(Y)\asymp(\log Y)^{-3/2}\). This follows from the
same mean theorem and its finite Euler factors; the last estimate uses the
prime mean and the convergent higher Euler-factor errors.
Pairs in this population sharing a prime greater than \(Y\) have total
weight at most

\[
 \sum_{Y<p\le D}\left(\sum_{\substack{d\le D\\p\mid d}}G(d)\right)^2
 \le9\sum_{Y<p\le D}\left(\sum_{u\le D/p}G(u)\right)^2
 \ll D^2\log D\sum_{p>Y}p^{-2}
 \ll\frac{D^2\log D}{Y}.
\]

Choose one sufficiently large fixed \(Y\) so this loss is smaller than
the positive pair mass; this is possible since
\(Yc(Y)^2\asymp Y/(\log Y)^3\to\infty\). Coprime pairs retain a
positive multiple of \(D^2\log D\), proving (9).

At \(D\ge\sqrt N\), this absolute rounding budget is not \(o(N)\).
At a smaller fixed-power cutoff, the budget can be small, but (4) shows that
the omitted divisors still have full-order absolute moments. Small-divisor
truncation followed only by this rounding bound and a global
Cauchy--Schwarz remainder estimate therefore does not give the required
\(o(N)\) binary correlation error. Additional signed cancellation is needed.

## A limited Hecke-recursion diagnostic

At an allowed prime, \(h(p^e)=(e+1)^2\) is not an ordinary unramified
\(\mathrm{GL}_2\) Hecke eigen-coefficient sequence. An order-two recurrence
with eigenvalue \(h(p)=4\) and \(h(1)=1\) would force its determinant
parameter to be \(16-h(p^2)=7\). It would then predict
\(h(p^3)=4\cdot9-7\cdot4=8\), whereas the actual value is 16.
Rescaling by \(p^{te}\) leaves the mismatch. This excludes only that direct
eigen-sequence model; it does not exclude a non-eigenfunction construction,
higher-rank automorphic methods or a different spectral model for the binary
correlation.

The imported mean and prime-counting statements above were checked in the
actual cited lemmas, with polynomial prime-power growth and divisor bounds
verifying their hypotheses. The conclusions are analytic obstructions to the
specified approximation steps. No pointwise binary lower bound, fixed-prime
marginal asymptotic, novelty determination, external mathematical review or
formal verification is claimed. No finite coverage search is used.
