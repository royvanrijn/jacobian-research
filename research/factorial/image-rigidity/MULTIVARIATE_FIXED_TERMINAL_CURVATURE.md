# Fixed-terminal curvature bounds for multivariate Strong Factorial

Take distinct nonnegative integer exponent vectors \(a,b\) in any
number of variables. Set \(d_i=b_i-a_i\) and
\[
 S=\sum_{i:d_i\ne0}\frac{d_i^2}{\max(a_i,b_i)} .
\]
This quantity is symmetric in the two monomials. The factorial
functional is the ordinary one, \(\Lambda(X^c)=\prod_i c_i!\).

## Small curvature at each fixed terminal

**Theorem.** For every integer \(m\ge1\) there is \(\varepsilon_m>0\)
such that if \(0<S<\varepsilon_m\), the consecutive factorial moments
of \(AX^a+BX^b\) at orders \(m,m+1\) cannot both vanish for any
\((A,B)\ne(0,0)\) over any field of characteristic zero.
The threshold is uniform in the number of variables and the exponents.

Normalize the two monomial coefficients by setting \(t=B/A\), and
write \(Q_n(t)=\Lambda((X^a+tX^b)^n)\).
The coordinates with \(a_i=b_i\) only contribute nonzero scalars
to \(Q_n\), so remove them. Whenever \(S<1\) we have
\(a_i,b_i>0\) in each remaining coordinate, since a private coordinate
would contribute a positive integer to \(S\).
With
\[
 V=\sum_i d_i^2/a_i,\qquad \sigma=\sqrt V,\qquad
 \Delta=\sum_i d_i ,
\]
one has \(S\le V\le S/(1-S)\). In fact
\(|d_i|/\max(a_i,b_i)\le d_i^2/\max(a_i,b_i)\le S\), so
\(a_i\ge(1-S)\max(a_i,b_i)\). We also use the bounds
\[
 \max_i\frac{|d_i|}{a_i}\le V,\qquad
 \max_i\frac1{a_i}\le V,\qquad
 \sum_i\frac{|d_i|}{a_i}\le V .
 \tag{1}
\]
All follow from the nonzero integer gaps \(|d_i|\ge1\).

### Gaussian limit of the centered roots

For fixed \(n\ge1\), take independent unit-rate Gamma variables \(Y_i\)
with shapes \(z_i=na_i+1\). This is the change of measure obtained
by weighting the ordinary factorial integral by \(X^{na}\).
Put \(T=\sum_i d_i\log Y_i\). Then
\[
 \frac{Q_n(t)}{Q_n(0)}=\mathbb E(1+t e^T)^n .
\]
Write \(\mu_n=\mathbb E T=\sum_i d_i\psi(z_i)\), where
\(\psi=(\log\Gamma)'\), and \(W=(T-\mu_n)/\sigma\).
Uniformly on any fixed complex disk \(|u|\le R\), its moment
generating function satisfies
\[
 \log\mathbb E e^{uW}
 = \frac{u^2}{2n}+O_{R,n}(\sigma) \quad\text{as }S\to0.
 \tag{2}
\]
For small enough \(V\), the segment from \(z_i\) to
\(z_i+u d_i/\sigma\) lies in \(\operatorname{Re}z\ge z_i/2\)
by (1). The Gamma moment formula gives the left side as
\[
 \sum_i\left(\log\Gamma(z_i+u d_i/\sigma)
      -\log\Gamma(z_i)-u d_i\psi(z_i)/\sigma\right),
\]
using the analytic logarithm of \(\Gamma\) in the right half-plane.
By the polygamma series
\(\psi'(z)=\sum_{j\ge0}(z+j)^{-2}\),
one has \(|\psi''(z)|\ll z_i^{-2}\) on that segment.
For positive real \(z_i\), comparison with an integral also gives
\(\psi'(z_i)=z_i^{-1}+O(z_i^{-2})\).
The quadratic term is therefore
\[
 \frac{u^2}{2V}\sum_i d_i^2\psi'(z_i)
 =\frac{u^2}{2n}+O_{R,n}(V),
\]
because \(\sum_i d_i^2/a_i^2\le V^2\).
The Taylor remainders are bounded by
\[
 O_{R,n}\!\left(V^{-3/2}\sum_i |d_i|^3/a_i^2\right)
 =O_{R,n}(\sigma).
\]
Here \(\sum_i |d_i|^3/a_i^2\le V\max_i |d_i|/a_i\le V^2\).
These estimates establish (2) uniformly over the supports.

Taking derivatives at zero in (2) gives convergence of the
moments of \(W\) to those of a centered normal variable \(Z_n\)
of variance \(1/n\). Indeed the moment generating functions
converge uniformly on a fixed disk and Cauchy's formula applies.
The same limit holds for the moments of
\(U=(e^{\sigma W}-1)/\sigma\).
For \(r\ge1\), the inequalities
\[
 |U|\le |W|e^{\sigma |W|},\qquad
 |U-W|\le\tfrac12\sigma |W|^2e^{\sigma |W|}
\]
give \(\mathbb E|U^r-W^r|=O_r(\sigma)\) for sufficiently small
\(\sigma\). To see this, bound the difference by
\(r|U-W|(|U|+|W|)^{r-1}\) and apply the uniform bound for
\(\mathbb E e^{2|W|}\) from (2) at \(u=\pm2\).
It follows that after setting \(t=e^{-\mu_n}(-1+\sigma x)\),
the scaled moment polynomial converges coefficientwise:
\[
 \sigma^{-n}\frac{Q_n(e^{-\mu_n}(-1+\sigma x))}{Q_n(0)}
 =\mathbb E\bigl(x+(\sigma x-1)U\bigr)^n
 \ \longrightarrow\ H_n(x):=\mathbb E(x-Z_n)^n.
 \tag{3}
\]
The limit is monic and hence all scaled roots are uniformly bounded.

### The consecutive limiting polynomials are coprime

Let \(\mathrm{He}_n\) be the monic probabilists' Hermite polynomial,
defined by \(\sum \mathrm{He}_n(y)u^n/n!=e^{uy-u^2/2}\).
Then
\[
 H_n(x)=(\mathrm i\sqrt n)^{-n}\mathrm{He}_n(\mathrm i\sqrt n\,x).
\]
For \(n=2r\), the monic square-variable polynomial
\(E_r(y)=(2r)^{-r}\mathrm{He}_{2r}(\sqrt{2r y})\) is irreducible
over \(\mathbb Q\). For \(r=1\) it is linear. For \(r\ge2\), use
Schur's even Hermite irreducibility theorem in the following precise
form. Theorem 1 of Allen and Filaseta, *Acta Arithmetica* 114 (2004),
p. 183, states that
\(\sum_{j=0}^r c_j X^{2j}/u_{2j}\) is irreducible for integer
\(c_j\) with \(|c_0|=1\) and \(0<|c_r|<2r-1\), where
\(u_{2j}=(2j-1)!!\) and \(u_0=1\).
The parametrization for \(\mathrm{He}_{2r}(X)/\mathrm{He}_{2r}(0)\)
is \(c_j=(-1)^j\binom rj\), with both endpoint coefficients units.
Irreducibility in \(X\) implies irreducibility in \(X^2\), and
rescaling the square variable by \(2r\) preserves it.
[Allen and Filaseta, Theorem 1](https://www.impan.pl/shop/publication/transaction/download/product/83056)

If \(H_{2r}\) shared a root with \(H_{2r-1}\) or \(H_{2r+1}\),
it would be nonzero since \(H_{2r}(0)\ne0\). Its negative
square would be a common root of \(E_r\) and the corresponding
monic square-variable polynomial of the odd neighbor, after
removing the odd factor. The lower neighbor has degree \(r-1\).
The upper neighbor has degree \(r\), but the coefficient of
\(y^{r-1}\) is \(-r\) there and \(-(2r-1)/2\) in \(E_r\).
These coefficients follow from the coefficient
\(-n(n-1)/2\) of \(X^{n-2}\) in \(\mathrm{He}_n(X)\).
Thus the irreducible \(E_r\) cannot divide either polynomial.

### Comparing the two centers

Suppose there is a sequence of counterexamples with \(S\to0\)
at this pair of fixed orders. At each common complex root \(t\)
write
\[
 x_n=(t e^{\mu_n}+1)/\sigma,\qquad n=m,m+1.
\]
Both coordinates are bounded by (3). Hence
\[
 e^{\mu_{m+1}-\mu_m}
   =\frac{-1+\sigma x_{m+1}}{-1+\sigma x_m}=1+O_m(\sigma).
 \tag{4}
\]
On the other hand the digamma recurrence gives
\[
 \begin{aligned}
 \mu_{m+1}-\mu_m
 &=\sum_i d_i \sum_{j=1}^{a_i}\frac1{ma_i+j}\\
 &=\Delta\log\frac{m+1}{m}+O_m(V).
 \end{aligned}
\]
For the second line, compare the inner sum to the integral
from \(ma_i\) to \((m+1)a_i\) of \(1/u\) with error at most
\(1/(ma_i)\), then use (1). Since \(\Delta\) is an integer,
(4) forces \(\Delta=0\) for sufficiently small \(S\).
The center difference is then \(O_m(V)\), and (4) gives
\(x_{m+1}-x_m=O_m(\sigma)\).
A subsequence has a common limit \(x\), which by (3) is a
common root of \(H_m,H_{m+1}\). This contradicts their
coprimality and proves the uniform existence of \(\varepsilon_m\).

Finally, a pair with a missing coefficient is a single monomial
and has nonzero factorial moments. For two nonzero coefficients,
the ordinary moment polynomials \(Q_m,Q_{m+1}\) lie in \(\mathbb Q[t]\).
Their complex coprimality implies a rational Bézout identity
and thus the conclusion over every characteristic-zero field.

## Large curvature and the endpoint radii

**Theorem.** Let \(N\ge4\), and let \(\alpha,\beta\) be distinct
nonnegative integer exponent vectors in any finite dimension. Set
\[
 S=\sum_{i:\,\max(\alpha_i,\beta_i)>0}
       \frac{(\beta_i-\alpha_i)^2}{\max(\alpha_i,\beta_i)},\qquad
 \Delta=\sum_i(\beta_i-\alpha_i).
\]
In the ordinary factorial problem at terminal \(N\), the two moment
polynomials are coprime if either
\[
 \begin{split}
 S&\ge4\log N+10,\\
 \text{or}\qquad N&\ge1024,\quad \Delta\ne0,\quad S\ge2\log N+8 .
 \end{split}
\]
Both criteria hold over every characteristic-zero field, with arbitrary
row orientations and no bound on the number of rows.

**Proof.** For \(m=N-1,N\), put
\[
 C_{m,k}=\prod_i((m-k)\alpha_i+k\beta_i)!,\qquad
 G_m=(C_{m,m}/C_{m,0})^{1/m},\qquad
 P_m(z)=\sum_{k=0}^m c_{m,k}z^k,\quad
 c_{m,k}=\binom mk\frac{C_{m,k}}{C_{m,0}G_m^k}.
\]
Here \(z=G_m t\) for the original coefficient ratio \(t\), and
\(c_{m,0}=c_{m,m}=1\). Within one row, order its exponents as \(a\le b\)
and put \(d=b-a\). The second difference of \(\log((ma+dk)!)\) at
\(1\le k<m\) is
\[
 \sum_{j=1}^d\log\!\left(1+\frac{d}{ma+d(k-1)+j}\right)
 \ \ge\ \frac{d^2}{bm}.
\]
This uses \(\log(1+u)\ge u/(1+u)\) and \(ma+dk+j\le bm\).
Reversing the row replaces \(k\) by \(m-k\) and keeps the bound.
After summing the second differences, subtract a quadratic of
second difference \(S/m\) and compare with the chord through the
endpoints. It follows that
\[
 c_{m,k}\le\binom mk\exp\!\left(-\frac{S\,k(m-k)}{2m}\right).
 \tag{A}
\]
Thus a single scalar bound applies to both orders.

Let \(\epsilon_m=\sum_{k=1}^{m-1}c_{m,k}\). When \(\epsilon_m<1\),
every root \(z\) of \(P_m\) satisfies
\[
 (1-\epsilon_m)^{1/m}\le |z|\le (1-\epsilon_m)^{-1/m},
 \qquad
 \bigl|m\arg z-(2j+1)\pi\bigr|\le\arcsin\epsilon_m
 \quad\text{for some }j\in\mathbb Z.                     \tag{B}
\]
For \(|z|\le1\), the root equation gives
\(|1+z^m|\le\epsilon_m\). The disk about \(-1\) of radius
\(\epsilon_m\) gives both the lower radius and the angle. If
\(|z|>1\), apply this argument to \(1/z\) for the reciprocal
polynomial and use symmetry of the odd multiples of \(\pi\).

For the first criterion use \(S_0=4\log N+10\) in (A). Since
\(2\log N\le N-1\) for \(N\ge4\), the upper bound for \(k=1\)
is at most \(e^{-7/3}/N<1/(8N)\). For \(2\le k\le m/2\),
write \(w=k/m\) and use \(\binom mk\le(e/w)^k\).
The logarithm of the resulting bound divided by \(k\) is
\[
 f(w)=1-\log w-S_0(1-w)/2.
\]
It is convex on \([2/m,1/2]\). At \(2/m\) use
\(S_0/m\le 2+10/4\) and \(m\le N\) to get
\(f(2/m)\le-\log(2N)+1/2<-\log N\). At \(1/2\) we have
\(f(1/2)=\log2-3/2-\log N<-\log N\). By symmetry of the bound
in (A), for \(N\ge5\) this gives
\[
 \epsilon_m<\frac{1}{4N}+2\sum_{k\ge2}N^{-k}
            \le\frac{3}{4N}.
\]
For \(N=4\) the same bound holds since the only possible additional
term in the first half is at \(k=2\), and \(2/N^2=1/(2N)\)
there. The bound also covers \(m=3\), when the interval for \(f\)
is unused. The elementary estimate \(\arcsin x\le2x\) on
\([0,1/2]\) now gives angle errors less than \(3/(2N)\).
A purported common root has angle \(\theta\in(0,\pi]\) after
conjugation. Its two normalized arguments agree, since the scalings
are positive. Subtract the two inequalities in (B).
The difference of the odd multiples must be zero, since their
difference is a multiple of \(2\pi\), \(\theta\le\pi\), and the
sum of the errors is less than \(\pi\). If the common odd multiple
is \((2j+1)\pi\), then
\[
 |(2j+1)\pi|
 \ \le\ N\arcsin\epsilon_{N-1}+(N-1)\arcsin\epsilon_N
 \ <\ 3,
\]
which is impossible. This proof does not ask the two endpoint
centers to move in the same direction.

For the second criterion assume \(2\log N+8\le S<4\log N+10\),
as the first criterion already handles larger \(S\).
Now put \(S_0=2\log N+8\) into (A).
For \(N\ge64\), \(S_0/(N-1)<1\), by monotonic decrease of the
quotient and \(\log64<5\). The first coefficient bound is at most
\(e^{-7/2}<1/25\). For the other coefficients in the first half,
the same convexity argument gives
\[
 f(2/m)\le-2-\log2,\qquad
 f(1/2)=\log2-1-\tfrac12\log N .
\]
Both exponentials are less than \(1/10\), since \(N\ge64\).
Hence \(\epsilon_m<2(1/25+1/90)<1/8\). These estimates use
only elementary exponential bounds, for example \(e^{7/2}>25\)
from the first six terms of its series and \(e>5/2\).

To compare the endpoint centers, extend \(\log G_m\) to real
\(x>0\) as
\[
 H(x)=\sum_i \frac{g(b_i x)-g(a_i x)}{x},
 \qquad g(t)=\log\Gamma(t+1),\quad a_i=\alpha_i,\quad b_i=\beta_i.
\]
Differentiating the quotient gives
\[
 H'(x)=x^{-2}\sum_i\int_{a_i x}^{b_i x}t g''(t)\,dt .
\]
Comparison of the trigamma series with its integrals gives
\(t/(1+t)\le t g''(t)\le1\) for \(t\ge0\). Thus \(H'(x)\)
differs from \(\Delta/x\) by at most
\(x^{-2}\sum_i|\log((1+b_i x)/(1+a_i x))|\).
For a sorted row \(a<b\), with \(d=b-a\), if \(d\ge b/2\) then
\(\log((1+b x)/(1+a x))\le d\log(1+x)\le2(d^2/b)\log(1+x)\).
The first bound is the product of \(d\) successive ratios, each
at most \(1+x\). If \(d<b/2\), use
\(\log((1+b x)/(1+a x))\le\log(b/a)\le d/a<2d/b\).
For \(x\ge3\) this too is at most \(2(d^2/b)\log(1+x)\).
It follows that on \(N-1\le x\le N\),
\[
 |xH'(x)-\Delta|\le \frac{2S\log(N+1)}{N-1}<\frac23.
\]
For the last inequality,
\((4\log N+10)\log(N+1)/(N-1)\) decreases for \(N\ge1024\):
the two numerator logarithmic derivatives sum to less than
\(1/N\), while the denominator has logarithmic derivative
\(1/(N-1)\). At \(1024\), the displayed bound is at most
\(2(4\cdot7+10)7/1023<2/3\), since \(\log1025<7\).

Because \(\Delta\) is a nonzero integer, \(H'\) has its sign
and \(|H'(x)|>1/(3x)\). In particular,
\[
 |\log G_N-\log G_{N-1}|>\tfrac13\log(N/(N-1))>\frac1{3N}.
\]
By (B), at a common root \(t\) the same difference
\(\log|G_N t|-\log|G_{N-1}t|\) has absolute value at most
\[
 \frac{-\log(1-\epsilon_N)}{N}
  +\frac{-\log(1-\epsilon_{N-1})}{N-1}
  <\frac17\left(\frac1N+\frac1{N-1}\right)<\frac1{3N}.
\]
The contradiction proves the second criterion. Rational Bézout
extends both conclusions from complex coprimality to every
characteristic-zero field. \(\square\)

## An all-order bound from private variables and shared curvature

**Theorem.** Suppose each monomial has at least one variable absent
from the other. Set
\[
 S_{\mathrm{sh}}=\sum_{i:\,a_i b_i>0}
                   \frac{(b_i-a_i)^2}{\max(a_i,b_i)}.
\]
If \(S_{\mathrm{sh}}\ge4\log8\), the two monomials satisfy the ordinary
Strong Factorial assertion at every terminal over every
characteristic-zero field. This permits private exponent sums one
on both sides.

**Proof.** The normalized private factorial contributions cancel
the binomial coefficient in the sense of the
[private-variable coefficient bound](MULTIVARIATE_TWO_MONOMIAL_EVENTUAL_BORDER.md#coefficient-estimates).
For shared rows the coefficient proof of (A) gives an additional
factor \(\exp[-S_{\mathrm{sh}}k(m-k)/(2m)]\). Thus
\(c_{m,k}\le q^{\min(k,m-k)}\), where
\(q=\exp(-S_{\mathrm{sh}}/4)\le1/8\).
By [the initial terminal proofs](MULTIVARIATE_TWO_MONOMIAL_NONDIVISIBILITY.md),
it suffices to consider lower order \(m\ge4\).

Split \(P_m(z)=L(z)+z^m U(1/z)\) at degree \(h=\lfloor m/2\rfloor\),
so \(L,U\) have degree at most \(h\), constant coefficient one
and coefficient bounds \(q^j\) in degree \(j\ge1\).
On the unit disk both moduli are at least \((1-2q)/(1-q)\)
and at most \(1/(1-q)\). Consequently at any root \(z\) of
\(P_m\), with \(r=|z|\), the equation \(L(z)=-z^m U(1/z)\)
gives
\[
 1-2q\le r^{m-h}\le (1-2q)^{-1}.
\]
For \(r\ge1\) use the lower bound on \(|U(1/z)|\) and
\(|L(z)|\le r^h/(1-q)\); for \(r\le1\) use the lower
bound on \(|L(z)|\) and \(|U(1/z)|\le r^{-h}/(1-q)\).
Since \(m-h\ge m/2\) and \(m\ge4\), we have
\(r^{\pm1}\le (1-2q)^{-2/m}\le\sqrt{4/3}<7/6\).
Thus in both \(L(rw)\) and \(U(r^{-1}w)\) the coefficient in
degree \(j\) is at most \(p^j\), where \(p=7/48\).

For \(0\le\theta\le\pi\) the argument of either polynomial
evaluated at \(w=e^{i\theta}\) is less in absolute value than
\(\theta/4\) when \(\theta>0\).
Indeed its real part is at least \(1-p/(1-p)>0\),
and the imaginary part is bounded in absolute value by
\(\theta\sum_{j\ge1}jp^j\). Therefore the absolute argument is
at most
\[
 \theta\frac{p}{(1-p)(1-2p)}<\frac{\theta}{4}.
\]
At a root, the phase equation is \(m\theta=(2j+1)\pi+\delta_m\)
where \(\delta_m=\arg L(re^{i\theta})+\arg U(r^{-1}e^{i\theta})\).
In particular \(|\delta_m|<\theta/2\) for \(0<\theta\le\pi\).
A common root of the consecutive moment polynomials has the same
argument in the two positive normalizations. Conjugate it if
necessary to put its argument in \((0,\pi]\); positivity excludes
argument zero. Subtracting the phase equations would place
\(\theta-\delta_{m+1}+\delta_m\) strictly between 0 and \(2\pi\)
while requiring it to be a multiple of \(2\pi\).
This is impossible, and rational Bézout completes the proof.
\(\square\)

## The limiting pair when each relative gap tends to zero

The two curvature bounds confine any counterexample at a fixed
terminal \(N=m+1\ge5\) to a compact positive interval of \(S\).
The next result describes accumulation inside that interval when
the gaps become small relative to their own rows.

**Theorem.** Fix \(m\ge1\), set \(N=m+1\), and put
\[
 \eta=\max_{i:d_i\ne0}\frac{|d_i|}{\max(a_i,b_i)},\qquad
 B_n(z;s)=\sum_{k=0}^n\binom nk
           e^{-s k(n-k)/(2n)}z^k.
\]
Suppose a sequence of ordinary two-monomial supports has bounded
\(S\), has \(\eta\to0\), and has a common root at this fixed
pair of orders for every support. After passing to a subsequence,
\(S\to s>0\), \(\Delta=\sum_i d_i\) is constant, and
\(L=\sum_{i:d_i\ne0}\log(b_i/a_i)\) tends to a limit \(\ell\).
Then \(B_m(z;s)\) and \(B_N(z;s)\) have a common root, and
\[
 \ell=2mN\Delta\log(N/m),\qquad |\ell|\le s.                \tag{C}
\]
For \(N\ge11\) this forces \(\Delta=\ell=0\).

In particular, there are at most finitely many possible \(s>0\)
at a fixed terminal. Define
\[
 R_m(u)=\operatorname{Res}_z\left(
       \sum_{k=0}^m\binom mk u^{N k(m-k)} z^k,\
       \sum_{k=0}^N\binom Nk u^{m k(N-k)} z^k\right).
 \tag{D}
\]
A necessary condition is \(R_m(e^{-s/(2mN)})=0\). If \(R_m\)
has no zeros in \((0,1)\), there exists \(\eta_m>0\) such that
all supports with \(\eta<\eta_m\) satisfy the two-moment statement
at this terminal, uniformly in the exponents and dimension.

**Proof.** Remove rows with equal exponents. As \(\eta\to0\),
all remaining exponents are positive. Put \(r_i=\min(a_i,b_i)\).
Then
\[
 \max_i |d_i|/r_i\le\eta/(1-\eta),\qquad
 S\le \sum_i d_i^2/r_i\le S/(1-\eta).
\]
The same first bound applies to \(\max_i1/r_i\) since the gaps are
nonzero integers. Take a subsequence with \(S\to s\).
Small-curvature separation proves \(s>0\).

For fixed \(n=m,N\), use the normalization of (A).
Taylor expansion of \(g(x)=\log\Gamma(x+1)\) about \(na_i\)
in each row, followed by subtracting the endpoint chord at \(k\),
gives
\[
 \log\frac{c_{n,k}}{\binom nk}
 =-\frac{k(n-k)}2\sum_i d_i^2 g''(na_i)
    +O_n\!\left(\sum_i\frac{|d_i|^3}{r_i^2}\right).
\]
The derivative bound used here is \(g'''(x)=O(x^{-2})\) on the
positive segment, by the differentiated trigamma series.
Also \(g''(na_i)=1/(na_i)+O_n(r_i^{-2})\) by series-integral
comparison. The total errors are \(O_n(S\eta/(1-\eta)^2)\).
In particular, \(P_n(z)\to B_n(z;s)\) coefficientwise. These
are Kabluchko's unitary Hermite polynomials in (1.5), with the
substitution \(-z\), scalar factor \((-1)^n\) and parameter
\(\sigma^2=s/n\). Thus their roots all have modulus one by
Lemma 2.1. [Kabluchko, (1.5) and Lemma 2.1](https://doi.org/10.5802/ahl.227)

The normalized roots at a common coefficient ratio \(t\) are
\(G_m t\) and \(G_N t\). Along a further subsequence both have
limits, since the normalized polynomials converge to monic
polynomials of fixed degree. Their limits have modulus one.
As \(G_N/G_m\) is a positive real number, it tends to one and
the two limits are the same. This proves the common-root
condition for \(B_m,B_N\).

We give the center comparison to account for a term that can
survive even when every relative gap vanishes. Since
\(\log G_n=\sum_i\int_{a_i}^{b_i}\psi(nx+1)\,dx\), we have
\[
 \log(G_N/G_m)
 =\Delta\log(N/m)-\frac{L}{2mN}
    + O_m\!\left(\sum_i\frac{|d_i|}{r_i^2}\right).        \tag{E}
\]
In detail, the series for \(\psi'(v+1)\) gives
\(\psi'(v+1)=1/v-1/(2v^2)+O(v^{-3})\). To obtain the second term,
compare the integral of \(f(y)=(v+y)^{-2}\) on \(y\ge0\) with its
trapezoidal sum. The error is bounded by
\(\tfrac18\int_0^\infty |f''(y)|\,dy\), by integrating the
second derivative against the interpolation kernel on each unit
interval. Thus \(\sum_{j\ge1}f(j)=\int_0^\infty f-\tfrac12f(0)
+O(v^{-3})\). Integrating this expansion from \(mx\) to \(Nx\)
gives \(\psi(Nx+1)-\psi(mx+1)=\log(N/m)-1/(2mNx)+O_m(x^{-2})\).
This proves (E), including in decreasing rows.

The error in (E) tends to zero since
\(\sum_i |d_i|/r_i^2\le \max_i(1/r_i)\sum_i d_i^2/r_i\).
Also \(|L|\le\sum_i|d_i|/r_i\le S/(1-\eta)\).
Now \(G_N/G_m\to1\), so (E) shows that the integers \(\Delta\)
are bounded. Choose a subsequence with constant \(\Delta\) and
take limits to obtain (C). At \(N\ge11\), the first large-curvature
criterion gives \(s\le4\log N+10<2(N-1)\). Indeed the strict
inequality holds at \(N=11\) by \(\log11<5/2\), and its margin
increases. Since \(2mN\log(N/m)>2m\), (C) gives
\(\Delta=0\) and hence \(\ell=0\).

At \(u=0\) the two inputs to the resultant in (D) are
\(1+z^m,1+z^N\). They have no common root since \(N-m=1\).
Thus \(R_m\) is a nonzero integer polynomial. The substitution
\(u=e^{-s/(2mN)}\) agrees with the two limiting polynomials above,
and injects \(s>0\) into \((0,1)\). This proves finiteness of the
candidate set. If it is empty, a sequence of counterexamples
with \(\eta\to0\) is impossible: \(S\) stays bounded by the first
large-curvature criterion when \(N\ge4\), and terminals 2 and 3
are already closed. Over other characteristic-zero fields use
rational Bézout as above. \(\square\)

### Exact exclusions through terminal thirteen in the limiting model

For each \(2\le m\le12\), the polynomial \(R_m\) has no zeros in
\((0,1)\). Consequently, at each of these fixed terminals there is
a uniform relative-gap threshold \(\eta_m>0\) as stated in the theorem.
The computation concerns the one-parameter limiting model in (D).
It does not certify all supports at those terminals.

Here is the exact algorithm. Treat both inputs to (D) as reciprocal
polynomials in \(z\), over \(\mathbb Z[u]\). Remove \(z+1\) from the
odd-degree member by exact division. Each remaining even reciprocal
polynomial of degree \(2r\) has the form \(z^r F(z+z^{-1})\).
To obtain \(F(x)\), convert the pairs \(z^j+z^{-j}\) using the
recurrence \(p_0=2\), \(p_1=x\), \(p_{j+1}=xp_j-p_{j-1}\).
The original polynomials have nonzero constant terms. Any common
root other than \(-1\) would therefore give a zero of the resultant
of the two angle polynomials in \(x\). The omitted root is checked
separately by substituting \(x=-2\) into the even member.

For both the angle resultant and the separate evaluation, remove all
powers of \(u\) and \(u-1\) by exact division. For the residual \(P(u)\)
of degree \(D\), certify that every coefficient of
\[
 (1+v)^D P(v/(1+v))
\]
is strictly of one sign. The substitution maps \(v>0\) onto
\(0<u<1\). For coefficient \(k\), the verifier computes
\(\sum_{i\le k}P_i(D-i)!D!/(k-i)!\) by integer convolution.
This is the transformed coefficient multiplied by the positive
integer \(D!(D-k)!\), so the check uses only exact integers.

Replay with python-flint 0.9.0:

```sh
.venv/bin/python research/factorial/image-rigidity/verify_multivariate_diffuse_resultants.py
```

All eleven orders passed both checks. The reduced angle-resultant
degrees for \(m=2,\ldots,12\) were
\(3,15,51,119,250,450,770,1210,1845,2665,3773\).
The verifier also accepts a nonvanishing polynomial with both
endpoint factors, and rejects both the zero polynomial and
\(2u-1\), which has an interior root.
