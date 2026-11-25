# Fractional-index transport for ordinary factorial coefficients

Let \(a,b\) be two distinct nonnegative integer exponent columns in
finite dimension. Write \(n\ge1,\ N=n+1,\ \alpha=N/n\), and put
\(d_i=b_i-a_i\). For \(s\in[0,1]\) define
\[
 X_i(s)=(1-s)a_i+s b_i,\qquad
 H_j(s)=\sum_i\left[\log\Gamma(1+jX_i(s))
   -(1-s)\log\Gamma(1+ja_i)-s\log\Gamma(1+jb_i)\right].
\]
Then \(\log F_{j,k}=H_j(k/j)\), where \(F_{j,k}\) is the
endpoint-normalized multiplier after removing the binomial
coefficient. Also define the shifted logarithmic Jensen gap
\[
 L_\delta(s)=\sum_i [\log(\delta+X_i(s))
          -(1-s)\log(\delta+a_i)-s\log(\delta+b_i)],
 \qquad \delta>0.
\]

**Theorem.** At every \(s\in[0,1]\),
\[
 \frac{L_{1/n}(s)}{2n}
   \ \le\ \alpha H_n(s)-H_N(s)\
   \le\ \frac{L_{1/(3N)}(s)}{2n}.
 \tag{1}
\]
The inequalities are strict in the interior. There exists a positive
random variable \(Y\), with mean one, such that
\[
 \exp(H_N(s)-\alpha H_n(s))=\mathbb E[Y^s]\qquad(0\le s\le1).
 \tag{2}
\]
Its logarithm has an explicit nonnegative jump measure. These claims
include profiles with private rows. They compare the two normalized
coefficient functions at the same fractional index, independently
of the radii used to normalize the two polynomials.

## The comparison kernel

Put
\[
 D_n(x)=\frac1{e^{x/N}-1}-\frac{N/n}{e^{x/n}-1}.
\]
For all \(x>0\),
\[
                    \frac{e^{-x/n}}{2n}
                   <D_n(x)<\frac{e^{-x/(3N)}}{2n}.
 \tag{3}
\]
To prove this, put \(g(v)=v/(e^v-1)\). Then
\[
 D_n(x)=\frac{N}{x}\int_{x/N}^{x/n}[-g'(v)]\,dv
       =\frac1n\operatorname{avg}_{[x/N,x/n]}[-g'(v)].
\]
The bounds used on the integrand are
\[
               \frac12 e^{-v}<-g'(v)<\frac12 e^{-v/3}
               \qquad(v>0).
\]
For the lower bound, clearing denominators gives the positivity
of \((2v-3)e^{2v}+4e^v-1\). This function vanishes at zero
and its derivative is \(4e^v((v-1)e^v+1)>0\).
For the upper bound clear denominators to compare
\[
 e^{2v}-2e^v+1-2(v-1)e^{4v/3}-2e^{v/3}
\]
to zero. For \(m\ge1\) let \(A_m\) be its Taylor coefficient at
\(v^m\), multiplied by \(m!3^m/2\). Direct expansion gives
\(A_m=6^m/2-3^m-(3m/4-1)4^m-1\), with \(A_1=A_2=A_3=0\).
The constant coefficient vanishes as well. For \(m\ge3\),
\[
 A_{m+1}-4A_m=6^m+3^m-3\cdot4^m+3>0
\]
since \(6^3>3\cdot4^3\). Hence the Taylor series is positive
for \(v>0\). Averaging the two strict bounds proves (3).

The [trigamma series](https://dlmf.nist.gov/5.15.E1) and the
integral of \(v e^{-t v}\) give
\[
 \psi'(1+t)=\int_0^\infty \frac{v e^{-t v}}{e^v-1}\,dv
 \qquad(t\ge0).
\]
The two derivatives in \(s\), followed by the zero endpoint
conditions, therefore give the exact chord integral
\[
 H_j(s)= \sum_i\int_0^\infty
 \frac{e^{-X_i(s)x}-(1-s)e^{-a_i x}-s e^{-b_i x}}
                  {x(e^{x/j}-1)}\,dx .
 \tag{4}
\]
In detail, the second derivative of (4) is
\(\sum_i j^2 d_i^2\psi'(1+jX_i(s))\) by the substitution
\(x=jv\), so it agrees with the second derivative of \(H_j\).
The quadratic chord error makes the integrand bounded at zero,
and its exponential denominator controls infinity. These facts
also justify differentiation on \(0\le s\le1\).
The numerator of (4) is strictly negative for every \(x>0\)
on a nonconstant row at an interior \(s\). In the difference
\(H_N-\alpha H_n\), replace the denominator kernel by \(D_n(x)/x\).
Subtracting exponentials in the integral with kernel \(e^{-\delta x}/x\)
gives \(-L_\delta(s)\), by integrating \(e^{-t x}\) first in \(x\)
and then between the two exponents in \(t\). Applying (3)
proves (1).

To prove (2), push \(e^{-a_i x}D_n(x)\,dx/x\) on each unequal
row to \(z=-d_i x\), and sum the resulting measures to form
\(\nu\) on \(\mathbb R\setminus\{0\}\). Then
\[
 H_N(s)-\alpha H_n(s)=\int
       [e^{sz}-1-s(e^z-1)]\,d\nu(z).
\]
The upper bound in (3) proves \(\int |z|\,d\nu<\infty\) and
\(\int_{|z|>1} e^z\,d\nu<\infty\), including the private orientations.
Take independent compound Poisson sums of jumps on the dyadic
annuli that partition \(\mathbb R\setminus\{0\}\). The expected
total absolute size of their jumps is \(\int |z|\,d\nu\), so their
sum converges absolutely almost surely. Put \(Z\) equal to that
sum minus \(\int(e^z-1)\,d\nu\), and take \(Y=e^Z\). For \(0\le
s\le1\), the positive jumps give convergence of the exponential
moments by monotone convergence and the negative jumps by bounded
convergence. Independence gives the displayed exponent for
\(\mathbb E e^{sZ}\), and at \(s=1\) this equals one.

## A simultaneous bound in terms of both adjacent coefficients

For \(1\le k\le n\) use \(s=k/N\), and write
\[
 m_i^\pm = \underset{t\in[(k-1)/n,k/n]}{\max/\min}\, X_i(t),
 \qquad
 C_-=\sum_i\frac{d_i^2}{m_i^++1/n},\qquad
 C_+=\sum_i\frac{d_i^2}{m_i^-+1/n}
                      \left(1+\frac{1}{1+n m_i^-}\right).
\]
The maximum defines \(m_i^+\) and the minimum \(m_i^-\).
The comparison of the trigamma series with its integrals gives
\(nC_-\le H_n''(t)\le n C_+\) on this interval. The quadratic
chord bounds around \(s\), combined with (1), then give
\[
 \begin{aligned}
 &\frac{k(N-k)}{2n^2N}C_- + \frac{L_{1/n}(s)}{2n}\\
 &\qquad\le
 \log\frac{F_{n,k-1}^{k/n}\,F_{n,k}^{(N-k)/n}}{F_{N,k}}\\
 &\qquad\le
 \frac{k(N-k)}{2n^2N}C_+
                 + \frac{L_{1/(3N)}(s)}{2n}.
 \end{aligned}
 \tag{5}
\]
Indeed the distances from \(s\) to the ends of the interval
are \((N-k)/(nN)\) and \(k/(nN)\). The chord error has
one half of their product times the respective second-derivative
bound. Multiplication by \(\alpha=N/n\) gives the coefficients
in (5) and the indicated exponents of the two lower multipliers.

Formula (5) keeps both neighboring coefficients in one inequality.
It complements the
[fixed-distance endpoint attenuation](MULTIVARIATE_ORDINARY_ADJACENT_COEFFICIENT_ATTENUATION.md).
Neither (2) nor (5) alone supplies a comparison of root arguments.

## Fractional transport can break circular root location

For a nontrivial mirrored profile the difference \(H_N-\alpha H_n\)
is invariant under \(s\mapsto1-s\). For \(0\le k\le N\), set
\[
 A_k=H_N(k/N)-\alpha H_n(k/N),\qquad
 D_\tau(t)=\sum_{k=0}^N\binom Nk e^{\tau A_k}t^k.
\]
For every \(N\ge4\), \(D_\tau\) has off-circle roots for all
sufficiently small \(\tau>0\), even though \(D_0=(1+t)^N\) is
circle-rooted. Thus the actual fractional transport does not
act as a circular-root-preserving semigroup on degree-\(N\)
polynomials.

To see this, write \(L(t)=\sum_k\binom Nk A_k t^k\), the first
variation of \(D_\tau\) at zero. For even \(N\), the jump integral
and the binomial formula give
\[
 L(-1)=\int (1-e^{z/N})^N\,d\nu(z)>0.
\]
The constant and linear terms in the jump integrand cancel under
the alternating sum. At least one unequal row supplies mass
away from \(z=0\). Thus, under \(t=-1+\tau^{1/N}\xi\),
\[
 \tau^{-1}D_\tau(-1+\tau^{1/N}\xi)\longrightarrow
                      \xi^N+L(-1)
\]
uniformly on compact \(\xi\)-sets. For odd \(N\), reciprocity
gives a structural factor \(1+t\), so divide it out first.
The first variation of the quotient at \(-1\) is
\[
 L'(-1)=N\int e^{z/N}(1-e^{z/N})^{N-1}\,d\nu(z)>0 .
\]
Here differentiation of the binomial sum supplies the displayed
factor and still cancels the constant and linear terms. Applying
\(t=-1+\tau^{1/(N-1)}\xi\) to the quotient gives the limit
\(\xi^{N-1}+L'(-1)\) after dividing by \(\tau\). These limits
have even degree at least four and have simple roots with
nonzero real parts. Rouché's theorem gives corresponding scaled
roots for small \(\tau>0\). Roots on the unit circle tending to
\(-1\) at either indicated scale must have purely imaginary scaled
limits, by expanding \(|-1+\tau^{1/m}\xi|^2=1\). The
transported polynomial therefore has roots off that circle.

This concerns the action of the moment quotient on a test
polynomial. It does not assert a shared root of any ordinary
factorial moment pair.

