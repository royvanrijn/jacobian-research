# A polynomial diffuse tolerance for ordinary factorial angular separation

Let \(a,b\) be distinct nonnegative integer exponent columns.
Remove equal rows and put
\[
 d_i=b_i-a_i,\quad L_i=\max(a_i,b_i),\quad
 S=\sum_i d_i^2/L_i,\quad \eta=\max_i |d_i|/L_i .
\]
For every terminal \(N=n+1\ge2\), the ordinary factorial moment
polynomials of orders \(n,N\) have disjoint root arguments under
independent positive rescalings if
\[
                      S\ge13,\qquad S\eta N\le10^{-6}.       \tag{1}
\]
In particular the two polynomials are coprime over every
characteristic-zero field. Together with the
[absolute high-curvature theorem](MULTIVARIATE_EFFECTIVE_ABSOLUTE_CURVATURE_BOUND.md),
the new range is \(13\le S<64\) under the diffuse bound in (1).

The proof starts with the two nonvanishing analytic halves of
the [Gaussian reference](MULTIVARIATE_DIFFUSE_EXPLICIT_S13_INTERLACING.md).
An analytic fractional-index multiplier transfers to each half
separately with a relative error. Thus cancellation between the
two halves remains available when comparing the actual root phases.

## The analytic fractional-index multiplier

Normalize the endpoints at each order by dividing each coefficient
by its positive endpoint chord. For complex \(v\) put \(x_i(v)=a_i+d_i v\)
and
\[
 \begin{split}
 E_j(v)=\sum_i\big[&
 \log\Gamma(jx_i(v)+1)
 -(1-v)\log\Gamma(ja_i+1)-v\log\Gamma(jb_i+1)\big]\\
 &+\frac{jS}{2}v(1-v).
 \end{split}
\]
The normalized actual polynomial is
\[
 P_j(z)=\sum_{k=0}^j ([z^k]B_j(z;S))\,e^{E_j(k/j)}z^k ,
\]
where \(B_j(z;S)=\sum_k\binom jk e^{-S k(j-k)/(2j)}z^k\).
The hypothesis implies \(\eta<1/60\). On \(|v|\le5\),
\(|x_i(v)-L_i|\le6|d_i|\) and
\(\Re x_i(v)\ge (1-6\eta)L_i\ge9L_i/10\).
Thus all Gamma logarithms use the branch analytic in the positive
real half-plane and real on the positive real axis.

For \(\Re q>0\), compare the trigamma series for \(\psi'(q+1)\)
with \(\int_0^\infty(q+t)^{-2}dt=1/q\). On each unit interval
the difference between the right-endpoint value and the integral
is bounded by the integral of \(2/(\Re q+t)^3\) on that interval.
Summing gives \(\left|\psi'(q+1)-1/q\right|\le(\Re q)^{-2}\).
[DLMF, equation 5.15.1](https://dlmf.nist.gov/5.15.E1)

Because \(d_i\) is a nonzero integer, \(1/L_i\le\eta\).
Differentiate \(E_j\) twice and compare the series first to
\(1/(jx_i(v))\) and then to \(1/(jL_i)\). This gives
\[
 \sup_{|v|\le5}|E_j''(v)|
 \le \frac{6jS\eta}{1-6\eta}+
                \frac{S\eta}{(1-6\eta)^2}\le10jS\eta .
\]
Since \(E_j(0)=E_j(1)=0\), integration over \([0,1]\) bounds
\(|E_j'(0)|\) by half this supremum. A second integration,
along the segment from zero to \(v\), gives
\(\sup_{|v|\le5}|E_j(v)|\le150jS\eta\).
For \(j=n,N\), the two functions \(e^{E_j(v)}-1\) and
\(e^{E_j(1-v)}-1\) therefore have modulus at most
\(300jS\eta\) on \(|v|\le4\). Their coefficients at degree \(m\)
are bounded by \(300jS\eta\,4^{-m}\).

## Relative transport on the two halves

Write
\[
 B_j(z;S)=I_j(z)+z^j I_j(1/z),\qquad H_j=\log I_j .
\]
For \(j\ge22\), use the polynomial half-Gaussian integral
from the \(s\ge13\) proof. Its coefficients are nonnegative
by expanding the integral on the real half-line.
On \(|z|\le2\) that proof gives
\[
 I_j=e^{jF_*}(1+\omega_j),\qquad
 |F_*|<1/250,\qquad |\omega_j|<r_0=21/100.
\]
Here \(F_*\) is common to both orders with the same \(S\).
Let \(C_0=r_0+r_0^2/(2(1-r_0))\).
Use the logarithmic series to define \(H_j\), so
\(\|H_j-jF_*\|_{|z|\le2}\le C_0\).

For \(j\le21\) use the half polynomial of the terms with \(k<j/2\),
with half the middle coefficient included at even \(j\).
Its coefficients are nonnegative. For \(S\ge13\) we have
\(\|I_j-1\|_{|z|\le2}<1/10\), by the finite rational bound below.
Thus \(\|H_j\|_{|z|\le2}\le1/9\), by the logarithmic series.
All the logarithms are real on the positive radius and have
conjugation symmetry. In particular, for \(j\ge2\),
\[
 \|H_j\|_{|z|\le2}\le j/18,\qquad
 \|H_N'-H_n'\|_{|z|\le1}<12/25.                            \tag{2}
\]
For the derivative difference use Cauchy's bound with radius one.
When both degrees are at least \(22\), its upper bound is
\(1/250+2C_0<12/25\). With the degrees \(21,22\) it is at
most \(1/9+22/250+C_0<12/25\).
For both degrees at most \(21\) it is at most \(2/9\).

Write \(D=z\,d/dz\). For every \(m\ge1\) and \(|z|\le3/2\),
\[
             \left|\frac{(D/j)^m I_j(z)}{I_j(z)}\right|
                              \le2(5/4)^m.                 \tag{3}
\]
If \(m\le j/4\), apply Cauchy's derivative inequality in \(h\)
to \(I_j(ze^h)/I_j(z)\) on \(|h|\le m/j\).
There \(|ze^h|\le(3/2)e^{1/4}<39/20\).
For \(j\ge22\), Cauchy on the \(F_*\) disk gives
\(\|F_*'\|_{|z|\le39/20}\le2/25\).
Since \(|ze^h-z|\le(39/20)|h|\), the quotient has modulus
at most \(2 e^{m/5}\) on the \(h\) disk.
For \(j\le21\), the quotient has modulus less than \(2\).
The derivative estimate then gives (3) using
\(m!/m^m\le1\) and \(e^{1/5}\le5/4\).
If \(m>j/4\), coefficient positivity and the degree bound give
\(\left|(D/j)^m I_j(z)\right|\le I_j(|z|)\).
The quotient \(I_j(|z|)/|I_j(z)|\) is at most \(2e^{j/125}\)
by the same half estimates, and \(j<4m\) gives (3).

Define positive-coefficient polynomials \(\widetilde I_{j,+}\)
and \(\widetilde I_{j,-}\) by multiplying the coefficients of
\(I_j\) at \(z^k\) by \(e^{E_j(k/j)}\) and \(e^{E_j(1-k/j)}\),
respectively. The actual polynomial satisfies exactly
\[
        P_j(z)=\widetilde I_{j,+}(z)+z^j\widetilde I_{j,-}(1/z).
\]
Expand the analytic multipliers at zero and substitute \(D/j\).
This acts on \(z^k\) by \(k/j\), so the convergent series
reproduces the stated coefficients. Its operator action on the
half is bounded by (3). On \(|z|\le3/2\) the relative half
errors are bounded by
\[
 300jS\eta\left(1+2\sum_{m\ge1}(5/16)^m\right)<600jS\eta .
\]
In particular the actual halves are nonvanishing there.
Their logarithms have the form \(H_j+e_{j,\pm}\), with
\(\|e_{j,\pm}\|_{|z|\le3/2}\le\kappa j\) for
\(\kappa=1200S\eta\), using the series and (1).

## Angular comparison of the ordinary roots

All roots of each actual polynomial of degree \(j\ge2\)
lie in the annulus \(|\log|z||<1/8\).
On the outer circle the reciprocal term in the decomposition
strictly dominates because \(1/8>2(1/18+\kappa)\).
That reciprocal polynomial has all its zeros inside the circle,
since its unreversed half is nonvanishing on \(|z|\le3/2\).
Rouché counts all roots inside. The same argument after
reciprocal reversal proves the inner bound.

Take any root \(z=e^{q+i\theta}\) with \(0<\theta\le\pi\).
Positive coefficients exclude \(\theta=0\). Put
\[
\begin{split}
 J_j^0(q,\theta)&=H_j(e^{q+i\theta})-H_j(e^{-q-i\theta}),\\
 \delta J_j(q,\theta)&=e_{j,+}(e^{q+i\theta})-e_{j,-}(e^{-q-i\theta}).
\end{split}
\]
Then \(j(q+i\theta)=J_j^0(q,\theta)+\delta J_j(q,\theta)
+(2k_j+1)i\pi\) for some integer \(k_j\).
The baseline has \(\Re J_j^0(0,\theta)=0\) by conjugation.
Throughout \(|q|\le1/8\) the evaluation points have modulus
less than \(6/5\). Cauchy applied to (2) on \(|z|\le2\) gives
\[
 |\partial_q J_j^0|\le j/6,\qquad
 |\partial_q\partial_\theta J_j^0|\le2j/3.
\]
Indeed, if \(M_j=\|H_j\|_{|z|\le2}\), its first and second
derivatives on \(|z|\le6/5\) are bounded by \(5M_j/4\) and
\(25M_j/8\). The two derivative bounds for \(J_j^0\) are
then \(3M_j\) and \(12M_j\).
The real part of the root equation gives \(|q|\le12\kappa/5\).

Cauchy on the actual half-error disk \(|z|\le3/2\) gives
\(|\partial_\theta\delta J_j|\le8\kappa j\).
All the imaginary parts vanish when \(\theta=0\) for real \(q\).
Integrating in \(\theta\), and using the mixed derivative for
the change in \(q\), gives
\[
 \left|\Im[J_j^0(q,\theta)+\delta J_j(q,\theta)-J_j^0(0,\theta)]\right|
                      < 10\kappa j\theta .
\]
If the two orders had a common root argument, conjugate
the roots to the same \(0<\theta\le\pi\).
At \(q=0\) the difference of the two baseline imaginary
parts is \(2\Im(H_N-H_n)(e^{i\theta})\), of modulus at most
\((24/25)\theta\) by (2) and reality at \(\theta=0\).
The two ordinary errors contribute at most
\(20\kappa N\theta\le(3/125)\theta\).
Subtracting the imaginary root equations would place an
integer multiple of \(2\pi\) between \(2\theta/125>0\)
and \(248\theta/125<2\pi\), which is impossible.

This proves the result for \(n\ge2\).
At \(n=1\) the only argument of the first moment is negative
real, and the second moment is positive on the real line by
the integral of the square of a nonzero binomial.
The rational coefficient Bézout identity transfers coprimality
to every characteristic-zero field.

## Exact low-degree half bound

For \(S=13\), the rational calculation below bounds the nonconstant
part of every low-degree half on \(|z|\le2\).
All its coefficients decrease as \(S\) increases.
Use the reciprocal of a positive Taylor partial sum to bound
each \(e^{-x}\) above.

~~~python
from fractions import Fraction as Q
from math import comb, factorial

for j in range(1, 22):
    upper = Q(0)
    for k in range(1, j//2+1):
        x = Q(13*k*(j-k), 2*j)
        exp_upper = 1/sum(x**m/factorial(m) for m in range(33))
        midpoint_weight = Q(1,2) if 2*k == j else 1
        upper += midpoint_weight * comb(j,k) * 2**k * exp_upper
    assert upper < Q(1,10)
~~~

