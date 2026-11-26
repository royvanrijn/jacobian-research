# All-degree interlacing with one private unit pair and small shared ratios

Consider the ordinary mirrored two-monomial profile
with one private unit pair \((1,0),(0,1)\) and
additional mirrored shared pairs
\((\ell_i+d_i,\ell_i),(\ell_i,\ell_i+d_i)\), where
\(\ell_i,d_i\ge1\). Set
\[
                            v=\sum_i d_i^2/\ell_i
\]
with one term per additional pair.
Equal rows contribute only a scalar.

**Theorem.** If \(v\le2\), the consecutive ordinary
moment polynomials strictly interlace on the unit
circle in every degree. This proves SFC(2) for these
profiles over every characteristic-zero field.
The weak restriction here is on the shared
ratios, while the private unit pair is exact.

The proof converts the matched likelihood-ratio
comparison of the even log-ratio laws into an
angular inequality at the actual lower roots.
A uniform fourth-moment ratio controls the lost
phase margin under the matched dilation.

## The phase margin from the matched dilation

Use \(R_j=A_j/B_j\) in the
[one-private angle formula](MULTIVARIATE_PRIVATE_UNIT_MIRRORED_ROOT_MONOTONICITY.md#an-even-log-ratio-law-and-a-two-term-root-equation).
In particular its positive-half angle boxes are
the
[geometric resolvent boxes](MULTIVARIATE_ONE_PRIVATE_MIRRORED_ANGLE_BOXES.md#the-exact-resolvent-and-angle-boxes).
If the shared core is empty, the moment polynomials
are geometric sums and already strictly alternate,
so assume the core is nonempty.

Write \(N=n+1\), fix \(0<\theta\le\pi/2\) and
\(c=\cos\theta\), and let \(\rho_n^+\) denote
the folded even law of the shared log ratio.
The
[matched folded likelihood-ratio sandwich](MULTIVARIATE_PRIVATE_UNIT_MIRRORED_ROOT_MONOTONICITY.md#the-folded-likelihood-ratio-comparison)
says that \((n/N)\rho_n^+\) is below \(\rho_N^+\).
The ratio defining \(R_N\) is the expectation
of the increasing function
\(\cosh((N/2+1)s)/\cosh(Ns/2)\) under the common
positive weight \(\cosh(Ns/2)/(\cosh s-c)\).
Therefore
\[
 \begin{split}
  R_N(\theta)&\ge Q(n/N),\\
  Q(t)&=
   \frac{\displaystyle\int
        \cosh((n/2+t)s)(\cosh(ts)-c)^{-1}\,d\rho_n^+(s)}
        {\displaystyle\int
        \cosh(ns/2)(\cosh(ts)-c)^{-1}\,d\rho_n^+(s)},
                  \qquad Q(1)=R_n(\theta).
 \end{split}                                                   \tag{1}
\]
The functions in (1) have finite integrals since
the shared exponents are positive.

Here is the differential estimate. Put
\(D_t(s)=\cosh(ts)-c\) and
\(f_t(s)=\cosh((n/2+t)s)-\cosh(ns/2)\).
Then \(Q(t)-1\) is the quotient of the integrals
of \(f_t/D_t\) and \(\cosh(ns/2)/D_t\).
The ratio of these integrands,
\[
  f_t/\cosh(ns/2)=\cosh(ts)-1+\tanh(ns/2)\sinh(ts),
\]
increases in \(s>0\).
Also \(\partial_t\log D_t=s\sinh(ts)/D_t\)
increases in \(s\). To see the latter, the
derivative in \(x=ts\) of
\(x\sinh x/(\cosh x-c)\) has numerator
\(\sinh x(\cosh x-c)+x(1-c\cosh x)\).
It is smallest at \(c=1\), where it equals
\((\cosh x-1)(\sinh x-x)\ge0\).
Thus the denominator-weight terms in the logarithmic
derivative of \(Q-1\) contribute nonpositively,
by monotone likelihood-ratio comparison between
the two positive integrands.

Factor
\(f_t(s)=2\sinh(ts/2)\sinh((n+t)s/2)\).
It follows that
\[
     \partial_t\log(Q(t)-1)\le
                   \frac1t+\frac1{n+t}+\mathbb E_{\nu_t}s,  \tag{2}
\]
where \(\nu_t\) is the probability with weight \(f_t/D_t\)
against \(\rho_n^+\). Here
\(\partial_t\log f_t\le s+1/t+1/(n+t)\),
using \(\coth x\le1+1/x\) for \(x>0\).
Differentiation is justified on \([n/N,1]\):
at infinity the numerator weight grows like
\(e^{ns/2}\), and the extra factor \(s\) from
differentiation is integrable against the endpoint
law. That law is the folded log ratio under the
ordinary tilt by the \(n\)-th power of a shared
monomial, and has an exponential moment at
positive rate because each smaller exponent is
positive.

Let \(H\) be the full, unfolded log ratio under
that endpoint tilt. Its folded density is
proportional to \(\cosh(ns/2)d\rho_n^+(s)\).
Relative to that folded law the \(\nu_t\)
weight is \(K_t=f_t/(\cosh(ns/2)D_t)\).
With \(z=\tanh(ts/2)\),
\[
 \frac{K_t(s)}{s^2}
   =\frac{2(z/s)^2+2(z/s)(\tanh(ns/2)/s)}
                  {(1-c)+(1+c)z^2}
\]
is decreasing: both \(\tanh\)-ratios decrease
and the denominator increases. Comparing with
the \(s^2\)-weighted folded law therefore gives
\[
                    \mathbb E_{\nu_t}s
                \le\frac{\mathbb E|H|^3}{\mathbb E H^2}
                \le\sqrt{\frac{\mathbb E H^4}{\mathbb E H^2}}. \tag{3}
\]

## Exact Gamma bounds for the endpoint tilt

The variables in
\(H=\sum_i d_i(\log Y_i-\log X_i)\)
are independent with shapes
\(n\ell_i+1\) for \(Y_i\) and
\(n(\ell_i+d_i)+1\) for \(X_i\).
Write \(a=|\mathbb EH|\), \(\sigma^2=\operatorname{Var}H\)
and \(\kappa_r\) for the cumulants of order \(r\).
The polygamma series gives
\[
 a\le\sum_i d_i\log((\ell_i+d_i)/\ell_i)\le v,\qquad
 \sigma^2\le 2v/n,\qquad
 |\kappa_3|\le 2(v/n)\sigma^2,\qquad
 0\le \kappa_4\le 6(v/n)^2\sigma^2 .
\]
Indeed \(\psi'(1+x)\le1/x\) gives the first
two bounds by integration and variance addition.
The series gives
\(|\psi''(z)|\le2\psi'(z)/z\) and
\(\psi'''(z)\le6\psi'(z)/z^2\), and
\(d_i/(n\ell_i)\le v/n\) because \(d_i\ge1\).

Using \(\mathbb EH^2=a^2+\sigma^2\), the
fourth moment expansion gives
\[
 \begin{split}
 \frac{\mathbb EH^4}{\mathbb EH^2}
 &\le a^2+5\sigma^2+4(v/n)\sigma+6(v/n)^2\\
 &\le v^2+10v/n+4(v/n)\sqrt{2v/n}+6v^2/n^2 .
 \end{split}
\]
The first line bounds
\(a^4+6a^2\sigma^2+3\sigma^4\)
by \((a^2+5\sigma^2)(a^2+\sigma^2)\);
the cubic contribution uses
\(a\sigma/(a^2+\sigma^2)\le1/2\).
When \(v\le2\), the square root in (3) is
at most \(2+5/n+4/n^{3/2}+6/n^2\) by concavity
of the square root. Since \(1/t\le1+1/n\)
and \(1/(n+t)\le1/n\), (2) gives
\[
             \partial_t\log(Q(t)-1)\le
                   C_n:=3+7/n+4/n^{3/2}+6/n^2.
\]
Integrating the differential bound from \(n/N\) to one yields
\[
            R_N(\theta)\ \ge\
                  1+e^{-C_n/N}(R_n(\theta)-1).               \tag{4}
\]

## Strict adjacent comparison

Let \(\theta\) be the \(j\)-th lower positive-half
angle. If \(\theta>\pi/2\), the boxes give
\(4j>n+1\), hence \(4j\ge n+2\), which already
places the upper \(j\)-th angle before it.
In the other case \(\pi/n<\theta\le\pi/2\).
The order-\(N\) equation evaluated at this
lower root has right side
\[
      \cos\theta-\tan(N\theta/2)\sin\theta
           = c+\frac{(1+c)(R_n(\theta)-1)}{R_n(\theta)+1}.
                                                               \tag{5}
\]
This follows by adding \(\theta/2\) to the
lower phase equation.
If the lower angle is already to the right
of the \(j\)-th upper box the desired inequality
holds. Otherwise it is in that box, since it
exceeds the left endpoint.
The right side of the upper root equation is
strictly decreasing there and \(R_N\) increases,
as in the fixed-degree monotonicity proof.
Thus \(R_N(\theta)\) strictly above (5)
places the upper angle before the lower.

Write \(r=R_n-1>0\), \(\lambda=1-e^{-C_n/N}\)
and \(\delta=1-c\).
Subtract (5) from the lower bound (4) and
multiply by \(2+r\). The result is
\[
                       (1-\lambda)r^2+2(\delta-\lambda)r
                                                      +2\delta.
\]
Its discriminant is \(4(\lambda^2-\sin^2\theta)\).
The quadratic is strictly positive for every \(r\) whenever
\[
                             1-e^{-C_n/N}<\sin(\pi/n).
                                                               \tag{6}
\]
For \(n\ge6\), \(\sin(\pi/n)\ge3/n\) by
concavity of sine up to \(\pi/6\).
The cubic lower bound for
\(-\log(1-3/n)\) shows that the right side
of \(C_n<-N\log(1-\sin(\pi/n))\) exceeds
\[
                      3+\frac{15}{2n}+\frac{27}{2n^2}
                                         +\frac9{n^3}.
\]
Its difference from \(C_n\), multiplied by \(n^2\),
is positive. With \(u=\sqrt n-4\), that difference
is \(u^2/2-1/2+9/(\sqrt n)^2\), bounded below by
\[
                        u^2/2-9u/32+1/16>0
\]
by the tangent line for \(9/x^2\) at four.
For \(n=4\), \(C_n/N=9/8\) and
\(e^{9/8}<10/3\) gives \(1-e^{-C_n/N}<7/10<\sin(\pi/4)\).
For \(n=5\), \(C_n/N<21/25\) and
\(e^{21/25}<50/21\) gives
\(1-e^{-C_n/N}<29/50<\sin(\pi/5)\).
The exponential upper bounds follow by the
Taylor sum through degree six with its positive
geometric tail bound. This proves (6).

The other side of alternation follows from the
boxes without a phase comparison.
Starting degrees \(n\le3\) follow from the
[direct low-degree identities](MULTIVARIATE_ONE_PRIVATE_MIRRORED_ANGLE_BOXES.md#consecutive-order-comparison).
Thus every adjacent pair strictly alternates.
The rational Bézout identity gives the factorial
assertion over characteristic-zero fields.
