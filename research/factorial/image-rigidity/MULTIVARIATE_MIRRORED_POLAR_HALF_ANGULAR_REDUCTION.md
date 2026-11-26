# A polar-half criterion for consecutive mirrored angular interlacing

For an ordinary mirrored two-monomial profile, write \(F_m\)
for the degree-\(m\) endpoint-normalized polynomial. Its
[simple unit-circle root property](MULTIVARIATE_MIRRORED_ROOT_CIRCLE.md#mirrored-profiles-at-every-fixed-order)
holds in every degree. The argument below applies to any
real reciprocal \(F_m\) with \(F_m(0)=1\), positive
coefficients, and simple unit-circle roots.
It reduces consecutive interlacing to an exact relative
sign bound on one linked coefficient multiplier at the
roots of the lower array.

## A stable half and its increasing phase

Define the polar half
\[
 H_m(z)=F_m(z)-\frac{z}{m}F_m'(z)
       =\frac{z^{m-1}}{m}F_m'(1/z).
\]
The second equality follows by differentiating reciprocity.
It gives \(F_m(z)=H_m(z)+z^m H_m(1/z)\).
The polynomial \(H_m\) is nonvanishing on the closed unit
disk. In fact, if \(\zeta\) runs over the roots of \(F_m\),
then
\[
              \frac{zF_m'(z)}{F_m(z)}
                         =\sum_\zeta\frac{1}{1-\zeta/z}.
\]
For \(|z|>1\) and for points of \(|z|=1\) that are not
roots, each summand has positive real part. At a root
\(F_m'\ne0\) by simplicity. Thus every zero of \(F_m'\)
lies strictly inside the circle, and the reciprocal
polynomial \(H_m\) has no zero on or inside it.

Take \(\log H_m(0)=0\) and define
\[
       \phi_m(\theta)=\Im\log H_m(e^{i\theta}),\qquad
       \psi_m(\theta)=m\theta/2-\phi_m(\theta).
\]
The centered polynomial \(e^{-im\theta/2}F_m(e^{i\theta})\)
is \(2|H_m(e^{i\theta})|\cos\psi_m(\theta)\).
For each zero \(\xi\) of \(H_m\), the contribution of its
factor to \(zH_m'(z)/H_m(z)\) is \(1-(1-z/\xi)^{-1}\).
On \(|z|=1\) its real part is less than \(1/2\), since
\(|z/\xi|<1\).
It follows that \(\psi_m'>1/2\) for \(m\ge2\), with
\(\psi_1'=1/2\) when \(H_1=1\).
The analytic log gives \(\psi_m(0)=0\) and
\(\psi_m(2\pi)=m\pi\).
This supplies increasing phases at every curvature.

## An exact lower-root test

Let \(N=n+1\), and denote the lower root arguments in
the ordering cut at \(1\) by
\(0<\theta_1<\cdots<\theta_n<2\pi\).
Then \(\psi_n(\theta_i)=(i-1/2)\pi\).
Strict interlacing holds exactly when
\[
          \theta_i/2-\pi <
                  \phi_N(\theta_i)-\phi_n(\theta_i)
                                < \theta_i/2
                             \qquad (1\le i\le n).              \tag{1}
\]
Indeed the correct upper interval at \(\theta_i\) is
\((i-1/2)\pi<\psi_N(\theta_i)<(i+1/2)\pi\).

There is an equivalent test without choosing a numerical
lift for either argument. Set \(z_i=e^{i\theta_i}\) and
\(T_i=H_N(z_i)/H_n(z_i)\).
For the real centered functions
\(f_m(\theta)=e^{-im\theta/2}F_m(e^{i\theta})\), we have
\[
         e^{-in\theta_i/2}H_n(z_i)=\frac{i}{n}f_n'(\theta_i),
         \qquad
     \frac{f_N(\theta_i)}{f_n'(\theta_i)}
                     =-\frac2n\Im\bigl(e^{-i\theta_i/2}T_i\bigr).
                                                               \tag{2}
\]
Thus strict interlacing is equivalent to
\(\Im(e^{-i\theta_i/2}T_i)<0\) for all lower roots.
Necessity follows from the alternating signs of
\(f_N(\theta_i)\) and \(f_n'(\theta_i)\).
For sufficiency, the stated signs force an upper root
in each interval between successive lower arguments,
as well as in both end intervals at the positive cut,
since the upper endpoint signs are \(+\) at zero and
\((-1)^N\) at \(2\pi\). These account for all \(N\) upper
roots. In particular (2) is equivalent to the global
collection of tests in (1).

## The ordinary tilt and exact multiplier

For positive monomials \(U,V\) let
\(G_m=(\mathbb E[V^m]/\mathbb E[U^m])^{1/m}\), so that the
endpoint-normalized \(F_m\) evaluates the physical moment
at \(z/G_m\). Its polar half is
\[
       H_m(z)=\frac{\mathbb E[U(U+(z/G_m)V)^{m-1}]}{\mathbb E[U^m]}.
\]
Under the fixed tilt \(d\mathbb P_U=U\,d\mathbb P/\mathbb E U\),
this is the degree-\(m-1\) moment polynomial evaluated at
\(z/G_m\) and normalized at zero with respect to
\(\mathbb P_U\). The two polar halves use the same tilt.
For ordinary mirrored rows \(U,V\) are exchangeable, so
\(G_m=1\) at every order.

For a direct computation from the lower coefficients,
write \(F_n(z)=\sum_{k=0}^n c_{n,k}z^k\).
Let \(A_m(s)\) be the endpoint-normalized factorial
log symbol of the
[holomorphic matched log transport](MULTIVARIATE_HOLOMORPHIC_MATCHED_LOG_TRANSPORT.md#a-paired-multiplier-on-the-same-degree).
In its notation the exact lower-to-upper polar multiplier is
\[
       H_N(z)=\sum_{k=0}^n c_{n,k} e^{g_0(k/n)}z^k,\qquad
       g_0(s)=N A_N(ns/N)-n A_n(s).
\]
This follows from \((1-k/N)\binom Nk=\binom nk\).
If a common normalized Gaussian strength \(\tau\) is
included, replace \(A_m(s)\) by
\(A_m(s)-\tau s(1-s)/2\), or equivalently multiply each
polar multiplier by \(e^{-\tau n s^2/(2N)}\).
The factorial transport controls its own \(A_m\) term as
in the cited formula.

Since \(F_n(z_i)=0\) and \(H_n(z_i)=-z_iF_n'(z_i)/n\),
the sign test (2) uses the exact quotient
\[
       T_i=\frac{-n}{z_iF_n'(z_i)}
          \sum_{k=0}^n c_{n,k}
                     \bigl(e^{g_0(k/n)}-C\bigr)z_i^k
                \qquad(C\in\mathbb R\ \text{arbitrary}).
\]
For example a sufficient relative estimate is
\[
   \frac{n}{|F_n'(z_i)|}
      \left|\sum_{k=0}^n c_{n,k}
                \bigl(e^{g_0(k/n)}-C+\lambda_i k/n\bigr)z_i^k\right|
                  <\lambda_i \sin(\theta_i/2),\qquad\lambda_i>0.
\]
The affine baseline contributes exactly \(\lambda_i\) to
\(T_i\). The estimate therefore puts
\(e^{-i\theta_i/2}T_i\) in the lower half-plane.

The stable polar half is available even when a truncated
symmetric coefficient half has zeros inside the disk.
The remaining obligation below the proved angular
ranges is the relative sign or remainder estimate at the
lower roots. The matched symbol and its derivatives alone
do not supply this estimate.
