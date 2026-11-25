# Angular interlacing under Euler operators

The [Euler contiguity relations](MULTIVARIATE_LINKED_GAMMA_ANGLE_COMPACTIFICATION.md#euler-operator-relations-for-integer-exponent-rows)
relate consecutive moment polynomials through forward Euler multipliers
and their inverses. Angular interlacing of general positive,
binomial-normalized log-convex arrays does not suffice to transfer
interlacing through these operations. Even strict interlacing of
reciprocal inputs can be lost. The examples here concern the operator
property and leave open the angular comparison for linked Gamma profiles.

Put \(D=t\,d/dt\), \(Q_c(t)=t^2+ct+1\), and \(T_\eta=1+\eta D\).
If \(z\ne0\) is a simple root of \(P\), the corresponding root
angles of \(T_\eta P\) and \(T_\eta^{-1}P\), for sufficiently
small real \(\eta\), are
\[
 \arg z+\eta^2\,\Im\frac{zP''(z)}{2P'(z)}+O(\eta^3),
 \qquad
 \arg z-\eta^2\,\Im\frac{zP''(z)}{2P'(z)}+O(\eta^3),        \tag{1}
\]
respectively. The branches follow by the implicit-function theorem.
Their root expansions through the quadratic term are
\(z-\eta z+\eta^2(z+z^2P''/(2P'))\) and
\(z+\eta z-\eta^2z^2P''/(2P')\). To obtain the second expansion,
use \(T_\eta^{-1}P=P-\eta DP+\eta^2D^2P+O(\eta^3)\).
If \(z=e^{\mathrm i\theta}\) and \(P\) has all roots on the unit
circle, the imaginary quantity in (1) is
\[
 \Im\frac{zP''(z)}{2P'(z)}
     =\tfrac12\sum_{w\ne z}\cot\!\left(\frac{\arg w-\theta}{2}\right),
\]
where the sum runs over the other roots with multiplicity. A factor
\(Q_c\) with \(c\ne -2\cos\theta\) contributes
\(-2\sin\theta/(-2\cos\theta-c)\).

## A common forward multiplier

Take
\[
 P=Q_{9/5}Q_{7/5},\qquad
 R_\alpha=(1+t)Q_{6/5}Q_\alpha,\qquad
 17/10\le\alpha<9/5.
\]
The input arguments strictly interlace, since the upper root
angle of \(Q_c\) strictly increases with \(c\). At the limiting
endpoint \(\alpha=9/5\), there is a shared largest upper angle
\(\theta\). In (1), the quintic's quadratic angle displacement
minus the quartic's has positive coefficient
\[
 2\sin\theta\left(\frac1{9/5-7/5}-\frac1{9/5-6/5}\right)
   +\tfrac12\cot((\pi-\theta)/2)>0.
\]
The terms from the common conjugate factor cancel. For every
sufficiently small positive \(\eta\) the transformed upper angle
of \(R_{9/5}\) thus exceeds that of \(P\), while the transformed upper
angle of \(R_{17/10}\) remains below it by continuity. The simple
endpoint roots continue uniformly for small \(\eta\) over this
compact \(\alpha\)-interval. The intermediate value theorem gives
an \(\alpha\) strictly below \(9/5\) at which \(T_\eta P\) and
\(T_\eta R_\alpha\) have equal root arguments.

The input quartic has binomial-normalized outer and inner
coefficients \(u=.8\) and \(v=113/150\), so \(u^2<v<u\).
The corresponding quintic values satisfy
\(u\in[.78,.80]\), \(v\in[.694,.716]\), which gives the
same strict inequalities uniformly. The transformed arrays
retain positivity and strict log convexity for small enough
\(\eta\) by continuity on the compact parameter interval.

## Inverse multipliers at arbitrarily large orders

Fix \(m\ge2\), put
\(\lambda_j=2(m-j)+3\) for \(1\le j<m\), and set
\(\lambda_m=1/(16m)\). For small \(\tau>0\) let
\[
 P(t)=\prod_{j=1}^m Q_{2-\tau\lambda_j}(t),\qquad
 R_\nu(t)=(1+t)Q_{2-\tau\nu}(t)
                       \prod_{j=1}^{m-1}Q_{2-\tau(\lambda_j+1)}(t).
\]
The order-\(2m\) and order-\((2m+1)\) roots strictly interlace
when \(\lambda_m<\nu<\lambda_m+1/2\). This follows by comparing
the quadratic parameters of the two polynomials.

Both polynomial families have strictly log-convex binomial-normalized
coefficients for small positive \(\tau\), uniformly in
\(\nu\in[\lambda_m,\lambda_m+1/2]\). Indeed, in either degree \(d\)
the polynomial is
\((1+t)^d-\tau L t(1+t)^{d-2}+O(\tau^2)\) coefficientwise, where
\(L>0\) is the sum of its quadratic parameters \(\lambda_j\)
or their counterparts in \(R_\nu\). Hence the log of the
coefficient at index \(k\), divided by \(\binom dk\) before
taking the logarithm, is
\[
 -\tau L\frac{k(d-k)}{d(d-1)}+O(\tau^2).
\]
Its second difference is \(2\tau L/(d(d-1))+O(\tau^2)>0\).
At fixed degree \(L\) is bounded away from zero uniformly on
the parameter interval.

At the largest upper angle \(\theta\) of \(P\), put
\(a=2-\tau\lambda_m\). The coefficient of the displacement
in (1), before the minus sign for the inverse multiplier, is
\[
 \frac{a}{4\sin\theta}
       -2\sin\theta\sum_{j<m}\frac{1}{a-(2-\tau\lambda_j)} .
\]
The contribution from the conjugate of the root is the first
term. This expression is positive for \(a>1\), since
\[
 8\sin^2\theta\sum_{j<m}\frac{1}{a-(2-\tau\lambda_j)}
 \le 8\lambda_m(m-1)<\tfrac12 .
\]
Here \(\sin^2\theta\le\tau\lambda_m\) and each denominator
is at least \(\tau\). Thus \(T_\eta^{-1}P\) has a largest
upper angle \(\theta_\eta<\theta\) for all sufficiently small
positive \(\eta\) at fixed \(\tau\). Set
\(\nu_\eta=(2+2\cos\theta_\eta)/\tau\). This lies in
\((\lambda_m,\lambda_m+1/2)\) for small \(\eta>0\).
The input \(R_{\nu_\eta}\) strictly interlaces \(P\), and its
largest upper angle coincides with that of \(T_\eta^{-1}P\).
The transformed coefficients retain positivity and strict log
convexity for small enough \(\eta\).

The inverse used here is
\[
 T_\eta^{-1}f(t)=\eta^{-1}\int_0^1 x^{1/\eta-1}f(tx)\,dx ,
\]
of the same form as an individual Euler resolvent in contiguity.
Consequently, even in arbitrarily large degrees, positivity,
strict log convexity, and strict input interlacing do not suffice
for a transfer that applies an inverse Euler factor only at the
lower order.
