# Common roots in the first sector for two comparable exponent rows

For an ordinary factorial binomial, delete every row on which its two
exponents agree. Suppose the remaining rows have gaps of one sign,
and at most two rows remain. Then at every consecutive pair of orders
\(n,N=n+1\), \(n\ge2\), the two moment polynomials have no common
root \(z\) with \(|\arg z|\le 2\pi/N\). The gap sizes may differ.
The same sector exclusion holds with any number of positive-gap rows
whose lower exponents are proportional to their gaps. The result
concerns the same complex point in both moments.

## A conditional moment comparison

Let \(U,V>0\) be random variables with all moments needed below, put
\(T=V/U\), and write \(P_m(z)=\mathbb E[(U+zV)^m]\). Use the tilt
\(d\mathbb P_n=U^n\,d\mathbb P/\mathbb E[U^n]\). Suppose a version of
\(g(T)=\mathbb E_n[U\mid T]\) extends to a positive \(C^2\) function
on \((0,\infty)\). Assume either \(g\) is constant or \(g'>0\) and
\(1/g\) is convex. Under these assumptions the sector conclusion
holds whenever \(T\) is nonconstant.

By conjugation it suffices to take \(z\) with \(0<\theta=\arg z
\le 2\pi/N<\pi\). Positive real \(z\) is excluded by positivity.
Let \(\phi(t)=\arg(1+zt)\), which increases from zero to \(\theta\).
If \(\theta\le\pi/n\), then
\(\operatorname{Im}(1+zt)^n>0\) for every finite \(t>0\).
Consequently a zero of \(P_n\) requires \(\theta>\pi/n\).

Define \(f(t)=(1+zt)^n\) and \(M(t)=g(t)(1+zt)\). The exact identities
are
\[
 \frac{P_n(z)}{\mathbb E[U^n]}=\mathbb E_n[f(T)],\qquad
 \frac{P_N(z)}{\mathbb E[U^n]}=\mathbb E_n[f(T)M(T)].
\]
In the nonconstant case,
\(M'(t)=g'(t)(1+z\lambda(t))\) where
\[
 \lambda(t)=t+\frac{g(t)}{g'(t)},\qquad
 \lambda'(t)=\frac{g(t)^3}{g'(t)^2}(1/g)''(t)\ge0.
\]
Thus the tangent argument \(\chi(t)=\arg M'(t)\) is nondecreasing
with \(0\le\chi(t)\le\theta\). A constant \(g\) gives
\(\chi(t)=\theta\).

The function \(n\phi(t)+\chi(t)\) increases strictly, starts below
\(\pi\) and ends above \(\pi\). Let \(t_*\) be its unique crossing.
For \(t<t_*\), the chord \(M(t_*)-M(t)\) is an integral of tangents
whose arguments lie in \([0,\chi(t_*)]\). The argument of its
product with \(f(t)\) therefore lies in \((0,\pi)\).
For \(t>t_*\), the chord \(M(t)-M(t_*)\) has argument in
\([\chi(t_*),\theta]\). The unwrapped argument of its product with
\(f(t)\) lies in \((\pi,2\pi)\), since \(n\phi(t)+\theta<N\theta\).
Both cases imply
\[
             \operatorname{Im}\bigl[f(t)(M(t)-M(t_*))\bigr]<0
                         \qquad(t\ne t_*).
\]
These expressions are integrable by the assumed moments. If
\(P_n(z)=0\), taking the expectation eliminates the subtracted
constant term. The strict sign on a set of positive probability
forces \(\operatorname{Im}P_N(z)<0\).

## Comparable factorial rows

Realize the factorial moments using independent standard
exponentials as the row variables. A row on which the exponents
agree contributes a positive constant factor to each moment
polynomial. If every row is deleted, the only root is negative
real. Otherwise reverse the two monomials if needed to make all
remaining gaps positive. This sends nonzero roots to their
reciprocals and preserves the sector.

For proportional lower exponents and gaps, \(U=T^\alpha\)
with \(\alpha\ge0\), so the condition on \(g\) holds. This covers
a single remaining row as well.

For two remaining rows, write
\[
                U=X^aY^b,\qquad V=X^{a+d}Y^{b+e},
                   \qquad d,e\ge1,\quad a,b\ge0.
\]
Exchange the row variables if needed so
\(r=e/d\) and \(\delta=b-ra\ge0\). Set \(t=T^{1/d}\) and
\(\nu=nb+1-r(na+1)\). Under the tilt, \(X,Y\) are independent
Gamma variables of shapes \(na+1,nb+1\) and rate one.
Changing from \((X,Y)\) to \((t,Y)\) gives \(X=tY^{-r}\).
Conditionally on \(t\), \(Y\) has density proportional to
\(y^{\nu-1}e^{-y-ty^{-r}}\). With
\[
      I_\mu(t)=\int_0^\infty y^{\mu-1}e^{-y-ty^{-r}}\,dy,\qquad
      g(T)=t^a\frac{I_{\nu+\delta}(t)}{I_\nu(t)},
\]
these integrals and their derivatives in \(t>0\) converge for every
real \(\mu\). On compact positive \(t\)-intervals the exponential
terms dominate the derivative integrands at both endpoints.

Write \(\mathbb E_{\mu,t}\) and \(\operatorname{Var}_{\mu,t}\) for
expectation and variance under the density defining \(I_\mu(t)\).
For \(W=Y^{-r}\),
\[
 \left(\log\frac{I_{\nu+\delta}}{I_\nu}\right)'
       = \mathbb E_{\nu,t}W-\mathbb E_{\nu+\delta,t}W,\qquad
 \left(\log\frac{I_{\nu+\delta}}{I_\nu}\right)''
       = \operatorname{Var}_{\nu+\delta,t}W
         -\operatorname{Var}_{\nu,t}W .
\]
The mean of \(W\) decreases strictly as \(\mu\) increases, since
raising \(\mu\) multiplies its density by a strictly decreasing
power of \(W\). The variance decreases as well. This implies
that \(\log g\) is increasing and concave as a function of \(t\).
The change \(t=T^{1/d}\) is increasing and concave, so \(\log g\)
is increasing and concave in \(T\). It follows that \(1/g\) is convex.
The derivative \(g'>0\) if \(a>0\) or \(\delta>0\), and otherwise
\(g\) is constant.

## The variance comparison

At fixed \(r,t>0\), \(W\) has density proportional to
\(w^{-\mu/r-1}\exp(-w^{-1/r}-tw)\). For two independent copies
\(W_1,W_2\), set \(u=\sqrt{W_1W_2}\) and
\(h=|W_1-W_2|/2\). The change
\(W_1=u e^x,\ W_2=u e^{-x},\ h=u\sinh|x|\) gives the conditional
density in \(h>0\), up to a constant depending on \(u\),
\[
  \frac{\exp\bigl[-2t\sqrt{u^2+h^2}
       -2u^{-1/r}\cosh(\operatorname{arsinh}(h/u)/r)\bigr]}
       {\sqrt{u^2+h^2}}.
\]
Its logarithmic derivative in \(h\) is
\[
 -\frac{h}{u^2+h^2}
 -\frac{2th}{\sqrt{u^2+h^2}}
 -\frac{2u^{-1/r}\sinh(\operatorname{arsinh}(h/u)/r)}
        {r\sqrt{u^2+h^2}}.
\]
All three terms increase strictly with \(u>0\) at fixed \(h>0\).
The conditional densities at larger \(u\) therefore have increasing
likelihood ratios as functions of \(h\). In particular,
\(\mathbb E[h^2\mid u]\) increases with \(u\).

The conditional density is independent of \(\mu\), and increasing
\(\mu\) by \(\epsilon>0\) multiplies the marginal density of \(u\)
by \(u^{-2\epsilon/r}\), up to normalization. This shifts its
expectations of increasing functions downward. Both ordering
statements follow by writing the covariance of a monotone function
and a likelihood ratio as half the expectation of their paired
differences. All moments used here are finite from the exponential
factors in the density of \(W\). Since
\(\operatorname{Var}(W)=\mathbb E[(W_1-W_2)^2]/2=2\mathbb E[h^2]\),
this proves the variance comparison.

## Scope of the conditional argument

The conditional assumptions alone allow a common root at a larger
argument. Take \(U=1\) and \(V=T\) with values \(1/2,1,3/2\)
of probabilities \(1/18,8/9,1/18\), respectively. Direct expansion
gives \(P_4(z)=(1+2z+13z^2/12)^2\). Each of its two nonreal zeros
is also a zero of \(P_3(z)\), by \(zP_4'(z)=4(P_4(z)-P_3(z))\).
This is a positive moment law with the required constant \(g\).
