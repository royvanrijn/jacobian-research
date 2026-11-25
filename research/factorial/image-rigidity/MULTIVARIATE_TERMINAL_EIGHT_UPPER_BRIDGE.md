# Shared terminal-eight separation from outer mean .88 to .95

Let \(a,b\) be distinct ordinary exponent columns whose unequal
rows are shared. Let \(u\) be the order-seven outer geometric mean
after normalizing the endpoints to one and removing the binomial
factors. If \(22/25\le u\le19/20\), their order-seven and order-eight
polynomial root arguments are disjoint under independent positive
rescalings. The dimension and positive exponents are unrestricted.

The [exact Cayley checker](verify_multivariate_terminal_eight_upper_bridge.py)
uses the quadratic coefficient model of the
[terminal-seven upper bridge](MULTIVARIATE_TERMINAL_SEVEN_UPPER_BRIDGE.md)
with the order-seven and order-eight coefficient links of the
[terminal-eight private-linear proof](MULTIVARIATE_TERMINAL_EIGHT_PRIVATE_LINEAR.md#reflected-means-and-skew-links).

## Ordinary row bounds

On a shared row of lower exponent \(\ell\) and gap \(d\), the
order-seven reciprocals are \(x_h=d/(7\ell+h)\), \(1\le h\le d\).
The row's squared outer mean bounds \(u^2\) from above:
\[
 u^2\le \prod_{h=1}^d\frac{7\ell+h}{7\ell+6d+h}.
\]
For \(d=1\) it follows that \((7\ell+1)/(7\ell+7)\ge u^2\).
For \(d\ge2\), put \(y=d/(7\ell)\). Then
\[
 u\le\left(\frac{1+y}{1+7y}\right)^{d/2}
 \le\frac{1+y}{1+7y},\qquad
 x_h\le y\le\frac{1-u}{7u-1}.
\]
Also \(x_h\le(u^{-2}-1)/6\) because \(u^{-2}\) is the
product of \(1+6x_h\) over the shared factors.

On an interval of lower endpoint \(u_-\), let \(\ell_-\) be the
least positive integer satisfying
\((7\ell_-+1)/(7\ell_-+7)\ge u_-^2\). The bound used for
the order-seven reciprocals is
\[
 t=\min\left(\frac{u_-^{-2}-1}{6},\
     \max\left(\frac1{7\ell_-+1},\frac{1-u_-}{7u_--1}\right)\right).
\]
Since \(u_-\ge22/25\), we have \(\ell_-\ge3\) and \(t\le1/22\).
The matching order-eight cap is \(7t/(8-t)\).
These caps also give the
[linked outer ratio bound](MULTIVARIATE_TERMINAL_EIGHT_PRIVATE_LINEAR.md#ordinary-reciprocal-caps)
on \(\rho=\log w/\log u\), with \(w\) the order-eight outer mean.

## Quadratic Cayley families

At a subinterval midpoint \(u_0\), let the symmetric reference
multipliers be \(u_0^{e_{7,k}}\) and \(u_0^{(49/48)e_{8,k}}\),
where \(e_{j,k}=k(j-k)/(j-1)\). The log correction \(\sigma_k\)
to the reference is a sum of the linked skew and mean defect
contributions and one outer-mean contribution. The latter uses
the vector \((e_{j,k})_k\) and scalar \(\log(u/u_0)\) in order seven,
and scalar \((49/48)\log(u/u_0)+(49/48-\rho)(-\log u)\) in order
eight. The cap-based bounds for all other scalars are those of
the [paired coefficient envelopes](MULTIVARIATE_TERMINAL_EIGHT_PRIVATE_LINEAR.md#reflected-means-and-skew-links).

If \(|\sigma_k|\le M_k<1\), expand its exponential to second order:
\[
 |\exp(\sigma_k)-1-\sigma_k-\sigma_k^2/2|
       \le M_k^3/[6(1-M_k)].
\]
The checker groups the linear scalars and each quadratic product,
and encloses the third-order remainder coefficientwise.

Use \(p=1-u=\lambda^2\) and \(t_{\mathrm{poly}}=(\lambda X-1)/
(\lambda X+1)\). If \(c_k\) are the normalized physical coefficients
then the transformed polynomial is
\(\sum_k c_k(\lambda X-1)^k(\lambda X+1)^{j-k}\).
The coefficient of \(X^h\) is \(\lambda^h\) times an integer
linear combination of the physical coefficients. The program
applies this map to each perturbation vector before multiplying
by its scalar interval. Since the reference coefficients at
\(k,j-k\) coincide, it first groups the reflected terms of the
integer map. This preserves their cancellations.

## Ray exclusions and root counts

At each accepted subinterval the transformed leading term
dominates on \(|X|=2\). Since \(4p<1\), the Cayley map has
nonzero denominator on this disk and its image lies in the
open left half-plane. The six ray parameters \(b\) are
strictly decreasing in their index and piecewise affine in \(u\).
Their physical arguments are \(\pi-2\arctan(\lambda b)\).
Writing \(X=x+iY\), use the ray equation
\[
 (Y-b)(1+pbY)=-pbx^2,\qquad 0<Y\le b .
\]
The checker encloses \(Y\) by three interval updates starting at
\([0,b]\). It checks \(0\le x\le2\), since \(X\mapsto-\overline X\)
corresponds to reversing coefficients on the same upper ray.
Reversal changes the signs of the skew scalars and preserves the
coefficient model. The projection signs are tested on
\(\Im(\widetilde P_j(X)\overline X^{\,h})\), with each scalar
introduced after projecting its complete coefficient contribution.
A rectangle is accepted only when some projection has a strict
sign, and an unresolved rectangle is bisected in \(u\) or \(x\).

At fixed \(u\), take the straight coefficient segment from the
actual arrays to the symmetric reference of ratio \(49/48\).
The coefficient model is convex as an enclosure in the
independent scalar, product and remainder intervals, and contains
both arrays. Then move \(u\) to \(22/25\) along the reference
arrays. The constant ratio \(49/48\) belongs to all the ratio
bounds. The root disk and ray exclusions hold along both paths.

At \(u=22/25\) put \(Z=t_{\mathrm{poly}}+t_{\mathrm{poly}}^{-1}\).
If \(a,b,c\) are the first three normalized interior coefficients
of the symmetric septic, its quotient after removing \(t_{\mathrm{poly}}+1\)
and dividing by the reciprocal power of \(t_{\mathrm{poly}}\) is
\[
 Z^3-Z^2-2Z+1+a(Z^2-Z-1)+b(Z-1)+c .
\]
If \(d,e,f,g\) are the first four interior coefficients of the
symmetric octic, its reciprocal quotient is
\[
 Z^4-4Z^2+2+d(Z^3-3Z)+e(Z^2-2)+fZ+g .
\]
The checker verifies their signs at \(\pm2\) and at the six
rational anchor abscissae \(Z=2(pb^2-1)/(1+pb^2)\).
The signs give one upper anchor root in each successive sector,
of orders \(8,7,8,7,8,7,8\), and the septic has the additional root
\(-1\). These counts continue along the deformations. In the
last sector joined to its lower reflection, the single septic
root stays real. The actual octic roots there are nonreal
because the eighth factorial moment is strictly positive on the
real line by its Gamma integral. Hence all arguments are disjoint
between the two orders, and rational Bézout proves coprimality.

The exact replay certifies 124,353 Cayley projection cells.

Replay from the repository root with

~~~sh
python3 research/factorial/image-rigidity/verify_multivariate_terminal_eight_upper_bridge.py
~~~
