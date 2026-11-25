# Shared terminal-seven separation on the upper outer-mean interval

Let \(a,b\) be distinct ordinary exponent columns whose unequal rows
are shared. Let \(u\) be the order-six outer geometric mean after
normalizing the two endpoints to one and removing the binomial factors.
If \(22/25\le u\le19/20\), the order-six and order-seven polynomial
root arguments are disjoint under independent positive rescalings.
The dimension and positive exponents are unrestricted.

The [exact checker](verify_multivariate_terminal_seven_upper_bridge.py)
uses outward dyadic interval arithmetic in the Cayley variable. The
interval joins the upper end of the
[paired light bands](MULTIVARIATE_TERMINAL_SEVEN_PAIRED_LIGHT.md)
and the lower end of the
[near-flat interval](MULTIVARIATE_TERMINAL_SEVEN_NEAR_FLAT.md).

## Ordinary row bounds

Write \(\ell\) for the lower positive exponent and \(d\) for the gap
of one shared unequal row. The Gauss reciprocals at order six are
\(x_h=d/(6\ell+h)\), \(1\le h\le d\). As in the
[shared outer-mean comparison](MULTIVARIATE_TERMINAL_SEVEN_PRIVATE_LINEAR.md#shared-reciprocal-caps),
\[
 u^2\le\prod_{h=1}^d\frac{6\ell+h}{6\ell+h+5d}.
\]
For \(d=1\) this requires \((6\ell+1)/(6\ell+6)\ge u^2\).
For \(d\ge2\) put \(v=d/(6\ell)\). The row bound gives
\[
 u\le\left(\frac{1+v}{1+6v}\right)^{d/2}
    \le\frac{1+v}{1+6v},
 \qquad
 x_h\le v\le\frac{1-u}{6u-1}.
\]
Also \(x_h\le(u^{-2}-1)/5\), since \(u^{-2}\) is the product
of the factors \(1+5x_h\) across all shared rows.

Thus for an interval with lower endpoint \(u_-\), let \(\ell_-\) be
the least positive integer such that
\((6\ell_-+1)/(6\ell_-+6)\ge u_-^2\). The reciprocal cap is
\[
 t=\min\left\{\frac{u_-^{-2}-1}{5},\
       \max\left(\frac1{6\ell_-+1},\frac{1-u_-}{6u_--1}\right)\right\}.
\]
The lower endpoint \(u_-\ge22/25\) gives \(\ell_-\ge3\) and the second
entry of the maximum is less than \(1/19\). Hence the cap is always
at most \(1/19\). This uses ordinary row data and allows a new bound
on each subinterval without enumerating all profiles.

## Keep the quadratic terms through the Cayley expansion

Use the linked skew, defect, and outer-log ratio bounds of the
[paired light proof](MULTIVARIATE_TERMINAL_SEVEN_PAIRED_LIGHT.md#coefficient-enclosures)
at the cap \(t\) just obtained. Write \(U=-\log u\), \(W=-\log w\)
with \(w\) the septic outer mean, and \(\rho=W/U\). At order \(j\)
the symmetric reference multiplier at subinterval midpoint \(u_0\) is
\(u_0^{\rho_j e_{j,k}}\), where
\[
 e_{j,k}=\frac{k(j-k)}{j-1},\qquad \rho_6=1,\quad \rho_7=36/35 .
\]
The log correction \(\sigma_k\) to that reference is a sum
\(\sum_i A_i v_{i,k}\), where \(v_i\) are the
[linked coefficient vectors](MULTIVARIATE_TERMINAL_SEVEN_PAIRED_LIGHT.md#projection-replay).
In each order one variable accounts for the common outer log
and has vector \((e_{j,k})_k\). Its scalar at order six is
\(\log(u/u_0)\). At order seven it is
\((36/35)\log(u/u_0)+(36/35-\rho)U\). The other scalar bounds
come from the cap-based estimates for the skew and the
symmetric defects.

If \(|\sigma_k|\le M_k<1\), then
\[
 \left|\exp(\sigma_k)-1-\sigma_k-\frac{\sigma_k^2}{2}\right|
 \le \frac{M_k^3}{6(1-M_k)} .
\]
This follows from the exponential series and \(n!\ge6\) for \(n\ge3\).
Keep the shared variables \(A_i\) in the linear terms and the
shared products \(A_i A_l\) in the quadratic terms. For \(i=l\)
use \(A_i^2/2\) as the scalar. The vector of a quadratic term
is the entrywise product of its two linear term vectors.
Bound the third-order remainder coefficientwise.

Put \(p=1-u=\lambda^2\) and \(t_{\mathrm{poly}}=(\lambda X-1)/
(\lambda X+1)\). If \(c_k\) is the original normalized coefficient,
the transformed polynomial is
\[
 \widetilde P_j(X)=
 \sum_{k=0}^j c_k(\lambda X-1)^k(\lambda X+1)^{j-k}.
\]
Each coefficient of \(X^h\) is therefore \(\lambda^h\) times an
integer linear combination of the original coefficients. Apply
this integer map to each perturbation vector before multiplying
by its scalar interval. Pair the contributions at \(k\) and
\(j-k\) before multiplying by their common reference coefficient.
This preserves the cancellations from the symmetric reference
array. The checker uses outward rational power comparisons for
the reference coefficients and \(\lambda\), and exact integer
Cayley coefficients.

## Ray exclusions and sector counts

The five ray parameters \(b\) are strictly decreasing in their
index and piecewise affine in \(u\), with rational knots.
The physical ray has argument
\(\theta=\pi-2\arctan(\lambda b)\).
On each certified subinterval the transformed leading term
strictly dominates at \(|X|=2\) by its lower coefficient bound
and the absolute upper bounds for lower degrees. Hence all roots
satisfy \(|X|<2\). Since \(4p<1\), the denominator of the Cayley map
is nonzero on this disk and its image lies in the open left
half-plane.

Write \(X=x+iY\) on a physical upper ray. The ray equation is
\[
 (Y-b)(1+pbY)=-pbx^2,\qquad 0<Y\le b .
\]
The checker encloses \(Y\) by starting with \([0,b]\) and
intersecting with three interval updates of the equation.
It suffices to check \(0\le x\le2\). Reversal of the coefficients
corresponds to \(X\mapsto-\overline X\) on the upper ray.
It negates the skew variables and leaves the outer and mean
defect variables unchanged, preserving their intervals.
This reflection also preserves the enclosures for their
quadratic products and the independent remainders.
The checker tries the projections
\(\Im(\widetilde P_j(X)\overline X^{\,h})\), \(0\le h\le j\).
It projects each grouped contribution before multiplying by
the corresponding scalar interval. Every accepted rectangle
has a strictly signed projection. Otherwise it bisects \(u\)
or \(x\), with the subdivision limit treated as an error.

At a fixed \(u\) deform each actual array along a straight
coefficient segment to the symmetric reference with outer-log
ratio \(\rho_7=36/35\) at that \(u\). Both arrays lie in the
certified coefficient model. The model also contains their
segment because it permits independent intervals for the linear
scalars, product scalars and remainders. Then move \(u\) to
\(22/25\) along the symmetric reference arrays. The ratio
\(36/35\) is in every cap-based ratio interval.
The root disk and ray exclusions apply to both paths.

At \(u=22/25\) write \(Z=t_{\mathrm{poly}}+t_{\mathrm{poly}}^{-1}\).
The symmetric sextic quotient and the symmetric septic quotient
after removing \(t_{\mathrm{poly}}+1\) are the two cubics in the
[paired light anchor](MULTIVARIATE_TERMINAL_SEVEN_PAIRED_LIGHT.md#projection-replay).
The anchor ray abscissae are exact rationals
\(Z=2(pb^2-1)/(1+pb^2)\). The checker encloses the cubics'
values at the five abscissae and at \(\pm2\).
The signs give one upper root per sector in the order
\(7,6,7,6,7,6\). The septic has the additional root \(-1\).

The sector counts are preserved with the first and last sectors
joined to their respective lower reflections. The last joined
sector contains only one septic root, which stays real by
conjugation. The actual sextic is strictly positive for every
real coefficient ratio by the Gamma integral of the sixth power
of a nonzero binomial, so its two roots in the last joined sector
are nonreal. Thus the actual root arguments are disjoint.
The original factorial moment polynomials have rational coefficients
and Bézout gives coprimality over every characteristic-zero field.

The exact replay certifies 25,096 Cayley projection cells.
Replay from the repository root with

~~~sh
python3 research/factorial/image-rigidity/verify_multivariate_terminal_seven_upper_bridge.py
~~~
