# Uniform closure of joint factorial curvature and center drift

Consider finite profiles of unequal ordinary exponent rows \((a_i,b_i)\).
Write \(A_i=\max(a_i,b_i)\), \(d_i=b_i-a_i\), and
\(s_i=d_i^2/A_i\). Assume \(\sum_i s_i\le S_0<\infty\).
The [joint curvature representation](MULTIVARIATE_JOINT_FACTORIAL_CURVATURE.md)
assigns to a row its actual factorial measure \(\mu_{a_i,b_i}\)
on the square, and to the profile their sum \(\mu\).

Under a curvature bound, the only diffuse measure that can be added
to the retained factorial row measures is a multiple of \(r\,dr\)
on the diagonal \((r,r)\), \(0<r<1\). This determines the joint
curvature arrays simultaneously at every order. A separate signed
diffuse scalar determines the possible renormalized center drift.
Below we give uniform truncation errors for both conclusions.

## All-order curvature error for an escaping row

Let \(u,v\ge0\) be real and set \(j=u+v+2\). For an increasing
row \(a\le b\), write \(d=b-a>0\). Its joint curvature moment
from the factorial measure is exactly
\[
 \begin{aligned}
 K_{a,b}(u,v)
 &=\int x^u y^v\,d\mu_{a,b}(x,y)\\
 &=d^2\int_0^\infty e^{-(aj+dv)t}
    \left(\frac{1-e^{-dt}}{dt}\right)^2
    \frac{t}{e^t-1}\,dt .
 \end{aligned}                                                       \tag{1}
\]
Reversal interchanges \(u,v\). This uses \(x=e^{-at}\),
\(y=e^{-bt}\) in the joint measure. The two bracketed scalar
factors before squaring are between \(e^{-dt}\) and \(1\),
and between \(e^{-t}\) and \(1\), respectively. Hence
\[
 \frac{d^2}{aj+d(v+2)+1}\le K_{a,b}(u,v)
       \le \frac{d^2}{aj+dv}
\]
when \(a>0\). Put \(s=d^2/b\) and \(\eta=d/b\). The weaker
orientation-independent bounds
\[
 \frac{s}{j+1/b}\le K_{a,b}(u,v)\le \frac{s}{j(1-\eta)}
 \tag{2}
\]
follow. The same reasoning also yields (2) for all real
moment indices \(u,v\) in its domain.

Choose an integer cutoff \(B>S_0\), and for rows with \(A_i>B\) put
\(\tau_B=\sum s_i\) and \(\eta_B=\sqrt{S_0/(B+1)}<1\).
Every such row is shared, with \(|d_i|/A_i\le\eta_B\).
Let \(K_{\rm tail}\) be their summed joint curvature.
Summing (2) gives the uniform bound
\[
 \sup_{u,v\ge0}
 \left|(u+v+2)K_{\rm tail}(u,v)-\tau_B\right|
 \le
 \tau_B\max\!\left(\frac1{2(B+1)},\frac{\eta_B}{1-\eta_B}\right).
 \tag{3}
\]
The measure \(r\,dr\) on \((r,r)\) has moment \(1/(u+v+2)\).
Thus (3) identifies and controls the diffuse diagonal measure
jointly across all orders.

For integer \(n\ge2\), curvature at the interior coefficient
indices is \(K(u,v)\) with \(u+v+2=n\).
If two curvature arrays differ by at most \(\epsilon/n\) at
those indices, their normalized coefficient logs differ by at
most \(\epsilon k(n-k)/(2n)\) at index \(k\).
Indeed both normalized arrays have zero endpoint logarithms;
comparison with the quadratic of constant second difference
\(\epsilon/n\), and its negative, gives the estimate.
Therefore (3) is also an all-order error bound for the normalized
log coefficient arrays.

## Exact closure of the measures

The sequence of joint moment functions is given the weighted
uniform metric
\[
 \|K_1-K_2\|_*=\sup_{u,v\ge0}(u+v+2)|K_1(u,v)-K_2(u,v)|.
\]
A shared row has finite norm by (2). A private row can have
infinite norm. For example, for a private row \((0,d)\), \(K(u,0)\) can be nonzero and
independent of \(u\). Consequently the metric is used on
differences of profiles with the same private rows, and private
multiplicities stabilize along the subsequences below.

For each fixed \(B\) only finitely many row types have \(A_i\le B\),
and their multiplicities are bounded because \(s_i\ge1/B\).
A diagonal subsequence retains fixed multiplicities of every
type and makes total curvature converge to some \(S_*\le S_0\).
Let \(S_{\rm ret}\) be the sum of the retained \(s_i\), and put
\(\tau=S_*-S_{\rm ret}\ge0\). Applying (3) at successively larger
cutoffs shows that the joint curvature converges in the weighted
metric, after stabilizing private rows, to the moment function of
\[
                 \mu_*=\sum_{\rm retained}\mu_{a_i,b_i}
                       +\tau\, (r,r)_*(r\,dr).
 \tag{4}
\]
In detail, at a fixed cutoff the retained contributions agree
eventually. The total tail curvature tends to the difference between
\(S_*\) and the cutoff-retained curvature. Apply (3) to both the
finite tails and the retained limit tails, then let \(B\) increase.
Boundedness of the total masses of all retained measures follows
by using (2) when \(|d_i|/A_i\le1/2\). For the remaining rows use
\(K(0,0)\le\log\binom{2|d_i|}{|d_i|}\le2|d_i|\log2
\le4s_i\log2\), directly from the factorial ratios.
Convergence of the nonnegative measures also holds weakly on
the square, since polynomial moments converge and the masses
are bounded.

Conversely, any countable sum of retained ordinary rows of finite
curvature with a nonnegative diffuse mass in (4) can be approximated
by finite profiles. Take longer finite truncations of the rows
and add many unit-gap rows of bases tending to infinity, with
total curvature approaching \(\tau\). Equation (3) proves
convergence throughout the same metric. For a prescribed curvature
band the resulting finite total curvatures converge to the desired
limiting value.

The retained rows in (4) and \(\tau\) are unique from all joint
moments. The diagonal measure is disjoint from the primitive
curve supports of the retained unequal rows outside their
measure-zero endpoints. On a fixed primitive curve there are only
finitely many scale factors \(g\) with retained positive
multiplicity, because the curvature of \(g(a,b)\) is linear in
\(g\). Hence the density uniqueness argument for the
[joint curvature representation](MULTIVARIATE_JOINT_FACTORIAL_CURVATURE.md#a-common-measure-for-all-orders)
also applies to this closure.

## Uniform signed center correction

For \(n\ge1\) real, let \(\mu_n\) denote the logarithm of the
positive endpoint center, and put \(\Delta=\sum_i d_i\).
For this paragraph use \(\mu_n\) for the logarithmic scalar
center, distinct from the curve measure \(\mu\).
Define the renormalized derivative
\[
 Z(n)=n^2\left(\partial_n\mu_n-\frac{\Delta}{n}\right)
     =n^2\sum_i\int_{a_i}^{b_i}
                    \left(x\psi'(1+nx)-\frac1n\right)\,dx .
 \tag{5}
\]
The integral follows by writing
\(\mu_n=\sum_i[\log\Gamma(1+n b_i)-\log\Gamma(1+n a_i)]/n\)
and differentiating its digamma integral.

The trigamma estimate from the
[fixed-terminal compactification](MULTIVARIATE_FIXED_TERMINAL_COMPACTIFICATION.md#uniform-estimates-for-rows-of-large-size)
is
\[
 \psi'(1+t)=\frac1t-\frac1{2t^2}+R(t),\qquad
 0\le R(t)\le \frac1{4t^3}\quad(t>0).
\]
For the tail at cutoff \(B\), write
\(L_B=\sum_{A_i>B}\log(b_i/a_i)\) in the actual orientation.
As all tail rows are shared, application of the estimate in (5)
gives, uniformly for \(n\ge1\),
\[
 \left|Z_{\rm tail}(n)+\frac{L_B}{2}\right|
 \le \frac{\tau_B}{4n(1-\eta_B)^2(B+1)},\qquad
 |L_B|\le \frac{\tau_B}{1-\eta_B}.
 \tag{6}
\]
The first bound uses
\(|d_i|/\min(a_i,b_i)^2\le s_i/[(1-\eta_B)^2(B+1)]\);
the second uses \(|\log(b_i/a_i)|\le |d_i|/\min(a_i,b_i)\).

As in the fixed-terminal extraction, pass to a subsequence on
which the total signed log ratio over shared rows converges and
subtract the retained contribution. By (6), the remaining signed
mass \(\lambda\) satisfies \(|\lambda|\le\tau\). Equations (5)
and (6) show that the whole renormalized center derivative
converges uniformly, as a difference from the limiting function,
to its retained-row function minus \(\lambda/2\). Any private
contributions are exactly retained after a finite cutoff, so they
cancel in the comparison. In particular, for every \(1\le m<n\),
the center difference after subtracting \(\Delta\log(n/m)\)
has diffuse contribution
\[
               -\frac{\lambda}{2}\left(\frac1m-\frac1n\right).
\]
The error divided by \(1/m-1/n\) is bounded by the uniform
error in (5), by integrating the derivative with weight \(1/t^2\).

Every \(\lambda\in[-\tau,\tau]\) is attainable independently
of the retained rows. At a large base \(A\), use unit-gap rows
with the two orientations in numbers asymptotic to
\(A(\tau+\lambda)/2\) and \(A(\tau-\lambda)/2\).
Their curvature and signed log ratio tend to \(\tau,\lambda\).
One can also prescribe the signed integer degree difference
\(\Delta\): balance any discrepancy by unit-gap rows of the
appropriate orientation with bases so large that their total
curvature and signed log ratio tend to zero. The bounds (3)
and (6) are uniform under this adjustment.

Thus the actual joint factorial curve densities impose a
specific diagonal diffuse law at every order, with one additional
signed parameter for the center drift. These restrictions
do not require a separate curvature approximation at each terminal.
