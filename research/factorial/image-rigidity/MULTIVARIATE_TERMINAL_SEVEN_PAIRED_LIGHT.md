# Paired Bernstein separation for light shared rows at terminal seven

Let \(a,b\) be distinct ordinary exponent columns whose unequal
rows are shared. Write \(\ell_i=\min(a_i,b_i)\) and \(d_i=|b_i-a_i|\)
on these rows. Let \(u\) be the order-six outer geometric mean
after endpoint normalization and removal of binomial coefficients.
The order-six and order-seven factorial moments have disjoint root
arguments under independent positive rescalings on these bands.

| Range of \(u\) | Maximum \(d_i/(6\ell_i+1)\) |
| --- | ---: |
| \([.10,.42]\) | \(1/10\) |
| \([.42,.48]\) | \(3/40\) |
| \([.48,.62]\) | \(1/20\) |
| \([.62,.72]\) | \(1/16\) |
| \([.72,.75]\) | \(7/200\) |
| \([.75,.88]\) | \(1/40\) |

The [checker](verify_multivariate_terminal_seven_paired_light.py)
replays every bound by exact rational and outward dyadic arithmetic.
It is independent of the [near-flat proof](MULTIVARIATE_TERMINAL_SEVEN_NEAR_FLAT.md).

## Coefficient enclosures

For one Gauss factor let \(x=d/(6\ell+h)\), \(x'=d/(7\ell+h)\)
for \(1\le h\le d\). Then \(x'=6x/(7-xh/d)\), so when \(x\le t\),
\[
 x'\le\frac{6t}{7-t},\qquad
 1+\tfrac{36}{7}x\le1+6x'\le\frac{1+5x}{1-x/7}.
\]
For \(U=-\log u\) and \(W=-\log w\), where \(w\) is the
order-seven outer mean, this gives
\[
 \frac{\log(1+(36/7)t)}{\log(1+5t)}
 \le \rho=\frac WU
 \le1-\frac{\log(1-t/7)}{\log(1+5t)} .
\]
The two quotient monotonicities are those of the
[paired terminal-six proof](MULTIVARIATE_TERMINAL_SIX_PAIRED_LIGHT_ROWS.md).
The interval always contains \(36/35\).

For each order use a local cap \(t_j\) on its reciprocal,
and write \(s_t=\frac12\log(1+(j-1)t_j)\). The total of the
local strengths \(s\) equals \(U\) or \(W\). For \(k\ge1\),
\(\sum x_j^k\le(t_j^k/s_t)\sum s\), by concavity of the
logarithm and \(x_j\le t_j\). The order-six linked variables
are as in the [terminal-six coefficient note](MULTIVARIATE_TERMINAL_SIX_COEFFICIENT_LINKS.md).
The checker bounds the outer skew \(\beta\) by evaluating
its single-factor skew-to-strength ratio at the cap.
For each pair of indices \(i,j-1-i\),
\[
 \log\frac{(1+i x_j)(1+(j-1-i)x_j)}{1+(j-1)x_j}
 =\log\big(1+2a(1-a)(\cosh(2s)-1)\big),
 \qquad a=i/(j-1).
\]
This is the logarithm of a positive linear combination of
three real exponentials, hence convex in \(s\). The outer
skew magnitude is a positive sum of such terms and vanishes
at zero, so the ratio increases. Use
\(|4\beta-5\eta|\le6\sum x_j^4\) and
\(J_6\le(6/5)\sum x_j^5\) for the other links. The
mean defect \(Y\) is bounded by
\((G_{3/5}(s_t)/s_t)\sum s\), with
\(G_a(s)=as-\operatorname{artanh}(a\tanh s)\).
Its derivative \(a(1-a^2)\tanh^2(s)/(1-a^2\tanh^2(s))\)
is nondecreasing, giving convexity and the ratio bound.

For order seven the two signed skew links are
\[
 L=\alpha-\beta=(2\Delta^4g_0+3\Delta^4g_1+2\Delta^4g_2)/14,\qquad
 K=\alpha-3\beta+5\gamma=\Delta^6g_0/14 ,
\]
where \(g_h=\log(1+h x')\) for one factor. The orientation
reverses their signs. Centered second differences give the
single-factor bound \(|\gamma|\le(x')^2\); the fourth and sixth
derivatives give \(|L|\le3(x')^4\) and \(|K|\le(60/7)(x')^6\).
The order-seven means satisfy \(H=\log(y/w^{5/3})=\sum G_{2/3}(s)\)
and \(\log(z/w^2)=(9/5)H+E\) with
\(E=\sum(G_{1/3}(s)-(4/5)G_{2/3}(s))\).
The difference of the derivatives is
\[
 -\frac{8\tanh^4(s)}
        {81(1-\tanh^2(s)/9)(1-4\tanh^2(s)/9)}.
\]
Its magnitude is at most \(s^4/5\), giving
\(-s_t^4 W/25\le E\le0\). Bound \(H\) using convexity
and the ratio \(G_{2/3}(s_t)/s_t\).
Each logarithm bound in the checker uses a finite
atanh series and its explicit geometric remainder.

## Projection replay

Use \(e_{j,k}=k(j-k)/(j-1)\). At each interval midpoint
\(u_0\), the reference multipliers are
\(u_0^{e_{6,k}}\) and \(u_0^{(36/35)e_{7,k}}\). The log
correction to the reference in order six consists of
the common mean term \(e_{6,k}\log(u/u_0)\) and the
linked skew and defect arrays from the paired terminal-six
proof. In order seven the mean term is
\(e_{7,k}[(36/35)\log(u/u_0)+(36/35-\rho)U]\).
The other log contributions in its seven-degree array are
the following shared scalars times their respective vectors.
\[
\begin{array}{c|rrrrrrrr}
 \gamma&0&5/2&5/2&1&-1&-5/2&-5/2&0\\
 L&0&3/2&1/2&0&0&-1/2&-3/2&0\\
 K&0&-1/2&-1/2&0&0&1/2&1/2&0\\
 H&0&0&1&9/5&9/5&1&0&0\\
 -E&0&0&0&-1&-1&0&0&0 .
\end{array}
\]
For a log correction \(\sigma\) with \(|\sigma|\le A<1\),
use \(0\le e^\sigma-1-\sigma\le A^2/[2(1-A)]\).
Keep each shared variable together in the linear term.
The nonlinear remainder uses independent coefficient bounds.

On an upper ray \(t=r e^{i\theta}\), each projection uses
\(\sin((k-h)\theta)/\sin\theta\), bounded by the exact
recurrence across the ray cosine interval. For each linked
variable compute the Bernstein coefficients of its projected
polynomial in \(r\) on \([1/10,10]\) before multiplying by
the variable interval. The Bernstein map on this positive
radius interval has nonnegative entries. Subdivision uses
outward rounded de Casteljau averages, and only strictly
signed Bernstein coefficient arrays certify a box. The
checker bisects the \(u\) interval if a radial check reaches
its subdivision limit.

The actual arrays, the symmetric arrays, and the straight
coefficient segments joining them have interior multipliers
at most one. Since \((11/10)^7<2\), the constant and leading
coefficients dominate on the boundaries of the tested annulus.
The rational affine ray-cosine knots match across bands.
At fixed \(u,\rho\), first deform to the symmetric arrays
with the same means. The linked linear and independent remainder
enclosures contain this straight segment. Then vary \(\rho\)
to \(36/35\), and move \(u\) to \(1/10\) along the symmetric
arrays. These lie in the enclosures for every intervening band.

At that anchor, use \(X=t+t^{-1}\). The sextic quotient is
\(X^3-3X+a(X^2-2)+bX+c\), with \(a,b,c\) its first three
interior coefficients. After removing \(t+1\), the septic
quotient is \(X^3-X^2-2X+1+d(X^2-X-1)+e(X-1)+f\), with
\(d,e,f\) its first three coefficients. The checker proves
the signs at \(X=2,-2\) and at the five ray abscissae, giving
successive upper roots of orders \(7,6,7,6,7,6\).
The remaining septic root is real. The same counts hold at
actual profiles by continuity, joining the last sector with
its conjugate. There the actual sextic pair is nonreal
by strict positivity on the real axis from the sixth-power
factorial moment representation. Thus the actual arguments
are disjoint. Rational Bézout gives coprimality over every
characteristic-zero field.

Replay from the repository root with

~~~sh
python3 -O research/factorial/image-rigidity/verify_multivariate_terminal_seven_paired_light.py
~~~
