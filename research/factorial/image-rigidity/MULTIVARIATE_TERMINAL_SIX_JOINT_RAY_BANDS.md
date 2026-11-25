# Coupled angular bands for fully shared terminal six

Let \(a,b\) be distinct ordinary exponent columns whose unequal rows
are fully shared. In each such row write \(\ell=\min(a_i,b_i)>0\) and
\(d=|b_i-a_i|\). Put \(x_*=\max d/(5\ell+1)\) and let \(u\) be the
order-five outer geometric mean after normalizing both endpoints
to one and removing the binomial coefficients. The following ranges
give disjoint root arguments for the order-five and order-six
polynomials under independent positive rescalings.

| Allowed \(u\) | Maximum \(x_*\) |
| --- | ---: |
| \((0,.28]\) | \(1/8\) |
| \([.28,.36]\) | \(2/25\) |
| \([.36,.44]\) | \(11/200\) |
| \([.44,.50]\) | \(7/200\) |
| \([.50,.55]\) | \(23/1000\) |
| \([.55,.60]\) | \(3/200\) |
| \([.60,.65]\) | \(9/1000\) |
| \([.65,.70]\) | \(1/200\) |
| \([.70,.72]\) | \(7/2000\) |

Every displayed endpoint is rational. The relative gap
\(\eta=\max d/(\ell+d)\) gives the sufficient condition
\(\eta\le5t/(1+5t)\) for a tabulated cap \(t\) on \(x_*\).
The direct \(x_*\) condition may be stronger. The test for \(u\)
also uses rational arithmetic after squaring the interval bounds,
because the product defining \(u^2\) is rational.

## Coupling the outer means

For one ordinary Gauss factor write
\(x=d/(5\ell+h)\), \(x'=d/(6\ell+h)\), \(1\le h\le d\).
If \(x_*\le t\), then
\[
 d/\ell\le5t/(1-t),\qquad
 x'\le5t/(6-t),\qquad
 \frac{25}{24}\le \frac{5x'}{4x}\le \frac{25}{4(6-t)} .
 \tag{1}
\]
Indeed \(x\ge d/(5\ell+d)\), and
\(x'/x=(5\ell+h)/(6\ell+h)\) lies between
\(5/6\) and \(5/(6-t)\). The first inequality in (1)
also follows directly from \(d/(5\ell+1)\le t\).

For the outer means \(u,w\) at orders five and six, set
\(U=-\log u,\ W=-\log w\). At one factor let
\(q=1+4x\), \(q'=1+5x'\).
The function \(s\mapsto\log(1+(25/24)(e^s-1))\)
is concave on \(s\ge0\) and vanishes at zero.
Therefore its quotient by \(s\) decreases, giving
\[
 L(t)\, U\le W\le \frac{25}{4(6-t)}\,U,\qquad
 L(t)=\frac{\log(1+(25/24)4t)}{\log(1+4t)} .
 \tag{2}
\]
For the upper bound use \(1+\lambda(q-1)\le q^\lambda\)
for \(\lambda\ge1\). Both inequalities sum exactly across
the factors. Thus the two coefficient boxes use the same \(u\)
and retain a narrower cross-order constraint than independent
bounds on the two curvatures.

## Skew and shape envelopes

Use the notation of the
[terminal-six coefficient links](MULTIVARIATE_TERMINAL_SIX_COEFFICIENT_LINKS.md)
and the formulas in the
[one-private moving-ray proof](MULTIVARIATE_TERMINAL_SIX_ONE_PRIVATE_SIDE.md#coefficient-envelopes).
In particular the shared arrays are
\[
 (1,u e^\alpha,v e^\chi,v e^{-\chi},u e^{-\alpha},1),\qquad
 (1,w e^\beta,y e^\eta,z,y e^{-\eta},w e^{-\beta},1),
\]
and their nonnegative defects are
\(V=\log(v/u^{3/2})\), \(Y=\log(y/w^{8/5})\),
\(Z=\log(z/(y w^{1/5}))\).
The local magnitudes of the skews and the defects are convex
in the local strength \(s=\frac12\log(1+(j-1)x_j)\)
and vanish at zero. Therefore the ratio of each magnitude to
\(s\) is nondecreasing. Evaluate at \(x_j=t\) for \(j=5\)
and at \(x_j=5t/(6-t)\) for \(j=6\), then sum across the rows.
This bounds \(|\alpha|,|\chi|,V\) by multiples \(A,C,D\)
of \(U\); it bounds \(|\beta|,|\eta|,Y,Z\) by multiples
\(B,E,H,J\) of \(W\). For \(\eta\) use
\(5|\eta|\le4|\beta|+|4\beta-5\eta|\).

With \(\ell_0\le L(t)\) the chosen rational lower bound in
(2), and \(\ell_1=25/[4(6-t)]\), the interval \([r,s]\)
in the following table means normalized multiplier
in \([u^s,u^r]\).
\[
\begin{array}{c|c|c}
j&\text{indices}&\text{exponent interval}\\ \hline
5&1,4&[1-A,1+A]\\
5&2,3&[3/2-D-C,\,3/2+C]\\
6&1,5&[\ell_0(1-B),\,\ell_1(1+B)]\\
6&2,4&[\ell_0(8/5-H-E),\,\ell_1(8/5+E)]\\
6&3&[\ell_0(9/5-H-J),\,9\ell_1/5].
\end{array}
\]
The checker proves each constant before using it. It encloses
\(\log q=2\operatorname{artanh}((q-1)/(q+1))\) by the
positive power sum through the term of degree \(31\) and
a geometric upper bound for the tail. It rounds each exponent
outwards to the nearest multiple of \(1/10000\).

## Moving rays and exact replay

The
[checker](verify_multivariate_terminal_six_joint_ray_bands.py)
gives the four continuous ray cosines by rational affine
interpolation on each subinterval of its displayed knot list.
The four rays stay ordered throughout. It bounds the powers
of \(u\) in the coefficient enclosures by rational intervals
of denominator \(10^8\), certifying each interval by integer
powers. A floating-point proposal can accelerate selection of
the endpoints; the same integer comparisons certify it, and
bisection provides an exact fallback.

The constant and leading coefficients dominate on radii
\(1/8\) and \(8\), respectively. On each two-dimensional
angle-parameter and radius box the checker uses the exact
projected sine ratios as in the
[fixed-ray proof](MULTIVARIATE_TERMINAL_SIX_PRIVATE_ANGLES.md#exact-ray-inequalities).
For a given parameter interval it encloses each sine ratio
by the polynomial recurrence, and bounds the projection from
below by a rational polynomial in the radius. Subdivision in
radius uses exact
[Bernstein coefficients](MULTIVARIATE_TERMINAL_SIX_DIFFUSE_BANDS.md#exact-ray-and-anchor-verification).
On each accepted radius interval, one of the projected
polynomials has all Bernstein coefficients strictly positive.
When the parameter interval is too wide to certify in this
way, it is bisected. The nine replays use \(1546,1398,2472,
2117,2248,3313,4495,6373,3541\) accepted radial interval
records, respectively.

At \(u=0\) the anchor polynomials are \(1+z^5,1+z^6\).
The anchor ray cosines lie in the four gaps between their
distinct nonreal upper arguments, and the checker proves
these inequalities exactly. To reach an actual profile,
first let \(u\) increase continuously from zero with the
canonical exponents \(1,3/2,3/2,1\) at order five and
\((25/24)(1,8/5,9/5,8/5,1)\) at order six.
These lie in the respective enclosures in every intervening
band, since \(L(t)\le25/24\le\ell_1\).
At the target \(u\), interpolate each coefficient to its
actual value inside its interval. The sector counts are
unchanged along both paths. The initial and final sectors
are counted with their conjugates. The final sector's one
quintic root remains real, and the actual even factorial
moment is strictly positive on the real axis. The argument
separation proves the characteristic-zero result by rational
Bézout.
