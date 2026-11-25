# Terminal six with large or two-sided private degree

Let \(a,b\) be distinct nonnegative integer exponent vectors. Write
\[
 S=\sum_{a_i\ne b_i}\frac{(a_i-b_i)^2}{\max(a_i,b_i)},\qquad
 P=\max\left\{\sum_{a_i=0}b_i,\ \sum_{b_i=0}a_i\right\}.
\]
The ordinary factorial moment polynomials of orders five and six are
coprime whenever \(P\ge3\) or \(S+P\ge34/5\). Both conditions are
uniform in the exponent sizes and the dimension. In fact these
polynomials have disjoint root arguments under independent positive
normalizations. The proof for \(P=3\) uses a comparison between
individual endpoint coefficients at consecutive orders.
When there are private variables on both sides, separation
also holds for private degrees \(2,2\) and \(2,1\).
For arbitrary profiles, a two-sided coefficient link
also gives separation at \(S\ge6\).

## Bounds at a fixed distance from an endpoint

For real \(n\ge k>0\) and positive Gamma shapes \(\gamma_i\), define
\[
 F_{n,k}
 =\prod_i
   \frac{\Gamma((n-k)a_i+kb_i+\gamma_i)}
        {\Gamma(na_i+\gamma_i)^{1-k/n}
         \Gamma(nb_i+\gamma_i)^{k/n}} .
\]
The ordinary endpoint-normalized coefficient is
\(c_{n,k}=\binom nk F_{n,k}\) at each integer order \(n\), with
\(\gamma_i=1\). We use the following comparison at each fixed \(k\):
\[
 F_{n,k}^{(n+1-k)/(n-k)}\le F_{n+1,k}\le F_{n,k}\qquad(n>k)
 \tag{1}
\]
when all shapes are at least \(1/2\). Reversing the rows gives
the same comparison at distance \(k\) from the other endpoint.

For one row put \(a=a_i,b=b_i,\gamma=\gamma_i\). Its contribution
to \(-\log F_{n,k}\) is
\[
 k\int_a^b\big[\psi(nx+\gamma)
            -\psi((n-k)a+kx+\gamma)\big]\,dx .
 \tag{2}
\]
The function \(t\mapsto t\psi'(t+\gamma)\) increases for \(t\ge0\).
Indeed, the trigamma integral, after scaling its integration
variable by \(t\), averages the kernel
\[
 K_\gamma(v)=\frac{v e^{-(\gamma-1/2)v}}{2\sinh(v/2)}
\]
against \(e^{-u}\,du\) at \(v=u/t\). The kernel decreases for
\(v>0\) when \(\gamma\ge1/2\), since
\(\sinh(v/2)/(v/2)\) increases. The value at \(t=0\) follows by
continuity. To differentiate (2), set \(y=(n-k)a+kx\).
When \(x\ge a\), we have \(nx\ge y\) and \(a/y\le1/n\) whenever
\(y>0\). Thus \(x\psi'(nx+\gamma)\ge a\psi'(y+\gamma)\).
When \(x\le a\) both comparisons reverse. The oriented integral
of the derivative is nonnegative in both cases. This gives the
upper bound in (1).

For the lower bound, rewrite (2) as
\[
 k(n-k)\int_a^b (x-a)\int_0^1
 \psi'\!\left(kx+(n-k)(a+u(x-a))+\gamma\right)\,du\,dx .
\]
The measure \((x-a)\,dx\) is positive in either orientation.
The argument of the trigamma increases with \(n\), and the trigamma
decreases. After division by \(n-k\) this expression therefore
decreases with \(n\). Summing the two comparisons over the rows
gives (1).

## Individual coefficient bounds

Reverse both columns if necessary so the larger private degree is
on the increasing side. Select its rows \((0,d)\), \(d>0\).
The endpoint-normalized contribution of one selected row at
order \(j\) and index \(k\) is
\[
 f_{j,k}(d)=\frac{(dk)!}{((dj)!)^{k/j}} .
\]
For \(d\ge1\),
\[
 f_{j,k}(d)\le f_{j,k}(1)^d .
 \tag{3}
\]
Indeed \((dk)!/(k!)^d\) is log-convex in integer \(k\), since its
ratio at adjacent indices is
\(\prod_{h=1}^d(dk+h)/(k+1)\), with each factor nondecreasing
in \(k\). The chord inequality gives (3).
Also \(f_{j,k}(d)\) decreases with real \(d>0\) after using the
Gamma interpolation, because its log derivative is
\(k[\psi(dk+1)-\psi(dj+1)]<0\) for \(0<k<j\).

The curvature of all nonselected rows is \(S-P\). Applying the
[discrete log-curvature bound](MULTIVARIATE_FIXED_TERMINAL_CURVATURE.md#large-curvature-and-the-endpoint-radii)
to these rows gives
\[
 c_{j,k}\le\binom jk f_{j,k}(1)^P
     \exp\!\left(-\frac{(S-P)k(j-k)}{2j}\right).
 \tag{4}
\]
In particular the region \(S+P\ge34/5\) is covered by the four
integer cases \(P=0,1,2,3\) with \(S-P\ge34/5-2P\),
and the case \(P\ge4\) with remaining curvature zero as a
lower bound. In the latter case cap the private power by four.
A single private gap at least three gives the additional cap
\(\binom jk f_{j,k}(3)\), even without remaining curvature.

## Exact ray inequalities

Let \(p_j(z)=1+\sum_{k=1}^{j-1} c_{j,k}z^k+z^j\).
For each cap \(B_{j,k}\) from (4), temporarily allow the coefficients
to vary independently in \([0,B_{j,k}]\).
Bounds on \(|z|=1/2\) and \(|z|=2\) give strict dominance
by the constant and leading terms, respectively, so every root
has \(1/2<|z|<2\).

For \(z=r e^{i\theta}\), \(x=\cos\theta\), and an integer rotation
\(h\in[0,j]\), consider
\[
 \frac{\Im(e^{-ih\theta}p_j(r e^{i\theta}))}{\sin\theta}
 =\sum_{k=0}^j c_{j,k}\,r^k\,\frac{\sin((k-h)\theta)}{\sin\theta}.
 \tag{5}
\]
The sine ratios obey the polynomial recurrence with first two
values \(0,1\) and next value \(2x\) times the previous value
minus the preceding one. The ratio for a negative index is the
negative of the positive-index ratio. Thus projection bounds
are exact rational calculations when \(x\) is rational.
To cap the exponential, use \(e^u\ge\sum_{\ell=0}^{64}u^\ell/\ell!\)
for \(u\ge0\) after raising (4) to the \(j\)-th power. Round the
algebraic upper coefficient bound upwards to a multiple of
\(1/10000\) by integer comparison of \(j\)-th powers.

The [exact checker](verify_multivariate_terminal_six_private_angles.py)
splits \(1/2\le r\le2\) into rational intervals, using
\(r^k\in[r_-^k,r_+^k]\) on each. It tests that some projection
(5) has a uniform strict sign on every interval. For all five
cases \(S+P\ge34/5\) the separating rays have cosines
\(81/100,2/5,-7/40,-17/25\). For a private gap at least three
use \(41/50\) as the first cosine and retain the other three.

Starting from \(1+z^j\), multiply all intermediate coefficients
by a parameter in \([0,1]\). Each selected ray remains zero-free.
The first ray lies between \(\pi/6\) and \(\pi/5\), the second
between \(\pi/5\) and \(\pi/2\), the third between \(\pi/2\) and
\(3\pi/5\), and the fourth between \(3\pi/5\) and \(5\pi/6\).
Thus the sector about the positive real axis contains one
sextic conjugate pair and no quintic root. The next three
upper sectors each contain one root, respectively quintic,
sextic and quintic. The final sector about the negative real
axis contains two sextic roots and one quintic root. The latter
is real by conjugation. The actual ordinary sextic is strictly
positive on the entire real axis by the factorial integral,
so no root argument agrees in that final sector either.

## Three unit degrees or a split of two and one

It remains to treat total private degree \(P=3\) without assuming
\(S+P\ge34/5\). Use (4) with three unit factors and no remaining
curvature. The checker validates the three separating rays with
cosines \(2/5,-1/5,-7/10\). By the same sector count as above,
only the first upper sector can contain an agreement of arguments.
It also validates a sextic ray with cosine \(39/50\), to the
left of the first such separator in angular order. The anchor
\(1+z^6\) shows that the first sextic upper root lies below this
angle. Every quintic root in the upper half-plane has angle
strictly greater than \(\pi/5\), since at smaller angles or at
\(\pi/5\) all imaginary terms are nonnegative and the linear
one is strictly positive. Hence a common argument would satisfy
\[
 \frac{\pi}{5}<\theta<\arccos(39/50).
 \tag{6}
\]

The checker excludes roots on this band with quintic radius
\(r\ge11/10\) or sextic radius \(s\le17/20\).
It does so by the same projection (5), now also enclosing
\(x\in[39/50,81/100]\) in the sine recurrence. This interval
contains all cosines in (6), since \(\cos(\pi/5)=(1+\sqrt5)/4\).

Write the quintic and sextic coefficients at index \(k\) as
\(c_k,d_k\). Remove the three private degrees with their actual
partition into rows, and apply (1) to all the remaining rows.
Using the comparison from the left for \(d_1,d_2,d_3\), and
from the right for \(d_4\), gives
\[
 d_1\le\tfrac45 c_1,\qquad
 d_2\le\tfrac23 c_2,\qquad
 d_3\le\tfrac35 c_3,\qquad
 d_4\le\tfrac{11}{10}c_3 .
 \tag{7}
\]
For each private partition \((1,1,1),(2,1),(3)\), the constants
are upper bounds for the binomial coefficient ratio times
the quotient of the corresponding exact private multipliers.
The checker verifies each by raising to the thirtieth power;
the normalized shared multipliers contribute a ratio at most
one by (1).

Put \(A_k=c_k/r^{5-k}\). For (6) all \(\sin(k\theta)\) for
\(1\le k\le4\) are positive, and \(\sin(4\theta)\le\sin(3\theta)\).
By \(r/s<22/17\), \(1/s<20/17\), and (7),
\[
 \sum_{k=1}^4\frac{d_k}{s^{6-k}}\sin(k\theta)
 \le \frac{31}{10}\sum_{k=1}^4 A_k\sin(k\theta).
 \tag{8}
\]
For the terms with \(A_1,A_2,A_3\), the separate factors before
rounding to \(31/10\) are respectively
\[
 \tfrac45(22/17)^4(20/17),\quad
 \tfrac23(22/17)^3(20/17),\quad
 (22/17)^2\big(\tfrac35(20/17)+\tfrac{11}{10}\big).
\]
The quintic imaginary equation gives the sum on the right of
(8), without \(31/10\), equal to \(-\sin(5\theta)\).
The sextic imaginary equation makes the left side at least
\(-\sin(6\theta)\), since \(d_5\ge0\) and
\(\sin(5\theta)<0\). But
\[
 -\sin(6\theta)>\frac{31}{10}\big(-\sin(5\theta)\big)
\]
on (6). The checker verifies this last inequality by dividing
by \(\sin\theta>0\) and enclosing the same sine recurrence over
one hundred rational subintervals of \([39/50,81/100]\).
This contradicts (8) and completes the separation for \(P=3\).

All normalizations of the variable used here are positive and
preserve arguments. At ordinary shapes the original polynomials
have rational coefficients, so rational Bézout identities
transfer coprimality to every characteristic-zero field.

## Private variables on both sides with one total degree two

Suppose the two private total degrees are \(2,2\) or
\(2,1\), up to reversing the columns. In the coefficient
box retain both sets of private rows. For private
degrees \(P_+,P_-\) the two contributions are capped by
\(f_{j,k}(1)^{P_+} f_{j,j-k}(1)^{P_-}\) using (3).
All other factors lie in \((0,1]\).

For \(2,2\), the same four rays with cosines
\(81/100,2/5,-7/40,-17/25\) give the original sector
count by homotopy to \(1+z^j\). The exact checker
verifies all four rays at both orders in 22 radial
intervals. This separates the two root arguments.

For \(2,1\), the last three rays give disjoint sectors
except possibly at the first upper roots, as in the
private degree-three comparison. The extra sextic ray
has cosine \(157/200\). The angle count from \(1+z^6\)
puts the first upper sextic root below that ray.
Thus a common argument must lie in
\(\pi/5<\theta<\arccos(157/200)\).
On this band the same rational projection tests give
\(r<109/100\) for a quintic root and \(s>9/10\) for
a sextic root. The two tests use coefficient caps for
both orientations, with all remaining rows unrestricted.

Keep the exact private coefficients in the adjacent link.
For either partition of the two increasing degrees, the
checker obtains
\[
 d_1\le\tfrac78 c_1,\qquad d_2\le\tfrac45 c_2,\qquad
 d_3\le\tfrac23 c_3,\qquad d_4\le\tfrac9{10}c_3 .
\]
These follow by the same thirtieth-power comparison as
(7), with the additional exact decreasing private
contribution at index \(j-k\).
Now \(r/s<109/90\) and \(1/s<10/9\). In the imaginary
equation (8), replacing the four link constants by
these values gives a multiplier at most \(5/2\).
The last sine comparison becomes
\(-\sin(6\theta)>(5/2)(-\sin(5\theta))\), verified
on the enclosing cosine band by the same interval
recurrence. The quintic and sextic equations
therefore contradict each other. The \(2,1\) replay
uses 69 radial or angle-radius intervals.

## Total curvature at least six

In this range use the coefficient cap
\(\binom jk \exp(-6k(j-k)/(2j))\) from total discrete
curvature. The checker verifies three separating rays
with cosines \(2/5,-11/50,-17/25\) at both orders,
and the ray with cosine \(3907/5000\) at order six.
The sector counts from \(1+z^j\) again leave only a
possible agreement in the first upper sector. Such
an agreement requires
\[
            \pi/5 <\theta<\arccos(3907/5000).
\]
The ray and radius projections on this band force
\(r<537/500\) for the quintic radius and \(s>459/500\)
for the sextic radius. All cosines in that band lie
between \(3907/5000\) and \(8091/10000\).

The two-sided link (1), with the binomial factors
restored, gives
\[
 d_1\le\tfrac65 c_1,\qquad d_2\le\tfrac32 c_2,\qquad
 d_3\le2\min(c_2,c_3),\qquad d_4\le\tfrac32c_3.
\]
The two bounds for \(d_3\) come from holding its distance
from the left or the right endpoint fixed. Set
\(A_k=c_k\sin(k\theta)/r^{5-k}\) for \(1\le k\le4\),
and let \(h=(537/500)/(459/500)\) and \(s_*=459/500\).
Since \(\sin(3\theta)\le\sin(2\theta)\) on the band and
\(\sin(4\theta)\le(621/1000)\sin(3\theta)\), the first
four terms of the sextic imaginary equation divided by
\(s^6\) are bounded by
\[
 L_1 A_1 + L_2 A_2 + L_4 A_3 + \min(M_2 A_2,M_3 A_3),
\]
where \(L_1=(6/5)h^4/s_*\), \(L_2=(3/2)h^3/s_*\),
\(L_4=(3/2)(621/1000)h^2\), \(M_2=2h^3\), and
\(M_3=2h^2/s_*\). The sine-ratio bound follows by checking
that \((621/1000)(4x^2-1)-(8x^3-4x)\) is positive on
the cosine band. The quotient
\[
 \frac{L_2 A_2+L_4 A_3+\min(M_2 A_2,M_3 A_3)}{A_2+A_3}
\]
is piecewise affine in \(A_2/(A_2+A_3)\). Its maximum is
at an endpoint or at \(M_2 A_2=M_3 A_3\).
Testing those three values and \(L_1\) gives a strict
upper bound of \(3469/1000\).

The quintic equation gives
\(\sum_{k=1}^4 A_k=-\sin(5\theta)\). In the sextic
equation the term \(d_5\sin(5\theta)/s\) is nonpositive,
so the sum of its first four divided terms is at least
\(-\sin(6\theta)\). But
\(-\sin(6\theta)>(3469/1000)(-\sin(5\theta))\) on the
whole band. The checker verifies positivity at the lower
cosine endpoint and a positive derivative bound throughout
for \(-U_5(x)+(3469/1000)U_4(x)\), where
\(U_j(\cos\theta)=\sin((j+1)\theta)/\sin\theta\).
This contradicts the two imaginary equations. The exact
replay uses 123 ray or angle-radius intervals in this
curvature case.
