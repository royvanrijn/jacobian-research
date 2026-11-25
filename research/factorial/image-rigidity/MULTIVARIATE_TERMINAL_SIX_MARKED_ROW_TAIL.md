# Shared terminal six with a sole large-reciprocal (1,2) row

Let the unequal rows of distinct ordinary exponent columns \(a,b\)
all be shared. Suppose one row is \((1,2)\) or \((2,1)\) and every
other row satisfies \(d/(5\ell+1)\le1/8\), where
\(\ell=\min(a_i,b_i)\) and \(d=|b_i-a_i|\). If
\[
 S=\sum_i\frac{(b_i-a_i)^2}{\max(a_i,b_i)}<4 ,
\]
their order-five and order-six factorial polynomials have
disjoint root arguments up to independent positive rescalings.
Their rational Bézout identity gives coprimality over every
characteristic-zero field.

In fact, the angular proof covers every profile of the same
form with residual order-five deficit \(p_R\le.86\), where
\(p_R=1-u/\sqrt{3/5}\) and \(u\) is the full order-five
outer geometric mean. The
[preceding marked-row proof](MULTIVARIATE_TERMINAL_SIX_MARKED_ROW_EXTENSION.md)
covers \(p_R\le.161\). The
[continuation checker](verify_multivariate_terminal_six_marked_row_tail.py)
covers \([.161,.86]\) with outward dyadic intervals.

## Curvature and reciprocal constraints

The marked row contributes \(1/2\) to \(S\).
Write \(S_R=S-1/2\) and \(U_R=-\log(1-p_R)\).
For a residual row put \(r=d/\ell\). Its reciprocal cap
implies \(r\le5/7\), since \(d/(5\ell+d)\le1/8\). Its
contribution to \(U_R\) is at most
\((d/2)\log(1+4d/(5\ell))\), by bounding all Gauss reciprocals
by \(d/(5\ell)\) in the outer-mean product. Dividing by
its curvature contribution gives
\[
 F(r)=\frac{1+r}{2r}\log(1+4r/5)\le\frac65\log(11/7).
\]
Indeed,
\[
 2r^2 F'(r)=\frac{(4/5)r(1+r)}{1+(4/5)r}-\log(1+(4/5)r)\ge0
\]
by \((1+r)/(1+4r/5)\ge1\) and \(\log(1+s)\le s\).
Summing gives
\(p_R\le1-(7/11)^{21/5}<.851\) for \(S<4\). The last
comparison is the integer-power inequality
\((149/1000)^5<(7/11)^{21}\).

For the first three continuation intervals the checker uses
the smaller residual reciprocal caps
\[
 \begin{array}{c|c}
 [.161,.250]&1/11\\
 [.250,.275]&11/100\\
 [.275,.295]&117/1000 .
 \end{array}
\]
Every residual unit-gap row has \(\ell\ge2\), because
\((\ell,d)=(1,1)\) would violate the tail cap. The reciprocal
is at most \(1/11\) in that case. For \(d\ge2\) use
\(x_5\le(1-u_R)/(5u_R-1)\), with \(u_R=1-p_R\), as in the
preceding marked-row proof. All remaining intervals use \(1/8\).

## Exact ray heights

The checker reuses the joint residual coefficient expansion and
outer-log ratio interval of the
[intermediate proof](MULTIVARIATE_TERMINAL_SIX_INTERMEDIATE_CAYLEY.md).
It fixes four moving rays by rational affine interpolation of their
cosines. The first knot at \(p_R=.161\) matches the preceding
marked-row proof. The normalized interior multipliers are at most
\(4/5\) on the actual and nominal arrays and on their straight
segments. The earlier dominance estimate therefore places all
roots in \(1/8<|t|<8\) throughout the homotopy.

With the unscaled Cayley transform \(X=(1+t)/(1-t)=x+iY\)
and \(t=R e^{i\theta}\), write \(c=\cos\theta\) and
\(s=\sin\theta>0\). Then
\[
 x=\frac{1-R^2}{1-2Rc+R^2},\qquad
 Y=\frac{2Rs}{1-2Rc+R^2}.
\]
The value of \(Y\) is unchanged by reciprocal radii and increases
with \(R\) through \(R=1\). Since \(-.9\le c\le .8\),
\[
 Y\ge\frac{16s}{65-16c}>\frac1{15}
\]
on the verified annulus. If \(c\le C\le1/5\), the bound
\(|x|\le63/(65-16C)\) follows from the decreasing function
\((1-R^2)/(1-2RC+R^2)\) for \(1/8\le R\le1\).
Its derivative has the sign of \(C(1+R^2)-2R\),
which is negative on that interval. For larger \(C\), use
\(|x|\le1/\sqrt{1-C^2}\), from the identity
\[
 (1-2Rc+R^2)^2-(1-c^2)(1-R^2)^2
     =(c(1+R^2)-2R)^2 .
\]
Here \(C\) is the upper endpoint of the cosine on a ray knot
interval.

On the earlier continuation intervals the checker encloses \(Y\)
using the fixed-point equation of the preceding marked-row proof
and the positive lower bound just proved. On the later intervals
it encloses both exact candidate heights
\[
       Y=\frac{c}{s}\ \pm\sqrt{\frac1{s^2}-x^2}.
\]
These follow from the ray equation
\(x^2+(Y-c/s)^2=1/s^2\). A candidate is discarded only when
its discriminant is strictly negative or it lies below the
positive lower bound. For every remaining candidate height,
one rotated imaginary projection of the Cayley polynomial
is strictly signed throughout the accepted box. In particular,
the checks include both heights when the ray crosses the right
half-plane and a fixed value of \(x\) permits two heights.

At a given \(p_R,\rho\), use the same affine homotopy to the
nominal marked array times the symmetric residual as in the
preceding proof. Then vary \(\rho\) to \(25/24\) within its
enclosure, and vary \(p_R\) to .161 along that nominal family.
The rays coincide with the preceding certificate there, so its
sector counts apply. The joined last sector's quintic root is
real and its two sextic roots are nonreal by the strict positivity
of the actual sextic factorial moment on the real line. This
proves disjoint arguments.

Run the exact checks from the repository root with

~~~sh
python3 -O research/factorial/image-rigidity/verify_multivariate_terminal_six_marked_row_tail.py
~~~

The interval endpoints for the replay are the dyadic grid points
used in the earlier checker. All 38 knot intervals use strict
outward enclosures. Optional integer arguments select the lower
deficit knots of the eleven cap ranges for parallel replays.
