# Phase separation on the two upper saddle sheets

Let \(a_i,b_i\) be positive integers, discard equal rows, and put
\(d_i=b_i-a_i\), \(p_i=-a_i/d_i\), \(x_i(q)=a_i+d_iq\), and
\(\Delta=\sum_i d_i\). Use the principal logarithms on the upper
half-plane and the stationary functions
\[
 \begin{split}
 X(q)&=\sum_i d_i\log x_i(q)+\log(1-q)-\log q,\\
 F_a(q)&=\sum_i a_i\log x_i(q)-\log(1-q),\\
 F_b(q)&=\sum_i b_i\log x_i(q)-\log q=F_a(q)+X(q).
 \end{split}
\]
The [stationary map](MULTIVARIATE_SADDLE_MAP_TWO_SHEETS.md)
has degree two over \(\mathcal T=\{w:-\pi<\Im w<0\}\).

**Theorem.** Suppose \(\Delta\ne0\). For any distinct upper
half-plane points \(q_1,q_2\) with \(X(q_1)=X(q_2)\in\mathcal T\)
and equal real parts of \(F_a-\Delta q\), one has
\(e^{F_a(q_1)}\ne e^{F_a(q_2)}\).
Consequently, for a fixed such profile, consecutive common
roots cannot accumulate in a regular complex region where the
two upper saddles are exactly the contributors and have the
uniform nonvanishing-amplitude expansion required in the
[consecutive comparison](MULTIVARIATE_SHARED_SUPPORT_CONDITIONAL_SADDLE.md#comparing-the-two-terminal-orders).
The theorem concerns the comparison of the two saddles. It
does not require a contour expansion to prove their phase
separation.

## A bound on fibers in the exterior quadrants

There is at most one solution in \(\Re q<0,\ \Im q>0\) of
\(X(q)=w\) for any \(w\in\mathcal T\), counted with multiplicity.
By swapping \(a,b\) and taking \(q\mapsto1-\overline q\),
the same bound holds in \(\Re q>1,\ \Im q>0\). Indeed the new
\(X\) at \(1-\overline q\) is \(-\overline{X(q)}\), which
preserves the target strip.

For the left quadrant, use a large quarter-disk closed on the
negative real and positive imaginary axes. Indent around zero
and the real poles. As in the upper half-plane winding proof,
count oriented crossings of the ray at height \(\eta=\Im w\)
to the right of \(\Re w\). Initially take generic targets so
that the image of the boundary avoids \(w\) and the crossings
are transverse. The real-axis pieces have heights that are
integer multiples of \(\pi\). Each indentation at a negative
real pole has real part tending to \(-\infty\), so contributes
zero. The indentation at zero travels from argument \(\pi\)
to \(\pi/2\) and contributes one upward crossing exactly
when \(-\pi<\eta<-\pi/2\).

On the outer arc \(q=R e^{i\theta}\), \(\pi/2\le\theta\le\pi\),
the imaginary part tends uniformly to
\(\Delta\theta+\pi(D_--1)\), where \(D_-=\sum_{d_i<0}|d_i|\),
and its angular derivative tends uniformly to \(\Delta\).
The real part is \(\Delta\log R+O(1)\). Set
\(D_+=\sum_{d_i>0}d_i\), \(D=D_++D_-\). If \(D\ge2\)
the initial height is \(\pi(D-2)/2\ge0\).
If \(\Delta>0\) it then increases and cannot hit the target
height. If \(\Delta<0\) the crossing is too far left to count,
and if \(\Delta=0\) the constant limiting height is nonnegative.
For \(D=1\) the only possible contribution is one upward
crossing for \(-\pi/2<\eta<0\) in the case of an increasing
unit row. The arc goes from height \(-\pi/2\) to zero and
has real part tending to \(+\infty\). These are all outer
contributions.

It remains to count the image of the imaginary axis, traced
upward. For \(y>0\), put \(\ell_i=a_i/|d_i|\). Its height is
\[
 B(y)=\Im X(iy)=-\frac\pi2+\sum_i |d_i|\arctan(y/\ell_i)
                                     -\arctan y .
\]
The derivative is
\[
 B'(y)=\sum_i\frac{a_i}{\ell_i^2+y^2}-\frac1{1+y^2}.
\]
If some \(\ell_i\le1\), its term alone is at least
\(1/(1+y^2)\), by \(a_i\ge1\). Thus the height is increasing
unless it is constant in the single-row \(a_i=d_i=1\)
case. If all \(\ell_i>1\), the ratio of the positive sum
to \(1/(1+y^2)\) strictly increases in \(y\).
Hence there can be at most one decreasing interval followed
by an increasing interval. The height starts at \(-\pi/2\)
and tends to \(\pi(D-2)/2\ge-\pi/2\). For \(\eta>-\pi/2\)
there is at most one upward crossing. If the outer arc also
contributes for this target, then \(D=1\), and the limiting
height on the axis is \(-\pi/2\), so there is no axis crossing.

Consider \(\eta<-\pi/2\), where it is necessary to compare
the axis crossings with the already counted crossing at zero.
If there are two axis crossings let \(0<y_1<y_2\) be their
locations, with downward and upward directions, respectively.
This case requires all \(\ell_i>1\). For \(\ell>0\) put
\[
 L(\ell)=\arctan(y_2/\ell)-\arctan(y_1/\ell),\qquad
 M(\ell)=\frac12\log\frac{\ell^2+y_2^2}{\ell^2+y_1^2}.
\]
Equality of heights gives \(\sum_i |d_i|L(\ell_i)=L(1)\).
In particular \(L(\ell_i)\le L(1)\). The tangent of
\(L(\ell)\) is \((y_2-y_1)\ell/(\ell^2+y_1y_2)\), and
\(0<L(\ell)<\pi/2\). Since \(\ell_i>1\), the inequality
implies \(\ell_i\ge\ell_*=\max(1,y_1y_2)\).
Also \(L(\ell_*)=L(1)\).
The ratio \(M(\ell)/L(\ell)\) decreases with \(\ell\):
it is the weighted mean of \(t\) with density
\((1+t^2)^{-1}\) on \([y_1/\ell,y_2/\ell]\), and the
mean is increasing in each endpoint of the interval.
As a result,
\[
 \sum_i d_i M(\ell_i)
 \le \sum_i |d_i|L(\ell_i)\frac{M(\ell_*)}{L(\ell_*)}
 = M(\ell_*).
\]
For \(y_1y_2\ge1\) one has
\(M(\ell_*)+M(1)=\log(y_2/y_1)\).
For \(y_1y_2\le1\) one has
\(M(\ell_*)+M(1)=2M(1)\le\log(y_2/y_1)\), since
\((y_2-y_1)(1-y_1y_2)\ge0\).
Consequently
\[
 \Re X(iy_2)-\Re X(iy_1)
 =\sum_i d_i M(\ell_i)+M(1)-\log(y_2/y_1)\le0 .
\]
Whenever the upward crossing counts, the downward crossing does
too. Thus the total axis contribution for \(\eta<-\pi/2\) is
nonpositive. If the height only touches \(\eta\), generic
perturbation avoids that case.

The total winding is at most one in both ranges for \(\eta\).
The arc asymptotics and the endpoint limits justify taking the
large radius and indentations uniformly for each generic target.
The bound extends to any target in the strip, since zeros in
the open quadrant, counted with multiplicity, persist in small
disks under a sufficiently small perturbation of the target. The
non-generic boundary images can be avoided by such perturbations:
the imaginary axis has a real-analytic parametrization, and the
exceptional heights of its extrema and limiting endpoints form
a finite set.

## Exclusion on a vertical fiber

If the theorem's two exponentials were equal, their log moduli
would agree, and equal real parts of \(F_a-\Delta q\)
would force \(\Re q_1=\Re q_2=x\), because \(\Delta\ne0\).
The quadrant bounds imply \(0\le x\le1\).
Write \(q_j=x+iy_j\) with \(0<y_1<y_2\), and now put
\(\ell_i=|x-p_i|\). For the same function \(M\) used
above, with \(M(0)\) defined by continuity, equality of
the moduli of \(e^{F_a},e^{F_b}\) gives
\[
 \sum_i a_i M(\ell_i)=M(1-x),\qquad
 \sum_i b_i M(\ell_i)=M(x).
\]
Here \(M\) is strictly decreasing on \([0,\infty)\).
Since all the row weights are integers at least one,
each \(\ell_i\) is at least \(\max(x,1-x)\).

For \(y_1\le y\le y_2\), consider
\[
 f(y)=\Re X'(x+iy)=
 \sum_i\frac{|d_i|\ell_i}{\ell_i^2+y^2}
    -\frac{x}{x^2+y^2}-\frac{1-x}{(1-x)^2+y^2}.
\]
The ratio of the positive sum to the two subtracted terms
is nondecreasing in \(y^2\). In fact, on differentiating
the ratio and multiplying by the square of its positive
denominator, each pair of one numerator term and one
denominator term contributes a nonnegative multiple of
\(\ell_i^2-x^2\) or \(\ell_i^2-(1-x)^2\), respectively.
Hence \(f\) has at most one sign change unless it is
identically zero.

The identities \(X(q_1)=X(q_2)\) and
\((F_a-\Delta q)'=-qX'\) give
\[
 \int_{y_1}^{y_2} f(y)\,dy=0,\qquad
 \int_{y_1}^{y_2} y f(y)\,dy=0 .
\]
A nonzero continuous function with at most one sign change
cannot satisfy both equations: when the first integral is
zero, subtract its multiple by the sign-change location
from the second integral to obtain a strict sign.
Thus only the identically-zero case remains.

For \(0<x<1\) that case is impossible. If \(x\ne1/2\),
the numerator-denominator ratio is strictly increasing
because the two denominator distances differ and are
both positive. If \(x=1/2\), every row pole lies strictly
outside \([0,1]\), so every row distance is strictly
larger than \(1/2\), again giving strictness. For \(x=0\),
an identically zero ratio difference requires every
\(\ell_i=1\) and \(\sum_i |d_i|=1\). This is the single
row \((a,b)=(1,2)\). On \(q=iy\), the imaginary part of
\(F_a\) is \(2\arctan y\), strictly increasing in \((0,\pi)\).
Thus the two distinct heights still cannot give equal
exponentials. The case \(x=1\) follows by reversal of
the columns and conjugation, which swaps \(F_a,F_b\).
This proves the theorem.

## Equal total degrees at and above the real cusp

**Theorem.** Suppose \(\Delta=0\) and
\[
 \max_{0<x<1}\rho(x)\ge1,\qquad
 \rho(x)=x(1-x)\sum_i\frac{d_i^2}{x_i(x)} .
\]
Then the same two-saddle phase separation holds on the whole
strip \(\mathcal T\). In this case the equal-modulus locus is
a single graph over the imaginary part of \(w\). Along it,
the phase difference decreases strictly from \(i\pi\) at
\(\Im w=-\pi\) to zero at \(\Im w=0\).

First there is a fixed divider for the two branches. Choose an
\(x_0\in(0,1)\) with \(\rho(x_0)\ge1\). For \(y>0\) write
\[
 \Im X(x_0+iy)=\sum_i |d_i|\arctan(y/\ell_i)-
                    \arctan(y/x_0)-\arctan(y/(1-x_0)),
 \qquad \ell_i=|x_0-p_i|.
\]
Every \(\ell_i>\min(x_0,1-x_0)\) by the pole positions.
Consider the list of reciprocals \(1/\ell_i\), with
\(|d_i|\) copies for each row. Its sum is at least
\(1/x_0+1/(1-x_0)\) by \(\rho(x_0)\ge1\), while every
entry is strictly smaller than \(1/\min(x_0,1-x_0)\).
For completeness, one can reduce entries while preserving
nonnegativity until the sum is exactly \(1/x_0+1/(1-x_0)\).
At fixed sum and with each entry at most \(1/\min(x_0,1-x_0)\),
one minimizes the sum of the strictly concave function
\(h\mapsto\arctan(yh)\) by combining masses into the
largest possible entries. Indeed transferring mass from a
smaller positive entry to a larger entry below the cap
decreases the sum strictly. The resulting list consists,
up to zero entries, of \(1/x_0,1/(1-x_0)\).
Therefore \(\Im X(x_0+iy)>0\) for all \(y>0\).
This uses integrality of the gap multiplicities.

The two inverse branches given by the proper covering can
be labeled \(q_L,q_R\) with \(\Re q_L<x_0<\Re q_R\).
To see that both sides are occupied, take a target just
above height \(-\pi\). On the corresponding real boundary
intervals \((p_L,0)\) and \((1,p_R)\), the derivative
\(X'=\sum_i d_i^2/x_i-1/q+1/(q-1)\) is positive,
the height is \(-\pi\), and the real parts each cover
the real line. Here \(p_L=\max_{d_i>0}p_i\) and
\(p_R=\min_{d_i<0}p_i\) are finite because \(\Delta=0\)
and some gaps are nonzero. The local inverses therefore
give one branch on each side near the boundary. They
cannot cross the divider. In particular, the phase difference
\[
 D(w)=F_a(q_R(w))-F_a(q_L(w)),\qquad
 D'(w)=q_L(w)-q_R(w)
\]
satisfies \(\Re D'<0\) everywhere in the strip.

For a fixed \(-\pi<\eta<0\), as \(\xi=\Re w\to+\infty\)
the inverse branches converge to \(0,p_R\), respectively.
At zero use the local coordinate for the singular part
\(-\log q\), and at \(p_R\) use the singular part
\(-D_R\log(p_R-q)-i\pi\), with
\(D_R=\sum_{p_i=p_R}|d_i|\). The regular parts are real
on the interval of approach, so the inverse branches given
by these coordinates are upper for the chosen height.
They account for the degree two of the map. As
\(\xi\to-\infty\) the same construction gives \(p_L,1\)
as the limits. The latter also follows by reversal of the
columns, which replaces \(w\) by \(-\overline w\).
At \(0\), \(F_a\) stays bounded. At \(p_R\) and \(p_L\)
its real part tends to \(-\infty\) because some \(x_i\)
vanishes with \(a_i>0\). At \(1\) its real part tends
to \(+\infty\) due to \(-\log(1-q)\). Hence for each
height \(\eta\), \(\Re D(\xi+i\eta)\) decreases from
\(+\infty\) to \(-\infty\). Its zero defines a unique
smooth graph \(\xi=\xi(\eta)\), and the Cauchy–Riemann
equations give
\[
 \frac{d}{d\eta}\Im D(\xi(\eta)+i\eta)
   =\frac{|D'|^2}{\Re D'}<0.
\]

On the lower boundary the two real inverse branches have
\(\Im F_a=0,\pi\), respectively, by the logarithm of
\(1-q\). The same limits at the two real ends and the
strictly negative derivative of the real phase difference
give a unique equal-modulus point there. Local inversion
and the implicit function theorem extend the graph to this
point. Thus the phase difference at the lower endpoint
is \(i\pi\).

On the upper boundary, when \(\max\rho>1\), the positive real
coexistence point of the
[conditional integral](MULTIVARIATE_SHARED_SUPPORT_CONDITIONAL_SADDLE.md#a-local-exclusion-around-the-positive-real-axis)
has two distinct real maxima on opposite sides of \(x_0\).
At each \(X'<0\), so their inverse branches approach
from the upper half-plane when \(\Im w<0\). Their
phases \(F_a\) agree there since \(\Delta=0\), and the
derivative of the real phase difference is nonzero.
The implicit function theorem extends the graph to
that boundary point with \(D=0\).

When \(\max\rho=1\), take \(x_0\) to be its unique maximizer
and \(w_c=X(x_0)\) real. Strict concavity of \(\rho\) and
\(X'=(\rho-1)/(x(1-x))\) on the real interval give
\[
 X(q)=w_c-c(q-x_0)^3+O((q-x_0)^4),\qquad c>0.
\]
This has a holomorphic cubic coordinate obtained by taking
a local analytic cube root of the nonvanishing factor.
For \(w=w_c+i\eta\) with negative small \(\eta\), two
local inverse branches are in the upper half-plane and
converge to \(x_0\), so \(D(w)\to0\).
For small nonzero real \(h\), as \(w=w_c+h+i\eta\) approaches
the boundary from below, one branch tends to the real local
inverse and the other to the upper nonreal inverse. The
real inverse continues from above because \(X'<0\) there.
The cubic roots have real parts separated by
\(3(|h|/c)^{1/3}/2+O(|h|^{2/3})\), so the limiting
\(\Re D'(w_c+h+i\eta)\) is strictly negative.
Integrate \(\Re D'\) along a short horizontal segment
from \(w_c+i\eta\), and let \(\eta\to0\). The local inverses
are bounded on the segment, so dominated convergence applies.
The limiting real parts of \(D\) have opposite strict
signs at \(w_c-h\) and \(w_c+h\) for small \(h>0\).
Since \(\Re D\) decreases, \(\xi(\eta)\to w_c\).
The two local branches then coalesce, giving \(D\to0\)
also along the graph.

The imaginary part on the open graph consequently lies
strictly between zero and \(\pi\). This excludes
\(\exp D=1\) and proves the equal-degree assertion.
