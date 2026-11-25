# Two sheets of the stationary map for shared supports

Let \(a,b\in\mathbb Z_{>0}^r\) be distinct, delete the rows
with \(d_i=b_i-a_i=0\), and put \(x_i(q)=a_i+d_iq\).
On the upper half-plane \(\mathbb H\), use principal logarithms in
\[
 X(q)=\sum_i d_i\log x_i(q)+\log(1-q)-\log q .
\]
The right-hand side is single-valued there because each logarithm
has its argument in one of the two open half-planes. In the
[conditional integral](MULTIVARIATE_SHARED_SUPPORT_CONDITIONAL_SADDLE.md),
\(e^{-X(q)}\) is the stationary value of the external parameter
\(\zeta\).

**Theorem.** Put \(\mathcal T=\{w:-\pi<\Im w<0\}\).
The restriction \(X:X^{-1}(\mathcal T)\cap\mathbb H\to\mathcal T\)
is a proper holomorphic map of degree two. For positive real \(x\)
between zero and one, define
\[
 \rho(x)=x(1-x)\sum_i\frac{d_i^2}{a_i+d_ix}.
\]
If \(\max_{0<x<1}\rho(x)<1\), the preimage in \(\mathbb H\) is
connected and has exactly one simple critical point. If the maximum
is at least one, there are two preimage components, each mapped
conformally onto \(\mathcal T\).

For an external parameter in the upper half-plane this classifies
the solutions of \(X(q)=-\operatorname{Log}\zeta\) in \(\mathbb H\).
Where the conditional contour stays in its strip and \(s\) has
positive imaginary part, these are the only possible stationary
points. Determining the contributing contours at these points is
a separate step.

## Argument principle and properness

Write \(p_i=-a_i/d_i\). Positive gaps have poles to the left of
zero and negative gaps have poles to the right of one. Let
\(p_L=\max_{d_i>0}p_i\) and \(p_R=\min_{d_i<0}p_i\), using
\(-\infty,+\infty\) when the corresponding lists are empty. Along
the real axis away from the listed poles and \(0,1\), the limiting
imaginary parts from above are
\[
 \Im X(x+i0)=
 \pi\left(\sum_{\substack{d_i>0\\x<p_i}}d_i+
           \sum_{\substack{d_i<0\\x>p_i}}|d_i|
           -\mathbf1_{\{x<0\text{ or }x>1\}}\right).
 \tag{1}
\]
Every term \(d_i\Im\log x_i(q)\) is positive in
\(\mathbb H\), and \(0<\arg q-\arg(1-q)<\pi\), since
\(q/(1-q)\) lies in \(\mathbb H\). Thus \(\Im X>-\pi\)
throughout \(\mathbb H\).

Fix \(w=\xi+i\eta\) in \(\mathcal T\). Use the positively oriented
boundary of a large upper half-disk with small upper half-disks
removed around the real singularities. The straight portions can
be interpreted as upper limits, or displaced slightly into
\(\mathbb H\). The winding number of \(X-w\) is the signed number
of crossings of the horizontal ray \(\Im X=\eta,\ \Re X>\xi\),
with upward crossings positive. By (1) the straight portions
do not cross. On each small semicircle about a singularity \(b\),
the real part of \(X\) uniformly tends to either \(+\infty\) or
\(-\infty\) from its nonzero logarithmic coefficient. Its imaginary
part is affine in the circle angle up to a term whose value and
angular derivative tend uniformly to zero. Thus a crossing
occurs exactly when the limiting endpoint heights straddle
\(\eta\), provided the real part tends to \(+\infty\).

At \(0\), \(-\log q\) gives \(+\infty\) with heights increasing
from \(-\pi\) to zero, hence one crossing. At \(1\) and at a
finite \(p_L\), the real part tends to \(-\infty\). At a finite
\(p_R\) it tends to \(+\infty\) with heights increasing from
\(-\pi\) to a nonnegative multiple of \(\pi\), hence one
crossing. All other poles have nonnegative endpoint heights.

For \(q=R e^{i\theta}\), \(0\le\theta\le\pi\), the large
semicircle expansion, with an angularly differentiable uniform
error of order \(R^{-1}\), is
\[
 X(q)=\Delta\log R+C+
 i\left(\Delta\theta+\pi(D_--1)\right)+O(R^{-1}),\qquad
 \Delta=\sum_i d_i,\quad D_-=\sum_{d_i<0}|d_i|,
 \tag{2}
\]
where \(C\) is real. This follows by factoring out \(q\) in
each logarithm, with \(\log(-q)=\log q-i\pi\) on the upper
semicircle. If both signs occur, both limiting endpoint heights
are nonnegative, and the affine imaginary part has no crossing.
If every gap is positive, the height increases from \(-\pi\)
to a nonnegative multiple of \(\pi\) and \(\Re X\to+\infty\).
This gives the second crossing in the absence of \(p_R\).
If every gap is negative, \(\Re X\to-\infty\), so the count
is already supplied by \(0,p_R\). The total winding is two
in all cases. The argument principle counts the preimages with
multiplicity.

To see properness, consider a sequence leaving every compact
subset of \(\mathbb H\), with \(X\) in a compact subset of
\(\mathcal T\). Any finite real limit away from the singularities
has one of the boundary heights in (1). At a singularity the
real part escapes. At infinity the real part escapes if
\(\Delta\ne0\). If \(\Delta=0\) then both signs occur and
(2) tends to a height \(\pi(D_--1)\ge0\). All cases contradict
compactness inside \(\mathcal T\).

Whenever \(q\) is the image of an \(s\) in the upper part of the
conditional strip, \(\Im s\in(0,\pi/2)\) and
\(\Im X=\Im s-\arg(q/(1-q))\in(-\pi,\pi/2)\).
If \(e^{-X}\) is in the upper half-plane, the only lift in
this range is \(X=-\operatorname{Log}\zeta\). This verifies
the application to the conditional integral.

## Connectedness and ramification

The [critical point count for the conditional map](MULTIVARIATE_SHARED_SUPPORT_CONDITIONAL_SADDLE.md#stationary-equation-and-its-critical-points)
shows that \(X'\) has exactly one simple zero in \(\mathbb H\) if
\(\max\rho<1\), and has no zero there if \(\max\rho\ge1\).
In the second case a proper unramified map onto a simply connected
strip restricts to a conformal map on every component.
The degree two gives exactly two components.

For \(\max\rho<1\), we show that both preimages of one target
belong to the same component. The real intervals \((0,1)\) and
\((1,p_R)\) have boundary heights zero and \(-\pi\), respectively.
On the first interval \(X'=(\rho-1)/(x(1-x))<0\) and
\(X\) decreases from \(+\infty\) to \(-\infty\). On the
second interval
\[
 X'=\sum_i d_i^2/x_i-1/x+1/(x-1)>0
\]
and \(\Re X\) increases from \(-\infty\) to \(+\infty\).
The far endpoint limit follows from the logarithm at \(p_R\),
or from (2) when \(p_R=+\infty\), in which case all gaps are
positive. The stated derivative signs give curves in
\(X^{-1}(\mathcal T)\) immediately above every compact
subinterval of each interval. These curves on either side of
one are in the same component. Indeed, for small \(\epsilon>0\)
on \(q=1+\epsilon e^{i\theta}\), \(0<\theta<\pi\), the
imaginary part is \(\theta-\pi+O(\epsilon\sin\theta)\).
The error bound follows from the analyticity and real Taylor
coefficients of \(X(q)-\log(1-q)\) near one. Hence the open
semicircle stays in the preimage and connects the two curves
near its endpoints.

For sufficiently large \(\xi>0\), one preimage of
\(\xi-i\pi/2\) lies near zero. In the expansion
\(X=-\log q+B(q)\), the holomorphic \(B\) is real on the
positive real axis near zero. The conformal coordinate
\(q e^{-B(q)}\) lifts the segment
\(\xi-i\pi/2\) to \(\xi\) to a path ending on \((0,1)\).
Its interior is in the upper half-plane since the exponential
coordinate is upper there and the inverse is real and increasing
on a small real interval. The other preimage has a corresponding
path to the boundary at height \(-\pi\) along
\(\xi-i\pi/2\) to \(\xi-i\pi\). At a finite right pole use
\[
 X=-D_R\log(p_R-q)+B_R(q)-i\pi,\qquad
 D_R=\sum_{p_i=p_R}|d_i|,
\]
where \(B_R\) is holomorphic near \(p_R\) and real on the
real interval to its left. The coordinate
\((p_R-q)e^{-B_R(q)/D_R}\) inverts the exponential of
\(-(X+i\pi)/D_R\) and gives the path to \((1,p_R)\).
For \(p_R=+\infty\), (2) gives
\(X=\Delta\log q+C-i\pi+B_\infty(1/q)\) with
\(B_\infty\) holomorphic and real near \(0\).
Inverting the conformal coordinate
\(q^{-1}e^{-B_\infty(1/q)/\Delta}\) gives the same path
to a large positive real \(q\).
In each case the coordinate on the right is real with
negative derivative in \(q\), and its lower half-plane values
lift to the upper half-plane. The two preimages near distinct
endpoints are therefore in the component already connected
through \(1\).

The restriction of a proper holomorphic map to a preimage
component is proper and surjective onto the strip: its image
is both open and closed. Both preimages of the selected
target lie in one component, so there is no other component.
A connected degree-two proper map to a simply connected
strip must ramify. The critical point count identifies the
ramification as just the unique simple zero of \(X'\) in
\(\mathbb H\).

## Equal-modulus arcs for the two branches

On any simply connected region of \(\mathcal T\) away from
ramification, label the two inverse branches \(q_1(w),q_2(w)\).
Put \(F(q)=\sum_i a_i\log x_i(q)-\log(1-q)\), and use the
difference of the stationary phases
\[
 D(w)=F(q_2(w))-\Delta q_2(w)-F(q_1(w))+\Delta q_1(w).
\]
Differentiating gives
\[
 (F(q)-\Delta q)'=-qX'(q),\qquad D'(w)=q_1(w)-q_2(w).
\]
Thus \(D'\) is nowhere zero in the region. The locus
\(\Re D=0\) consists locally of smooth real-analytic arcs,
and \(\Im D\) is a strict coordinate on each connected
component in a region with fixed branch labels. This is the
equal-modulus locus when both saddles are in the contour
domain and contribute to the asymptotic expansion.

In a region with an exhaustive two-saddle expansion, the
[consecutive phase condition](MULTIVARIATE_SHARED_SUPPORT_CONDITIONAL_SADDLE.md#comparing-the-two-terminal-orders)
becomes
\[
 \exp(D(w)-\Delta D'(w))=1
\]
at a limiting common zero. For equal total degrees this
requires \(\Im D\in2\pi\mathbb Z\) along the equal-modulus
arc. For unequal total degrees it also requires
\(\Re D'=0\). At such a point the equal-modulus arc has a
horizontal tangent in the \(w\)-plane, since its defining
gradient is \((\Re D',-\Im D')\). These conditions only use
the two inverse branches and do not rely on a proposed
choice of contributing contours.
