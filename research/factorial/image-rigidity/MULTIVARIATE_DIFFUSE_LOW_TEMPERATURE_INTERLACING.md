# Consecutive interlacing in the low-temperature Curie–Weiss model

For \(s>0\), write
\[
 B_n(z;s)=\sum_{k=0}^n\binom nk e^{-s k(n-k)/(2n)}z^k.
\]
These are the polynomials in the [diffuse limit of the multivariate
problem](MULTIVARIATE_FIXED_TERMINAL_CURVATURE.md#the-limiting-pair-when-each-relative-gap-tends-to-zero).

**Theorem.** If \(s\ge64\), every consecutive pair \(B_n,B_{n+1}\)
has simple, strictly interlacing unit-circle roots. For every
\(\delta>0\) there is an integer \(n_\delta\) such that the same
conclusion holds for \(n\ge n_\delta\) and \(s\ge4+\delta\),
uniformly in \(s\).

Order the arguments in \((0,2\pi)\) as \(\theta_{n,1}<\cdots<
\theta_{n,n}\). Strict interlacing means
\(\theta_{n+1,k}<\theta_{n,k}<\theta_{n+1,k+1}\) for every \(k\).

## The phase comparison

We use a decomposition
\[
 B_n(z;s)=I_n(z)+z^n I_n(1/z),\qquad I_n=e^{H_n},
\]
where \(H_n\) is single-valued and analytic near the unit circle,
and is real at \(z=1\) and satisfies conjugation symmetry. At a
unit-circle point, put
\[
 T_n(\theta)=n\theta-2\Im H_n(e^{i\theta}).
\]
Then \(B_n(e^{i\theta};s)=0\) exactly when \(T_n(\theta)\)
is an odd multiple of \(\pi\). If \(T_n\) is strictly increasing,
its endpoints \(0,2n\pi\) give all \(n\) roots of \(B_n\).
The inequality \(0<T_{n+1}(\theta)-T_n(\theta)<2\pi\) for
\(0<\theta<2\pi\) then gives the required interlacing.

## A uniform estimate at every degree for \(s\ge64\)

Put \(a=s/2\). The Gaussian integral, expanded by the binomial
formula, gives
\[
 B_n(z;s)=\sqrt{\frac n{4\pi a}}\int_{\mathbb R}
       \bigl[e^{-x^2/(4a)}(1+ze^{x-a})\bigr]^n\,dx.
\]
Indeed the integral of the term at index \(k\), divided by
\(\binom nk z^k\), is \(e^{-ak+ak^2/n}\). Define \(I_n(z)\)
as the same integral over \((-\infty,a]\), with the same
prefactor. The substitution \(x\mapsto2a-x\) on the other
half gives the exact decomposition.

Use \(|z|\le2\), \(a\ge32\), and the phase
\(F(x,z)=-x^2/(4a)+\log(1+ze^{x-a})\) on
\(\Re x\le a/2\), with its logarithmic series. The saddle
equation is
\[
 x_*=\frac{2aze^{x_*-a}}{1+ze^{x_*-a}}.
\]
On \(|x|=1\) the right-hand side has modulus at most
\(8a e^{-a/2}<1\). Rouché's theorem gives a unique simple
saddle \(|x_*|<1\), analytic in \(z\) on a neighborhood of
the closed disk. The saddle value \(F_*=F(x_*,z)\) satisfies
\[
 |F_*|\le\frac1{4a}+4e^{-a/2}<\frac1{100}.
\]
These bounds use \(|ze^{x-a}|\le2e^{-a/2}\) on the saddle
disk and \(e^{16}>2^{16}\) at \(a=32\).
The functions \(a e^{-a/2}\) and \(e^{-a/2}\) decrease
throughout this range.

Move the contour to the horizontal line through \(x_*\), and
include the vertical connector at \(a\). The integrand before
introducing the logarithm is entire in \(x\). The left closing
side tends to zero by Gaussian decay. On the horizontal part up
to \(\Re x=a/2\), the inequality
\[
 \left|F_{xx}(x,z)+\frac1{2a}\right|
  =\left|\frac{ze^{x-a}}{(1+ze^{x-a})^2}\right|
  \le4e^{-a/2}
\]
gives
\[
 F(x_*+t,z)=F_*-\frac{t^2}{4a}+E(t),\qquad
 |E(t)|\le\frac{t^2}{400a}.
\]
Here \(t\) is real and the last inequality follows from
\(800a e^{-a/2}<1\) at \(a\ge32\). Using
\(e^v-1=v\int_0^1e^{uv}\,du\) and integrating the Gaussian
majorant, the normalized error on this portion after
division by \(e^{nF_*}\) is at most
\[
 \sqrt{\frac n{4\pi a}}\int_{\mathbb R}
 \frac{nt^2}{400a}e^{-99nt^2/(400a)}\,dt
 =\frac{1/100}{2(99/100)^{3/2}}<\frac1{100}.
\]
The right Gaussian completion starts at \(t\ge a/2-1\).
After substituting \(v=\sqrt{n/(4a)}\,t\) the lower endpoint
is at least \(b=15\sqrt{an}/64>5/4\), with \(b^2>7/4\).
Its normalized mass is less than
\(e^{-b^2}/(2b\sqrt{\pi})<(4/15)e^{-7/4}<1/20\).
The last bound follows from
\(e^{7/4}>(8/3)(1+3/4+(3/4)^2/2)=65/12\).

On the remaining horizontal part and the connector,
\(a/2\le\Re x\le a\) and \(|\Im x|<1\). The modulus of
the base of the integrand is at most
\[
 e^{1/(4a)}\bigl(e^{-a/16}+2e^{-a/4}\bigr).
\]
The first term bounds \(e^{-(\Re x)^2/(4a)}\).
The second follows since the exponent
\(\Re x-a-(\Re x)^2/(4a)\) increases to \(-a/4\) at \(x=a\).
After division by \(|e^{F_*}|\) this bound is at most
\(e^{-a/17}\), since
\[
 \frac1{4a}+\frac1{100}+2e^{-3a/16}
 \le\frac1{128}+\frac1{100}+\frac1{32}<\frac{2}{17}
 \le\frac a{272}.
\]
The combined length is at most \(a/2+1\le17a/32\), so the
normalized cost is at most
\[
 \frac{17}{32}\sqrt{\frac{an}{4\pi}}e^{-an/17}
 \le\frac{85}{96}e^{-32/17}<\frac7{50}.
\]
For the middle inequality, \(\sqrt{x}e^{-x/17}\) decreases
for \(x\ge32\) and \(\sqrt{32/(4\pi)}<5/3\).
For the last comparison,
\(32/17>15/8\) and
\[
 e^{15/8}>
 \frac{27}{10}\left(1+\frac78+\frac{(7/8)^2}{2}
                               +\frac{(7/8)^3}{6}\right)
 >\frac{85\cdot50}{96\cdot7}.
\]

Combining the three errors gives
\(I_n=e^{nF_*}(1+\eta_n)\) with \(|\eta_n|<1/5\) for
\(|z|\le2\). Therefore
\(H_n=nF_*+\log(1+\eta_n)\) is analytic there and has
conjugation symmetry. The logarithmic series gives
\(|\log(1+\eta_n)|\le\log(5/4)<23/100\). Cauchy's estimate
on radius-one disks gives, for \(|z|\le1\),
\[
 |H_n'|<\frac n{100}+\frac{23}{100}<\frac n4,\qquad
 |H_{n+1}'-H_n'|<\frac{47}{100}.
\]
Thus \(T_n'>0\). The imaginary part of \(H_{n+1}-H_n\)
vanishes at \(\theta=0,2\pi\), and its magnitude at
\(e^{i\theta}\) is less than
\((47/100)\min(\theta,2\pi-\theta)\). It follows that
\(0<T_{n+1}(\theta)-T_n(\theta)<2\pi\). This proves the
first assertion.

## Uniform eventual interlacing for \(s\ge4+\delta\)

It suffices to work on a compact interval \(2<a_-\le a\le a_+\)
for \(a=s/2\), because the first assertion covers \(a\ge32\).
Use the centered coordinate \(t=x-a\) and the entire base
\[
 G_a(t,z)=e^{-a/4-t^2/(4a)}(e^{-t/2}+ze^{t/2}).
\]
Then \(B_n(z;s)=\sqrt{n/(4\pi a)}\int_{\mathbb R}G_a(t,z)^n\,dt\).

We construct the left saddle on \(|z|=1\). For \(z=e^{i\theta}\)
with \(0\le\theta\le\pi\), write it as \(t_L=-r+iy\).
Set \(u=\theta+y\), \(C=\cos u\). The stationary equations are
\[
 P_C(r):=\frac{r(\cosh r+C)}{\sinh r}=a,\qquad
 y=\frac r{\sinh r}\sin u .
\]
For fixed \(C\in[-1,1]\), \(P_C\) strictly increases from
\(1+C\) to infinity on \(r>0\). Its derivative is affine
in \(C\). The endpoint functions \(r\tanh(r/2)\) and
\(r\coth(r/2)\) are strictly increasing. For the latter,
the derivative has the sign of \(\sinh r-r>0\).
Thus the first equation has a unique positive solution
for each \(0\le u\le\pi\) and \(a>2\).
This solution increases with \(u\), since \(P_C(r)\)
increases in \(C\). The function \(r/\sinh r\) is decreasing
and is strictly less than one. Consequently
\(u-(r/\sinh r)\sin u\) has strictly positive derivative
and increases from \(0\) to \(\pi\). This constructs \(t_L\)
for every \(\theta\), with negative real part uniformly
separated from zero on the compact \(a\)-interval. For
negative \(\theta\) use conjugation.

On the horizontal line through \(t_L\), the real logarithm
of the modulus is, for \(x=\Re t\), a constant plus
\[
 -\frac{x^2}{4a}+\frac12\log(\cosh x+C).
\]
On the left half-line \(x=-r'\) with \(r'>0\), its derivative
in \(r'\) is \(r'(1/P_C(r')-1/a)/2\).
It has a unique strict maximum at \(r'=r\), with second
derivative \(-r P_C'(r)/(2a^2)<0\). The endpoint \(x=0\)
has smaller modulus, including when \(C=-1\) and that
modulus is zero. The reflected saddle \(t_R=-\overline{t_L}\)
lies on the same horizontal line and has the mirror
maximum on the right half-line.

The nonzero stationary value and strictly negative real
curvature imply that the saddles are simple and extend
analytically in \(z\) to a common annulus about \(|z|=1\).
For detail, the analytic implicit function theorem gives
a locally unique root of the stationary equation near
each left saddle. Compactness in \(a,\theta\) gives uniform
neighborhood sizes after shrinking to keep the root isolated
from the other stationary points. Local continuations agree
on overlaps by this isolation and the continuity of the
chosen left saddle on the circle. The annulus can be chosen
invariant under \(z\mapsto1/z\), and we set
\(t_R(z)=-t_L(1/z)\) there. The symmetry
\(G_a(-t,1/z)=G_a(t,z)/z\) shows that it is the right saddle.
On the circle conjugation makes this definition agree
with the reflected saddle above.

Let \(M(z)=(t_L(z)+t_R(z))/2\). Define \(I_n(z)\) as the
integral of \(G_a(t,z)^n\) over a contour from the negative
real infinity to \(M(z)\), with the original prefactor.
The finite part is path independent because the integrand
is entire. Since \(M(1/z)=-M(z)\), the piece from \(M(z)\)
to positive real infinity is exactly \(z^n I_n(1/z)\).
This gives the exact decomposition with analytic \(I_n\).
On the circle \(M(z)\) lies on the imaginary axis at the
height of both saddles.

For the estimate on \(I_n\), use the horizontal line through
\(t_L\) up to \(\Re M\), with the final vertical connector
to \(M\). On the circle the connector has zero length and
the horizontal part has the strict maximum established above.
Outside any fixed small saddle neighborhood, the remaining
contour has a strict gap in log modulus on bounded pieces.
For the unbounded left tail, use, for bounded \(|z|\) and
\(|\Im t|\), the uniform bound
\[
 |G_a(t,z)|\ \le\ C_0
       \exp\!\left(-\frac{(\Re t)^2}{4a_+}+\frac{|\Re t|}{2}\right)
\]
with a constant \(C_0\) depending only on the compact parameter
range. The saddle modulus is bounded away from zero. After
division by it the far-left bound can thus be reduced to
\(e^{-c(\Re t)^2}\) for some \(c>0\). These bounded-piece
and tail estimates persist uniformly on a smaller annulus
by continuity. In particular the final vertical connector,
which starts at the strict-loss endpoint on the circle,
lies entirely in a strict-loss region on the smaller annulus.

Write \(F_*=\log G_a(t_L,z)\), using
\[
 F_*=-a/4-t_L^2/(4a)-t_L/2+\log(1+ze^{t_L}).
\]
The last logarithm has a single analytic branch by its series,
since \(|ze^{t_L}|<1\) on the circle uniformly in the parameters.
Let \(A_*=-\partial_t^2\log G_a(t_L,z)\), with
\(\Re A_*>0\) on this annulus after shrinking it.
On a fixed small horizontal segment through \(t_L\) the phase
equals \(F_*-A_*t^2/2+O(|t|^3)\) in the displacement \(t\).
Uniform bounded derivatives and the positive real curvature
give a Gaussian comparison error bounded by
\[
 C_1\sqrt n\int_{\mathbb R}n|t|^3 e^{-c_1nt^2}\,dt=O(n^{-1/2}),
\]
in units of \(e^{nF_*}\), with constants \(C_1,c_1>0\).
Gaussian completion costs an exponentially small error.
The finite pieces away from the saddle contribute
\(O(\sqrt n e^{-c_2 n})\), and the far-left Gaussian bound
gives the same estimate or better. All constants can be
chosen uniformly on the smaller annulus and compact \(a\)-range.
We obtain
\[
 I_n(z)=e^{nF_*}(2a A_*)^{-1/2}(1+\eta_n(z)),\qquad
 \eta_n=O(n^{-1/2}).
\]
The amplitude uses the square root with positive real part.
The quotient defining \(\eta_n\) is analytic, so its
\(z\)-derivative is also \(O(n^{-1/2})\) on a still smaller
annulus by Cauchy's estimate.

For large \(n\) this gives a single-valued logarithm
\(H_n=nF_*-\tfrac12\log(2a A_*)+\log(1+\eta_n)\).
All three terms have conjugation symmetry and are real at
\(z=1\) by their chosen branches. By stationarity,
\[
 z\partial_z F_*=\frac{ze^{t_L}}{1+ze^{t_L}}.
\]
In particular
\[
 \alpha(a,\theta):=1-2\Re(z\partial_z F_*)
      =\frac{1-|ze^{t_L}|^2}{|1+ze^{t_L}|^2}>0
       \qquad (z=e^{i\theta}).
\]
It has a uniform positive minimum on the compact parameter
range. The phase derivatives therefore satisfy
\[
 T_n'(\theta)=n\alpha(a,\theta)+O(1),\qquad
 (T_{n+1}-T_n)'(\theta)=\alpha(a,\theta)+O(n^{-1/2}).
\]
The amplitude logarithm cancels in the difference. Both
derivatives are positive for large \(n\), uniformly in \(a\).
The phase difference has endpoints \(0,2\pi\). The phase
comparison now proves strict interlacing throughout the
compact range, and the first assertion supplies \(a\ge32\).
\(\square\)

## Localization of possible common roots to the edge

Let \(a=s/2\) and define
\[
 \theta_*(a)=
 \begin{cases}
 \arccos(a-1)-\sqrt{1-(a-1)^2},&0<a\le2,\\
 0,&a\ge2 .
 \end{cases}
\]
Fix \(0<a_-\le a_+<\infty\) and \(\varepsilon>0\). There is
an integer \(n_0\) such that for \(a\in[a_-,a_+]\) and
\(n\ge n_0\), any common root of \(B_n,B_{n+1}\) has
principal absolute argument \(\theta\in[0,\pi]\) satisfying
\[
 |\theta-\theta_*(a)|<\varepsilon .
\]
In this assertion roots are on the unit circle by the
[unitary Hermite circle theorem](https://doi.org/10.5802/ahl.227).
The argument estimates below concern the two closed regions
\(\theta\le\theta_*(a)-\varepsilon\) and
\(\theta\ge\theta_*(a)+\varepsilon\) inside the compact parameter
rectangle \([a_-,a_+]\times[0,\pi]\).

On the first region one has \(a<2\). Put \(u_c=\arccos(a-1)\).
The function \(u\mapsto u-a\tan(u/2)\) increases strictly from
zero to \(\theta_*(a)\) on \(0\le u\le u_c\), since its
derivative is \(1-a/(1+\cos u)\), positive up to the last
endpoint. Solve \(\theta=u-a\tan(u/2)\) and put
\(t_0=i a\tan(u/2)\). This is a stationary point,
with \(u=\theta+\Im t_0\). On the horizontal line through it
the log modulus again has derivative in \(r'=|\Re t|>0\)
equal to \(r'(1/P_C(r')-1/a)/2\), with \(C=\cos u\).
Here \(P_C(0)=1+C>a\), so the modulus decreases strictly
away from \(t_0\). The real logarithmic curvature there is
\[
 -\frac{1}{2a}+\frac{1}{2(1+C)}<0.
\]
The saddle value is nonzero. On the compact first region
the curvature and saddle modulus have uniform positive
margins in absolute value. Shift the entire real line to
this horizontal line, with both closing sides vanishing by
Gaussian decay. The same local Taylor estimate as above
and strict loss off the saddle neighborhood, now on both
sides, give a nonzero Gaussian leading term with relative
error \(O(n^{-1/2})\) uniformly in the region. In particular
\(B_n(e^{i\theta};s)\) itself is nonzero there for large \(n\).

On the second region, when \(a\le2\) use \(u>u_c\) in the left
saddle construction. One has \(1+\cos u<a\), so \(P_C(r)=a\)
again has a unique solution \(r>0\). As \(u\) decreases to \(u_c\)
this solution tends to zero. The map
\[
 u\longmapsto u-\frac{r}{\sinh r}\sin u
\]
increases strictly from \(\theta_*(a)\) to \(\pi\). Its
derivative is positive by the same monotonicity of \(r\)
and \(r/\sinh r\) used above. This includes \(a=2\) with
\(u_c=0\). For \(a>2\) the earlier left saddle construction
applies on the whole arc. On the compact second region,
\(r\) has a positive minimum and the horizontal maxima
have uniformly nonzero curvature. The right saddle is
the reflected one at the same height. The horizontal
maximum argument on each half-line and the strict loss
at their midpoint therefore hold with uniform margins.

Use analytic continuations on neighborhoods of these
arcs and their conjugates, and the same split at the
analytic midpoint \(M(z)\). The compactness and implicit
function argument above gives a uniform neighborhood
width, with \(t_R(z)=-t_L(1/z)\). At the real endpoint
\(z=-1\), the two arcs have consistent continuations by
conjugation and simplicity of the saddles. On each such
neighborhood the same horizontal and final vertical
contours give the two-packet estimate above. The gap
on the finite contour pieces persists by continuity, and
the Gaussian bound on the left tail is uniform. All
logarithms in
\(H_n=nF_*-\tfrac12\log(2a A_*)+\log(1+\eta_n)\) use the
same branches as above. In particular the difference of
phase functions satisfies, on the compact second region,
\[
 T_{n+1}(\theta)-T_n(\theta)
   =\theta-2\Im F_*(e^{i\theta})+O(n^{-1/2})
\]
uniformly in \(a,\theta\).

For \(a\le2\), as \(\theta\) decreases to \(\theta_*(a)\)
the left saddle becomes \(i\sin u_c\). At that point
\(G_a\) is \(e^{i\theta_*(a)/2}\) times a positive real
factor, and the chosen logarithm has imaginary part
\(\theta_*(a)/2\). Indeed
\(1+e^{i u_c}=2e^{iu_c/2}\cos(u_c/2)\) with
\(\cos(u_c/2)>0\), and the same determination is the
limit from \(|ze^{t_L}|<1\). For \(a>2\) the phase begins
at zero because \(t_L\) and \(F_*\) are real at \(z=1\).
The function \(\theta-2\Im F_*\) thus starts at zero
at \(\theta_*(a)\) in both cases. Its derivative is
\(\alpha(a,\theta)>0\) for \(\theta>\theta_*(a)\)
by the stationary-value formula. It ends at \(\pi\)
since \(F_*\) is real at \(z=-1\). In particular it
has a positive minimum on the compact second region
and is bounded above by \(\pi\). For large \(n\) the
phase difference is strictly between zero and \(2\pi\)
there, which excludes a common root.

These estimates leave the region where the stationary
points coalesce at the edge. At \(a=2,\theta=0\) that region
contains the critical cusp.

## The bounded Airy window in the subcritical range

Fix a compact interval \([a_-,a_+]\subset(0,2)\) and a constant
\(L>0\). Uniformly for \(a\) in this interval and for all sufficiently
large \(n\), \(B_n\) and \(B_{n+1}\) have no common unit-circle root
with principal absolute argument in
\[
 \bigl|\theta-\theta_*(a)\bigr|\le L n^{-2/3}.
\]
The proof tracks the simple Airy roots using one differentiable
interpolation in the order. This distinguishes the corresponding
roots at consecutive orders even though their distance is smaller
than the leading Airy displacement.

Use \(u_c=\arccos(a-1)\), \(t_c=i\sin u_c\), and
\(\eta=\theta-\theta_*(a)\). For the real normalized function
\(e^{-in\theta/2}B_n(e^{i\theta};s)\) the phase is
\[
 f_a(\zeta,\eta)=
 -\frac a4-\frac{(t_c+\zeta)^2}{4a}
 +\log\!\left(2\cosh\frac{t_c+\zeta+i(\theta_*(a)+\eta)}2\right).
\]
The logarithm uses the branch real at \(\zeta=0\) for real \(\eta\)
near zero. In the strip where the imaginary part of
\(t_c+\zeta+i(\theta_*+\eta)\) lies in \((-\pi,\pi)\), the argument
of this logarithm has positive real part. The original real line
in \(t=t_c+\zeta\) and the horizontal line through \(t_c\) are
in the strip, with a uniform margin for small real \(|\eta|\)
on the compact \(a\)-interval. Their integrals, with the Gaussian
prefactor, agree for every positive real order in this branch.
This follows by deformation, with the closing sides vanishing
by Gaussian decay. At integer orders the integral is exactly the
normalized polynomial. For real \(\eta\) the horizontal integrand
has conjugation symmetry under \(\zeta\mapsto-\bar\zeta\)
and its integral is real.

At \(\zeta=\eta=0\), the first two \(\zeta\)-derivatives vanish.
With \(b=\tan(u_c/2)>0\), differentiation gives
\[
 f_{\zeta\eta}=\frac{i}{2a},\qquad
 f_{\zeta\zeta\zeta}=-\frac{ib}{2a},\qquad
 f_{\zeta\zeta\zeta\zeta}=-\frac{b^2}{2a}-\frac1{2a^2}.
\]
In particular, for real \(\eta=0\),
\(f_a(\zeta,0)-f_a(0,0)=-iA\zeta^3-D\zeta^4+O(|\zeta|^5)\)
with \(A=b/(12a)>0\) and \(D>0\), both bounded away from
zero uniformly in \(a\).

In the horizontal contour replace the segment \([-r_0,r_0]\)
by the rays of arguments \(-5\pi/6,-\pi/6\) between the
origin and radius \(r_0\), with circular arcs of radius \(r_0\)
joining their outer endpoints to the real line. Choose \(r_0>0\)
small enough uniformly that the contour stays in the logarithm
strip. On the right arc of argument \(\phi\in[-\pi/6,0]\)
the leading real terms at \(\eta=0\) are
\(A r_0^3\sin(3\phi)-D r_0^4\cos(4\phi)\).
For \(-\pi/12\le\phi\le0\) the first term is nonpositive
and \(\cos(4\phi)\ge1/2\). For \(-\pi/6\le\phi\le-\pi/12\)
the cubic term is strictly negative uniformly, dominating
the quartic remainder. This gives strict loss on the whole
arc after decreasing \(r_0\), and reflection gives the same
on the left arc. Along both rays the real part drops from
the value at zero by at least \(c|\zeta|^3\) at \(\eta=0\)
for some uniform \(c>0\). The retained horizontal pieces
also have strict loss away from the origin, since
\(P_{\cos u_c}(r)>a\) for \(r>0\). Infinite tails have the
Gaussian majorant from the regular-arc proof.

Set \(\lambda=n^{-1/3}\), and interpolate by taking \(\lambda>0\)
real, \(\eta=\lambda^2 v\). The parameter \(v\) is kept in
a fixed bounded real interval. Divide the contour integral
by \(\lambda e^{f_a(0,\lambda^2v)/\lambda^3}\), omitting the
nonzero Gaussian prefactor. Denote the resulting real function
by \(J(a,v,\lambda)\). On the rays, with \(\zeta=\lambda \xi\),
its exponent is
\[
 \Phi(a,v,\lambda,\xi)=
   \lambda^{-3}\bigl(f_a(\lambda\xi,\lambda^2v)
                   -f_a(0,\lambda^2v)\bigr).
\]
A mixed derivative bound gives an \(O(|\eta||\zeta|)\)
perturbation of the difference in the exponent from its
value at \(\eta=0\). Thus on the scaled rays
\(\Re\Phi\le-c|\xi|^3+C|v||\xi|\) uniformly.
The fixed-radius arcs and retained horizontal pieces
remain in a strict-loss region for small \(\lambda\).
After this normalization their integrals and their
derivatives in \(v,\lambda\) are exponentially small.
Indeed differentiation gives only powers of \(\lambda^{-1}\)
and at most polynomial growth in \(|\Re\zeta|\) on the
far tails, where the uniform Gaussian bound applies.

For completeness, the \(\lambda\)-derivative on the
scaled rays is uniformly integrable. Expand first in
\(\zeta\) at fixed \(\eta\) and use the three displayed
derivatives at the origin. One obtains
\[
 \begin{split}
 \Phi={}&i\frac{v\xi}{2a}-iA\xi^3\\
 &+\lambda\bigl(v\xi^2R_2(\lambda^2 v)
                  +\xi^4R_4(\lambda\xi,\lambda^2 v)\bigr)\\
 &+\lambda^2\bigl(v^2\xi R_1(\lambda^2 v)
                     +v\xi^3R_3(\lambda^2 v)\bigr),
 \end{split}
\]
where the four remainders and their first derivatives
are bounded uniformly for \(|\lambda\xi|\le r_0\) and
small \(|\lambda^2 v|\). This bounds
\(\partial_\lambda\Phi\) and \(\partial_v\Phi\) by a
constant times \(1+|\xi|^4\), using \(\lambda|\xi|\le r_0\)
when differentiating \(R_4\). The endpoint terms from
varying the ray lengths \(r_0/\lambda\) are exponentially
small. Dominated convergence and differentiation now give
uniform convergence of \(J,\partial_v J\) to their limiting
values, and a uniform bound on \(\partial_\lambda J\).

With \(c_a=(3A)^{1/3}\), the limiting ray integral is
\[
 J(a,v,0)=\int_{\Gamma}
       e^{\,i v\xi/(2a)-iA\xi^3}\,d\xi
       =\frac{2\pi}{c_a}\,
           \mathrm{Ai}\!\left(-\frac{v}{2a c_a}\right),
\]
where \(\Gamma\) consists of the two full rays in their
left-to-right orientation. To get the last equality use
\(w=i c_a\xi\), giving the usual Airy contour from
argument \(-\pi/3\) to \(\pi/3\).
Its real zeros are at \(v=2a c_a|a_k|>0\), where \(a_k<0\)
are the real Airy zeros. They are simple by the Airy
differential equation and uniqueness for its initial data.
[DLMF Airy contour integral](https://dlmf.nist.gov/9.5.E4).
[DLMF real Airy zeros](https://dlmf.nist.gov/9.9.i).

For every bounded \(v\)-interval only finitely many of
these Airy zero branches can enter it, uniformly in \(a\).
Their positive locations and pairwise spacings are bounded
away from zero on the compact parameter range. In disjoint
small neighborhoods of those branches, uniform convergence
of \(J,\partial_v J\) gives one real zero \(v_k(a,\lambda)\)
of \(J\), with \(|\partial_v J|\) bounded below. On the
remaining bounded \(v\)-region \(J\) is separated from zero
for small \(\lambda\). The implicit function theorem
and the bound on \(\partial_\lambda J\) give uniformly
bounded \(\partial_\lambda v_k\). Thus along each branch
\[
 \frac{d}{d\lambda}
  \left(\theta_*(a)+\lambda^2v_k(a,\lambda)\right)
  =2\lambda v_k(a,\lambda)+\lambda^2\partial_\lambda v_k(a,\lambda)
  >0
\]
for small positive \(\lambda\). The angles at \(n\) and \(n+1\)
cannot agree on the same branch. They cannot agree on
different branches in this bounded window either, since a
common angle would give two \(v\)-coordinates differing by
a factor \((1+1/n)^{2/3}\), while the scaled branches have a
uniform positive separation. Taking the bounded interval for
the convergence larger than \([-2L,2L]\) accounts for both
orders. This proves the claimed exclusion. \(\square\)

A one-saddle estimate controls the intermediate region below
the subcritical edge. Write \(h=u_c-u>0\) for the stationary
parameter below the edge. The derivative of \(u-a\tan(u/2)\)
vanishes at \(u_c\), and its second derivative there is
\(-\sin(u_c)/a\), bounded away from zero on the compact
subcritical interval. Thus, as \(h\to0^+\),
\(\theta_*(a)-\theta\) is comparable to \(h^2\), uniformly in \(a\).
The real negative curvature \(d\) of the phase at the horizontal
saddle \(t_0=ia\tan(u/2)\) is comparable to \(h\), by its explicit
formula
\[
 d=\frac{1}{2a}-\frac1{2(1+\cos u)}.
\]
The real part of the phase on this horizontal is even in the
real displacement \(x\), and its fourth derivative converges
to the strictly negative edge derivative computed above. Therefore
on some fixed \(|x|\le r_0\) its difference from the saddle value
has real part at most \(-c(hx^2+x^4)\) for a uniform \(c>0\).
The complex difference from its quadratic approximation is
\(O(|x|^3)\). Comparing the exponentials by
\(e^v-1=v\int_0^1e^{uv}\,du\), the relative Gaussian
replacement error is at most a constant times
\[
 \sqrt{nh}\int_{\mathbb R} n|x|^3e^{-cnh x^2}\,dx
   =O\bigl((nh^3)^{-1/2}\bigr).
\]
Here the normalized Gaussian amplitude is comparable to
\(h^{-1/2}\). Completion of the Gaussian has an exponentially
small relative error in \(nh\). The rest of the actual horizontal
has a uniform strict loss at fixed distance from the saddle,
including at the limiting edge, and the far tails have the same
quadratic majorant as before. This proves nonvanishing below
the edge when \(n^{2/3}(\theta_*(a)-\theta)\to+\infty\).
Together with the bounded Airy estimate and the regular arcs,
it excludes all arguments below the edge for large \(n\) on
compact subcritical intervals.

The next argument covers distances above the edge larger
than the bounded Airy window.

## Intermediate arguments above the subcritical edge

**Theorem.** On every compact interval \(0<s_-\le s\le s_+<4\)
there is an order \(n_0\) such that \(B_n,B_{n+1}\) are coprime
for all \(n\ge n_0\), uniformly in \(s\).

To complete the subcritical range, let \(a=s/2\) lie in a compact
subinterval of \((0,2)\) and \(\theta>\theta_*(a)\) be close to
the edge. In the two-saddle parametrization,
\[
 \cos u=\frac{a\sinh r}{r}-\cosh r
       =a-1-\frac{3-a}{6}r^2+O(r^4),
 \qquad y=\frac r{\sinh r}\sin u .
\]
Since \(a-1\) stays in a compact subset of \((-1,1)\), \(u,y\)
are analytic functions of \(r^2\) near zero. Expansion of
\(\theta=u-y\) gives
\[
 \theta-\theta_*(a)=\frac{2-a}{2\sin u_c}r^2+O(r^4).
\]
In particular this distance is comparable to \(r^2\) uniformly
on the compact interval.

After removing the phase \(e^{in\theta/2}\), center at \(iy\)
and use the phase
\[
 F_r(\zeta)=-\frac a4-\frac{(iy+\zeta)^2}{4a}
               +\log\!\left(2\cosh\frac{\zeta+iu}{2}\right).
\]
The logarithm is defined in the same strip as in the Airy
argument, with \(F_r(0)\) real. Let \(F_s=F_r(-r)\) denote
the left saddle value in this section. Write
\(A=\tan(u_c/2)/(12a)\) and
\[
 M_r(\zeta)=iA(3r^2\zeta-\zeta^3),\qquad
 E_r(\zeta)=F_r(\zeta)-M_r(\zeta).
\]
Because \(F_r\) depends analytically on \(r^2\), the cubic
coefficient at the edge and the vanishing quadratic term there
give
\[
 |E_r''(\zeta)|\le C(r^2+|\zeta|^2)
\]
on a fixed small disk about zero, uniformly in \(a\).
Both \(F_r'\) and \(M_r'\) vanish at \(-r\). Integrating the
last estimate on straight segments in the disk gives
\[
 |E_r(\zeta)-E_r(-r)|\le C|\zeta+r|^2(r^2+|\zeta|^2),
 \qquad
 |E_r(\zeta)-E_r(0)|\le C r^3|\zeta|
           \quad (|\zeta|\le r).
\]
For the second estimate use \(E_r'(-r)=0\) to bound \(E_r'\)
first on \(|\zeta|\le r\). In particular
\(\Im F_s=-2Ar^3+O(r^4)\), uniformly in \(a\).

Shift the full contour to the horizontal line through \(iy\)
and split at \(\zeta=0\). At integral order \(n\) the normalized
polynomial is twice the real part of the left integral, by
reflection. On the whole contour the same logarithm gives a
continuous interpolation for real \(n>0\). Deform the left
half as follows, using \(w=\zeta/r\).

- Approach \(w=-1\) along \(w=-1-(1+i)\tau\), decreasing
  \(\tau\) from its intersection with a fixed circle
  \(|\zeta|=r_0\).
- Leave along \(w=-1+(1+i)\tau\) for \(0\le\tau\le1/2\).
- Approach the endpoint zero on \(w=(-1+i)\tau\), decreasing
  \(\tau\) from \(1/2\) to zero.

Join the beginning of the first segment to the old horizontal
at \(-r_0\) by a circular arc. The orientation is from the
negative real infinity to zero. The contour stays in the
logarithm strip for small fixed \(r_0\) and sufficiently small
\(r\), uniformly in \(a\). The old far-left part stays horizontal.

On the first two segments the real model losses from the saddle
are \(-Ar^3(6\tau^2+2\tau^3)\) and
\(-Ar^3(6\tau^2-2\tau^3)\), respectively.
The first perturbation estimate above adds at most
\(C r^4\tau^2(1+\tau^2)\). On the left segment its ratio to
the loss is bounded by a constant times \(r(1+\tau)\le r_0\),
and on the other segment by a constant times \(r\).
For \(r_0\) and \(r\) sufficiently small the real losses are
therefore at least \(c r^3(\tau^2+\tau^3)\) on the approach
and \(c r^3\tau^2\) on the departure, with \(c>0\) uniform.
On the endpoint segment the model's real part relative to
zero is \(-Ar^3\tau(3-2\tau^2)\).
The second perturbation estimate gives an error \(O(r^4\tau)\).
The exact phase at zero has smaller real part than \(F_s\) by
the horizontal maximum argument. Thus the loss relative to
\(F_s\) on the endpoint segment is at least \(c r^3\tau\).

On the circular connector the angles lie in
\([-\pi,-3\pi/4]\). At the limiting edge, the cubic and quartic
terms give strict loss uniformly on this arc. The proof is the
reflection of the Airy arc estimate, with the right arc range
enlarged to \([-\pi/4,0]\). Its portion away from zero argument
still has \(\sin(3\phi)\le-\sqrt{2}/2\) on
\([-\pi/4,-\pi/12]\), while the quartic term has negative real
part on \([-\pi/12,0]\). By continuity the loss persists at
small \(r\). The retained horizontal also has uniform strict
loss off a fixed saddle neighborhood, and the same quadratic
majorant at the far-left infinity.

Let
\[
 m=nr^3,\qquad
 \Psi_r(w)=\frac{F_r(rw)-F_s}{r^3},\qquad
 K(m)=\int_{\Gamma_r}e^{m\Psi_r(w)}\,dw,
\]
where \(\Gamma_r\) is the preceding contour in \(w\). For fixed
\(a,\theta\) this contour is independent of \(m\). Including the
prefactor, the left integral is
\(\sqrt{n/(4\pi a)}\,r e^{nF_s}K(m)\).
Along the local line \(w=-1+(1+i)\tau\), \(\Psi_r\) has a
quadratic term \(-L_r\tau^2/2\) with
\(L_r=12A+O(r)\), since \(F_r''(-r)=6iAr+O(r^2)\).
Its third derivative on a fixed small local interval is
uniformly bounded. The Gaussian comparison there gives
\[
 \sqrt m\,K(m)=C_r+O(m^{-1/2}),\qquad
 C_r=(1+i)\sqrt{\frac{2\pi}{L_r}},
                                                        \tag{6}
\]
uniformly for \(m\ge1\) as \(r\to0\). In this estimate the endpoint
segment contributes \(O(m^{-1})\) before multiplication by
\(\sqrt m\), by its linear loss. The nonlocal diagonals have
exponential loss in \(m\). The connector and the original
horizontal tails have exponential loss in \(m/r^3\), which
absorbs their rescaled lengths. The amplitudes \(C_r\) are
bounded and bounded away from zero.

We also need the order derivative of the same normalized
integral. Differentiation with respect to \(m\) gives
\[
 \frac{d}{dm}(\sqrt m\,K(m))
    =\frac{K(m)}{2\sqrt m}
        +\sqrt m\int_{\Gamma_r}\Psi_r(w)e^{m\Psi_r(w)}\,dw
    =O(m^{-1/2}).                                           \tag{7}
\]
For the local integral, \(|\Psi_r(w)|\le C|\tau|^2\) gives
\(O(m^{-1})\) after normalization. On the endpoint segment
\(\Psi_r\) is bounded and the linear loss gives \(O(m^{-1/2})\)
after normalization. On the nonlocal diagonals the model
and perturbation bound \(|\Psi_r|\) by a constant times
\(\tau^2+\tau^3\), integrable with the exponential loss.
On the connector and old horizontal, any powers of
\(r^{-1}\) and the quadratic growth of \(F_r\) at infinity
are absorbed by the exponential loss in \(m/r^3\).
These same bounds justify the differentiation.

Write \(\sqrt m K(m)=C_r(1+\epsilon_r(m))\). By (6) and (7),
\(\epsilon_r=O(m^{-1/2})\) and
\(\epsilon_r'=O(m^{-1/2})\). For \(m\) sufficiently large
the packet has a nonzero analytic-in-order amplitude, and
its change in argument from \(n\) to \(n+1\) can be computed
from the logarithm of \(1+\epsilon_r\). The real prefactor
\(\sqrt{n/(4\pi a)}\,r/\sqrt m\) is independent of \(n\).
The change is therefore
\[
 \Im F_s+O(r^3 m^{-1/2})=-2Ar^3+O(r^4+r^3m^{-1/2}).
\]
This lies strictly between \(-\pi\) and zero for small
\(r>0\) and large \(m\), uniformly in \(a\). If both
normalized polynomials vanished, their left packets would
both be nonzero and purely imaginary, requiring the change
in argument to be a multiple of \(\pi\). Thus the consecutive
pair has no common root here.

Since \(r^2\) is comparable to the distance above the edge,
bounded \(m=nr^3\) is included in a fixed bounded Airy window.
Large \(m\) is covered by the preceding packet estimate.
Below the edge the one-packet and Airy bounds already exclude
roots, and fixed distances from the edge are covered by the
regular-arc estimates. This completes the subcritical compact
range. The theorem does not pass through \(s=4\).
\(\square\)
