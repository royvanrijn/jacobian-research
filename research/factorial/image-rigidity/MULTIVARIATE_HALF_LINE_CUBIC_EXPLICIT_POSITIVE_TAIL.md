# Explicit positive-temperature bounds for the half-line cubic

Let
\[
 F(T,H)=\int_0^\infty e^{-y^3/3-Ty^2/2+Hy}\,dy .
\]
At a simple zero, let \(H_T\) denote the temperature derivative
of its zero branch.

**Theorem.** If \(T\ge10^6\), every zero is simple and every
upper zero satisfies \(\Im(2H-T H_T)>0\).
If \(T\ge10^{11}\), distinct zeros also have distinct imaginary
heights. These bounds are explicit, with no finite root count
required. This quantifies the positive tail in the
[half-line zero-separation theorem](MULTIVARIATE_HALF_LINE_CUBIC_ZERO_SEPARATION.md).

## Exterior fields

Take \(T\ge1\), put \(a=T/2\), and let \(R=|H|\ge20\sqrt T\).
At an upper zero \(\Re H>0\) by the decreasing-density sine
transform argument, and there are no real zeros.
Write \(Z=H+a^2\) and \(\omega=e^{2\pi i/3}\). The
[exact Airy connection](MULTIVARIATE_HALF_LINE_CUBIC_ZERO_SEPARATION.md#an-exact-airy-decomposition)
is \(F=E+S\), where
\[
 E=\omega F(T\omega^2,\omega H),\qquad
 S=2\pi i\omega^2e^{-aH-2a^3/3}\operatorname{Ai}(\omega^2 Z).
\]
The endpoint term is the integral of the original cubic
integrand along the ray of angle \(7\pi/9\), by rotation
from the ray of angle \(2\pi/3\). On this ray, the real
parts of the quadratic and cubic exponents are nonpositive.
The real linear exponent is at most \(-R u/2\), where
\(u\) is distance along the ray. Comparison with the linear
integral and direct differentiation thus give
\[
 \left|E+\frac1H\right|\le\frac1R
            \left(\frac{32}{R^3}+\frac{8T}{R^2}\right),
 \qquad |E_H|\le4/R^2,\qquad |E_T|\le8/R^3 .
\]
The relative error in this endpoint comparison is below
\(1/20\) for the stated range. In particular
\(|E_H/E|\le5/R\), \(|E_T/E|\le9/R^2\).

Set \(P=\sqrt Z-a\), and let \(\Xi=\log(S/E)\) locally.
The [explicit Airy-sector estimates](EXPLICIT_AIRY_CONTOUR_ENVELOPE.md#3-explicit-primary-airy-error-bounds-when-xige4)
give
\[
 \left|\omega^2\frac{\operatorname{Ai}'(\omega^2 Z)}
                         {\operatorname{Ai}(\omega^2 Z)}
                          -\sqrt Z\right|\le3/|Z|.
\]
Indeed \(Z\) is in the right quadrant and \(|Z|\ge R\ge20\),
so the relative errors of the Airy value and derivative
are bounded by \((25/24)|Z|^{-3/2}\) and
\((49/48)|Z|^{-3/2}\). In particular \(S\) is nonzero.
Differentiating the exact quotient gives
\[
             |\Xi_H-P|\le8/R,qquad
             |\Xi_T+P^2/2|\le3a/|Z|+9/R^2.              \tag{1}
\]

To localize the exterior zero angles, put
\(V=R/(a+\sqrt R)^{1/2}\). Then
\[
 V^2=\frac{R^2}{a+\sqrt R}\ge\frac{400}{1/2+\sqrt{20}}>80 ,
\]
using monotonicity in \(R\) and the lower bound
\(R\ge20\sqrt T\). The saddle exponent has the exact
representation
\[
 \mathcal S=\tfrac23 Z^{3/2}-aH-\tfrac23a^3
     =H^2\int_0^1\frac{1-u}{2\sqrt{a^2+uH}}\,du.
\]
If \(\theta=\arg H\), the integral after multiplication
by \(H^2\) has integrand arguments between
\(3\theta/2\) and \(2\theta\).
Thus \(\Re\mathcal S\ge V^2/8\) for \(\theta\le\pi/6\),
and \(\Re\mathcal S\le-V^2/12\) for \(\theta\ge5\pi/12\).
The Airy and endpoint bounds give
\[
       (V/2)e^{\Re\mathcal S}\le |S/E|
                             \le4V e^{\Re\mathcal S}.
\]
Here \(|Z|^{1/4}\) is between
\((a+\sqrt R)^{1/2}/\sqrt2\) and \((a+\sqrt R)^{1/2}\).
Since \(V>8\) and \(e^{V^2/12}>4V\) on that range,
every exterior zero lies in the open cone
\(\pi/6<\theta<5\pi/12\).

On the closed cone,
\(P=H/(\sqrt Z+a)\) gives
\(\Re P\ge V^2/(8R)\), since
\(|\sqrt Z+a|\le2(a+\sqrt R)\) and \(\cos(5\pi/12)>1/4\).
Thus \(\Re\Xi_H>0\) by (1).
At a zero this gives simplicity. On the segment between
two exterior zeros of the same imaginary height, the
modulus does not leave the exterior range and the angle
stays in the cone. There \(\log|S/E|\) is strictly increasing
horizontally, excluding the pair.

The same estimates give a strict imaginary order derivative.
At an exterior zero \(H_T=-\Xi_T/\Xi_H\). From (1),
\[
 \left|(2H-T H_T)-(2H-aP)\right|
 =\left|\frac{T(\Xi_T+P^2/2)+aP(\Xi_H-P)}{\Xi_H}\right|
 \le \frac{11R}{V^2/2-8}.
\]
Here \(T|\Xi_T+P^2/2|\le6+9T/R^2<7\),
\(a|P|\le R/2\), and
\(|\Xi_H|\ge |P|-8/R\ge(V^2/2-8)/R\).
Writing \(y=\Im H\), one has
\(\Im P=y/(2\Re\sqrt Z)\) and \(\Re\sqrt Z\ge a\).
Thus \(\Im(2H-aP)\ge3y/2>3R/4\).
The displayed error is less than \(11R/32\), proving
the strict positive sign in the exterior.

## Half-Gaussian comparison at finite scaled field

Write \(s=T^{-3/2}\), \(z=H/\sqrt T=x+iy\), and
\[
 G_s(z)=\sqrt T F(T,\sqrt T z)
       =\int_0^\infty e^{-u^2/2-su^3/3+zu}\,du,qquad
 G=G_0 .
\]
Use \(x,y\ge0\). For \(|z|\le K\), deform the contour to the
vertical segment from zero to \(iy\), followed by the
horizontal from \(iy\) to \(+\infty+iy\).
On the vertical the Gaussian real exponent is
\(v^2/2-yv\le0\) and the cubic is imaginary.
On the horizontal, at \(u=t+iy\), put \(v=t-x\).
The real Gaussian exponent is
\(\Re(z^2)/2-v^2/2\).
The cubic real exponent is at most \(s t y^2\).
When \(s K^3\le1\), its positive part is at most
\(1+|v|\). The closing segments at infinity vanish
for the integral and its derivatives.

In these estimates use \(|u|\le |v|+K\) on the horizontal
and \(e^{-v^2/2+1+|v|}\le e^2e^{-v^2/4}\).
For \(1\le d\le4\),
\[
 \int_{\mathbb R}(|v|+K)^d e^{-v^2/2+1+|v|}\,dv
                  \le 8\cdot2^{d-1}(48+4K^d).
\]
The vertical integral of \(|u|^d\) is at most
\(K^{d+1}/(d+1)\). Consequently, for \(K=20\),
\[
 \begin{split}
 |G_s-G|&\le 10^7s(1+|e^{z^2/2}|),\
 |\partial_s G_s|&\le10^7(1+|e^{z^2/2}|),\
 |G_s''|&\le10^7(1+|e^{z^2/2}|).
 \end{split}                                                   \tag{2}
\]
The first bound uses
\(|e^{-su^3/3}-1|\le (s/3)|u|^3e^{\max(0,-s\Re u^3/3)}\).
For \(K=42\) the same calculation gives both
\[
 |G_s-G|,quad |G_s'-G'|
           \le10^9s(1+|e^{z^2/2}|).                         \tag{3}
\]
Each bound is individual.

At a zero of \(G_s\) with \(|z|\le20\) and \(T\ge10^6\),
we have \(10^7s\le1/100\).
Since \(G(z)+G(-z)=\sqrt{2\pi}e^{z^2/2}\) and \(|G(-z)|<2\),
(2) shows \(|e^{z^2/2}|<2\) and \(|G(z)|<3/100\).
Integration by parts gives the exact relation
\(G_s'=1+zG_s-sG_s''\). Thus the zero is simple,
\(|G_s'|>1/2\), and its derivative in \(s\) satisfies
\(|s z_s|\le6\cdot10^7s\le 0.06\).

Also \(y>1/10\) at such a zero. Indeed \(G'=zG+1\);
integrating this equation along the segment from \(z\)
to \(x\) would give
\[
 |G(x)|\le e^{y^2/2}(|G(z)|+y)<1/2\qquad(0\le y\le1/10),
\]
contrary to \(G(x)>1/2\) for \(x\ge0\).
On the physical branch \(H=\sqrt T z(s)\) the order
derivative is \((3/2)\sqrt T(z+s z_s)\).
Its imaginary part is positive by the preceding bounds.
This proves the first assertion of the theorem together
with the exterior argument.

## Quantitative height gaps in the inner fields

For \(T\ge10^{11}\) we have \(s<10^{-16}\).
At any actual zero \(|z|\le40\), (3) and the Gaussian
connection give \(|e^{z^2/2}|<2\) and \(|G(z)|<3\cdot10^{-7}\).
There is a Gaussian zero within radius \(10^{-6}\) of
that point. To verify this by Rouché, use
\[
 G(w)=e^{(w^2-z^2)/2}G(z)+\int_z^w e^{(w^2-u^2)/2}\,du .
]
On the circle \(|w-z|=10^{-6}\) the difference of the
right side from \(w-z\) is less than half the radius.
The nearby Gaussian zero has modulus below \(41\) and is
still in the upper half-plane: actual zeros in the inner
range have \(y>1/10\), and exterior ones are in the proved
cone.

For completeness, all upper Gaussian zeros have real
part greater than one. They lie in the right half-plane
by the Gaussian Stieltjes formula in the
[half-line separation proof](MULTIVARIATE_HALF_LINE_CUBIC_ZERO_SEPARATION.md#all-zeros-at-large-positive-temperature).
For \(0<x\le1\) and \(0<y\le3\) the imaginary part of
\(G(-z)\) is negative by the decreasing sine transform,
whereas that of \(e^{z^2/2}\) is positive. For \(y\ge3\),
rotate the endpoint contour for \(G(-z)\) to angle
\(-\pi/4\). This gives
\[
                |G(-z)-1/z|<3/|z|^3 .
]
The endpoint modulus is then at least \(2/[3(y+1)]\)
when \(x\le1\), whereas the Gaussian exponential has
modulus at most \(3e^{(1-y^2)/2}\), strictly smaller for
\(y\ge3\). Both cases exclude a zero at \(x\le1\).

Two distinct Gaussian zeros of modulus at most \(41\)
are at complex distance at least \(1/200\), by the same
integrating-factor formula centered at one zero.
On the rectangle \(1\le x\le41\), \(0\le y\le41\), write
\(U=\log|G(-z)|-\Re(z^2)/2\).
The Gaussian Stieltjes formula and Cauchy–Schwarz give
\(|G(-z)|^2\le \Re G(-z)/x\).
Thus \(U_x=-\Re(1/G(-z))\le-1\).
Moreover \(|U_y|\le|1/G(-z)|\le100\):
for \(|z|\ge3\) use the rotated endpoint estimate, and
for \(|z|<3\) use the real part of the Stieltjes formula
restricted to \(-1\le v\le1\), giving
\(\Re G(-z)\ge1/100\).
All these zeros have the same value of \(U\).
Hence two distinct such zeros satisfy
\(|x_1-x_2|\le100|y_1-y_2|\), and their height difference
is at least \(1/20200\).

Distinct actual roots cannot correspond to the same
Gaussian zero in the \(10^{-6}\) neighborhood.
On that disk the local integrating-factor bound gives
\(G'=1+O(10^{-4})\), and (3) gives
\(|G_s'-G'|<3\cdot10^{-7}\). The latter estimate uses
\(|e^{w^2/2}|<2\) on the disk, since the exponential
has modulus below one at the Gaussian zero.
The derivative stays within a disk of radius less
than one about one, so integrating on any segment
gives injectivity. For different corresponding Gaussian
zeros, the certified height gap exceeds the two
\(10^{-6}\) displacements. Thus the actual roots of
modulus \(|z|\le40\) have distinct heights.

Finally, an exterior root sharing a height with one
satisfying \(|z|\le20\) must itself have \(|z|\le40\)
by the exterior cone. A pair of exterior roots is
already excluded by the horizontal derivative there.
This proves distinct heights for every \(T\ge10^{11}\).
\(\square\)
