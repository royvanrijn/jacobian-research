# Zero separation at large negative temperature for the half-line cubic

For real \(T\) and complex \(H\), write
\[
 F(T,H)=\int_0^\infty e^{-y^3/3-Ty^2/2+Hy}\,dy.
\]
At a simple zero let \(H_T\) denote the derivative of the
local zero branch in \(T\).

**Theorem.** For all sufficiently large negative \(T\), the
zeros of \(F(T,\cdot)\) are simple, distinct zeros have
distinct imaginary heights, and every zero satisfies
\(2H-T H_T\ne0\). Together with the
[local finiteness and the positive-temperature theorem](MULTIVARIATE_HALF_LINE_CUBIC_ZERO_SEPARATION.md),
this confines the union of the three exceptional temperature
sets to a finite subset of \(\mathbb R\). This set is unevaluated.

Put \(T=-A\) with \(A\to\infty\). There are two relevant
rescalings. For bounded \(H/A^2\), put
\[
              y=A z,\qquad b=H/A^2,\qquad M=A^3.
]
The exact scaled model, up to the factor \(A\), is the
integral of \(e^{M P_{a,b}(z)}\) over \(z\ge0\), with
\[
             P_{a,b}(z)=bz-a z^2/2-z^3/3,qquad a=-1 .
 \tag{1}
\]
For \(|H|/A^2\to\infty\), put \(R=\sqrt{|H|}\) and
\[
              y=R z,qquad b=H/R^2,qquad
              M=R^3,\qquad a=-A/R\to0 .
]
This gives (1) with \(|b|=1\); the value of \(R\) can be
held fixed when differentiating in \(H,T\) nearby.
The real axis itself has no zeros by positivity, and
conjugation permits the upper half-plane throughout.

## Endpoint and saddle packets

For clarity we give the features of the
[cubic contours](MULTIVARIATE_PRIVATE_ENDPOINT_UNIFORM_THERMAL_TAILS.md#the-contours-and-the-limiting-equal-modulus-fields)
needed here. At a nonreal limiting \(b\) in the upper
half-plane and \(a\in\{-1,0\}\), the right saddle is
\(z_s=x+iy=(-a+\sqrt{a^2+4b})/2\), with \(x,y>0\)
and \(2x+a>0\). Its value is
\(D=P_{a,b}(z_s)=a z_s^2/2+2z_s^3/3\).
Shift along the finite vertical segment from zero to
\(iy\), then along the rightward horizontal.
On the vertical segment the real phase is
\(-(Im b)v+a v^2/2\), strictly negative for \(0<v\le y\).
The real derivative along the horizontal at \(u+iy\) is
\((x-u)(x+u+a)\). Hence its only possible maxima are
the starting point and the saddle. The starting
point is strictly below zero and the saddle has
negative real curvature \(-2x-a\).
Thus a zero limit requires \(\Re D=0\), with precisely
an endpoint and a simple saddle active there.
The shifted path is an exact deformation by cubic
decay at positive real infinity. These estimates
persist for nearby actual parameters.

For \(a=-1\) on the real \(b\)-axis the sole equality
is \(b=-3/16\), where the real saddle is \(x=3/4\).
Indeed any strictly positive local maximum is at
\(x>1/2\) with value \(x^2(2x/3-1/2)\).
At equality split off a short real endpoint interval,
leaving the finite starting connector of the saddle
contour in a strict loss region. When \(a=0\) and
\(|b|=1\), the negative real ray has only an endpoint
maximum, and the positive real ray has a strictly
dominant saddle. In the nonreal upper half-plane
\(a=0\) gives \(D=2 b^{3/2}/3\), so the only possible
zero angle is \(\arg b=\pi/3\).
For the thermal model \(a=-1\), each imaginary height
has at most one limiting equality. This follows
from \(D_b=z_s\) and the fact that the right-saddle
branch is analytic along every strictly positive
height, with \(\Re z_s>0\). At height zero the
comparison of the real maxima gives the stated
unique equality.

Fix a limiting equality, including the real thermal
one. On a small parameter disk use the two exact
integrals \(E,S\) on the split contour. Use a fixed
split point from the limiting contour, so the
integrals are analytic in the parameters and also
in \(M\) in a fixed narrow sector of the positive
ray. For estimates one may move the second integral
through the nearby actual saddle, with short
connectors below its strict maximum. By the linear
loss at the endpoint and the Gaussian maximum at
the saddle,
\[
 \begin{split}
 E(M,a,b)&=M^{-1} e(a,b)(1+O(M^{-1})),\
 S(M,a,b)&=M^{-1/2}s(a,b)e^{M D(a,b)}(1+O(M^{-1})).
 \end{split}                                                  \tag{2}
\]
Here \(e=-1/b\) and \(s\) is the nonzero analytic
Gaussian coefficient. The bounds are uniform on
smaller fixed parameter disks and a fixed order
sector. In particular the endpoint coefficient
comes from its linear exponential; the quadratic
local error integrates to a relative \(O(M^{-1})\).
At the saddle expand to cubic order on a small
symmetric interval. The first cubic term is odd
and integrates to zero, and the quartic and
squared-cubic terms have relative size \(O(M^{-1})\).
The cubic tail at large positive real part controls
both estimates on the complex order sector, and
all bounded complementary pieces have a strict loss.
Thus the analytic logarithm of the packet ratio is
\[
 W(M,a,b)=M D(a,b)+\tfrac12\log M+L(a,b)+R(M,a,b),qquad
 L=\log(s/e),qquad R=O(M^{-1}).                               \tag{3}
\]
At a zero, \(e^W=-1\). The choice of the additive
imaginary constant in \(L\) has no effect on the
derivatives. Cauchy's estimate gives
\(R_M=O(M^{-2})\) and \(R_b,R_a=O(M^{-1})\) on smaller
disks. Stationarity gives
\[
                          D_b=z_s,qquad D_a=-z_s^2/2 .
 \tag{4}
\]
In particular \(\Re W_b\ge c M\) for real large \(M\).
This implies simplicity, as well as uniqueness of
the zero at a given imaginary height in the
parameter disk. Indeed a zero requires \(\Re W=0\),
and \(\Re W\) strictly increases on a horizontal
segment in that disk.

## Bounded thermal fields

At a zero near a regular thermal equality,
\(\Re W=0\) and (3) imply
\[
 \Re D=-\frac{\log M}{2M}+O(M^{-1}),\qquad
 \Re W_M=-\frac{\log M}{2M}+O(M^{-1}).
]
Thus \(b_M=-W_M/W_b\ne0\) on its simple branch,
with \(a=-1\) fixed. In physical variables
\(H=A^2 b(A^3)\) and \(T=-A\), so
\[
                     2H-T H_T=-3A^5 b_M\ne0 .
]
If two bounded thermal fields at the same imaginary
height were zeros along a sequence of temperatures,
pass to their limiting fields. Both would have to be
at an equality at the same limiting height, which
is unique. They would then lie in the same parameter
disk in (3), contradicting the strict horizontal
monotonicity there.

## Escaping fields and mixed pairs

At an escaping zero its field angle tends to
\(\pi/3\), by the localization for \(a=0\).
Use (3) and (4) with the real \(R\) appropriate to
the root held fixed when differentiating. In these
coordinates
\[
 (2H\partial_H+T\partial_T)W
     = M(2b z_s-a z_s^2/2)+O(1)
     = M(2z_s^3+3a z_s^2/2)+O(1),
]
since the saddle satisfies \(b=z_s^2+a z_s\).
The leading coefficient tends to \(2z_s^3\ne0\).
Together with \(W_H=R^{-2}W_b\ne0\), implicit
differentiation proves \(2H-T H_T\ne0\).
For two escaping roots at the same height their
arguments both tend to \(\pi/3\), so their moduli
are asymptotically equal. Use the same \(R\) for
both, taken from the first root. Both scaled fields
then tend to the same \(e^{i\pi/3}\) and are
excluded by the strict horizontal monotonicity of
(3). A bounded thermal field cannot share a height
with an escaping one, since the latter has
imaginary part comparable to its modulus and
that modulus divided by \(A^2\) diverges.

If any of the three obstructions persisted at
arbitrarily large negative temperatures, a sequence
of those temperatures would have a subsequence in
one of the two field regimes or the mixed case
just treated. All are excluded, which proves the
theorem. \(\square\)
