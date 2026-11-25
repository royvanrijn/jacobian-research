# A compact interval certificate for positive-temperature cubic speed

Write
\[
 G(s,z)=\int_0^\infty e^{P(u;s,z)}\,du,\qquad
 P(u;s,z)=-u^2/2-su^3/3+zu,\qquad
 D(s,z)=4zG_z(s,z)-G_{zz}(s,z).
\]
The [compact field reduction](MULTIVARIATE_HALF_LINE_CUBIC_SPEED_COMPACT_FIELD.md)
shows that a speed zero of a simple upper zero at temperature \(T\ge2\)
would give \(G(s,z)=D(s,z)=0\), where \(s=T^{-3/2}\), \(z=H/\sqrt T\),
\(\Re z>0\), \(\Im z>\sqrt2\), and \(|z|<4\).

**Proposition.** The functions \(G,D\) have no common zero on the closed
domain
\[
 0\le s\le 3/8,\quad
 0\le x=\Re z\le4,\quad 7/5\le y=\Im z\le4,\quad x^2+y^2\le16.
 \tag{1}
\]
In particular \(2H-T H_T\ne0\) at every simple upper zero for \(T\ge2\).
Complex conjugation gives the corresponding lower-half-plane assertion.

## Center values and bound over a cell

The [verifier](verify_half_line_cubic_positive_speed.py) covers (1) by
exact rational closed boxes. On each box let \((s_0,x_0,y_0)\) be its
center and \(z_0=x_0+iy_0\). The verifier chooses one exact dyadic
\(b>0\) for the whole box, using an Arb midpoint of
\[
        \Im\frac{2z_0}{1+\sqrt{1+4s_0z_0}} .
\]
Only the location of the contour is affected by the midpoint rounding.
The integrals for \(G\) and its first three field derivatives are taken
along \(0\to ib\to\infty+ib\). For each fixed \(b\), the deformation
from the nonnegative real axis holds by Cauchy's theorem. On the right
vertical connector with real part tending to infinity, the term
\(-u^2/2\) controls the integrand uniformly for \(s\ge0\) in a box.
The cubic makes an additional negative contribution for sufficiently
large real part.

Arb complex quadrature encloses each integral of \(u^k e^{P(u;s_0,z_0)}\)
on the first two finite segments \(0\to ib\to L+ib\), where \(L=16\).
The integrand is entire. On the horizontal tail the real exponent
has derivative \(x_0-t-s_0(t^2-b^2)\), with \(t=\Re u\).
The verifier checks \(d=L-x_0+s_0(L^2-b^2)>0\).
Since the derivative decreases for \(t\ge L\), the omitted tail has
modulus at most
\[
 e^{\Re P(L+ib;s_0,z_0)}
 \sum_{m=0}^k {k\choose m}(L+b)^{k-m}\frac{m!}{d^{m+1}}.
 \tag{2}
\]
Here \(|t+ib|\le t+b\) and the terms are elementary exponential
moments in \(t-L\). The code adds bound (2) to the radius of each
coordinate of the quadrature result.

For each \(0\le k\le6\), let \(B_k\) bound the integral of
\(|u|^k e^{\Re P(u;s,z)}\,|du|\) on the displaced contour for every
point in the box. One can obtain such bounds as follows. On the vertical
segment the real exponent at \(u=iv\) is \(v^2/2-yv\), \(0\le v\le b\).
If \(d_v=y_{\min}-b/2>0\), the vertical integral is at most
\[
       \min\left(\frac{b^{k+1}}{k+1},\frac{k!}{d_v^{k+1}}\right).
\]
Otherwise convexity bounds the exponent by the larger of its endpoint
values at \(y=y_{\min}\). The absolute value of \(b^2/2-y_{\min}b\)
gives a valid upper bound for this maximum.

On the horizontal ray \(u=t+ib\), \(t\ge0\), the real exponent is
\[
 -t^2/2+(x+sb^2)t-st^3/3+b^2/2-yb.
\]
Majorize it by replacing \(x,s,y\) by the appropriate box endpoints.
Writing \(A=x_{\max}+s_{\max}b^2\) and
\(a=2A/(1+\sqrt{1+4s_{\min} A})\), the resulting majorant is maximized
at \(a\) and has second derivative at most \(-1\) on \(t\ge0\).
It is therefore bounded by \(C-(t-a)^2/2\), where
\[
 C=-a^2/2+Aa-s_{\min}a^3/3+b^2/2-y_{\min}b.
\]
Using \(|t+ib|\le\sqrt{a^2+b^2}+|t-a|\) and extending the integrals in
\(t-a\) to the full real line gives the horizontal bound
\[
  e^C\sum_{m=0}^k{k\choose m}(a^2+b^2)^{(k-m)/2}
       2^{(m+1)/2}\Gamma((m+1)/2).
 \tag{3}
\]
All moment bounds are combined using outward-rounded Arb arithmetic.

## Excluding a common zero on a box

Let \(h_s,h_x,h_y\) be the coordinate radii and
\(h_z=(h_x^2+h_y^2)^{1/2}\). The affine part of the change in \(G\)
from the center is \((z-z_0)G_z-(s-s_0)G_{zzz}/3\), with derivatives
evaluated at the center. The integral remainder along the straight
segment to any point in the box is bounded in modulus by
\[
            \frac{h_z^2 B_2}{2}
           +\frac{h_z h_s B_4}{3}
           +\frac{h_s^2 B_6}{18}.
 \tag{4}
\]
The second directional derivative introduces the factor
\(((z-z_0)u-(s-s_0)u^3/3)^2\) inside the integral.
The center value of \(G\) excludes zero when the lower bound for
the absolute value of either coordinate exceeds the sum of the
corresponding affine coordinate bound and (4).

If that test fails, bounds on the changes in \(G_z,G_{zz}\) are
\[
    \delta_1=h_zB_2+h_sB_4/3,\qquad
    \delta_2=h_zB_3+h_sB_5/3 .
\]
Writing \(p=G_z(s_0,z_0)\), these give
\[
       |D(s,z)-D(s_0,z_0)|
           \le 4h_z|p|+4(|z_0|+h_z)\delta_1+\delta_2.        \tag{5}
\]
The box is certified if the lower bound for \(|D(s_0,z_0)|\)
strictly exceeds the upper bound in (5).

Boxes whose lower \(x,y\) corners satisfy \(x_{\min}^2+y_{\min}^2>16\)
are outside (1). Each other box is tested as above or bisected exactly
across a widest coordinate after multiplying the \(s\)-width by eight.
The verifier requires termination by depth 33.
The replay at 128-bit precision checked 6,526 non-skipped boxes and
finished with 225 outside-disk leaves, 2,988 leaves certified with
\(G\), and 163 with \(D\). The split boxes cover their parent boxes,
so these leaf tests prove the proposition.

From the repository root the standalone replay command, with
`python-flint` installed, is

```sh
python research/factorial/image-rigidity/verify_half_line_cubic_positive_speed.py
```
