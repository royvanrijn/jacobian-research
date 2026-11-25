# A compact field region for positive-temperature speed zeros

Let
\[
 F(T,H)=\int_0^\infty e^{-u^3/3-Tu^2/2+Hu}\,du .
\]
At a simple zero, \(H_T\) denotes the local derivative of the
zero branch with respect to \(T\).

**Theorem.** For \(T\ge2\) a simple upper zero with
\(2H-T H_T=0\) must satisfy
\[
            \Re H>0,\qquad
            \Im H>\sqrt{2T},\qquad |H|<4\sqrt T.
                                                               \tag{1}
\]
In scaled coordinates \(s=T^{-3/2}\), \(z=H/\sqrt T\), the
remaining condition is that
\[
 G_s(z)=0,\qquad 4zG_s'(z)-G_s''(z)=0,\qquad
 G_s(z)=\int_0^\infty e^{-u^2/2-su^3/3+zu}\,du                 \tag{2}
\]
with \(0<s\le2^{-3/2}\) and \(z\) in the region given by (1).
These equations are nonsingular at the Gaussian limit \(s=0\).

## Exterior speed exclusion

For \(\Re H\le0\) and \(T\ge0\) the function
\(u\mapsto e^{-u^3/3-Tu^2/2+(\Re H)u}\) is strictly decreasing
on \(u>0\). Its sine transform is positive at every positive
frequency. To see this, pair consecutive intervals of length
\(\pi/\Im H\) in the sine integral. Thus \(\Re H>0\) at an
upper zero.

Write \(a=T/2\), \(Z=H+a^2\), \(R=|H|\), and
\(\omega=e^{2\pi i/3}\). Use the same
[exact Airy identity](MULTIVARIATE_HALF_LINE_CUBIC_ZERO_SEPARATION.md#an-exact-airy-decomposition)
as in the [positive tail comparison](MULTIVARIATE_HALF_LINE_CUBIC_EXPLICIT_POSITIVE_TAIL.md).
It writes \(F=E+S\) with
\[
 E=\omega F(T\omega^2,\omega H),\qquad
 S=2\pi i\omega^2 e^{-aH-2a^3/3}\operatorname{Ai}(\omega^2 Z).
\]
On the endpoint ray of angle \(3\pi/4\), the quadratic
exponent is purely imaginary and the cubic exponent has
negative real part. The ray is obtained from the original
angle \(2\pi/3\) by deformation in the cubic decay sector.
Throughout the closed upper right field quadrant the
real linear exponent on this ray is at most \(-Ru/\sqrt2\),
with \(u\) the positive ray coordinate. Comparing to the
linear integral gives
\[
 \begin{aligned}
 |E+1/H| &\le R^{-1}(8/R^3+2\sqrt2\,T/R^2),\\
 |E_H|&\le 2/R^2,\qquad |E_T|\le 2\sqrt2/R^3 .
 \end{aligned}
\]
For example the value bound follows from
\(\lvert e^v-1\rvert\le |v|e^{\max(0,\Re v)}\) applied
to the quadratic and cubic exponents together. The
three estimates then integrate elementary exponential
moments. If \(R\ge4\sqrt T\), \(T\ge2\), the relative
value error is at most \(5\sqrt2/32<1/4\). Hence
\[
             |E_H/E|\le 8/(3R),\qquad
             |E_T/E|\le 8\sqrt2/(3R^2).                      \tag{3}
\]

In the same range \(|Z|\ge4\sqrt2\) and
\(\lvert\arg(\omega^2 Z)\rvert\le2\pi/3\).
The [explicit Ai value and derivative bounds](EXPLICIT_AIRY_CONTOUR_ENVELOPE.md#3-explicit-primary-airy-error-bounds-when-xige4)
give relative errors \((25/24)|Z|^{-3/2}\) and
\((49/48)|Z|^{-3/2}\) respectively in their leading terms.
In particular \(S\ne0\), and
\[
  \left|\omega^2\frac{\operatorname{Ai}'(\omega^2 Z)}
                 {\operatorname{Ai}(\omega^2 Z)}-\sqrt Z\right|
   \le \frac{99}{48|Z|(1-(25/24)|Z|^{-3/2})}
   <\frac{9}{4|Z|}.
\]
The last inequality uses \((4\sqrt2)^{3/2}>25/2\).

Set \(\Xi=\log(S/E)\) locally and \(P=\sqrt Z-a\).
The leading terms of \(\Xi_H,\Xi_T\) are \(P,-P^2/2\).
By (3) and the Airy bound, the Euler numerator
\(2H\Xi_H+T\Xi_T\) differs from \(P(2H-aP)\) by at most
\[
 \frac{9(2R+Ta)}{4|Z|}+\frac{16}{3}
                +\frac{8\sqrt2\,T}{3R^2}
 \le \frac{9\sqrt2}{2}+\frac{16}{3}+\frac{\sqrt2}{6}<12 .
                                                               \tag{4}
\]
Here \(Ta=2a^2\) and \(|Z|^2\ge R^2+a^4\); the last
comparison follows from \(\sqrt2<10/7\).
On the other hand \(P=H/(\sqrt Z+a)\) and
\(\Re\sqrt Z\ge a\), so \(a|P|\le R/2\). It follows that
\[
 |P(2H-aP)|\ge\frac{3R^2}{2(a+\sqrt{a^2+R})}
       \ge \frac{24}{1/2+\sqrt{1/4+\sqrt2}}>13 .                \tag{5}
\]
For the second inequality the preceding ratio increases
with \(R\), since \(R^2/(a+\sqrt{a^2+R})=
R(\sqrt{a^2+R}-a)\). At \(R=4\sqrt T\) the ratio on
the right of (5) is minimized at \(T=2\).
For the last comparison use
\(\sqrt{1/4+\sqrt2}<4/3\).
At a simple zero \(H_T=-\Xi_T/\Xi_H\).
The strict inequalities (4) and (5) exclude zero order
speed in the entire exterior \(R\ge4\sqrt T\).

## A strip without roots

Let \(z=x+iy\) and regard \(G_s(x+iy)/G_s(x)\) as the
characteristic function of a positive random variable \(U\)
with density proportional to \(e^{-V(u)}\) on \(u\ge0\),
where \(V(u)=u^2/2+su^3/3-xu\) and \(s\ge0\).
Integration by parts gives
\[
  \mathbb E[U V'(U)]=1,\qquad
  \mathbb E[V'(U)]=1/G_s(x)>0.
\]
If \(U'\) is an independent copy, then
\((U-U')(V'(U)-V'(U'))\ge (U-U')^2\) because \(V''\ge1\).
Consequently
\(\operatorname{Var}(U)\le\operatorname{Cov}(U,V'(U))<1\).
The real part of the centered characteristic function is
at least \(1-y^2\operatorname{Var}(U)/2\), using
\(\cos v\ge1-v^2/2\). It is therefore positive for
\(0\le y\le\sqrt2\), proving the lower height in (1).

Finally \(sG_s''+G_s'-zG_s=1\) by integration by parts.
Differentiating in \(z\) and using
\(\partial_s G_s=-G_s'''/3\) shows that at a simple zero
\[
          2H-T H_T
            =\frac{\sqrt T}{2G_s'(z)}
                        (4zG_s'(z)-G_s''(z)).
\]
This identity gives (2) without dividing by \(s\).
The same equations at \(s=0\) cannot share a zero:
\(G_0'=1+zG_0\) gives \(G_0'=1\) and \(G_0''=z\)
there, so the second expression in (2) is \(3z\ne0\).
