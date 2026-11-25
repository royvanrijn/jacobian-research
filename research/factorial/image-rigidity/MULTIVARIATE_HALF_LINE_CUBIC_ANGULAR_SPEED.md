# Positive angular speed for all half-line cubic zeros

For real \(T\), let
\[
 F(T,H)=\int_0^\infty \exp(-u^3/3-Tu^2/2+Hu)\,du.
\]
Let \(H_T\) denote the temperature derivative on a simple zero branch.

**Theorem.** For every real \(T\) and every upper-half-plane zero \(H\) of
\(F(T,\cdot)\),
\[
                     \Im(2H-T H_T)>0.
\]
Lower-half-plane zeros have the opposite sign. The real axis has no zeros.

The [Riccati sign theorem](MULTIVARIATE_HALF_LINE_CUBIC_RICCATI_POSITIVITY.md)
proves simplicity, \(\Im H_T>0\) in the upper half-plane, and the asserted
speed sign for \(T\le0\).
The [extended temperature certificate](MULTIVARIATE_HALF_LINE_CUBIC_EXTENDED_TEMPERATURES.md)
gives the sign on \(0\le T\le2\).
The following exterior comparison and compact interval replay cover \(T\ge2\).

## Exterior fields

Let \(a=T/2\), \(R=|H|\), \(Z=a^2+H\), \(P=\sqrt Z-a\), and
\(\omega=e^{2\pi i/3}\). An upper zero has \(\Re H>0\), by the
decreasing-density sine-transform comparison in the
[compact speed reduction](MULTIVARIATE_HALF_LINE_CUBIC_SPEED_COMPACT_FIELD.md#exterior-speed-exclusion).
Use its Airy decomposition \(F=E+S\) and the local logarithm
\(\Xi=\log(S/E)\), with
\[
 E=\omega F(T\omega^2,\omega H),\qquad
 S=2\pi i\omega^2e^{-aH-2a^3/3}\operatorname{Ai}(\omega^2Z).
\]
The endpoint ray \(3\pi/4\) and the Airy derivative estimate from that
proof apply for \(T\ge2\), \(R\ge4\sqrt T\).
In particular \(|E+1/H|<1/(4R)\), \(|E_H/E|\le8/(3R)\),
\(|E_T/E|\le 8\sqrt2/(3R^2)\), and the Airy logarithmic derivative
differs from \(\sqrt Z\) in the \(Z\) coordinate by less than \(9/(4|Z|)\).
Both terms are nonzero. The derivatives of their log ratio satisfy
\[
 |\Xi_H-P|<\frac5R,\qquad
 \left|\Xi_T+\frac{P^2}{2}\right|
   \le\frac{9a}{4|Z|}+\frac{8\sqrt2}{3R^2}.                 \tag{1}
\]

Suppose now \(R\ge6\sqrt T\). We first exclude small angles at a zero.
Write \(V=R/(a+\sqrt R)^{1/2}\). The exact saddle exponent is
\[
 \mathcal S=\frac23Z^{3/2}-aH-\frac23a^3
       =H^2\int_0^1\frac{1-u}{2\sqrt{a^2+uH}}\,du .
\]
For \(0\le\arg H\le\pi/6\), the integrand after multiplication by
\(H^2\) has argument between \(0\) and \(\pi/3\).
Since \(|\sqrt{a^2+uH}|\le a+\sqrt R\), this gives
\(\Re\mathcal S\ge V^2/8\).
The Airy value's relative error is bounded by
\((25/24)|Z|^{-3/2}<1/4\), as in the
[explicit Airy bounds](EXPLICIT_AIRY_CONTOUR_ENVELOPE.md#3-explicit-primary-airy-error-bounds-when-xige4).
Together with \(|E|\le5/(4R)\) and
\(|Z|^{1/4}\le(a+\sqrt R)^{1/2}\), it yields
\[
          |S/E|\ge (V/2)\exp(\Re\mathcal S).
\]
The function \(R^2/(2a+\sqrt R)\) is increasing in \(R>0\)
for fixed \(a\ge0\). At \(R\ge6\sqrt T\), \(T\ge2\),
\[
 V^2\ge\frac{R^2}{2a+\sqrt R}
      \ge\frac{36}{1+\sqrt6/2^{3/4}}
      >\frac{6048}{413}>\frac{73}{5},                      \tag{2}
\]
where \(\sqrt6<49/20\) and \(2^{3/4}>42/25\).
Thus \(|S/E|>1\) at angles at most \(\pi/6\).
At a zero \(|S/E|=1\), so \(\arg H>\pi/6\).

The identity \(P=H/(\sqrt Z+a)\) gives
\[
 |P|\ge\frac{R}{2a+\sqrt R},\qquad a|P|\le R/2,
\]
since \(\Re\sqrt Z\ge a\).
Consequently (1) and (2) show that \(\Xi_H\ne0\).
Differentiating \(S/E=-1\) on the local zero branch gives
\(H_T=-\Xi_T/\Xi_H\), hence
\[
 \left|(2H-T H_T)-(2H-aP)\right|
 \le
 \frac{T|\Xi_T+P^2/2|+a|P||\Xi_H-P|}{|\Xi_H|}
 <\frac{(36/5)R}{R^2/(2a+\sqrt R)-5}<\frac{3R}{4}.
\]
Indeed \(|Z|\ge a^2\), \(R^2\ge36T\), and the numerator is
at most \(9/2+8\sqrt2/108+5/2<36/5\).
For \(Y=\Im H>0\) we have \(\Im P=Y/(2\Re\sqrt Z)\).
This gives \(\Im(2H-aP)\ge3Y/2>3R/4\), proving the required sign.

## Compact interval replay

For \(T\ge2\) put \(s=T^{-3/2}\), \(z=H/\sqrt T=x+iy\), and
\[
   G(s,z)=\int_0^\infty e^{-u^2/2-su^3/3+zu}\,du.
\]
By the [root-free strip](MULTIVARIATE_HALF_LINE_CUBIC_SPEED_COMPACT_FIELD.md#a-strip-without-roots)
and the exterior proof, it is enough to consider
\[
 0\le s\le3/8,\qquad 0\le x\le6,\qquad 7/5\le y\le6,\qquad
 x^2+y^2\le36.                                             \tag{3}
\]
The larger closed box includes the limiting value \(s=0\).
At a zero, integration by parts gives \(G_z+sG_{zz}=1\).
Thus
\[
 \Im\bigl((4zG_z-G_{zz})\overline{G_z}\bigr)
      =4y|G_z|^2-\Im G_{zz}.
\]
For \(s>0\), the [scaled speed identity](MULTIVARIATE_HALF_LINE_CUBIC_SPEED_COMPACT_FIELD.md#a-strip-without-roots)
gives
\[
 \Im(2H-T H_T)
      =\frac{\sqrt T}{2|G_z|^2}\bigl(4y|G_z|^2-\Im G_{zz}\bigr).
                                                               \tag{4}
\]
Simplicity ensures that the denominator is positive.

The [interval checker](verify_half_line_cubic_angular_speed.py) uses the
same exact dyadic contour choice, complex integration on two finite
segments, exponential tail, and Gaussian moment bounds as the
[compact speed certificate](MULTIVARIATE_HALF_LINE_CUBIC_POSITIVE_SPEED_CERTIFICATE.md).
For a box with radii \(h_s,h_x,h_y\) and \(h_z=(h_x^2+h_y^2)^{1/2}\),
the bounds \(B_k\) there give
\[
 |G_z-G_{z,0}|\le h_z B_2+h_sB_4/3=\delta_1,\qquad
 |G_{zz}-G_{zz,0}|\le h_z B_3+h_sB_5/3=\delta_2,
\]
where the subscript zero means the box center. These follow by
integrating the directional derivative along the segment in parameter
space from the center to the evaluation point. The same second-order
bound for \(G\) as in the compact speed certificate excludes a root
whenever a center coordinate is strictly separated from zero.
The new sign test is the strict enclosure
\[
       4y_{\min}\max(0,|G_{z,0}|-\delta_1)^2
                > \Im G_{zz,0}+\delta_2 .
\]
The checker rounds the left side down and the right side up.
Thus even when the root exclusion fails, the sign test certifies
(4) at every possible root in the box. A box is skipped only when
its lower \(x,y\) corner has squared norm greater than \(36\).
Exact bisections and those two strict certificates cover (3).
The 128-bit replay checked 46,040 non-skipped boxes, with
829 outside-disk leaves, 22,096 root exclusions,
and 510 strict speed certifications.

Replay with `python-flint==0.9.0` from the repository root.

~~~sh
python research/factorial/image-rigidity/verify_half_line_cubic_angular_speed.py
~~~

## Effective angular exclusion in bounded ordinary private windows

Consider the consecutive ordinary factorial moment polynomials of a
two-monomial profile with exactly one private unequal row, oriented
\((0,1)\), and suppose \(3/4\le\sigma\le5/4\). Use the symbols
\(\sigma,\ell,T_j,H_j,\epsilon_j,\xi_j\) of the
[effective bounded-window proof](MULTIVARIATE_EFFECTIVE_BOUNDED_PRIVATE_CUBIC.md).
For every pair of rational bounds \(B_T,B_H\ge1\) there is a computable
terminal cutoff excluding equal arguments for zeros at orders \(n\)
and \(N=n+1\) whenever both have \(|T_j|\le B_T\) and \(|H_j|\le B_H\).
Positive rescalings of either polynomial preserve the conclusion.

For a proof, let \(\lambda=\log(1+1/n)\), \(T=T_n\), \(H=H_n\), and
\(\epsilon=\epsilon_n\).
Equal arguments imply that the two \(\xi_j\) have the same imaginary
part. Indeed both are on the positive-endpoint chart for large \(n\),
and all changes in their logarithm normalizations are real.
Consequently
\[
 T_N=e^{\lambda/3}T,\qquad
 H_N=e^{2\lambda/3}H+v,\qquad
 \epsilon_N=e^{-\lambda/3}\epsilon,\qquad v\in\mathbb R.
\]
The bounded fields give \(|v|\le2B_H+1\) for large \(n\).
Let \(D=B_T+2B_H\) and choose a rational \(U\ge1\) bounding the remainder
constant and all derivatives of order at most two of \(F,C_1,C_2,C_3\)
on \(|T|\le B_T+1,\ |H|\le4B_H+4\), as in the real-lattice expansion
\[
 \mathcal Q_j=F(T_j,H_j)+\sum_{p=1}^3\epsilon_j^p C_p(T_j,H_j)
                          +\mathcal R_j,\qquad
 |\mathcal R_j|\le U\epsilon_j^4 .
\]
Here \(\mathcal Q_j\) is a nonvanishing normalization of the exact
moment polynomial at its own field and \(C_p\) are common to the
two orders. In particular \(|F(T,H)|\le4U\epsilon\) at a root.

Define
\[
 B(T,H,v)=\int_0^1 F_H(T,H+t v)\,dt,\qquad
 A(T,H)=T F_T(T,H)+2H F_H(T,H).
\]
On the compact parameter rectangle \(|T|\le B_T\),
\(|\Re H|,|\Im H|\le B_H\), \(|v|\le2B_H+1\), we can
compute positive rational lower bounds \(\mu_1,\mu_0\) for
\[
          |F|+|B|,\qquad |F|+|\Im(A\overline{F_H})|,
                                                               \tag{5}
\]
respectively. At a model zero, \(B\ne0\) by the
[Riccati sign theorem](MULTIVARIATE_HALF_LINE_CUBIC_RICCATI_POSITIVITY.md)
and its same-height exclusion, using simplicity when \(v=0\).
The new speed theorem gives the second strict positivity because
\(A=F_H(2H-T H_T)\) there, with the sign changing under conjugation.
There are no real model zeros. Compactness then gives positive
infima. Rational cell quadrature and derivative enclosures compute
positive lower bounds as in the existing bounded-window proof.
For the second expression this uses only the displayed derivatives
of \(F\) and ordinary complex multiplication.

For \(n\ge1\), the parameter changes satisfy
\[
 |T_N-T|+|H_N-H-v|\le D\lambda,\qquad
 |\epsilon_N^p-\epsilon^p|\le\lambda\epsilon^p,\qquad
 \lambda\ge1/(2n).
\]
At the proposed zeros the linked difference of the exact lattice
expansions gives
\[
 |vB(T,H,v)|
 \le U D\lambda+3U\epsilon|v|
      +3U(D+1)\epsilon\lambda+2U\epsilon^4
 \le (4D+7)U\lambda+3U\epsilon|v|.
\]
The first inequality follows by replacing the leading difference
by \(vB\) and bounding the change of each correction by
\(U\epsilon^p(|v|+(D+1)\lambda)\).
If \(16U\epsilon\le\mu_1\), inequality (5) gives
\(|B|-3U\epsilon\ge\mu_1/2\), and hence
\[
                 |v|\le K\lambda,\qquad K=(8D+14)U/\mu_1 .
\]

Now expand the linked difference at \((T,H)\) to first order.
Its leading model term is \((\lambda/3)A+vF_H\)
with error at most \(U(D+(D+K)^2)\lambda^2\).
Indeed the field increment differs from
\((\lambda T/3,\,2\lambda H/3+v)\) by at most
\(D\lambda^2\) in the sum norm, and the sum norm of the
actual increment is at most \((D+K)\lambda\).
After division by \(\lambda/3\), the correction sum has
modulus at most \(9U(K+D+1)\epsilon\), and the two
remainders cost at most \(12U\epsilon\).
With \(M=3(D+(D+K)^2)+9(K+D+1)+12\), this yields
\[
                |A+3(v/\lambda)F_H|\le UM\epsilon .
\]
Multiplication by \(\overline{F_H}\) and taking the imaginary part
eliminate the \(v\) term because \(v/\lambda\) is real.
Thus \(|\Im(A\overline{F_H})|\le U^2 M\epsilon\), contrary to
(5) when \((4U+U^2M)\epsilon<\mu_0\).
The tests on \(\epsilon\) and the original lattice cutoff are
computable and uniform over the stated profiles.
