# Priced slice repair and the variance required at a nonsimplex tie

Use the Gaussian cone profiles \(h_{n,\omega}\), their mean \(\bar h_n\),
the ratio marginal \(\mu_n\), and the radial lift form \(B_n\) from
(AS1). Suppose balanced three-cell cones \(F_n\) converge in Gaussian
label distance to the cylinder of a fixed planar global maximizer
\(A\) at correlation \(0<r<1\). Put \(V=S_r(A)\), write \(C(A)\)
for the radial loss coefficient in (CH2), and set
\[
\begin{split}
 W_n&=\int\sum_i\|T_{\sqrt r}(h_{n,\omega,i}-\bar h_{n,i})\|_2^2
                                                          \,d\sigma(\omega),\\
 Q_n&=\int |m_n(\omega)-m|^2\,d\sigma(\omega),\qquad
 m_n(\omega)=\int h_{n,\omega}\,d\mu_n .
\end{split}
\]
Here \(m\) is the common balanced mass vector and the \(L^2\) norms
use standard planar Gaussian measure. For every \(r<\theta<1\) there
are finite constants \(C_1,C_2\) such that
\[
 V-S_r(F_n)\ \ge\
 \frac{C(A)}n+(1-\theta)W_n
            -C_1\sqrt{W_nQ_n}-C_2Q_n+o(n^{-1}).              \tag{PR1}
\]
For spectator-even cones the result holds for every
\(r^2<\theta<1\). In particular a nonsimplex limit \(A\) and a
deficit \(o(n^{-1})\) imply \(\liminf nQ_n>0\).
[Angular slice comparison](GAUSSIAN_CONE_ANGULAR_SLICE_COMPARISON.md)

## Use a common translation price

For a planar simplex-valued profile \(h\), let \(h^v(x)=h(x-v)\).
Its Gaussian mass derivative \(M'_h(0)\) and stability derivative
\(S'_h(0)\), as functions of the translation vector \(v\), are
\[
 M'_h(0)v=-\int (x\cdot v)h(x)\,d\gamma_2(x),\qquad
 S'_h(0)v=-r^{-1/2}\int (x\cdot v)\sum_i(T_{\sqrt r}h_i)^2\,d\gamma_2 .
\]
The second formula follows by translating the Gaussian argument
after applying \(T_{\sqrt r}\). Both maps of the translation
parameter have uniformly bounded second derivatives on a fixed
small parameter ball, by differentiating the Gaussian density.
For any two bounded simplex-valued profiles, the differences of
their first derivatives at zero are bounded by a constant times
their distance in the norm
\(\|h\|_r^2=\sum_i\|T_{\sqrt r}h_i\|_2^2\). For masses use that
linear functions are eigenfunctions of \(T_{\sqrt r}\). For
stability use the boundedness of the smoothed profiles and
Cauchy–Schwarz with \(x\cdot v\).

Write \(\beta_n\) for the centered vector determined by
\(S'_{\bar h_n}(0)=2\beta_n\cdot M'_{\bar h_n}(0)\). This equation
has a unique solution for large \(n\). In fact \(\bar h_n\to A\)
in Gaussian \(L^1\) by the profile convergence argument for (AS4),
and the limiting mass derivative is invertible onto the
two-dimensional mass tangent space by the full centroid rank
of \(A\). Thus \(\beta_n\to\beta(A)\), by the first variation
identity at \(A\).
[Translation rank](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#one-active-dimension-is-impossible)
[Volume derivative and score prices](THREE_CELL_VOLUME_VALUE_DERIVATIVES.md#ambient-transport-and-uniform-mass-correction)

Given \(h\), put \(e=\|h-\bar h_n\|_r\) and
\(\delta_G=\int h\,d\gamma_2-m\). For \(e\) below a fixed
threshold and \(n\) large, the mass map of \(h\) has invertible
derivative at zero with uniformly bounded inverse, and
\(\delta_G\) is small. Uniform second derivatives give an exact
correcting translation \(v\) with \(|v|\le C|\delta_G|\).
Subtract twice the common price from its stability derivative.
The derivative at zero of the resulting priced objective is
bounded by \(Ce\), since the corresponding derivative for
\(\bar h_n\) vanishes. The corrected profile is feasible for \(V\),
which yields
\[
 S_r(h)\le V+2\beta_n\cdot\delta_G+
                             C(e|\delta_G|+|\delta_G|^2).    \tag{PR2}
\]
The same inequality holds for \(e\) above the threshold after
increasing \(C\). In that case direct mass transfers give
\(S_r(h)\le V+C|\delta_G|\), and \(e|\delta_G|\) absorbs the
transfer and the bounded price. The argument also applies to
fractional profiles. Their mass corrections can move fractions
between labels, and their balanced stability is at most the
ordinary optimum by positivity and rounding.

## Average the corrected bounds

Write \(\delta=\int h\,d\mu_n-m\). The law expansions (AX2) give
\(\delta_G=\delta-n^{-1}\mathfrak m(h)+O(n^{-3/2})\) and
\(B_n(h)=S_r(h)+n^{-1}\mathfrak s(h)+O(n^{-3/2})\), uniformly
over profiles. The difference \(|\delta-\delta_G|\) is \(O(n^{-1})\)
and \(|\delta|\le e+O(n^{-1})\) since \(\bar h_n\) is balanced
under \(\mu_n\). Applying (PR2) thus gives
\[
 B_n(h)\le V+2\beta_n\cdot\delta
       +n^{-1}\bigl(\mathfrak s(h)-2\beta_n\cdot\mathfrak m(h)\bigr)
       +C(e|\delta|+|\delta|^2)+O(e/n+n^{-3/2}).
\]
The deviations of the slices satisfy \(\int\delta\,d\sigma=0\)
exactly, so their linear prices cancel. The mean Gaussian \(L^1\)
distance of the slices to \(A\) tends to zero by (AS4)'s
profile convergence argument. Integrability of the coefficient
densities in (AX2) and boundedness of the profiles imply that
the mean of the coefficient in parentheses tends to
\(\mathfrak s(A)-2\beta(A)\cdot\mathfrak m(A)=-C(A)\).
Therefore
\[
 D_n:=\int B_n(h_{n,\omega})\,d\sigma(\omega)
 \ \le\ V-C(A)/n+C_1\sqrt{W_nQ_n}+C_2Q_n+o(n^{-1}).
\]
[Ratio law coefficients](GAUSSIAN_CONE_AXIAL_INSTABILITY.md#a-uniform-expansion-for-radial-spectator-blocks)
[Identification of the radial loss](GAUSSIAN_CONE_AXIAL_INSTABILITY.md#exact-leading-loss-for-the-entire-radial-profile-class)

By the quadratic identity \(D_n-B_n(\bar h_n)=\int
B_n(h_{n,\omega}-\bar h_n)\,d\sigma(\omega)\), with \(B_n\)
extended to signed tuples in this formula. The expansion (AX2)
gives \(D_n-B_n(\bar h_n)=W_n+o(n^{-1})\). Indeed the coefficient
in the difference is the integral of \(J_r(x,y)\) against the
averaged products of the two centered profiles at \(x,y\).
Its absolute value is bounded by the integral of \(|J_r(x,y)|\)
against the mean absolute centered profile at \(x\). The latter
tends to zero in Gaussian \(L^1\) and is uniformly bounded,
so this integral tends to zero. The uniform \(O(n^{-3/2})\)
remainder also applies to these bounded signed tuples.
Substitute both bounds into (AS1) to obtain (PR1).

By Young's inequality, the mixed term in (PR1) can be absorbed
into half of \((1-\theta)W_n\) and a fixed multiple of \(Q_n\).
The strict positivity of \(C(A)\) for nonsimplex \(A\) then gives
the asserted necessary slice mass variance for deficit \(o(n^{-1})\).
This estimate leaves open the case of slice mass fluctuations
at or above the \(n^{-1/2}\) scale in root mean square.
The separate first-order comparison excludes an \(o(n^{-1})\)
deficit at a nonsimplex limit for every slice mass profile.
[First-order comparison for unrestricted cones](GAUSSIAN_CONE_DECONING_ASYMPTOTICS.md)
