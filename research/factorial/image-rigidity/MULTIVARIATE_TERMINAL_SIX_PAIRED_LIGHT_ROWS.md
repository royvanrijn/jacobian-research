# Paired ray projections for shared rows at terminal six

Suppose that the two ordinary exponent columns are distinct and every
unequal row has positive exponents. Write \(\ell_i=\min(a_i,b_i)\) and
\(d_i=|a_i-b_i|\) in each unequal row. Let \(u\) be the order-five
outer geometric mean after removing the binomial coefficients and
normalizing the constant and leading coefficients to one. If
\[
 \frac7{25}\le u\le\frac35,
 \qquad \max_i\frac{d_i}{5\ell_i+1}\le\frac18 ,
\]
then the order-five and order-six polynomial root arguments are
disjoint under independent positive rescalings. In particular the
ordinary factorial moment polynomials are coprime over every
characteristic-zero field. The
[exact rational checker](verify_multivariate_terminal_six_paired_light_rows.py)
certifies the full closed interval using the same cap \(1/8\) throughout.

## Coefficient intervals

Use the ordinary Gauss factors from the
[terminal-five normalization](MULTIVARIATE_TERMINAL_FIVE_ANGULAR_SEPARATION.md#the-normalized-coefficient-ranges)
and the [terminal-six coefficient notation](MULTIVARIATE_TERMINAL_SIX_COEFFICIENT_LINKS.md).
The order-five and order-six normalized multiplier arrays are
\[
 (1,u e^\alpha,v e^\chi,v e^{-\chi},u e^{-\alpha},1),\qquad
 (1,w e^\beta,y e^\eta,z,y e^{-\eta},w e^{-\beta},1).
\]
Put \(U=-\log u,\ W=-\log w,\ \rho=W/U\) and
\[
 V=\log(v/u^{3/2}),\quad H=\log(y/w^{8/5}),\quad
 J=-\log z+\tfrac32\log y-\tfrac35\log w .
\]

In an ordinary Gauss factor indexed by \(1\le h\le d\), write
\(x=d/(5\ell+h)\), \(x'=d/(6\ell+h)\). Then
\[
 x'=\frac{5x}{6-xh/d}\le\frac{5x}{6-x},
 \qquad
 1+\frac{25}{6}x\le 1+5x'\le\frac{1+4x}{1-x/6}.
\]
The cap on \(x\) is \(t_5=1/8\), and the cap on \(x'\) is
\(t_6=5/47\). For the lower bound on \(\rho\), write
\(\log(1+(25/6)x)=f(\log(1+4x))\), where
\(f(s)=\log(1+(25/24)(e^s-1))\) is concave and \(f(0)=0\).
Its quotient \(f(s)/s\) decreases, so its value at \(t_5\) gives
a lower bound on the ratio of the two outer logarithms. The upper
bound comes from \(1-\log(1-x/6)/\log(1+4x)\), which increases
because \(-\log(1-x/6)\) is convex, \(\log(1+4x)\) is concave,
and both vanish at zero. These bounds sum over the factors.

The other intervals use the single-factor convexity formulas in the
[quintic skew proof](MULTIVARIATE_TERMINAL_FIVE_ANGULAR_SEPARATION.md#the-normalized-coefficient-ranges),
the [sextic skew and defect proof](MULTIVARIATE_TERMINAL_SIX_ONE_PRIVATE_SIDE.md#coefficient-envelopes),
and the [sextic finite differences](MULTIVARIATE_TERMINAL_SIX_COEFFICIENT_LINKS.md#within-order-bounds).
With \(g_h=\log(1+hx_j)\), the local quantities used for the
nonnegative magnitudes at the two ends of each reflected pair are
\[
\begin{array}{c|c}
 |\alpha| & [2(g_1+g_2+g_3)-3g_4]/10\\
 |\chi| & (4g_2-g_1-g_3-g_4)/10\\
 |\beta| & (g_1+g_2+g_3+g_4-2g_5)/6\\
 |4\beta-5\eta| & [3(g_1+g_4)-2(g_2+g_3)-g_5]/2 .
\end{array}
\]
These magnitudes, \(V\), and \(H\) are convex in the local outer
logarithm \(s=\tfrac12\log(1+(j-1)x_j)\) and vanish at zero.
Thus the ratio to \(s\) increases, and evaluating it at \(t_j\)
bounds the sum by a multiple of the total outer logarithm. Reversing
a row changes only the signs of its skews. The triangle inequality
gives the bound for \(\eta\) from those for \(\beta\) and \(4\beta-5\eta\).

Two links use the derivative bounds directly. A single quintic factor
has \(\alpha-2\chi=\Delta^4 g_0/10\), so its absolute value is bounded
by \(3x^4/5\). A single sextic factor has \(J=\Delta^5g_0/20\), so
\(0\le J\le6(x')^5/5\). These bounds follow by integrating the fourth
and fifth derivatives of \(\log(1+tx_j)\) over unit cubes.
For any \(k\ge1\), concavity of the logarithm gives
\[
 \sum x_j^k\le
 \frac{2t_j^k}{\log(1+(j-1)t_j)}\sum s .
\]
Consequently the checker verifies all of the following rational
bounds, with lower bounds zero for the defects and magnitudes.
\[
\begin{array}{c|c@{\qquad}c|c}
 \rho &[1.034,1.052] & &\\ \hline
 |\alpha|/U&.05050&|\beta|/W&.05666\\
 |\chi|/U&.02510&|\eta|/W&.04563\\
 |\alpha-2\chi|/U&.00080&|4\beta-5\eta|/W&.00147\\
 V/U&.00509&H/W&.00578\\
 &&J/W&.00009
\end{array}
\]
The logarithm comparisons use the positive power series for
\(2\operatorname{artanh}((s-1)/(s+1))\) and a geometric tail bound.

## Preserve the paired contributions in each projection

For an interval of \(u\), let \(u_0\) be its midpoint and write
\(e_{j,k}=k(j-k)/(j-1)\). The reference binomial coefficient
multipliers are \(u_0^{e_{5,k}}\) and \(u_0^{(25/24)e_{6,k}}\).
Relative to these reference arrays, the log corrections are
\[
\begin{aligned}
 \sigma_5={}&m_5(e_{5,k})_{k=0}^5
   +(0,\alpha,V+\chi,V-\chi,-\alpha,0),\\
 \sigma_6={}&m_6(e_{6,k})_{k=0}^6
   +(0,\beta,H+\eta,\tfrac32 H-J,H-\eta,-\beta,0),\\
 m_5={}&\log(u/u_0),\\
 m_6={}&\tfrac{25}{24}\log(u/u_0)+(\rho-\tfrac{25}{24})\log u .
\end{aligned}
\]
The weighted outer logarithm contributes one shared scalar per
order. In the quintic use \(\alpha=2\chi+(\alpha-2\chi)\) so that
one signed variable \(\chi\) controls both reflected coefficient
pairs. In the sextic use
\(\beta=\tfrac54\eta+(4\beta-5\eta)/4\), and keep the same \(H\)
in the paired and central coefficients.

If \(|\sigma_{j,k}|\le A_{j,k}<1\), write
\[
 \exp(\sigma_{j,k})=1+\sigma_{j,k}+E_{j,k},
 \qquad 0\le E_{j,k}\le\frac{A_{j,k}^2}{2(1-A_{j,k})}.
\]
The lower bound is the tangent inequality and the upper bound
comes from the exponential series. Enclose the remainder
coefficientwise. In the linear term collect each shared variable's
contributions before evaluating its interval.

The four ray cosines are continuous piecewise affine functions
of \(u\) with the rational knots in the checker. For each ray
\(t=r e^{i\theta}\), the checker tries the projections
\[
 \frac{\Im(e^{-ih\theta}P_j(r e^{i\theta}))}{\sin\theta}
 =\sum_{k=0}^j c_{j,k}r^k\frac{\sin((k-h)\theta)}{\sin\theta},
 \qquad 0\le h\le j .
\]
The ratios on the right are bounded by the sine recurrence on
each cosine interval. For each shared variable, compute lower
and upper Bernstein coefficient arrays of its projection
polynomial in \(r\) on \([1/8,8]\). The change of basis has
nonnegative entries for intervals with nonnegative left endpoint,
as does subdivision. At each Bernstein index first multiply that
group's interval by the interval of its shared scalar, then add
the contributions. If one projection is certified strictly
positive or negative at every Bernstein index, that radial cell
contains no root. Otherwise bisect in \(r\), or bisect the
\(u\) interval if the radial test reaches its depth limit.
Each accepted cell has a strict rational sign certificate.
The depth limit on \(u\) is an error condition.

## Root counts

Each actual normalized array is log convex before restoring the
binomial coefficients. Hence its interior multipliers are at
most \(\max(ue^{|\alpha|},we^{|\beta|})\). The preceding bounds give
both skews less than \(2/25\), while \(0<w\le u\le3/5\).
Since \(e^x\le(1-x)^{-1}\) for \(0\le x<1\), the maximum is at
most \(15/23\). The reference arrays and their straight
coefficient segments satisfy the same bound. In both orders,
\[
 \frac{15}{23}\big((1+\tfrac18)^6-1\big)+(\tfrac18)^5<1.
\]
The constant and leading coefficients therefore dominate at
\(1/8\) and \(8\), respectively. All roots are in the checked
annulus.

At fixed \(u,\rho\) deform each actual array by a straight
coefficient segment to the symmetric array with those same
outer means and \(V=H=J=\alpha=\beta=\chi=\eta=0\). The coefficient
models are convex enlargements of these arrays, since the
certified expressions are linear in all shared variables and
independent remainder intervals. These segments stay in the
certified models. Next move \(\rho\) to \(25/24\), and then
move \(u\) to \(7/25\) using the symmetric arrays and the
continuous rays. The ray exclusions and annulus bounds preserve
root counts along the entire deformation.

At the anchor put \(X=t+t^{-1}\). The symmetric quintic,
after removing \(t+1\) and the reciprocal power, has quotient
\(B(X)=X^2+(5u-1)X+10v-5u-1\). The sextic quotient is
\(C(X)=X^3+6wX^2+(15y-3)X+20z-12w\).
The checker evaluates their signs at \(X=2\cos\theta\) on the
four anchor rays. In increasing upper angle the signs are
\((+,-,-,+)\) for \(B\) and \((-,-,+,+)\) for \(C\).
Also \(C(2)>0>C(-2)\). The quotient degrees give the unique
sector counts, with upper roots in the order
\[
 \text{sextic},\ \text{quintic},\ \text{sextic},\
 \text{quintic},\ \text{sextic}.
\]
The quintic has the additional root \(-1\). When continuing the
counts, join the first and last sectors with their respective lower
reflections. The sole quintic root in the last joined sector
stays real by conjugation. The actual sextic is strictly positive
on the real line by the ordinary factorial Gamma integral of an
even power of a nonzero real binomial. Its two roots in the last
joined sector are therefore nonreal. This proves disjointness of
arguments for the actual arrays. Positive normalizations preserve
arguments and rational Bézout transfers the coprimality conclusion.

Run the checker from the repository root with
```sh
python3 research/factorial/image-rigidity/verify_multivariate_terminal_six_paired_light_rows.py
```
