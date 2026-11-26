# Sharp stationary exclusion and unfolding for finite ordinary cores

Fix adjacent degrees \(n,N=n+1\), \(n\ge2\), a finite ordinary
two-monomial core, and the branch pair at a simple physical contact
with parameters \(W_*>0\), \(|\lambda_*|<W_*\), and integer tail
gap \(h\). Use the midpoint and root-response notation of the
[positive diffuse jet tests](MULTIVARIATE_POSITIVE_DIFFUSE_EXCEPTION_SIGNED_JETS.md).
Suppose the Gaussian angular difference \(g\) vanishes to even
order \(2r\), with leading coefficient \(a\ne0\), and that the
skew response \(A=\Im\Delta\kappa\) vanishes to order at least
\(r\), with coefficient \(b\) of order \(r\).
Put \(B=\Im\Delta\gamma(W_*)\),
\(C=\Im\Delta\sigma(W_*)\), \(D=\Im\Delta\beta^H(W_*)\),
and set \(M=W_*-|\lambda_*|>0\). Define
\[
 \mathcal T_M=\operatorname{conv}\{(0,0),(M/2,0),(M,M^2)\},
 \qquad
 E(u,v)=\operatorname{sgn}(a)
          \left[B+Cu+\left(D-\frac{b^2}{4a}\right)v\right].
\]

**Theorem.** If the minimum of \(E\) on \(\mathcal T_M\) is
strictly positive, the indicated branch pair cannot give an
ordinary common root in a neighborhood of the reference for any
sufficiently diffuse tail, with no bound on the tail row count.
If the minimum is strictly negative, the reference is a limit of
simple transverse interior physical Gaussian contacts of finite
ordinary cores. The combined references converge in normalized
coefficients and physical comparison scales. The adjacent orders
and the total integer degree difference are preserved. Any strict upper bound
on the total limiting curvature persists.

For a finite ordinary core, if \(B\ne0\) one of the two strict
tests always applies. All three signs are then effectively
decidable by evaluating the affine function at the three vertices.
The unfolding branch produces Gaussian references with ordinary
integer anchors. The exact common-root question for arbitrary
ordinary completions near those transverse contacts remains open.

## Uniform exclusion from the signed curvature allowance

Use the [signed curvature triangle](MULTIVARIATE_STATIONARY_BIAS_CURVATURE_TRIANGLE.md)
with \(Z=\sum_i(d_i/m_i)^2\), \(k=K/\sqrt Z\), and \(u=B_4/Z\).
For any candidate contact sequence approaching \((W_*,\lambda_*)\),
the possible limit points of \((u,k^2)\) belong to \(\mathcal T_M\).
The trimmed weights in that theorem give explicit errors
\(\eta\) and \(\sqrt Z\) for the two normalized moments. In
particular the bound is uniform over all row counts.

By the imaginary part of the
[exact root expansion](MULTIVARIATE_POSITIVE_DIFFUSE_EXCEPTION_SIGNED_JETS.md#root-coefficients-at-an-algebraic-exceptional-clock),
the bounded variable \(Y=(W-W_*)^r/\sqrt Z\) of any such
ordinary contact satisfies
\[
                     0=aY^2+bkY+B+Cu+Dk^2+o(1).
\]
The boundedness follows by absorbing the Taylor errors in
\(g(W)=a(W-W_*)^{2r}+O((W-W_*)^{2r+1})\) and
\(A(W)=O((W-W_*)^r)\), and using bounded \(k,u\).
The remainder is uniform with the same root isolations.
Completing the square yields the first assertion by the strict
margin on \(\mathcal T_M\). The coefficient enclosures and the
displayed inequalities make this a computable local criterion.

## Exact integer anchors and transverse unfolding

If the minimum is negative, continuity permits a point with
\(E(u,v)<0\) on \(\mathcal T_L\) for some \(0<L<M\).
Indeed all vertices at \(L\) tend to their counterparts at \(M\).
Set \(k^2=v\), choosing the sign so that \(bk/a\le0\).
The [ordinary realization of the triangle](MULTIVARIATE_STATIONARY_BIAS_CURVATURE_TRIANGLE.md)
gives finite shared anchors \(A_\varepsilon\) at mass \(L\) with
\[
 W_A\to L,\quad \lambda_A\to0,\quad
 Z_A/\varepsilon^2\to1,\quad K_A/\sqrt{Z_A}\to k,\quad
 B_{4,A}/Z_A\to u,\quad \eta_A\to0 .
\]
Here one may omit the signed-gap correction in the realization.
The signed anchor gap \(h_A\) is an integer and may grow.

Keep the anchor rows exact and vary the remaining Gaussian
strength \(T=W-W_A\) in a small neighborhood of \(W_*-L>0\).
The lower normalization includes the anchor endpoint factor.
The angular root functions come from multiplying the original
core coefficients at strength \(W\) by the exact anchor corrections
after its quadratic strength \(W_A\) is removed.
From the matched factorial expansion, these corrections in
log coefficients are
\(K_A H_j+Z_A G_j+B_{4,A}I_j+O_j(\eta_A Z_A)\).
The root perturbation expansion holds also with one \(W\)
derivative of its remainder. For this, isolate the reference
simple roots on a fixed complex \(W\)-disk about \(W_*\), and
apply the analytic root map to the exact coefficient
multipliers, which are independent of \(W\) after removal
of the quadratic strength. The Taylor remainder is uniformly
\(O(\eta_A Z_A)\) on a smaller disk by the same parameter
bounds as for real \(W\). Cauchy's estimate gives that bound
for the derivative on a further smaller disk.

Set \(W-W_*=Z_A^{1/(2r)}X\) and divide the angular
difference by \(Z_A\). On bounded real \(X\) intervals the
function and its first derivative converge to
\[
                   aX^{2r}+bkX^r+B+Cu+Dk^2 .                 \tag{1}
\]
Indeed the Gaussian and linear skew Taylor terms give the
first two terms of (1), and the other three Taylor coefficients
are evaluated at \(W_*\) up to errors tending to zero in both
the function and its \(X\) derivative. The completed-square
minimum of (1) after multiplying by \(\operatorname{sgn}(a)\)
is negative and occurs at \(X^r=-bk/(2a)\ge0\).
Its larger zero as a quadratic in \(X^r\) is strictly positive
and simple. Consequently (1) has a simple positive zero in
\(X\), and the exact angular difference changes sign with
nonzero derivative at nearby points for all sufficiently small
\(\varepsilon\).

At each such angular zero, solve the radial equality by
choosing the remaining bias and use remaining signed gap
\(h-h_A\). In the anchor scale, the combined physical logarithmic
argument multiplier is
\[
 h\log(N/n)-\frac{\lambda_A+\lambda_{\mathrm{rem}}}{2nN}
              +O_{n,N}(Z_A).
\]
This follows from the exact center expansion, since
\(|L_{2,A}|\le Z_A\) and the remaining part is exactly Gaussian.
The ratio of the root moduli tends to that of the original
contact. The equality of that ratio with the physical multiplier
then gives \(\lambda_{\mathrm{rem}}	o\lambda_*\). Its strict feasibility follows from
\(T\to W_*-L>|\lambda_*|\). This proves the physical and
integer-gap assertions. The total curvature of the anchor
using larger exponents differs from \(W_A\) by \(O(\eta_A W_A)\).
Thus the limiting total curvature agrees with the original
one and every strict cap persists.

Finally, because \(v_*=e^{-\lambda_* /(2nN)}\),
\[
          M=-2nN\log\!\left(q_*e^{|\lambda_*|/(2nN)}\right).
\]
The argument is a positive algebraic number different from one
by the interior condition and the algebraicity of \(q_*,v_*\).
The [algebraic clock and bias test](MULTIVARIATE_DIFFUSE_PHYSICAL_EXCEPTION_POINTS.md)
and Hermite–Lindemann make \(M\) transcendental. The three
vertex values are evaluations of algebraic polynomials in \(M\)
having nonzero constant term if \(B\ne0\). Each such value
is nonzero and rational interval refinement decides its sign.

