# Exact curvature signs for ordinary factorial row mergers

An ordinary row merger replaces several exponent rows
\((a_i,b_i)\) by their entrywise sum. The normalized moment coefficients
can move in either direction. There is an exact pointwise condition
for the change in log curvature, and an explicit order beyond which
every nonproportional shared merger flattens the normalized array.

## The pointwise sign

At any fixed fractional coefficient index, write the factorial
arguments in the rows being merged as \(x_i\ge0\), and their signed
slopes as \(d_i\). Put
\[
 f(x)=\psi'(1+x),\qquad h(x)=1/f(x),\qquad
 A_i=h(x_i),\qquad B=\sum_i A_i,\qquad
 C=h\!\left(\sum_i x_i\right),\qquad D=B-C .
\]
For a merger of at least two rows, \(D>0\). With \(d=\sum_i d_i\),
the change in second derivative of the log moment multiplier is
\[
 \begin{aligned}
 \Delta
 &=d^2 f\!\left(\sum_i x_i\right)-\sum_i d_i^2 f(x_i)\\
 &=\frac{d^2 D}{BC}
      -\sum_i\frac{(d_i-A_i d/B)^2}{A_i}.
 \end{aligned}                                                      \tag{1}
\]
If the factorial arguments are \(x_i=n(a_i+d_i t)\), with
\(d_i=b_i-a_i\) the exponent gaps, multiply (1) by \(n^2\) to
obtain the second derivative in \(t\).
Endpoint normalization removes an affine function of \(t\),
so the same second-derivative formula holds for normalized
log coefficients. Equation (1) follows by expanding the squares.

To prove \(D>0\), first \(h(x)>x\), by the integral comparison
\(\psi'(1+x)<1/x\) for \(x>0\), with the case \(x=0\) immediate.
Also \(h'(x)<1\). Here is a short proof of the latter fact.
For \(g(t)=t/(1-e^{-t})\), extended by \(g(0)=1\), we have
\[
 (\log g)''(t)=-\frac1{t^2}+\frac1{(2\sinh(t/2))^2}<0 .
\]
Thus \(g(s)g(t-s)>g(t)\) for \(0<s<t\), by strict concavity at
the endpoints \(0,t\). Integrating in \(s\) gives
\((g*g)(t)>t g(t)\). Taking Laplace transforms in the
polygamma integrals proves
\(\psi'(z)^2+\psi''(z)>0\) for \(z>0\), which gives \(h'(x)<1\).
The same trigamma inequality appears in stronger form in
[Qi's completely monotonic degree theorem](https://doi.org/10.3934/math.2020219).
It follows that \(h(x+y)-h(x)<y<h(y)\) for \(y>0\);
at \(y=0\), use \(h(0)>0\). Iterating proves \(C<B\).

If two rows have opposite signed gaps or one has gap zero,
then \((d_1+d_2)^2\le d_1^2+d_2^2\).
Since \(f(x_1+x_2)\le\min(f(x_1),f(x_2))\), (1) is
nonpositive, and is strictly negative when both arguments
are positive and at least one gap is nonzero. Any balanced
batch with \(d=0\) also has strictly negative \(\Delta\) if
one gap is nonzero.

In contrast, for proportional shared rows all the ratios
\(d_i/x_i\) agree at each point. In that case (1) is positive
when the common gap direction is nonzero. Indeed
\(x f(x)\) increases strictly on \((0,\infty)\).
With \(K(v)=v/(e^v-1)\), this follows by writing
\[
 x f(x)=\int_0^\infty e^{-u}K(u/x)\,du
\]
and using strict decrease of \(K\).
Consequently \(x^2 f(x)\) is strictly superadditive on positive
arguments, as needed after substituting \(d_i=\lambda x_i\).
The [joint factorial curvature representation](MULTIVARIATE_JOINT_FACTORIAL_CURVATURE.md)
gives a complementary comparison across proportional operations.

For two rows with \(d_1+d_2\ne0\), the positive-change region
can equivalently be written as the explicit interval
\[
 \left|\frac{d_1}{d_1+d_2}-\frac{A_1}{A_1+A_2}\right|
 <\frac{\sqrt{A_1A_2 D/C}}{A_1+A_2}.
 \tag{2}
\]
For example, since \(0<h(x)-x\le h(0)<1\), its radius is
at most \(1/\sqrt{2(x_1+x_2)}\) when \(x_1+x_2>0\).
Its center differs from \(x_1/(x_1+x_2)\) by less than
\(1/(x_1+x_2)\). These bounds exhibit the narrowing of the
positive-change interval at larger factorial arguments.

## An explicit flattening threshold for nonproportional rows

Now suppose all merged rows are shared. Write
\[
 u_i(t)=a_i+d_i t,\qquad
 L_i=\min(a_i,b_i)>0,\qquad U_i=\max(a_i,b_i),\qquad
 D_{ij}=d_i a_j-d_j a_i .
\]
If the rows are nonproportional, the quantities
\[
 \mathcal D=\sum_{i<j}
      \frac{D_{ij}^2}{U_iU_j\sum_\ell U_\ell},
 \qquad
 \mathcal E=\sum_i\frac{d_i^2}{L_i^2}
\]
satisfy \(\mathcal D>0\). For every real order
\(n>\mathcal E/\mathcal D\), merging these rows strictly raises
every interior endpoint-normalized coefficient multiplier at
all fractional indices \(0<t<1\).

Indeed integral comparison gives
\[
 \frac1{v+1}\le f(v)\le\frac1v\qquad(v>0).
\]
Using these inequalities in the first line of (1) with
\(x_i=n u_i(t)\) and multiplying by \(n^2\) gives
\[
 \begin{aligned}
 n^2\Delta
 &\le n\left(\frac{(\sum_i d_i)^2}{\sum_i u_i}
                  -\sum_i\frac{d_i^2}{u_i}\right)
       +\sum_i\frac{d_i^2}{u_i^2}\\
 &=-n\sum_{i<j}
     \frac{D_{ij}^2}{u_i u_j\sum_\ell u_\ell}
       +\sum_i\frac{d_i^2}{u_i^2}
 \le -n\mathcal D+\mathcal E.
 \end{aligned}                                                       \tag{3}
\]
The equality expands the weighted variance, using
\(d_i u_j(t)-d_j u_i(t)=D_{ij}\).
The difference of the normalized log multipliers is then
strictly concave and vanishes at both endpoints, so is
strictly positive inside the interval. For proportional
shared rows the pointwise sign above is reversed and the
merged interior multipliers are strictly smaller whenever
the gap direction is nonzero.

The same weighted-variance identity gives an exact sign at
the level of the leading stationary curvature
\(\rho(t)=t(1-t)\sum_i d_i^2/u_i(t)\). A nonproportional
shared merger reduces \(\rho(t)\) at every \(0<t<1\); a
proportional merger preserves it. In particular, for fully
shared supports a nonproportional merger strictly reduces
the maximum of \(\rho\) on \([0,1]\), including when the
merged rows have zero total gap. This gives a one-direction
comparison for the large-order saddle transition parameter.

These comparisons describe coefficient curvature and the
leading stationary geometry. The corresponding transfer of
coprimality for arbitrary ordinary row mergers requires
additional information about the roots.
