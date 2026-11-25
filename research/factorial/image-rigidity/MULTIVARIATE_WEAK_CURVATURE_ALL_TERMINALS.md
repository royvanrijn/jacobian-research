# An explicit weak-curvature bound at every terminal

Let \(a,b\) be distinct nonnegative integer exponent vectors, and let
\(X_i\) be independent unit-rate Gamma variables of positive shapes
\(\gamma_i\). Write \(U=X^a,\ V=X^b\), \(d_i=b_i-a_i\), and
\(P_j(t)=\mathbb E(U+tV)^j\). Indices with \(d_i=0\) can be removed
since they contribute a nonzero constant to each polynomial. For
\(N=n+1\ge3\) put
\[
 \lambda_i(j)=\gamma_i+\tfrac j2(a_i+b_i),\qquad
 V_j=\sum_i\frac{d_i^2}{\lambda_i(j)},\qquad
 \sigma=\sqrt{V_n}.
\]

**Theorem.** If
\[
 \sigma\le\frac{1}{3\cdot10^7N^{7/2}},
 \tag{1}
\]
then \(P_n\) and \(P_N\) have disjoint sets of root arguments.
Each of the two polynomials has simple roots. In particular, for
ordinary factorial moments in arbitrary dimension, the terminal
pair is coprime over every characteristic-zero field if
\[
 S=\sum_i\frac{d_i^2}{\max(a_i,b_i)}
 \ \le\ \frac{1}{10^{16}N^6}.
 \tag{2}
\]
The constants apply at every terminal. They depend neither on the
dimension nor on the exponent profile or the Gamma shapes.

This strengthens the fixed-terminal weak-curvature limit by making its
dependence on \(N\) explicit. Centering at the midpoint tilt and using the
logarithm of the coefficient ratio make the angular comparison independent
of the two radial centers. [Fixed-terminal curvature proof](MULTIVARIATE_FIXED_TERMINAL_CURVATURE.md)

## Quantitative comparison of the Hermite pairs

Write \(\mathrm{He}_j\) for the monic probabilists' Hermite polynomial
defined by \(\exp(xt-t^2/2)=\sum_j\mathrm{He}_j(x)t^j/j!\).
For any \(\rho\in[n/N,1]\) consider the roots of
\(\mathrm{He}_n\) and the roots of \(\mathrm{He}_N\) scaled by
\(\sqrt{\rho}\). The two sets have distance at least
\[
 \delta=\frac{1}{200 N^{5/2}}.
 \tag{3}
\]
Within each set every two distinct roots have distance at least
\(b_0=1/(20\sqrt N)\).

The recurrence \(\mathrm{He}_{j+1}=x\mathrm{He}_j-j\mathrm{He}_{j-1}\)
identifies the roots with the eigenvalues of the size-\(j\)
tridiagonal matrix with zero diagonal and off-diagonal entries
\(\sqrt{k}\) for \(1\le k<j\). They are simple because an eigenvector
is determined by its first entry, and \(|x|\le2\sqrt j\) by the row
sum bound. For spacing, put \(x=x_0+t/\sqrt j\) in
\(\mathrm{He}_j''-x\mathrm{He}_j'+j\mathrm{He}_j=0\) and normalize
the solution so that its initial data are \((0,1)\) at a root.
The scaled equation is \(g''=(x_0/\sqrt j+t/j)g'-g\).
On \(|t|\le1/10\) with \(|x_0|\le2\sqrt j\) the norm
\(|g|+|g'|\) is at most \(e^{4|t|}<3/2\), by integrating the
first-order system. The integral equation for \(g'\) then gives
\(|g'(t)-1|\le(31/10)(3/2)|t|<1/2\). Thus a second root is more
than \(1/(10\sqrt j)\) away. With initial data \((1,0)\), the
same bound gives \(|g(t)-1|\le(3/2)|t|\). This excludes roots at
distance at most \(1/(10\sqrt j)\) from any stationary point
in \([-2\sqrt j,2\sqrt j]\) with nonzero value. These two estimates
also give the stated spacing after scaling the larger root set.

Order the positive roots decreasingly as \(x_k\) for the lower
degree and \(y_k(\rho)\) for the scaled higher degree, using zero
or the relevant negative root if an upper neighbor is nonpositive.
The leading principal block of the larger matrix is the lower
matrix scaled by \(\sqrt\rho\), so the inner neighbor
\(y_{k+1}(\rho)\) lies below \(x_k\).
If this neighbor is positive, its distance below \(x_k\) is at
least the distance for \(\rho=1\). Each \(x_k\) is a stationary
point of \(\mathrm{He}_N\) because
\(\mathrm{He}_N'=N\mathrm{He}_n\). Its value there is nonzero
by principal interlacing and the recurrence. The preceding
stationary-point bound applies.

For the outer neighbor, the trailing principal block has
off-diagonal entries \(\sqrt{\rho(k+1)}\). Relative to the lower
matrix their ratios decrease with \(k\) and are bounded below by
\[
 c=\sqrt{1+\frac{1}{N(n-1)}}\ \ge\ 1+\frac1{3N^2}.
\]
Dividing these ratios by \(c\), the zero-diagonal Jacobi
comparison bounds each positive eigenvalue of the trailing
block below by \(c x_k\). That comparison follows by writing
a decreasing sequence of edge ratios at least one as products
of adjacent vertex weights at least one and applying the
singular value comparison to the bipartite matrix.
Principal interlacing therefore gives \(y_k(\rho)-x_k\ge
x_k/(3N^2)\). Symmetry and same-set spacing give \(x_k\ge
b_0/2\). This proves (3) for the outer neighbors. For an inner neighbor at zero or on the other side of zero,
the gap is at least \(x_k\ge b_0/2\). If the smaller degree
has the zero root, the larger root set is even, with its least
positive root at least \(b_0/2\) by symmetry and spacing.
Reflection covers the negative roots. [Zero-diagonal Jacobi comparison](MULTIVARIATE_MIRRORED_ROOT_CIRCLE.md#consecutive-interlacing-for-a-single-mirrored-unit-gap)

## Centering the exact moments

Let \(\mathbb E_j\) denote expectation with independent shapes
\(\lambda_i(j)\). Under this measure put
\[
 L_j=\sum_i d_i\log X_i,\qquad
 \mu_j=\mathbb E_j L_j,\qquad W_j=(L_j-\mu_j)/\sigma,\qquad
 \rho_j=V_j/\sigma^2 .
\]
The reference variances satisfy \(\rho_n=1\) and
\(n/N\le\rho_N\le1\), because
\(\lambda_i(n)\le\lambda_i(N)\le (N/n)\lambda_i(n)\).
For \(z=-\exp(-\mu_j+\sigma x)\), the identity
\[
 U+zV=-2(UV)^{1/2}e^{(-\mu_j+\sigma x)/2}
                  \sinh\bigl(\sigma(x+W_j)/2\bigr)
\]
gives, after removing a nonvanishing factor, the entire function
\[
 F_j(x)=\mathbb E_j\!\left[(x+W_j)^j
        D(\sigma(x+W_j))^j\right],\qquad
 D(y)=\frac{\sinh(y/2)}{y/2}.
 \tag{4}
\]
Indeed the original weight multiplied by \((UV)^{j/2}\)
is a constant times the midpoint Gamma weight.

The function \(D(y)\) is the moment generating function of a
uniform variable on \([-1/2,1/2]\). Write \(Y_j\) for the sum
of \(j\) such independent variables, independent of \(W_j\).
With \(h=\sigma Y_j\) and
\[
 Q_{j,h}(x)=\mathbb E_j\!\left[e^{h W_j}(x+W_j)^j\right],
\]
equation (4) is exactly \(F_j(x)=\mathbb E[e^{hx}Q_{j,h}(x)]\).
Here \(|h|\le\sigma N/2\). The expectations can be interchanged
since the Gamma moments with bounded tilts are finite in the
open strip of their shapes. The bounds below place the entire
closed range of \(h\) inside that strip.

## A relative error bound on Hermite circles

Let \(K_j(z)=\log \mathbb E_j e^{zW_j}\) with \(K_j(0)=0\).
Since the integer gaps satisfy \(|d_i|\ge1\),
\[
 \frac{|d_i|}{\lambda_i(j)\sigma}\le\sigma,\qquad
 \frac{1}{\lambda_i(j)}\le\sigma^2 .
\]
This follows by bounding each \(d_i^2/\lambda_i(j)\) by
\(\sigma^2\). The Gamma moment formula gives Taylor coefficients
of \(K_j\) from the polygamma series. The integral bounds on
\(\psi'(\lambda)=\sum_{k\ge0}(\lambda+k)^{-2}\) give
\(1/\lambda\le\psi'(\lambda)\le1/\lambda+1/\lambda^2\).
The quadratic coefficient of \(K_j\) therefore differs from
\(\rho_j/2\) by at most \(\sigma^2/2\).
For degree \(s\ge3\), differentiating the same series bounds
the coefficient in absolute value by
\[
 \frac{1}{s}\sum_i\left|\frac{d_i}{\sigma}\right|^s
                         \sum_{k\ge0}(\lambda_i(j)+k)^{-s}
 \ \le\ \frac{2}{s}\sigma^{s-2}.
\]
We used \(\sum_i d_i^2\psi'(\lambda_i(j))/\sigma^2\le1+\sigma^2\le2\).
The Taylor series is analytic for \(\sigma|z|<1\), where all
shifted Gamma shapes have positive real parts. Consequently
\[
 \left|K_j(z)-\tfrac12\rho_j z^2\right|
 \le \tfrac12\sigma^2|z|^2+\tfrac43\sigma|z|^3
 \qquad(\sigma|z|\le1/2).
 \tag{5}
\]
The second inequality follows by summing the geometric
majorant from degree three. [DLMF, trigamma series](https://dlmf.nist.gov/5.15.E1)

For \(Z\) standard real Gaussian define
\(G_j(x)=\mathbb E(x+\sqrt{\rho_j}Z)^j\).
Its roots are the imaginary rotations of the two Hermite root
sets compared above. Use circles of radius \(r=\delta/4\)
around these roots. For a point \(x\) on a circle of its own
order, put \(A=1/r=800N^{5/2}\). The sum of the reciprocal
distances to all roots other than the center is at most
\(80\sqrt N(1+\log N)\) by the spacing bound. Indeed
a root \(k\) places away on the imaginary axis is at distance
at least \(k b_0/2\), and the harmonic sum on each side is at
most \(1+\log N\). The logarithmic derivative of \(G_j\)
therefore satisfies
\[
 \left|\frac{G_j'(x)}{G_j(x)}\right|
 \le A+80\sqrt N(1+\log N).
\]
To bound higher derivatives, differentiate the Hermite equation
\(\rho_j G_j''+xG_j'-jG_j=0\), which also follows by comparing
coefficients of the generating function. At derivative order \(k\),
\[
 \rho_j G_j^{(k+2)}+xG_j^{(k+1)}-(j-k)G_j^{(k)}=0.
\]
On the root circles \(|x|\le3\sqrt N\) and \(\rho_j\ge2/3\).
Taking \(D_0=6\sqrt N\), the normalized coefficients of the
two lower derivatives in this recurrence have absolute values
at most \(3/4\) and \(1/24\). Since \(1+\log N\le N^2\), the
bounds at \(k=0,1\) and induction thus give
\[
 \left|\frac{G_j^{(k)}(x)}{G_j(x)}\right|
 \le 150N^2 D_0^k\qquad(0\le k\le j).
 \tag{6}
\]

Take \(R=2D_0=12\sqrt N\). For \(|z|\le R\) and \(|h|\le\sigma N/2\)
we have \(|z+h|\le2R\) and \(\sigma|z+h|\le1/2\) by (1).
The points on the root circles satisfy \(|x|\le3\sqrt N<R\).
Using (5) and \(N\le R^2\),
\[
 \left|hx+K_j(z+h)-\tfrac12\rho_j z^2\right|
 \le 24\sigma R^3\le1/2 .
\]
For example the error term in (5) is at most \(16\sigma R^3\),
and \(|hx|+\rho_j|hz+h^2/2|\le3|h|R\le2\sigma R^3\).
If \(c_k\) are the Taylor coefficients in \(z\) of
\(\exp(hx+K_j(z+h)-\rho_j z^2/2)-1\) at this fixed \(x,h\),
Cauchy's estimate gives \(|c_k|\le48\sigma R^3/R^k\).
Extracting the coefficient at degree \(j\) of
\(\exp(xz+\rho_j z^2/2)\) times this series gives
\[
 e^{hx}Q_{j,h}(x)-G_j(x)=\sum_{k=0}^j c_k G_j^{(k)}(x).
\]
Equations (6) and \(R=2D_0\) now give, also after averaging in \(h\),
\[
 \left|\frac{F_j(x)}{G_j(x)}-1\right|
 \le48\sigma R^3(1+150N^2)
 \le48\cdot12^3\cdot151\,\sigma N^{7/2}<\tfrac12 .
\]
The initial term uses \(G_j/G_j=1\) and the remaining terms
are bounded by (6), since \(\sum_{k\ge1}(D_0/R)^k=1\).

Rouché's theorem places exactly one zero of \(F_j\) in each
circle's interior, counted with multiplicity. All these zeros
satisfy \(|x|<3\sqrt N\), so the exponential map in (4) is
injective on their union because \(6\sigma\sqrt N<2\pi\).
They thus give all \(j\) roots of \(P_j\), which are simple.
Since both orders use the same scale \(\sigma\), their arguments
are distinct by (3). The midpoint means \(\mu_j\) are real and
contribute only radial factors.

Finally, \(\sigma^2\le2S/n\) since
\(\lambda_i(n)\ge(n/2)\max(a_i,b_i)\). With \(n\ge 2N/3\), (2)
implies (1). For ordinary shapes the moment polynomials have
integer coefficients, so rational Bézout transfers coprimality
to every characteristic-zero field.
