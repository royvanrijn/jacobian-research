# First-sector radial comparison for comparable exponent columns

Let \(a_i,b_i\) be distinct nonnegative integer exponent columns in
any finite dimension, and suppose \(a_i\le b_i\) on every row. Write
\[
 P_m(z)=\sum_{k=0}^{m}\binom{m}{k}
                  \prod_i((m-k)a_i+kb_i)!\,z^k .
\]
For \(n\ge2\), \(N=n+1\), and every real \(\kappa\ge1\), the two
equalities \(P_n(z)=P_N(\kappa z)=0\) are impossible in
\(|\arg z|\le 2\pi/N\). In particular, comparable ordinary profiles
have no common physical root in the first sector, in every degree.
If they have roots on the same ray in that sector, the upper
root's radius must be strictly smaller in the increasing-gap
orientation. For decreasing gaps, reversing the columns
reciprocates the roots and reverses this radius comparison.

The proof verifies the
[conditional chord criterion](MULTIVARIATE_COMPARABLE_ROWS_FIRST_SECTOR_COMMON_ROOTS.md#a-conditional-moment-comparison)
without a restriction on the number of rows. It uses a
completely monotone ratio of the two densities of the monomial
ratio.

## Dilation of a hyperbolically completely monotone density

A positive function \(f\) is hyperbolically completely monotone
if \(f(uv)f(u/v)\) is completely monotone as a function of
\(v+v^{-1}\) for every \(u>0\). Such a function is a positive
pointwise limit of finite products of the form
\(Cx^\gamma\prod_j(1+c_jx)^{-b_j}\), with \(c_j,b_j>0\).
The class of distributions with these densities is closed
under independent products and powers of absolute value at
least one. These are the approximation and closure results
of Bondesson, Chapter 5, including Theorem 5.1.1 and the
power statement on p. 69. Section 2 of Bondesson and Simon
also states the approximation and convolution properties.
[Bondesson's chapter](https://link.springer.com/chapter/10.1007/978-1-4612-2948-3_5),
[Bondesson and Simon](https://doi.org/10.1090/tran/7123).

For such \(f\) and \(c\ge1\), the ordinary function
\(x\mapsto f(cx)/f(x)\) is completely monotone. Indeed on an
elementary factor in the approximating product, up to a
constant, the quotient is \(L(x)=((x+A)/(x+B))^b\) with
\(A\ge B>0\). Its negative logarithmic derivative is
\[
  -(\log L)'(x)=b\big[(x+B)^{-1}-(x+A)^{-1}\big]\ge0 ,
\]
and this derivative is completely monotone by the elementary
Laplace representation of the reciprocal. Induction in
\(L'=-L[-(\log L)']\) using the product rule proves complete
monotonicity of \(L\). The powers \(x^\gamma\) in the
approximation contribute only a positive constant \(c^\gamma\).
Products and finite positive pointwise limits of completely
monotone functions preserve complete monotonicity. The desired
quotient follows by taking the approximating limit in both
numerator and denominator.

Suppose \(X\) has a positive smooth hyperbolically completely
monotone density \(f\), and \(B\in(0,1]\) is independent of \(X\).
Assume \(BX\) has a positive smooth density \(f_0\). Then
\[
                R(x)=\frac{f_0(x)}{f(x)}
                  =\mathbb E\!\left[\frac{f(x/B)}{Bf(x)}\right]
\]
is completely monotone. The displayed equality holds first
almost everywhere by conditioning. Its right side is
nonincreasing, by the dilation result, and finite everywhere
since it is finite almost everywhere. It must agree everywhere
with the continuous left side. The alternating finite
differences of each integrand in \(x\) have the required signs.
Taking the expectation and then the difference-quotient
limits gives the derivative signs for the smooth function
\(R\). If both laws have first moments and \(\mathbb P(B<1)>0\),
the first moments differ, so \(R\) is nonconstant. Bernstein's
Laplace representation then gives \(R'(x)<0\) everywhere.

## The two ordinary exponential tilts

Represent the factorials using independent standard exponential
variables \(Y_i\). Put \(U=\prod_iY_i^{a_i}\),
\(V=\prod_iY_i^{b_i}\), \(d_i=b_i-a_i\) and \(T=V/U\).
Equal rows contribute positive factors and are dropped. Under
the probability tilt by \(U^j\), the \(Y_i\) are independent
unit-rate gamma variables of shapes \(j a_i+1\). Thus the
probability law of \(T\) under either tilt has a positive
smooth density on \((0,\infty)\), by convolution in the
logarithmic coordinates.

If \(a_i>0\), take \(Z_i\) gamma of shape \(Na_i+1\) and
\(B_i\) beta of parameters \(na_i+1,a_i\), independently.
Then \(B_iZ_i\) is gamma of shape \(na_i+1\). The identity
follows by changing from two independent same-rate gamma
variables to their sum and their proportion. If \(a_i=0\),
use \(B_i=1\). Make the pairs independent across rows.
Consequently the lower law is that of \(B X\), where
\[
     X=\kappa\prod_i Z_i^{d_i},\qquad
     B=\kappa^{-1}\prod_i B_i^{d_i}\le1 ,
\]
and \(X,B\) are independent. The law of \(X\) is the law
of \(\kappa T\) under the upper tilt. It has a hyperbolically
completely monotone density \(f\) by the product and power
closure properties, since gamma densities themselves are
in that class.

After division by positive constants, the two polynomials
under comparison now have the form
\[
      \int_0^\infty(1+zx)^n f_0(x)\,dx,\qquad
      \int_0^\infty(1+zx)^N g(x) f_0(x)\,dx,\qquad
      g(x)=f(x)/f_0(x).
\]
The reciprocal \(1/g\) is completely monotone. Either \(g\)
is constant or \(g'(x)>0\) and \(1/g\) is convex. The proof
of the conditional chord criterion applies to these two
integrals, taking its multiplier to be \(g(x)(1+zx)\).
This establishes the sector and radius comparison. At
\(\kappa=1\) the function \(g\) is also the conditional mean
\(\mathbb E_n[U\mid T=x]\) divided by \(\mathbb E_n[U]\).

## Mixed gaps and the conditional monotonicity condition

With a negative gap present, the density ratio of the two
tilts cannot be nondecreasing, even after any finite positive
dilation of the upper law. Indeed, the positive moment of
order \(s>0\) of \(T\) under the order-\(j\) tilt is finite
exactly when
\[
                 s<\min_{d_i<0}\frac{ja_i+1}{|d_i|}.
\]
Every \(a_i\) in the minimum is strictly positive, so the
threshold for \(j=N\) exceeds the one for \(j=n\).
Choose \(s\) strictly between the thresholds. The lower
\(s\)-moment is infinite and the upper one finite, a relation
unchanged by scaling. A nondecreasing density ratio would
order the upper law above the lower one for all increasing
functions, by the covariance identity for a monotone
likelihood ratio. Applying that order to the truncated
\(s\)-powers gives a contradiction. A profile with both
signs has this obstruction in either orientation.
