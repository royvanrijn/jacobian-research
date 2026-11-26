# Absolute diffuse angular separation at every curvature and degree

For distinct nonnegative integer exponent columns \(a,b\), omit equal
rows and set
\[
       S=\sum_i\frac{(b_i-a_i)^2}{\max(a_i,b_i)},\qquad
       \eta=\max_i\frac{|b_i-a_i|}{\max(a_i,b_i)} .
\]

**Theorem.** There is one computable positive rational number
\(\eta_0\le 1/2\) such that every ordinary profile with
\(\eta\le\eta_0\) has disjoint consecutive moment root arguments
at every degree, at arbitrary \(S>0\) and with arbitrary row
orientations. The tolerance is absolute and unevaluated.
Consequently these ordinary profiles satisfy SFC(2) over every
characteristic-zero field.

There is also an absolute positive unevaluated tolerance \(\eta_1\)
with the same conclusion for independent blends on the two
consecutive grids. More precisely, if \(L_j\) is the
endpoint-normalized factorial log symbol per degree, \(N=n+1\),
and \(\lambda_n,\lambda_N\in[0,1]\) are independent, use
\[
   \widehat L_j=(1-\lambda_j)L_n+\lambda_j L_N,\qquad
   P_j(z)=\sum_{k=0}^j\binom jk e^{j\widehat L_j(k/j)}z^k,
                \quad j=n,N .
\]
When \(\eta\le\eta_1\), these two polynomials have disjoint root
arguments for every pair of orders and every such blend.
Taking \(\eta_1\le1/2\) makes these ordinary profiles fully shared.

## A uniform Gaussian gap on each compact positive curvature interval

Choose an integer \(J\ge4\) larger than the eventual terminal cutoff
for the comparisons under consideration. For the actual pair use
the [effective absolute ordinary angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
For the independent blends use the
[absolute fully shared blend theorem](MULTIVARIATE_NEAR_CRITICAL_SHARED_COMMON_SYMBOL_ANGLE_SEPARATION.md).
The transfer below needs only the finite range \(3\le N<J\).
The first order has its normalized root at \(-1\). The normalized
second polynomial, also for every blend, is \(1+v z+z^2\)
with \(0<v<2\), by strict convexity of both endpoint-normalized
symbols. Thus their arguments are already separated.

Take the compact interval \(s_0\le S\le27\), where
\(0<s_0\le1\) is rational, and put \(t=S/2\).
Compare the sharp Gaussian references
\[
 G_j(z;t)=\sum_{k=0}^j\binom jk e^{-t k(1-k/j)}z^k ,
                  \qquad j=n,N.
\]
The [all-coupling phase reserve](MULTIVARIATE_GAUSSIAN_SUPERCRITICAL_PHASE_RESERVE.md)
excludes a collision everywhere on this interval for all positive
degrees. The [polynomial physical gap](MULTIVARIATE_GAUSSIAN_POLYNOMIAL_PHYSICAL_GAPS.md)
gives circular distance greater than \(t/[36nN(1+t)]\).
Since \(t/(1+t)\ge s_0/3\), and the chord is at least half
the circular distance, all the reference gaps in the finite
range \(3\le N<J\) have the common chord lower bound
\[
                               g=\frac{s_0}{216J^2}.
\]
This bound uses the absolute eventual cutoff without requiring
resultants or checks at individual Gaussian degrees.

## A chord-error bound uniform in the row list

Assume \(\eta\le1/2\). Write \(d_i=b_i-a_i\),
\(M_i=\max(a_i,b_i)\), and \(x_i(s)=a_i+d_i s\).
Then \(x_i(s)\in[(1-\eta)M_i,M_i]\) for \(s\in[0,1]\).
The endpoint-normalized symbol satisfies
\[
                  L_j''(s)=\sum_i d_i^2 j\psi'(1+jx_i(s)).
\]
Integral comparison for the trigamma series gives
\(0\le1/x-j\psi'(1+jx)\le1/(jx^2)\) for \(x>0\).
Also \(1/M_i\le\eta\), since the unequal integer exponents
have \(|d_i|\ge1\). Consequently at every \(j\ge2\)
\[
     \|L_j''-S\|_{[0,1]}
       \le\frac{S\eta}{1-\eta}
         +\frac{S\eta}{j(1-\eta)^2}\le4S\eta.
\]
The same bound holds for every blend of \(L_n,L_N\).
Subtract the quadratic of curvature \(S\) with zero endpoints.
The Green kernel for the interval with zero endpoint values
has integral \(s(1-s)/2\). Therefore
\[
              |j[\widehat L_j(s)+Ss(1-s)/2]|\le j S\eta/2.
                                                               \tag{1}
\]
At the fractional grid \(s=k/j\) the quantity in (1) is the
exact log ratio of the blended coefficient to the corresponding
Gaussian coefficient. Both exponents are nonpositive by
convexity and endpoint normalization, so the exponential
has Lipschitz constant one in this comparison. For either
an exact or a blended factorial polynomial,
\[
                    |P_j(z)-G_j(z;q)|
                        \le (jS\eta/2)(1+|z|)^j .             \tag{2}
\]

## Uniform transfer, including independent blends

If \(z=r\omega\) is a root of \(P_j\), then \(r>0\) and
\(|\omega|=1\). For any unit-circle root \(\zeta\) of the
reference,
\[
                     |r\omega-\zeta|/(1+r)
                                 \ge |\omega-\zeta|/2,
\]
as follows by squaring and using \(4r\le(1+r)^2\).
Factor the monic \(G_j\) and use (2). Some reference root
is within \(2(jS\eta/2)^{1/j}\) of \(\omega\).
On \(K=[s_0,27]\), choose
\[
             \eta\le
       \min\!\left\{\frac12,\frac{2}{27J}(g/8)^{J-1}\right\}.
\]
Then for \(j<J\) every such displacement is at most \(g/4\).
Two actual roots of equal argument would put their reference
neighbors within \(g/2\), contradicting the Gaussian gap.

Choose \(s_0>0\) rational in the
[ordinary weak-curvature range](MULTIVARIATE_UNIFORM_WEAK_CURVATURE_ANGLE_SEPARATION.md)
for the actual pair, and sufficiently small for the
[weak independent-blend theorem](MULTIVARIATE_WEAK_ASYMMETRIC_COMMON_SYMBOL_ANGLE_SEPARATION.md)
in the blended case. The latter is given in terms of
\(W=\sum_i d_i^2/((a_i+b_i)/2)\), where \(S\le W\le2S\).
For \(S\ge27\) apply the
[absolute large-curvature theorem](MULTIVARIATE_ABSOLUTE_CURVATURE_27.md)
or its
[independent-blend extension](MULTIVARIATE_ABSOLUTE_CURVATURE_27_COMMON_SYMBOL.md),
respectively. Together these cover the full curvature range.
The high-degree thresholds handle all \(N\ge J\) in
the remaining compact interval. This proves both statements.

In particular, any ordinary counterexample must have
\(S<27\) and contain a row with \(|d_i|/M_i>\eta_0\).
For that row,
\[
                    |d_i|<27/\eta_0,\qquad
                    M_i<27/\eta_0^2 ,
\]
since its strength \(d_i^2/M_i\) is at most \(S\).
Thus every possible counterexample has a retained row in one
finite absolute set of exponent pairs. The conclusion leaves
the interaction of these rows with the remaining profile open.
