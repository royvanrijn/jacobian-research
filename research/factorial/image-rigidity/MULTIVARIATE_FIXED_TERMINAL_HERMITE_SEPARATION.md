# A linked Hermite comparison at every fixed terminal

For ordinary distinct fully shared exponent columns, let
\(a_i\) be the smaller entry of each unequal row and \(d_i>0\)
its gap. This quantifies the
[fixed-terminal small-curvature separation](MULTIVARIATE_FIXED_TERMINAL_CURVATURE.md#angular-separation-at-small-curvature)
using a linked Cayley expansion. At lower order \(n\) define
\[
              X_n=\sum_i\sum_{h=1}^{d_i}\frac{d_i}{na_i+h}.
\]
It is the sum of the Gauss reciprocals in the lower coefficient
array. Set \(N=n+1\) and \(S=\sum_i d_i^2/(a_i+d_i)\).

**Theorem.** For each fixed \(n\ge1\) there is a computable rational
\(\lambda_*(n)>0\) such that \(X_n\le\lambda_*(n)^2\) gives
disjoint root arguments of the two factorial moments, also after
independent positive rescalings. In particular, it suffices that
\[
                    S\le\frac{n\lambda_*^2}{1+n\lambda_*^2}.
\]
For \(n\ge2\), it also suffices that the order-\(n\) outer
geometric mean is at least \((1+(n-1)\lambda_*^2)^{-1/2}\).
The result uses the actual common exponent rows to determine both
limiting Hermite polynomials.

## Exact cancellation in the Cayley transform

Normalize each factorial array at its two endpoints. The
[Gauss-factor formula](MULTIVARIATE_LINKED_GAMMA_ANGLE_COMPACTIFICATION.md#linked-factors-and-comparison-across-orders)
writes the normalized multiplier at \(k\) from an increasing
factor of reciprocal \(x\) as the exponential of
\[
      \sum_{l=0}^{k-1}\log(1+l x)
             -\frac{k}{j}\sum_{l=0}^{j-1}\log(1+l x).
\tag{1}
\]
A decreasing factor replaces \(k\) by \(j-k\). Put \(X=X_n>0\)
and write \(x_{j,i,h}=X\xi_{j,i,h}\) for both orders.
Then \(R_j=\sum_{i,h}\xi_{j,i,h}\) satisfies \(R_n=1\) and
\[
                 \frac nN\le R_N\le\frac nN+\frac{X}{N},
                 \qquad R_N\le1 .
\tag{2}
\]
Indeed \(x_{N,i,h}/x_{n,i,h}=n/N+h/[N(Na_i+h)]\) and
\(h/(Na_i+h)\le x_{N,i,h}\le X\).
In particular all \(\xi_{j,i,h}\) lie in \((0,1]\) and have
sum at most one at each order.

For the fixed normalized factor lists, replace \(X\) in (1)
by a complex parameter \(s\), and write \(F_{j,k}(s)\) for the
product of these deformed multipliers. On \(|s|\le(4N^2)^{-1}\)
the logarithms converge by their Taylor series and
\(|\log F_{j,k}(s)|\le2j(j-1)|s|\).
The coefficient of \(s^r\) in the logarithm, at either
orientation, is a polynomial in \(k\) of degree at most
\(r+1\), by summing \(l^r\) and subtracting the endpoint term.
For \(r=1\) their sum is \(R_j(k^2-jk)/2\).
Consequently the coefficient of \(s^p\) in \(F_{j,k}(s)\)
has degree at most \(2p\) and leading \(k^{2p}\)-coefficient
\((R_j/2)^p/p!\). The leading degree can only come from the
\(r=1\) terms, since \(r+1<2r\) for \(r>1\).

Set \(s=\lambda^2\) and apply the Cayley map
\(t=(-1+\lambda z/2)/(1+\lambda z/2)\).
If \(P_j\) denotes the endpoint-normalized polynomial at this
parameter value, consider its transformed polynomial
\[
       C_j(z,\lambda)=\lambda^{-j}
                   (1+\lambda z/2)^j\,P_j(t).
\]
For the binomial factor at \(k\), the Cayley multiplier before
the division is \((-1)^k(1+\lambda z/2)^{j-k}
(1-\lambda z/2)^k\).
The coefficient of \(w^l\) in \((1+w)^{j-k}(1-w)^k\) is a
polynomial in \(k\) of degree at most \(l\) with leading
coefficient \((-2)^l/l!\). This follows by expanding
\((1+w)^j\exp(k[\log(1-w)-\log(1+w)])\).
The binomial alternating sum in \(k\) vanishes on polynomials
of degree less than \(j\), and takes \(k^j\) to \((-1)^j j!\).
Thus every power \(\lambda^{2p+l}\) of degree below \(j\)
cancels before the division. At degree \(j\) the only surviving
contributions are the highest-degree terms in \(k\). We obtain
\[
 C_j(z,0)=\mathcal H_j(z;R_j)
     =j!\sum_{p=0}^{\lfloor j/2\rfloor}
                 \frac{(R_j/2)^p z^{j-2p}}{p!(j-2p)!}.
\tag{3}
\]
This is an analytic continuation in \(\lambda\) for the
deformed family. The actual pair is at \(\lambda=\sqrt X\).

## Quantitative comparison with the two Hermite anchors

Put \(\lambda_0=1/(2N)\). On \(|\lambda|=\lambda_0\),
\(|F_{j,k}(\lambda^2)|<2\) by (1) and the preceding
logarithm bound. The modulus of the coefficient of \(z^l\)
in \(C_j(z,\lambda)\) is at most
\[
          B_{j,l}=2^{j-l+1}\binom jl\lambda_0^{l-j}.
\]
Indeed the coefficient of \(w^l\) in the Cayley factor
has modulus at most \(\binom jl\) at integer \(k\).
Cauchy's formula at \(\lambda_0\) therefore bounds the
difference of the coefficient from its value at zero by
\(2B_{j,l}\lambda/\lambda_0\) for \(0\le\lambda\le\lambda_0/2\).
Use the two fixed anchors
\(P_n^0(z)=\mathcal H_n(z;1)\),
\(P_N^0(z)=\mathcal H_N(z;n/N)\).
By (2), the corresponding difference of the parameter
\(R_N\) in (3) costs at most \(\lambda^2/N\) times the
coefficient bounds for \(\partial_R\mathcal H_N\) on
\(n/N\le R\le1\). Both together give a computable rational
\(C_R\) such that \(|C_j(z,\lambda)-P_j^0(z)|\le C_R\lambda\)
for \(|z|\le R\), every profile, and \(0<\lambda\le\lambda_0/2\).

For the anchors, \(i^{-j}\mathcal H_j(ix;\alpha)\) is the
monic probabilists' Hermite polynomial with variance \(\alpha\).
It is the determinant of the symmetric tridiagonal Jacobi
matrix with off-diagonal entries
\(\sqrt{\alpha},\sqrt{2\alpha},\ldots,\sqrt{(j-1)\alpha}\),
subtracted from \(xI\). This follows from the three-term
determinant recurrence. Its eigenvalues are real and simple
because an eigenvector is determined by its first entry.
Thus all anchor roots are simple and purely imaginary.

The anchors have no common root, even after the variance
scaling specified above. To verify this for every order,
write their real polynomials as \(x^{j\bmod 2} B_j(x^2)\).
At even order \(j=2m\) the squared polynomial is irreducible
over the rationals. For \(m=1\) this is immediate. For
\(m>1\), its constant-normalized expression in \(y=x^2/\alpha\)
has coefficients \((-1)^k\binom mk/u_{2k}\), where \(u_{2k}\)
is the product of the odd integers at most \(2k\).
Apply the even Schur irreducibility result to the expression
in \(x^2\), with constant and leading numerator of modulus
one. Its irreducibility as a polynomial in \(x\) implies
irreducibility as a polynomial in \(x^2\); a rational nonzero
scaling of \(y\) preserves the latter. [Allen and Filaseta, Theorem 1](https://www.impan.pl/shop/publication/transaction/download/product/83056)

Thus when \(n\) is odd the upper even squared polynomial
has larger degree and cannot share a root with the lower
squared polynomial. When \(n\) is even the lower squared
polynomial is irreducible and both squared polynomials are
monic of the same degree. Their subleading coefficients
are \(-n(n-1)/2\) and \(-n^2/2\), so they are distinct.
Only the odd anchor has a zero at the origin.
This uses irreducibility of the even squared polynomial only.

Isolate all imaginary heights of these two rational anchors
by real polynomial root isolation. Their pairwise separation
has a positive rational lower bound \(g_0\).
Choose a rational \(d\le\min(1/4,g_0/8)\) and a radius \(R\)
covering all disks of radius \(d\) about the anchor roots.
On the boundary of each disk about an anchor of order \(j\)
the modulus of that anchor is at least
\(d(g_0-d)^{j-1}\). The checker can improve this bound by
multiplying the individual separation lower bounds minus
\(d\). Decrease \(\lambda_*\le\lambda_0/2\) so that
\(C_R\lambda_*\) is smaller than all these rational margins.
Rouché then places exactly one root of \(C_j\) in each
corresponding disk. Since there are \(j\) disks, all its
roots are accounted for.

## Transfer to root arguments and a reproducible bound

For \(\lambda R\le1\), the signed clockwise angle from the negative real axis
under the Cayley map is
\[
 \operatorname{Im}\log\frac{1+\lambda z/2}{1-\lambda z/2}.
\]
Dividing by \(\lambda\) differs from \(\operatorname{Im} z\)
by at most \(\lambda^2R^3/3\), using the convergent odd-power
logarithm series. Require \(\lambda_*^2R^3/3<d\).
Each exact root argument then belongs to the angular interval
whose divided height is within \(2d\) of its anchor height.
The intervals for both orders are disjoint. The transformation
misses \(t=1\), which cannot be a root since the coefficients
are positive. Positive normalizations of the physical parameter
preserve angles, so the two original polynomials are coprime.
A rational Bézout identity transfers the result to every
characteristic-zero field.

If \(S<1\), each \(d_i/(a_i+d_i)\le S\) and
\(X_n\le\sum_i d_i^2/(na_i)\le S/[n(1-S)]\).
For the outer-mean variant,
\(u^{-2}=\prod_{i,h}(1+(n-1)x_{n,i,h})\ge1+(n-1)X_n\).
These give the stated alternative sufficient conditions.

The [rational checker](verify_multivariate_fixed_terminal_hermite.py)
computes \(\lambda_*\) from exact Hermite coefficients, a rational
gcd, real root isolation, and the explicit bounds above. It can
be run with any lower order \(n\). Its certificate gives a
small, uniform near-flat region independently of any precomputed
coefficient envelopes for that degree.

At terminal nine the replay certifies \(\lambda_*(8)=2^{-52}\)
from 17 disjoint anchor intervals. Thus \(X_8\le2^{-104}\)
is one exact sufficient near-flat bound. This radius is
conservative because it uses a Cauchy majorant for each
transformed coefficient.
