# Uniform weak same-grid transfer for mirrored factorial symbols

Take a fully shared mirrored profile with positive smaller
exponents \(\ell_i\) and integer gaps \(d_i\ge1\). Count both
orientations of each mirrored pair. Write
\[
      a_i=\ell_i/d_i+\tfrac12,\qquad
      W=\sum_i d_i/a_i,\qquad L=\min_i\ell_i/d_i .
\]
At each order \(j\), let \(L_j(s)\) be the log-factorial coefficient
symbol per degree, with the binomial coefficient removed and
the affine endpoint interpolation subtracted. For \(n\ge2\)
let \(N=n+1\).

**Theorem.** If \(W\le10^{-7}\), the polynomial
\[
      \sum_{k=0}^n\binom nk e^{nL_N(k/n)}z^k
\]
has simple unit-circle roots. They pair with the order-\(n\)
factorial roots in endpoint normalization with angular
displacement at most
\[
                    13000 W^{3/2}\frac{1+\log n}{n^2}.        \tag{1}
\]
The threshold and bound are uniform in the degree, the
number of rows and the smaller exponents. The result
compares two factorial symbols on the same order-\(n\)
grid.

## A mesh bound and the shifted product

The [mirrored circle-root factorization](MULTIVARIATE_MIRRORED_ROOT_CIRCLE.md)
applies to positive real \(\ell_i\) as well as integers, since
its Gauss factorization only uses the integrality of \(d_i\).
Use the Gegenbauer factor from a row attaining \(L\), with
parameter \(\alpha=(n\ell_i+1)/d_i\). Then
\(nL\le\alpha\le nL+1\). For \(\alpha\ge1\), its Dirichlet
equation for \(v(t)=(\sin t)^\alpha C_n^\alpha(\cos t)\)
has coefficient at most
\[
           (n+\alpha)^2-\alpha(\alpha-1)
                           =n^2+2n\alpha+\alpha .
\]
The Poincaré inequality between consecutive zeros, including
the left endpoint and the first zero, gives gap at least
\(\pi/\sqrt{n^2+2n\alpha+\alpha}\). At zero the boundary term
is \(v v'=O(t^{2\alpha-1})\) and the potential times \(v^2\)
is integrable because \(\alpha\ge1\). In the circular variable,
reflection makes the wrap-around gap four times the first
positive \(t\) root, so the same circular lower bound applies.
For \(n\ge2\) and \(L\ge1600\) the denominator is at most
\(n\sqrt{1+3L}\). The transformation \(z=e^{2it}\) and
the [Schur–Szegő mesh comparison](MULTIVARIATE_MATCHED_LAPLACE_NONLINEAR_TRANSFER.md#actual-mirrored-reference-profiles)
therefore give the actual order-\(n\) polynomial angular mesh
at least \(\gamma/n\) with
\[
                    \gamma=2\pi/\sqrt{1+3L}.
\]
Here \(W\ge2/(L+1/2)\) since the row attaining \(L\)
has its mirrored partner and \(d_i\ge1\). The assumed weak
cap in particular gives \(L\ge1600\).

Fix an actual reference root and use the gauged function
\(\widetilde F\) and real phase coefficient \(b\) from the
[full sine-mesh estimate](MULTIVARIATE_MATCHED_LAPLACE_MESH_SHARPENING.md#the-exponential-rate-from-the-full-circular-mesh).
Use the following constants:
\[
 \begin{gathered}
 M=\pi/\gamma,\qquad D=2+4M^2,\qquad a_0=1/(32M),\\
 K_0=M+\frac1{2n}+\frac{a_0D}{2},\qquad
 K_1=K_0+\frac1{2a_0},\\
 B=\exp(4a_0^2D+a_0),\qquad
 b_0=\tfrac12\exp(-4a_0^2D).
 \end{gathered}
\]
Then \(|b|\le M(1+\log n)\), and the shifted-product bounds are
\[
       |\widetilde F(h+y)|\le B(|h|+2a_0)e^{K_0|h|},
       \qquad |\widetilde F(y)|\ge b_0|y|
                      \qquad (h\in\mathbb R,\ |y|\le2a_0).
\]
For \(|y|\le a_0\) Cauchy's estimate at radius \(a_0\) bounds
\(|\widetilde F'(y)|\) by \(2B\). Its second-derivative bound at
\(y+h\) is \(2B(|h|+2a_0)e^{K_0|h|}/a_0^2\).
Taylor's formula along the real segment therefore gives
\[
 |\widetilde F(y+h)-\widetilde F(y)-h\widetilde F'(y)|
                  \le \frac{2B}{a_0}h^2e^{K_1|h|}.
                                                               \tag{2}
\]
Use \(x+2a_0\le2a_0 e^{x/(2a_0)}\) for \(x\ge0\)
in the second-derivative remainder.

## Symmetric jumps and the uniform bound

The [positive rowwise Laplace increment](MULTIVARIATE_MATCHED_LAPLACE_NONLINEAR_TRANSFER.md#positive-curvature-increment-for-ordinary-shared-rows)
for \(L_N-L_n\) gives measure \(\nu(dh)=h^{-2}\mu(dh)\).
It is symmetric because each increasing row has its reflected
decreasing row. Let
\[
                     T=n\int h^2 e^{K_1|h|}\,\nu(dh).
\]
First truncate \(\nu\) to \(|h|\ge\epsilon>0\).
Let \(H\) be the compound Poisson sum with jump intensity
\(n\nu\) on this truncation and let
\(\chi=\mathbb E e^{-ibH}\). Its logarithmic moment generating
function is \(n\int(\cosh(th)-1)\,d\nu\) by symmetry.
Differentiating the characteristic function and using
\(1-\cos x\le x^2/2\) gives
\[
         |\chi|\ge e^{-b^2T/2},\qquad
         |\mathbb E[H e^{-ibH}]/\chi|\le |b|T .
\]
At \(t=K_1\) the logarithmic moment and its first and second
derivatives are bounded in magnitude by
\(K_1^2 T/2,\ K_1T,\ T\), respectively.
Thus
\[
             \mathbb E[H^2 e^{K_1|H|}]
                      \le 2e^{K_1^2T/2}(T+K_1^2T^2).
\]
The factor two uses the sum of the exponential moments at
\(K_1\) and \(-K_1\).

After dividing the phase-adjusted shifted polynomial by \(\chi\),
the linear term from \(\widetilde F'(y)\) is bounded separately.
By (2) the difference from \(\widetilde F(y)\), throughout
\(|y|\le a_0\), is at most
\[
 E=2B|b|T+\frac{4B}{a_0}e^{(K_1^2+b^2)T/2}(T+K_1^2T^2).
                                                               \tag{3}
\]
All estimates are uniform in the truncation. The limit is the
exact multiplier \(e^{n(L_N-L_n)}\), up to the affine term in the
[positive jump representation](MULTIVARIATE_MATCHED_LAPLACE_NONLINEAR_TRANSFER.md#positive-curvature-increment-for-ordinary-shared-rows).
Here its linear part vanishes by reflection.

For the constants above, \(1\le M\le\sqrt L\), \(B\le2\),
\(b_0\ge1/3\), and \(K_1\le18\sqrt L\le L/2\).
For example \(a_0^2D\le6/1024\) and
\(K_1=(273/16)M+1/(32M)+1/(2n)\).
The exact row bound \(\mu_i(du)\le u e^{-a_i u}\,du/(2nN)\)
and \(a_i-K_1\ge a_i/2\) imply
\[
                    T\le\frac{2}{N}\sum_i a_i^{-2}
                        \le\frac{2W}{Ln}.
\]
Consequently \(K_1^2T\le648W/n\) and
\(b^2T\le2W(1+\log n)^2/n\).
For \(n\ge2\), the maximum of \((1+\log n)^2/n\) is \(4/e<2\),
so the exponent
in (3) is at most \(164W\).
In particular \(e^{(K_1^2+b^2)T/2}\le2\) and
\(T+K_1^2T^2\le2T\). Using \(1/a_0=32M\), (3) becomes
\[
                     E\le2056(1+\log n)\frac{W}{n\sqrt L}.
\]
If \(E>0\), the radius \(\rho=2E/b_0\) is less than \(a_0\):
indeed \(\rho/a_0\le394752(1+\log n)W/n\le394752W\).
On that circle the strict Rouché comparison holds by the lower
bound \(b_0\rho\). The uniform strict inequality passes to the
truncation limit. The disjoint log disks each contain one root.
Reciprocal conjugation preserves each disk and the transformed
polynomial, so each root is on the unit circle. When \(E=0\)
the compared arrays agree. Finally
\(\rho/n\le12336 W(1+\log n)/(n^2\sqrt L)\) and \(L\ge1/W\),
which give (1).
