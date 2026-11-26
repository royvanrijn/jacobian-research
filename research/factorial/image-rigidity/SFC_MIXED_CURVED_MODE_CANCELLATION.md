# Two consecutive cancellations in a repeated-degree saddle block

Consider the exponents of \(X^D,Y^D,X^{2D}\), with integer \(D\)
fixed and sufficiently large. As the starting order \(m\) grows,
coefficients can be chosen so that the two degree-\(D\) saddle
contributions cancel to leading order at two consecutive shifts.
Their scaled limiting sum is a nonzero multiple of
\(\binom{\ell}{2}\), where \(\ell\) is the shift. This is an
obstruction to assigning only one cancellation to a two-mode
block when a foreign monomial survives at a repeated saddle.
The conclusion concerns a limiting profile and gives no exact
ordinary three-window collision.

Use the normalized peaks of the
[separated-support theorem](SFC_DILATED_SUPPORT_SADDLE_ENVELOPE.md)
at \(t=m+\ell\) and take raw weights at \(m\) equal to
\(1,q_2,q_3\) respectively. Peel the common order-dependent
weight \((t/m)^D\) in the degree-\(D\) block. Write
\(P(x)=xe^{1-x}\), \(r_D=(e/4)^D\) and \(\lambda=r_Dq_3\).
On the \(X\) box, the \(X^{2D}\) term relative to \(X^D\)
is exactly \(\lambda(t/m)^D x^D\); on the \(Y\) box both
\(X\) terms vanish to order at least \(D\) at the inactive
endpoint. For a fixed small \(\delta>0\) define
\[
 R_t(\lambda)=
 \frac{\int_{1-\delta}^{1+\delta} P(x)^{Dt}
                          (1+\lambda (t/m)^D x^D)^t\,dx}
      {\int_{1-\delta}^{1+\delta}P(x)^{Dt}\,dx}.
\]
The inactive endpoint factor is common. The sum of the
two reduced local terms, after division by the pure
local integral, is \(R_t+q_2^t\). The omitted foreign terms
on the two corresponding local boxes lose \(O(t^{1-D})\)
relative to the respective saddle sizes by the inactive-endpoint
estimate of the
[incomparable-coordinate theorem](SFC_MIXED_DEGREE_INCOMPARABLE_COORDINATES.md#two-local-estimates).

For small complex \(B\), let \(g(B)\) be the saddle value of
\(D\log P(x)+\log(1+Bx^D)\) near \(x=1\). This is analytic,
\(g'(B)=1+O(D|B|)\) on \(|B|\le 2r_D\), and its second
derivative is \(O(D)\) there. The critical point solves
\(x=1+Bx^D/(1+Bx^D)\), which is a contraction near \(1\)
for \(|B|\le 2r_D\); differentiating the critical value
gives \(g'(B)=x(B)^D/(1+Bx(B)^D)\).
On a fixed complex neighborhood of \(\zeta=1\) for fixed
\(D\), the translated saddle estimate yields
\[
   \log R_{m\zeta}(\lambda)
                 =m\zeta g(\lambda\zeta^D)+B_m(\lambda,\zeta),
\]
where \(B_m\) is uniformly bounded and holomorphic locally
for \(|\lambda|<r_D\). Indeed the numerator uses the
Gaussian saddle with analytic nonzero leading amplitude
near its pure value, and the denominator uses the pure
Gaussian saddle. The peeled errors are \(O(m^{-1/2})\)
uniform in a small complex neighborhood, so the ratio
has a bounded logarithm after the exact exponential is
removed. Cauchy bounds give the same boundedness for the
first three \(\zeta\)-derivatives on a smaller neighborhood.

The analytic functions
\[
 \Lambda_m(\lambda)=\log R_{m+1}(\lambda)
                         -(1+1/m)\log R_m(\lambda)
\]
converge locally uniformly to \(D\lambda g'(\lambda)\).
There is a simple solution
\(\lambda_*=i r_D/2+O(D r_D^2)\) of
\(D\lambda_*g'(\lambda_*)=i D r_D/2\), since
\(\lambda g'(\lambda)=\lambda+O(D\lambda^2)\) with
derivative \(1+O(D|\lambda|)\). Choose odd integers
\(2k_m+1\) such that \((2k_m+1)\pi/m\to D r_D/2\).
Rouché's theorem near \(\lambda_*\) gives
\(\lambda_m\to\lambda_*\) with
\(m\Lambda_m(\lambda_m)=i(2k_m+1)\pi\).
Set \(q_3=\lambda_m/r_D\) and
\(q_2=R_{m+1}(\lambda_m)/R_m(\lambda_m)\).
Then \(q_2^m=-R_m\) and \(q_2^{m+1}=-R_{m+1}\) exactly.

Put \(F(\zeta)=\zeta g(\lambda_*\zeta^D)\).
For each fixed \(\ell\), Taylor expansion at \(\zeta=1\),
using the bounded derivatives of \(B_m\), shows that
\[
  m\,\frac{R_{m+\ell}+q_2^{m+\ell}}{R_m q_2^\ell}
                       \longrightarrow F''(1)\binom{\ell}{2}.
\]
Here
\(F''(1)=D(D+1)\lambda_*g'(\lambda_*)+
D^2\lambda_*^2g''(\lambda_*)\ne0\) for large \(D\), since the
second term is smaller by \(O(D r_D)\).
Because \(|q_3|\to|\lambda_*/r_D|=1/2+O(D r_D)\),
the third peak and the exterior integrals are exponentially
small relative to the displayed normalization, by the peak
separation estimates. The two inactive-endpoint errors
are \(o(1/m)\) for \(D>2\). Consequently the same quadratic
leading limit holds for the ordinary moments after the
common scaling and division by the pure local integral.

