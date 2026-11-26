# Exact soft-tail exclusion at transverse Gaussian contacts

Fix an ordinary two-monomial core with raw moment polynomials
\(Q_j\), and take adjacent orders \(n\ge2\), \(N=n+1\).
For strength \(t\ge0\) let
\[
 F_j(t,z)=\sum_{k=0}^j [z^k]Q_j(z)\,
                e^{-t k(j-k)/(2j)}z^k,\qquad \Lambda=N/n .
\]
Fix a positive strength \(\tau_*\) and integer \(h\). A reference
contact at zero limiting signed bias is a pair of roots
\[
      F_n(\tau_*,\zeta_n)=F_N(\tau_*,\zeta_N)=0,\qquad
                             \zeta_N=\Lambda^h\zeta_n .
\]
For roots simple in their respective polynomials, write
\[
  \alpha_j=\left.\partial_t\log\zeta_j(t)\right|_{t=\tau_*}
      =\left.\left(\frac12-
          \frac{\mathcal D^2F_j}{2j\mathcal D F_j}\right)
                               \right|_{(\tau_*,\zeta_j)},
      \qquad \mathcal D=z\partial_z .
\]

For a tail of shared unequal rows, use the signed gaps \(d_i\)
and midpoints \(M_i\), and put
\[
 \tau=\sum_i\frac{d_i^2}{M_i},\qquad
 \kappa=\frac12\sum_i\frac{d_i^2}{M_i^2},\qquad
 \varepsilon=\max_i\frac{d_i^2}{M_i}.
\]
A reflected pair is counted as two rows in these sums.

**Theorem.** Suppose all the reference contacts at \(\tau_*\)
are simple in both polynomials. Each of the following tests
gives a neighborhood of \(\tau_*\) and a positive tolerance
for \(\varepsilon\) in which the indicated ordinary tails
have no physical common root at the pair of orders.

1. The tail consists of reflected pairs, \(h=0\), and at every
   reference contact \(\operatorname{Im}(\overline{\alpha_N}\alpha_n)
   \ne0\).
2. The tail rows have a common upper exponent, signed total
   gap \(h\), and
   \[
        \left|\sum_i d_i|d_i|\right|+\left|\sum_i d_i^3\right|
            \le\varepsilon\sum_i d_i^2 .
   \]
   At every contact require
   \[
   \operatorname{Im}\left[
     \overline{(\alpha_N-\alpha_n)}
        \left(\frac{\alpha_n}{n}-\frac{\alpha_N}{N}
                                      +\frac{h}{nN\tau_*}\right)\right]
      \ne0 .
   \]
   For \(h\ne0\), it suffices that
   \(\operatorname{Im}(\alpha_N-\alpha_n)\ne0\)
   at every contact.

The thresholds can depend on the fixed core, orders, contact
strength and \(h\). The number of tail rows, their gaps, and
their exponents are unrestricted subject to the given inequalities.
The proof also permits an independent small tolerance in the
signed power bound in the second test.

## Ordinary tail logarithms through the next order

Let \(L_j(s)\) be the endpoint-normalized tail log coefficient
per degree, with the binomial coefficient stripped, and put
\(s=k/j\) on the sampling grid. It is the sum of the normalized
log-Gamma chords. It vanishes at \(s=0,1\) and has exact curvature
\[
 L_j''(s)=\sum_i j d_i^2
             \psi_1\bigl(1+j(M_i+d_i(s-1/2))\bigr).
\]
Define
\[
 S_1=\sum_i\frac{d_i}{M_i},\quad
 S_2=\sum_i\frac{d_i}{M_i^2},\quad
 S_3=\sum_i\frac{d_i^3}{M_i^2},\qquad
 V=\sum_i\frac{d_i^4+d_i^2}{M_i^3}.
\]
For \(\varepsilon\le1/2\), one has \(1/M_i\le\varepsilon\)
and \(|d_i|/M_i\le\varepsilon\), by the integer gaps.
The trigamma integral, or its Taylor expansion as
in the [signed radial response](MULTIVARIATE_ALGEBRAIC_PROFILE_TRAJECTORIES.md#signed-radial-response-in-the-vanishing-tail),
gives
\[
                 L_j''(s)=\tau-\frac{\kappa}{j}
                                 -(s-1/2)S_3+O_J(V)
                            \qquad(0\le s\le1).
\]
Indeed each term uses \(1/(M_i+x)-1/[2j(M_i+x)^2]\),
where \(x=d_i(s-1/2)\), with error \(O(M_i^{-3})\) before
multiplication by \(d_i^2\). The linear term of the reciprocal
is \(-x/M_i^2\); its remainder and the quadratic-denominator
remainder are \(O((d_i^2+1)/M_i^3)\).
The constants in the estimates are independent of the row count.
In particular \(V\le4\varepsilon\kappa\).
Comparing endpoint-zero chords gives an error
\(O_J(V+|S_3|)\) in each normalized log coefficient, relative
to the Gaussian clock \(\tau-\kappa/j\).

Write \(D_j\) for the tail endpoint argument scale, so that
the two normalized variables at a common physical root have
ratio \(D_N/D_n\). At order \(j\) a row contributes
\(\int_{a_i}^{b_i}\psi(1+jx)\,dx\) to \(\log D_j\).
To compare orders use
\[
 x\psi_1(1+jx)
      =j^{-1}-(2j^2x)^{-1}+(6j^3x^2)^{-1}+O_J(x^{-3}).
\]
For completeness, the integral for \(\psi_1(1+t)\)
is the Laplace transform of \(u/(e^u-1)\), which differs
from \(1-u/2+u^2/12\) by \(O(u^3)\) on \(u\ge0\).
Integrating first in \(j\) and then along the exponent interval
gives
\[
 \log(D_N/D_n)=h\log\Lambda-\frac{S_1}{2nN}
                   +\frac{n^{-2}-N^{-2}}{12}S_2+O_J(V).
\]
Here the signed integrals of \(x^{-1}\) and \(x^{-2}\)
are \(d_i/M_i+O(|d_i|^3/M_i^3)\) and
\(d_i/M_i^2+O(|d_i|^3/M_i^4)\) by reflection about \(M_i\).
The errors and the last integral are bounded by \(O(V)\).

For a reflected tail, \(S_3=0\) and the endpoint scales are
both exactly one. For a common upper exponent \(U\) in the
second test, \(M_i=U-|d_i|/2\). On a fixed positive compact
strength interval,
\[
 \frac1U=\frac{2\kappa}{\tau}+O(\varepsilon\kappa).
\]
This follows by writing \(2\kappa/\tau\) as the weighted mean
of \(1/M_i\), with weights \(d_i^2/M_i\). Expanding each
midpoint reciprocal gives
\[
 \begin{split}
 S_1&=\frac hU+\frac{\sum_i d_i|d_i|}{2U^2}
                           +O(\varepsilon\kappa),\\
 S_2&=O_{h,\tau_*}(\varepsilon\kappa),\\
 S_3&=\frac{\sum_i d_i^3}{U^2}
                           +O(\varepsilon\kappa).
 \end{split}
\]
Indeed the remainders after either signed numerator in
\(S_1,S_3\) are bounded by a constant times
\(U^{-3}\sum_i(|d_i|^3+|d_i|^4)\le O(\varepsilon\kappa)\).
For \(S_2\) the leading term is \(h/U^2\) and the
remaining bound is \(O(U^{-3}\sum_i d_i^2)\).
The signed power bound in the theorem therefore implies
\[
          \log(D_N/D_n)=h\log\Lambda-
                  \frac{h\kappa}{nN\tau}+O(\varepsilon\kappa),
       \qquad S_3=O(\varepsilon\kappa).
\]
In both cases the exact normalized moment arrays are within
\(O(\varepsilon\kappa)\) coefficient error of the two Gaussian
arrays at clocks \(\tau-\kappa/n,\tau-\kappa/N\).

## The transverse clock test

A possible common root for a sufficiently soft tail near
\(\tau_*\) must approach one of the reference contacts.
This follows by compactness of their nonzero root annuli,
coefficient convergence, and convergence of the positive
endpoint-scale ratio to \(\Lambda^h\). Work in fixed
neighborhoods of these finitely many contacts.
Near each simple root, the two Gaussian root branches are
analytic in their separate strengths. The perturbation of a
simple root under the coefficient error above is
\(O(\varepsilon\kappa)\), uniformly in a small neighborhood.
For example this follows from Rouché on a circle around the
reference root with radius a sufficiently large constant
times the coefficient error.

Choose the logarithm branches at each contact so their
difference is \(h\log\Lambda\) there. In the second test set
\[
 G(t,k)=\log\zeta_N(t-k/N)-\log\zeta_n(t-k/n)
                      -h\log\Lambda+\frac{hk}{nN t}.
\]
In the first test use the same expression with \(h=0\).
An actual common root would give
\(|G(\tau,\kappa)|=O(\varepsilon\kappa)\).
The two complex derivatives of \(G\), in the two real variables
at \((\tau_*,0)\), are precisely the vectors in the indicated
determinant. For \(h=0\) its nonvanishing is equivalent to
\(\operatorname{Im}(\overline{\alpha_N}\alpha_n)\ne0\).
An invertible real derivative and a Taylor bound give
\[
             |G(\tau,\kappa)|\ge c(|\tau-\tau_*|+\kappa)
\]
throughout a small neighborhood, for some \(c>0\).
Since \(\kappa>0\) for a nonempty tail of strength near
\(\tau_*>0\), reducing the tolerance contradicts the
coefficient-error estimate.

At a positive reference contact the
[angular eliminant](MULTIVARIATE_DIFFUSE_PHYSICAL_EXCEPTION_POINTS.md#the-finite-set-of-physical-reference-contacts)
makes \(q_*=\exp[-\tau_* /(2nN)]\) algebraic. Thus the
contact roots, the \(\alpha_j\), and their real and imaginary
parts are algebraic. The positive strength \(\tau_*\) is
transcendental by Hermite–Lindemann. The determinant in the
second test is
\[
       (1/n-1/N)\operatorname{Im}(\overline{\alpha_N}\alpha_n)
       +\frac{h}{nN\tau_*}
                    \operatorname{Im}\bigl(\overline{\alpha_N-\alpha_n}\bigr).
\]
If the latter numerator is nonzero, the two terms cannot
cancel. This proves the exact sufficient criterion for \(h\ne0\).

For a fixed core and \(h\), only finitely many positive reference
strengths on a compact interval can give a physical contact
at the finitely many terminals below the
[absolute angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
The contact and derivative tests above can be performed at
algebraic \(q_*\) using exact root isolation.
If each contact on the interval passes the applicable test,
ordinary weak tails in the stated class satisfy SFC(2) at every
degree there. Compactness covers the strengths away from
contact neighborhoods. The comparison explicitly depends on
the signed shape and transverse derivative hypotheses.

## Diffuse reflected tails with large individual strengths

The reflected case admits an exact test using the maximum relative
gap \(\eta=\max_i |d_i|/M_i\) in place of \(\varepsilon\). Set
\[
 Z=\sum_i\frac{d_i^4}{M_i^3},\qquad
 H(s)=\frac{(s-1/2)^4-1/16}{12},\qquad
 F_j(t,u,z)=\sum_k[z^k]F_j(t,z)\,e^{ju H(k/j)}z^k .
\]
At a simple reference contact, let \(\gamma_j\) be the derivative
of the log root of \(F_j(\tau_*,u,z)\) at \(u=0\). Implicit
differentiation gives \(\gamma_j=-\partial_uF_j/\mathcal D F_j\)
at the reference root. Write
\[
 A=\alpha_N-\alpha_n,\qquad C=\gamma_N-\gamma_n,\qquad
 B=\frac{\alpha_n}{n}-\frac{\alpha_N}{N}+\frac{h}{nN\tau_*},
 \qquad I(r)=\operatorname{Im}\bigl(\overline A(B+rC)\bigr).
\]

**Theorem.** Suppose the reference contacts at \(\tau_*\) are simple.

1. For reflected tails take \(h=0\). If \(I(0)\) and
   \(I(\tau_*)\) have the same strict sign at each reference contact,
   then sufficiently small \(\eta\) excludes exact common roots
   for strengths in a neighborhood of \(\tau_*\).
2. For tails with common upper exponent \(U\) and signed total gap
   \(h\), put
   \[
       X=\frac{\sum_i d_i|d_i|}{2U^2},\qquad
       Y=\sum_i\frac{d_i^3}{M_i^2}.
   \]
   If \(I(0)\) and \(I(2\tau_*)\) have the same strict sign at each
   reference contact, then sufficiently small \(\eta\) and
   \(\delta\) exclude exact common roots for strengths in a
   neighborhood of \(\tau_*\), provided
   \(|X|+|Y|\le\delta\kappa\).

All the bounds are uniform in the number of rows and their individual
strengths. The allowable ratio \(Z/\kappa\) can also be restricted
to any closed subinterval on which the same strict sign holds,
uniformly up to a small enlargement of that interval.

To prove the estimates, use \(M_i^{-1}\le\eta\) by integrality
and write \(g=s-1/2\). Expanding the reciprocal in the trigamma
formula through its quadratic term gives
\[
  L_j''(s)=\tau-\kappa/j-gY+g^2Z+O_J(\eta\kappa).
\]
The cubic reciprocal remainder is \(O(\eta Z)\). The error from
the second trigamma term is bounded by a constant times
\(\sum_i |d_i|^3/M_i^3\le2\eta\kappa\). The remaining trigamma
error is \(O_J(\sum_i d_i^2/M_i^3)=O_J(\eta\kappa)\).
Here \(Z/\kappa\) is twice the weighted mean of the individual
strengths \(d_i^2/M_i\) with weights \(d_i^2/M_i^2\), so
\(0\le Z/\kappa\le2\tau\). In a reflected pair the two individual
strengths agree, which gives the stronger bound \(Z/\kappa\le\tau\);
reflection also gives \(Y=0\) exactly. The upper endpoint
\(Z/\kappa=\tau\) is attained by a single reflected pair. It can
remain at positive limiting strength while \(\eta\to0\), by
taking growing integer gaps and midpoints proportional to their squares.

Since \(H''(s)=g^2\) and \(H(0)=H(1)=0\), the exact normalized
array is within \(O_J((\eta+\delta)\kappa)\) of
\(F_j(\tau-\kappa/j,Z,z)\). For reflected tails use \(\delta=0\).
The order-dependent scale ratio is exactly one for a reflected
tail. In the common-upper case the difference of the two exact
endpoint log scales has the sharper remainder
\[
 \log(D_N/D_n)
 =h\log\Lambda-\frac{S_1}{2nN}
                 +\frac{n^{-2}-N^{-2}}{12}S_2+O_J(\eta\kappa).
\]
Indeed the three errors after the displayed signed terms in the
earlier endpoint calculation are bounded by a constant times
\(\sum_i(|d_i|^3/M_i^3+|d_i|^3/M_i^4+|d_i|/M_i^3)\),
which is \(O(\eta\kappa)\).
Since \(M_i=U-|d_i|/2\), direct reciprocal expansion gives
\[
 \frac1U=\frac{2\kappa}{\tau}+O(\eta\kappa),\qquad
 S_1=\frac hU+X+O(\eta\kappa),\qquad
 S_2=\frac h{U^2}+O(\eta\kappa)=O(\eta\kappa).
\]
For the first formula, average \(M_i^{-1}-U^{-1}\) with weights
\(d_i^2/M_i\). For the second, the error is bounded by
\(\sum_i |d_i|^3/(4U^2M_i)\le\eta\kappa/2\).
For the third, \(M_i^{-2}-U^{-2}\le|d_i|/(M_i^2U)\)
and \(U^{-1}\le\eta\). These estimates are on a fixed compact
positive strength interval and \(h\) is fixed. In particular
\[
 \log(D_N/D_n)=h\log\Lambda-\frac{h\kappa}{nN\tau}
                           -\frac{X}{2nN}+O(\eta\kappa).
\]

As before, any possible common root near \(\tau_*\) must lie
near one of the finitely many reference contacts. Expanding the
log-root difference of the two simple branches and accounting
for the endpoint ratio at that contact yields
\[
 0=A(\tau-\tau_*)+B\kappa+CZ+
       O\bigl((\eta+\delta)\kappa+(\tau-\tau_*)^2+\kappa^2\bigr).
\]
The hypothesis implies \(A\ne0\), so shrinking the neighborhood
first gives \(|\tau-\tau_*|=O(\kappa)\) by the triangle inequality.
Also \(\kappa\le\eta\tau/2\).
Taking the imaginary part after multiplying by \(\overline A\)
then gives
\[
           0=\kappa I(Z/\kappa)+O((\eta+\delta)\kappa).
\]
The affine function \(I\) has a strict margin on the prescribed
closed interval by the endpoint test, and retains that margin
on a slightly enlarged interval around \(\tau_*\).
Reducing the tolerances for \(\eta,\delta\) gives the contradiction.

The new response \(\gamma_j\) is algebraic at the isolated
algebraic \(q_*\), as are the original responses. The signs in
the theorem admit exact decisions: multiplying by the positive
\(\tau_*\) gives polynomials in the transcendental \(\tau_*\)
with real algebraic coefficients. A zero polynomial is decided
algebraically and any nonzero evaluation is separated by certified
approximations. As in the soft-tail test, verifying these
conditions on the finite reference set in a compact positive
strength interval gives an all-degree exclusion for tails in
the stated classes.
