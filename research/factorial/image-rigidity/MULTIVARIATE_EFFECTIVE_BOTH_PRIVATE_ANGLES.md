# Effective angular separation with private variables on both sides

Suppose each of two distinct ordinary monomials has a private
variable, meaning its exponent in the other monomial is zero.

**Theorem.** There is a computable absolute terminal cutoff above which
the consecutive factorial moment polynomials have disjoint root
arguments, including under independent positive rescalings.
The bound is uniform over the two exponent columns and their dimension.

If the total private degree on each side is at least two, the
[private-degree angle bounds](MULTIVARIATE_TWO_MONOMIAL_NONDIVISIBILITY.md#all-orders-when-both-private-exponent-sums-are-at-least-two)
use only the root arguments at lower orders \(n\ge3\).
If the two private degrees are one and at least two, respectively,
the [mixed private-degree phase proof](MULTIVARIATE_ONE_PRIVATE_LINEAR_TAIL.md#4-consecutive-separation)
uses only the root arguments at lower orders \(n\ge360\).
Thus those estimates work for the angular claim with independently
varying normalized radii. We treat the remaining case of one private
linear row on each side.

Use the notation of the
[effective two-private linear proof](MULTIVARIATE_EFFECTIVE_BOTH_PRIVATE_LINEAR.md).
In particular \(W=\sum_i d_i^2/m_i\) is the shared curvature and
\(P_j(z)=\sum_{k=0}^j c(j,k)z^k\) is the endpoint-normalized moment.
If \(W=0\), these are consecutive geometric sums, whose nontrivial
root arguments are distinct between adjacent orders. If \(W\ge128\),
the [high-curvature angular theorem](MULTIVARIATE_EFFECTIVE_ABSOLUTE_CURVATURE_BOUND.md)
applies. In what follows \(0<W<128\).

## Geometric halves

For each prescribed \(w_0>0\), the first geometric-tail argument
of the [effective proof](MULTIVARIATE_EFFECTIVE_BOTH_PRIVATE_LINEAR.md#the-two-geometric-tail-ranges)
already excludes equal arguments for \(W\ge w_0\). Indeed it compares
the angular logarithmic derivatives on a common disk, at the two
actual radii, and its root radius bounds are single-order estimates.
The [phase subtraction](MULTIVARIATE_BOTH_PRIVATE_UNIFORM_EVENTUAL.md#an-angular-comparison-for-split-polynomials)
only requires equality of the two arguments.

For \(W\le w_0\) and \(\mathcal L=nW\) large, use the smooth
cutoff halves \(L_j^*,U_j^*\) and the tilded halves in
[estimate (8)](MULTIVARIATE_BOTH_PRIVATE_UNIFORM_EVENTUAL.md#small-curvature-with-nw-unbounded).
The nonvanishing disks, root annulus and derivative bounds are
single-order estimates. Let \(\lambda_j=\log|z_j|\) at a proposed
root of argument \(0<\theta\le\pi\) and set
\[
 K_j(\lambda,\theta)=
      \log|L_j^*(e^{\lambda+i\theta})|
           -\log|U_j^*(e^{-\lambda-i\theta})|.
\]
The modulus equation at the root is \(j\lambda_j=K_j(\lambda_j,\theta)\).
Estimate (8) gives, throughout the comparison annulus,
\[
 |\partial_\lambda K_j|\le C/W,\qquad
             |K_N(\lambda,\theta)-K_n(\lambda,\theta)|
                       \le C/(n^2 W).
\]
For the radial derivative, the unchanged geometric factors
contribute \(O(W^{-1})\) and the tilded factors
\(O((jW^2)^{-1})\). At fixed radius the geometric factors cancel
and the real-order change of the tilded log moduli is
\(O((n^2W)^{-1})\) by (8). Subtracting the modulus equations
therefore gives
\[
       (N-C/W)|\lambda_N-\lambda_n|
                            \le|\lambda_n|+C/(n^2 W).
\]
The single-order radius bound \(n|\lambda_n|\le
C((nW)^{-1}+W)\) from the effective small-curvature proof now
yields
\[
                   |\lambda_N-\lambda_n|
                       \le C\bigl((n^3 W)^{-1}+W/n^2\bigr).
\]
Absorb \(C/W\) on enlarging the cutoff for \(\mathcal L\).
At fixed radius, (8) bounds the adjacent-order change of the
logged angular derivatives by \(C(nW)^{-2}\).
The radial move changes them by at most
\(C W^{-2}|\lambda_N-\lambda_n|\), using (8) through
two field derivatives and the geometric factors.
Thus the total angular slope difference is at most
\[
                       C(\mathcal L^{-2}+\mathcal L^{-3}
                                                +(n\mathcal L)^{-1}).
\]
This is less than one for sufficiently large computable
\(\mathcal L\) and \(n\). Phase subtraction gives the angular
exclusion in this range without a constraint on the relative
normalized radii.

## Bounded \(nW\) and escaping scaled angles

Fix the computable cap \(\mathcal L_0\) left by the geometric
estimates. For \(nW\le\mathcal L_0\), the
[bounded-curvature estimates](MULTIVARIATE_EFFECTIVE_BOTH_PRIVATE_LINEAR.md#bounded-nw)
give the annulus \(|\lambda_j|\le B/j\) at all roots,
and
\[
 (1-z)P_j(z)=1-z^{j+1}+D_j(z),\qquad
 D_j(z)=\sum_{k=1}^j [c(j,k)-c(j,k-1)]z^k.
\]
On the enlarged annulus \(|\log|z||\le2B/n\) the same first
and second coefficient differences give
\[
 |(1-z)D_j(z)|\le C/n,\qquad
 |(1-z) z D'_j(z)|\le C,\qquad
 |D_N(z)-D_n(z)|\le C/n.                                   \tag{1}
\]
For the second bound, apply one first difference to the
coefficients \(k(c(j,k)-c(j,k-1))\). Their endpoint sizes
are \(O(1)\) and the differences are \(O(1/n)\).
The annulus bounds every \(|z|^k\) by a fixed constant.
For the last bound, at common indices the two order differences
of the first coefficient differences are \(O(n^{-2})\) by
the mixed real-order estimate, while the extra final difference
is \(O(n^{-1})\). These are the coefficient comparisons in the
effective two-private proof and do not use an equality of radii.

Suppose \(n\theta\ge H\), with \(H\) a fixed sufficiently large
computable number. Throughout the horizontal radial comparison,
\(|1-z|\ge (2/\pi)|z|^{1/2}\theta\), by the chord bound for
\(\sin(\theta/2)\). Thus (1) gives \(|D_j(z)|\le C/H\).
Choose \(H\) so this is less than \(1/4\).
The analytic logarithm of \(1+D_j(z)\) near one has radial
derivative bounded by \(Cn/H\), using the second bound in (1).
The exact root moduli satisfy
\[
 (j+1)\lambda_j=\Re\log(1+D_j(z_j)).
\]
The log difference at a fixed radius is \(O(1/n)\) by (1).
Therefore
\[
       (N+1-Cn/H)|\lambda_N-\lambda_n|
                              \le|\lambda_n|+C/n.
\]
Increasing \(H\) to absorb the derivative term gives
\(|\lambda_N-\lambda_n|\le C'/n^2\), with a bound independent
of further increases in \(H\). The root arguments likewise
satisfy \((j+1)\theta=2\pi m_j+\Im\log(1+D_j(z_j))\).
Their subtraction gives
\[
               |\theta-2\pi(m_N-m_n)|\le C''/n .
\]
Here the log difference costs \(O(1/n)\) at a fixed radius
and \(O(n|\lambda_N-\lambda_n|/H)=O(1/n)\) on moving the radius.
The constant \(C''\) is uniform for all larger \(H\).
Taking \(H>C''\) excludes the angular equation, because the
distance of \(\theta\in(0,\pi]\) to the nearest multiple of
\(2\pi\) is \(\theta\).

## Bounded scaled angles

In the remaining range \(n\theta\le H\), use \(h=j\log z_j\)
at each order, with the logarithm near one. These fields are
bounded by the annulus and the angle cap. Write
\[
 h=h_n,\qquad h_N=h+h/n+v,\qquad
                  v=N(\lambda_N-\lambda_n)\in\mathbb R .
\]
The annulus bounds \(v\) on a fixed compact interval.
With \(\Lambda=nW\), the model is
\[
 F_\Lambda(h)=\int_0^1 e^{-\Lambda u(1-u)/2+hu}\,du.
\]
The two trapezoid estimates and the changes of \(\Lambda\) and
\(h/n\) imply \(|F_\Lambda(h)|+|F_\Lambda(h+v)|=O(n^{-1})\)
at a root pair. The
[model zero theorem](MULTIVARIATE_BOTH_PRIVATE_UNIFORM_EVENTUAL.md#small-curvature-with-bounded-nw)
shows that every zero of \(F_\Lambda\) is imaginary and simple.
Consequently, on the compact parameter and field caps,
\[
       |F_\Lambda(h)|
          +\left|\int_0^1\partial_h F_\Lambda(h+t v)\,dt\right|
\]
has positive infimum. At a zero of the first term, vanishing of
the divided difference would give distinct real-displaced zeros
if \(v\ne0\), or a multiple zero if \(v=0\).
This positive bound is computable by interval subdivision.
The equation \(F_\Lambda(h+v)-F_\Lambda(h)=O(n^{-1})\)
therefore yields \(v=O(n^{-1})\).

At the two orders the integral in the trapezoid expansion has
exponents \(L_n(u)+hu\) and \(L_N(u)+h_Nu\).
The real-order coefficient estimates give
\[
 L_n(u)=-\tfrac{\Lambda}{2}u(1-u)+O(n^{-1}),\qquad
 L_N(u)-L_n(u)=-\tfrac{\Lambda}{2n}u(1-u)+O(n^{-2}).
\]
For the latter identity use \(L_N-L_n=G+R_N-R_n\) and the
small \(W\) expansion of \(G\); \(H(u)\) is unchanged in this
difference in the notation of that estimate. The explicit
trapezoidal endpoint terms change by \(O(n^{-2})\), since
\(h_N-h=O(n^{-1})\). Expanding the linked integral difference
to first order and multiplying by \(n\) thus give
\[
 \left|\mathcal D F_\Lambda(h)+nv\,\partial_h F_\Lambda(h)\right|
                      \le C/n,\qquad
                  \mathcal D=h\partial_h+\Lambda\partial_\Lambda.
\]
Multiplication by \(\overline{\partial_h F_\Lambda}\) and
taking imaginary parts remove the real term \(nv\).
The compact residual
\[
   |F_\Lambda(h)|+
          \left|\Im\bigl(\mathcal D F_\Lambda(h)
                                 \overline{\partial_h F_\Lambda(h)}\bigr)\right|
\]
has positive infimum. To see this at a model zero \(h=iy\), use
\(F_\Lambda=e^{h/2-\Lambda/8}G_\Lambda\) from the model zero theorem.
At a zero, the factors \(\partial_h F_\Lambda\) and
\(\mathcal D F_\Lambda\) equal the same nonzero prefactor times
\(\partial_h G_\Lambda\) and \(\mathcal D G_\Lambda\).
Even symmetry and real coefficients make the first of the
\(G_\Lambda\) expressions purely imaginary and the second real.
Simplicity and equation (12) of that theorem make both nonzero.
The imaginary product residual is therefore nonzero there.
Interval quadrature gives computable positive lower bounds on
these compact residuals as in the effective linear proof.
The \(O(n^{-1})\) bounds contradict the last residual for all
sufficiently large effective \(n\).

Taking the maximum of the computable thresholds in the three
ranges gives an absolute angular cutoff for two private linear
rows, and the original mixed and higher-private angle estimates
complete the theorem.
