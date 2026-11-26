# Second-order contact tests for weak ordinary factorial tails

Fix an ordinary two-monomial core \(C\), and write \(Q_j\)
for its factorial moment at degree \(j\). Append ordinary
unequal shared rows. For an appended row put
\(m_i=(a_i+b_i)/2\) and \(d_i=b_i-a_i\), with signed \(d_i\).
Use the tail aggregates
\[
 W=\sum_i d_i^2/m_i,\qquad K=\sum_i d_i^3/m_i^2,\qquad
 Z=\sum_i d_i^2/m_i^2,\qquad h=\sum_i d_i .
\]
The first-order radial response at a limiting integer-twist
contact is described in the
[algebraic trajectory theorem](MULTIVARIATE_ALGEBRAIC_PROFILE_TRAJECTORIES.md#signed-radial-response-in-the-vanishing-tail).
The following result retains the second-order angular response
of the ordinary rows and applies to arbitrary weak-tail sequences.

## Exact second-order response of a simple core root

Let \(L_j(s)\) be the endpoint-normalized tail log coefficient
per degree, with the binomial coefficient stripped. Define
polynomials of the coefficient index \(k\) by
\[
 \begin{split}
 F_j(k)&=\frac{k^2-jk}{2j},\qquad G_j(k)=-\frac{F_j(k)}{2j},\\
 H_j(k)&=\left(\frac{k}{j}-\frac12\right)\frac{k(j-k)}{6j}.
 \end{split}
\]
Then as \(W\to0\),
\[
 jL_j(k/j)=WF_j(k)+ZG_j(k)+KH_j(k)+O_j(W^3).                \tag{1}
\]
The error is uniform over all tails, without a row-count or
gap bound.

To verify this, put \(t=s-1/2\) and \(x_i=m_i+d_i t\).
The exact symbol curvature is
\(L_j''(s)=\sum_i d_i^2 j\psi'(1+jx_i)\).
The trigamma integral with
\(0\le u/(e^u-1)-1+u/2\le u^2/12\) gives
\[
 0\le j\psi'(1+jx)-\frac1x+\frac1{2jx^2}
                         \le\frac1{6j^2x^3}\qquad(x>0).
\]
The cited [first-order proof](MULTIVARIATE_ALGEBRAIC_PROFILE_TRAJECTORIES.md#signed-radial-response-in-the-vanishing-tail)
gives the scalar integral bound. For \(W\le1/2\), the
integer gaps give \(|d_i|/m_i\le W\), \(1/m_i\le W\), and
\(x_i\ge3m_i/4\) on \(0\le s\le1\). Expand
\(x_i^{-1}=m_i^{-1}-d_it/m_i^2+O(d_i^2/m_i^3)\),
and \(x_i^{-2}=m_i^{-2}+O(|d_i|/m_i^3)\).
The first error has constant at most \(1/3\), by the
exact geometric remainder and \(|t|\le1/2\).
The second has constant at most \(5/2\), by differentiating
\((1+y)^{-2}\) on \(|y|\le1/4\).
Since \(\sum_i |d_i|^a/m_i^3\le W^3\) for \(a=2,3,4\),
these bounds give
\[
          L_j''(s)=W-K(s-1/2)-Z/(2j)+O(W^3)
\]
with error at most \(2W^3\) for \(j\ge2\).
Integrate twice with the two zero endpoint values.
The affine curvature term uses
\(((s-1/2)s(s-1))''=6(s-1/2)\).
This proves (1).

For a simple root \(\zeta_j\ne0\) of \(Q_j\), let
\(\mathcal D=w\,d/dw\) and use the notation
\([T]=T(\mathcal D)Q_j(\zeta_j)\) for polynomial operators
\(T\). Put \(q_1=[\mathcal D]\), \(q_2=[\mathcal D^2]\), and
\[
 \begin{split}
 \alpha_j&=-[F_j]/q_1,\\
 \beta_j&=-\bigl(q_2\alpha_j^2/2+
               \alpha_j[\mathcal D F_j]+[F_j^2]/2\bigr)/q_1,\\
 \gamma_j&=-[G_j]/q_1=-\alpha_j/(2j),\qquad
 \kappa_j=-[H_j]/q_1 .
 \end{split}                                                 \tag{2}
\]
In these brackets \(\mathcal D\) in an operator denotes its
polynomial argument. All constants are algebraic in the exact
core coefficients and its isolated root. The corresponding root
\(\zeta_j^{\rm tail}\) after normalizing the tail alone satisfies
\[
 \log(\zeta_j^{\rm tail}/\zeta_j)
     =W\alpha_j+W^2\beta_j+Z\gamma_j+K\kappa_j+O_{C,j,\zeta_j}(W^3).
                                                               \tag{3}
\]
In fact \(Z,|K|\le W^2\). Expand the coefficient exponentials
in (1), insert \(w=\zeta_j e^u\), and compare through degree
two in \(W\), keeping the linear terms in \(K,Z\).
The coefficient of \(u\) is \(q_1\ne0\); solving gives (2).
For a uniform effective remainder, isolate the simple root in
a rational disk excluding all other core roots and zero, with
radius less than \(|\zeta_j|/2\). Rouché on its boundary gives
a holomorphic single root for all coefficient vectors in a
small complex polydisk about the core, using the product lower
bound on \(|Q_j|\) there. Its logarithmic displacement has
bounded derivatives on a smaller polydisk by Cauchy's formula.
Apply this analytic map to the coefficient expansion (1).
This supplies effective constants in (3).

## A finite algebraic gate at a limiting twist

Take adjacent \(n,N=n+1\) and suppose that along a sequence
with \(W\to0\) the signed tail gap is the same integer \(h\).
At a limiting physical contact of simple roots write
\(\zeta_N=(N/n)^h\zeta_n\), and let
\(\Delta\alpha=\alpha_N-\alpha_n\), and similarly for the
other three responses. Use \(\xi_i=1/d_i\) and the normalized
tail curvature weights \(w_i=d_i^2/(m_i W)\). Then
\[
 \mu=\sum_i\xi_i w_i,\qquad
 X=K/W^2=\sum_i\xi_i w_i^2,\qquad
 Y=Z/W^2=\sum_i\xi_i^2w_i^2.
\]
They are bounded real numbers. Along any convergent subsequence
their limits must satisfy
\[
 2nN\Delta\alpha=-\mu_*,\qquad
 \Im\!\left(\Delta\beta+Y_*\Delta\gamma+X_*\Delta\kappa\right)=0.
                                                               \tag{4}
\]
The first equation is the cited signed radial response.
For the second, the ratio of the two tail-normalized roots at a
physical contact is positive real. Subtract (3), use the first
equation to see that \(\Im\Delta\alpha=0\), and divide the
remaining imaginary equation by \(W^2\).
The error tends to zero uniformly.

If the allowed tail has a bound on \(\sum_i|d_i|\), enumerate
the finite signed gap lists with that bound and sum \(h\).
For each list let \(w_i\ge0\), \(\sum_iw_i=1\), and test (4)
using the displayed polynomials in \(w_i\). These are finite
real polynomial tests over the algebraic numbers. This retains
the simultaneous signed moment constraints of the same rows.
If infeasible for every list at a fixed simple contact, a
computable weak-tail tolerance excludes that contact. To obtain
a tolerance, bound the sum of squares of the two residuals in
(4) below on each compact simplex by real quantifier elimination,
and use the uniform remainders in (1)–(3) and in the first-order
radial response. At fixed core and gap budget there are finitely
many twists, fixed terminals below the
[absolute angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md),
and algebraic core roots to inspect.

There is also a finite outer gate for arbitrary tail sizes.
Put \(b_i=\xi_iw_i\), let \(P_+=\sum_{b_i>0}b_i\),
\(P_-=\sum_{b_i<0}|b_i|\), and let \(Y_+,Y_-\) sum
\(b_i^2\) on those two signs. In every limit
\[
 \begin{gathered}
 P_\pm\ge0,\quad P_++P_-\le1,\quad P_+-P_-=\mu_*,
 \quad 0\le Y_\pm\le P_\pm^2,\quad Y_++Y_-=Y_*,\\
 Y_+-Y_- -P_-(1-P_+-P_-)\le X_*
       \le Y_+-Y_- +P_+(1-P_+-P_-).
 \end{gathered}                                             \tag{5}
\]
Indeed \(w_i\ge|b_i|\) and
\(X=\sum_i b_i|b_i|+\sum_i b_i(w_i-|b_i|)\).
The second sum is bounded using the remaining mass
\(1-P_+-P_-\) and \(-P_-\le b_i\le P_+\).
Thus infeasibility of (4)–(5) excludes all weak ordinary tails
at the contact, without a gap budget.

At a reciprocal core with a simple unit-circle root, the
Gaussian and quadratic perturbations are reciprocal, so
\(\alpha_j,\beta_j,\gamma_j\) are purely imaginary.
The skew operator \(H_j\) reverses sign under \(k\mapsto j-k\).
Reciprocal conjugation therefore sends the nearby root at a
real skew parameter \(u\) to the nearby root at \(-u\).
Differentiating the logarithmic root position gives
\(\kappa_j\) real. In this case (4) requires
\(\Delta\alpha=0,\ \mu_*=0\) and reduces the remaining angular
test to
\(\Im\Delta\beta+Y_*\Im\Delta\gamma=0\).
This is a simplification of the exact gate. When both
normalized core arrays have all roots on the unit circle and
their raw endpoint scales agree, a nonzero physical integer
twist already has different root moduli.

## A Gaussian jet criterion without simple roots

There is also a higher-order sufficient test for weak tails
using only one algebraic polynomial in the coefficient arrays.
Fix \(n,N\), and use the
[off-diagonal root quotient resultant](MULTIVARIATE_DIFFUSE_PHYSICAL_EXCEPTION_POINTS.md#the-finite-set-of-physical-reference-contacts),
denoted here by \(\mathcal A(c)\), on the pair of coefficient
arrays \(c\). Work on the affine space where the endpoints
are the fixed core endpoints. Equal root arguments of the
actual pair force \(\mathcal A(c)=0\): the even ordinary moment
excludes real roots and a common nonreal ray gives a common
quotient of a root and its conjugate. The criterion is invariant
under the tail's positive normalizations.

Let \(c_0\) be the exact core coefficient vector and let \(r\)
be the least nonzero total Taylor degree of \(\mathcal A\)
at \(c_0\) in these affine coordinates. In the Gaussian pair
the coefficient multipliers are \(q^{(nN/j)k(j-k)}\), where
\(q=e^{-W/(2nN)}\). Let \(\nu\) be the multiplicity at \(q=1\)
of the corresponding polynomial in \(q\). Its value at \(q=0\)
is nonzero, since the endpoint binomials have disjoint
nontrivial root-of-unity quotients. Hence \(r,\nu\) are both
finite and exactly computable. In particular \(r\le\nu\).

For every sufficiently small positive \(W\) there is an effective
\(c_1>0\), depending only on the core and orders, such that the
Gaussian resultant has modulus at least \(c_1 W^\nu\).
Factor out \((q-1)^\nu\), bound the remaining polynomial away
from zero near one, and use \(1-e^{-W/(2nN)}\ge W/(4nN)\).
On the other hand let \(d_*=\min_i|d_i|\) in a nonempty tail.
The curvature remainder above improves to \(2W^3/d_*^2\),
since \(\sum_i|d_i|^a/m_i^3\le W^3/d_*^2\) for \(a\le4\).
Also \(|K|\le W^2/d_*\), \(Z\le W^2/d_*^2\).
It follows from the Green chord bound that the actual
coefficient vector differs from the Gaussian one by
\(O_{C,n,N}(W^2/d_*)\). Both vectors have distance \(O_{C,n,N}(W)\)
from \(c_0\). For \(r\ge1\), the gradient of \(\mathcal A\) on
the segment between them has size \(O_{C,n,N}(W^{r-1})\)
by its Taylor expansion. For \(r=0\) use a bounded gradient.
In either case the two resultant values differ by at most
\(c_2 W^{r+1}/d_*\) for an effective constant \(c_2\).

Consequently
\[
                    d_*>(c_2/c_1)\,W^{r+1-\nu}
\]
excludes angular contact for all sufficiently small \(W\).
In particular if \(\nu=r\) then every sufficiently weak
ordinary tail is angularly separated at these orders,
including when the core has multiple roots. All constants
and orders in this criterion concern a fixed exact core.
Together with the absolute angular cutoff it gives
all-degree separation when the finitely many remaining
pairs meet the stated condition.
