# Effective first-sector stability under mixed tails

Start with an ordinary two-monomial core whose exponent columns satisfy
\(a_i\le b_i\) on every row. Append shared ordinary rows \((a_i,b_i)\),
\(a_i,b_i\ge1\), of either orientation. Define the appended curvature
and signed degree change by
\[
       \tau=\sum_{\rm tail}\frac{(b_i-a_i)^2}{\max(a_i,b_i)},
       \qquad D=\sum_{\rm tail}(b_i-a_i).
\]
Assume \(D\ge0\). There is an absolute, computable positive
\(\tau_0\) such that if \(\tau<\tau_0\), the full ordinary
moment polynomials cannot satisfy \(P_n(z)=P_{n+1}(\lambda z)=0\)
with \(\lambda\ge1\) and \(|\arg z|\le 2\pi/(n+1)\), for any
\(n\ge2\). This includes comparison at the same complex point.
The condition on \(\lambda\) uses the orientation in which the
core gaps are nonnegative.

The bound is given by a terminating real-algebraic procedure in
terms of the [absolute eventual angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
Its numerical value is unevaluated. The procedure only tests
compact families of continuous parameters and uses the
[first-sector comparison for comparable columns](MULTIVARIATE_COMPARABLE_FIRST_SECTOR_HCM_RADIAL_COMPARISON.md)
to certify positive separation throughout those families.

## Bounds at a proposed common ray

Choose \(J\ge4\) at least the absolute eventual angular cutoff,
and put \(R=1+2^J\). The eventual theorem excludes every candidate
with terminal \(N=n+1\ge J\), and the
[absolute curvature theorem](MULTIVARIATE_ABSOLUTE_CURVATURE_27.md)
excludes every candidate with full curvature at least \(27\).
For distinct full columns, the first and second moments cannot
share a ray, since the quadratic moment is strictly positive
on the real axis. If \(\tau=0\), the tail is diagonal and the
first-sector comparison for the comparable core applies directly.
If that core is also diagonal, its only root is negative real
and outside the indicated sector. Hence we may assume
\(3\le N<J\), full curvature below \(27\), and \(0<\tau<1/10\).

For order \(j\), let \(G_j\) and \(D_j\) be the positive endpoint
scales of the core and the tail, respectively. The \(j\)-th powers
of these scales are the corresponding leading-to-constant
coefficient ratios. Let \(\Pi_j\) denote the core polynomial
in its endpoint normalization. The full normalized polynomial
has coefficients \([w^k]\Pi_j(w)\) multiplied by normalized
tail factors \(t_{j,k}\). Its endpoints are one and every
coefficient is bounded by \(\binom jk\) by Hölder in the positive
moment representation. Cauchy's root bound and reversal put
each full normalized root \(w_j\) in
\[
                           R^{-1}\le |w_j|\le R.                \tag{1}
\]
For the assumed zeros at \(z,\lambda z\), put
\(w=w_n=G_nD_nz\), so
\(w_N=\lambda (G_N/G_n)(D_N/D_n)w\).

The
[order derivative and chord estimates for the tail](MULTIVARIATE_COMPARABLE_CORE_EFFECTIVE_TAIL_STABILITY.md#one-physical-comparison-point)
give
\[
  0\le 1-t_{j,k}\le \tau j/4,\qquad
  \left|\log\frac{D_N}{D_n}-D\log\frac{N}{n}\right|
                          \le\frac{\tau}{nN(1-\tau)}<\tau .      \tag{2}
\]
For clarity, the endpoint scale is
\(\log D_j=\sum d_i\int_0^1\psi(1+j(a_i+u d_i))\,du\).
The bound \(0\le j^{-1}-x\psi'(1+jx)\le j^{-2}/x\) follows
from comparison of the trigamma series with integrals.
Integration of the order derivative from \(n\) to \(N\) gives
(2), because
\(\sum |d_i|/\min(a_i,b_i)\le\tau/(1-\tau)\).
The normalized tail log array is the endpoint chord difference
of a function with second derivative in the coefficient index
between zero and \(\tau/[j(1-\tau)]\). Its absolute value is
bounded by \(\tau j/[8(1-\tau)]\), which gives the first
bound in (2).

The
[shared row scale estimate](MULTIVARIATE_UNIFORM_ROOT_ANNULUS.md#shared-rows-and-large-degree-imbalance)
writes the shared part of the core scale ratio in logarithmic
coordinates as \(\Delta_{\rm sh}\log(N/n)\) with error at most
\(783/(nN)\). Here \(\Delta_{\rm sh}\) is the sum of the shared
core gaps. To see the constant, write \(s=d^2/(a+d)\) for
an increasing shared row with \(a\ge1\). Then \(d/a\le1+s\)
and \(d^2/a\le s(2+s)\). Summing at curvature below \(27\)
gives \(783\). Each private core row has a nondecreasing
endpoint scale in the order. In fact for any real \(a\ge0\),
integer \(d>0\), the row's log scale is
\(\int_0^d\psi(1+j(a+x))\,dx\), which increases with \(j\).
The sum of private core gaps is at most \(26\), since they
are positive integers and contribute their full gaps to curvature.

From (1), \(\lambda\ge1\), and (2) we obtain
\[
   (\Delta_{\rm sh}+D)\log(N/n)
       \le 2\log R+783/(nN)+\tau < 2(J+1)+132 .
\]
Since \(\log(N/n)>1/N\) and \(J\ge4\), the complete list
\(d_1,\ldots,d_r\) of positive core gaps and the tail twist
\(D\) satisfy
\[
                       d_1+\cdots+d_r+D\le 2J(J+80).            \tag{3}
\]
Diagonal core rows can be dropped. Define
\[
       \widehat w_N=\lambda(G_N/G_n)(N/n)^D w .
\]
By (2), \(\widehat w_N\) differs from \(w_N\) by a positive
rescaling with log of modulus less than \(\tau\).
In particular \(|\widehat w_N|\le2R\) and
\(1\le\lambda\le2R^2\), using \(G_N/G_n\ge1\) and \(D\ge0\).
At the candidate zero the residuals obey
\[
  |\Pi_n(w)|\le C_J\tau,\qquad
  |\Pi_N(\widehat w_N)|\le C_J\tau,\qquad
                     C_J=2J(1+2R)^J .                          \tag{4}
\]
Indeed the tail coefficients cost at most
\(\tau J(1+R)^J/4\) at each full normalized root by (2).
Rescaling the upper root costs at most
\(\tau J(1+2R)^J\) by integrating the derivative in log
radius along the positive radial segment. The derivative
bound follows directly from the binomial bounds for
the normalized core coefficients.

## A compact algebraic bound for the core residuals

For each \(3\le N<J\), and each finite list of positive gaps
and integer \(D\ge0\) satisfying (3), replace each actual
baseline \(a_i\) by the coordinate \(x_i=1/(1+a_i)\). Allow
the whole interval \(0\le x_i\le1\). For a gap \(d\) define
\[
 F_{j,k}^{(d)}(x)=\prod_{h=1}^{kd}\big[j(1-x)+hx\big].
\]
This is strictly positive on the closed interval, including
the empty product at \(k=0\). Its normalized row factor at
index \(k\) is
\[
 \frac{F_{j,k}^{(d)}(x)}{F_{j,j}^{(d)}(x)^{k/j}},
\]
and its cross-order scale ratio is
\[
                 \frac{F_{N,N}^{(d)}(x)^{1/N}}{
                          F_{n,n}^{(d)}(x)^{1/n}}.
\]
The positive real roots in these formulas have unique values.
At \(x=0\) the normalized array is one and the scale ratio
is \((N/n)^d\). The polynomials \(\Pi_j\) and the core scale
ratio \(A=G_N/G_n\) therefore extend as continuous
semialgebraic functions of the entire baseline box. Each
\(x_i>0\) represents a real nonnegative baseline and an
integer positive gap by the same Gamma formulas as before.
The beta coupling and density proof of the first-sector
comparison use only \(a_i\ge0\), \(d_i\ge1\) and the
corresponding Gamma shapes, so they apply to these real
baselines as well.

For each list take the closed domain with \(x_i\in[0,1]\),
\(1\le\lambda\le2R^2\), \(R^{-1}\le|w|\le R\),
\(\operatorname{Re} w\ge |w|\cos(2\pi/N)\), and
\(|\lambda A(N/n)^D w|\le2R\). On this compact domain the
function
\[
     |\Pi_n(w)|^2+
             \big|\Pi_N\big(\lambda A(N/n)^D w\big)\big|^2      \tag{5}
\]
is strictly positive. To check this, remove the coordinates
with \(x_i=0\). Their only cross-order effect is to multiply
the additional dilation by the product of their factors
\((N/n)^{d_i}\). If the remaining core is nonempty, (5) can
vanish only when its lower and upper physical moments vanish
at two arguments on the same ray with ratio
\(\lambda(N/n)^{D+\sum_{x_i=0}d_i}\ge1\).
The first-sector comparison excludes that point. If every
coordinate is removed, the normalized lower polynomial is
\((1+w)^n\), which has no zero in the domain.

All these domains and residuals have real-algebraic data.
For the angle boundary use the designated real algebraic
number \(\cos(2\pi/N)\), and encode every positive root in
the scales by its power equation and positivity. Real
quantifier elimination can therefore test any rational
lower bound. There are finitely many gap lists and degrees,
so compactness and strict positivity ensure that a search over
dyadic \(b>0\) terminates with (5) strictly above \(2b^2\) on
each domain, including the vacuous conditions for empty domains.
Set
\[
                     \tau_0=\min\{1/10,\ b/C_J\}.
\]
A candidate at \(\tau<\tau_0\) would satisfy (4) and put (5)
below \(2b^2\), a contradiction. This establishes the
effective first-sector tolerance uniformly in the ordinary
core and appended row data.
