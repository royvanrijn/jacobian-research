# Support-size factorial windows on separated homogeneous supports

Let \(b_1,\ldots,b_s\) be distinct nonnegative integer columns
in \(r\) variables. Assume they have a common total degree
and a common number of positive entries. For an integer \(D>0\),
put \(F_D(X)=\sum_i c_i X^{Db_i}\).

**Theorem.** For all sufficiently large \(D\), every positive
window of \(s'\) ordinary factorial moments of \(F_D\)
contains a nonzero moment, for every coefficient vector with
\(s'\) nonzero entries. The threshold is computable from the
exponent columns alone. The assertion holds in every
characteristic-zero field.

There is a uniform version for actual homogeneous integer
columns \(a_i\). Write \(b_i=a_i/D\) for a real scale \(D>0\).
Fix bounds \(r,s,0<\epsilon\le M\), and \(g>0\).
If every positive entry of \(b_i\) lies in \([\epsilon,M]\),
the same number of entries of each column is positive, and
\(\|b_i-b_j\|_\infty\ge g\) for \(i\ne j\), then the conclusion
holds for the original support \(X^{a_i}\) once \(D\) exceeds
a computable threshold in those bounds. No common divisor of
the entries is required.

For equal total degrees, the separate monomial saddle rates
can have equal magnitudes at all orders. The proof samples
a full window and applies a polynomial in the shift operator.
Nearby saddle frequencies whose separation is smaller than
the reciprocal order contribute with the same phase. The
shift polynomial removes the main terms of the other
frequencies. Derivatives of the saddle amplitudes control
the errors after this removal.

## A uniform analytic saddle amplitude

Discard zero coefficients. The one-term case is immediate.
Use the monomial peaks \(u_i(y)\) and \(P_i(y)\), disjoint
boxes \(V_i\), and logarithmic phases of
[the separated-support saddle proof](SFC_DILATED_SUPPORT_SADDLE_ENVELOPE.md#exact-integral-and-separation-of-the-peaks).
The exact integral after scaling \(X=mD\,y\) now has a
common factor in all the coefficients because all \(d_i\)
equal the same \(d\). Write
\[
                    q_i=c_i u_i(b_i)^D
\]
and rescale the polynomial by a common positive scalar
so that \(\max_i|q_i|=1\). Call \(i\) eligible if
\(|q_i|\ge e^{-\gamma D}\), with the fixed positive \(\gamma\)
from the peak estimate. Eligibility is independent of \(m\).
Write \(\lambda=mD\). If each column has \(r_0\) positive
entries, the common amplitude exponent is
\(\rho=r_0/2+r-r_0\). After nonzero scalar normalization,
the order-\(m\) moment is
\[
  M_m=\lambda^\rho\int_{\mathbb R_{\ge0}^r}
                   \left(\sum_i q_iP_i(y)^D\right)^m\,dy .
\]
In fact the scalar removed from the original factorial
moment is \(\lambda^{r-\rho+ m D d}\).

For an eligible \(i\), let \(\Phi_i\) be the phase on its
box and \((\theta_i,0)\) its local saddle, with inactive
coordinates evaluated at zero. Both are independent of
\(m\). Put
\[
  U_i=q_i e^{D\Phi_i(\theta_i,0)},\qquad
  B_i(\lambda)=\lambda^\rho e^{-\lambda\Phi_i(\theta_i,0)}
                      \int_{V_i}e^{\lambda\Phi_i(y)}\,dy .
\]
The original saddle estimates give
\[
 M_m=\sum_{i\ {\rm eligible}}U_i^m B_i(mD)+R_m,\qquad
 |R_m|\le C(mD)^\rho e^{-c mD},                              \tag{1}
\]
with constants independent of the coefficient vector and
\(m\). The omitted boxes and the outside integral were
bounded before any saddle deformation.

For use with finite differences, the saddle estimate also
holds throughout a fixed complex sector \(|\arg\lambda|<\theta_0\),
\(|\lambda|\ge1\):
\[
                 B_i(\lambda)=K_i+O(|\lambda|^{-1/2}),        \tag{2}
\]
where \(K_i\) is the Gaussian and endpoint amplitude of the
separated-support proof. Here \(K_i\) approaches the positive
pure-phase amplitude uniformly as \(D\) increases. To
verify sector uniformity, deform the active box by the same
translation as for real \(\lambda\). On each side contour the
phase difference from its saddle is close to a negative real
pure-phase value bounded away from zero. On the translated
box use the active-coordinate path from the saddle followed by
the inactive-coordinate path. The active Hessians and inactive
slopes are uniformly close to the respective real pure-phase
ones. For every \(|\alpha|\le\theta_0\), after a fixed reduction
of \(\theta_0\) their real parts after multiplication by
\(e^{i\alpha}\) retain a negative definite Hessian and strictly
negative inactive slopes. Integration along the two paths
therefore bounds the rotated real part of the phase difference
by \(-c'(|x|^2+\sum y_\nu)\). The quadratic-linear Taylor
polynomial has the same bound. The Taylor error
\(O(|x|^3+|x|\sum y_\nu+(\sum y_\nu)^2)\) gives integrated
error \(O(|\lambda|^{-\rho-1/2})\), using the segment formula
for the exponential difference and scaling by \(|\lambda|\).
The side and Gaussian-tail errors are exponentially small.
The leading Gaussian and endpoint integrals are
\(K_i\lambda^{-\rho}\) by analytic continuation from real
positive \(\lambda\), using the sector branch. This proves
(2) with a uniform bound in a closed smaller sector.

Let \(A_i(t)=B_i(tD)\). Cauchy estimates on discs of radius a
fixed multiple of \(t\), contained in this smaller sector,
now give, for \(t\ge1\) and \(1\le p\le s\),
\[
 |A_i(t)-K_i|\le C(tD)^{-1/2},\qquad
 |A_i^{(p)}(t)|\le C_p(tD)^{-1/2} t^{-p}.                     \tag{3}
\]
The saddle displacement, bounds on the phase and compact
parameter estimates used here are those in the preceding
saddle proof. For the uniform version with real \(b_i\), its
local holomorphic ratio argument applies because the actual
powers \(D(b_j-b_i)=a_j-a_i\) are integral. Its compactness
argument needs the lower separation of centers to be
\(g>0\), as assumed here, rather than a lower separation
of their total degrees.

## Filtering a positive window

Fix a window beginning at \(m\ge1\). All eligible \(|U_i|\)
are bounded by an absolute constant for the given support,
and at least one eligible \(|U_i|\) is at least
\(\exp(-C e^{-3\gamma D})\). This is because the saddle
values in the definition of \(U_i\) are
\(O(D^{-1}e^{-3\gamma D})\) and one raw \(|q_i|=1\).
Choose an eligible \(p\) maximizing
\[
         |U_p|^m\prod_{j\ne p}
                     \bigl(|U_p-U_j|+m^{-1}\bigr),             \tag{4}
\]
where all products and groups below use the eligible set.
Then \(|U_p|\) is bounded below by a fixed positive constant.
Indeed the product of any index in (4) is at least
\(m^{-(s-1)}\) and is uniformly bounded above, so comparison
with a raw maximizing index gives
\[
              |U_p|^m\ge c_0 m^{-(s-1)}
                                  \exp(-C m e^{-3\gamma D}).
\]
Taking \(m\)-th roots gives the bound uniformly for \(m\ge1\).

Choose a small positive \(\delta<1/4\), independently of
\(D,m\), as quantified below. There is a \(q\) with
\(1\le q\le 1+\binom s2\) such that none of the values
\(m|U_i-U_j|\) lies in \((\delta^{q+1},\delta^q]\).
This follows from the disjointness of these intervals and the
number of pairs. The relation \(m|U_i-U_j|\le\delta^{q+1}\)
is an equivalence: two such inequalities imply distance at
most \(2\delta^{q+1}<\delta^q\), and the empty interval gives
the required smaller bound. Let \(I\) be the class of \(p\),
and put \(P(z)=\prod_{j\notin I}(z-U_j)\). Its degree is at
most \(s'-1\), including the possibility \(P=1\). For \(i\in I\),
\[
 \left|(U_i/U_p)^m-1\right|\le C\delta,\qquad
 \left|P(U_i)/P(U_p)-1\right|\le C\delta.                      \tag{5}
\]
The first estimate uses the uniform lower bound on \(|U_p|\)
and \(m|U_i-U_p|\le\delta^{q+1}\). For the second, every
factor \(U_p-U_j\) with \(j\notin I\) has modulus more than
\(\delta^q/m\), so the relative change in each factor is at
most \(\delta\). By choosing \(\delta\) small compared with
the uniform positive bounds for the pure amplitudes \(K_i\),
and then \(D\) large, (2) and (5) imply
\[
       \left|\sum_{i\in I} U_i^m P(U_i)A_i(m)\right|
                             \ge c_1 |U_p|^m |P(U_p)|.       \tag{6}
\]
Indeed division by \(U_p^mP(U_p)\) leaves a sum of numbers
uniformly close to the positive pure amplitudes. Thus there is
no internal cancellation at the group selected by (4).

Let \(E\) denote the forward shift on sequences. In applying
\(P(E)\) to the \(i\)-th term of (1), the leading part with
\(A_i(m)\) frozen is \(U_i^mP(U_i)A_i(m)\). The error is bounded by
\[
            C(mD)^{-1/2}|U_i|^m
                    \prod_{j\notin I}(|U_i-U_j|+m^{-1}).     \tag{7}
\]
To prove this, write \(P(z)=\sum_l p_l z^l\) and Taylor expand
\(A_i(m+l)\) around \(m\) through degree \(s-1\). By (3) the
\(h\)-th Taylor coefficient, \(h>0\), is bounded by
\(C(mD)^{-1/2}m^{-h}\). The corresponding sum over the shift
coefficients is \((z\partial_z)^hP(U_i)\). Differentiating the
product for \(P\) multiplies the bound by at most a fixed
constant when each differentiated factor is replaced with
\(m^{-1}\) and the undifferentiated factors are bounded
by \(|U_i-U_j|+m^{-1}\). This gives (7) for every retained
Taylor term. The remainder is \(O((mD)^{-1/2}m^{-s})\)
because \(l\le s-1\), and the \(p_l\) and \(U_i\) are bounded.
The product in (7) is at least \(m^{-(s-1)}\), so the
remainder obeys that bound as well.

For clarity, let \(D_i=\prod_{j\notin I}(|U_i-U_j|+m^{-1})\)
and \(D'_i=\prod_{j\ne i}(|U_i-U_j|+m^{-1})\). For every \(i\),
\[
                   D_i\le m^{|I|-1}D'_i,\qquad
                   D_p\ge c_2m^{|I|-1}D'_p .
\]
For \(i\in I\), these comparisons remove \(|I|-1\) factors
bounded below by \(m^{-1}\); for \(p\) they are bounded above
by \((1+\delta^{q+1})m^{-1}\). For \(i\notin I\), the
comparison removes \(|I|\) such factors and inserts the
factor \(m^{-1}\) from \(j=i\). Hence (4) gives
\(|U_i|^mD_i\le c_2^{-1}|U_p|^mD_p\).
Also \(D_p\le C_\delta |P(U_p)|\), by the lower separation
\(\delta^q/m\) from every outside frequency and the fixed
upper bound on \(q\). Thus the sum of (7) is
\(O((mD)^{-1/2}|U_p|^m|P(U_p)|)\).

The filtered remainder in (1) is bounded by
\(C((m+s)D)^\rho e^{-c mD}\) since \(P\) has bounded
coefficients. From (4), \(D_p\le C_\delta |P(U_p)|\),
and the estimate at one raw maximizing index,
\[
              |U_p|^m|P(U_p)|\ge
                      c_3 m^{-(s-1)} e^{-C m e^{-3\gamma D}} .
\]
The filtered remainder is uniformly smaller than a fixed
fraction of this lower bound for all \(m\ge1\) as \(D\) grows.
Combining with (6) proves \(P(E)M_m\ne0\). The shift
polynomial has degree at most \(s'-1\), so the moments in the
window cannot all vanish.

Every bound uses fixed compact peak margins, a fixed
\(\delta\) and polynomial or exponential inequalities in \(D\).
The sector estimates and the rational-box constructions in
the separated-support proof make these thresholds computable.
As usual one embeds the subfield generated by the finitely many
coefficients into \(\mathbb C\) to obtain the
characteristic-zero statement. The common number of
positive coordinates gives the common normalization
\(\lambda^\rho\) in (1).

## Density one on homogeneous supports

Fix \(r,s\) and a zero pattern with \(r_0\ge2\) positive entries
in each column. Choose the positive entries of each column
independently and uniformly among the positive compositions
of \(H\) into \(r_0\) parts. As \(H\) tends to infinity, the
fraction of supports satisfying the actual support-size
window at every order tends to one, uniformly over all
coefficient vectors. This includes fully shared homogeneous
supports, and therefore a density-one fraction of all
homogeneous exponent arrays for fixed dimension \(r\ge2\)
and fixed support size \(s\). The thresholds for any
prescribed tolerance in the fraction are computable.

For any part the probability of taking a specified integer
value is at most
\[
                 \frac{\binom{H-2}{r_0-2}}
                      {\binom{H-1}{r_0-1}}
                       =\frac{r_0-1}{H-1}.
\]
Consequently a fraction at most
\(rs(r_0-1)(\epsilon H+1)/(H-1)\) has an entry in the
interval \((0,\epsilon H)\). For each pair of columns,
condition on the second column. Choose any active coordinate
of the first and use the same marginal bound to estimate
the probability that the two column vectors have sup
distance less than \(gH\). This is at most
\((2gH+2)(r_0-1)/(H-1)\), even if the second column has zero
at the chosen coordinate. Sum over pairs. The remaining
supports satisfy the uniform theorem with \(D=H\) and \(M=1\).
Choose \(\epsilon,g\) small and then \(H\) large. The
estimate on pairs includes all identical columns.
Among weak compositions of \(H\) into \(r\ge2\) parts,
the fraction with a zero in any fixed part equals
\((r-1)/(H+r-1)\) for that part. Thus the fraction of all
arrays outside the fully shared pattern tends to zero.
For \(r_0=1\), a pattern with distinct active coordinates
for all columns already has pairwise distance \(H\),
so its only possible array at each \(H\) satisfies the
uniform theorem for sufficiently large \(H\).
