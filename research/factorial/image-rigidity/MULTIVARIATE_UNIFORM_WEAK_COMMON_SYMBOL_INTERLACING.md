# Common-symbol interlacing at weak mirrored curvature

Consider a distinct fully shared mirrored ordinary profile and set
\[
    W=\sum_i\frac{d_i^2}{\ell_i+d_i/2},
\]
where the sum counts both orientations of each mirrored row.
Let \(L_j(s)\) denote the order-\(j\) normalized log-factorial
symbol per degree. There is a computable absolute positive
curvature cap \(W_*\) such that, for every \(N=n+1\ge2\) and
\(0<W\le W_*\), the two polynomials
\[
      Q_m(z)=\sum_{k=0}^m\binom mk e^{m L_N(k/m)}z^k,
                      \qquad m=n,N ,
\]
have simple unit-circle roots and strictly interlace across
degrees. The exact upper symbol is used at both degrees.

The proof retains a lower bound on the actual adjacent angular
gap in the weak saddle comparison. The
[uniform symmetric jump estimate](MULTIVARIATE_UNIFORM_WEAK_MIRRORED_LAPLACE_TRANSFER.md)
then makes the change from the lower actual symbol to the upper
one on the lower grid too small to cross any upper root argument.
This argument gives a common-symbol theorem on the weak band.

## A robust adjacent gap for the actual mirrored profiles

Write \(\varepsilon=\sqrt{W/2}\). In the
[uniform weak saddle comparison](MULTIVARIATE_UNIFORM_WEAK_CURVATURE_ANGLE_SEPARATION.md),
a root variable has the form \(z=-e^{\varepsilon\xi}\) for the
mirrored profile. The ordinary radial center in that proof is
one because the signed row contributions cancel. The
[mirrored circle-root theorem](MULTIVARIATE_MIRRORED_ROOT_CIRCLE.md)
makes \(\xi=iy\) purely imaginary. Choose a sufficiently small
fixed weak cap so the large-field and one-saddle exclusions of
the weak comparison cover all roots outside a fixed strip
\(|y|\le Y\), with \(\varepsilon Y<\pi/2\).

There are absolute computable \(g>0\), \(n_0\ge2\), and a
possibly smaller absolute weak cap, such that all adjacent
actual root heights obey
\[
                       |y_N-y_n|>g n^{-5/3}
                                  \qquad (n\ge n_0).          \tag{1}
\]
To verify this, suppose \(|y_N-y_n|\le g n^{-5/3}\).
It is enough to consider nonnegative heights: if the signs are
opposite, reflection gives a pair of nonnegative heights with
distance no larger. For large \(n\) the two heights lie in one
common admitted chart, by the fixed overlap of the regular
and fold charts in the weak saddle proof.

In a common regular chart, that proof gives the two root
equations with the same analytic phase difference \(D\) and
nonzero amplitude \(\Gamma\):
\[
    \Gamma(iy_j)\exp(jD(iy_j))=-(1+e_j),\qquad
                            |e_j|\le K/\sqrt n,\quad j=n,N.
\]
The compact chart supplies a fixed positive lower bound for
\(|e^{D(iy_n)}-1|\) and fixed upper bounds for the derivatives
of \(D\) and \(\log\Gamma\) on the joining segment. Dividing the
two equations, the log error from the two arguments is
\(O(N|y_N-y_n|)=O(g n^{-2/3})\). The separate relative errors
are \(O(n^{-1/2})\). Thus \(e^{D(iy_n)}\) tends to one uniformly
under the proposed near coincidence, contradicting the margin.

In the fold chart let \(\Lambda\) be the admitted analytic
Airy parameter. The
[Airy label estimate](MULTIVARIATE_UNIFORM_WEAK_CURVATURE_ANGLE_SEPARATION.md#airy-labels-at-the-same-imaginary-field)
places \(j^{2/3}\Lambda(iy_j)\) within \(O(n^{-1/3})\) of an
Airy zero. The near coincidence gives
\[
 |N^{2/3}\Lambda(iy_N)-n^{2/3}\Lambda(iy_n)|
     \le N^{2/3}\|\Lambda'\|\,g n^{-5/3}
            +(N^{2/3}-n^{2/3})\|\Lambda\| .
\]
The first term is \(O(g/n)\). The same fixed fold shrinkage
and Airy spacing inequality as in that proof force the same
Airy label. For that label the analytic order branch has
\[
  \xi'(j)=(\Lambda^{-1})'
        \left(-\tfrac23 a_k j^{-5/3}+O(n^{-2})\right)
                        \qquad(n\le j\le N).
\]
The imaginary part of the inverse derivative is uniformly
positive and \(-a_k\) is uniformly positive. Integration gives
an actual height gap at least a positive fixed multiple of
\(n^{-5/3}\), since the exact roots differ from their analytic
branch by exponentially small errors. Choose \(g\) smaller
than that multiple and increase \(n_0\). Both comparisons use
strict effective margins on fixed charts, which proves (1).

This argument also applies when the smaller exponents \(\ell_i\)
are positive reals at least one and the gaps \(d_i\) are
positive integers. The relevant integrality bounds in the
[effective density and saddle estimates](MULTIVARIATE_EFFECTIVE_VANISHING_CURVATURE.md)
use \(|d_i|\ge1\) for unequal rows and the order \(j\) being
integer at the actual moments. The signed-degree difference
vanishes by reflection. The Gamma identities and estimates
allow real smaller exponents, and the annulus bound in the
imaginary field follows directly from the mirrored circle
location.

## Replacing the lower symbol without crossing

The physical angular gap from (1) exceeds
\(g\sqrt{W/2}\,n^{-5/3}\). For \(W\le10^{-7}\), the uniform
jump estimate moves the lower-order roots by at most
\(13000W^{3/2}(1+\log n)/n^2\) when \(L_n\) is replaced by
\(L_N\) on the lower grid. Since
\((1+\log n)/n^{1/3}\le3e^{-2/3}<2\), taking
\(W\le g/100000\) makes the movement less than half the actual
adjacent gap. This gives disjoint arguments for the common
symbol whenever \(n\ge n_0\).

## The bounded degrees and the orientation of interlacing

The same conclusion for the finitely many remaining degree
pairs uses the ordinary mirrored midpoint hierarchy. To see
its uniformity, freeze an order \(j\ge2\), put
\(c=\tfrac12L_j''(1/2)\) and write \(x=s-1/2\). By the
[positive row density representation](MULTIVARIATE_MATCHED_LAPLACE_NONLINEAR_TRANSFER.md#positive-curvature-increment-for-ordinary-shared-rows),
grouping each mirrored pair once gives
\[
 c=\sum_i d_i \int_0^\infty e^{-a_i u} f(u/(d_i j))\,du,
 \quad
 L_j''(1/2+x)=2\sum_i d_i \int_0^\infty
                      e^{-a_i u} f(u/(d_i j))\cosh(ux)\,du,
\]
where \(a_i=\ell_i/d_i+1/2\) and \(f(t)=t/(e^t-1)\).
Use \(e^{-t}\le f(t)\le1\), \(j\ge2\), and \(a_i\ge1/2\)
to obtain \(\sum_i d_i/a_i\le2c\) and
\(\max_i a_i^{-1}\le2c\). Also \(c\le W/2\).
Integrating the even Taylor expansion gives
\[
 L_j(1/2+x)-L_j(1/2)=c x^2+\sum_{r\ge1}\lambda_r x^{2r+2},
 \quad 0\le\lambda_r\le
               \frac{4\cdot4^r c^{2r+1}}{(2r+1)(2r+2)}.     \tag{2}
\]
The estimate follows by integrating the \(2r\)-th Laplace
moment and then integrating twice in \(x\).

For a fixed degree \(m\), set \(x_k=k/m-1/2\) and introduce
an auxiliary complex variable \(\tau\). Define
\[
 A_\tau(x)=\tau^2x^2+
       \sum_{r\ge1}\frac{\lambda_r}{c^{2r+1}}
                       \tau^{4r+2}x^{2r+2},\qquad
 F_m(\tau,y)=\sum_{k=0}^m(-1)^k\binom mk
                 e^{m A_\tau(x_k)+i\tau my x_k}.
\]
At \(\tau=\sqrt c\) this is the order-\(m\) common-symbol
polynomial in angle \(z=-e^{i\sqrt c\,y}\) after removal of
nonzero factors. The coefficient of \(\tau^w\) in the
exponential is a polynomial in \(x_k\) of degree at most \(w\).
Every occurrence of a nonquadratic term reduces the maximal
degree by at least two. The alternating binomial sum kills
all terms of degree less than \(m\). At weight \(m\) only the
quadratic and angle terms survive. Reflection in \(k\) makes
the coefficient of weight \(m+1\) vanish as well. Since
\(\sum_k(-1)^k\binom mk x_k^m=(-1)^m m!/m^m\), the leading
coefficient is
\[
             (-i)^m(2/m)^{m/2}\mathrm{He}_m(\sqrt{m/2}\,y).
\]
For \(|\tau|\le1/2,\ |x|\le1/2\), (2) gives
\(|A_\tau(x)|\le1/16+1/60<1/10\).
On \(|y|\le R\) the Cauchy estimate for the remaining weights
then bounds the error in \(F_m(\sqrt c,y)/c^{m/2}\) by
\(8\cdot4^m e^{m(1/10+R/4)}c\) for \(\sqrt c\le1/4\).
All these bounds are independent of the profile and frozen
order.

The [scaled Hermite interlacing theorem](MULTIVARIATE_SCALED_HERMITE_INTERLACING.md)
applies to these two limit arrays, since the upper variance
relative to the lower is \(n/(n+1)\). For every fixed
\(2\le n<n_0\), isolating disks around the limiting Hermite
roots and the error bound above give an effective uniform
positive cap on \(c\), and hence on \(W\), for strict
interlacing. The case \(n=1\) follows from the linear factor
\(1+z\) and the reciprocal quadratic with middle coefficient
strictly between zero and two.

Finally, for any covered profile, continuously replace
\(\ell_i\) by \(t\ell_i\) for \(t\ge1\), retaining each gap
and its mirror. The hypotheses remain valid and \(W\)
decreases to zero. The common-symbol roots stay simple and
on the circle by the uniform same-grid and mirrored circle
theorems. The common-symbol pairs have disjoint arguments
by the preceding large- and bounded-degree comparisons.
Their cyclic order can therefore never change.
The small-\(W\) Hermite limit identifies that order as strict
interlacing, completing the proof.

## Sharp-clock Gaussian corollary

There is a computable absolute \(c_*>0\) such that the Gaussian
polynomials
\[
         R_m(z;c/m)=\sum_{k=0}^m\binom mk
                          e^{-c k(m-k)/m}z^k,\qquad m=n,n+1,
\]
strictly interlace for every \(n\ge1\) when \(0<c\le c_*\).
In particular the lower clock ratio \(n/(n+1)\) is admitted
throughout a fixed weak interval, independent of the degree.

For the [diffuse mirrored unit-gap sequence](MULTIVARIATE_GAUSSIAN_PHYSICAL_CLOCK_COMPARISON.md#connection-with-the-matched-ordinary-tanh-moments),
take \(\lfloor c\ell\rfloor\) copies of the pair with smaller
exponent \(\ell\) and gap one, and let \(\ell\) tend to infinity
through positive integers. Then
\(W_\ell=2\lfloor c\ell\rfloor/(\ell+1/2)\) tends to \(2c\)
and is at most \(2c\). At each fixed degree \(m\) the actual
polynomial converges to \(R_m(z;c/m)\). The same-grid transfer
also shows that the upper-symbol polynomial on the lower grid
has the same limiting ordered root locations: before replacing
\(1/\sqrt L\) by \(\sqrt W\), its angular bound is
\(12336W(1+\log n)/(n^2\sqrt L)\), and \(L=\ell\) here.

For \(n\ge n_0\), use \(2c\le W_*\).
The common-symbol pair in this sequence interlaces and has
angular gap at least
\((g/2)\sqrt{W_\ell/2}\,n^{-5/3}\), by (1) and the
same-grid perturbation bound. This lower bound tends to
\((g/2)\sqrt c\,n^{-5/3}>0\), so its interlacing remains
strict in the Gaussian limit. For the finitely many smaller
\(n\), the Gaussian Taylor expansion is the case
\(\lambda_r=0\) of (2). The same explicit Hermite disk
comparison gives a positive finite-degree cap on \(c\).
The minimum of these caps and \(W_*/2\) proves the assertion.
