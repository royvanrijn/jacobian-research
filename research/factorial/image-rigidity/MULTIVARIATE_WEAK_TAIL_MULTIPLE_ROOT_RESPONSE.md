# Multiple-root separation under arbitrary weak ordinary tails

Fix an ordinary two-monomial exponent profile, called the core.
Let \(Q_m\) be its factorial moment polynomial in order \(m\).
Append a finite nonempty list of unequal shared exponent rows.
For an appended row write \(d_i=b_i-a_i\),
\(M_i=(a_i+b_i)/2\), and put
\[
 W=\sum_i d_i^2/M_i .
\]
No bound on the number, individual gaps or signed total gap of
the rows is imposed. The results concern \(W>0\) sufficiently small.

At a core root \(\alpha\), write \(r=r_m(\alpha)\) for its
multiplicity and \(\mathcal D=z\,d/dz\). Define
\[
          c_m(\alpha)=\frac{\mathcal D^{r+1}Q_m(\alpha)}
                             {(r+1)\mathcal D^r Q_m(\alpha)} .
\]

**Theorem.** Take adjacent orders \(n,N=n+1\), \(n\ge2\), and
a pair of core roots \(\alpha,\gamma\) with equal argument.
All the nearby augmented roots in these two clusters have
different arguments for sufficiently small \(W>0\) if
one of the following holds.

- At least one of \(r_n(\alpha),r_N(\gamma)\) is even.
- Both multiplicities are odd and
  \(\operatorname{Im}(c_n(\alpha)/n)\ne
    \operatorname{Im}(c_N(\gamma)/N)\).

If both multiplicities are odd, every possible contact for
small \(W\) must use the unique central branch of each cluster.
Their angular changes are
\[
             -\frac{W}{m}\operatorname{Im}c_m(\alpha)
                           +O_{Q_m,m}(W^{3/2}),                 \tag{1}
\]
with \(\gamma\) used at the upper order. The error and the
threshold are uniform over the tail.

The second condition has an exact arithmetic certificate.
Use a common nonarchimedean embedding for the core roots and
a valuation \(v_p\) normalized by \(v_p(p)=1\). Let
\[
 \delta=v_p(\alpha-\overline\alpha)-v_p(\alpha),\qquad
 R=\max\left(0,\
       \left\{v_p(\xi-\zeta)-v_p(\xi):
       (m,\xi)=(n,\alpha),(N,\gamma),\
       Q_m(\zeta)=0,\ \zeta\notin\{\xi,\overline\xi\}\right\}\right).
\]
Set \(e=v_p(2)\) and \(I=N r_n(\alpha)-n r_N(\gamma)\).
If
\[
                        \delta>e,\qquad 2(\delta-R)>e+v_p(I),     \tag{2}
\]
the two central angular changes differ. For an odd prime
not dividing \(I\), (2) says that the conjugate pair at each
order is strictly closer \(p\)-adically than the other roots
are to that pair, relative to its root valuations.
This permits repeated first residual roots and multiple
core roots. For two simple roots \(I=1\).

## Uniform Hermite expansion of an entire cluster

Independently normalize the tail endpoints at order \(m\) by
positive scalars and a positive variable dilation. The new
multiplier \(T_{m,k}\) on coefficient \(k\) of \(Q_m\) has
\(T_{m,0}=T_{m,m}=1\). Its logarithm \(g(k)\) is the sum,
over the appended rows, of the Gamma chords
\[
 \log\Gamma(1+ma_i+kd_i)
 -(1-k/m)\log\Gamma(1+ma_i)-(k/m)\log\Gamma(1+mb_i).
\]
For \(W<1\) and \(0\le k\le m\), the Gamma argument minus one
is at least \(mM_i/2\). Also \(1/M_i\le W\) and
\(|d_i|/M_i\le W\), since every gap has integer magnitude
at least one. The series for the polygammas gives, for \(x>0\),
\[
 0\le x^{-1}-\psi'(1+x)\le x^{-2},\qquad
 |\psi^{(j-1)}(1+x)|\le (j-2)!x^{-(j-1)} \quad(j\ge2).
\]
These estimates follow by comparing the decreasing summands
of the polygamma series to their integrals. Hence uniformly
in the tail,
\[
   g''(k)=W/m+O_m(W^2),\qquad
   \|g^{(j)}\|_{[0,m]}\le C_{m,j}W^{j-1}\quad(j\ge2).
\]
For the second inequality sum
\(|d_i|^j/M_i^{j-1}
\le (d_i^2/M_i)W^{j-2}\).
For the first, replacing \(ma_i+kd_i\) by \(mM_i\) incurs
a bounded relative \(O(W)\) error, and the trigamma error
has the same size. These bounds also give
\(g(k)=W(k^2/m-k)/2+O_m(W^2)\) by the zero endpoint values.

Interpolate the values \(g(0),\ldots,g(m)\) by a polynomial
\(G(k)=\sum_j g_j k^j\) of degree at most \(m\).
The finite difference formula with basis \(\binom{k}{j}\)
and the derivative bounds show that
\[
 \begin{split}
  g_0&=0,\qquad g_1=-W/2+O_m(W^2),\\
  g_2&=W/(2m)+O_m(W^2),\qquad
  g_j=O_m(W^{j-1})\quad (j\ge3).
 \end{split}                                                     \tag{3}
\]
Indeed every order-\(j\) finite difference is the integral
of \(g^{(j)}\) over a \(j\)-dimensional unit cube.
The first two coefficients of the interpolation of the
quadratic approximation are its actual coefficients, and
the interpolation operator on this fixed grid is bounded.

In the local coordinate \(y=\log(z/\alpha)\) write
\(Q_m(\alpha e^y)=b y^r(1+c y+O(y^2))\), with
\(c=c_m(\alpha)\). The augmented, tail-normalized polynomial
in this coordinate is exactly
\(\exp(G(\partial_y))Q_m(\alpha e^y)\): on the finite
span of the exponentials \(e^{ky}\), the operator multiplies
by \(\exp(G(k))=T_{m,k}\).
Set \(t=\sqrt W\), \(y=t u\), and write
\[
     K_{r,m}(u)=\exp(\partial_u^2/(2m))u^r .
\]
The expansion on every fixed compact \(u\)-disk, after
division by \(b t^r\), is
\[
 K_{r,m}(u)
 +t\bigl[cK_{r+1,m}(u)-\tfrac12 K'_{r,m}(u)
                         +\mu_m K'''_{r,m}(u)\bigr]+O(t^2),     \tag{4}
\]
where \(\mu_m=g_3/W^2\) is real and uniformly bounded
and is zero when \(m<3\).

For the error control in (4), truncate the exponential on
its original finite spectrum after order \(K\), with
\(2(K+1)\ge r+2\) and \(2K\ge r+1\). Its residual is
\(O(t^{2K+2})\) before division. In the finite expansion
each monomial operator of derivative order \(j\) has a
coefficient \(O(t^j)\) by (3), so replacing the Taylor
series of \(Q_m(\alpha e^y)\) through degree \(r+1\) costs
\(O(t^{r+2})\) in total. Terms differentiating the remainder
above degree \(r+1\) have the same bound. After scaling,
the second derivative coefficient is \(1/(2m)+O(t^2)\);
the first and third are respectively \(-t/2+O(t^3)\)
and \(t\mu_m\); every higher term is \(O(t^2)\).
The differential operators commute and their expansions
on the degree-\((r+1)\) Taylor polynomial are finite.
This proves (4) with uniform errors, including the errors
needed for root perturbation on smaller disks.

The zeros of \(K_{r,m}\) are
\(i h/\sqrt m\) for \(h\) a zero of the ordinary monic
probabilists' Hermite polynomial \(\mathrm{He}_r\).
They are all simple and imaginary. Rouché's theorem on
disjoint small disks around them gives all \(r\) branches
of the root cluster of the augmented polynomial. For
odd \(r\) the central branch has \(u\to0\) and the
implicit simple-root expansion gives
\[
          u=-t\left(\frac{c}{m}-\frac12+
                            \mu_m m(r-1)\right)+O(t^2).
\]
Here \(K_{r+1,m}(0)/K'_{r,m}(0)=1/m\) and
\(K'''_{r,m}(0)/K'_{r,m}(0)=m(r-1)\), directly from
the exponential series. Taking the imaginary part of
\(y=tu\) proves (1).

## Why the noncentral branches cannot match

The square-variable polynomials obtained from
\(\mathrm{He}_r(x)\) upon removing a factor \(x\) in the
odd case are irreducible over \(\mathbb Q\).
This is the classical Hermite irreducibility of Schur
in the square variable, or follows from the corresponding
assertion for the even polynomial and the odd quotient
in \(x\), followed by the rational rescaling of the square
variable from the classical to the probabilists' convention.
[Schur, section II](https://doi.org/10.1515/crll.1931.165.52)

If two nonzero leading angular changes were equal at
adjacent orders, their squares would be a common root
of these irreducible polynomials scaled by \(1/n\) and
\(1/N\). Their degrees must therefore agree, say both
are \(a\ge1\). Their multiplicities have the form \(2a\)
or \(2a+1\). The mean of the squared nonzero Hermite
roots is \(2a-1\) in the even case and \(2a+1\) in the
odd case, from the coefficient \(-r(r-1)/2\) of
\(x^{r-2}\). These means after division by \(n,N\)
never coincide. For equal parity the numerator is
the same and positive. For distinct parity either the
larger numerator has the smaller denominator, or equality
would require \(2a-1=2n\). This last equation has opposite
parities on its two sides.

A noncentral root therefore cannot match any branch in
the other cluster to first order in \(t\). An even
multiplicity has no central branch. When both
multiplicities are odd, (1) separates the central
branches under the second hypothesis. The omitted
endpoint normalizations are positive and do not affect
arguments. This proves the theorem.

## Local separation at a close conjugate pair

Factor at a root \(\xi\) of multiplicity \(r\) to obtain
\[
        c_m(\xi)=\frac r2+
                      \sum_{\zeta\ne\xi}\frac{\xi}{\xi-\zeta}, \tag{5}
\]
with the other roots counted with multiplicity.
This follows by expanding \((e^y-1)^r\) and the remaining
factor near \(\xi\). Because \(\alpha,\gamma\) have equal
argument, \(\gamma/\alpha=\overline\gamma/\overline\alpha\).
Thus the relative depth \(\delta\) and the ratios
\((\xi+\overline\xi)/(\xi-\overline\xi)\) in the two orders
are algebraically the same.

Consider \(N(c_n(\alpha)-c_n(\overline\alpha))
-n(c_N(\gamma)-c_N(\overline\gamma))\).
The contribution of the repeated conjugate pair in (5)
is \(I(\alpha+\overline\alpha)/(\alpha-\overline\alpha)\).
The integer \(I\) is nonzero, since \(\gcd(n,N)=1\) and the
multiplicities of nonreal roots are strictly below the
polynomial degrees. Under \(\delta>e\) this term has
valuation \(v_p(I)+e-\delta\).

For any remaining root \(\zeta\) at either order, put
\(t_0=v_p(\xi)\), \(u=v_p(\zeta)\) and
\(s=v_p(\xi-\zeta)-t_0\). Since \(s<\delta\), the distance from
\(\overline\xi\) to \(\zeta\) has the same valuation. The
difference of its two fractions in (5) has valuation
\[
                         \delta+(u-t_0)-2s \ \ge\ \delta-2R.
\]
Indeed for \(u\ge t_0\) this follows from \(s\le R\), and
for \(u<t_0\) one has \(s=u-t_0<0\). Multiplication by
\(n\) or \(N\) cannot lower this bound. Condition (2)
makes the conjugate-pair term the unique term of least
valuation, so the displayed difference is nonzero.
In the complex embedding this is exactly nonzero
difference of the imaginary responses.

For a simple quadratic contact in a common rational
chart, a convenient sufficient test at an odd prime is
that the monic shared quadratic has integral coefficients,
unit constant coefficient and discriminant divisible by
\(p\). Require that the two cofactors have integral
coefficients in their respective scalar normalizations
and nonzero reduction at the repeated residual root.
All other roots then satisfy \(R=0\) by integral
polynomial reduction at that unit residue. Here
\(\delta=v_p(\operatorname{disc} H)/2>0\) and \(I=1\).
At \(p=2\) a discriminant valuation at least three
also suffices under the same cofactor test.
For exact algebraic controls, take \(Q_3=H(z)(z+2)\) and
\(Q_4=H(z)(z^2+3z+3)\) with \(H=z^2+2z+1+p\) at an
odd prime, or \(H=z^2+2z+3\) at \(2\). In these controls
the shared conjugates lie on the same valuation face and
have the same residue. The controls have positive rational
coefficients and positive even parent on the real axis.

## All-degree consequence and the bounded-gap reduction

Take \(J\) at least the
[absolute angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
If every equal-argument core pair at \(2\le n<N<J\)
satisfies one of the two conditions, take the minimum
of the finitely many thresholds. Core pairs with unequal
arguments cannot meet under sufficiently small perturbations.
The pair starting at order one is separated by positivity
of the even moment, and the absolute theorem handles the
higher orders. This yields all-degree angular separation
for every sufficiently weak nonempty ordinary shared tail,
independent of its signed degree difference. Thresholds
are effective using polynomial root isolation and the
uniform analytic estimates in (3) and (4) with positive
margins for the finite comparisons.

For original physical coprimality one can use a weaker
hypothesis. For a prescribed signed tail gap \(h\), test
only those pairs which are actual common roots after
twisting the upper core polynomial by \(((n+1)/n)^h\).
Indeed the endpoint ratio of the actual tail differs
from that twist by \(1+O(W)\), uniformly, by the
[core twist comparison](MULTIVARIATE_ARBITRARY_CORE_TWISTED_TAIL_STABILITY.md#proof-of-stability).
The finitely many different limiting physical roots
cannot coincide for sufficiently small \(W\), and the
new argument separates every tested contact by angle.
The possible integer twists for each fixed core are
already effectively finite by the same core comparison.

In particular, this gives a finite branching procedure
for the original bounded-total-gap property whenever the
escaping-core tests succeed. Fix a total gap budget \(L\).
Order the unequal rows by their height. The
[uniform weak-curvature theorem](MULTIVARIATE_UNIFORM_WEAK_CURVATURE_ANGLE_SEPARATION.md)
excludes profiles whose smallest unequal height exceeds a
computable bound, since \(W\le2L^2/\min_i\max(a_i,b_i)\).
Enumerate the first row below that bound. At a retained prefix
\(C\), test its original moments exactly and apply the
physical criterion to each integer \(h\) with
\(|h|\le L-\sum_{i\in C}|d_i|\).
If every test succeeds, take the minimum positive weak-tail
tolerance. When all further rows have height above a
sufficiently large cutoff depending on \(C\), their
\(W\) is below this tolerance. Otherwise branch over the
next row below that cutoff. There are finitely many
children at each such step and the depth is at most \(L\).
A node with an unresolved odd central contact remains
an open branch. Compared to the
[bounded-gap box test](MULTIVARIATE_BOUNDED_GAP_FINITE_DECISION.md#the-bounded-gap-physical-limit),
this can certify a branch even when the limiting affordable
twist itself has an exact common root.

## Stability when the body only approximates a reference

At fixed adjacent degrees the local expansion also applies
to a fixed positive-coefficient reference pair, including
a pair obtained by a Gaussian deformation. The body
polynomials in an actual profile can approach that pair
after separate positive variable scalings and scalar
normalizations. Suppose the reference pair satisfies the
multiple-root separation criterion for all its coincident
arguments. Let \(r_*\) be the maximum multiplicity over
its roots. If the coefficient error of the body is
\(o(W^{(r_*+1)/2})\), an additional shared tail of curvature
\(W\downarrow0\) gives angular separation for sufficiently
small \(W\).

Indeed, multiplying by the normalized tail coefficients
has operator norm at most one on the maximum absolute
coefficient norm, since \(0<T_{m,k}\le1\). On every fixed
scaled \(u\)-disk the body error is still of its original
order before division by \(t^r\) in (4). Afterwards it is
\(o(t)\) for every cluster. Thus the simple-root displacements
on the Hermite scale change by \(o(t)\), by the same
analytic simple-root estimate, so their logarithmic
displacements change by \(o(W)\). Both the noncentral
separation and the central velocity separation persist.
This gives a condition on two scales near a fixed reference.
The comparison concerns profiles with the stated body
approximation relative to the additional curvature.
