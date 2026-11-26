# Exact radial separation for ordered opposite factorial pairs

Start with any ordinary mirrored core, allowing private pairs
and repeated pairs. Add a finite collection of opposed shared
pairs. In each pair the two integer gaps have the same absolute
value \(d\ge1\), and the increasing row has midpoint \(m\)
while the decreasing row has midpoint \(M\). Assume
\[
             m\le M\quad\text{in each pair},\qquad
             m<M\quad\text{in at least one pair}.
\]
The gaps and both midpoints are otherwise independent across
the pairs. One may reverse the inequality throughout by
reversing the monomials. Put
\[
                       \sigma=\max_{\text{pairs}} d^2/m .
\]

**Theorem.** There is a computable absolute \(\sigma_0>0\)
such that \(\sigma\le\sigma_0\) gives ordinary SFC(2) in every
degree, with physical coefficient ratio. There is no bound on
the number of rows or the total curvature of the core or
of the paired tail. Exact mirrored pairs can be put in the
core and have no effect on \(\sigma\). The proof gives a
radial exclusion even when the signed center mismatch tends
to zero faster than the tail curvature.

Equivalently, for each fixed gap the increasing and decreasing
midpoint lists can be matched in sorted order with one list
coordinatewise no greater than the other, with the same ordering
direction at every gap. The hypothesis of strict inequality
in at least one pair leaves the fully mirrored case outside
this radial criterion.

## An exact skew test for any diffuse tail of a mirrored core

For a shared unequal tail row put \(d_i=b_i-a_i\),
\(m_i=(a_i+b_i)/2\),
\(w_i=d_i^2/m_i\), and \(\eta_i=|d_i|/m_i\).
Write \(W=\sum_i w_i\) and \(\eta=\max_i\eta_i\).
For \(j\ge2\) and \(s\in[0,1]\) define the row's log-factorial
chord per degree by
\[
 \mathcal L_{j,i}(s)=\frac1j
   \left[g(j(a_i+d_i s))-(1-s)g(ja_i)-s g(jb_i)\right],
        \qquad g(x)=\log\Gamma(1+x).
\]
The total tail log multiplier at grid index \(k\) is
\(J_{j,k}=j\sum_i\mathcal L_{j,i}(k/j)\).
Define its odd part and maximum size by
\[
             J^o_{j,k}=\tfrac12(J_{j,k}-J_{j,j-k}),
             \qquad \delta_j=\max_{0\le k\le j}|J^o_{j,k}|.
\]
Let \(G_j>0\) be the exact endpoint argument scale of the
full profile, so \(G_j^j\) is the ratio of its leading
and constant coefficients. Write \(C=\log(G_N/G_n)\) at
adjacent orders \(n,N=n+1\).

The following test also holds for ordinary shared tails
without any pairing. At bounded curvature and sufficiently
small relative gaps, a physical common root requires
\[
                             |C|\le 2B(\delta_n+\delta_N)       \tag{1}
\]
with one effective absolute \(B\) at all relevant orders.
The signed scale \(C\) and the odd multiplier in (1) are
both exact. The even tail perturbation is used only to
secure the simple circle reference.

Here are uniform choices for (1). Let \(J\ge4\) be an
integer at least the
[absolute eventual angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md),
and choose a rational \(0<\eta_*\le1/2\) from the
[absolute diffuse angular theorem](MULTIVARIATE_DIFFUSE_ALL_ORDER_RESULTANT_TRANSFER.md).
At a possible collision the
[absolute curvature bound](MULTIVARIATE_ABSOLUTE_CURVATURE_12.md)
gives total curvature less than \(12\), and thus \(W<24\).
Orders \(N\ge J\) are already excluded by the eventual
theorem. Order two is separated by positivity of the
even moment on the real axis. If all core rows and tail
rows have relative gaps, measured against their larger
exponents, at most \(\eta_*\), the absolute diffuse theorem
excludes the collision. The case left for the argument has
\(3\le N<J\) and at least one mirrored core row of relative
gap greater than \(\eta_*\).

Endpoint-normalize the exact core and multiply its
\(k\)-th coefficient by \(e^{-jW(k/j)(1-k/j)/2}\).
The resulting reference \(F_j\) has simple unit-circle
roots. These are the mirrored Gaussian references of the
[uniform mesh comparison](MULTIVARIATE_MIRRORED_CORE_SIGNED_DIFFUSE_RADIAL_SEPARATION.md#a-mesh-bound-uniform-over-the-core).
Its distinguished nondiffuse mirror pair gives chord
separation at least
\[
                         g_0=\min\left(1,\frac4{J(2+1/\eta_*)}\right).
\]
These bounds use the actual mirrored core in both orders.

For \(\eta\le1/2\), the integral comparison of the
trigamma series gives
\[
 \left|\sum_i\mathcal L_{j,i}''(s)-W\right|\le 2W\eta .
\]
Indeed, at \(x=m_i+d_i(s-1/2)>0\) compare \(j g''(jx)\) to
\(1/x\) with error at most \(1/(jx^2)\), and then compare
\(1/x\) to \(1/m_i\). Use \(1/m_i\le\eta_i\) and
\(x\ge m_i(1-\eta/2)\). Integrating the difference with
zero endpoint values shows that
\(J_{j,k}+jW(k/j)(1-k/j)/2\) has modulus at most \(jW\eta/4\).

Replace \(J_{j,k}\) in the exact endpoint-normalized
polynomial by its even part
\((J_{j,k}+J_{j,j-k})/2\), and call the resulting
reciprocal polynomial \(F_j^{\rm sym}\).
All these normalized stripped log multipliers are
nonpositive by convexity. Therefore
\[
            |F_j^{\rm sym}(z)-F_j(z)|
                          \le (jW\eta/4)(1+|z|)^j,
      \qquad
            |P_j(z)-F_j^{\rm sym}(z)|
                          \le \delta_j(1+|z|)^j,
\]
where \(P_j\) is the actual endpoint-normalized parent.
Also \(\delta_j\le jW\eta/4\) because the Gaussian
multiplier is even.

Set
\[
         B=J3^J(4/g_0)^{J-1},
         \qquad \eta_1=\min\left(\eta_*,\frac{g_0}{192JB}\right).
\]
For \(\eta\le\eta_1\) the radius \(r=B W\eta\) is at most
\(g_0/8\). Rouché on both radii \(r,2r\)
around each root of \(F_j\) places exactly one root of
\(F_j^{\rm sym}\) in each disk at those radii.
On either circle the monic product is bounded below by
the radius times \((g_0/2)^{j-1}\), which dominates the
displayed error by the choice of \(B\).
Reciprocal conjugation maps a root in the smaller disk
to a root in the larger disk of the same center, since
\(r/(1-r)<2r\). Uniqueness proves it lies on the unit
circle. The roots of \(F_j^{\rm sym}\) thus have chord
separation at least \(3g_0/4\).

Rouché around these new roots with radius \(B\delta_j\)
then gives \(|\log|\zeta||\le2B\delta_j\) for every root
\(\zeta\) of \(P_j\). Indeed \(B\delta_j\le g_0/8\),
the same monic product lower bound applies, and
\(-\log(1-u)\le2u\) for \(0\le u\le1/2\).
If \(\delta_j=0\), \(P_j=F_j^{\rm sym}\) already.
At a common physical root the two normalized log radii
differ by \(C\). This proves (1).

## Same-gap pairs give a center larger than their skew

For an increasing row of gap \(d\ge1\) and midpoint \(s\),
write
\[
 H_j(s)=\frac{g(j(s+d/2))-g(j(s-d/2))}{j},
 \qquad E(s)=H_N(s)-H_n(s)-d\log(N/n).
\]
These functions give the row's exact adjacent endpoint
scale. Reversal changes their signs. Suppose \(s\ge2d\).
Differentiating in \(s\) gives
\[
 E'(s)=\int_{-d/2}^{d/2}
         [N g''(N(s+y))-n g''(n(s+y))]\,dy .
\]
For \(x>0\) the bracket at \(x=s+y\) is
\[
          \int_0^\infty e^{-xu}
               [q(u/N)-q(u/n)]\,du,\qquad
                     q(t)=\frac{t}{e^t-1}.
\]
This follows from the same trigamma integral used for the
[exact factorial center](MULTIVARIATE_FIXED_TERMINAL_COMPACTIFICATION.md#uniform-estimates-for-rows-of-large-size).
Since \(q\) is decreasing, the integrand is positive.
For \(0\le t\le1\) one has \(-q'(t)\ge1/18\).
For \(t>0\) the numerator of \(-q'\) is
\(\int_0^t u e^u\,du\ge t^2/2\), and its denominator
is at most \(e^2t^2<9t^2\). Restricting the Laplace
integral to \(0<u<1/x\) gives a lower bound
\(1/(108nN x^2)\). Indeed \(u/n\le1\) there,
the difference of the \(q\)'s is at least \(u/(18nN)\),
and \(e^{-xu}\ge1/3\). As \(x\le5s/4\), this proves
\[
                         E'(s)\ge\frac{d}{200nN s^2}>0 .        \tag{2}
\]

Let \(\mathcal L^o_{j,d,s}(v)\) be the odd part,
under \(v\mapsto1-v\), of the chord per degree
of this one increasing row. Put \(y=v-1/2\).
Its curvature is
\[
  (\mathcal L^o_{j,d,s})''(v)
     =\tfrac12 j d^2
         [g''(j(s+dy))-g''(j(s-dy))].
\]
After differentiating in \(s\), the mean value bound
using \(0<g^{(4)}(x)\le2/x^3\) gives
\[
       |\partial_s(\mathcal L^o_{j,d,s})''(v)|
                              \le \frac{d^3}{(s-d/2)^3}.
\]
The bound on \(g^{(4)}\) follows by differentiating the
trigamma series twice and integrating its positive terms.
Integrating the chord with zero endpoint values costs at
most \(1/8\). Since \(s-d/2\ge3s/4\), it follows that
\[
                 \|\partial_s(j\mathcal L^o_{j,d,s})\|_\infty
                                      \le \frac{j d^3}{3s^3}.   \tag{3}
\]

In one opposed pair the signed center contribution is
\(E(m)-E(M)\) and its odd log contribution at order \(j\)
is \(j\mathcal L^o_{j,d,m}-j\mathcal L^o_{j,d,M}\).
The term \(d\log(N/n)\) cancels and mirroring reverses
the odd log part. Integrate (2) and (3) from \(m\) to \(M\).
If \(m<M\), the center contribution is negative, and
the supremum of the odd contribution is at most
\[
          \frac{200}{3}\,j nN\,\frac{d^2}{m}\bigl[E(M)-E(m)\bigr].
                                                                    \tag{4}
\]
Both contributions vanish when \(m=M\).
Since the centers of the mirror core cancel, all nonzero
pair contributions have the same sign, and \(C<0\).
Adding (4) gives
\(\delta_j\le70 J^3\sigma |C|\) for \(j=n,N<J\).

Choose
\[
                      \sigma_0=\min\left(\eta_1,\frac1{560 B J^3}\right).
\]
Because \(d\ge1\), each row of the paired tail has
\(\eta_i\le d^2/m\le\sigma_0\), and \(m\ge2d\).
Thus (1) and (4) would require
\(|C|\le280 B J^3\sigma |C|\le |C|/2\).
Its strict positivity in modulus gives a contradiction.
For orders outside the comparison range the absolute
angular and curvature theorems or real positivity already
give separation. Rational Bézout extends the result to
characteristic-zero coefficient fields.

The same comparison can use a weighted bound when the
pairs do not all have the same ordering. At one pair of
orders put \(C_\ell=E(m_\ell)-E(M_\ell)\) and
\(\sigma_\ell=d_\ell^2/\min(m_\ell,M_\ell)\).
If \(\eta\le\eta_1\) and \(C=\sum_\ell C_\ell\ne0\),
the tested orders are separated whenever
\[
             \sum_\ell \sigma_\ell |C_\ell|
                          \le \frac{|C|}{560 B J^3}.
\]
Indeed (4) with the smaller midpoint gives the bound
on the sum of the absolute odd contributions even
without the common center sign. This permits some larger
individual curvatures if their exact center contributions
are small. Both sides are computable directly from the
factorial endpoint scales.

For example add the same positive number of copies of
\((A,A+1)\) and \((A+2,A+1)\) to an arbitrary mirrored
core. For sufficiently large \(A\) the theorem gives
all-degree ordinary coprimality uniformly in both the
multiplicity and the core. The total tail gap is zero,
and the ratio of the signed log bias to the tail curvature
is asymptotic to \(1/(2A)\). It therefore tends to zero,
unlike the fixed positive fractional bias in the earlier
[signed diffuse radial criterion](MULTIVARIATE_MIRRORED_CORE_SIGNED_DIFFUSE_RADIAL_SEPARATION.md).
The same proof permits independent gap sizes and midpoint
ratios within the stated pairing and per-pair curvature
bound. It leaves tails with canceling adjacent center
contributions and large odd chord changes outside the
ordered-pair criterion.
