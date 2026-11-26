# Root monotonicity under mirrored row changes with one private unit pair

Consider the ordinary two-monomial profile consisting of a private
unit pair \((1,0),(0,1)\) and any finite collection of additional
mirrored pairs \((\ell_i+d_i,\ell_i),(\ell_i,\ell_i+d_i)\) with
\(d_i\ge1\) and \(\ell_i\ge0\). Write \(\theta_{j,k}\in(0,\pi)\)
for the positive root arguments in increasing order at degree \(j\).

**Theorem.** At every degree \(j\ge2\), the roots are simple on
the unit circle. Increasing a lower exponent \(\ell_i\) at fixed
gap shifts every \(\theta_{j,k}\) strictly towards \(\pi\).
Increasing a gap at fixed lower exponent, inflating both exponents
of a mirrored pair by the same integer factor greater than one, or adding a
mirrored pair shifts every \(\theta_{j,k}\) strictly towards zero.
The fixed-order comparisons also apply to real exponents and
arbitrary equal positive Gamma shapes in each additional pair
whenever the indicated parameter comparisons hold. The private
unit pair retains ordinary shape one.

The same comparisons hold after adding the Gaussian symbol
\(-\tau q(1-q)/2\) at every degree. At fixed degree, increasing
\(\tau\ge0\) strictly shifts all positive arguments towards zero.

There is an exact matched-order comparison of the even log-ratio
laws used in the proof. Let \(\rho_j\) be the law of the core
log ratio under the geometric mean moment tilt in degree \(j\).
In likelihood-ratio order of absolute values,
\[
                  (j/(j+1))\,\rho_j\ \preceq\
                          \rho_{j+1}\ \preceq\ \rho_j .
                                                               \tag{1}
\]
Dilation acts on the variable of the law. Relation (1) compares
the two positive integral kernels before the angular clock
changes.

## An even log-ratio law and a two-term root equation

Let \(U,V\) be the monomials from the additional pairs. Integration
over the private unit pair cancels the binomial coefficient.
Thus up to the positive factor \(j!\) the order-\(j\)
polynomial has coefficient at \(k\) equal to
\(\mathbb E[U^{j-k}V^k]\).
Put \(s=\log(V/U)\) and weight the original law by
\((UV)^{j/2}\). After division by the positive total mass the
law of \(s\) is the even probability measure \(\rho_j\).
If there are no additional pairs then \(\rho_j=\delta_0\).

Within a pair the geometric tilt gives both independent Gamma
variables the same shape
\(\alpha_{j,i}=\gamma_i+j(\ell_i+d_i/2)\).
Integrating their ratio gives density proportional to
\[
                    \cosh\bigl(s/(2d_i)\bigr)^{-2\alpha_{j,i}}.
                                                               \tag{2}
\]
For example, the ratio of two independent Gamma variables of
shape \(\alpha\) has density proportional to
\(r^{\alpha-1}/(1+r)^{2\alpha}\) on \(r>0\).
The substitution \(s=d_i\log r\) gives (2).
The full \(\rho_j\) is the convolution of these even densities.

Write \(\phi=\theta/2\), \(c=\cos\theta\), and consider the
real centered circle value. Summing the finite geometric
series inside the integral gives it as a positive scalar
times
\[
               A_j(\theta)\cos(j\phi)
                    -B_j(\theta)\cos((j+2)\phi),              \tag{3}
\]
where, using the law of \(|s|\) under \(\rho_j\),
\[
 A_j(\theta)=\int\frac{\cosh((j/2+1)s)}{\cosh s-c}\,d\rho_j^+(s),
 \qquad
 B_j(\theta)=\int\frac{\cosh(js/2)}{\cosh s-c}\,d\rho_j^+(s).
\]
Here \(\rho_j^+\) denotes the pushforward to the positive
half-line including zero. To verify (3), set \(a=e^{-s}\) for
\(s\ge0\) and sum the two series with coefficients
\(e^{js/2}(a^k+a^{j-k})\). The real numerator is proportional
to \((1+a^{j+2})\cos(j\phi)-(a+a^{j+1})\cos((j+2)\phi)\).
The denominator is \(1-2ac+a^2=2a(\cosh s-c)\), which yields
(3). The moments in this calculation are finite.

For a nonempty core \(A_j>B_j>0\). At
\(\theta=(2k-1)\pi/j\), (3) has sign \((-1)^{k+1}\).
At \(\theta=2k\pi/(j+1)\) the two cosines agree and
the sign is \((-1)^k\). Thus each interval
\[
          \left((2k-1)\pi/j,\ 2k\pi/(j+1)\right),
                    \qquad 1\le k\le\lfloor j/2\rfloor,
                                                               \tag{4}
\]
has a root. These intervals and their reflected intervals
contain the entire degree count, including the structural
root \(-1\) in odd degree. This proves location and simplicity.
For an empty core the polynomial is geometric and its
positive arguments equal the right endpoints of (4).

## The folded likelihood-ratio comparison

For laws with positive densities on the positive half-line,
write \(\mu\preceq\nu\) when the ratio of their densities
\(d\nu/d\mu\) is nondecreasing. This order increases the
expectation of each increasing function whenever the
expectations exist. It also survives a common positive
reweighting.

For \(f(x)=\cosh(c_0x)^{-\alpha}\), \(c_0,\alpha>0\),
consider the folded kernel
\(K_f(x,y)=f(x-y)+f(x+y)\) for \(x,y\ge0\).
It has positive order-two minors at strictly ordered
arguments. To see this set \(p=\tanh(c_0x)\) and
\(q=\tanh(c_0y)\). Removing positive factors in each
variable gives
\[
 (1-pq)^{-\alpha}+(1+pq)^{-\alpha}
            =2\sum_{r\ge0}\frac{(\alpha)_{2r}}{(2r)!}
                                       p^{2r}q^{2r}.
\]
Expanding the two-by-two determinant of a finite partial
sum gives only nonnegative terms, from pairs of equal
or increasing degrees. The pair of degrees zero and
two contributes a strictly positive term when the two
arguments of each variable are distinct. Passing to
the convergent series proves the claim.

Folded kernels compose under convolution of even densities.
The two-by-two determinant formula for an integral kernel
composition expresses its minors as integrals over ordered
pairs of products of minors of the two kernels. This follows
by expanding the two products and swapping the integration
coordinates on half of their domain. Thus convolutions of
the densities (2) have the same strict order-two property.
The same identity with one kernel column replaced by each
of two densities proves that convolution by such an even
density preserves the likelihood-ratio order of the folded
laws: the integrand is a nonnegative kernel minor times the
nonnegative density cross product from the likelihood ratio.

At fixed gap, increasing the lower exponent increases
\(2\alpha_{j,i}\), so the new-to-old density ratio is strictly
decreasing in \(|s|\). At fixed lower exponent, increasing
the gap decreases both the coefficient and the argument in
the magnitude of the log derivative
\[
 \frac{\alpha_{j,i}}{d_i}
          \tanh\!\left(\frac{s}{2d_i}\right)
       =\left[\frac j2+\frac{j\ell_i+\gamma_i}{d_i}\right]
          \tanh\!\left(\frac{s}{2d_i}\right),\qquad s>0 .
\]
Hence the new-to-old ratio is strictly increasing in \(|s|\).
Inflating both exponents together also decreases that log
derivative magnitude. Increasing the common Gamma shape
has the same order effect as increasing overlap.
Replacing factors successively and applying the folded
convolution result gives each comparison for the full law.

For addition of a pair to a nonempty core, the folded kernel
of the old law satisfies
\(K(x,y)/K(x,0)\) strictly increasing in \(x\) for \(y>0\).
Integrating against the new folded law proves the increased
likelihood-ratio order relative to the old folded density.

For (1), applying the log derivative formula to each factor
shows that the order-\((j+1)\) density is more concentrated
than the order-\(j\) density. Dilating the latter factor by
\(c_0=j/(j+1)\) instead changes its log derivative magnitude
to
\(\alpha_{j,i}/(c_0d_i)\tanh(s/(2c_0d_i))\).
Both factors here are greater than their counterparts at
order \(j+1\), since
\(\alpha_{j,i}/c_0=\alpha_{j+1,i}+\gamma_i/j\).
Thus the dilated order-\(j\) factor is more concentrated
than the order-\((j+1)\) factor. Dilation commutes with
convolution, so the two comparisons prove (1).

## From likelihood ratios to movement of every root

Let \(R_j=A_j/B_j\). It is the weighted expectation of
\(\cosh((j/2+1)s)/\cosh(js/2)\), a strictly increasing
function of \(s>0\). The weight
\(\cosh(js/2)/(\cosh s-\cos\theta)\) is positive.
Thus comparisons of the folded laws transfer to \(R_j\)
at every fixed \(\theta\). Also \(R_j\) increases with
\(\theta\in(0,\pi)\): for \(\theta_2>\theta_1\), the ratio
of those weights is
\[
                  \frac{\cosh s-\cos\theta_1}
                       {\cosh s-\cos\theta_2},
\]
strictly increasing in \(s>0\).

Inside each interval (4) the root equation is
\[
            R_j(\theta)=g_j(\theta)
                   :=\cos\theta-\tan(j\theta/2)\sin\theta .
\]
The right side strictly decreases there. Indeed let
\(y=-\tan(j\theta/2)>0\). Its derivative is
\(-[1+(j/2)(1+y^2)]\sin\theta+y\cos\theta\).
This is negative when \(\cos\theta\le0\). Otherwise use
\((1+y^2)/2\ge y\) and \(j\tan\theta>\pi\), since
\(\theta>\pi/j\).
Increasing \(R_j\) pointwise therefore decreases the unique
root argument in each interval. Decreasing \(R_j\) increases
that argument. The strict orders already proved give the
claimed row operations. The empty-core comparison for
addition follows from (4).

## Adding the Gaussian symbol

Multiplying the degree-\(j\) coefficients by
\(\exp[-\tau k(j-k)/(2j)]\) is, up to a common scalar, exactly
convolution of \(\rho_j\) with the centered real Gaussian of
variance \(\tau/j\). Indeed its moment generating function
at \(k-j/2\) is
\(\exp[\tau(k-j/2)^2/(2j)]\).
The folded Gaussian kernel, after removal of the separate
positive factors in its two variables, is
\(2\cosh(xy/(\tau/j))\). Its even series gives the same strict
order-two determinant property as above.

At fixed degree, increasing the variance adds an even Gaussian
convolution and therefore increases the folded likelihood-ratio
law. This yields the asserted shift towards zero in (4).
The sandwich (1) still holds for the Gaussian-augmented laws.
The variance at degree \(j+1\) is \(j/(j+1)\) times the variance
at degree \(j\), between its undilated value and its value
after dilation by \(j/(j+1)\). The Gaussian density ratios
are ordered on the positive half-line, and the folded
convolution preserves those comparisons.

The consecutive angular comparison also changes the two
cosines in (3) and the weighting in \(R_j\). The
[linked phase comparison](MULTIVARIATE_ONE_PRIVATE_MIRRORED_LINKED_PHASE_INTERLACING.md)
uses the lower half of the likelihood-ratio sandwich to
control both changes and prove strict interlacing.
