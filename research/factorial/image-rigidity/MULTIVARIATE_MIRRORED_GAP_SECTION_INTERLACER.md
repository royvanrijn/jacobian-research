# An exact section interlacer for the mirrored gap-two pair

Let \(P_{n,d,\ell}\) be the order-\(n\) factorial polynomial for
the shared mirrored pair \((\ell,\ell+d),(\ell+d,\ell)\).
For \(d=2\) there is an explicit all-degree comparator for the
actual upper polynomial. It comes from padding the *underlying*
degree-\(2n\) Gegenbauer polynomial by one conjugate pair and
then taking its exact weighted even section. The padded section
strictly interlaces \(P_{n,2,\ell}\), and gives an exact sign
benchmark for the unpadded upper section.

## The weighted section from polarization

For \(A=1+n\ell\), define
\[
 G_{dn,A}(w)=\sum_{j=0}^{dn}\binom{dn}{j}
                         (A)_j(A)_{dn-j}w^j,\qquad
 {\cal S}_{n,d}\left(\sum_{j=0}^{dn} c_jw^j\right)
           =\sum_{k=0}^{n}\frac{\binom nk}{\binom{dn}{dk}}c_{dk} z^k .
\]
Then \(P_{n,d,\ell}=(n\ell)!^2{\cal S}_{n,d}G_{dn,A}\).
Equivalently, the section is the symmetric multi-affine
polarization of \(G_{dn,A}\) evaluated on \(n\) copies of
one \(d\)-th root-of-unity orbit. If its elements are
\(\omega^h w\), the elementary symmetric weights follow
from
\(\prod_{h=0}^{d-1}(1+\omega^h w t)=1+(-1)^{d-1}w^dt^d\).
This gives the section variable \(z=(-1)^{d-1}w^d\).
The circle-root property also follows by convolving \(G_{dn,A}\)
with \((1+w^d)^n\) in degree \(dn\) and reading the result
as a polynomial in \(w^d\). The two circle factors are the
Gegenbauer factor and the indicated binomial power.
[Mirrored roots and circular composition](MULTIVARIATE_MIRRORED_ROOT_CIRCLE.md#mirrored-profiles-at-every-fixed-order)

## Padding by a conjugate pair

Fix \(d=2\), \(F={\cal S}_{n,2}G_{2n,A}\), and write
\(E=z\,d/dz\). The exact output of the padded lower input is
\[
 \begin{aligned}
 \widehat F
   &={\cal S}_{n+1,2}\big((1+w^2)G_{2n,A}(w)\big)\\
   &=\frac{z(2E+1)F+(2n+1-2E)F}{2n+1}.
 \end{aligned}                                                   \tag{1}
\]
Indeed the coefficient indexed \(2k\) in the padded input is
\(c_{2k}+c_{2k-2}\). For the first summand, the ratio of the
new section weight to the old one at the same index is
\((2n+1-2k)/(2n+1)\). For the second it is
\((2k-1)/(2n+1)\) relative to the old weight at \(k-1\).

The operator in (1) is a strict degree-raising interlacer
on a broader class. Take any real reciprocal polynomial \(F\)
of degree \(n\) with simple unit-circle roots and \(F(1)>0\).
Put \(R=(2E+1)F\) and \(B=(2n+1-2E)F\). Then \(B=R^*\)
under reversal in degree \(n\). All roots of \(R\) lie
strictly in the disk. For \(|z|\ge1\) away from the roots of
\(F\), this follows from
\(\Re(zF'(z)/F(z))\ge n/2\) by the unit-circle factorization.
At any root of \(F\), simplicity gives \(R(z)=2zF'(z)\ne0\).
Thus \(R/B\) is a finite Blaschke product up to a unit constant.
It has increasing lifted boundary argument starting at zero,
since \(R(1)=B(1)=(n+1)F(1)>0\).
The equation \(R/B=-1\) is exactly \(F=0\).
Meanwhile \(\widehat F=0\) on the circle is the equation
\(zR/B=-1\). Between consecutive positive-cut arguments
\(\theta_i,\theta_{i+1}\) of \(F\), the lifted argument of
\(zR/B\) advances from \((2i-1)\pi+\theta_i\) to
\((2i+1)\pi+\theta_{i+1}\). It crosses exactly one odd
multiple of \(\pi\). The preceding and final intervals
at the cut give the two remaining roots. Thus \(\widehat F\)
strictly interlaces \(F\).

There is also a direct sign benchmark. Write
\(f(\theta)=e^{-in\theta/2}F(e^{i\theta})\), which is real.
At any zero \(\theta_i\) of \(f\), (1) gives
\[
 e^{-i(n+1)\theta_i/2}\widehat F(e^{i\theta_i})
            =\frac{4\sin(\theta_i/2)}{2n+1} f'(\theta_i).          \tag{2}
\]
Consequently, for the exact upper section
\({\cal S}_{n+1,2}G_{2n+2,A+\ell}\), the missing comparison
is a sign agreement with (2) at the lower zeros.
Its centered value is positive at \(\theta=0\), so these
signs would force one zero in each strict interlacing
interval of the lower polynomial, including the two cut
intervals. The padded identity provides the benchmark.
A principal-block comparison of the underlying Gegenbauer
matrices gives inequalities at the underlying degree
difference two, leaving this comparison of the sections open.

The section map itself has no general monotonicity under a
symmetric paired angle shift. For example
\[
 {\cal S}_{2,2}\big((w^2-2cw+1)(w^2-2ew+1)\big)
            =1+\frac{2+4ce}{3}z+z^2 .
\]
Vary \(c\) in \((-3/5,-2/5)\). The resulting positive-half
root angle increases with \(c\) for \(e=1/5\), and decreases
with \(c\) for \(e=-1/5\), by the cosine formula for the
quadratic. All inputs in these two families have positive
coefficients and simple circular roots. Thus a comparison
based on monotonicity in just the underlying paired angles
would need an additional Gegenbauer constraint.
