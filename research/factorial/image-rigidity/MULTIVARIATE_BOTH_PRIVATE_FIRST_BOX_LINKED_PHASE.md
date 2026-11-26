# First angle box separation for linked both private moments

Consider an ordinary two-monomial profile with at least one
private variable in each orientation and arbitrary shared rows.
Let \(Q_m\) be its degree-\(m\) moment polynomial and let
\(\Delta=\sum_i(b_i-a_i)\) be the total signed degree difference.
At \(N=n+1\), compare \(Q_n(z)\) with \(Q_N((n/N)^\Delta z)\).

**Theorem.** These polynomials have no common root whose argument
lies in the first positive angle box of both orders, or in its
complex conjugate, for any \(n\ge2\). In particular the ordinary
physical polynomials \(Q_n(z),Q_N(z)\) have this property whenever
\(\Delta=0\). The bound has no curvature or degree cutoff, and
the private exponents can be arbitrary positive integers.
For a common root, the general angle boxes already force the same
positive label at both orders, so the remaining possible collisions
for the compensated pair have label at least two.
[General private-unit angle boxes](MULTIVARIATE_PRIVATE_UNIT_ASYMMETRIC_ANGLE_BOXES.md)

The compensation uses an exact signed matched tanh statistic.
Each unequal row contributes one oriented Gauss slot per unit of
its gap. The difference of their counts in the two orientations
is \(\Delta\). The proof below uses the actual linked coefficient
quotient in both orientations.
[Joint ordinary Gauss grid](MULTIVARIATE_MATCHED_TANH_INTEGER_MULTIPLICITY.md)

## Coefficient quotients at the same test radius

Reserve one private row in each orientation. The Gauss formula
extracts one unit factorial from each reserved row, so the two
factorials cancel the binomial coefficient up to the scalar \(m!\).
The remaining factors of each reserved exponent \(d\) have
unit gaps in the same orientation and Gamma offsets \(h/d\),
\(1\le h<d\). The sum of the signed gaps is still \(\Delta\),
since one gap in each orientation was canceled. The additional
positive exponential scales are independent of the order and
their \(k\)-th powers cancel in linked coefficient ratios.

Divide each polynomial by its constant coefficient. In one
orientation the resulting coefficient \(A_{m,k}\) is, up to
these same scales, a product of
\(\Gamma(ma_i+d_i k+\gamma_i)/\Gamma(ma_i+\gamma_i)\).
Here unreserved rows have \(\gamma_i=1\), and for the residual
private factors \((a_i,d_i)=(0,1)\) or \((1,-1)\)
with the positive Gauss offset above.
At fixed integer \(k\le n\), real-order differentiation gives
\(\partial_m\log A_{m,k}\le\Delta k/m\).
For an unreserved increasing row, subtract the derivative of
each factor \(ma_i+h\), \(1\le h\le d_i k\), from one share
\(1/m\). The loss is \(h/[m(ma_i+h)]\ge0\).
For an unreserved decreasing row the factors are reciprocals
of \(ma_i-h\), \(0\le h<|d_i|k\), and the loss from the share
\(-1/m\) is \(h/[m(ma_i-h)]\ge0\).
An increasing residual factor is independent of \(m\).
For a decreasing residual factor the denominators are
\(m+\gamma_i-j\), \(1\le j\le k\), all positive and at most
\(m\), giving derivative at most \(-k/m\).
Thus \(A_{N,k}(n/N)^{\Delta k}\le A_{n,k}\).

The actual linked ratio has a stronger order. For each factor
with \(a_i>0\), take independent
\(X_i\sim\mathrm{Beta}(na_i+\gamma_i,a_i)\). A factor with
\(a_i=0\) contributes ratio one before the compensation.
The beta integral gives
\[
        \frac{A_{n,k}}{A_{N,k}(n/N)^{\Delta k}}
                =\mathbb E\left[
                   \left((N/n)^\Delta\prod_{i:a_i>0}X_i^{d_i}
                   \right)^k\right] .
\]
All moments used here are finite since the shifted gamma arguments
are positive. The moments are log-convex in \(k\), start at one,
and are at least one at \(k=1\). The ratios are therefore
nondecreasing. In the reversed orientation the compensation
reciprocates and \(\Delta\) changes sign, so the same conclusion
holds after reversal and normalization by the leading coefficient.
Evaluation at a fixed positive radius leaves both coefficient
ratio comparisons intact.
[Positive Mellin coefficient quotients](MULTIVARIATE_BOTH_PRIVATE_LINKED_COEFFICIENT_ORDER.md#positive-mellin-representation-for-the-quotient)

## A first-coefficient radius gate

Write \(A_{m,1}^+\) and \(A_{m,1}^-\) for the first coefficients
of the unscaled moment polynomials after separate constant
normalization in the two orientations. The conclusion also
excludes roots of the original unscaled polynomials at
\(r_n e^{i\theta}\) and \(r_N e^{i\theta}\) in both first boxes
whenever
\[
                 \frac{A_{N,1}^-}{A_{n,1}^-}
                \ \le\ \frac{r_N}{r_n}
                \ \le\ \frac{A_{n,1}^+}{A_{N,1}^+}.           \tag{1}
\]
Indeed the linked quotient at these two radii is the same
positive moment ratio with the deterministic scale \(r_n/r_N\)
in the forward orientation and \(r_N/r_n\) in reverse.
Condition (1) makes each first moment at least one, giving
the same nondecreasing quotient in both orientations.
The compensated case has \(r_N/r_n=(n/N)^\Delta\), which is
in the interval by the derivative inequality. For arbitrary
\(\Delta\), the physical same-argument exclusion applies if
the interval contains one.

## Phase comparison except for the central terms

Suppose there is a common angle \(\theta\) in the first positive
angle box of both orders, with radii satisfying (1). Then
\(\pi/n<\theta<2\pi/(n+1)\).
The case consisting of a lone private unit pair and equal
rows gives consecutive geometric sums and is immediate.

Write \(a_k=A_{n,k}r_n^k\) at the lower order. These weights
are positive and log-convex. Split the lower polynomial
into the left and right halves with half weight on the
center if \(n\) is even. Use the same splitting rule for
the upper polynomial at its proposed radius. In either order \(m\)
write \(\phi_{m,L},\phi_{m,R}\) for the arguments of the two
half sums with nonnegative angular exponents, using the
reversed order for the right half. Every half lies in an
angular sector of width less than \(\pi\). The original
root equation says
\[
                    m\theta=\pi+\phi_{m,L}+\phi_{m,R}.          \tag{2}
\]
Indeed \(m\theta-\phi_{m,L}-\phi_{m,R}\) is positive and
strictly less than \(2\pi\), and the half sums at a root
have opposite phases after restoring the right half.

Augment the lower half sums to the upper index ranges and
weights, using the corresponding lower coefficients for each
additional term. If \(n=2p\), this adds the other half of
\(a_p\) to each half. If \(n=2p+1\), add half the coefficient
at index \(p+1\) in each orientation of the lower array.
The ratio of each upper coefficient to the corresponding
augmented lower coefficient is nonincreasing with the index,
after removal of an overall positive scalar. This follows
from the quotient comparison above in each orientation.

In a sector of width less than \(\pi\), such coefficient
ratios can only decrease the phase. In fact, for
\(F=\sum_k b_k e^{ik\theta}\) and
\(G=\sum_k u_k b_k e^{ik\theta}\), \(b_k>0\),
\[
 \operatorname{Im}(G\overline F)=
     \sum_{j>k} b_j b_k(u_j-u_k)\sin((j-k)\theta)\le0
\]
when \(u_k\) is nonincreasing and the index range has angular
width less than \(\pi\). The two arguments lie in that
sector, fixing the sign of their difference. Thus we only
need to bound the phase increases caused by the lower
central additions.

## The two central additions

Rotate both original physical lower halves by \(e^{-in\theta/2}\).
At a lower root they now equal \(x-iy\) and \(-x+iy\)
with \(y>0\). The common \(y\) is the sine sum of either
half, with positive angles measured from the center.

If \(n=2p\), set \(S=\sum_{h=1}^p\sin(h\theta)\).
Because \(a_{p-h}a_{p+h}\ge a_p^2\) by log convexity,
Cauchy–Schwarz in the two equal sine sums gives
\(y\ge a_p S\). Both added terms have modulus \(a_p/2\).
Writing \(u=\theta/2\), one has
\[
 2\sin u\, S=\cos u-\cos((p+1/2)\theta)
                  >\cos u+\sin u>1 .
\]
Here \(p\theta>\pi/2\) and \((p+1/2)\theta<\pi\).

If \(n=2p+1\), set \(S=\sum_{h=0}^p\sin((h+1/2)\theta)\).
If \(a_{p+1}\ge a_p\), log convexity shows that each
coefficient in the right half is at least \(a_{p+1}\).
If \(a_p\ge a_{p+1}\), it shows that each coefficient in
the left half is at least \(a_p\). The common sine sum
therefore gives \(y\ge\max(a_p,a_{p+1})S\).
The moduli of both central additions are at most
\(\max(a_p,a_{p+1})/2\). In this case
\[
     2\sin u\, S=1-\cos((p+1)\theta)>1 ,
\]
since \(\pi/2<(p+1)\theta<\pi\).

In either parity the ratio of the modulus of each addition
to that of its original half is at most \(1/(2S)<\sin u\).
The disk about the original half with this relative radius
subtends an angle at most \(\arcsin(1/(2S))<u\) on either
side of its argument. The absolute phase change in each half
is therefore less than \(u\). Rotations and conjugation used
to describe the right half do not affect this bound.
Combining with the monotone coefficient comparison gives
\[
       \phi_{N,L}+\phi_{N,R}-\phi_{n,L}-\phi_{n,R}<2u=\theta .
\]
Subtracting the two required equations (2) contradicts this
bound. Conjugation gives the claim in the lower half-plane.

For comparable shared rows the
[conditional monotonicity theorem](MULTIVARIATE_BOTH_PRIVATE_COMPARABLE_SHARED_INTERLACING.md)
already excludes all common physical roots without an assumption
on \(\Delta\). The first-box comparison here applies to arbitrary
mixed orientations and the radius gate (1). When
\(\Delta=0\), it uses the original physical argument.

Since a common label must satisfy \(4k<n+2\), the
first-box theorem excludes all common roots of the compensated
pairs through terminal seven. The first and second moments
have no common root by even-order positivity.
