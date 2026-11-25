# Shared terminal seven at curvature at least 17/4

Let the unequal rows of two distinct ordinary exponent columns
be fully shared, with lower exponent \(a_i\ge1\), gap \(d_i>0\)
and either orientation. Put \(S=\sum_i d_i^2/(a_i+d_i)\).
If \(S\ge17/4\), the normalized order-six and order-seven factorial
moments have disjoint root arguments. Hence the original moments
are coprime over every characteristic-zero field.

Use the [high-curvature certificate](MULTIVARIATE_TERMINAL_SEVEN_HIGH_CURVATURE.md)
when \(S\ge73/10\). For \(17/4\le S<73/10\), the
[exact checker](verify_multivariate_terminal_seven_shared_high_curvature.py)
uses a common order-six outer geometric mean \(u\), the minimum
total curvature, and the discrete possibilities for a row with a
large Gauss reciprocal. Its ray calculations use rational arithmetic.

## A row restriction from the total outer loss

At any lower order \(n\ge2\), put
\(s_i=d_i^2/(a_i+d_i)\). The outer mean on row \(i\) satisfies
\[
 v_i^2=\prod_{h=1}^{d_i}
             \frac{na_i+h}{na_i+(n-1)d_i+h}.
\]
The inequality \(\log(1+y)\ge y/(1+y)\) gives
\(-\log v_i\ge(n-1)s_i/(2n)\), since every denominator
\(na_i+(n-1)d_i+h\) is at most \(n(a_i+d_i)\).
If \(S\ge S_0\), it follows that
\[
 u\le v_i e^{-(n-1)\max(S_0-s_i,0)/(2n)} .
 \tag{1}
\]
This improves \(u\le v_i\) whenever the other rows must supply
positive curvature.

There are finitely many rows with \(d_i/(n a_i+1)>t_0\) and
\(s_i\le B\). Indeed \(a_i<d_i/(nt_0)-1/n\) and
\(d_i<B(1+1/(nt_0))\). The checker takes \(n=6,\ t_0=1/6\),
and \(B=73/10\), enumerates this rectangle and applies the exact
predicates. For each candidate it bounds the right side of (1),
squared, by dividing the exact \(v_i^2\) by the positive Taylor
lower bound for the corresponding exponential. If this upper
bound is below the current interval's lower endpoint squared,
the row cannot appear. The remaining rows and the floor \(t_0\)
give a reciprocal cap at order six. At order seven use the
maximum of the remaining candidates' reciprocals and
\(6t_0/(7-t_0)\) for rows at or below the floor.
The total outer loss also gives \(u\le e^{-5S_0/12}<.171\).

## Joint coefficient bounds

On a shared Gauss factor with lower-order reciprocal \(x\), the
upper-order reciprocal is \(x'=n x/(n+1-\gamma x)\) with
\(0<\gamma\le1\). Consequently
\[
 1+\frac{n^2}{n+1}x\le1+n x'
       \le\frac{1+(n-1)x}{1-x/(n+1)}.
\]
If \(0<x\le t\), the ratio between the two outer logarithms lies
between
\[
 \frac{\log(1+n^2 t/(n+1))}{\log(1+(n-1)t)}
 \quad\text{and}\quad
 1-\frac{\log(1-t/(n+1))}{\log(1+(n-1)t)} .
 \tag{2}
\]
For the lower comparison,
\(\log(1+c(e^z-1))/z\) decreases for \(c\ge1\) by concavity
of its numerator. For the upper comparison
\(-\log(1-x/(n+1))/\log(1+(n-1)x)\) increases since its
numerator is convex, the denominator is concave and both vanish
at zero. Summing across the factors gives (2) for the whole profile.

Let \(j=6,7\), let \(x\) be a reciprocal at that order, and
write \(g_m=\log(1+mx)\) and \(s=\tfrac12 g_{j-1}\).
For one increasing factor the normalized log multiplier at \(k\) is
\[
 f_k=\sum_{m<k}g_m-\frac{k}{j}\sum_{m<j}g_m .
\]
For \(k\le j/2\) the mean loss of the reflected pair is
\[
 -\tfrac12(f_k+f_{j-k})
   =\tfrac12\sum_{m=0}^{k-1}(g_{j-1-m}-g_m)
   =\sum_{m=0}^{k-1}\operatorname{artanh}(c_m\tanh s),
 \qquad c_m=\frac{j-1-2m}{j-1}.
\]
For \(0\le c\le1\) the summand is increasing and concave in \(s\).
Its derivative \(c(1-\tanh^2 s)/(1-c^2\tanh^2 s)\) decreases,
and its derivative at zero is \(c\). Dividing by \(s\) gives a
lower ratio at the reciprocal cap and the upper ratio
\(\sum_m c_m=k(j-k)/(j-1)\).

The skew is bounded in absolute value by \((f_{j-k}-f_k)/2\),
whose log weights are
\[
     \tfrac12\left(\mathbf 1_{\{k\le m<j-k\}}-\frac{j-2k}{j}\right).
\]
For the noncentral pairs at orders six and seven its second
derivative in \(s\) is nonnegative on \(0\le x\le2\). The
checker forms the derivative from these weights with
\(d/ds=2(1+(j-1)x)d/dx/(j-1)\), substitutes \(x=2z/(1+z)\),
and checks the numerator and denominator coefficient signs.
The denominator has a positive constant. The value and first
derivative at zero vanish. The ratio of this skew magnitude to
\(s\) increases, so its value at the reciprocal cap bounds the
aggregate signed skew in either orientation. The middle coefficient
at \(j=6\) has skew zero. Rational logarithm intervals evaluate
these ratios and (2), enclosing each multiplier in powers of \(u\).

Intersect the envelopes with the
[total curvature cap](MULTIVARIATE_FIXED_TERMINAL_CURVATURE.md#large-curvature-and-the-endpoint-radii)
and the [fixed-distance endpoint attenuation](MULTIVARIATE_ORDINARY_ADJACENT_COEFFICIENT_ATTENUATION.md#an-all-order-endpoint-inequality).
At distance \(k\) from an endpoint the consecutive attenuation
factor is at most \((90/91)^{S_0 k^2}\) before the binomial coefficient.
The checker uses the larger \(e^{-S_0 k^2/91}\), by
\(\log(91/90)\ge1/91\), and applies the bound from both ends.

## Rays and sector counts

The two Gaussian anchor arrays, after removing binomial coefficients,
use the multipliers
\[
 u^{k(6-k)/5},\qquad u^{(36/35)k(7-k)/6}.
\]
Both lie within the envelopes. In fact \(u\le e^{-5S_0/12}\)
gives the required total-curvature caps, and the order-seven
exponent exceeds the order-six exponent at the same endpoint distance
\(k\) by \(k^2/35\). The relative multiplier is at most
\(e^{-S_0 k^2/84}\), which fits the broader attenuation cap
used by the checker. The right-end comparison follows by reflection.
Thus the coefficient boxes on each interval of \(u\) contain the
actual arrays and the Gaussian arrays whenever \(u\le e^{-5S_0/12}\),
and their straight segments. The curvature caps put all roots
between radii \(1/8\) and \(8\).

The checker uses five piecewise affine ray cosines as functions of
\(u\). It tests each ray for both orders by
[rational Bernstein projections](MULTIVARIATE_TERMINAL_SIX_JOINT_RAY_BANDS.md#moving-rays-and-exact-replay),
bisecting the intervals of \(u\) or radius where necessary. It also
checks that the rays stay ordered. At \(u=0\) the anchors are
\(1+t^6\) and \(1+t^7\). Cubic sign tests in the cosine of the
argument place their upper roots in alternating sectors, starting
with order seven. Continuous variation through the Gaussian arrays
preserves these sector counts, as does each segment to the actual
array. The last joined sector across the negative real axis has a
single septic root and two sextic roots. The septic root there
is real by conjugacy. The even factorial moment is positive on
the real axis, so the sextic roots there are nonreal. This gives
disjoint arguments, and the rational Bézout identity transfers
coprimality to any characteristic-zero field.
