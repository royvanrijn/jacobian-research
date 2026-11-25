# Matched shared rays from outer mean .48 to .76 at terminal nine

Let all unequal exponent rows of two distinct ordinary exponent columns
be shared, and write \(u\) for the normalized order-eight outer reflected
mean. The eighth and ninth factorial moment polynomials have disjoint
root arguments under independent positive rescalings under any
of the following conditions.

| Outer mean \(u\) | Upper bound on all \(d/(8\ell+1)\) |
| --- | ---: |
| \([.48,.55]\) | \(1/8\) |
| \([.55,.62]\) | \(1/10\) |
| \([.62,.70]\) | \(3/40\) |
| \([.70,.76]\) | \(11/200\) |

In particular the moment polynomials are coprime over every
characteristic-zero field. This extends the
[preceding matched-ray bands](MULTIVARIATE_TERMINAL_NINE_SHARED_MATCHED_EXTENSION.md).

## Ordinary reciprocal bounds

An individual row with gap \(d\) and lower base \(\ell\)
has squared outer mean
\[
 u_{\ell,d}^2=\prod_{h=1}^d\frac{8\ell+h}{8\ell+7d+h}
              \le\left(\frac{8\ell+d}{8\ell+8d}\right)^d.
\]
The full mean is the product of the row means, each at most one.
Thus \(u\ge u_-\) implies \(u_{\ell,d}^2\ge u_-^2\) for every row.
If \(d/(8\ell+1)>t\), then \(8\ell+1<d/t\), which gives
\[
 u_{\ell,d}^2<\left(1+\frac{7t}{1+t}\right)^{-d}.
\]
For the floor \(t=1/40\), these inequalities bound the candidate
gaps and bases above the floor by finite rational tests. Taking the
maximum of the floor and each reciprocal allowed by the exact row
mean and declared cap gives the effective caps \(2/17,2/25,2/33,2/41\),
respectively. An order-eight cap \(c\) gives the order-nine cap
\(8c/(9-c)\).

## Common coefficient and ray displacement

The [exact checker](verify_multivariate_terminal_nine_shared_matched_upper.py)
uses the [matched skew and mean model](MULTIVARIATE_TERMINAL_NINE_SHARED_MATCHED_LIGHT.md#matched-skew-parameter)
with these effective caps. Fix a subdivision interval with midpoint
\(u_0\) and set \(h=(u-u_0)/u_0\). If \(|h|\le H<1\), then
\[
 -\frac{H^2}{2(1-H)}\le\log(u/u_0)-h\le0 .
\]
The upper bound is \(\log(1+h)\le h\), and the lower bound
follows by summing the absolute terms of the logarithm series
of degrees at least two. Separate the linear contribution
\(h e_{8,k}\) at order eight or \(h(64/63)e_{9,k}\) at order nine
from the rest of the mean correction. Include the resulting
\(h\)-direction and its products in the quadratic exponential
expansion. Thus the coefficient enclosure at index \(k\)
can be written with interval parameters as
\[
 A_k+hG_k+\sum_{l\ge2}p_l C_{l,k}.
\]
The base vector \(A\) includes the exponential remainder and the
omitted log-series tails. The vector \(G\) is the linear contribution
from the Gaussian coefficient shift just described.

Each ray cosine is affine in \(u\), so it has the form
\(x=x_0+b h\) with explicitly rational \(x_0,b\).
For the signed sine-ratio polynomial \(S_{k-i}(x)\) at projection
pivot \(i\), denote the actual value by \(D\). The enclosure
\[
 D=D_0+hD_1+R_D,\qquad
 D_0=S_{k-i}(x_0),\quad D_1=b S'_{k-i}(x_0)
\]
uses the Bernstein bounds for \(S''_{k-i}\) on the cosine interval.
The remainder belongs to their product with
\([0,b^2H^2/2]\). Multiply the coefficient model by this direction
and group the common linear displacement as
\[
 A_k(D_0+R_D)+h(A_kD_1+G_kD)+\sum_{l\ge2}p_l C_{l,k}D.
\]
This is an exact decomposition before taking interval enclosures.
The grouped \(h\)-coefficient preserves the cancellation between
the coefficient displacement and the ray displacement. Applying the
positive radial Bernstein transform and its dyadic subdivisions to
each group verifies strict projection signs.

The checker verifies the annulus \(1/12<|z|<12\).
The [same logarithmic homotopies](MULTIVARIATE_TERMINAL_NINE_SHARED_MATCHED_LIGHT.md#exact-projection-and-continuation)
connect actual arrays to the Gaussian reference arrays at the lower
endpoint of each band. Exact Rouché disks independently certify each
anchor and the alternating root orders. The single ninth-order root
in the sector joined across the negative axis is real. The actual
eighth moment is strictly positive on the real axis by its Gamma
representation, so its roots there are nonreal.

## Reproduction

Run from the repository root with

~~~sh
python3 research/factorial/image-rigidity/verify_multivariate_terminal_nine_shared_matched_upper.py --native
~~~

The optional --native flag compiles the
[exact Bernstein implementation](verify_multivariate_terminal_nine_shared_matched_upper.cpp)
in a temporary directory using a C++17 compiler with signed 128-bit
integers. Without the flag the checker uses the same dyadic interval
operations in Python. The native path divides each interval endpoint
product into quotient and remainder at denominator \(2^{48}\).
It checks input magnitudes that ensure all intermediate values fit
in a signed 128-bit integer. The wrapper uses the Python replay
for inputs outside those bounds. Acceptance requires strict signs
in outward interval arithmetic.

The replay verifies 54,443, 18,768, 5,211, and 5,559 radial cells
for the four respective bands, including their independently
certified Gaussian anchors.
