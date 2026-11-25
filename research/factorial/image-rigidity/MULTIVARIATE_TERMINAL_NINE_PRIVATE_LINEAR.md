# Terminal-nine separation for one private linear row

Let \(a,b\) be distinct ordinary exponent columns with total private
degree one. The eighth and ninth factorial moment polynomials have
disjoint root arguments under independent positive rescalings and
are coprime over every characteristic-zero field. Together with the
[larger private-degree theorem](MULTIVARIATE_TERMINAL_NINE_PRIVATE_ROWS.md),
this covers every ordinary profile with a private row at terminal nine.

Orient the unique private row as \((0,1)\). Its normalized coefficients,
including binomial factors, are
\(B_{j,k}=\binom jk\, k!/(j!)^{k/j}\).
For shared curvature \(S_0\ge21/4\), the coefficient caps are
\(B_{j,k}\exp[-S_0 k(j-k)/(2j)]\). At \(S_0=21/4\) the
[exact checker](verify_multivariate_terminal_nine_private_linear.py)
compares these caps on seven fixed upper-half-plane rays. The same
boxes enclose all larger strengths, with the binomial array as anchor.

Below \(21/4\), use the shared means \(u,w\) at orders eight and nine.
The reciprocal enumeration, linked range of \(\rho=(-\log w)/(-\log u)\),
and independent exponent boxes are proved in the
[terminal-nine outer-mean comparison](MULTIVARIATE_TERMINAL_NINE_PRIVATE_ROWS.md#shared-outer-mean-boxes).
Those arguments use only the given upper bound on \(S_0\), and the
checker replays the row enumeration with \(S_0\le21/4\).
It subdivides the \(u\)-intervals and can also use the following grouped
model on the intermediate range \(24/100\le u\le88/100\).

## Grouping the shared logarithmic corrections

For one Gauss factor at order \(j\), write \(m=j-1\),
\(s=\frac12\log(1+mx)\), and \(T=\tanh s=mx/(2+mx)\), where \(x\)
denotes that order's reciprocal. Put \(p_h=1-2h/m\).
For a coefficient at distance \(r=\min(k,j-k)\) from the
endpoints, its reflected mean loss is
\(\sum_{h=0}^{r-1}\operatorname{artanh}(p_hT)\).
Since \(\sum_{h<r}p_h=e_{j,k}=k(j-k)/(j-1)\), the defect
from the reference mean loss \(e_{j,k}s\) is
\[
 \sum_{l\ge1}\frac{T^{2l+1}}{2l+1}
           \sum_{h<r}(p_h-p_h^{2l+1}).
\]
Every inner sum is nonnegative, since \(0<p_h\le1\) in this range.

For skew use
\(\omega_h=(\mathbf 1_{\{r\le h<j-r\}}-(j-2r)/j)/2\).
The weights sum to zero and satisfy \(\omega_{m-h}=\omega_h\).
As \(1+hx=(1+mx/2)(1-p_h T)\), its log gives the local
skew up to orientation as
\[
 \sum_h\omega_h\log(1+hx)
      =-\sum_{l\ge1}\frac{T^{2l}}{2l}\sum_h\omega_h p_h^{2l}.
\]
The contributions for indices \(k\) and \(j-k\) have opposite
signs. Following the reference convention for one orientation
and allowing the coefficient of each skew vector to have either
sign also covers every mixture of row orientations.

The checker retains the first two mean-defect vectors and the
first three skew vectors displayed by these expansions.
Let \(T_*\) be the upper \(T\) from the reciprocal cap and let
\(U_j\) be the sum of \(s\)'s in the order. Since \(s\ge T\), the
scalar on the \(l\)-th mean vector is in
\([0,U_jT_*^{2l}/(2l+1)]\). The absolute scalar on the \(l\)-th
skew vector is at most \(U_jT_*^{2l-1}/(2l)\).
The two omitted tails at coefficient \(k\) are bounded in
absolute sum by
\[
 E_k=U_j\left[
  \frac{rT_*^6}{7(1-T_*^2)}
  +\frac{T_*^7\sum_h|\omega_h|}{8(1-T_*^2)}
 \right].
\]
This follows by bounding the remaining series geometrically
and using \(|p_h|\le1\). All formulas remain valid when there
are no shared unequal rows.

In a subinterval centered at \(u_0>0\), the reference coefficients
are \(B_{8,k}u_0^{e_{8,k}}\) and
\(B_{9,k}u_0^{(64/63)e_{9,k}}\).
Relative to this center the outer mean contributes the Gaussian
vector \(e_{j,k}\) with scalar \(\log(u/u_0)\) at order eight
and scalar
\[
 \tfrac{64}{63}\log(u/u_0)+(\tfrac{64}{63}-\rho)(-\log u)
\]
at order nine. The other scalar bounds use
\(U_8=-\log u\) and \(U_9=-\rho\log u\).
If \(\tau_k\) is the resulting retained logarithmic correction
and \(\sigma_k\) the exact one, then
\(|\sigma_k-\tau_k|\le E_k\). The checker obtains \(M_k<1\)
bounding \(|\tau_k|\) and \(|\sigma_k|\), and uses
\[
 \left|e^{\sigma_k}-1-\tau_k-\frac{\tau_k^2}{2}\right|
 \le \frac{M_k^3}{6(1-M_k)}+E_k(1+M_k).
\]
The cubic bound uses \(n!\ge6\) for \(n\ge3\).
The square difference uses
\(|\sigma_k^2-\tau_k^2|\le2M_k|\sigma_k-\tau_k|\).
Multiplying the reference by this model keeps each linear
vector and each product of two vectors grouped in the ray
projection. This yields a convex coefficient enclosure around
the reference sequence.

## Rays and continuation

The seven ray cosines are piecewise affine in \(u\), with rational
knots. The checker uses exact rational bounds or outward dyadic
bounds at every acceptance test. It verifies constant and leading
coefficient dominance on the two annulus boundaries and certifies
strict projected ray signs on all radial cells by the
[paired Bernstein test](MULTIVARIATE_TERMINAL_SEVEN_PAIRED_LIGHT.md#projection-replay)
or the rational independent test.

At the actual \(u\), the chosen enclosures contain both the true
arrays and the reference arrays with factors \(u^{e_{8,k}}\) and
\(u^{(64/63)e_{9,k}}\). In the grouped model the latter have the
same mean displacement with mean defects, skew, and transport
error set to zero, using \(\rho=64/63\). Because the enclosures
are convex, the straight coefficient segments exclude all seven
rays. Then move \(u\) to zero along the reference arrays, using
continuous ordered rays. The resulting binomial anchors have upper
sector orders \(9,8,9,8,9,8,9,8\), as checked by exact anchor
signs. The same counts hold in the high-strength branch.

Join the first and last sectors to their conjugates. The single
ninth-order root in the last joined sector stays real and negative.
The actual eighth-order polynomial has no real root by the strict
positivity of its Gamma integral. The actual root arguments are
therefore disjoint, and rational Bézout transfers coprimality to
every characteristic-zero field.

Replay from the repository root with

~~~sh
python3 research/factorial/image-rigidity/verify_multivariate_terminal_nine_private_linear.py
~~~

The replay verifies 7 skew convexity inequalities, 92 high-strength
radial cells, and 13,399 low-strength radial cells. The reciprocal-cap
enumeration retains 700 candidate row types.
