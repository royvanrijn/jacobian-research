# A profile-uniform root annulus for multivariate binomial moments

For arbitrary distinct exponent columns and positive Gamma shapes,
the outer geometric mean of the normalized coefficient array gives a
degree-independent bound on the moduli of all normalized roots.
For fully shared ordinary profiles of bounded curvature, this puts
both normalized scales of a possible common root into one compact
annulus, uniformly when the exponents and dimension vary with the
terminal. It also excludes profiles with a large total-degree difference.

## Product radii with both orientations

Let \(a_i,b_i\) be the original exponent columns. Delete equal rows,
put \(d_i=|b_i-a_i|\), \(\ell_i=\min(a_i,b_i)\), and define
\(\sigma_i=1\) for \(b_i>a_i\), \(-1\) for \(b_i<a_i\).
For Gamma shape \(\gamma_i>0\), set
\[
 l_{n,i}=(\gamma_i+n\ell_i)_{d_i},\qquad
 h_{n,i}=(\gamma_i+n\ell_i+(n-1)d_i)_{d_i}.
\]
For the actual scalar coefficient ratio \(t\), every reciprocal
root modulus lies in
\[
 \frac{\prod_{\sigma_i=1}l_{n,i}}
      {\prod_{\sigma_i=-1}h_{n,i}}
 \ \le\ |1/t|\ \le\
 \frac{\prod_{\sigma_i=1}h_{n,i}}
      {\prod_{\sigma_i=-1}l_{n,i}}.                          \tag{1}
\]
Both bounds are strict if \(n\ge2\). In particular two orders
are separated whenever the numerical intervals (1) at those
orders are disjoint or just touch at one endpoint.

For one increasing row these are the
[sharp Gamma-product annulus bounds](GAMMA_PRODUCT_SHARP_LOWER_ANNULUS_AND_SHAPE_SEPARATION.md#2-every-positive-shape-a-sharper-product-annulus).
Their single-factor lower endpoint is the classical Bessel circle
bound, also obtained from the cardioid theorem after substituting
\(y_n(2t;T-n+1)\) for a rising factor of base \(T>0\).
[de Bruin, Saff and Varga, definition (2.1) and Theorem 3.1](https://www.math.kent.edu/~varga/pub/paper_127.pdf)
Reversing the row polynomial reciprocates its root moduli.
The multivariate coefficients are the binomial-basis coefficient
products of the individual rows. The
[closed-disc composition lemma](GAP_TWO_GRACE_GLOBAL_LOCALIZATION.md#1-a-self-contained-closed-disc-convolution-lemma)
multiplies the upper reciprocal-root bounds. Reverse all factors
and the product to multiply the lower bounds. Strictness propagates
because at every finite order the individual strict annuli can be
replaced by slightly smaller closed annuli. This proves (1).

Write \(C_{n,k}\) for the coefficient divided by \(\binom nk\).
Let \(G_n=(C_{n,n}/C_{n,0})^{1/n}\), so the endpoint-normalized
variable is \(z=G_n t\), and put
\[
 u_n^2=\frac{C_{n,1}C_{n,n-1}}{C_{n,0}C_{n,n}} .
\]
In one increasing rising factor of base \(T\), the center
\(((T)_n)^{1/n}\) is between \(T\) and \(T+n-1\).
Reversal inverts the center and the radial bounds.
Thus \(G_n\) lies between the two product endpoints of (1).
Their ratio is \(u_n^{-2}\), so every normalized root satisfies
\[
 u_n^2\le |z|\le u_n^{-2}.                                   \tag{2}
\]
The endpoint ratio follows from
\(\prod_i h_{n,i}/l_{n,i}\), or by telescoping the outer
coefficient ratio in each row. This bound needs no total-curvature
hypothesis. If the lower-order \(u_n\) stays bounded away from
zero, the upper-order \(u_{n+1}\) does also, by the
[linked outer-mean comparison](MULTIVARIATE_LINKED_GAMMA_ANGLE_COMPACTIFICATION.md#linked-factors-and-comparison-across-orders).
Thus such sequences have both normalized scales in a fixed annulus
even when the number or orientation of rows changes.

## Shared rows and large degree imbalance

Suppose now that \(\ell_i\ge1\) for every unequal row and the shapes
are ordinary, \(\gamma_i=1\). Put
\[
 V=\sum_i \frac{d_i^2}{\ell_i},\qquad
 \Delta=\sum_i(b_i-a_i).
\]
For every normalized root at order \(n\ge2\),
\[
 |\log |z||<\frac{n-1}{2n}V .                                \tag{3}
\]
Indeed split a row into rising factors of bases
\(T_j=(1+n\ell_i+j)/d_i\), \(0\le j<d_i\).
The normalized center of such a factor has log
\(n^{-1}\sum_{k=0}^{n-1}\log(T_j+k)\).
Its distance to either logarithmic endpoint \(\log T_j\) or
\(\log(T_j+n-1)\) is bounded by \((n-1)/(2T_j)\), by
\(\log(1+x)\le x\) and summation over \(k\).
Orientation reversal changes the sign of these logarithmic
coordinates and preserves the bound. Adding over factors and using
\(T_j\ge n\ell_i/d_i\) gives (3).

For consecutive orders \(n,N=n+1\), put \(\mu_n=\log G_n\).
The same integral comparison as in the
[varying-profile phase estimate](MULTIVARIATE_VARYING_PROFILE_PHASE.md#uniform-estimates-with-an-order-derivative)
gives the explicit bound
\[
 |\mu_N-\mu_n-\Delta\log(N/n)|
   \le\frac{1}{nN}\sum_i|\log(b_i/a_i)|\le\frac{V}{nN}.
\]
Here \(|\log(b_i/a_i)|\le d_i/\ell_i\le d_i^2/\ell_i\).
If there is a common root, its two normalized moduli differ by
\(e^{\mu_N-\mu_n}\), so (3) gives
\[
 |\Delta|\log(N/n)
 <\left(1-\frac{2n-1}{2nN}\right)V,\qquad
 |\Delta|<\left(n+\frac1{2n}\right)V.                        \tag{4}
\]
The second inequality uses \(\log(N/n)>1/N\).
Failure of (4) is a sufficient condition for coprimality.

To express this in the symmetric curvature, put
\(m_i=(a_i+b_i)/2\),
\(W=\sum_i d_i^2/m_i\), and
\(R=\max_i d_i/(a_i+b_i)\). Then \(V\le W/(1-R)\).
On a bounded-\(W\) band of integer shared rows, \(R\) stays
uniformly below one by the
[varying-profile estimate](MULTIVARIATE_VARYING_PROFILE_PHASE.md).
Consequently (3) holds in a fixed annulus and a common root
requires \(|\Delta|=O(n)\), even along arbitrary escaping profiles.
A simpler small-curvature consequence is that \(W\le1/(n+1)\)
with \(\Delta\ne0\) suffices for separation. Indeed
\(2R\le W\), so \(V\le W/(1-W/2)\), and the right side of (4)
is then at most one. The proof excludes a nonzero integer \(\Delta\)
in that range.

Combining (3) with the
[uniform subcritical saddle](MULTIVARIATE_VARYING_PROFILE_SUBCRITICAL_SADDLE.md)
also gives a fixed zero-free angle about the positive real axis
at all sufficiently large orders for varying shared profiles with
bounded \(W\) and a fixed subcritical margin. No extra assumption
on the normalized scale is needed, since the annulus supplies it.
The interior complex angles for bounded degree imbalance remain
outside these conclusions. For ordinary moments, all the stated
coprimality results transfer by rational Bézout to every
characteristic-zero field.
