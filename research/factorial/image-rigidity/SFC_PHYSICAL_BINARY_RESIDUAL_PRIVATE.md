# A binary physical merger with private degree left on both sides

Let two monomials be built from arbitrary integer exponent
rows. In a Gaussian row-merger path, correlate a selected
pair of coordinates and leave all other coordinates independent.
Suppose the untouched rows supply total private exponent
at least two on each side. Then the two consecutive moment
forms have no common nonzero projective zero at any positive
starting order and any squared correlation \(q\in[0,1]\).

The same statement holds for several correlated coordinates
confined to the base, provided all the reserved private
coordinates stay independent of the base and of each other.
The ordinary endpoints therefore permit the row merger
and any further one-sided private unit-row inflation
that leaves the reserved private variables in place.
There is no balance or exchangeability hypothesis.

## A distribution-independent coefficient bound

In fact, let \(U_0,V_0>0\) be any jointly distributed base
features with finite positive mixed moments at every relevant
order. Take independent unit-rate exponentials \(X_i,Y_j\),
also independent of the base, and let
\[
 U=U_0\prod_i X_i^{a_i},\qquad
 V=V_0\prod_j Y_j^{b_j},\qquad
 a_i,b_j\in\mathbb Z_{>0},\quad \sum_i a_i,\sum_j b_j\ge 2 .
\]
At order \(m\) write
\[
 C_{m,k}=\mathbb E(U^{m-k}V^k),\qquad
 P_m(t)=\sum_{k=0}^m\binom mk C_{m,k}t^k.
\]
After dividing by \(C_{m,0}\) and putting
\(z=t(C_{m,m}/C_{m,0})^{1/m}\), the normalized polynomial
has the form
\[
                 T_m(z)=1+z^m+\sum_{k=1}^{m-1} c_{m,k}z^k .
\]
The base factor has normalized coefficients at most one,
since Hölder gives
\[
 \mathbb E(U_0^{m-k}V_0^k)\le
   (\mathbb E U_0^m)^{1-k/m}(\mathbb E V_0^m)^{k/m}.
\]
For a private exponent \(a\), the sequence
\((aj)!/(j!)^a\) starts at one and is log-convex in \(j\).
Its consecutive ratio is
\(\prod_{\ell=1}^a (a-(a-\ell)/(j+1))\), which increases in
\(j\). Endpoint normalization therefore bounds the private
factor belonging to \(k\) by \((k!/(m!)^{k/m})^a\),
and the other side uses \(m-k\) in place of \(k\).
Since both private degrees are at least two, the product
of the normalized private factors is at most
\((k!(m-k)!/m!)^2\). Hence
\[
                   0<c_{m,k}\le \binom mk^{-1}.           \tag{1}
\]
This is the coefficient suppression of the
[ordinary private-degree proof](MULTIVARIATE_TWO_MONOMIAL_NONDIVISIBILITY.md#all-orders-when-both-private-exponent-sums-are-at-least-two).
Hölder makes it valid for the entire positive base.

## Angular separation at every order

For \(m\ge3\) set
\(\epsilon_m=\sum_{k=1}^{m-1}\binom m k^{-1}\) and
\(e_m=\arcsin\epsilon_m\). Counting the two outer terms
exactly and bounding the others by
\(2/(m(m-1))\) gives
\[
            \epsilon_m\le\frac{4(m-2)}{m(m-1)}\le2/3 .
\]
For a root of \(T_m\) of modulus at most one, (1) gives
\(|1+z^m|\le\epsilon_m\). At larger modulus apply the
same argument to the reciprocal polynomial and use the
symmetry of the odd angle grid. By conjugation a root
angle is taken in \(0<\theta\le\pi\), with zero excluded
by coefficient positivity. It satisfies
\[
                 |m\theta-(2j_m+1)\pi|\le e_m
                            \quad (j_m\ge0).             \tag{2}
\]
Positive normalizations preserve the argument of a proposed
common physical root of the two consecutive orders.

If \(m\ge3\) and the same angle occurs at both orders,
subtract (2). Since \(e_m+e_{m+1}<\pi\) and
\(0<\theta\le\pi\), the two grid labels agree and
\(\theta\le e_m+e_{m+1}\). On \([0,2/3]\),
\(\arcsin x\le (11/10)x\), by convexity and
\(\sin(11/15)>2/3\). Thus
\[
 (2j_m+1)\pi\le(m+1)e_m+m e_{m+1}<44/5<3\pi .
\]
Here each term is less than \(22/5\) by the bound on
\(\epsilon_m\), so \(j_m=j_{m+1}=0\).
For \(n\ge5\), count the first two and last two terms
exactly and bound the rest by \(\binom n3^{-1}\). This gives
\[
 (n-1)\epsilon_n
   \le 2+\frac{8n-34}{n(n-2)}\le\frac83 .
\]
The same bound holds at \(n=4\) directly. Consequently
\(e_n<\pi/(n-1)\) for all \(n\ge4\).
With \(n=m+1\), (2) would force \(\theta<\pi/m\).
Every term of \(T_m\) has nonnegative imaginary part
when \(0<\theta\le\pi/m\), and the linear term is
strictly positive in imaginary part. This is impossible.

For \(m=1\) the sole root is negative real, whereas
\(T_2(z)=1+c_{2,1}z+z^2\), \(0<c_{2,1}\le1/2\), has
two nonreal roots. For \(m=2\), its upper root has angle
\(\pi/2<\theta\le\pi/2+\arcsin(1/4)<2\pi/3\).
Thus \(3\theta\) has distance greater than \(\pi/2\) from
every odd multiple of \(\pi\), contrary to (2) for
order three. This covers every positive starting order.

For the physical Gaussian application, all rows involved
in the correlation are included in \((U_0,V_0)\).
Its features are positive almost surely, and their
mixed moments are finite, also at \(q=1\). No continuity
or monotonicity in \(q\) is used. If one projective
coefficient vanishes, the remaining pure moments
are positive. Over the ordinary endpoints the
unnormalized moments have integer coefficients, so a
rational Bézout identity transfers the complex
coprimality to every characteristic-zero field.
