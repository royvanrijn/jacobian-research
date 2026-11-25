# Terminal-six separation with one small private side

Suppose the two ordinary exponent columns have private rows on just
one side. Remove equal rows and orient the private rows to increase.
For each shared unequal row let \(\ell_i=\min(a_i,b_i)>0\) and
\(d_i=|a_i-b_i|\), and set
\[
 V=\sum_{\rm shared}\frac{d_i^2}{\ell_i},\qquad
 S_0=\sum_{\rm shared}\frac{d_i^2}{\ell_i+d_i}.
\]
The ordinary order-five and order-six polynomials are coprime in
these ranges, uniformly in dimension and the exponents.

| Positive gaps of the private rows | Allowed \(V\) | Sufficient bound on \(S_0\) |
| --- | ---: | ---: |
| \(1\) | \(9/100\) | \(9/109\) |
| \(1,1\) | \(7/8\) | \(1/2\) |
| \(2\) | \(13/10\) | \(13/20\) |

The angular proof covers the \(V\) bounds. The endpoint
\(S_0=1/2\) in the second row also uses an Eisenstein test.

## Shared coefficient lower bounds

After removing the binomial coefficient and normalizing both
endpoints, let \(F_{j,k}^{\rm sh}\) be the contribution of the
shared rows to the order-\(j\) coefficient. In a shared row the
second difference at index \(k\) is nonnegative by log-convexity,
and for an increasing orientation has upper bound
\[
 \sum_{h=1}^{d_i}
 \log\left(1+\frac{d_i}{j\ell_i+(k-1)d_i+h}\right)
 \le\frac{d_i^2}{j\ell_i}.
\]
Reflection preserves the bound. After summing and comparing to
the quadratic with the same endpoints and constant second
difference \(V/j\), we obtain
\[
 \exp\left(-\frac{Vk(j-k)}{2j}\right)
 \le F_{j,k}^{\rm sh}\le1 .
 \tag{1}
\]
The private contribution for each private gap \(e\) is exactly
\((ek)!/((ej)!)^{k/j}\). Thus (1) encloses every coefficient by
explicit lower and upper bounds around the polynomial of the
private rows alone.

## Exact angle certificates

The
[checker](verify_multivariate_terminal_six_small_shared.py)
uses the ray projections (5) in the
[private degree-three proof](MULTIVARIATE_TERMINAL_SIX_PRIVATE_ANGLES.md#exact-ray-inequalities).
The four decreasing ray cosines for each case are respectively
\[
 \begin{array}{c|rrrr}
 (1)&.146&-.196&-.640&-.771\\
 (1,1)&.616&.160&-.376&-.720\\
 (2)&.688&.200&-.322&-.710
 \end{array}
\]
where all decimals stand for exact rationals.
Raise the coefficient bounds to the order \(j\) to leave only
a rational bound on the exponential. For \(0\le u<66\),
\[
 \sum_{h=0}^{64}\frac{u^h}{h!}\ \le\ e^u\ \le\
 \sum_{h=0}^{64}\frac{u^h}{h!}+\frac{u^{65}}{65!(1-u/66)} .
\]
The upper estimate bounds subsequent Taylor terms by a geometric
series. Invert that upper estimate to get the lower coefficient
bound. Enclose the resulting algebraic bounds with rational
numbers of denominator \(10000\) by integer powers.

On the circles of radii \(1/5\) and \(5\), the constant and
leading terms respectively dominate for all admissible
coefficients. Thus only the closed radius interval
\([1/5,5]\) needs testing along each ray. The checker divides
that interval into rational subintervals and verifies that
one projection has strictly positive or strictly negative
sign uniformly on each. The coefficient bounds remain
valid while linearly interpolating each coefficient from its
private-only value to its actual value.

To verify the private-only sector counts, the checker uses
the unnormalized polynomial
\(\sum_k\binom jk\prod_e(ek)!\,t^k\), whose variable differs
only by a positive scaling. For each explicitly given
rational complex center \(z_0\) it expands this polynomial
as \(\sum_m b_m(t-z_0)^m\). On the circle of radius \(1/10000\)
centered at \(z_0\) the lower bound on \(|b_1(t-z_0)|\)
strictly dominates the upper bound on the other terms.
Both estimates use rational real and imaginary parts.
These disjoint disks account for all roots. The checker
places each disk away from the separating rays using the
rational cosines and rational bounds on their positive sines.

In each case the three sextic upper roots lie respectively
before the first ray, between the second and third rays,
and after the fourth ray. The two quintic upper roots lie
respectively between the first and second rays and between
the third and fourth rays. There is also a negative real
quintic root. Along the interpolation, the two sectors about
the real axes are counted as full conjugate-symmetric sectors.
The final sector then still contains precisely one quintic
root, which must be real by conjugation. The actual ordinary
sextic is strictly positive on the real axis. The argument
separation proves coprimality and descends by rational Bézout.

## Conversion to \(S_0\)

For \(S_0<1\), each shared row has
\(d_i/(\ell_i+d_i)\le S_0\), so \(V\le S_0/(1-S_0)\).
This proves the first sufficient \(S_0\) bound.

When \(S_0\le1/2\), a shared row with \(\ell_i=1,d_i=1\)
uses the full curvature budget, and is then the only unequal
shared row. In the case of two private unit rows, apply the
[marked-row test](MULTIVARIATE_MARKED_ROW_VARIABLE_PRIME.md)
to \((1,2)\) or its reversal with \(N=6\) and \(p=7\).
All other exponents are at most one. If there is no such shared
row, then \(d_i/(\ell_i+d_i)\le1/3\) for every row. This holds
for a unit gap since \(\ell_i\ge2\), and for \(d_i\ge2\) since
\(d_i/(\ell_i+d_i)\le S_0/d_i\le1/4\). Hence
\(V\le(3/2)S_0\le3/4<7/8\).

Finally, if \(S_0\le13/20\) then \(\ell_i\ge d_i\) in every
shared row. Otherwise \(d_i\ge2\) and
\(d_i^2/(\ell_i+d_i)>d_i/2\ge1\).
Therefore \(V\le2S_0\le13/10\), as required.
