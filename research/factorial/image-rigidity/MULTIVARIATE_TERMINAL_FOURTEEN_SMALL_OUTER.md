# Terminal-fourteen exact separation near zero outer mean

For fully shared ordinary profiles, use the endpoint-normalized
order-thirteen reflected outer mean \(u\). The order-thirteen and
order-fourteen factorial polynomials have disjoint arguments for
\(0<u\le7/400\) if all lower-order reciprocals
\(d/(13\ell+1)\) are at most \(1/8\). Together with the
[retained profile result](MULTIVARIATE_TERMINAL_FOURTEEN_RETAINED_PROFILES.md),
this covers every fully shared profile on that interval.

At consecutive orders, the
[independent exponent estimates](MULTIVARIATE_TERMINAL_THIRTEEN_SHARED_LIGHT.md#matched-logarithmic-enclosure)
at reciprocal cap \(1/8\) apply to the updated orders 13 and 14.
For each coefficient they give a positive lower and upper exponent
with respect to the lower order reflected outer mean. The Gaussian
nominal exponent lies between the two. The normalized endpoint entries
are one. On \(u=0\), the coefficient array has only the two
endpoint entries. The checker tests the annulus and the expected
initial sector counts directly.

The twelve rays have rational cosines affine in \(u\) on each
interval between consecutive knots \(j/400\). The checker
bounds every polynomial coefficient and each ray cosine
on a cell. Exact radial Bernstein projection then excludes roots
on the rays in the closed annulus. The homotopy to the binomial
anchors can be taken continuously by first interpolating
coefficients to their Gaussian values at fixed \(u\) and then
varying \(u\) to zero. The annulus and ray checks apply to
every point of the homotopy.

Join the sectors adjacent to the two real rays with their conjugate
sectors. Their root counts follow from the exact initial counts
of \(1+z^{13}\) and \(1+z^{14}\). The odd negative root remains
real by conjugation and uniqueness in its sector. The even
moment is strictly positive on the real axis by the Gamma
integral. The two actual polynomials consequently have disjoint
arguments throughout both real-adjacent sectors as well as the
interior sectors.

Replay the seven closed intervals with the
[exact checker](verify_multivariate_terminal_fourteen_small_outer.py).
There are 13,731 radial cells, with 3,080 at order thirteen and
10,651 at order fourteen.
