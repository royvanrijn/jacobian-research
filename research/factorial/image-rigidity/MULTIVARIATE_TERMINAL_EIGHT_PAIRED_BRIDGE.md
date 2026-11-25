# Paired terminal-eight separation on the intermediate outer-mean bands

Let all unequal rows of distinct ordinary exponent columns be
shared. Put \(\ell_i=\min(a_i,b_i)\) and \(d_i=|b_i-a_i|\) on
these rows, and let \(u\) be the normalized order-seven outer
geometric mean. The order-seven and order-eight factorial moments
have disjoint root arguments under independent positive
rescalings on the following bands.

| Range of \(u\) | Maximum \(d_i/(7\ell_i+1)\) |
| --- | ---: |
| \([.48,.52]\) | \(1/8\) |
| \([.52,.60]\) | \(2/25\) |
| \([.60,.66]\) | \(13/200\) |
| \([.66,.70]\) | \(9/200\) |

Rational Bézout gives coprimality over every characteristic-zero
field. The model and orientation reduction are those of the
[paired light proof](MULTIVARIATE_TERMINAL_EIGHT_PAIRED_LIGHT.md).
In particular each reciprocal cap is at most \(1/8\), so both
the exact single-factor skew transfer and the proof of convexity
of the local skew magnitudes apply.

## Paired projections

The [bridge checker](verify_multivariate_terminal_eight_paired_bridge.py)
uses the same correlated Gaussian coefficient model and annulus
\(1/12<|t|<12\). It keeps the new ray cosines affine in \(u\), and
subdivides both \(u\) and the normalized order-seven inner skew.
For the sine-ratio term in a projection use the polynomial
\(S_m(x)=\sin(m\theta)/\sin\theta\), \(x=\cos\theta\).
The integer coefficient recurrence is \(S_0=0\), \(S_1=1\),
\(S_{m+1}=2xS_m-S_{m-1}\), and \(S_{-m}=-S_m\).
The checker takes the exact Bernstein coefficient bounds of
these polynomials on each cosine interval. This retains the
sine-ratio cancellations near \(x=-1\).

As in the paired light proof, the Bernstein coefficients
in the radial variable preserve each shared perturbation
as one contribution. Strict outward bounds on the grouped
projections certify all accepted radial cells. If the bounds at
the midpoint of a radial cell do not certify any projection,
the checker subdivides the parameter box. That midpoint test
uses the nonnegative weights \(\binom jm/2^j\) in the Bernstein
basis of order \(j\); it never certifies a radial cell on its own.

At fixed \(u,\rho\), first scale all skew and defect corrections
in the logarithms of the normalized multipliers to zero. This
preserves the linear coefficient constraints and the cross-order
skew relation. Then vary \(\rho\) to \(49/48\) and move \(u\) to
.48 along the Gaussian model. The ray cosines at .48 match those
of the preceding paired-light proof exactly. Its Gaussian
anchor counts consequently apply. As there, the sole septic root
in the last sector joined to its conjugate is real, and the octic
pair there is nonreal by positivity of the actual eighth factorial
moment on the real axis. Hence the root arguments are disjoint.

Replay from the repository root with

~~~sh
python3 -O research/factorial/image-rigidity/verify_multivariate_terminal_eight_paired_bridge.py
~~~
