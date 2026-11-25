# Terminal-eleven fixed-ray separation for shared loss and private rows

Let \(a,b\) be distinct ordinary exponent columns. Write
\(S_0=\sum_i d_i^2/(\ell_i+d_i)\) over their shared unequal rows,
where \(\ell_i\) is the smaller exponent and \(d_i\) the gap.
Let \(p,q\) be the private degrees on the two sides and put
\(S=S_0+p+q\). If \(S\ge8\), or if a case in the table below
applies, the order-ten and order-eleven factorial moments have
disjoint root arguments under independent positive rescalings.
In particular they are coprime over every characteristic-zero field.
For fully shared profiles \(S=S_0\).

| Increasing private gaps | Decreasing private gaps | Required shared loss \(S_0\) |
| --- | --- | ---: |
| none | none | \(8\) |
| \(1\) | none | \(6\) |
| \(1,1,1,1\) | none | \(0\) |
| \(1,1\) | \(1,1\) | \(0\) |
| \(1,1,1\) | \(1\) | \(0\) |
| \(3\) | none | \(1/4\) |
| \(1,1,1\) | none | \(9/8\) |
| \(1,2\) | none | \(5/8\) |
| \(1,1\) | none | \(27/8\) |
| \(2\) | none | \(11/4\) |
| \(1\) | \(1\) | \(13/4\) |
| \(1,1\) | \(1\) | \(7/8\) |
| \(2\) | \(1\) | \(3/8\) |

Each box also covers profiles obtained by adding private factors
or increasing their gaps, and column reversal preserves the
conclusion. In particular the zero-loss cases cover a private
degree at least four on either side, private degrees at least
two on both sides, or a private degree at least three with a
positive opposite private degree. These three envelopes exhaust
total private degree at least four.

## Coefficient boxes

Normalize the endpoints to one by positive rescaling. An
increasing private gap \(d\) contributes
\[
 f_{j,k}(d)=\frac{(dk)!}{((dj)!)^{k/j}},
\]
and a decreasing gap contributes \(f_{j,j-k}(d)\). The
[private-factor and shared loss bounds](MULTIVARIATE_TERMINAL_SEVEN_PRIVATE_ROWS.md#private-factors-and-shared-loss)
give \(f_{j,k}(d)\le1\), monotonicity in \(d\), and
\(f_{j,k}(d)\le f_{j,k}(1)^d\). For the last bound,
the successive ratios of \((dk)!/(k!)^d\) are
\(\prod_{h=1}^d(dk+h)/(k+1)\), whose individual factors
are nondecreasing in \(k\). Applying its logarithmic chord gives
the stated inequality. Hence \(p\) units of private degree
are bounded by \(p\) unit factors in that orientation.

The second logarithmic differences for a shared row are bounded
below by \(d^2/(j(\ell+d))\), by \(\log(1+t)\ge t/(1+t)\).
The chord-normalized shared factors are therefore bounded by
\(\exp[-S_0 k(j-k)/(2j)]\). Together with the private factors
and \(\binom jk\), these give upper bounds for the interior
coefficients. The same second-difference calculation applies
when \(\ell=0\). At an interior index \(k\), each summand is
\(\log(1+d/((k-1)d+h))\ge d/(kd+h)\ge1/j\).
Summing over \(1\le h\le d\) gives the lower bound \(d/j\)
for the private row. Hence the coefficient caps
\(\binom jk\exp[-S k(j-k)/(2j)]\) apply to all profiles and
the first fixed box certifies \(S\ge8\).
The checker encloses the exponential reciprocals with finite
positive Taylor sums and rounds the resulting roots outward by
integer-power comparisons. It also checks constant
and leading coefficient domination at radii \(1/8\) and \(8\).

## Exact ray and anchor replay

The [checker](verify_multivariate_terminal_eleven_high_curvature.py)
uses nine fixed rational ray cosines for each case. Its box allows
each interior coefficient to range from zero to its cap independently,
with the endpoints fixed to one. The reference arrays \(1+z^{10}\),
\(1+z^{11}\) belong to these boxes. The
[exact Bernstein projection](MULTIVARIATE_TERMINAL_SIX_JOINT_RAY_BANDS.md#moving-rays-and-exact-replay)
excludes zeros on every ray of the annulus throughout the straight
coefficient segments joining the references to the actual arrays.

At the binomial references the angle-\(\theta\) sign for even order
\(2m\) is the sign of \(\cos(m\theta)\). For odd order \(2m+1\)
it is the sign of \(\cos((2m+1)\theta/2)/\cos(\theta/2)\) on
\(0\le\theta<\pi\), with its limiting sign at \(\pi\).
These are polynomials in the ray cosine, so their signs are
checked exactly at all sector boundaries. In increasing upper
argument the ten sectors contain one root each, of orders
\(11,10,11,10,11,10,11,10,11,10\). The odd reference also has
one negative real root. Join the first and last sectors to
their reflected lower sectors. The sole order-eleven root in
the last joined sector stays negative real by conjugacy.
The actual order-ten factorial moment is strictly positive on
the real axis by its Gamma integral, so the order-ten pair in
that sector is nonreal. The arguments are disjoint there too.
Rational Bézout gives the characteristic-zero claim.

Replay from the repository root with

~~~sh
python3 research/factorial/image-rigidity/verify_multivariate_terminal_eleven_high_curvature.py
~~~

The replay checks 13 coefficient boxes at both orders, 286 exact
anchor signs and 1,326 radial cells. The ray table is frozen.
