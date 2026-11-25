# Terminal-twelve Cayley separation on the band .66–.70

Let two distinct ordinary exponent columns have only shared unequal
rows. Let \(u\) be their endpoint-normalized, binomial-stripped
order-eleven outer mean. If \(33/50\le u\le7/10\) and every
order-eleven Gauss reciprocal \(d/(11\ell+h)\) is at most \(1/25\),
the order-eleven and order-twelve factorial moment polynomials
have disjoint root arguments under independent positive rescalings.
Their coprimality over characteristic-zero fields follows by
rational Bézout.

## Exact reduction of the reciprocal cap

For an unequal shared row let \(x=d/(11\ell+1)\) be its largest
reciprocal. If \(d=1\), the initial cap forces \(\ell\ge3\)
and hence \(x\le1/34\). If \(d=2\), it forces \(\ell\ge5\)
and \(x\le1/28\). If \(d\ge3\), each row reciprocal is at
least \(x/(1+x)\), and
\[
 u^{-2}
 \ge \left(1+\frac{10x}{1+x}\right)^d.
\]
For \(x\ge1/28\) the right side is at least
\((39/29)^3>(50/33)^2\), which contradicts \(u\ge33/50\).
Thus \(1/28\) is a common cap on the full ordinary band.

## Comparison and exact checks

Apply the [leading-skew Cayley model](MULTIVARIATE_TERMINAL_TWELVE_CAYLEY_LOWER.md)
at the reduced cap. The 10 rays have affine levels in the
order-eleven deficit \(p=1-u\) on each five-thousandth subband.
The checker proves leading dominance in the scaled disk of
radius three and verifies that the physical rays are in
the left half-plane. The closest ray uses the eight signed
second-moment subintervals, and the other rays use the unconditioned
grouped model. The same logarithmic-moment continuation at fixed
deficits connects each actual polynomial to its pure Gaussian array.
Following the Gaussian comparison ratio then brings the arrays to
the checked alternating anchor at \(p=.300\). The sector counts and
the real-axis conclusion follow as in the neighboring band.

The [exact checker](verify_multivariate_terminal_twelve_cayley_bridge.py)
verifies 61,105 projection cells across 160 order, ray, and
subband checks.

~~~sh
python3 research/factorial/image-rigidity/verify_multivariate_terminal_twelve_cayley_bridge.py
python3 research/factorial/image-rigidity/verify_multivariate_terminal_twelve_cayley_bridge.py 12 0 300 305
~~~
