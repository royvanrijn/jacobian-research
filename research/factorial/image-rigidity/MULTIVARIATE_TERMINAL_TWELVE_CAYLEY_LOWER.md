# Terminal-twelve Cayley separation on the band .70–.76

Let two distinct ordinary exponent columns have only shared unequal
rows. Let \(u\) be their endpoint-normalized, binomial-stripped
order-eleven outer mean. If \(7/10\le u\le19/25\) and every
order-eleven Gauss reciprocal \(d/(11\ell+h)\) is at most \(3/100\),
the order-eleven and order-twelve factorial moment polynomials
have disjoint root arguments under independent positive rescalings.
Their coprimality follows by rational Bézout over every
characteristic-zero field.

Use the [terminal-twelve Cayley model](MULTIVARIATE_TERMINAL_TWELVE_CAYLEY.md)
with weight nine, the same comparison scales \(p=1-u\), \(q=1-w\)
and \(r=q/p\), and the scaled root disk of radius three.
All ten rays lie in the physical left half-plane. As in the
upper bands their positive Cayley height is relevant only for
\(0\le\Re Z\le1/\lambda\), where \(\lambda\) is the square root
of the appropriate order deficit. The checker verifies the leading
coefficient dominance and the left half-plane ray inequalities
in every subband. For the closest ray it refines the reciprocal
cap by the integer row argument in the upper-band proof.

## Condition the logarithmic moments on the leading skew

At either order write \(h\) for its deficit and \(a_i=T_i/h\),
with \(T_i\) its normalized row factors. In the Cayley model
\(0\le a_i\le A\) and \(\sum_i a_i\le C\). Write
\[
 P=\frac12\sum_i\epsilon_i a_i^2
\]
for the signed logarithmic moment of weight two, where
\(\epsilon_i\in\{-1,1\}\) is its row orientation. The leading
ray checker partitions \([-CA/2,CA/2]\) into eight equal
subintervals. In a subinterval \([P_-,P_+]\) let \(P_0\) be
the minimum of \(|P|\). For odd weights \(d\ge3\) the
unsigned logarithmic moment lies in
\[
 \left[\frac{(2P_0)^{d-1}}{d C^{d-2}},\;
       \frac{C A^{d-1}}{d}\right].
\]
In fact \(\sum_i a_i^2\ge2P_0\), and Jensen's inequality for
the weights \(a_i/\sum_i a_i\) gives the lower bound on
\(\sum_i a_i^d\). The upper bound follows termwise.
For even weights \(d\ge4\) the signed logarithmic moment lies in
\[
 \left[-\frac{(CA-2P_-) A^{d-2}}{2d},\;
         \frac{(CA+2P_+) A^{d-2}}{2d}\right].
\]
Indeed the sums of the squared \(a_i\) of either orientation
have total at most \(CA\) and difference \(2P\). Bound
each signed weight-\(d\) contribution by \(A^{d-2}\) times
the squared \(a_i\) of its orientation. For the exponential
monomials, every power of a weight parameter is bounded as a
power of its chosen interval, including the square powers
of the leading signed moment. The other rays use the
unconditioned grouped enclosure.

## Continuation and exact replay

Keep the actual \(p,r\) fixed and multiply every logarithmic
correction to the Gaussian coefficient sequence by \(s\in[0,1]\).
In the conditioned enclosure this replaces \(P\) by \(sP\).
For odd weights the moment itself is multiplied by \(s\),
whereas its lower bound from the new \(P\) subinterval is
bounded by \(s^{d-1}\) times the lower bound associated to the
actual \(P\). The upper and signed moment bounds also hold
along the path: the sign-mass bounds follow from
\(s(CA\pm2P)\le CA\pm2sP\). The omitted logarithmic remainder
is multiplied by \(s\), so the same weight-nine remainder
estimate applies along this path. In particular every
intermediate polynomial is covered by the subinterval
containing \(sP\). This continuation does not require the
union of the conditioned enclosures to be convex.

At \(s=0\) the coefficients are Gaussian. Move the reference
ratio at fixed \(p\) to the comparison value
\(R(1-(R-1)p/2)\), with \(R=121/120\), and follow that
reference down to the anchor \(p=.24\). This path lies in the
ratio enclosures used for both sets of rays. The anchor sign
checks and the sector counts are the same as in the upper band,
with orders twelve and eleven alternating and with the negative
odd-order root lying in the final joined sector. The order-twelve
polynomial is strictly positive on the real line, so the two even
roots in that sector are nonreal.

The [exact checker](verify_multivariate_terminal_twelve_cayley_lower.py)
certifies 68,762 projection cells across 240 order, ray, and
subband checks. Each leading ray check covers all eight
leading-skew subintervals.

~~~sh
python3 research/factorial/image-rigidity/verify_multivariate_terminal_twelve_cayley_lower.py
python3 research/factorial/image-rigidity/verify_multivariate_terminal_twelve_cayley_lower.py 12 0 240 245
~~~
