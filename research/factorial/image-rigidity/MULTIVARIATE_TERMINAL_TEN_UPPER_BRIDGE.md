# Shared terminal-ten separation from outer mean .88 to .93

Let two ordinary exponent columns be distinct and have only shared
unequal rows. Let \(u\) be their order-nine outer geometric mean
after positive endpoint normalization and removal of the binomial
factors. If \(22/25\le u\le93/100\), their ninth and tenth
factorial moment polynomials have disjoint root arguments under
independent positive rescalings. The dimension and the positive
exponents are arbitrary. Rational Bézout gives coprimality in
characteristic zero.

The marked-unit-row theorem covers rows \((3,4)\), \((4,5)\)
and \((5,6)\) in either orientation on this band. Exclude those
rows for the remainder of the proof.
[Marked unit rows](MULTIVARIATE_TERMINAL_TEN_MARKED_UNIT_ROWS.md).

## Reciprocal and normalized moment bounds

Use \(p=1-u\), order-ten deficit \(q=1-w\), and \(r=q/p\).
An order-nine unit row with smaller exponent \(\ell\) contributes
\((9\ell+1)/(9\ell+9)\) to the square of \(u\). Thus \(\ell=1,2\)
are impossible here, and after the marked exclusions every unit
reciprocal is at most \(1/55\). For a nonunit row let
\(y=d/(9\ell)\). The squared row product bound gives
\(u\le((1+y)/(1+9y))^{d/2}\), so \(x_h\le y\le(1-u)/(9u-1)\).
At \(u\ge.88\) this gives \(x_h\le3/173<1/55\).

On a band of left endpoint \(u_-\), take the minimum of
\(1/55\), \((u_-^{-2}-1)/8\), and the maximum of the unit
and nonunit caps just proved. For the unit cap use
\(1/(9\ell_-+1)\), where \(\ell_-\) is the smallest integer
with \((9\ell_-+1)/(9\ell_-+9)\ge u_-^2\).
Write the resulting lower-order cap as \(c\).
The matched upper-order cap is \(9c/(10-c)\). The same outer
ratio comparison as in the near-flat proof gives
\[
  1+\frac{u(1+u)}{160(1+(1-u^2)/80)^{3/2}}
                  \le r\le\frac{81}{8(10-c)}.
\]
The checker verifies the rational range
\([5051/5000,5072/5000]\), the scale bound
\(\sqrt p\le3464102/10^7\), and \(q\le q_0=61/500\).
[Matched product comparison](MULTIVARIATE_TERMINAL_TEN_NEAR_FLAT.md#outer-deficits-and-common-angular-scale).

Use the weighted Taylor proof with the sharper normalized
tanh budgets supplied by these caps. In an order with deficit
\(h\), let \(h_-\) be a positive lower bound for the deficit on
the current subband. If \(T_*\) is the cap for that order's
tanh factors, put \(A=\min(C,T_*/h_-)\) and \(C=1/(1-q_0)\).
For order ten use \(h_-=r_-p_-\) from the lower ratio bound.
Then \(T_i\le A h\) and \(\sum_i T_i\le C h\).
For every \(d\ge2\) this gives
\[
               \sum_i T_i^d\le C A^{d-1}h^d .
\]
The same bound applies to signed moments in absolute value.
Thus replace each coefficient majorant \(C^d/d\) by
\(CA^{d-1}/d\) in the near-flat proof. Retain the same
six logarithmic directions and all products through weight
seven. The two omitted geometric tails use numerators
\(CA^7/8\) and \(CA^8q_0/9\), with denominator
\(1-(Aq_0)^2\), in place of the corresponding powers of
\(C\). The same rational exponential and Cayley tail bounds
then give the coefficient enclosure.
[Stable Taylor coefficient enclosure](MULTIVARIATE_TERMINAL_TEN_NEAR_FLAT.md#a-grouped-taylor-enclosure-through-deficit-weight-seven).

## Ray certification and continuation

The checker proves strict leading dominance in the same scaled
disk \(|Z|<2\) and excludes eight rays with piecewise affine
levels. The order-ten scale uses \(b/\sqrt r\).
The centered complex Taylor projection and the separate ratio
subdivision are the same as for the near-flat checker.
The disk maps into the left half-plane since \(4q_0<1\).

For each fixed \(p,r\), continue the actual arrays to their
Gaussian references by straight coefficient interpolation.
The affine enclosure in grouped scalars, products, and coefficient
remainders is convex and includes both arrays. Then move \(p\)
to \(.07\) through the reference arrays with the same \(r\).
These references belong to the enclosure throughout the band,
including the cap transitions. At \(p=.07\), the checker
evaluates the two reference polynomials on the scaled imaginary
axis and subdivides the ratio interval to certify the anchor signs.
They give nine upper sectors ordered
\(10,9,10,9,10,9,10,9,10\) and one additional real nonic root.
No root crosses the rays during the two continuations. The
nonic root in the last joined sector is real by conjugation
and the actual decic has no real zero by positivity of its
even factorial moment. This proves the claimed separation.

The exact replay certifies 203,802 projection cells over
160 order, ray, and subband checks. It verifies every expected
cell count in addition to the rational anchor bounds.
[Exact replay checker](verify_multivariate_terminal_ten_upper_bridge.py).

~~~sh
python3 research/factorial/image-rigidity/verify_multivariate_terminal_ten_upper_bridge.py
python3 research/factorial/image-rigidity/verify_multivariate_terminal_ten_upper_bridge.py 10 -1 70 90
~~~
