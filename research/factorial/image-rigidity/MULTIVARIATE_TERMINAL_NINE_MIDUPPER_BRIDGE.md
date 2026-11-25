# Shared terminal-nine angular separation from outer mean .76 to .88

Let two distinct ordinary exponent columns have only shared unequal
rows. Let \(u\) be the order-eight outer geometric mean after positive
endpoint normalization and removal of the binomial factors. For
\(19/25\le u\le22/25\), the order-eight and order-nine polynomials
have disjoint root arguments under independent positive rescalings.
The dimension and positive exponents are arbitrary. Rational Bézout
gives coprimality in characteristic zero.

Profiles containing the unit row \((2,3)\) in either orientation
are covered throughout this range. When \(u\ge21/25\), profiles
containing the unit row \((3,4)\) are also covered. The remainder of
the proof handles the profiles without those rows in the respective
ranges.
[Marked (2,3) theorem](MULTIVARIATE_TERMINAL_NINE_MARKED_TWO_THREE.md).
[Marked (3,4) theorem](MULTIVARIATE_TERMINAL_NINE_UPPER_MARKED_ROWS.md).

## Ordinary reciprocal caps

For a row with lower exponent \(\ell\) and gap \(d\), put
\(x_h=d/(8\ell+h)\), \(1\le h\le d\), and \(y=d/(8\ell)\).
The squared outer mean of the row gives
\[
 u^2\le \prod_{h=1}^d \frac{8\ell+h}{8\ell+7d+h}
        \le \left(\frac{1+y}{1+8y}\right)^d.
\]
For \(d=1\), the row \((1,2)\) is excluded by \(u>3/4\).
After the marked-row exclusions the unit reciprocals are at most
\(1/25\) below \(.84\) and at most \(1/33\) from \(.84\) onward.

For \(d\ge3\) and \(x_1>2/49\), we have \(y>2/49\) and thus
\(u^2<(51/65)^3<(19/25)^2\). For \(d=2\), the condition
\(x_1>2/49\) forces \(\ell\le5\). The product then gives
\(u^2\le41\cdot42/(55\cdot56)<(19/25)^2\).
This proves the \(2/49\) cap below \(.78\).

For \(d\ge3\) and \(x_1>1/25\), the same bound gives
\(u^2<(26/33)^3<(39/50)^2\). For \(d=2\), the condition
\(x_1>1/25\) forces \(\ell\le6\) and gives
\(u^2\le49\cdot50/(63\cdot64)<(39/50)^2\).
This proves the \(1/25\) cap from \(.78\) to \(.84\).

For \(u\ge21/25\) and \(d\ge2\), the general row bound gives
\[
        x_h\le y\le\frac{1-u}{8u-1}
                     \le\frac{1-21/25}{8(21/25)-1}<1/33.
\]
The checker uses exactly these three caps on the order-eight
reciprocals. For an order-eight cap \(c\), the linked order-nine
cap is \(8c/(9-c)\).

## Coefficient enclosure and ray continuity

At each midpoint \(u_0\) the two reference multiplier arrays are
\(u_0^{e_{8,k}}\) and \(u_0^{(64/63)e_{9,k}}\), where
\(e_{j,k}=k(j-k)/(j-1)\). Use the linked outer logarithm ratio
and the six grouped logarithmic directions. The coefficient enclosure
retains their linear and quadratic contributions and bounds the
remaining terms coefficientwise. For \(p=1-u=\lambda^2\),
apply the Cayley map \(t=(\lambda X-1)/(\lambda X+1)\).
For each direction, combine reflected pairs in the integer
Cayley matrix before introducing the symmetric reference
coefficient intervals. The interval projection and binomial quotient
anchor are those of the upper bridge.
[Grouped logarithmic coefficient model](MULTIVARIATE_TERMINAL_NINE_PRIVATE_LINEAR.md#grouping-the-shared-logarithmic-corrections).
[Exact Cayley method and anchor](MULTIVARIATE_TERMINAL_NINE_UPPER_BRIDGE.md#interval-coefficient-model-and-ray-exclusion).

The checker proves leading-term dominance on \(|X|=2\) and strict
exclusion of the seven shared rays. Since \(4(1-u)<1\), the
image of the root disk lies in the open left half-plane.
The real part of \(X\) need only be tested on \([0,2]\) because
reflection \(X\mapsto-\bar X\) reverses the coefficient arrays
and preserves the enclosure. For each ray level \(b\), the
height \(y\) at real coordinate \(x\) satisfies
\((y-b)(1+pby)=-pbx^2\). Three interval updates enclose the height.
Continue the actual arrays along their straight coefficient
segments to the reference arrays at the same \(u\), and continue
the reference arrays along the band to \(u=.76\). The convex
coefficient enclosure and the continuous ray table certify this
deformation even where the cap changes between subintervals.

At the anchor, the reciprocal quotients in \(Z=t+t^{-1}\)
have strict alternating signs at the ray abscissae and \(Z=\pm2\).
These place the upper roots in sectors of orders
\(9,8,9,8,9,8,9,8\), with one additional real nonic root.
The counts persist along the deformation. Join the last upper
sector to its lower reflection across the negative real axis.
This joined sector contains two octic roots and the sole real
nonic root. The nonic root stays real by conjugation. The
actual octic has no real root because it is the eighth
factorial moment of a nonzero real linear combination of distinct
monomials and that even moment is strictly positive.
Thus the arguments are disjoint throughout the band.

The exact replays certify 168 order, ray, and subband checks
with 3,785,361 projection cells. The two accelerated nonic
subbands \([.83,.84]\) and \([.87,.88]\) use tighter exact
power brackets, additional monomial-normal projections, and a
different subdivision ratio. The remaining checks use the
original enclosures and subdivision ratio. The checker verifies
the expected cell count for every replayed ray.
[Exact replay checker](verify_multivariate_terminal_nine_midupper_bridge.py).

~~~sh
python3 research/factorial/image-rigidity/verify_multivariate_terminal_nine_midupper_bridge.py
python3 research/factorial/image-rigidity/verify_multivariate_terminal_nine_midupper_bridge.py 830 840 9
~~~
