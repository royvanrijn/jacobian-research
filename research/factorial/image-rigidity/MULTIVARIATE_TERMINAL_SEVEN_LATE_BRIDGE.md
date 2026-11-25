# Shared terminal-seven separation from outer mean .75 to .88

Let \(a,b\) be distinct ordinary exponent columns whose unequal rows
are shared. Let \(u\) be the order-six outer geometric mean after
normalizing both endpoints to one and removing the binomial factors.
For \(3/4\le u\le22/25\), the order-six and order-seven polynomial
root arguments are disjoint under independent positive rescalings.
The dimension and exponents are unrestricted.

The [exact replay](verify_multivariate_terminal_seven_late_bridge.py)
uses the quadratic coefficient model of the
[upper Cayley bridge](MULTIVARIATE_TERMINAL_SEVEN_UPPER_BRIDGE.md).
The changes are a row bound, a tighter Cayley root disk, and a new
table of ray knots. The sectors are anchored at the same symmetric
arrays at \(u=22/25\).

## The row \((1,2)\)

A row with unordered exponents \((1,2)\) contributes \(7/12\) to
\(u^2\). There can be at most one such row since
\((7/12)^2<(3/4)^2\). Another row with an order-six reciprocal
exceeding \(1/8\) would contribute less than \(8/13\) to \(u^2\),
using just that Gauss factor. Hence a profile containing the marked
row has all remaining reciprocals at most \(1/8\), because
\((7/12)(8/13)<(3/4)^2\).

With \(S=\sum_i d_i^2/(\ell_i+d_i)\), the
[outer loss bound](MULTIVARIATE_TERMINAL_SEVEN_SHARED_HIGH_CURVATURE.md#a-row-restriction-from-the-total-outer-loss)
gives \(5S/12\le-\log u\). Thus \(u\ge3/4\) gives
\(S\le(12/5)\log(4/3)\le4/5<17/4\). The
[marked-row theorem](MULTIVARIATE_TERMINAL_SEVEN_MARKED_LIGHT_TAIL.md)
then proves the claim for profiles containing the row \((1,2)\)
in either orientation.

## Cap and ray model for the remaining rows

When there is no marked row, a unit-gap row has \(\ell\ge2\), so
its order-six reciprocal is at most \(1/13\).
For a gap \(d\ge2\), the
[upper bridge row bound](MULTIVARIATE_TERMINAL_SEVEN_UPPER_BRIDGE.md#ordinary-row-bounds)
gives \(d/(6\ell+h)\le(1-u)/(6u-1)\). For each interval take the
smaller of the upper bridge cap and
\(\max(1/13,(1-u_-)/(6u_--1))\), where \(u_-\) is the
interval's lower endpoint. This controls all its Gauss factors.
The checker supplies this cap to the existing Cayley coefficient
routine and leaves the linked linear and quadratic terms intact.

The new ray parameters \(b\) interpolate affinely between the
rational knot values. With \(p=1-u\) and \(\lambda^2=p\), the
physical ray has argument \(\pi-2\arctan(\lambda b)\).
The checker proves leading-term dominance at scaled radius
\(19/10\) for each coefficient box used in the ray tests.
Since \(p\le1/4\), this disk lies in the domain where the
physical Cayley image is in the open left half-plane.
For the real part \(x\) of the scaled variable it suffices to
check \(0\le x\le19/10\), by the reflection symmetry of the model.
The height is enclosed by the same equation and three interval
updates as in the upper bridge. Every accepted box has a
strictly signed projection; unresolved boxes are bisected.

At fixed \(u\) deform the actual arrays by a straight coefficient
segment to the symmetric arrays with order ratio \(36/35\).
Then increase \(u\) to \(22/25\) along those symmetric arrays.
The certified coefficient model contains these paths, and the
upper endpoint ray knots coincide with the upper bridge's
Gaussian anchor. The preserved upper sector counts occur in the
order \(7,6,7,6,7,6\). As in the upper bridge, joining the last
sector to its reflection isolates one real septic root and a
nonreal sextic pair. The latter is nonreal by the strict
real-axis positivity of the ordinary sixth moment. This proves
disjointness of arguments and rational Bézout gives coprimality.

The exact replay certifies 107,978 projection cells. Run it from the
repository root with

~~~sh
python3 research/factorial/image-rigidity/verify_multivariate_terminal_seven_late_bridge.py
~~~
