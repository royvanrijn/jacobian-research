# Terminal-thirteen shared light bands by grouped root disks

Let \(u\) be the normalized order-twelve outer geometric mean of a
fully shared ordinary profile. If \(\ell\) and \(d\) are the lower
exponent and gap on a row, its largest order-twelve reciprocal is
\(d/(12\ell+1)\). The orders twelve and thirteen have disjoint root
arguments on the following bands, provided every reciprocal is at
most the stated cap.

| Band in \(u\) | Reciprocal cap | Accepted rectangles |
| --- | ---: | ---: |
| \([.65,.66]\) | \(1/40\) | 249 |
| \([.69,.70]\) | \(1/50\) | 128 |
| \([.70,.71]\) | \(1/50\) | 128 |
| \([.72,.73]\) | \(1/50\) | 155 |
| \([.73,.74]\) | \(1/50\) | 256 |
| \([.74,.75]\) | \(1/50\) | 248 |
| \([.75,.76]\) | \(1/50\) | 304 |

## Enclosure

Use the signed-skew and outer-mean coordinates of the
[shared light proof](MULTIVARIATE_TERMINAL_THIRTEEN_SHARED_LIGHT.md#matched-logarithmic-enclosure).
Reversing both exponent columns if needed makes the lower signed
coordinate nonnegative and preserves angular separation. The
coordinates lie in the product of the listed outer interval with
the closed signed interval \([0,1]\). The same enclosure computes
the upper signed coordinate by matched transport. It also applies
the rational upper bound for the outer-loss ratio from the feasible
integer rows at the lower endpoint of each rectangle.

The grouped cubic expansion gives an interval-polynomial sum
\(\sum_\nu I_\nu P_\nu(z)\) for each order, where a real parameter
in \(I_\nu\) multiplies the entire corresponding vector. The
checker converts the outward dyadic intervals to outward Arb balls
and keeps this grouping when translating a polynomial. An empty
feasible parameter rectangle may be discarded, using the checks
in the shared light enclosure. Otherwise the following tests
certify it or the algorithm bisects it.

## Root disks and angular order

At a proposed center \(c\), translate each group by \(z=c+w\)
before applying its scalar interval. Write the resulting interval
coefficients as \(a_k\), so that they enclose the actual Taylor
coefficients at \(c\). A disk test for radius \(r>0\) is
\[
 \inf|a_1|\,r>\sup|a_0|+\sum_{k=2}^{j}\sup|a_k|r^k.
\]
This proves exactly one root in the open disk by Rouché's
theorem. Approximate roots of a midpoint polynomial supply
candidate centers. The test uses their exact floating-point
values and outward enclosures of each Taylor coefficient.

For two successive upper disks with centers \(a,b\) and radii
\(r,s\), require the strict inequality
\[
 \operatorname{Im}(\overline a b)>|a|s+|b|r+rs.
\]
Every point in either disk lies in the upper half-plane, as verified
by a separate bound on its center. The inequality gives strict
angle ordering for every point in the first disk and every point
in the second. Tests on all adjacent pairs thus impose a strict
angular order on all the upper roots. Conjugation gives the lower
roots. An additional order-thirteen disk centered on the negative
real axis contains a unique root, so that root is real by
conjugation. These are distinct roots and exhaust both degree
counts. The two orders have disjoint arguments throughout each
accepted rectangle.

Replay all 1,468 accepted rectangles from the repository root with
Python and python-flint.

~~~sh
python3 research/factorial/image-rigidity/verify_multivariate_terminal_thirteen_shared_disks.py
~~~
