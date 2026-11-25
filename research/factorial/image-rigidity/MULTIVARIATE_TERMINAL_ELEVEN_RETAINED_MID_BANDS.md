# Retained shared profiles in the terminal-eleven middle bands

Let two distinct ordinary exponent columns have all unequal rows
shared. Write \(u\) for their endpoint-normalized order-ten reflected
outer mean, the square root of the product of the multipliers
at indices one and nine after removing the binomial coefficients.
Suppose that
every reciprocal \(d/(10\ell+1)\) is at most \(1/8\). In any row of the
following table assume \(u\) lies in the specified band and at least
one reciprocal exceeds the listed cap.

| Outer mean \(u\) | Light cap \(c\) | Retained profiles up to column reversal |
| --- | ---: | --- |
| \([.44,.62]\) | \(9/100\) | One \((1,2)\), two \((1,2)\) in the same orientation, two \((1,2)\) in opposite orientations, or one \((2,4)\) |
| \([.62,.70]\) | \(3/40\) | One \((1,2)\) |
| \([.70,.76]\) | \(11/200\) | One \((1,2)\) |

Then the tenth and eleventh factorial moment polynomials have
disjoint root arguments under independent positive rescalings.
Consequently they are coprime over every characteristic-zero field.

## Exhaustion and residual ranges

Retain each row with largest reciprocal exceeding \(c\). A shared
row with smaller exponent \(\ell\), gap \(d\), and this property has
\(10\ell+1<d/c\) and squared outer mean
\[
 u_{\ell,d}^{\,2}
  =\prod_{h=1}^{d}\frac{10\ell+h}{10\ell+h+9d}
  <\left(1+\frac{9c}{1+c}\right)^{-d}.
\]
Indeed \(10\ell+h\le 10\ell+d<d(1+1/c)\). All row means are
at most one, so every retained row and every retained product has
mean at least \(u_-\), the band's lower endpoint. These inequalities
bound the candidate gaps and bases. The checker enumerates both
orientations and all multiplicities until the product of squared
means falls below \(u_-^2\). It identifies a signed list with its
global reversal. The exact lists are those in the table.

Write \(v_C\) for the retained product mean. The residual mean
\(v=u/v_C\) belongs to
\([u_-/v_C,\min(1,u_+/v_C)]\) and all its reciprocals are at most
\(c\). Comparing squares against \(v_C^2\) encloses the endpoints
outward on the rational grid with denominator 10,000. This
bracket and the fixed retained factorial multipliers are the
inputs to the [paired ray comparison](MULTIVARIATE_TERMINAL_ELEVEN_RETAINED_PROFILES.md#joint-coefficient-comparison).
Its coefficient bounds and adjacent attenuation apply here
because the residual cap is at most \(1/8\).

The comparison keeps the core fixed while the residual log
skew, mean defects, and Gaussian outer ratio move as in the
retained comparison. The checker bounds the annulus on this
residual band by the upper coefficient bounds including the
core. Outward complex balls isolate the midpoint roots, and the
adaptive rational rays have strict Bernstein projection signs
throughout each accepted mean cell. The order-eleven negative
real root and the last nonreal order-ten pair remain separated
as in the retained comparison.

## Reproduction

From the repository root run the
[exact retained complement checker](verify_multivariate_terminal_eleven_retained_mid_bands.py).

~~~sh
python3 research/factorial/image-rigidity/verify_multivariate_terminal_eleven_retained_mid_bands.py --native
~~~

The replay checks all six lists, with 207 outer cells and 28,506
radial cells. The three bands account for 20,423, 3,993, and
4,090 radial cells respectively.
