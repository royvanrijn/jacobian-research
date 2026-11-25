# Smaller reciprocal caps for terminal thirteen

For a fully shared ordinary profile, let \(u\) be the order-twelve
reflected outer mean after endpoint normalization. The largest
reciprocal is the maximum of \(d/(12\ell+1)\) over its unequal rows,
where \(\ell\) is the smaller exponent and \(d\) is the gap.

**Theorem.** The order-twelve and order-thirteen factorial moment
polynomials have disjoint root arguments under independent positive
rescalings on each listed band when the largest reciprocal is greater
than the lower cap and at most the upper cap. In particular the
polynomials are coprime over every characteristic-zero field.

| \(u\) band | Lower cap | Upper cap | Lists | Outer cells | Radial cells |
| --- | ---: | ---: | ---: | ---: | ---: |
| \([.53,.54]\) | \(1/30\) | \(1/22\) | 10 | 40 | 6,597 |
| \([.54,.55]\) | \(1/30\) | \(1/22\) | 10 | 40 | 6,952 |
| \([.56,.57]\) | \(1/40\) | \(1/20\) | 67 | 246 | 36,806 |
| \([.57,.60]\) | \(1/40\) | \(1/20\) | 49 | 618 | 103,132 |
| \([.60,.62]\) | \(1/40\) | \(1/30\) | 14 | 112 | 28,363 |
| \([.62,.63]\) | \(1/40\) | \(1/30\) | 13 | 52 | 13,359 |
| \([.66,.67]\) | \(1/50\) | \(1/30\) | 51 | 184 | 34,263 |
| \([.67,.68]\) | \(1/50\) | \(1/30\) | 38 | 146 | 28,465 |

Apply the [retained band enumeration](MULTIVARIATE_TERMINAL_THIRTEEN_RETAINED_BANDS.md)
with the upper and lower caps in the table. The retained list is
nonempty. Its square contribution to the outer mean is exact, and it
leaves a residual mean in the interval obtained by dividing the
outer-band endpoints by its mean and intersecting with \((0,1]\).
The [retained-profile proof](MULTIVARIATE_TERMINAL_THIRTEEN_RETAINED_PROFILES.md)
uses the refined model on these residual intervals to certify the
sectors at both orders.

Run the [checker](verify_multivariate_terminal_thirteen_retained_fallbacks.py)
from the repository root with Python 3, python-flint, and a C++17 compiler.

~~~sh
python3 research/factorial/image-rigidity/verify_multivariate_terminal_thirteen_retained_fallbacks.py
~~~
