# Terminal-fourteen retained complements of the middle light band

Let \(u\) be the normalized order-thirteen reflected outer mean of a
fully shared ordinary profile. Write \(c_{\max}\) for the largest
reciprocal \(d/(13\ell+1)\) of its unequal rows, with smaller exponent
\(\ell\) and gap \(d\).

**Theorem.** On the bands in the table the order-thirteen and
order-fourteen factorial moment polynomials have disjoint root
arguments under independent positive rescalings whenever
\(1/60<c_{\max}\le1/8\).

| \(u\) band | Retained lists | Outer cells | Radial cells |
| --- | ---: | ---: | ---: |
| \([.64,.72]\) | 126 | 2,073 | 485,724 |
| \([.72,.74]\) | 29 | 219 | 66,515 |
| \([.74,.76]\) | 18 | 132 | 46,323 |
| \([.76,.81]\) | 13 | 218 | 65,795 |

Consequently, separation holds for every fully shared profile on
\([.64,.76]\). The [direct middle bound](MULTIVARIATE_TERMINAL_FOURTEEN_SHARED_DISKS.md)
covers \(c_{\max}\le1/60\), and the
[large-reciprocal theorem](MULTIVARIATE_TERMINAL_FOURTEEN_RETAINED_PROFILES.md)
covers \(c_{\max}>1/8\). Rational Bézout transfers coprimality
to every characteristic-zero field.

Retain every row above \(1/60\). Use the
[exact outer-mean product enumeration](MULTIVARIATE_TERMINAL_THIRTEEN_RETAINED_BANDS.md#retaining-rows-on-a-band)
at lower order thirteen to list the nonempty retained cores, up to
column reversal. Divide the band by the retained mean contribution and
bracket the remaining mean interval. Each residual cap is at most
\(1/60\). The exact row-cap refinement enumerates individual rows
above the floor \(1/150\) using the residual mean.

The [paired ray test](MULTIVARIATE_TERMINAL_FOURTEEN_RETAINED_PROFILES.md#paired-ray-checks)
applies to those residual intervals and their fixed cores. It uses
the refined matched skew enclosure, the same order-thirteen and
order-fourteen Gaussian references, and the exact endpoint attenuation
and annulus bounds. The orientation enumeration, all interval
endpoints, and the Bernstein projection are exact.

Run the [checker](verify_multivariate_terminal_fourteen_retained_middle.py)
from the repository root with Python 3, python-flint, and a C++17
compiler.

~~~sh
python3 research/factorial/image-rigidity/verify_multivariate_terminal_fourteen_retained_middle.py
~~~
