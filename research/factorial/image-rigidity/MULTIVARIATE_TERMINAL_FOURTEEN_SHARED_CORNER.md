# Terminal-fourteen separation at outer means from .02 to .04

Let \(u\) be the normalized order-thirteen reflected outer mean of
a fully shared ordinary exponent profile with distinct columns.
For \(u\in[.02,.04]\), the order-thirteen and order-fourteen
factorial moment polynomials have disjoint root arguments under
independent positive rescalings. In particular they are coprime
over every characteristic-zero field.

For an unequal row with smaller exponent \(\ell\) and gap \(d\),
write \(x=d/(13\ell+1)\). The proof splits at \(x=2/25\).
Profiles with a row above \(1/8\) are handled by the
[large-reciprocal result](MULTIVARIATE_TERMINAL_FOURTEEN_RETAINED_PROFILES.md).

## Profiles with all reciprocals at most \(2/25\)

The first three bands use the
[paired ray check](MULTIVARIATE_TERMINAL_FOURTEEN_RETAINED_PROFILES.md#paired-ray-checks)
with an empty retained core. The rays separate the Gaussian
anchors and the refined matched envelopes certify that these
same rays avoid every root through the homotopy, with the
annulus bounds checked on each band. The empty core has both
endpoint factors one.

| Mean band | Outer ray cells | Radial ray cells |
|---|---:|---:|
| \([.020,.024]\) | 1 | 111 |
| \([.024,.028]\) | 1 | 116 |
| \([.028,.032]\) | 1 | 124 |

The last two bands use the
[grouped projected Rouché comparison](MULTIVARIATE_TERMINAL_FOURTEEN_SHARED_LOW_U.md#angular-bounds-from-grouped-root-disks).
The declared cap remains \(2/25\). Enumeration raises the
integer reciprocal floor only for bounding outer ratios as in
that comparison.

| Mean band | Grouped root rectangles |
|---|---:|
| \([.032,.036]\) | 18 |
| \([.036,.040]\) | 9 |

## Profiles with a reciprocal in \((2/25,1/8]\)

Retain every row above \(2/25\). The
[exact row-square enumeration and tail bracketing](MULTIVARIATE_TERMINAL_THIRTEEN_RETAINED_BANDS.md#retaining-rows-on-a-band)
applied at lower order thirteen gives 43 distinct nonempty
oriented lists up to simultaneous reversal. In the enumeration
every product of squared retained-row means must be at least
\((.02)^2\). The residual means are outwardly bracketed using
the exact squared retained products and the endpoints .02 and
.04. The remaining rows all have reciprocals at most \(2/25\).

Keep all factorial coefficient factors of these lists. Use the
same refined paired ray check with each list and its full
residual mean band. The residual ordinary ratio bounds use
the reciprocal floor \(1/150\). The replay verifies 66 outer
and 6,903 radial cells for the 43 lists.

## Exact replay

The [checker](verify_multivariate_terminal_fourteen_shared_corner.py)
binds the shared and retained models to orders thirteen and
fourteen before populating their order-dependent caches. It
replays 69 outer and 7,254 radial cells and 27 grouped root
rectangles. It uses outward ball arithmetic for the root
comparisons and the existing exact native Bernstein projector.

~~~sh
.venv/bin/python research/factorial/image-rigidity/verify_multivariate_terminal_fourteen_shared_corner.py
~~~
