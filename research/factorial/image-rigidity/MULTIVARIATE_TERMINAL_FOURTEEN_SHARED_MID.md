# Terminal-fourteen angular separation from outer mean .30 to .45

Let \(u\) be the normalized order-thirteen reflected outer mean
of a fully shared ordinary exponent profile with distinct columns.
For \(u\in[.30,.45]\), its order-thirteen and order-fourteen factorial
moment polynomials have disjoint root arguments under independent
positive rescalings. In particular, they are coprime over every
characteristic-zero field by rational Bézout.

For an unequal row of smaller exponent \(\ell\) and gap \(d\),
write \(x=d/(13\ell+1)\). The direct test handles profiles with every row at or below the
indicated cap. The retained test covers profiles with some
\(x\) greater than that cap and all \(x\le1/8\).
The [large-reciprocal theorem](MULTIVARIATE_TERMINAL_FOURTEEN_RETAINED_PROFILES.md)
covers the remaining profiles.

| \(u\) band | Direct cap | Direct rectangles | Retained lists | Outer ray cells | Radial ray cells |
|---|---:|---:|---:|---:|---:|
| \([.30,.35]\) | \(1/22\) | 1,696 | 25 | 768 | 115,196 |
| \([.35,.40]\) | \(1/25\) | 1,700 | 21 | 616 | 98,164 |
| \([.40,.45]\) | \(1/25\) | 1,827 | 10 | 294 | 51,701 |

The direct test uses the same grouped logarithmic and
[Rouché disk estimates](MULTIVARIATE_TERMINAL_FOURTEEN_SHARED_LOW_MID.md)
at lower order thirteen and upper order fourteen. The integer
row enumeration uses floor \(1/80\) and is recomputed on each
closed subinterval of length \(1/100\). The disk comparison
matches the signed second moment across the orders and retains
the higher directions in the coefficient groupings.

For each retained list keep the exact factorial powers of
every row above the cap and bracket the residual mean by exact
outer products, as in the
[retained-band enumeration](MULTIVARIATE_TERMINAL_THIRTEEN_RETAINED_BANDS.md#retaining-rows-on-a-band).
The residual integer row enumeration uses floor \(1/150\).
The [paired ray test](MULTIVARIATE_TERMINAL_FOURTEEN_RETAINED_PROFILES.md#paired-ray-checks)
uses the refined matched skew enclosure on each retained list.
Simultaneous column reversal accounts for the other orientations.

The [checker](verify_multivariate_terminal_fourteen_shared_mid.py)
verifies 5,223 direct rectangles and 56 retained list occurrences,
with 1,678 outer and 265,061 radial retained ray cells.

~~~sh
.venv/bin/python research/factorial/image-rigidity/verify_multivariate_terminal_fourteen_shared_mid.py
~~~

A positional value \(30\), \(35\), or \(40\) selects a band.
The optional argument `--part direct` or `--part retained`
selects one certificate.
