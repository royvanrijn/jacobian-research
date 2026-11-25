# Terminal-fourteen retained complements between outer means .50 and .64

Let \(u\) be the normalized order-thirteen reflected outer mean of a
fully shared ordinary profile. Write \(c_{\max}\) for the largest
reciprocal \(d/(13\ell+1)\) among unequal rows, with smaller exponent
\(\ell\) and gap \(d\).

**Theorem.** The moment polynomials at orders thirteen and fourteen
have disjoint root arguments under independent positive rescalings
on each of the following bands when the indicated reciprocal cap is
strictly less than \(c_{\max}\le 1/8\).

| \(u\) band | Lower reciprocal cap | Retained lists | Outer cells | Radial cells | Disk covers |
|---|---:|---:|---:|---:|---:|
| \([.50,.55]\) | \(1/40\) | 66 | 1,207 | 213,668 | 0 |
| \([.55,.60]\) | \(1/40\) | 32 | 553 | 116,072 | 0 |
| \([.60,.63]\) | \(1/40\) | 17 | 253 | 51,093 | 64 |
| \([.63,.64]\) | \(1/44\) | 28 | 116 | 28,151 | 0 |

The complementary ranges \(c_{\max}\) at or below the lower cap
and \(c_{\max}>1/8\) are covered, respectively, by the
[direct root-disk estimate](MULTIVARIATE_TERMINAL_FOURTEEN_SHARED_LOW_MID.md)
and the
[large-reciprocal theorem](MULTIVARIATE_TERMINAL_FOURTEEN_RETAINED_PROFILES.md).
Thus the theorem gives full separation throughout \([.50,.64]\).
Together with the [middle](MULTIVARIATE_TERMINAL_FOURTEEN_RETAINED_MIDDLE.md),
[bridge](MULTIVARIATE_TERMINAL_FOURTEEN_RETAINED_BRIDGE.md), and
[upper](MULTIVARIATE_TERMINAL_FOURTEEN_RETAINED_UPPER.md) certificates,
this covers every fully shared profile with \(u\in[.50,1)\).
Rational Bézout yields coprimality over every characteristic-zero field.

Retain every row above the lower cap. The
[exact outer-mean product enumeration](MULTIVARIATE_TERMINAL_THIRTEEN_RETAINED_BANDS.md#retaining-rows-on-a-band)
at lower order thirteen determines the retained oriented multisets,
modulo simultaneous reversal. For each list, divide the outer mean
by the retained contribution, bracket the resulting interval, and
use the lower reciprocal cap for the residual factors.
The exact residual row enumeration uses reciprocal floor \(1/150\).
The ray certificates use the
[paired ray test](MULTIVARIATE_TERMINAL_FOURTEEN_RETAINED_PROFILES.md#paired-ray-checks)
with the refined matched skew enclosure and the order-thirteen and
order-fourteen coefficients of the retained rows.

The single retained row \((1,2)\) in the band \([.60,.63]\)
leaves the residual mean in \([511/625,4293/5000]\). Apply the
[centered signed disk certificate](MULTIVARIATE_TERMINAL_FOURTEEN_RETAINED_BRIDGE.md#centered-signed-disks-for-the-remaining-list)
on the two subintervals separated by \(21/25\). There are 32
accepted cells on each side. The centered inequalities and order
transport depend on the residual profile and mean bounds. The
fixed factorial powers in this use of the certificate are those
of the row \((1,2)\).

The [checker](verify_multivariate_terminal_fourteen_retained_low_mid.py)
uses outward ball arithmetic for the reference and rational interval
bounds and integer Bernstein projection for the certified rays.

~~~sh
.venv/bin/python research/factorial/image-rigidity/verify_multivariate_terminal_fourteen_retained_low_mid.py
~~~
