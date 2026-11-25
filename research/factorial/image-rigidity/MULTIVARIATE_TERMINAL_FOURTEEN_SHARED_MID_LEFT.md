# Terminal-fourteen angular separation from outer mean .25 to .30

Let \(u\) be the normalized order-thirteen reflected outer mean
of a fully shared ordinary exponent profile with distinct columns.
For \(u\in[.25,.30]\), the order-thirteen and order-fourteen factorial
moment polynomials have disjoint root arguments under independent
positive rescalings.

Use the [middle-band direct and retained estimates](MULTIVARIATE_TERMINAL_FOURTEEN_SHARED_MID.md).
The direct certificate uses cap \(1/18\) on \([.25,.28]\) and
cap \(1/20\) on \([.28,.30]\), producing 1,744 and 889 accepted
rectangles, respectively. The
[large-reciprocal theorem](MULTIVARIATE_TERMINAL_FOURTEEN_RETAINED_PROFILES.md)
handles profiles with a reciprocal above \(1/8\).

To connect these cases, on \([.25,.30]\) first retain every row
whose reciprocal exceeds \(1/18\), restricting to profiles whose
largest reciprocal is at most \(1/8\). This gives 31 nonempty
lists modulo simultaneous reversal, 790 outer cells and 111,384
radial cells. On \([.28,.30]\) the remaining profiles have every
reciprocal at most \(1/18\). Retain rows with reciprocal above
\(1/20\) there. This gives two more nonempty lists, 24 outer
cells and 3,380 radial cells. The direct bound includes equality
at its cap in either interval.

Together with the
[adjacent middle range](MULTIVARIATE_TERMINAL_FOURTEEN_SHARED_MID_FLANKS.md),
this gives separation for the full shared interval \([.25,1)\).
Rational Bézout gives coprimality over every characteristic-zero field.

~~~sh
.venv/bin/python research/factorial/image-rigidity/verify_multivariate_terminal_fourteen_shared_mid_left.py
~~~
