# Terminal-fourteen angular separation from outer mean .45 to .50

Let \(u\) be the normalized order-thirteen reflected outer mean
of a fully shared ordinary exponent profile with distinct columns.
For \(u\in[.45,.50]\), the order-thirteen and order-fourteen factorial
moment polynomials have disjoint root arguments under independent
positive rescalings.

Use the [middle-band direct and retained tests](MULTIVARIATE_TERMINAL_FOURTEEN_SHARED_MID.md)
with reciprocal cap \(1/30\). The direct certificate consists
of 2,640 rectangles on five closed subintervals of length \(1/100\).
For the retained complement with largest reciprocal at most \(1/8\),
the checker certifies 39 oriented lists modulo simultaneous reversal,
915 outer cells and 157,316 radial cells. The
[large-reciprocal theorem](MULTIVARIATE_TERMINAL_FOURTEEN_RETAINED_PROFILES.md)
covers profiles with a reciprocal greater than \(1/8\).
The adjacent
[middle](MULTIVARIATE_TERMINAL_FOURTEEN_SHARED_MID.md) and
[upper](MULTIVARIATE_TERMINAL_FOURTEEN_RETAINED_LOW_MID.md)
theorems extend the full interval to \([.30,1)\).
Rational Bézout gives coprimality over every characteristic-zero field.

~~~sh
.venv/bin/python research/factorial/image-rigidity/verify_multivariate_terminal_fourteen_shared_mid_flanks.py
~~~
