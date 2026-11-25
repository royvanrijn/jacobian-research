# Terminal-fourteen shared profiles from outer mean .50 to .64

Let \(u\) be the normalized order-thirteen reflected outer mean
of a fully shared ordinary exponent profile. The order-thirteen and
order-fourteen factorial moments have disjoint root arguments under
independent positive rescalings under either condition below.

| \(u\) | Cap on \(d/(13\ell+1)\) |
|---|---:|
| \([.50,.63]\) | \(1/40\) |
| \([.63,.64]\) | \(1/44\) |

Use the [grouped logarithmic enclosure](MULTIVARIATE_TERMINAL_THIRTEEN_SHARED_LIGHT.md#matched-logarithmic-enclosure)
with orders thirteen and fourteen and lower ratio \(R=169/168\).
The formulas are configured before the first order-dependent
cache evaluation. The integer row enumeration uses reciprocal
floor \(1/80\). For every row above this floor the product mean
bound limits the gap and smaller exponent. The finite list
bounds both the maximum legal reciprocal and the ratio of the
upper and lower outer losses on the current mean interval.
Below the floor the corresponding termwise bounds apply.

The [terminal-fourteen grouped disk test](MULTIVARIATE_TERMINAL_FOURTEEN_SHARED_DISKS.md)
uses the same matched skew parameter at both orders. It bounds
the shifted polynomial with outward complex arithmetic on each
parameter rectangle. The Rouché test isolates one root in
each proposed disk. Strict oriented-area inequalities order the
upper disks by argument. Conjugation gives the lower roots, and
the real disk at odd order gives one negative real root. The
one-root certificates exhaust both degrees.

The [checker](verify_multivariate_terminal_fourteen_shared_low_mid.py)
covers fourteen closed subintervals of length \(1/100\), with
3,663 accepted parameter rectangles. The original
polynomials are coprime over every characteristic-zero field
by rational Bézout.

~~~sh
.venv/bin/python research/factorial/image-rigidity/verify_multivariate_terminal_fourteen_shared_low_mid.py
~~~
