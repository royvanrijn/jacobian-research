# Terminal-fourteen shared middle band by grouped root disks

**Theorem.** Let \(u\) be the normalized order-thirteen reflected
outer mean of a fully shared ordinary profile, and suppose
\(u\in[.64,.76]\). If every shared row with smaller exponent
\(\ell\) and gap \(d\) satisfies \(d/(13\ell+1)\le1/60\), then
the moment polynomials at orders thirteen and fourteen have
disjoint root arguments under independent positive rescalings.

Use the [grouped logarithmic enclosure](MULTIVARIATE_TERMINAL_THIRTEEN_SHARED_LIGHT.md#matched-logarithmic-enclosure)
with lower order \(13\), upper order \(14\), and lower ratio
\(R=169/168\). The formulas for the reciprocal cap, the exact
upper ratio and the higher signed moments are valid with these
orders substituted. Their estimates apply to the same integer
row factors at the two orders, as in the
[matched transport theorem](MULTIVARIATE_TANH_EXCESS_LOSS_TRANSPORT.md).
The [terminal-fourteen retained checker](verify_multivariate_terminal_fourteen_retained_profiles.py)
configures the common formulas before any order-dependent cache
is populated.

The [shared disk test](MULTIVARIATE_TERMINAL_THIRTEEN_SHARED_DISKS.md)
then applies on each parameter rectangle. Propose centers from
the two midpoint polynomials and translate the grouped interval
vectors with outward complex arithmetic. A Rouché inequality
certifies one root in each upper disk. Adjacent disks have
strictly ordered arguments by the oriented-area inequality.
Conjugate disks give the lower roots. There is one additional
disk for order thirteen on the negative real axis. Conjugation
and the one-root certificate make that root real. Counting the
disks exhausts the degrees at both orders.

The [checker](verify_multivariate_terminal_fourteen_shared_disks.py)
bisects the normalized signed coordinate and \(u\), separately
on twelve closed intervals of length \(1/100\). It verifies 2,306
accepted rectangles. The rational Bézout identity transfers
coprimality over the complex numbers to every characteristic-zero
field.

~~~sh
python3 research/factorial/image-rigidity/verify_multivariate_terminal_fourteen_shared_disks.py
~~~
