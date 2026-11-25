# Terminal-seven separation with one marked (2,4) row

Let all unequal rows of two distinct ordinary exponent columns
be shared. Suppose one row is \((2,4)\) or \((4,2)\), and on
every other unequal row \(d/(6\ell+1)\le1/8\), where
\(\ell=\min(a_i,b_i)\) and \(d=|b_i-a_i|\). If
\[
              S=\sum_i\frac{(b_i-a_i)^2}{\max(a_i,b_i)}<17/4,
\]
the order-six and order-seven factorial moment polynomials
have disjoint root arguments up to independent positive
rescalings. Rational Bézout gives coprimality over every
characteristic-zero field.

Having exchanged the columns if necessary, take the marked row
to be \((2,4)\). It contributes \(1\) to \(S\), so the unmarked
curvature is below \(13/4\). The
[marked-tail curvature estimate](MULTIVARIATE_TERMINAL_SEVEN_MARKED_LIGHT_TAIL.md#residual-outer-mean-bounds)
then gives residual order-six outer mean
\(u_R>(7/12)^{169/48}>.149\). This is the bound for residual
curvature below \(13/4\), independent of the marked row.
The same residual reciprocal caps apply on \([.149,.620]\),
\([.620,.662]\), \([.662,.900]\), and \([.900,1]\), with values
\(1/8,2/19,2/25,1/20\), respectively.

The normalized order-\(j\) multiplier of the marked row is
\[
       \frac{(2j+2k)!}{((2j)!)^{\,1-k/j}((4j)!)^{\,k/j}}
\]
at index \(k\). The [exact replay](verify_multivariate_terminal_seven_marked_24_tail.py)
obtains it by integer comparisons for a rational \(j\)-th root.
It multiplies each nominal and perturbation coefficient of
the linked residual array by this same factor. In the code
the existing marked-tail coefficient models already include
one \((1,2)\) weight, so multiplication by the ratio of the
two positive row factors replaces that weight. Interval
operations round outward. All subsequent Bernstein
projection inequalities are the same as in the
[paired terminal-seven proof](MULTIVARIATE_TERMINAL_SEVEN_PAIRED_LIGHT.md).
The new rational ray knots and their exact projection checks
cover the whole residual range \([.149,1]\).

At fixed \(u_R,\rho\), deform to the marked array times the
symmetric residual with those means, then vary \(\rho\) to
\(36/35\) and \(u_R\) to one. The linked projections preserve
the sector counts throughout these paths. The limiting
polynomials are the normalized moments of the sole marked
\((2,4)\) row. At that anchor, the checker isolates their
roots in rational disks using the same Rouché inequality
and signed-distance tests as in the
[marked \((1,2)\) proof](MULTIVARIATE_TERMINAL_SEVEN_MARKED_LIGHT_TAIL.md#projections-and-pure-marked-anchors).
The six upper sectors have roots of orders \(7,6,7,6,7,6\).
The sole other septic root is real and negative. The
constant and leading terms dominate at radii \(1/10,10\)
throughout the homotopies by \((11/10)^7<2\) and the unit
upper bound on the interior multipliers. In the final
sector joined with its conjugate the actual sextic pair
is nonreal by real-axis positivity of the sixth factorial
moment, while the single septic root is real by conjugacy.
The arguments are therefore disjoint at the original
profile.

Replay from the repository root with

~~~sh
python3 -O research/factorial/image-rigidity/verify_multivariate_terminal_seven_marked_24_tail.py
~~~
