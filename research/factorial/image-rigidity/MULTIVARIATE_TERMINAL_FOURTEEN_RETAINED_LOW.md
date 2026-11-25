# Terminal-fourteen retained complements at low outer means

Let \(u\) be the normalized order-thirteen reflected outer mean of
a distinct fully shared ordinary exponent profile. On the listed
intervals the order-thirteen and order-fourteen factorial moment
polynomials have disjoint root arguments under independent positive
rescalings if the largest lower reciprocal \(d/(13\ell+1)\) exceeds
the indicated cap and is at most \(1/8\).

| \(u\) band | Reciprocal cap | Retained lists | Outer cells | Radial cells |
|---|---:|---:|---:|---:|
| \([.10,.11]\) | \(1/9\) | 1 | 4 | 489 |
| \([.13,.15]\) | \(1/10\) | 1 | 16 | 2,064 |
| \([.15,.16]\) | \(1/11\) | 3 | 32 | 4,239 |
| \([.16,.17]\) | \(3/40\) | 4 | 48 | 6,474 |
| \([.18,.20]\) | \(3/40\) | 4 | 96 | 12,805 |

The [middle-band retained method](MULTIVARIATE_TERMINAL_FOURTEEN_SHARED_MID.md)
applies to each band. The exact outer products determine the retained
oriented lists, modulo simultaneous reversal, and bound the residual
mean intervals. The calculator keeps the factorial powers of all
retained rows and uses the refined matched-moment ray test at the
indicated cap for the residual rows. The ratio enumeration uses
reciprocal floor \(1/150\).

There are 13 list occurrences, 196 outer cells and 26,071 radial
cells. The [checker](verify_multivariate_terminal_fourteen_retained_low.py)
accepts an optional lower endpoint in hundredths to replay a single
band.

~~~sh
.venv/bin/python research/factorial/image-rigidity/verify_multivariate_terminal_fourteen_retained_low.py
~~~
