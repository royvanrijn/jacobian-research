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
| \([.12,.13]\) | \(1/10\) | 1 | 8 | 995 |
| \([.13,.15]\) | \(1/10\) | 1 | 16 | 2,064 |
| \([.15,.16]\) | \(1/11\) | 3 | 32 | 4,239 |
| \([.16,.17]\) | \(3/40\) | 4 | 48 | 6,474 |
| \([.17,.18]\) | \(3/40\) | 4 | 48 | 6,446 |
| \([.18,.20]\) | \(3/40\) | 4 | 96 | 12,805 |
| \([.20,.21]\) | \(3/40\) | 4 | 40 | 5,414 |
| \([.21,.22]\) | \(3/40\) | 3 | 24 | 3,240 |
| \([.22,.25]\) | \(1/14\) | 6 | 160 | 22,948 |

The [middle-band retained method](MULTIVARIATE_TERMINAL_FOURTEEN_SHARED_MID.md)
applies to each band. The exact outer products determine the retained
oriented lists, modulo simultaneous reversal, and bound the residual
mean intervals. The calculator keeps the factorial powers of all
retained rows and uses the refined matched-moment ray test at the
indicated cap for the residual rows. The ratio enumeration uses
reciprocal floor \(1/150\).

The table has 31 list occurrences, 476 outer cells and 65,114 radial
cells.

## The intermediate reciprocal strip at \([.21,.22]\)

For \(.21\le u\le.22\), the same test also separates every profile with
\(1/14<c_{\max}\le3/40\), where \(c_{\max}=\max d/(13\ell+1)\).
It retains all rows above \(1/14\) and uses residual cap \(1/14\).
The exact enumeration, up to simultaneous reversal, is as follows.
A triple gives the smaller exponent, the gap, and the orientation.

| Retained rows | Residual mean enclosure | Outer cells | Radial cells |
|---|---:|---:|---:|
| \((2,2,-1)\) | \([.3933,.4121]\) | 4 | 1,283 |
| \((2,2,-1),(2,2,-1)\) | \([.7366,.7718]\) | 8 | 1,045 |
| \((2,2,-1),(2,2,1)\) | \([.7366,.7718]\) | 8 | 1,054 |
| \((3,3,-1)\) | \([.5405,.5664]\) | 8 | 1,232 |

These are 28 outer and 4,614 radial cells. The exact row-product
bounds bracket the residual means outward. Only these four
lists meet the row-square and profile-square constraints in the strip.
Combining the strip with the table, the
[direct small-reciprocal check](MULTIVARIATE_TERMINAL_FOURTEEN_SHARED_LOW_U.md),
and the [large-reciprocal result](MULTIVARIATE_TERMINAL_FOURTEEN_RETAINED_PROFILES.md)
covers \([.21,.22]\) for every fully shared ordinary profile with
distinct columns.

The [checker](verify_multivariate_terminal_fourteen_retained_low.py)
accepts an optional lower endpoint in hundredths to replay a single
band. Endpoint 21 replays both reciprocal parts of this band.

~~~sh
.venv/bin/python research/factorial/image-rigidity/verify_multivariate_terminal_fourteen_retained_low.py
~~~
