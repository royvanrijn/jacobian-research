# Retained complements for the terminal-thirteen lower bridge

For a fully shared ordinary profile, use \(u\) for the order-twelve
reflected outer mean after endpoint normalization. The largest
reciprocal is the maximum of \(d/(12\ell+1)\) over the unequal rows,
where \(\ell\) is the smaller exponent and \(d\) is the gap.

**Theorem.** At each listed band, the order-twelve and order-thirteen
factorial moment polynomials have disjoint root arguments under
independent positive rescalings if the largest reciprocal exceeds
the lower cap and is at most the upper cap. Their coprimality
follows over every characteristic-zero field.

| \(u\) band | Lower cap | Upper cap | Lists | Outer cells | Radial cells |
| --- | ---: | ---: | ---: | ---: | ---: |
| \([.63,.64]\) | \(1/40\) | \(1/30\) | 13 | 52 | 12,961 |
| \([.64,.66]\) | \(1/40\) | \(1/30\) | 13 | 92 | 21,807 |
| \([.63,.66]\) | \(1/50\) | \(1/40\) | 38 | 429 | 78,844 |
| \([.68,.72]\) | \(1/50\) | \(1/36\) | 32 | 448 | 106,829 |
| \([.72,.76]\) | \(1/50\) | \(1/34\) | 16 | 216 | 62,727 |
| \([.68,.72]\) | \(1/60\) | \(1/50\) | 7 | 68 | 13,580 |
| \([.72,.76]\) | \(1/60\) | \(1/50\) | 3 | 48 | 14,052 |

Apply the [retained band enumeration](MULTIVARIATE_TERMINAL_THIRTEEN_RETAINED_BANDS.md)
with the displayed caps. Each list retains every row above the lower
cap, so it is nonempty. The product of their exact row squares gives
the retained contribution to \(u^2\). The residual mean interval is
obtained by dividing the band endpoints by the retained mean and
intersecting with \((0,1]\). The
[retained-profile proof](MULTIVARIATE_TERMINAL_THIRTEEN_RETAINED_PROFILES.md)
uses the refined model to verify the sectors at both orders on those
intervals. Column reversal reduces the lists as in the retained band
enumeration.

The table joins the preceding retained proofs at cap \(1/30\) on
\([.63,.66]\), cap \(1/36\) on \([.68,.72]\), and cap \(1/34\) on
\([.72,.76]\). A direct result under the lower cap covers the
remaining profiles in the corresponding band.

The [checker](verify_multivariate_terminal_thirteen_retained_bridges.py)
uses the [smaller-cap verifier](verify_multivariate_terminal_thirteen_retained_fallbacks.py)
with the bridge band table. Run it from the repository root with
Python 3, python-flint, and a C++17 compiler.

~~~sh
python3 research/factorial/image-rigidity/verify_multivariate_terminal_thirteen_retained_bridges.py
~~~

The replay covers 122 list occurrences, 1,353 outer cells, and
310,800 radial cells.
