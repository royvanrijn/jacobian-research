# Terminal-thirteen separation above shared reciprocal caps

Let \(u\) be the order-twelve outer geometric mean of a fully shared
ordinary exponent profile, after dividing the factorial multipliers by
their endpoint chords. On a row with smaller exponent \(\ell\) and gap
\(d\), the largest order-twelve reciprocal is \(d/(12\ell+1)\).

**Theorem.** The order-twelve and order-thirteen factorial moment
polynomials have disjoint root arguments under independent positive
rescalings on each listed band whenever the largest reciprocal is
greater than the lower cap and at most the upper cap.

| \(u\) band | Lower cap | Upper cap | Lists | Outer cells | Radial cells | Disk covers |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| \([.18,.27]\) | \(3/40\) | \(1/8\) | 42 | 2,176 | 270,329 | 0 |
| \([.27,.43]\) | \(3/50\) | \(1/8\) | 25 | 1,216 | 167,621 | 0 |
| \([.43,.60]\) | \(1/20\) | \(1/8\) | 8 | 706 | 105,157 | 0 |
| \([.47,.50]\) | \(1/22\) | \(1/20\) | 1 | 16 | 2,212 | 0 |
| \([.52,.55]\) | \(1/22\) | \(1/20\) | 1 | 4 | 571 | 0 |
| \([.60,.66]\) | \(1/30\) | \(1/8\) | 9 | 200 | 48,818 | 0 |
| \([.66,.72]\) | \(1/30\) | \(1/8\) | 5 | 120 | 26,080 | 0 |
| \([.66,.72]\) | \(1/36\) | \(1/30\) | 1 | 32 | 5,565 | 0 |
| \([.72,.76]\) | \(1/34\) | \(1/8\) | 3 | 47 | 13,208 | 0 |
| \([.76,.855]\) | \(1/55\) | \(1/8\) | 13 | 289 | 82,932 | 0 |
| \([.855,.88]\) | \(1/60\) | \(1/8\) | 2 | 24 | 5,799 | 0 |
| \([.88,.92]\) | \(1/70\) | \(1/8\) | 2 | 48 | 11,863 | 0 |
| \([.92,.95]\) | \(1/90\) | \(1/8\) | 3 | 71 | 2,137 | 65 |

Rational Bézout gives coprimality in characteristic zero.


## Retaining rows on a band

For a positive lower outer mean \(u_-\), retain each row above a
positive cap \(c\). Restrict the largest reciprocal to an upper cap
\(C\). The row contributes the exact factor
\[
 q(\ell,d)=\prod_{h=1}^d\frac{12\ell+h}{12\ell+11d+h}
\]
to \(u^2\). If \(d/(12\ell+1)>c\), every Gauss reciprocal on that
row exceeds \(c/(1+c)\). Hence
\[
 q(\ell,d)<(1+11c/(1+c))^{-d}.
\]
Since \(q(\ell,d)\ge u_-^2\) is necessary, the inequality gives a
finite gap bound. At a fixed gap, the same reciprocal threshold
bounds \(\ell\). The checker checks each candidate's exact square,
then enumerates the oriented multisets whose products are at least
\(u_-^2\). All factors are less than one, so each multiplicity is
bounded. The empty list is omitted because this is the complement
above the cap \(c\). Sort candidates lexicographically before comparing
each retained multiset with its column reversal. Every orbit has
one retained representative.

For a fixed list with square product \(q\) and an upper outer mean
\(u_+\), the residual mean \(v\) lies in
\([u_-/\sqrt q,\min(1,u_+/\sqrt q)]\). The checker brackets both
endpoints by squared rational comparisons. Its residual cap is \(c\),
with further reductions from the residual mean on each cell.
It then applies the independent, grouped, or direct-disk tests of
the [retained-profile proof](MULTIVARIATE_TERMINAL_THIRTEEN_RETAINED_PROFILES.md)
with exactly that core. The residual mean ranges and reciprocal cap
suffice for the coefficient bounds and the argument separation
steps. A total strength budget is unnecessary for these band tests.

For any additional retention step the same enumeration applies to
the residual mean relative to the first core. Both orientations of
every new row are retained because the orientation of the first core
has already been chosen. The empty additional list is included with
the new lower cap, and each nonempty list uses a fresh residual
bracket and that same cap. The \([.43,.60]\) band applies this
second retention to the singleton unit row \((1,2)\). It lowers the
residual cap from \(1/20\) to \(1/30\), and verifies both orientations
of any additional rows above \(1/30\).


## Replay

Run the [checker](verify_multivariate_terminal_thirteen_retained_bands.py)
with Python 3, python-flint, and a C++17 compiler. For example, the
highest band is

~~~sh
python3 research/factorial/image-rigidity/verify_multivariate_terminal_thirteen_retained_bands.py --native cayley-upper-tight
~~~

The disk covers are included among the outer cells. The radial
count covers the outer cells using the ray and annulus test.
