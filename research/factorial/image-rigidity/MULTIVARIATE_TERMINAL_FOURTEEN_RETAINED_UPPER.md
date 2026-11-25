# Terminal-fourteen upper retained complements

Let \(u\) be the order-thirteen reflected outer mean of a fully shared
ordinary profile, with endpoint normalization and binomial factors
removed. The largest reciprocal is the maximum of \(d/(13\ell+1)\)
over the unequal rows, with smaller exponent \(\ell\) and gap \(d\).

**Theorem.** The order-thirteen and order-fourteen factorial moment
polynomials have disjoint root arguments under independent positive
rescalings on each listed band when the largest reciprocal exceeds
the lower cap and is at most \(1/8\).

| \(u\) band | Lower cap | Retained lists | Outer mean cells |
| --- | ---: | ---: | ---: |
| \([.855,.88]\) | \(1/70\) | 3 | 145 |
| \([.88,.92]\) | \(1/80\) | 3 | 178 |
| \([.92,.95]\) | \(1/100\) | 2 | 48 |

In particular every fully shared profile with \(.855\le u<1\) has
disjoint root arguments. On each of the three bands, the
[high-mean theorem](MULTIVARIATE_TERMINAL_FOURTEEN_HIGH_U.md)
covers all reciprocals at or below the listed lower cap, and the
complement here covers those above it through \(1/8\). The
high-mean theorem also gives the uncapped result on \([.95,1)\).
The [retained-profile theorem](MULTIVARIATE_TERMINAL_FOURTEEN_RETAINED_PROFILES.md)
covers largest reciprocals above \(1/8\). Rational Bézout transfers
coprimality over the complex numbers to every characteristic-zero field.

## Exact row enumeration

The [retained band enumeration](MULTIVARIATE_TERMINAL_THIRTEEN_RETAINED_BANDS.md#retaining-rows-on-a-band)
applies with lower order thirteen. A row contributes the exact
factor
\[
 q(\ell,d)=\prod_{h=1}^d\frac{13\ell+h}{13\ell+12d+h}
\]
to \(u^2\). If its largest reciprocal exceeds \(c\), every reciprocal
on the row exceeds \(c/(1+c)\). Its square contribution is thus
strictly smaller than \((1+12c/(1+c))^{-d}\), which bounds \(d\)
when \(q(\ell,d)\) must reach the square of the lower endpoint.
The checker enumerates the remaining bases and multiplicities
and keeps one sorted list from each column reversal orbit.
It divides the outer band by the exact retained mean contribution
and brackets the resulting tail interval.

At \([.855,.88]\) the singleton row \((5,6)\) requires a second
retention at cap \(1/85\) within its tail cap \(1/70\). The exact
product enumeration gives the empty additional list and the row
\((6,7)\) in either orientation. The residual cap is refined with
the tail mean on each cell. Individual rows are enumerated above
reciprocal floor \(1/150\) in the first band and \(1/40\) elsewhere.
Reciprocals at or below the floor are covered by it.

## Matched skew disks

Use \(v\) for the order-thirteen tail mean and \(U=-\log v\). With a
tail reciprocal cap \(c\), put \(T_*=12c/(2+12c)\),
\(R=169/168\), and \(\eta=c/(28+167c)\). The first signed skew
parameter \(P\) is \(M_2/2\) up to the common orientation convention.
The inequality \(\sum_i T_i\le U\) gives \(P=\sigma T_*U/2\) for
\(\sigma\in[-1,1]\), with the zero-loss case included. Squaring the
[matched tanh ratio bounds](MULTIVARIATE_MATCHED_TANH_MOMENTS.md)
and centering their interval gives
\[
 P'=R^2(1+\eta^2)P+E,\qquad |E|\le 2R^2\eta\, T_*U/2.
\]
On a signed subinterval, center the first skew shape at the
corresponding midpoint estimate and use this transport enclosure
in the other order. The other log terms and their tails use the
absolute bounds in the
[retained root-disk proof](MULTIVARIATE_TERMINAL_THIRTEEN_RETAINED_PROFILES.md#direct-root-disks-with-a-small-residual-loss).
The ratio and cap at order fourteen obey the
[terminal-fourteen residual bounds](MULTIVARIATE_TERMINAL_FOURTEEN_RETAINED_PROFILES.md#retention-and-residual-bounds).

If \(L<1\) bounds the full logarithm and \(E_{\rm log}\) its excluded
terms, the quartic exponential remainder is bounded by
\((L^5/5!+E_{\rm log})/(1-L)\).
Each scalar monomial before this remainder multiplies a complete
coefficient direction. An included direction with a coefficientwise
absolute bound \(10^{-23}\) may instead be absorbed into the
independent error. The retained factorial multipliers stay exact.

Translate complete coefficient groups at every proposed center.
For the shifted polynomial \(\sum_k a_kz^k\) use a radius \(r\)
satisfying \(r\inf|a_1|>\sup|a_0|+\sum_{k\ge2}r^k\sup|a_k|\).
Rouché's theorem places one root in each verified disk. At adjacent
upper centers \(z,w\) with radii \(r,s\), the strict bound
\[
 \operatorname{Im}(\overline z w)>r|w|+s|z|+rs
\]
orders the arguments everywhere in those disks. The upper order
pattern is \(14,13\) repeated six times, followed by \(14\).
Conjugate disks give the lower roots. The remaining disk at order
thirteen has a negative real center and its unique root is real by
conjugation. These disks exhaust both degrees.

The [checker](verify_multivariate_terminal_fourteen_retained_upper.py)
uses Python 3 and python-flint and replays eight retained list
occurrences in 371 outer mean cells, subdividing the signed parameter
within a cell where needed.

~~~sh
python3 research/factorial/image-rigidity/verify_multivariate_terminal_fourteen_retained_upper.py
~~~
