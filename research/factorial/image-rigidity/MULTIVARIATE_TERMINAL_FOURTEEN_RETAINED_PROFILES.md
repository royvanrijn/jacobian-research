# Terminal-fourteen retained profiles

Let \(a,b\) be distinct ordinary exponent columns, with equal rows
omitted. A shared unequal row with smaller exponent \(\ell\) and gap
\(d\) has strength \(d^2/(\ell+d)\) and largest order-thirteen
reciprocal \(d/(13\ell+1)\).

**Theorem.** The factorial moment polynomials at orders thirteen and
fourteen have disjoint root arguments under independent positive
rescalings if there is a private row, or if all unequal rows are
shared and some reciprocal \(d/(13\ell+1)\) exceeds \(1/8\). Their
coprimality follows over every characteristic-zero field.

## Retention and residual bounds

Use the shared-strength cutoffs in the
[terminal-fourteen fixed-ray theorem](MULTIVARIATE_TERMINAL_FOURTEEN_HIGH_CURVATURE.md).
That theorem also covers every total private degree at least four.
Below a cutoff, retain all private rows and all shared rows above the
reciprocal cap \(1/8\). In a strength budget \(B\), a row above a cap
\(c\) satisfies \(\ell<(d/c-1)/13\), hence \(d<B(1+1/(13c))\).
The checker enumerates the finitely many permitted rows with their
multiplicities. If the private gap multisets are equal, it compares
each sorted list with its reversal and keeps one representative.
Otherwise it keeps both shared orientations for the fixed private
orientation. The fully shared family requires a nonempty retained
list.

Let \(B\) now denote the strength budget left after retaining the
rows. Write \(F_{j,k}\) for the normalized residual factorial
multiplier, with endpoints one and the binomial factor removed.
Let \(v=\sqrt{F_{13,1}F_{13,12}}\). As in the
[terminal-thirteen residual mean bound](MULTIVARIATE_TERMINAL_THIRTEEN_RETAINED_PROFILES.md#residual-bounds),
a shared residual row with largest reciprocal \(x\) contributes at
most \((d/2)\log(1+12x)\) to \(-\log v\). With \(\sigma=d^2/(\ell+d)\)
one has \(d/\sigma<1+1/(13x)\). The expression
\((1+1/(13x))\log(1+12x)\) increases with \(x\). On putting \(z=12x\),
the derivative has the sign of
\(z(1+13z/12)-(1+z)\log(1+z)\), whose derivative is
\(13z/6-\log(1+z)>0\) for \(z>0\). At \(x\le 1/8\) this gives
\[
 -\log v\le \frac{21}{26}\log(5/2) B .
\]
Thus the initial mean interval extends from an outward rounding of
\((2/5)^{21B/26}\) to \(1\). At a cell's lower endpoint \(v_-\),
the product identity for \(v^{-2}\) bounds each reciprocal by
\((v_-^{-2}-1)/12\). The checker refines this bound by enumerating
admissible individual rows, using their exact mean products.

The [consecutive endpoint attenuation](MULTIVARIATE_ORDINARY_ADJACENT_COEFFICIENT_ATTENUATION.md#an-all-order-endpoint-inequality)
with \(M\ge2\) gives multiplier \((377/378)^{k^2S_{\rm res}}\)
at distance \(k\) from either endpoint. The exact tests
\(\log(378/377)\ge1/378\) and
\(378(21/26)\log(5/2)<280\) therefore give
\[
 F_{14,k}\le F_{13,k}v^{k^2/280},\qquad
 F_{14,k}\le F_{13,k-1}v^{(14-k)^2/280}
 \qquad(1\le k\le13).
\]
For the linked Gaussian array, let \(e_{j,k}=k(j-k)/(j-1)\).
As the lower reciprocal \(x\) increases to \(1/8\), the lower
outer-loss quotient
\(\log(1+(169/14)x)/\log(1+12x)\) decreases by concavity
of \(y\mapsto\log(1+(169/168)(e^y-1))\). Thus the lower quotient
is at least \(\log(281/112)/\log(5/2)>1+1/280\), with the
last strict inequality verified by rational logarithm bounds.
The Gaussian endpoint inequalities hold since
\[
 (1+1/280)e_{14,k}-e_{13,k}-k^2/280
 =k(k-1)\left(\frac1{13\cdot12}-\frac{14}{13\cdot280}\right)\ge0 .
\]
The reference ratio \(169/168\) and the residual interpolation
use the same [matched moment
enclosures](MULTIVARIATE_TERMINAL_THIRTEEN_RETAINED_PROFILES.md#residual-bounds)
as the terminal-thirteen proof.

## Empty one-private list

When the only private row has gap one and no common row was retained,
the checker subdivides the initial residual mean range as follows.
On the middle bands it retains the additional shared rows with
\(1/20<d/(13\ell+1)\le1/8\).

| \(v\) band | Residual cap | Lists | Outer cells | Radial cells |
| --- | ---: | ---: | ---: | ---: |
| Initial lower endpoint to \(.30\) | \(1/8\) | 1 | 48 | 6,187 |
| \([.30,.40]\) | \(1/20\) | 28 | 731 | 94,600 |
| \([.40,.50]\) | \(1/20\) | 10 | 464 | 60,593 |
| \([.50,.60]\) | \(1/20\) | 8 | 176 | 23,312 |
| \([.60,.75]\) | \(1/20\) | 3 | 166 | 22,247 |
| \([.75,1]\) | \(1/8\) | 1 | 64 | 8,748 |

The initial lower endpoint is the bound at \(B=23/4\). In each
subdivision the [exact \(u^2\) product
enumeration](MULTIVARIATE_TERMINAL_THIRTEEN_RETAINED_BANDS.md#retaining-rows-on-a-band)
lists all additional retained row multisets at the specified band
and caps. It includes the empty choice and keeps both orientations
because the private row is fixed. For a nonempty choice, divide the
mean interval by the exact mean contribution of the additional
rows and intersect with \((0,1]\). The narrowed cap applies to that
new residual array.

## Paired ray checks

The checker uses the terminal-thirteen sector verification at orders
thirteen and fourteen. The upper roots of the Gaussian midpoint
anchors have verified order pattern \(14,13\) repeated six times,
followed by \(14\). The odd order has one additional root on the
negative ray. The ray and annulus exclusions preserve the full-sector
root counts throughout the interpolation. In particular the final
sector joined to its reflection contains one order-thirteen root,
which is real by conjugation. The even order's two roots in that
sector are nonreal by strict real-axis positivity of its actual
moment. The other sectors separate the arguments of the remaining
roots. The original moment polynomials have integer coefficients,
so the coprimality transfers by rational Bézout.

The native projector is built in a temporary directory from the
existing source. The checker checks each of the three source
substitutions once and sets the maximum to 15 coefficient slots.
For this size it uses an outward dyadic shift of 17 on the Bernstein
transform, leaving entries bounded by \(2^{82}\) in the
48-bit grid. The input guards bound the projected sine and
derivative products by \(2^{62}+2^{66}\) on the same grid.
Thus each outward Bernstein control is bounded by
\(15(2^{62}+2^{66})2^{34}+1\). The binomial weights for a midpoint
sum to \(2^{14}\). Each scalar group parameter is bounded by \(2\),
and there are at most 256 groups. Including the product rounding,
the checker verifies that their signed accumulation is less
than \(2^{127}\).

Run the [checker](verify_multivariate_terminal_fourteen_retained_profiles.py)
with Python 3, python-flint, and a C++17 compiler from the repository
root.

~~~sh
python3 research/factorial/image-rigidity/verify_multivariate_terminal_fourteen_retained_profiles.py
~~~

The three replay families are \(\texttt{shared}\),
\(\texttt{one-private}\), and \(\texttt{other-private}\). The
fully shared family above cap \(1/8\) checks 278 lists, 682 outer
cells, and 67,984 radial cells. The 88 nonempty one-private lists
check 156 outer and 13,799 radial cells. The empty one-private
list uses the subdivisions in the table. The eight degree-two
and degree-three private allocations below their fixed-ray
thresholds check a further 33 lists, 97 outer cells, and
10,620 radial cells.
