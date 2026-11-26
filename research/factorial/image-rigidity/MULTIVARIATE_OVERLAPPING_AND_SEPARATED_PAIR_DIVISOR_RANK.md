# Full divisor rank for overlapping and scale-separated factorial pairs

Fix a finite ordinary two-monomial core and \(d\ge2\).
Consider the four appended types
\[
      (U,U+d),\qquad (U+1,U),\qquad
      (V,V+1),\qquad (V+d,V).                                \tag{1}
\]
There are two conclusions about exact coprimality.

- For \(U=M,\ V=M+K\), with any fixed integer \(K\ge1\), give
  all four types a common positive multiplicity \(r\).
  At any positive limiting curvature, the exact original moments
  are coprime in every degree for sufficiently large \(M\),
  uniformly over the actual integer \(r\) in any fixed
  positive compact range of \(r/M\). This includes overlapping
  ramps \(1\le K<d\) and the balanced obstruction \(K=d\).
- For \(U=M^a,\ V=M^b+K\), with integers \(1\le a<b\),
  \(v_2(a)\ne v_2(b)\), and \(K\ge0\), assign the four counts
  \[
                       r,\quad dr,\quad dr M^{b-a},\quad rM^{b-a}.
  \]
  The same exclusion holds when \(r/M^a\) ranges over a fixed
  positive compact interval. Both scale families contribute
  positive curvature.

In either case, for an empty core or for a fixed core whose own
ordinary moments are already coprime in consecutive orders, there
is one computable height threshold, uniform over **all** positive
integer \(r\), above which the original binomial satisfies
physical SFC(2) in every degree. When the core itself has a
collision, the conclusion here excludes positive limiting
strengths, and vanishing strength still requires a weak-tail
comparison. Rational coprimality extends to every
characteristic-zero field.

The proof gives full divisor row rank at every adjacent pair.
It uses the
[affine divisor-rank theorem](MULTIVARIATE_AFFINE_DIVISOR_RANK_DIFFUSE_EXCLUSION.md)
for the first case and the
[polynomial divisor-rank theorem](MULTIVARIATE_POLYNOMIAL_SCALE_DIFFUSE_EXCLUSION.md)
for the second. The affine source includes the disjoint-window
example \(K=d+1\). The rank proof below uses only selected
divisors at the two ends of the two windows, including when their
interiors overlap. In the polynomial case no single exponent
scale carries all four types of a full-rank template.

## Raw counts and the normalization map

At order \(j\) strip the smaller-end factorial of each appended
row. The row contributes its initial rising factors of length
\(dk\) or \(d(j-k)\) for gap \(d\), with the corresponding
length \(k\) or \(j-k\) for gap one. Scalar-normalize each
order at \(k=0\). For raw divisor columns this amounts to
subtracting the contribution at \(k=0\) from each coordinate.
Write \(S_\ell\) for the vector on \(1\le k\le j\) with
\(S_\ell(k)=1\) if \(k\ge\ell\) and \(S_\ell(k)=0\) otherwise.
The \(S_\ell\), \(1\le\ell\le j\), are a basis since
\(S_\ell-S_{\ell+1}\) is the unit vector at \(\ell\), taking
\(S_{j+1}=0\).

At adjacent \(n,N=n+1\), concatenate the raw coordinates of
the two orders. For a raw divisor column \(v\), the normalized
column of the divisor theorem has coordinates
\[
             (\Pi v)_{j,k}=n v_{j,k}-k v_{n,n},
                 \qquad (j,k)\ne(n,n).
\]
The omitted index is the lower leading coefficient, since
this normalization makes its multiplier one. The map is
surjective onto the \(n+N-1\) retained coordinates: prescribe
\(v_{n,n}=0\) and divide each desired retained coordinate by
\(n\). Therefore separate spanning of the two raw blocks by
order-exclusive divisors suffices for full rank. The common
leading-factor powers removed in either divisor theorem cancel
under this map, and nonzero rational constants do not contribute
divisor exponents.

## Spanning with overlapping affine windows

Consider any order \(j\ge2\) in the first template.
The proportional class of the factor \(j(M+q)\), with
\(q=t/j\) and \(1\le t\le j-1\), lies in the open low
interval \(0<q<1\). Thus it belongs only to the first two
rows of (1), since \(K\ge1\). The long increasing row first
contains its factor at \(\lceil t/d\rceil\).
The short decreasing row contains the factor precisely up
to \(k=j-t\), including \(k=0\). The raw column is consequently
\[
                   S_{\lceil t/d\rceil}-S_{j-t+1}.           \tag{2}
\]
Put \(L_j=\lceil(j-1)/d\rceil\). For each
\(1\le\ell\le L_j\), use instead the factor of normalized
offset
\[
                         q=K+d-\frac{(\ell-1)d+1}{j}.        \tag{3}
\]
The integer \((\ell-1)d+1\) is between \(1\) and \(j-1\).
Consequently \(q>d\) and \(q-K>d-1\ge1\). These inequalities
place the factor outside both low ramps and above the upper
short ramp. The long decreasing row in the upper window uses
this factor at \(k=0,\ldots,\ell-1\): its local threshold is
\[
                jd-((\ell-1)d+1)=(j-\ell+1)d-1.
\]
After subtracting the index-zero count, the column is exactly
\(-S_\ell\). Thus all \(S_\ell\) with \(\ell\le L_j\) are in
the column span. Adding them back in (2) produces
\(S_{j-t+1}\) for \(1\le t\le j-1\), hence the entire
raw block. Every selected \(q\) is nonintegral.
A nonintegral rational cannot have both \(nq\) and \(Nq\)
integral, since the orders are coprime. These columns are
therefore order-exclusive, and the map \(\Pi\) gives full rank
in every adjacent pair with \(n\ge2\).

The analytic hypotheses on counts and gaps for the affine
divisor-rank theorem hold with common slope one, no additive
counts, zero signed gap, and leading curvature parameter
\(2(d^2+1)r/M\). The stated exclusion on positive compact
strength ranges now follows from that theorem and the absolute
terminal cutoff.

## Spanning at two separated polynomial scales

In the second template write \(e=\gcd(a,b)\). Lower divisors
have the form \(M^a+q\) and upper divisors the form \(M^b+K+q'\),
with \(q,q'>0\) in the stripped ramps. They have no common
complex root. Otherwise raising the two equalities for
the powers of that root to \(b/e\) and \(a/e\), respectively,
would give the same power with opposite real signs: precisely
one of \(b/e,a/e\) is even. Within either scale, different
normalized rational offsets have coprime divisors since their
polynomials differ by a nonzero constant. Each such polynomial
has only simple roots in characteristic zero, so any one of
its irreducible factors supplies the same count column.

Use the lower offsets \(q=t/j\), \(1\le t\le j-1\), as
before. With the new multiplicity weights their columns are
\[
                     S_{\lceil t/d\rceil}-dS_{j-t+1}.
\]
At the upper scale choose the local offset
\(q'=d-((\ell-1)d+1)/j>1\) for
\(1\le\ell\le L_j\). Its column is
\(-M^{b-a}S_\ell\), from the long decreasing upper row alone.
There is no overlap with a lower divisor by the preceding
disjointness argument. Both sets of offsets are nonintegral,
so the columns are again order-exclusive.
Over \(\mathbb Q(M)\), their span contains all the
\(S_\ell\) by the same subtraction argument. This proves the
full rank required for polynomial templates, after applying
\(\Pi\).

Within each of the two scales the signed gap is zero with the
indicated count ratio, so the total signed gap is the zero
polynomial. The smallest exponent-degree difference in the
polynomial theorem is \(a>0\) for all four types:
the lower rows have exponent degree \(a\) and constant
main count polynomials, and the upper rows have degree \(b\)
and main count polynomials of degree \(b-a\).
Thus both scales are active. The leading curvature parameter
is \(2(d^2+d)r/M^a\), half from each scale. The hypotheses
of the polynomial divisor-rank theorem are satisfied.

## Passing to the full range of multiplicities

For each case, the positive-strength conclusions hold at
every possible Gaussian contact, independently of simplicity
of the limiting roots, by the cited divisor-rank theorems.
The relevant transcendence input is that the exponential
of a nonzero algebraic number is transcendental.
[Hermite–Lindemann](https://isa-afp.org/entries/Hermite_Lindemann.html)

The [absolute curvature bound](MULTIVARIATE_ABSOLUTE_CURVATURE_12.md)
bounds the limiting tail-strength parameter at any collision,
uniformly over the dimension. Concretely each of those parameters
is the actual tail curvature times \(1+O(M^{-1})\) in the affine
case, and times \(1+O(M^{-a})\) in the polynomial case.
Hence outside a bounded range of \(r/M\) or \(r/M^a\), the
absolute curvature theorem already certifies separation.
For strengths bounded above and away from zero the cited
divisor-rank results give computable neighborhoods of the
finitely many positive reference contacts; coefficient
continuity on the compact remainder gives a uniform cutoff.
The [absolute ordinary cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL.md)
leaves only finitely many paired orders to treat.

For a coprime fixed core, the signed twist is zero in both
cases, so the strength-zero moment pair in the same physical
normalization is exactly the fixed core pair. Ordinary
coefficient continuity handles strengths close to zero.
For the empty core, the
[uniform small-curvature theorem](MULTIVARIATE_UNIFORM_WEAK_CURVATURE_ANGLE_SEPARATION.md)
applies directly to the nonempty appended profile.
This proves the claimed all-degree threshold over all \(r>0\)
in these two core cases.
