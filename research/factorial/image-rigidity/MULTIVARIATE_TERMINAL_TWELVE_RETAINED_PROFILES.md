# Terminal-twelve separation with retained shared rows

Let \(a,b\) be distinct nonnegative integer exponent columns.
Delete equal rows, and write \(\ell,d\) for the smaller exponent and
positive gap in a shared unequal row. Its strength is
\(s=d^2/(\ell+d)\). Let \(S_0\) be the sum of these strengths.

**Theorem.** The eleventh and twelfth factorial moment polynomials
have disjoint root arguments under independent positive rescalings
if there is exactly one private degree-one row, or if all rows are
shared and some reciprocal \(d/(11\ell+1)\) exceeds \(1/8\).
The two polynomials are therefore coprime in characteristic zero.

The [fixed-ray theorem](MULTIVARIATE_TERMINAL_TWELVE_HIGH_CURVATURE.md)
covers the first family when \(S_0\ge6\) and the second when
\(S_0\ge8\). Below those cutoffs retain each shared row with
\(d/(11\ell+1)>1/8\) and its orientation. The retained list is
nonempty in the second family. In the first family the single
private row fixes the orientation of the whole profile, so both
orientations of each shared row are enumerated, including the
possibility of an empty retained list.

## Retained lists and the residual interval

For a row above a positive cap \(c\), the smaller exponent satisfies
\(\ell<(d/c-1)/11\). In a strength budget \(B\) this gives
\(d<B(1+1/(11c))\). These inequalities bound the integer candidate
rows before enumeration. Their positive strengths bound their
multiplicities. In the fully shared case a whole-column reversal
does not change angular disjointness. The checker puts the candidates
in lexicographic order of \((\ell,d,\mathrm{orientation})\), then
compares the ordered tuple of selected rows with the sorted tuple
of reversed rows. At least one of each pair is included.

Write \(B\) for the strength budget left after retaining the rows.
Let \(F_{j,k}\) be the normalized residual factorial multiplier,
without the binomial coefficient, and let
\(v=(F_{11,1}F_{11,10})^{1/2}\). For every residual row the Gauss
product gives a contribution to \(-\log v\) bounded above by
\((d/2)\log(1+10x)\), where \(x=d/(11\ell+1)\). Since
\(d/s<1+1/(11x)\) and \((1+1/(11x))\log(1+10x)\) increases with
\(x\), one obtains
\[
        -\log v \le \frac{19}{22}\log(9/4) B.
\]
To verify the monotonicity, substitute \(z=10x\). The derivative
has the sign of \(z(1+11z/10)-(1+z)\log(1+z)\), and the latter
increases from zero by \(11z/5-\log(1+z)>0\).
The replay thus encloses \(v\) between an outward rounding of
\((4/9)^{19B/22}\) and \(1\). On a subinterval starting at \(v_-\),
every reciprocal is at most \((v_-^{-2}-1)/10\), as
\(v^{-2}\) is the product of the positive factors \(1+10x_h\).
The [ordinary-cap enumeration](verify_multivariate_terminal_twelve_shared_light.py)
improves this bound by discarding rows that cannot occur at that
mean.

## Coupled coefficient bounds

The [all-order endpoint attenuation](MULTIVARIATE_ORDINARY_ADJACENT_COEFFICIENT_ATTENUATION.md#an-all-order-endpoint-inequality)
applies to the residual shared factors with \(M\ge2\). Its multiplier
at orders eleven and twelve is at most \((275/276)^{k^2 S_{\rm res}}\)
at distance \(k\) from either endpoint, where \(S_{\rm res}\) is the
residual strength. By \(\log(276/275)\ge1/276\) and the exact
comparison \(276(19/22)\log(9/4)<194\), this gives
\[
 F_{12,k} \le F_{11,k}\,v^{k^2/194},\qquad
 F_{12,k} \le F_{11,k-1}\,v^{(12-k)^2/194}
 \qquad (1\le k\le11).
\]

The endpoint-normalized Gaussian comparison factors are
\(v^{e_{11,k}}\) and \(v^{\rho e_{12,k}}\), respectively, with
\(e_{j,k}=k(j-k)/(j-1)\). For a lower reciprocal \(x\) and its
matched upper reciprocal \(x'=11x/(12-\gamma x)\), \(0<\gamma\le1\),
\[
 1+\frac{121}{12}x \le 1+11x'
    \le \frac{1+10x}{1-x/12}.
\]
The lower logarithmic quotient is nonincreasing in \(x\), by
concavity of \(y\mapsto\log(1+(121/120)(e^y-1))\) on \(y\ge0\).
The upper quotient is nondecreasing since
\(-\log(1-x/12)/x\) increases and \(\log(1+10x)/x\) decreases.
Thus at cap \(1/8\) the lower allowed ratio is
\(\log(217/96)/\log(9/4)>1+1/194\), verified by rational logarithm
enclosures. The two endpoint inequalities also hold on the
Gaussian families because
\[
 (1+1/194)e_{12,k}-e_{11,k}-k^2/194
  = k(k-1)\left(\frac1{110}-\frac{12}{11\cdot194}\right)\ge0
\]
and reflection treats the other endpoint. The reference ratio
\(121/120\) belongs to the allowed interval.

The independent coefficient bounds are intersected with the
endpoint inequalities. When necessary the checker switches to the
[grouped coefficient bounds](MULTIVARIATE_TERMINAL_TEN_SHARED_LIGHT.md#grouped-moment-and-moving-ray-comparison).
It uses the order-eleven and order-twelve scales, with the exact
matched skew constant from the
[tanh comparison](MULTIVARIATE_MATCHED_TANH_MOMENTS.md).
Each coordinate is multiplied by an outward enclosure of the
corresponding retained coefficient multiplier.

At actual \(v,\rho\) scale the logarithmic skew and mean-defect
moments jointly to zero. The Jensen lower bounds hold throughout
because \(t^{2l}\le t\) for \(t\in[0,1]\) and \(l\ge1\).
The logarithmic coefficient inequalities also persist, since
they are affine in the interpolated logarithms and hold at both
ends. Afterwards move \(\rho\) to \(121/120\), then move \(v\)
to its subinterval midpoint. These Gaussian paths obey the same
bounds.

## Exact sector separation

The [checker](verify_multivariate_terminal_twelve_retained_profiles.py)
bounds all roots in the annulus between radii \(1/12\) and \(12\).
It verifies constant and leading coefficient domination or
strict Eneström–Kakeya inequalities for its independent boxes.
Those boxes contain the complete interpolation, whichever
projection proof a ray uses. The two Gaussian midpoint arrays
with their retained core have isolating complex balls whose upper
root cosines are strictly ordered. For each gap the checker
selects a rational cosine and verifies a strict outward Bernstein
projection throughout the annulus. Floating-point estimates only
propose rays. All inclusions and separation tests use exact
rational bounds or outward complex ball bounds.

The midpoint upper root orders are \(12,11,12,11,12,11,12,11,12,11,12\)
in every accepted cell. The nonreal roots cannot cross a separating
ray during interpolation. The last sector joined to its reflected
sector contains one odd-order root, so conjugacy keeps that root
negative real. It contains two even-order roots, which are
nonreal for the actual moment because its Gamma integral is
strictly positive on the real axis. This separates the arguments
in all sectors. Rational Bézout gives characteristic-zero
coprimality.

Run from the repository root.

~~~sh
python3 research/factorial/image-rigidity/verify_multivariate_terminal_twelve_retained_profiles.py --native private
python3 research/factorial/image-rigidity/verify_multivariate_terminal_twelve_retained_profiles.py --native shared
~~~

The checker verifies 244 one-private lists with 401 outer cells
and 26,899 radial cells. It verifies 603 fully shared lists with
911 outer cells and 56,666 radial cells. Either command also accepts specific list indices. The
native path uses the same outward Bernstein routines as the
terminal-eleven proof, with one additional coefficient slot and
an outward dyadic rescaling of the degree-twelve transform.
