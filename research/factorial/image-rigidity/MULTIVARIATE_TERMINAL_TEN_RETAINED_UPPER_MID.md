# Terminal-ten separation between outer means .76 and .88

Let two distinct ordinary exponent columns have only shared unequal
rows. On each row let \(\ell\) and \(d>0\) be the smaller exponent
and gap. Write \(u\) for the normalized order-nine reflected outer
mean and \(S=\sum d^2/(\ell+d)\).

**Theorem.** If \(19/25\le u\le22/25\), the ninth and tenth
factorial moment polynomials have disjoint root arguments under
independent positive rescalings. This holds in arbitrary finite
dimension and transfers to coprimality in characteristic zero.

## The seventeen retained lists

The [larger reciprocal certificate](MULTIVARIATE_TERMINAL_TEN_RETAINED_PROFILES.md)
covers any \(d/(9\ell+1)>1/8\). For the remaining profiles retain all
rows with \(d/(9\ell+1)>1/45\), in their actual orientations. Since
\(-\log u\ge4S/9\) by the outer-loss bound in that proof,
\[
                       S\le(9/4)\log(25/19).
\]
Each retained row satisfies
\(\ell<d/(9/45)-1/9\). This and the curvature budget bound the gap
and base in the exact enumeration. The checker excludes any row
whose squared outer mean is less than \((19/25)^2\). It applies the
same test to the product of squared means for each list.
The [finite replay](verify_multivariate_terminal_ten_retained_upper_mid.py)
finds seventeen lists up to simultaneous reversal, including the empty
list, by integer and rational comparisons.

For a retained list with squared order-nine outer mean \(q\), the
residual mean \(v\) satisfies \(u=\sqrt q\,v\). The checker brackets
\(\sqrt q\) outward and intersects the induced interval for \(v\)
with the residual curvature bound \(v\ge 2^{-(17/18)B}\), where \(B\)
is the remaining curvature budget. It uses the residual reciprocal
cap \(\min(1/45,(v_-^{-2}-1)/8)\) on intervals starting at \(v_-\),
with the second bound rounded outward on a rational grid.

## Paired rays with an exact annulus

The [retained-profile comparison](MULTIVARIATE_TERMINAL_TEN_RETAINED_PROFILES.md#tail-strength-and-the-two-linked-polynomials)
gives independent coefficient intervals with both endpoint attenuation
constraints, and the same log-linear homotopy to Gaussian residual
multipliers. The reference midpoint here is the center of each residual
outer interval. The [terminal-ten cubic light model](MULTIVARIATE_TERMINAL_TEN_SHARED_LIGHT.md#grouped-moment-and-moving-ray-comparison)
also encloses this homotopy after multiplying each coefficient group
by the exact factorial core. The cubic model allows the full signed
tail interval because the retained core can be asymmetric. Since the
rays for a residual cell are fixed, the shared ray displacement has
zero slope in the fused projection.

For each residual cell the checker verifies the annulus
\(1/12<|z|<12\) by comparing consecutive coefficient intervals.
Indeed if \(0<a_k\) and \(a_k/a_{k+1}<R\), the coefficients of
\(P(Rz)\) are increasing. Multiplication by \(1-z\) makes its leading
term strictly dominant for \(|z|>1\), since all other coefficients
are nonnegative and sum to the absolute leading coefficient.
Apply the same argument to the reciprocal. The strict interval
checks give radii strictly inside the displayed annulus for the
actual array and every point on the coefficient homotopy.

At each residual midpoint the comparison roots are isolated in
outward complex balls. Rounding proposes rational cosines only.
The checker verifies that the rays lie strictly between successive
anchor arguments, and that the two orders alternate. Exact projection
of either the independent boxes or the cubic grouped model excludes
zeros on these rays throughout the annulus. Subdivision in the
residual outer mean or signed parameter is used when a strict
projection does not yet follow. Both types of enclosure contain
the same log-linear path and the Gaussian ratio adjustment of the
retained-profile proof.

The sector joined across the negative real axis contains the single
real nonic root and the two last decic roots. Conjugation preserves
the reality of the nonic root. The decic roots are nonreal since
the actual even moment is positive at every real argument.
The other sectors separate the arguments directly. Rational Bézout
gives coprimality in characteristic zero.

The checked native replay uses 264 outer cells and 38,291 radial
cells. The default Python path uses the same outward arithmetic
without the optional compiled projection.

~~~sh
python3 research/factorial/image-rigidity/verify_multivariate_terminal_ten_retained_upper_mid.py --native
~~~
