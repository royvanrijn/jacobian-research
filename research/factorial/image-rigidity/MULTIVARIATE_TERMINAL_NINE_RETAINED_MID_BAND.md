# Terminal-nine separation on the fully shared .55–.62 band

Let two distinct ordinary exponent columns have all unequal rows shared.
For an unequal row let \(\ell\) be the lower exponent and \(d\) the gap.
Write \(u\) for the normalized order-eight outer reflected mean. If
\[
       .55\le u\le .62,\qquad
       \max_{\rm rows}\frac{d}{8\ell+1}\le\frac18,
\]
the eighth and ninth moment polynomials have disjoint root arguments
under independent positive rescalings. In particular they are coprime
over every characteristic-zero field.

## Exact retained lists

The case where every \(d/(8\ell+1)\) is at most \(1/10\) is the
[light-row band](MULTIVARIATE_TERMINAL_NINE_SHARED_MATCHED_UPPER.md).
Otherwise retain every row whose reciprocal exceeds \(1/10\).
The [fully shared high-curvature theorem](MULTIVARIATE_TERMINAL_NINE_SHARED_HIGH_CURVATURE.md)
covers \(S\ge53/10\).
For the remaining profiles the
[retained enumeration](MULTIVARIATE_TERMINAL_NINE_RETAINED_BANDS.md#finite-enumeration-and-tighter-tail-cap)
uses the additional budget \((16/7)\log(1/u_-)\), with \(u_-=11/20\).
It prunes rows and lists by their exact squared order-eight
outer means, and identifies each list with its column reversal.

The [checker](verify_multivariate_terminal_nine_retained_mid_band.py)
confirms the following four lists. Here an exponent pair records a
retained row with both its base and its orientation.

| List up to column reversal | Exact verification |
| --- | --- |
| One \((2,4)\) row | Retained ray replay |
| Two \((1,2)\) rows, same orientation | Retained ray replay |
| Two \((1,2)\) rows, opposite orientations | Retained ray replay |
| One \((1,2)\) row | Matched moving rays with fixed core |

The first three lists use the existing
[retained annulus and paired rays](MULTIVARIATE_TERMINAL_NINE_RETAINED_BANDS.md#annuli-and-paired-rays)
with the residual reciprocal cap \(1/10\). Their replay tightens this
cap on each residual outer interval by the separate unit and larger
gap bounds. The final list uses the following fixed-core calculation.

## The marked unit row

Exchange columns if needed so the retained pair is \((1,2)\). The
exchange inverts roots up to a positive scaling, and real
coefficients preserve their upper arguments after conjugation.
The exact normalized marked coefficients, including binomial
factors, are
\[
 C_{j,k}=\binom jk
          \frac{(j+k)!}{(j!)^{1-k/j}((2j)!)^{k/j}},\qquad j=8,9.
\]
The squared order-eight reflected mean of this row is \(9/16\).
The residual mean \(v=u/(3/4)\) therefore lies in
\([11/15,62/75]\). Every residual reciprocal is at most \(1/10\).
The [ordinary row pruning calculation](MULTIVARIATE_TERMINAL_NINE_SHARED_MATCHED_UPPER.md#ordinary-reciprocal-bounds)
gives the smaller cap \(1/17\) with this lower residual endpoint.

Use the
[common coefficient and ray displacement model](MULTIVARIATE_TERMINAL_NINE_SHARED_MATCHED_UPPER.md#common-coefficient-and-ray-displacement)
for the residual profile, multiplying each grouped coefficient
of degree \(k\) by the positive exact factor \(C_{j,k}/\binom jk\).
This multiplies both parts of the grouped displacement at degree
\(k\) by the same factor. The marked row fixes the orientation,
so the signed residual parameter ranges over all of \([-1,1]\).
The checker uses outward rational enclosures of the marked factors.

Independent coefficient caps verify the annulus \(1/12<|z|<12\).
At fixed \(v\) the same logarithmic homotopy as for light rows
scales the residual defects and skews to zero, while preserving
the marked factors. Then move the outer-mean ratio to \(64/63\)
and slide the residual Gaussian references to \(v=11/15\).
The ray cosines are ordered, continuous and piecewise affine in \(v\).
The verified strict Bernstein projection signs permit the
homotopy along every ray cell.

Exact rational Rouché disks at that endpoint give alternating
upper sector orders \(9,8,9,8,9,8,9,8\), with one extra negative
real nonic root. The latter remains real in the final sector
joined to its conjugate. Strict positivity of the actual eighth
moment on the real axis follows from the Gamma integral. This
excludes any equality of root arguments in this sector as well.

## Reproduction

The complete checker replays the light-row band, the exact list
enumeration, all three retained ray computations and the marked
unit row. Run from the repository root with

~~~sh
python3 research/factorial/image-rigidity/verify_multivariate_terminal_nine_retained_mid_band.py --native
~~~

The optional --native flag uses the [same exact C++ Bernstein
implementation](verify_multivariate_terminal_nine_shared_matched_upper.cpp)
as the light-row checker. The replay verifies 18,768 light-row
radial cells. The three other retained lists use 12 outer cells
and 951 radial cells in total. The marked unit row uses 6,461
radial cells.
