# Terminal-ten angular separation with a marked shared unit row

Let two distinct ordinary exponent columns have only shared unequal
rows. Let \(u\) be their order-nine outer geometric mean after
positive endpoint normalization and removal of the binomial factors.
Their ninth and tenth factorial moment polynomials have disjoint
root arguments under independent positive rescalings whenever one
of the following rows is present in either orientation.

| Marked row | Lower bound on \(u\) |
|---|---:|
| \((2,3)\) | \(19/25\) |
| \((3,4)\) | \(21/25\) |
| \((4,5)\) | \(22/25\) |
| \((5,6)\) | \(22/25\) |

The dimension and all remaining shared exponents are arbitrary.
Rational Bézout gives coprimality in characteristic zero.

## Removing the marked row

For a marked row \((a,a+1)\), its contribution to the square of
the order-nine outer mean is \(u_a^2=(9a+1)/(9a+9)\).
Two copies would give \(u\le u_a^2\), below the displayed lower
bound, so the marked row is unique. Write \(v=u/u_a\) for the
residual order-nine outer mean. The checker verifies
\(u_{min}^2>u_a^2v_{min}^2\) with the following thresholds.

| \(a\) | \(v_{min}\) | Order-nine reciprocal cap \(c_0\) |
|---:|---:|---:|
| 2 | \(113/125\) | \(1/37\) |
| 3 | \(119/125\) | \(1/82\) |
| 4 | \(97/100\) | \(1/136\) |
| 5 | \(119/125\) | \(1/82\) |

For a residual row with lower exponent \(\ell\) and gap \(d\),
put \(x_h=d/(9\ell+h)\) and \(y=d/(9\ell)\).
Its contribution to the square of \(v\) gives
\[
 v^2\le\prod_{h=1}^d\frac{9\ell+h}{9\ell+8d+h}
      \le\left(\frac{1+y}{1+9y}\right)^d .
\]
When \(d=1\), the exact ratio \((9\ell+1)/(9\ell+9)\)
and \(v\ge v_{min}\) rule out every \(\ell\) below the
first unit row with \(1/(9\ell+1)\le c_0\). For \(d\ge2\),
the product bound gives \(x_h\le y\le(1-v)/(9v-1)<c_0\).
The checker verifies all these comparisons rationally.
Also \(v^{-2}\) is the product of \(1+8x_h\) over the
residual factors. Thus the sharper cap on a residual
subinterval with left endpoint \(v_-\) is
\(c=min(c_0,(v_-^{-2}-1)/8)\). The matched upper-order
reciprocal cap is \(9c/(10-c)\).

Swap both columns if needed to make the marked row increase.
The normalized marked core is
\[
       C_{j,k}=\binom jk
          \frac{(ja+k)!}{((ja)!)^{1-k/j}((j(a+1))!)^{k/j}},
          \qquad j=9,10 .
\]
Keep the core fixed and apply the grouped logarithmic coefficient
model to the residual factors. The midpoint references multiply
\(C_{9,k}\) by \(v_0^{e_{9,k}}\) and \(C_{10,k}\) by
\(v_0^{(81/80)e_{10,k}}\). The ratio of adjacent residual
outer logarithms is between
\[
  \frac{\log(1+(81/10)c)}{\log(1+8c)}
     \quad\text{and}\quad
  1+\frac{\log(1/(1-c/10))}{\log(1+8c)} .
\]
Indeed the one-factor matched quotient \(x'/x=9/(10-(h/d)x)\)
lies between \(9/10\) and \(9/(10-x)\). The bounding
logarithm ratios are respectively decreasing and increasing
in \(x\), by comparison of the ratios of their derivatives.
A weighted average with weights \(\log(1+8x)\) gives the
bounds for all residual factors. The checker keeps the
linked outer, three skew and two mean-defect directions
and their quadratic products. It bounds the log and
exponential tails per coefficient, and brackets each marked
factorial power and reference power by integer powers.
[Grouped coefficient enclosure](MULTIVARIATE_TERMINAL_NINE_PRIVATE_LINEAR.md#grouping-the-shared-logarithmic-corrections).

## Wide Cayley contours and the marked-core anchor

Put \(p=1/8\), \(\lambda=\sqrt p\), and
\(t=(\lambda X-1)/(\lambda X+1)\). Apply the Cayley matrix to
each complete coefficient direction. Use the circle
\(|X|=11/4\), which maps into the physical left half-plane
since \(p(11/4)^2=121/128<1\).
For each residual subinterval, the checker tests either strict
leading-term dominance on that circle or a strict projected
sign on each of its upper semicircle boxes. Conjugation covers
the lower semicircle. The boundary exclusion does not require
leading-term dominance in the subintervals where the triangle
bound is too wide.

The eight separating ray angles are \(\pi-2\arctan(\lambda b)\).
Their levels \(b\) are ordered and piecewise affine in \(v\).
For \(X=x+iy\) on such a ray the equation is
\((y-b)(1+pby)=-pbx^2\). Enclose \(y\) by three interval updates
from \([0,b]\) and test both signs of \(x\). For these rays
and for the circle, retain the grouped coefficient contributions
when evaluating the complex signs through centered Taylor
enclosures. All accepted projections have strictly signed
rational interval bounds.
[Centered Taylor projection bounds](MULTIVARIATE_TERMINAL_TEN_NEAR_FLAT.md#exact-rays-and-sector-counts).

The residual fractions can be deformed to zero while preserving
these verified bounds. Add the same real nonnegative shift
to both exponents in every residual row, leaving the marked
row fixed. The gaps remain integers, so the finite Gamma
product identities still hold. The matched reciprocals decrease,
preserving their initial \(c_0\) cap, while the residual mean
increases to one. The product bound for \(v^{-2}\) still holds
at every step. The normalized residual coefficients tend to one,
giving a continuous path to the marked core.

At that core, the checker verifies disjoint rational Rouché
disks containing all roots. Each disk is inside the circle and
in its stated sector. The upper sectors in order of increasing
argument have counts \(10,9,10,9,10,9,10,9,10\), and there is
one additional real nonic root. The circle and ray exclusions
preserve these counts along the deformation. In the last sector
joined to its conjugate, the nonic root is real by conjugacy.
The two decic roots there are nonreal because the actual tenth
factorial moment is strictly positive for every real linear
combination with nonzero constant coefficient of the two
distinct monomials. Thus the root arguments are disjoint.

The exact replay certifies 71,155 ray cells across 896
order, ray, and subinterval checks, in addition to verifying
the circle boundaries and all eight core anchors. The checker
validates every expected ray count.
[Exact replay checker](verify_multivariate_terminal_ten_marked_unit_rows.py).

~~~sh
python3 research/factorial/image-rigidity/verify_multivariate_terminal_ten_marked_unit_rows.py
python3 research/factorial/image-rigidity/verify_multivariate_terminal_ten_marked_unit_rows.py 2 10
~~~
