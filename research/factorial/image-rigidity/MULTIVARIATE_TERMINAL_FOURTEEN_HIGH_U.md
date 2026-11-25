# Terminal-fourteen angular separation near the flat limit

Let two distinct ordinary exponent columns have only shared unequal
rows. Write \(u\) for their order-thirteen reflected outer mean after
endpoint normalization and removal of the binomial coefficients.
The order-thirteen and order-fourteen factorial moment polynomials
have disjoint root arguments under independent positive rescalings
in the following ranges.

| \(u\) | Cap on \(d/(13\ell+h)\) |
|---|---:|
| \([19/20,1)\) | Unrestricted |
| \([23/25,19/20]\) | \(1/100\) |
| \([22/25,23/25]\) | \(1/80\) |
| \([171/200,22/25]\) | \(1/70\) |
| \([19/25,171/200]\) | \(1/60\) |

The original rational polynomials are coprime over every
characteristic-zero field by rational Bézout.

## Deficits and row constraints

Put \(p=1-u,\ q=1-w,\ r=q/p\), with \(w\) the order-fourteen outer
mean. Write \(U=-\log u,\ W=-\log w\) and \(R=169/168\).
The [matched outer-ratio inequality](MULTIVARIATE_MATCHED_OUTER_RATIO_MONOTONICITY.md)
and the individual upper comparison give
\[
 R\le W/U\le L=1-\frac{\log(1-c/14)}{\log(1+12c)},
\]
where \(c\) bounds the lower reciprocals. The upper comparison follows
from \(x'=13x/(14-(h/d)x)\) and
\(1+13x'\le(1+12x)/(1-x/14)\). The quotient of the two logarithms
increases with \(x\), by convexity of the numerator and concavity of
the denominator after taking both terms positive. Its value at
\(c=0\) is interpreted by the limit \(L=R\).

The product for \(u^{-2}\) gives \(c\le(u^{-2}-1)/12\).
The tighter lattice estimate keeps the largest \(1/(13\ell+1)\)
allowed by the cap for \(d=1\). For \(d\ge2\), the other reciprocals
in the row are at least \(x_1/(1+x_1)\) if \(x_1\) is the largest.
Using two factors in the product gives \(x_1\le(1-u)/(13u-1)\).
The maximum of the two case bounds is valid for the whole list.
For \(a\ge1\), the secant \(\big(1-(1-p)^a\big)/p\) decreases
with \(p\) by concavity of \(1-(1-p)^a\). On a deficit slab
\([p_-,p_+]\), use the outward bounds
\[
         \frac{1-(1-p_+)^R}{p_+}\le r
              \le \frac{1-(1-p_-)^L}{p_-}.
\]
The upper bound at \(p_-=0\) is \(L\). At a zero-width slab
at \(p=0\) the limiting ratio is \(R\). The checker takes dyadic
enclosures of these powers using the convergent exponential series.

At order \(j\) use the scale \(\lambda=\sqrt h\), where \(h=p\)
or \(q\). In the [grouped Cayley model](MULTIVARIATE_TERMINAL_TWELVE_CAYLEY.md#weighted-cayley-model)
the normalized tanh factors \(a_i=T_i/h\) satisfy
\(a_i\le A\) and \(\sum_i a_i\le C\), with
\(C=(1-q_0/2)/(1-q_0)\) for an outward deficit bound \(q_0\).
Bound \(A\) by \(C\) and, when the lower deficit endpoint is
positive, by the tanh cap divided by that endpoint.
For the upper tanh cap use \(x'\le13c/(14-c)\).

## Rays for \(u\ge .98\)

Use the [factorial-to-Gaussian annulus path](MULTIVARIATE_FACTORIAL_GAUSSIAN_ANNULUS_HOMOTOPY.md)
at fixed \(p\). Any initial lattice-refined bound survives, since each
surviving reciprocal decreases and the replacement Gaussian introduces
no tanh factors. The path ends at the Gaussian ratio
\(r_G(p)=(1-(1-p)^R)/p\). Continue this Gaussian curve to the limiting
scaled Cayley polynomials at \(p=0\).

There are twelve positive separator heights \(b\), on piecewise affine
rational curves in the lower Cayley scale. The upper scale divides
\(b\) by \(\sqrt r\). On the side of a ray with original modulus
\(\varrho\le1\), put
\(a=(1-\varrho)/(\lambda(1+\varrho))\). The inverse ray is
\[
                    Z=\frac{a+i b}{1+i h b a}
\]
in the order's own scale. If the maximum surviving tanh factor is
bounded by \(\alpha-1\), the annulus gives
\(\varrho\ge(1-h)^\alpha\ge 1-\alpha h\), hence
\[
                    0\le a\le\frac{\alpha\lambda}{2-\alpha h}.
\]
This cutoff has a continuous value at \(p=0\). A column reversal
and conjugation cover the side with \(\varrho\ge1\).

The checker uses weight-nine grouped coefficients. It verifies leading
dominance on the radius-three disk and excludes zeros on the inverse
ray cells. At the limit \(p=0\), exact signs at the reference heights
account for alternating positive heights of orders fourteen and
thirteen, starting and ending at order fourteen. The order-thirteen
structural real root lies at zero in the scaled variable. The joined
sector at the negative ray contains one odd-order root and two
even-order roots. Conjugation makes the odd-order root real. The
even-order factorial moment is strictly positive on the real axis
by its Gamma integral, so the even-order roots are nonreal.
The sector counts then give disjoint root arguments.

## Disks for \(.76\le u\le .98\)

The disk checker keeps all log and exponential terms of total weight
at most eleven. It uses the same grouped expansion and records the
first omitted even and odd logarithmic weights as twelve and
thirteen. The exponential remainder starts at the terms of weight
at least twelve among exponential powers up to five. The terms of
higher powers are bounded by the geometric estimate starting with
\(V^6/6!\), if \(V\) bounds the retained logarithm. For each retained
power of that logarithm, a truncated convolution keeps each weight
below twelve and the sum at weights twelve and above. Multiplication
by a nonnegative-weight factor can never return that last sum to the
retained weights, so its update gives exactly the omitted sum.
Every retained
direction is transformed as a whole before division by \(\lambda^j\).

Use \(P=M_2/(2p^2)\) at the lower order and its transported upper
interval from
\[
            P'=(R/r)^2 P+\Delta .
\]
The [matched tanh transport](MULTIVARIATE_MATCHED_TANH_MOMENTS.md#exact-tanh-transport-and-moment-errors-at-every-order)
bounds \(|\Delta|\) by
\(13^4 C^3 p_+/((168-T_*)^3r_-^2)\), where \(T_*\) is the
lower tanh cap and \(r_-\) is the lower deficit-ratio endpoint.
Here \(\sum T_i^3\le C^3 p^3\). For a positive lower deficit
the [excess-loss bound](MULTIVARIATE_TANH_EXCESS_LOSS_TRANSPORT.md)
gives the additional upper estimate
\[
  \frac{R^2 (56/12) B E^2 (L-R) U_+}{(p_-r_-)^2},
 \qquad
 B=\frac{168c}{26+168c},\quad E=\frac{12(1+14c)}{26+168c},
\]
where \(U_+\) is an outward upper endpoint for \(U\).
The checker uses the smaller of the two estimates.
The [conditional higher moment inequalities](MULTIVARIATE_SKEW_CONDITIONED_TANH_MOMENTS.md)
then carry the chosen signed skew intervals into the weight-eleven
enclosures. Column reversal allows the lower interval to be
restricted to \(P\ge0\).

For each interval cell, use roots of a representative transformed
polynomial to propose disk centers. For a proposed center \(z_0\)
write the actual shifted polynomial as \(\sum_k c_k z^k\).
Outward complex ball arithmetic verifies
\[
        \eta\,|c_1|_{\min}
             > |c_0|_{\max}+\sum_{k\ge2}\eta^k |c_k|_{\max}.
\]
Rouché's theorem then places exactly one root in the disk of
radius \(\eta\). The grouped polynomials are shifted before applying
the uncertain scalar parameters, preserving their cancellations.

For a disk in the upper half-plane, enclose
\[
                  \frac{\sqrt{h/p}\,\Im Z}{1-h|Z|^2}.
\]
The checker verifies that the denominator and imaginary coordinate
are positive. The ratio is the tangent of the displacement of the
original root argument from the negative ray, divided by the common
positive factor \(2\sqrt p\). Disjoint intervals therefore separate
the upper arguments. Conjugation accounts for the lower roots.
The one disk centered on the real axis at order thirteen contains
a real root by uniqueness and conjugation. All roots have now been
accounted for, so this real root has no even-order counterpart.

## Exact replay

The [ray checker](verify_multivariate_terminal_fourteen_cayley_flat.py)
uses 58,064 projection cells across 96 order, ray and subband checks.
The [disk checker](verify_multivariate_terminal_fourteen_cayley_disks.py)
uses 20,626 outward interval cells.

~~~sh
python3 research/factorial/image-rigidity/verify_multivariate_terminal_fourteen_cayley_flat.py
.venv/bin/python research/factorial/image-rigidity/verify_multivariate_terminal_fourteen_cayley_disks.py
~~~
