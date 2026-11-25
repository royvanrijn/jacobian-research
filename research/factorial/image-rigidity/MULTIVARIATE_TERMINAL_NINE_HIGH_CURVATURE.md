# Terminal nine at total curvature at least 1493/200

For distinct ordinary exponent columns \(a,b\), put
\(S=\sum_{i:a_i\ne b_i}(b_i-a_i)^2/\max(a_i,b_i)\).
If \(S\ge1493/200\), the order-eight and order-nine factorial
moments have disjoint complex root arguments. Consequently the
two polynomials are coprime over every characteristic-zero field.
The degree difference, the signs of the row gaps, and private
rows are unrestricted.

For \(j=8,9\), normalize both endpoint coefficients to one.
The [total-curvature estimate](MULTIVARIATE_FIXED_TERMINAL_CURVATURE.md#large-curvature-and-the-endpoint-radii)
bounds each interior coefficient by
\[
               \binom jk \exp\bigl(-1493 k(j-k)/(400j)\bigr).
\]
The [rational checker](verify_multivariate_terminal_nine_high_curvature.py)
encloses this by raising it to order \(j\), splitting the inverse
exponential into two equal factors, and applying the positive
Taylor lower bound to each factor. The checker then takes an
outward rational \(j\)-th-root bound.
The resulting coefficient boxes contain the actual normalized
arrays and every straight segment from \(1+t^j\) to an actual array.
The constant term and the leading term strictly dominate on the
respective radii \(1/8\) and \(8\) throughout the boxes.

The ray cosines, in increasing angular order, are
\[
 \frac{23097}{25000},\quad \frac{3741}{5000},\quad
 \frac{4147}{10000},\quad \frac{199}{2000},\quad
 -\frac{377}{1250},\quad-\frac{317}{500},\quad
 -\frac{8729}{10000}.
\]
The checker applies the
[rational Bernstein projection tests](MULTIVARIATE_TERMINAL_SIX_JOINT_RAY_BANDS.md#moving-rays-and-exact-replay)
on each ray throughout the annulus and for both orders.
It proves strict nonvanishing in 92 radial cells.

For \(x=\cos\theta\), the sign of the binomial anchor at the
order-eight ray after removing its central phase is the sign of
\(8x^4-8x^2+1\). At order nine, also divide by
\(2\cos(\theta/2)>0\) for \(0<\theta<\pi\). The sign polynomial
is \(16x^4-8x^3-12x^2+4x+1\).
In the reusable checker these polynomials come from the common
cosine recurrence, starting with \(1,x\) in even degree and
\(1,2x-1\) in odd degree. Exact signs at the rays and at \(x=\pm1\)
give one upper anchor root per angular sector, in the order
\(9,8,9,8,9,8,9,8\). The remaining ninth-order root is negative real.

Continue the root counts along the coefficient segments.
The first and last sectors are joined to their lower reflections
for counting across the real axes. The last joined sector contains
one ninth-order root and two eighth-order roots.
Conjugation makes the single ninth-order root there negative real.
For the actual even moment, positivity of
\(\mathbb E[(X^a+tX^b)^8]\) at every real \(t\) excludes a real root,
so its two eighth-order roots occupy one upper and one lower position.
Thus the argument sets are disjoint. The monomial normalization
uses positive real scalars, and the rational Bézout identity
transfers coprimality to every characteristic-zero field.
