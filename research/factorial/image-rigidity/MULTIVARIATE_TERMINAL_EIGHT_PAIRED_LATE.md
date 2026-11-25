# Shared terminal-eight separation from outer mean .72 to .76

Let all unequal rows of distinct ordinary exponent columns be shared.
Let \(u\) denote their normalized order-seven outer geometric mean
and put \(\ell_i=\min(a_i,b_i)\), \(d_i=|b_i-a_i|\) on unequal rows.
When
\[
          .72\le u\le .76,\qquad\max_i d_i/(7\ell_i+1)\le1/25 ,
\]
the order-seven and order-eight factorial moment polynomials have
disjoint root arguments under independent positive rescalings.
Their rational Bézout identity gives coprimality over every
characteristic-zero field.

## Additional linked coefficient terms

Use the [shared skew transfer](MULTIVARIATE_TERMINAL_EIGHT_PAIRED_LIGHT.md)
and the same Gaussian midpoint models. In the notation of the
[octic linked means](MULTIVARIATE_TERMINAL_EIGHT_NEAR_FLAT.md#octic-linked-means-and-skews),
retain the combination
\[
       J_4=\tfrac85J_3+R,\qquad
       R=\sum\big(G_{1/7}-\tfrac35G_{3/7}+\tfrac15G_{5/7}\big)(s).
\]
For one factor write \(v=\tanh s\) and \(D_a=1-a^2v^2\).
The derivative of the summand is
\[
       \frac{48\cdot24\cdot8}{343\cdot49^2}
         \frac{v^6}{D_{5/7}D_{3/7}D_{1/7}} .
\]
It is nonnegative. Use \(v\le s\le s_t\), with \(s_t\) from the
reciprocal cap, and divide the seventh-power bound by \(s_t\)
to sum across the factors. This keeps \(J_3\) as one shared
signed contribution to both the paired and central means.

At each reference coefficient separate the log correction as
\(\sigma=P+B\), with \(P\) the sum of the primary inner skew and
the common log-mean shift terms at that coefficient. The
remaining terms are in \(B\). Expand
\[
 e^\sigma=1+\sigma+\tfrac12 P^2+
              P B+\tfrac12 B^2+T,\qquad
 |T|\le\frac{|\sigma|^3}{6(1-|\sigma|/4)} .
\]
The square of either shared scalar in \(P\), and their product,
remain shared across all coefficients. Bound \(PB+B^2/2+T\)
in each coefficient using its interval. The exponential
remainder bound follows from the Taylor series and the
geometric bound on its tail.

## Exact projections and endpoint

The [late checker](verify_multivariate_terminal_eight_paired_late.py)
uses the exact cosine Bernstein evaluation and grouped radial
projections of the [preceding bridge](MULTIVARIATE_TERMINAL_EIGHT_PAIRED_BRIDGE.md),
with rational affine ray cosines. All actual and nominal interior
multipliers are at most one, so the same annulus
\(1/12<|t|<12\) applies. The logarithmic coefficient homotopy
scales down all skews and defects at fixed means. It preserves
the linear skew and defect constraints. Then vary the
outer-log ratio to \(49/48\) and the Gaussian outer mean to .72.
The six ray cosines at .72 match the preceding checker exactly,
so its sector counts apply. The last joined sector has one real
septic root and one nonreal octic pair, by conjugacy and the
strict positivity of the ordinary eighth factorial moment.
Thus the arguments are disjoint.

Replay from the repository root with

~~~sh
python3 -O research/factorial/image-rigidity/verify_multivariate_terminal_eight_paired_late.py
~~~
