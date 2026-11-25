# Shared terminal-eight separation from outer mean .70 to .72

Let all unequal rows of two distinct ordinary exponent columns be
shared. Put \(\ell_i=\min(a_i,b_i)\) and \(d_i=|b_i-a_i|\) on
those rows, and let \(u\) be the normalized order-seven outer
geometric mean. When
\[
           .70\le u\le .72,\qquad
           \max_i d_i/(7\ell_i+1)\le9/200 ,
\]
the order-seven and order-eight factorial moment polynomials
have disjoint root arguments under independent positive
rescalings. Their rational Bézout identity proves coprimality
over every characteristic-zero field.

The [exact checker](verify_multivariate_terminal_eight_paired_bridge_72.py)
uses the coefficient model and cross-order skew constraint from the
[paired light proof](MULTIVARIATE_TERMINAL_EIGHT_PAIRED_LIGHT.md)
and the cosine and radial Bernstein projections of the
[preceding bridge](MULTIVARIATE_TERMINAL_EIGHT_PAIRED_BRIDGE.md).
The six ray cosines are rational affine functions of \(u\),
starting at the preceding bridge's values at .70 and ending
at the values displayed in the checker at .72.
The same root annulus and logarithmic interpolation give
root counts by continuation through the .70 boundary. The last
joined sector contains a real septic root and a nonreal octic
pair, by conjugacy and real-axis positivity of the ordinary
eighth factorial moment. The root arguments are disjoint.

The twelve strict ray checks give cell counts
\(4575,243,4048,433,1288,108\) at order seven and
\(4010,456,22420,3154,18358,2619\) at order eight.
Replay from the repository root with

~~~sh
python3 -O research/factorial/image-rigidity/verify_multivariate_terminal_eight_paired_bridge_72.py
~~~
