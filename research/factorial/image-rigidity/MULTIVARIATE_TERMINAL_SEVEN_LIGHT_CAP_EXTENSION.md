# Terminal-seven separation on the light reciprocal cap 1/8 through outer mean .48

Let all unequal rows of two distinct ordinary exponent columns
be shared, and on these rows put
\(\ell_i=\min(a_i,b_i)\), \(d_i=|b_i-a_i|\).
Let \(u\) be the order-six normalized outer geometric mean.
When
\[
       .10\le u\le.48,\qquad \max_i d_i/(6\ell_i+1)\le1/8,
\]
the order-six and order-seven factorial polynomials have disjoint
root arguments up to independent positive rescalings. Rational
Bézout gives coprimality over every characteristic-zero field.
This widens the reciprocal caps of the
[paired light-row proof](MULTIVARIATE_TERMINAL_SEVEN_PAIRED_LIGHT.md)
on the indicated outer-mean range. The
[small-outer-mean proof](MULTIVARIATE_TERMINAL_SEVEN_SMALL_OUTER.md)
covers \(0<u\le.10\) at the same cap.

Apply the paired model with \(t_6=1/8\) and
\(t_7=6t_6/(7-t_6)\). Use a sharper bound on the
single order-seven inner skew \(\gamma\).
For local \(x'\) put \(g_h=\log(1+h x')\), \(s=g_6/2\).
The magnitude of that skew is
\[
 |\gamma_{\rm local}|=\frac{6g_3-g_1-g_2-g_4-g_5-g_6}{14}.
\]
It is nonnegative by the centered concavity of the logarithm.
The exact skew-convexity check of the
[one-private-linear proof](MULTIVARIATE_TERMINAL_SEVEN_PRIVATE_LINEAR.md)
verifies that this magnitude is convex in \(s\) for
\(x'\in[0,1]\). It vanishes at zero, so its ratio to \(s\)
increases. Evaluating the ratio at \(x'=t_7\) and summing
bounds \(|\gamma|/W\), where \(W=-\log w\).
The formula is \((6g_3-g_1-g_2-g_4-g_5-g_6)/(7g_6)\)
at that endpoint. The checker computes its enclosure by the
same outward rational logarithm bounds as the paired proof.

Keep the other linked coefficient enclosures of the paired
proof unchanged. The [extension checker](verify_multivariate_terminal_seven_light_cap_extension.py)
verifies the new ray exclusions across the larger reciprocal
cap. Only the first ray has modified knots, all at \(u\ge.36\).
The annulus and affine coefficient homotopy of the paired proof
are valid under these enclosures. Along the symmetric arrays,
move first to \(W/U=36/35\), then to \(u=.10\), where the rays
match the earlier Gaussian anchor exactly. The sector counts
and the real-axis positivity argument there apply without
change.

Replay from the repository root with

~~~sh
python3 -O research/factorial/image-rigidity/verify_multivariate_terminal_seven_light_cap_extension.py
~~~
