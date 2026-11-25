# Matched shared rays from outer mean .44 to .48 at terminal nine

For distinct ordinary exponent columns with all unequal rows shared, let
\(u\) be the normalized order-eight outer reflected mean. Suppose
\[
 .44\le u\le .48,\qquad
 \max_i\frac{d_i}{8\ell_i+1}\le\frac18.
\]
Then the eighth and ninth factorial moment polynomials have disjoint root
arguments under independent positive rescalings and are coprime over
every characteristic-zero field.

The [exact extension checker](verify_multivariate_terminal_nine_shared_matched_extension.py)
uses the coefficient model, common-skew bounds, and homotopies from
the [preceding matched-light comparison](MULTIVARIATE_TERMINAL_NINE_SHARED_MATCHED_LIGHT.md).
The only new ray data are rational cosines at .46 and .48.
At .44 the rays coincide with the preceding checker by construction,
so the Gaussian reference arrays have the same sector counts there.
The two orders continue to the .44 references along the same ordered
rays and Gaussian homotopies. The checker establishes the annulus
\(1/12<|z|<12\) for the enlarged range and the strict Bernstein
projection signs. Positivity of the actual eighth moment excludes
real roots in the last sector joined to its conjugate.

Replay from the repository root with

~~~sh
python3 research/factorial/image-rigidity/verify_multivariate_terminal_nine_shared_matched_extension.py
~~~
