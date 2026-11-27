# Initial windows at all heights for homogeneous trinomials of short span

**Theorem.** A nonzero homogeneous bivariate trinomial
whose exponents in one row have span at most 12 cannot
have all of its first three factorial moments zero,
over any characteristic-zero field. All independent
common monomial translations are included.

Write the offsets as \((0,d),(a,d-a),(d,0)\),
with \(0<a<d\). The least exponents in the two rows are
independent nonnegative integers \(L,H\). Reflecting the
two coordinates sends \(a\) to \(d-a\) and exchanges
\(L,H\). Thus it suffices to consider
\(1\le a\le \lfloor d/2\rfloor\).

At each order \(m\) strip \((mL)!(mH)!\) and form the
integer homogeneous moment \(P_m\). The coefficient
of the count vector \(\alpha\) is
\[
            \binom m\alpha (mL+1)_{U}(mH+1)_{md-U},
              \qquad U=a\alpha_1+d\alpha_2 .
\]
Here \((t)_k\) is the rising factorial. Let \(p_0,p_1,p_2\)
be the first-order coefficients. Substitute
\((p_2 X,p_2 Y,-p_0 X-p_1Y)\) for the column coefficients,
which spans the first zero plane since \(p_2>0\).
Write the new second and third forms as
\[
       A X^2+BXY+CY^2,\qquad
             D X^3+E X^2Y+F XY^2+G Y^3 .
\]
The second is positive definite over the reals for all
ordinary heights. A real nonzero coefficient direction
gives a nonzero polynomial and hence has a strictly
positive square integral in the product-exponential
measure. If there is a common projective zero over
\(\mathbb C\), the cubic is divisible by the quadratic
since its two nonreal projective roots are conjugate.
Comparison of coefficients forces
\[
                          ACF-C^2D-ABG=0.                    \tag{1}
\]
The same exclusion applies over any characteristic-zero
field by the Nullstellensatz over the rational
coefficient chart.

The [exact checker](verify_sfc_homogeneous_trinomial_short_spans.py)
forms (1) in \(\mathbb Z[L,H]\) for the 36 reflected
profiles with \(2\le d\le12\). Every nonzero coefficient is
negative and the constant coefficient is negative in
each case. This excludes every nonnegative height pair
without a height cutoff or a numerical search over heights.
The same test can be applied to further spans. The
coefficient sign for unrestricted spans is a sufficient
target that remains unproved.

Replay using Python 3 and SymPy from the repository root.

~~~sh
python3 research/factorial/image-rigidity/verify_sfc_homogeneous_trinomial_short_spans.py
~~~
