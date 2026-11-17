# Regular nilpotent B with its linear vector in ker B: retained partial reduction

17 November2026. **Classification remains OPEN.** This is a checkpoint
inside the full Image goal, not a closed boundary theorem. Normalize
\(H_2=zt\) and take \(H_1=az\), equivalently \(Bv=0\). The parameter
\(a\) is arbitrary and all seven cubic coefficients are retained in
[the exact generator](nilpotent_kernel_linear.py). Its input is extracted
from the unrestricted cyclic generator using the uniquely recoverable
linear-\(w\) exponent, not by deleting that fixed coefficient from an
already specialized polynomial. Direct spatial moments through5 were
checked during generation.

Under the stabilizing torus combined with amplitude, \(a\) is fixed
and \(\gamma_j\) has character \(j-1\). The complete nonpositive-charge
chart split is by the lowest nonzero cubic charge. The central-only bad
chart is excluded by \(\mu_2\). Positive-charge-only cubics give the
already safe one-sided family.

For lowest charge\(-1\), normalize \(\gamma_{-1}=1\), with lower
coefficients zero. The first two moments eliminate \(\gamma_0,\gamma_1\).
The remaining variables are \((a,\gamma_2,\gamma_3)\). The primitive
fourth and fifth equations have the constant-pivot form
\(143a^2-2\gamma_3+39\) and \(\gamma_2\). Eliminating them leaves

\[
 P_6=a(841a^2+285),\qquad
 P_8=44235a^4+39114a^2+6991.
\]

They are coprime. The [exported rational identity](nilpotent_kernel_lowest_minus1_unit.json)
contains their exact Bezout multipliers and the lifted seven-row
polynomial unit identity in moments4..10. Only indices through8 have
nonzero multipliers. The identity was directly multiplied exactly.

For lowest charge\(-2\), normalize \(\gamma_{-2}=1\) and
\(\gamma_{-3}=0\). Constant pivots give

\[
 \gamma_1=-\gamma_0\gamma_{-1}-6a,\qquad
 \gamma_2=-\frac{15\gamma_0^2-5\gamma_0\gamma_{-1}^2-30a\gamma_{-1}}4.
\]

Four variables remain. The [input](nilpotent_kernel_linear_charts.json),
[rational multiplier export](nilpotent_kernel_lowest_minus2_unit_witness.log)
and [direct FLINT checker](nilpotent_kernel_unit_check.py) certify a unit
identity through10. Although the coefficients have thousands of decimal
digits, there are only261 multiplier monomials. Exact sparse rational
multiplication takes less than a second and yields1. The checker parses
polynomial data without evaluating a giant expression. The proof does
not depend on trusting the standard-basis oracle that originally found
the multipliers.

The **remaining lowest\(-3\) chart** is normalized by
\(\gamma_{-3}=1\), with

\[
\begin{aligned}
 \gamma_2&=-\bigl(\gamma_0\gamma_{-1}+(\gamma_1+6a)\gamma_{-2}\bigr)/2,\\
 \gamma_3&=-\bigl(15\gamma_0^2+5\gamma_1\gamma_{-1}
       -2(\gamma_0\gamma_{-1}+(\gamma_1+6a)\gamma_{-2})\gamma_{-2}\bigr)/12.
\end{aligned}
\]

Its five variables \((a,\gamma_0,\gamma_1,\gamma_{-1},\gamma_{-2})\)
remain unrestricted. Exact input through12 is retained in
[the extension](nilpotent_kernel_linear_charts_prefix12.json). A modular
prefix10 basis has length18 and forces
\(\gamma_0=\gamma_{-1}=\gamma_{-2}=0\) in that special fibre. The
resulting **selected** two-variable family
\(az+zt+w^3+\gamma_1z(t^2-zw/2)\) is excluded by its exact rational
moments4,8,12. That selected exclusion does not exhaust the five-variable
chart over characteristic zero.

Full QQ std/slimgb runs stopped at35/60/90 seconds. A tracked modular
prefix12 calculation at1,073,741,827 returned a unit in about6 seconds,
but its multipliers have about183,000 monomials. The complete source,
7.8MB output and run receipt are preserved. There is no rational
reconstruction or rational unit identity for this last chart. Neither
the modular unit nor its prefix10 support reduction proves a rational
exclusion. Do not blindly restart the stopped rational algorithms.

A subsequent bounded weighted-homogeneous modular reconstruction also
stopped at90 seconds, without a reconstructed basis. Its source, output
and receipt are retained as
`nilpotent_kernel_lowest_minus3_prefix10_weighted_modstd.*`. It adds no
rational exclusion. No process from these experiments remains live.

The [new contraction-zero theorem](NILPOTENT_CONTRACTION_ZERO_AND_THREE_GLOBAL_TARGETS.md)
closes this boundary under the declared additional condition
Delta^2(H2*H3)=0, as part of an exhaustive arbitrary-linear nilpotent
classification. It does not assert that an unrestricted pure phase
automatically satisfies that condition.

If the final chart is excluded, the entire ker\(B\) family would be
one-sided whenever pure. That conclusion is **not yet proved**. The
[complete highest-cubic finite-module theorem](HIGHEST_CUBIC_GLOBALIZATION_AND_GRADED_GATE.md)
now supplies an alternative global framework including this boundary.
