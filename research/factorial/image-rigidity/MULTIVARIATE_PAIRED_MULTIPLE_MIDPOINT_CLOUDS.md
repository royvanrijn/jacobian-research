# Paired multiple Gaussian contacts of infinite ordinary cores

A finite integer increase of a single retained midpoint gives an
exact rising product. Two ordinary row clouds cancel selected
terms of its logarithm, allowing adjacent Gaussian clocks to
split at different strengths. The selected row can have unit gap.

**Theorem.** In the
[fixed-terminal reference](MULTIVARIATE_FIXED_TERMINAL_COMPACTIFICATION.md#exact-compactified-data),
fix consecutive orders \(n,N=n+1\), \(n\ge2\), and a physical interior
contact of positive strength \(T_*\) and bias \(\lambda_*\).
Assume both contacting nonreal roots have multiplicity at least two
and the ordinary retained core contains infinitely many unequal
rows with summable curvature. Then the contact is a limit of simple
transverse interior Gaussian contacts of finite ordinary cores.
The orders and total signed degree difference agree, and the
strengths, biases, and relative physical root scales converge.
Every strict total curvature cap persists.

## The integer midpoint operator

Select weak retained rows of absolute gap \(d>0\) and midpoint \(M\),
with \(d^2/M\to0\). Their existence follows from summability.
After a subsequence, reverse both columns if necessary so that the
selected gaps are positive. Set \(\eta=d/M\). Fix an integer \(L>0\)
and increase both endpoints of the selected row by \(L\).
For frequency \(v=k-j/2\) at order \(j=n,N\), the exact ratio
of its unnormalized factorial coefficients is
\[
       \prod_{h=1}^{jL}\bigl(jM(1+\zeta)+h\bigr),
               \qquad \zeta=\eta v/j .
\]
After subtracting the constant \(jL\log(jM)\), its logarithm is
\(jL\log(1+\zeta)+E_j(\zeta)\). Here \(E_j\) is analytic and
\(O_{j,L}(M^{-1})\) on \(|\zeta|\le1/4\), by expanding the factors
\(1+h/[jM(1+\zeta)]\) there. Endpoint normalization removes the
affine chord in \(\zeta\), which multiplies the polynomial by a scalar
and effects a real radial change with logarithmic slope \(O_{j,L}(\eta)\).
The latter can be set aside when comparing angles.

Write \(\xi_j\) for the logarithm of the original normalized root,
with the two imaginary parts equal, and \(m_j\) for its multiplicity.
In the local coordinate \(X=(x-\xi_j)/\eta\), the limiting
operator on \(X^{m_j}\) is
\[
              e^{y\partial_X^2/(2j)}(1+\partial_X/j)^{jL},
                \qquad T=T_*+\eta^2y .                         \tag{1}
\]
The common strength increment acts as the indicated heat operator
on \(e^{-jx/2}F_j(T_*,e^x)\), up to a nonzero scalar.
The rising product acts as a polynomial in differentiation on the
same finite frequency set. To justify the limit, truncate its logarithm
on \(|\zeta|\le1/4\) and the heat multiplier on the finite frequency
set at orders high enough for their errors to be \(o(\eta^{m_j})\).
Under differentiation in \(X\), the fixed Taylor coefficients of
the rescaled \(E_j\) term are \(O(M^{-1})\) by Cauchy's estimate.
The unperturbed function divided by \(\eta^{m_j}\) tends to a
nonzero scalar multiple of \(X^{m_j}\) on compact disks, with
derivatives. These facts prove (1).

Translation in \(X\) by a real constant is immaterial for angle
comparisons. Suppressing this translation in (1), the additional
logarithmic terms from the midpoint shift start with
\[
         -\frac{L}{2j}\partial_X^2
          +\frac{L}{3j^2}\partial_X^3
          -\frac{L}{4j^3}\partial_X^4 .
                                                               \tag{2}
\]
For every fixed \(h\ge5\) the coefficient of \(\partial_X^h\)
has size \(O_j(L)\).

## Two exact ordinary clouds

Let \(w_C,w_R>0\) be fixed and small. For a required nonzero cubic
coefficient \(\beta/j^2\), add \(\lfloor w_C H_C\rfloor\) unit-gap rows
of orientation \(-\operatorname{sgn}\beta\), with half-integer
midpoint \(H_C\sim w_C/(6|\beta|\eta^3)\).
Transfer their actual variance \(w_C(\eta)=\lfloor w_C H_C\rfloor/H_C\)
from the Gaussian strength to the rows. The rescaled log
multiplier then contributes \(\beta\partial_X^3/j^2\).

Indeed, at inverse midpoint \(u\) the compensated coefficient
of \(v^3\) is \(-u^2/(6j^2)\) for an increasing unit row
and has the opposite sign for a decreasing row. The other
terms of that layer have degree at most two in \(v\).
At a layer \(u^h\), with \(h\ge3\), the log coefficient has
degree at most \(h+1\) in \(v\), by summing powers in the
finite factorial product and subtracting its endpoint chord.
Thus the rescaled coefficients from this layer are
\(O(\eta^{3(h-1)-(h+1)})\), which vanish.
The quadratic and lower terms at \(h=2\) also vanish.

For \(A>0\), one may additionally transfer a variance tending to
\(w_R\) to reflected unit-gap pairs. Use \(\lfloor w_R H_R/2\rfloor\)
pairs at half-integer midpoint \(H_R\sim 1/(A\eta^2)\).
Their actual variance is \(w_R(\eta)=2\lfloor w_R H_R/2\rfloor/H_R\).
The limiting log multiplier from these pairs is
\[
             -\frac{w_R A}{4j^2}\partial_X^2
              +\frac{w_R A^2}{12j^3}\partial_X^4 .             \tag{3}
\]
To check (3), the compensated log multiplier per pair is
\[
 u^2\phi(k/j)+u^3\bigl(2j\Xi(k/j)-\phi(k/j)/(3j)\bigr)+O_j(u^4),
 \quad
 \phi(s)=s(1-s)/2,\quad
 \Xi(s)=((s-1/2)^4-1/16)/12 ,
\]
as in the [reflected chord calculation](MULTIVARIATE_MULTIPLE_ROOT_INTERIOR_HEAT_UNFOLDING.md#exact-reflected-row-compensation).
The top coefficients in the displayed layers are
\(-v^2/(2j^2)\) and \(v^4/(6j^3)\).
The degree at layer \(u^h\) is even and at most
\(2\lfloor(h+1)/2\rfloor\). For \(h\ge4\) the exponent
\(2(h-1)-2\lfloor(h+1)/2\rfloor\) is positive.
Both cloud calculations use analytic log factorial products
for fixed unit gaps. Truncations through sufficiently high
finite order have remainder \(o(\eta^{m_j})\) on the
original finite frequencies, so the coefficient limits give
local operator limits as in (1). A cloud can be omitted by
leaving all its corresponding terms zero.

## The two clock choices

In either choice take \(\beta=-L/3\), canceling the midpoint
cubic at both orders. First include (3) and set
\(A=\sqrt{3L/w_R}\), canceling the midpoint quartic.
Let
\[
                y=L+\sqrt L\,t,\qquad X=L^{1/4}Y .
\]
After ignoring radial translations, the normalized model as
\(L\to\infty\) is
\[
 \exp\left(\frac{t-\sqrt{3w_R}/(2j)}{2j}\,
                      \partial_Y^2\right)Y^{m_j}.             \tag{4}
\]
Indeed the uncanceled terms of degree \(h\ge5\) from (2)
and its higher analogues are \(O_j(L^{1-h/4})\) after dilation.
The limit is in the finite-dimensional space of operators
acting on \(Y^{m_j}\).

For \(m\ge2\), let \(h_m^-,h_m^+\) be the smallest and largest
positive roots of \(\mathrm{He}_m\). Suppose
\[
                 h_{m_n}^+/\sqrt n>h_{m_N}^-/\sqrt N .         \tag{5}
\]
At a fixed \(t\) between \(\sqrt{3w_R}/(2N)\) and
\(\sqrt{3w_R}/(2n)\), the lower roots of (4) are all real, while
the upper inner positive imaginary rank is strictly positive.
Comparing the lower outer positive rank of the positive-heat
model with the upper inner positive rank gives a strictly
negative difference at this clock. At large positive \(t\)
the same difference is positive by (5).
Choosing a sufficiently large integer \(L\) preserves the
strict signs in the models before taking the limit in \(L\).

The rank classification in the
[unbounded-gap proof](MULTIVARIATE_PAIRED_MULTIPLE_UNBOUNDED_GAPS.md#the-outstanding-hermite-configurations)
shows that (5) can only fail for the lower double and upper odd
multiplicities \(3,5,7,9\). This follows from Hermite and
same-parity Laguerre interlacing and \(h_{11}^-<1\). Equality
in these four comparisons would make \(N/n\) a rational
nonintegral squared Hermite root, contrary to the monic
integer polynomial \(\mathrm{He}_m(X)/X\) in \(X^2\).
The double/triple case is already covered on infinite cores
by the [reverse heat shift](MULTIVARIATE_REVERSE_CLOCK_ROW_REPLACEMENT.md).

For the cases \(m_n=2\), \(m_N\in\{5,7,9\}\) with (5) reversed,
omit the reflected cloud. The limiting model at the same
\(y=L+\sqrt L\,t\) and \(X=L^{1/4}Y\) is
\[
             \exp\left(\frac{t}{2j}\partial_Y^2
                         -\frac{1}{4j^3}\partial_Y^4\right)
                                      Y^{m_j}.                 \tag{6}
\]
At \(t=0\) the upper polynomials, after a real dilation
and a nonzero scalar, are respectively
\[
        Z^5-30Z,\qquad
        Z^7-210Z^3,\qquad
        Z^9-756Z^5+11340Z .
\]
Each has at least three real sign-changing zeros.
For the last polynomial, the quadratic
\(q^2-756q+11340\) has two different positive roots.
At a sufficiently small fixed \(t>0\), the upper model
therefore retains at least three real zeros, while the
lower double model has a positive imaginary zero. Thus
the difference of the same lower and upper ranks is
positive. At sufficiently large positive \(t\), dilating
by \(\sqrt t\) leaves pure heat in (6), which makes this
rank difference negative by the reversed inequality in (5).
Once again some fixed sufficiently large integer \(L\)
preserves these strict inequalities.

## Angular transfer and physical parameters

For each fixed choice of the clocks, \(L\), and cloud
parameters, select the weak retained row sufficiently far
out so that \(\eta\) is small. Truncate the original retained
core after that row with coefficient errors \(o(\eta^{\max(m_n,m_N)})\)
and center errors tending to zero. The split limits proved
above persist on both selected clocks. In the two isolating
root disks the sorted imaginary parts of the logarithmic
roots are continuous. The disk boundaries can be kept free
of roots throughout the clock interval, since the changes
in the coefficient arrays tend to zero. Apply the
[ordinary regularization lemma](MULTIVARIATE_MULTIPLE_CONTACT_ANGULAR_RANK_REVERSAL.md#ordinary-perturbations-with-simple-transverse-angular-contacts)
with perturbations small enough to preserve the two strict
rank signs. It gives a simple transverse angular contact
between the clocks.

The midpoint shift keeps the row's signed gap and changes
its curvature and center by \(o(1)\). Maintain the total
signed degree difference by the residual Gaussian data
as in the fixed-terminal reference. The transferred actual
strength is \(T_*+\eta^2 y-w_C(\eta)-w_R(\eta)\), with the
omitted cloud assigned zero variance. The reflected center
contributions cancel. In orientation \(\epsilon\), the
unit cloud center tends to \(-\epsilon w_C/(2nN)\),
by the [center expansion](MULTIVARIATE_THIRD_ORDER_BOUNDARY_AND_INTERIOR_REDUCTION.md#exact-row-expansions-and-the-projected-coefficients).
The moduli of the two normalized roots tend to the original
moduli, so recovering the bias from their ratio gives
\(\lambda_*-\epsilon w_C+o(1)\). For
\(2w_C+w_R<T_*-|\lambda_*|\) the resulting contact remains
interior. The curvature transferred to the cloud rows
cancels the Gaussian decrease. Finally let \(w_C,w_R\to0\)
diagonally, choosing \(L\) and then \(\eta\) for each fixed
choice. This proves the convergence and the curvature
assertion in the theorem.
