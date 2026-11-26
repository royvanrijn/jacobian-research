# Shifted factorial windows on common rays with doubled degree blocks

Take \(s\) distinct nonnegative integer columns \(b_i\) and add the
same positive integral vector \(T a\) to all of them, where each
\(a_\rho>0\). Partition them by the offset total degrees \(d_i=|b_i|\).
Suppose each degree block has at most two columns. Let \(b\) be
the number of degree blocks and \(h=s-b\) the number of doubled
blocks. In a doubled block \(j\), orient the difference of its
two offset columns as \(v_j\), and denote its degree by \(d_j\).
For \(B=(B_1,\ldots,B_h)\) define the quadratic form
\[
 V_t(B)=\sum_\rho \frac1{a_\rho}
                 \left(\sum_{j=1}^h B_j t^{d_j} v_{j,\rho}\right)^2 .
\]

**Theorem.** If \(h=1\), every prescribed set of \(s\) distinct
positive orders containing an even order has no common nonzero
coefficient direction for all sufficiently large \(T\). In
particular this covers each fixed consecutive support-size
window at an arbitrary start. The height threshold can depend
on the chosen orders.

For arbitrary \(h\ge1\), let \(I\) consist of \(s\) distinct
positive orders. Suppose it contains \(h\) even orders \(E\)
satisfying \(t/2>\max(4,2^{s-1})\) and such that \(V_t(B)=0\)
for \(t\in E\) forces \(B=0\). Then the same conclusion holds
for \(I\), with the other orders unrestricted. The condition on
the variance forms holds for any \(h\) distinct positive
orders if the pair differences are pairwise orthogonal in the
real metric with weights \(1/a_\rho\), if they are all collinear,
or if \(h\le2\). This gives arbitrary support size with arbitrarily
many doubled degree blocks.

For a fixed \(I\) in the theorem, the stripped homogeneous
resultant is a nonzero rational polynomial in \(T\), so the
finite exceptional heights are effectively bounded. The
statement holds over characteristic-zero fields.

## The common-height Gaussian reduction

Use the exact [tilted Gamma identity](SFC_TERNARY_COMMON_TRANSLATION_CUBIC_OBSTRUCTION.md#exact-tilted-gamma-normalization).
Divide out the scalar factorials at the common height and change
coefficient coordinates, independently of order, by
\[
              C_i=c_i T^{d_i}\prod_\rho a_\rho^{b_{\rho i}}.
\]
At order \(t\), the normalized moment is the \(t\)-th moment
of \(\sum_i C_i t^{d_i} Y^{b_i}\), with independent coordinates
\(Y_\rho\) equal to a unit-rate Gamma variable of shape
\(tTa_\rho+1\) divided by \(tTa_\rho\). Joint moments of
\(\sqrt T(Y-\mathbf 1)\) converge to those of centered real
Gaussians with variances \(1/(t a_\rho)\), by the elementary
cumulant formula for Gamma variables. All orders in \(I\) are
fixed in this argument.

Suppose a sequence \(T\to\infty\) has common projective zeros.
Set \(\max_i|C_i|=1\) and extract a limit. Write \(S_j\)
for the sum of the coefficients in each degree block \(j\).
Let \(p_t=\sum_j S_j t^{d_j}\) be the value of the normalized
phase at \(\mathbf 1\). In each fixed \(L^k\) norm the
difference from \(p_t\) is \(O(T^{-1/2})\), by the Taylor
expansion and the bounded joint moments above. If
\(\sqrt T |p_t|\) were unbounded, division of the vanishing
moment by \(p_t^t\) along a subsequence would give a limit
of one. Thus \(p_t=O(T^{-1/2})\) for every \(t\in I\).
The evaluation matrix \((t^{d_j})\) has full column rank,
as any \(b\) distinct positive orders give an invertible
matrix by Descartes' sign rule applied to a nonzero real
or imaginary part. Therefore \(S_j=O(T^{-1/2})\).

Pass to limits \(A_j=\lim\sqrt T S_j\). Each singleton
coefficient tends to zero; the two coefficients in doubled
block \(j\) tend to \((B_j,-B_j)\), with at least one
\(B_j\ne0\). The normalized phase multiplied by \(\sqrt T\)
now converges in joint moments at order \(t\) to a complex
Gaussian linear phase of mean
\(P_t=\sum_j A_j t^{d_j}\) and bilinear variance \(V_t(B)/t\).
Indeed the Taylor terms of order at least two vanish in this
scaling, and the gradient only retains the paired differences.
The limiting \(t\)-th moment must vanish.

Write \(\mathrm{He}_t\) for the monic probabilists' Hermite
polynomial with generating function \(e^{xz-z^2/2}\).
When \(V_t\ne0\), the limiting equation is
\(\mathrm{He}_t(P_t/\sqrt{-V_t/t})=0\). When \(V_t=0\)
the moment equals \(P_t^t\), so \(P_t=0\). In both cases
there is a squared root \(\eta_t\) of \(\mathrm{He}_t\)
such that
\[
                         t P_t^2+\eta_t V_t=0.                 \tag{1}
\]
Only odd orders can have \(\eta_t=0\). If \(V_t=P_t=0\) and
\(t\ge2\), either a zero or a nonzero permitted root can
be chosen.

## A single doubled block

Let its offset degree be \(d_0\). Then
\(V_t=\kappa B_1^2 t^{2d_0}\) with
\(\kappa=\sum_\rho v_{1,\rho}^2/a_\rho>0\).
It follows that
\[
              x_t=\frac{P_t}{i\sqrt\kappa\,B_1t^{d_0}}
\]
is a root of \(\mathrm{He}_t(\sqrt t\,x)\) at each order.
The evaluation matrix with entries \(t^{d_j-d_0}\) has
\(s\) rows and \(s-1\) columns, and each maximal minor is
nonzero by the same positive-order interpolation argument.
Thus a left null vector exists over \(\mathbb Q\) with
every coordinate nonzero, giving a rational linear relation
on the \(x_t\). An even order contributes a nonzero root.
The [scaled Hermite isolation theorem](SFC_SCALED_HERMITE_ROOT_ISOLATION.md)
excludes the relation. This proves the first claim.

## Several doubled blocks and a polynomial isolation argument

Consider the homogeneous quadrics
\[
                 t P_t(A)^2+X_t V_t(B) \qquad(t\in I)
\]
in the \(s=b+h\) variables \((A,B)\). Let \(\mathcal R(X)\)
be their homogeneous resultant. By multihomogeneity its
degree in each \(X_t\) is at most \(2^{s-1}\).
This is the resultant multidegree formula of
[Sturmfels, Theorem 4.4](https://www.math.ucdavis.edu/~deloera/MISC/LA-BIBLIO/trunk/Sturmfels/cbms.pdf).

Keep the actual values \(\eta_t\) fixed at all indices outside
\(E\), including any odd zero roots. The specialized resultant
in the \(X_t\) for \(t\in E\) is a nonzero polynomial over
the field of those other values. To see this, write the
equations at \(t\in E\) as \(V_t(B)+Y_t\,tP_t(A)^2=0\),
where \(Y_t=1/X_t\). At \(Y_t=0\) for every \(t\in E\),
variance regularity forces \(B=0\). The other \(b\) original
equations then give \(P_t(A)=0\) at distinct positive orders,
forcing \(A=0\). The resultant in the \(Y_t\) is nonzero at
zero. Multihomogeneity relates it to the specialized resultant
in the \(X_t\) by multiplication by
\(\prod_{t\in E}X_t^{2^{s-1}}\) and substitution \(Y_t=1/X_t\).

The squared Hermite values \(\eta_t\) for \(t\in E\) are
nonzero since \(t\) is even. The
[squared-root field isolation](SFC_SCALED_HERMITE_ROOT_ISOLATION.md)
gives degree \(t/2\) for each one over the field generated
by all the other chosen values. In particular its degree
exceeds \(2^{s-1}\). A nonzero polynomial over the field
of values outside \(E\), of separate degrees at most
\(2^{s-1}\), cannot vanish on the values at \(E\). Regard
it as a polynomial in one variable over the field of the
other values. Vanishing would force every coefficient
evaluation to be zero, by the degree at that variable.
Induct on the number of variables, using the same lower
bound over smaller fields. This contradicts vanishing of
the specialized resultant at (1).

For orthogonal differences \(V_t\) is a sum of
\(\kappa_j B_j^2 t^{2d_j}\) with \(\kappa_j>0\), so \(h\)
distinct orders force all \(B_j=0\). For collinear differences
\(V_t\) is a nonzero positive scalar times the square of a
generalized polynomial with coefficients nonzero multiples
of \(B_j\). The same interpolation applies. When \(h=2\)
and both variables are nonzero, a tested zero requires
\(t^{d_1-d_2} B_1/B_2\) to equal one of the roots of the
fixed real Gram quadratic. The Gram matrix is positive
semidefinite and has both diagonal entries positive.
Its roots are either a conjugate pair of equal nonzero
modulus or a repeated real root of nonzero modulus. Two
distinct positive \(t\) give distinct moduli, so they
cannot both vanish. If one variable is zero there is no
tested zero unless both are zero.

Finally, each factorial moment at a fixed order, after
stripping the factorials at the common height, has
polynomial coefficients in \(T\). The regularity just proved
at large ordinary \(T\) shows its window resultant is not
identically zero. Its integer roots bound the exceptional
heights and a standard root bound, or exact factorization,
makes this effective for the prescribed window. A finitely
generated field of coefficients over \(\mathbb Q\) embeds
in \(\mathbb C\), preserving moment vanishing.
