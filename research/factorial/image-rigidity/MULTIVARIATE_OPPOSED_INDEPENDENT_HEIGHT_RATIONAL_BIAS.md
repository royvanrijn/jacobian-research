# Opposed independent heights at transverse rational-bias contacts

Fix an ordinary finite binary core \(C\), positive integer gaps
\(d_+,d_-\), and adjacent orders \(n,N=n+1\) with \(n\ge2\).
Append \(r_+\) copies of the increasing row \((L_+,L_++d_+)\)
and \(r_-\) copies of the decreasing row \((L_-+d_-,L_-)\).
The ordinary smaller exponents \(L_+,L_-\) and the counts may
vary independently. Write \(M_\pm=L_\pm+d_\pm/2\),
\[
  x=\frac{r_+d_+^2}{M_+},\qquad y=\frac{r_-d_-^2}{M_-},\qquad
  T=x+y,\qquad \lambda=x/d_+-y/d_-,\qquad
  h=d_+r_+-d_-r_- .
\]
These midpoint parameters agree with the compactified parameters
at an escaping limit. Fix \(h\in\mathbb Z\).

**Theorem.** Let \((T_*,\lambda_*)\) be a positive physical
Gaussian reference contact of the core with twist \(h\).
Assume
\[
   -1/d_-<\lambda_*/T_*<1/d_+,\qquad
                    \lambda_*/T_*\in\mathbb Q .
\]
Suppose the two specified contacting roots are simple and the
derivative in \(T\) of their angular difference is nonzero.
Then there is no sequence of ordinary collisions continuing
this pair with \(L_\pm\to\infty\) and
\((T,\lambda)\to(T_*,\lambda_*)\). In particular there is a
neighborhood of that contact excluding the original physical
collisions of all such sufficiently high two-type tails.

This result places no initial relation between the two heights
or the two counts beyond the fixed integer signed gap. For
\(h\ne0\), the first non-Gaussian correction already excludes a
collision by integrality. For \(h=0\), the same correction reduces
any collision sequence to an exact affine relation between the
heights. An isolated first factorial factor or a shared
first-factor chord excludes that affine family.

**Affine consequence.** For any fixed two opposed affine smaller
heights with positive rational slopes and rational offsets,
there is a computable height threshold above which every nonempty
tail with arbitrary independent counts is physically coprime in
every degree. Here the parameter is integral and only ordinary
heights are used. The conclusion allows unequal gaps, including
at a shared affine smaller height.

## Analytic equations in the independent curvatures

Use the
[fixed-terminal endpoint and center normalization](MULTIVARIATE_FIXED_TERMINAL_COMPACTIFICATION.md#normalized-polynomials-and-row-data).
For the core alone let \(G_j^C\) be the positive \(j\)-th root
of its endpoint coefficient ratio, and let
\[
 F_j(T,z)=\sum_{k=0}^j a_{j,k}
             \exp[-T k(j-k)/(2j)]z^k,\qquad j=n,N
\]
with \(a_{j,k}\) its endpoint-normalized coefficients.
The scale ratio of the reference is
\[
 \exp(C_R-\lambda/(2nN)),\qquad
       C_R=\log(G_N^C/G_n^C)+h\log(N/n).
\]
At the specified contact continue the two simple nonzero roots
in their respective endpoint variables, with logarithms whose
difference is real at the contact. For a perturbation vector
\(V\), multiply the coefficients of \(F_j\) by the exponentials
of the coordinates \(V_{j,k}\), and write the new root branches
as \(z_j(T,V)\). The last coordinate \(V_E\) accounts for the
residual center change. The exact contact equation is
\[
 H(T,\lambda,V)=\log(z_N(T,V)/z_n(T,V))
                   -C_R+\lambda/(2nN)-V_E=0 .                 \tag{1}
\]
The two real equations in (1) are nonsingular in \(T,\lambda\)
at the reference. In fact \(\partial_\lambda H=1/(2nN)\)
is real and the imaginary part of \(\partial_T H\) is nonzero
by hypothesis.

For a row of signed gap \(d\), inverse midpoint \(u\), and
order \(j\), define
\[
 R_{j,k}^{(d)}(u)
      =\sum_{t=1}^{kd}\log\bigl(1+(t/j-d/2)u\bigr).
\]
When \(kd<0\), the sum means the negative of the sum from
\(kd+1\) to zero. The finite factorial quotient from index
zero to index \(k\) is
\((j/u)^{kd}\exp R_{j,k}^{(d)}(u)\). Thus the per-row
compensated vector \(V^{(d)}(u)\) has entries
\[
\begin{split}
 V_{j,k}^{(d)}(u)
   &=R_{j,k}^{(d)}(u)-\frac{k}{j}R_{j,j}^{(d)}(u)
                         +\frac{d^2u\,k(j-k)}{2j},\
 V_E^{(d)}(u)
   &=\frac1N R_{N,N}^{(d)}(u)-\frac1n R_{n,n}^{(d)}(u)
                         +\frac{du}{2nN}.
\end{split}
\]
The signed-sum identity \(\sum_{t=1}^{kd}t=kd(kd+1)/2\)
gives
\[
                    V^{(d)}(u)=A^{(d)}u^2+O(u^3)               \tag{2}
\]
with a rational vector \(A^{(d)}\). The analytic formulas
hold also when one signed sum is negative. They show
that (1) is exact for the original ordinary pair on
substituting the actual \(T,\lambda\) and the vector
\[
                 V=r_+V^{(d_+)}(1/M_+)
                         +r_-V^{(-d_-)}(1/M_-).
\]
Indeed the Gaussian term cancels the compensation in the
shapes, and the \(V_E\) entry restores the original relative
center.

## The affine relation forced by the integer heights

Write \(g=\gcd(d_+,d_-)\), \(c_+=d_-/g\), and \(c_-=d_+/g\).
The fixed gap equation gives
\[
                 r_+=c_+r+s_+,\qquad r_-=c_-r+s_-
\]
for fixed integers \(s_\pm\) with \(d_+s_+-d_-s_-=h\) and
a single integer \(r\). There are no admissible counts if
\(g\nmid h\). Let \(x_*,y_*>0\) be the two curvatures
specified by \(T_*,\lambda_*\). Along a proposed sequence
\(r\to\infty\) since both midpoint heights diverge.
Set \(v=1/r\). Then \(1/M_+=xv/(d_+^2(c_++s_+v))\),
and similarly for the minus row. Formula (2) expresses \(V\)
as an analytic function of \(x,y,v\) near the reference and
gives the first-order term
\[
 V=v\left(\frac{x^2 A^{(d_+)}}{c_+d_+^4}
                  +\frac{y^2 A^{(-d_-)}}{c_-d_-^4}\right)+O(v^2).
\]
The invertible rational transformation from \(x,y\) to
\(T,\lambda\) and the nonsingularity of (1) give unique
analytic functions \(x(v),y(v)\) on the collision branch.
Write \(X_1=x'(0)\) and \(Y_1=y'(0)\).

The [algebraic contact result](MULTIVARIATE_DIFFUSE_PHYSICAL_EXCEPTION_POINTS.md)
makes the reference clock \(q_*=\exp[-T_* /(2nN)]\)
algebraic. Both normalized reference roots and all first
derivatives in (1) are algebraic, by differentiating the
simple root equations. Real and imaginary parts are
algebraic as well. Since the bias ratio is rational,
\(x_*,y_*\) are rational multiples of \(T_*\). Solving
the first-order implicit equations therefore gives
\[
                     X_1/T_*^2,\quad Y_1/T_*^2
                                           \in\overline{\mathbb Q}.
\]
Also \(T_*\) is transcendental by Hermite--Lindemann,
as used in the
[algebraic-strength core theorem](MULTIVARIATE_ARBITRARY_CORE_DIFFUSE_EXTENSION.md#an-angular-eliminant-for-the-exact-core).

The two midpoint heights satisfy
\[
\begin{split}
 M_+&=\frac{c_+d_+^2}{x_*}\,r+
             \frac{s_+d_+^2}{x_*}
                -\frac{c_+d_+^2 X_1}{x_*^2}+O(r^{-1}),\
 M_-&=\frac{c_-d_-^2}{y_*}\,r+
             \frac{s_-d_-^2}{y_*}
                -\frac{c_-d_-^2 Y_1}{y_*^2}+O(r^{-1}).
\end{split}
\]
Let \(\rho=c_-d_-^2 x_*/(c_+d_+^2 y_*)\in\mathbb Q_{>0}\).
Subtraction and \(c_-/c_+=d_+/d_-\) show that
\[
 L_- -\rho L_+ \ \longrightarrow           -\frac{h d_-}{y_*}+\kappa,\qquad
                        \kappa\in\overline{\mathbb Q}\cap\mathbb R. \tag{3}
\]
Indeed the count offsets give the first term, while the
scaled root-response derivatives and the half-gap shifts
give \(\kappa\). The left side belongs to a fixed rational
lattice. A convergent sequence in that lattice is eventually
constant. If \(h\ne0\), the limit in (3) is transcendental
and cannot be in the lattice. If \(h=0\) and the limit is
in the lattice, the exact affine relation
\(L_-=\rho L_++\kappa\) holds eventually. This conclusion
uses integrality of the actual heights.

## First marks exclude the forced affine family

The following argument applies at the fixed pair to any
two affine heights \(L_i=c_i' A+b_i'\) with positive rational
slopes, rational offsets and integer values at admissible
integer \(A\to\infty\). First allow an arbitrary fixed bound
on the absolute signed gap of the appended rows. If both types
are present, each count is at least a fixed positive fraction
of \(1+r_++r_-\), because \(|d_+r_+-d_-r_-|\) is bounded
and a nonzero count is at least one. When only one type is
present, its first upper factor is unique on its own.
At order \(j=n,N\), the factors in row \(i\)
are \(j L_i+t\) for \(1\le t\le j d_i\). Their normalized
affine offsets are
\[
                        (b_i'+t/j)/c_i' .
\]
If the two first offsets in the upper order are distinct,
the smaller is the unique minimum over all factors at both
orders. If they tie and the slopes differ, choose the
lower-order first factor of the larger slope. Its normalized
offset exceeds the tied upper offset by exactly
\(1/(nN\max c_i')\). The other lower first offset is
strictly larger, and any upper nonfirst offset exceeds the
tied value by at least \(1/(N\max c_i')\).
Thus the chosen factor is again unique among all factors.
The remaining case has equal slopes and equal offsets, so
the two smaller heights coincide.

For a unique factor \(F=m L_i+1\), clear the fixed rational
affine coefficients. For each competing factor \(U\),
linear elimination gives a fixed nonzero integer divisible by
\(\gcd(F,U)\). Let \(B\) be a fixed positive integer divisible
by all such integers, and by the binomial and retained
factorial data at these orders. At a prime with \(E=v_p(F)>v_p(B)\)
every competing factor has valuation at most \(v_p(B)\).
Strip the smaller-end factorials and, if necessary, reverse the
columns so the marked type increases. The marked part of the
ordinate array is zero at index zero and \(r_i E\) at every
positive index. In the other order it is identically zero.
All remaining ordinate errors are nonnegative and bounded
by \(C(1+r_++r_-)v_p(B)\), where \(C\) is fixed.
The lower convex minorants differ by at most this bound
at each index by their chord formula. Their successive
slopes therefore differ by at most twice this bound.
The marked polygon has all slopes \(r_iE/m\) before the
error, and the other polygon has only zero slopes.
The uniform count fraction then implies that a common
factor forces \(E\le C'v_p(B)\) at every prime, for a fixed
integer \(C'\). This uses the
[Newton polygon product rule](MULTIVARIATE_RETAINED_COMMON_SLOPE_BUDGET.md).
It would give \(F\mid B^{C'}\) while \(F\to\infty\).

For equal smaller heights \(L\) and unequal gaps, mark
\(F=mL+1\) at either order \(m\). It is confined to that order
since two adjacent orders have distinct first fractional
offsets and \(1/m\) cannot be \(t/m'\) in the other order.
In the marked order it occurs in both interior blocks
and in one row type at each endpoint. Its literal count
array is
\[
                     (r_-,r_-+r_+,\ldots,r_-+r_+,r_+).
\]
The minorant of its valuation contribution is the endpoint
chord of slope \((r_+-r_-)E/m\), and the opposite minorant
has only zero slopes. If the two counts differ, the fixed
signed-gap bound implies that
\(|r_+-r_-|/(1+r_++r_-)\) is bounded below by a positive
constant. Indeed, with \(d_+\ne d_-\) the identity
\[
 d_+r_+-d_-r_-=(d_+-d_-)r_-+d_+(r_+-r_-)
\]
bounds \(r_-\) by a constant times \(1+|r_+-r_-|\).
The same depth argument forces \(F\mid B^{C'}\).
If instead the counts are equal, their common value is bounded
by the signed-gap bound divided by \(|d_+-d_-|\).
The finitely many fixed count choices are covered by the
[all-directions affine-ray theorem](MULTIVARIATE_AFFINE_BASE_RAY_ALL_DIRECTIONS.md).
For rational affine coefficients, partition admissible integral
parameters into residue classes modulo a common denominator
and translate each class far enough to have nonnegative integral
intercepts. Equal smaller heights and equal gaps have a common
upper exponent and are covered by the
[common-height bounded-gap theorem](MULTIVARIATE_COMMON_HEIGHT_BOUNDED_GAP_TAIL_CLOSURE.md).
This rules out the forced affine family in every case and
completes the contact exclusion.

For the unconditional affine consequence, the
[absolute angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md)
leaves finitely many pairs of orders. The first pair is coprime
by strict positivity of the quadratic on the real axis. At each
remaining possible collision the
[signed imbalance bound](MULTIVARIATE_ALL_DEGREE_SIGNED_IMBALANCE.md)
bounds the tail's absolute signed gap in terms of the cutoff
and the fixed core. All the depth constants and the remaining
bounded-count thresholds can therefore be chosen independent
of the counts. Each selected affine factor grows to infinity
with the parameter. This gives the claimed computable all-degree
threshold.

The proof treats each specified transverse branch, including a
nonzero signed-gap twist. It leaves nonrational bias, tangential
or multiple contacts, and limits where only one curvature is
positive as separate obligations.
