# Anisotropic rational monomial congruences by triangular descent

Let \(d\ge2\), let \(C_i\ge1\) for \(1\le i\le d\), and
write \(P=\prod_i C_i\). There is an absolute constant
\(B_d\) depending only on \(d\) such that, if \(k\ge B_d P\)
and \(m\ge2\) is \(k\)-full, the number of distinct exact
rational outputs
\[
 \Phi(\mathbf x)=\prod_i x_i^{a_i}\equiv K\pmod m,\qquad
 H(x_i)\le m^{C_i/k},
 \tag{ATC1}
\]
from positive rational \(m\)-unit coordinates is
\[
 O_{d,\mathbf C,k,\mathbf a,\varepsilon}
              \left(m^{d(d-1)/(2k)+\varepsilon}\right).     \tag{ATC2}
\]
The \(a_i\) are fixed nonzero signed integers, and the
bound is uniform in the unit target \(K\).
The fullness threshold is independent of the orders
\(a_i\). Chart-complexity constants may depend on them. The count is for exact outputs, so choose
one original representative per output.

This strengthens the fixed-height-constant threshold
in (FC1). The anisotropic threshold uses rectangular
monomial supports. A derivative of a rank-defect
polynomial yields a triangular next equation directly,
so no coefficient-height bound from general
elimination or projective covering is needed.
[Signed monomial congruences (FC1)](ABC_OMITTED_COORDINATE_GRAPH.md#signed-monomial-congruences-over-full-moduli)

## Clauses retained along the descent

At a stage choose \(j\) remaining original coordinates
indexed by \(I\), with \(s=d-j\) already removed
coordinates in their order of removal, \(i_1,\ldots,i_s\).
There are equations \(G_1=\cdots=G_s=0\) on the
original sample in this branch. Equation \(G_t\)
depends only on the original coordinates remaining
just before \(i_t\) was removed, and
\(\partial_{i_t}G_t\ne0\) on the branch.
The gradient submatrix on the removed coordinates
is upper triangular with nonzero diagonal.
In addition each equation has bounded degree at
fixed parameters, weighted degree
\(\sum_i C_i\deg_{x_i}(G_t)=O_d(P)\), and
coefficient height \(m^{O_d(P/k)+o(1)}\).
The initial stage \(j=d\) has no equations.

On the open unit locus with this triangular
minor nonzero, the preceding equations define
a smooth variety of dimension \(j\) over
characteristic zero. Append the logarithmic
monomial row
\[
                   (a_i\prod_{h\ne i}x_h)_{i=1}^d
\]
to their gradients. The locus where all the
enlarged maximal minors vanish has finitely
many geometric irreducible components, of
bounded degree at the fixed parameters.
On each its exact monomial value is constant:
the logarithmic differential annihilates
the tangent space to the smooth ambient
locus, hence annihilates the tangent space
at a generic smooth point of the critical
component. Characteristic zero identifies a
rational function with vanishing differential
as a constant over the algebraic closure.
Taking at most one chosen representative per
exact output, this contributes boundedly many
points. Thus no algebraic dimension recursion
on a critical component is needed.

For the other points choose an enlarged minor
\(J\ne0\), branching among finitely many choices.
The evaluated numerator of \(J\) has height
\[
                    m^{O_d(P/k)+o(1)}.
\]
Indeed each equation has the weighted degree
and coefficient-height caps stated above.
Coordinate denominator clearing adds only
\(O_d(P/k)\) to the exponent in a fixed
dimensional determinant. The integers \(a_i\)
are fixed with respect to \(m\), so their sizes
occur only in the fixed multiplicative constants.
All the clearing denominators are units at
the primes dividing \(m\).

For \(p^e\Vert m\), mark \(p\) bad if
\(v_p(J)>e/5\). For sufficiently large \(B_d\)
the combined bad log weight is at most
\((1/10)\log m\), increasing the lower
height threshold in the fixed parameters
if necessary. At a good prime write the
polynomial congruence as
\[
       F_0=\prod_i x_i^{a_i^+}
                    -K\prod_i x_i^{a_i^-}=0\pmod {p^e}.
\]
The gradient at the original unit point
equals a unit multiple of the displayed
logarithmic row modulo \(p^e\). With all but
the minor's coordinates fixed, multivariable
Hensel lifts to \(G_1=\cdots=G_s=F_0=0\)
at precision at least \(4e/5\), with the
minor remaining nonzero. The reference
set restricted to its full rank locus
has local dimension \(j-1\) in the unit box.

## A rectangular determinant on the current coordinates

When \(j\ge2\), take
\[
   P_I=\prod_{i\in I} C_i,\qquad
   D_i=\lceil A_d P_I/C_i\rceil\quad(i\in I)
\]
for sufficiently large \(A_d\). The rectangular
support has size \(N=\prod_{i\in I}(D_i+1)\),
at least \(A_d^jP_I^{j-1}\).
Set \(r_0=\lfloor N/2\rfloor\).
Clear evaluation rows by the product of the
current coordinate denominators to their
respective powers \(D_i\). Each row is integral
and has height at most \(m^{H_I/k}\), where
\[
                  H_I=\sum_{i\in I} D_i C_i=O_d(A_d P_I).
\]

At the good primes apply the uniform
\(T_r\) charts of (FC1) to the exact smooth
reference sets. Their defining equations
and degrees belong to a bounded-complexity
family for fixed parameters. Large primes
use uniform parametrization, and the
finitely many smaller primes use the
fixed-field parametrization. Divide the
chart parameters into classes modulo
\(p^{b_p}\) with \(b_p=\lfloor e/k\rfloor\).
The choices of chart, parameter class and
bad-prime flags cost \(O(m^{(j-1)/k+\varepsilon})\).
Write \(S_{r_0}\) for the sum of the first
\(r_0\) orders in the \(j-1\) local Taylor
parameters, in nondecreasing order with
multiplicity. Then
\[
        S_{r_0}/r_0\asymp_j r_0^{1/(j-1)}
                   \gg_d A_d^{j/(j-1)} P_I .
\]
Take \(A_d\) large so that \(S_{r_0}>8r_0 H_I\).
The last Taylor order is \(O_{d,A_d}(P_I)\).
Choose the approximation order \(r\) just
above it. Increasing \(B_d\) if needed, use
\(k\ge B_d P\) and \(P_I\le P\) to make
\(r b_p<4e/5\).

Each \(r_0\)-minor of a residue class is
divisible by \(p^{b_p S_{r_0}}\) at good
primes by Taylor expansion. Repeated
Taylor columns give zero, and lift and
truncation errors have order at least
\(r b_p\). Row multipliers are units in
the local determinant. Since the good
primes carry most of the total log
weight and \(b_p\ge e/(2k)\), the minor
is divisible by an integer at least
\(m^{S_{r_0}/(4k)}\).
Its absolute value is
\(O_N(m^{r_0 H_I/k})\) and so it vanishes.
A rank-half evaluation matrix yields
a nonzero polynomial \(F\) on the
current coordinate support, of
coefficient height \(O_N(m^{H_I/k})\).
For this short null vector, map an
integer coefficient box of side \(T\)
to the evaluations on independent
rows. There are \(O_N((T m^{H_I/k})^{N/2})\)
possible images, fewer than the domain size
when \(T\) is a suitable fixed multiple
of \(m^{H_I/k}\).

At each original point of \(F=0\)
choose a multi-index \(\gamma\) of
minimal order with
\(\partial^\gamma F\ne0\) at the point.
Its order is at least one, since
\(F\) vanishes, and at most the degree.
Choose \(i\in I\) with \(\gamma_i>0\).
The polynomial
\[
                    G_{s+1}=\partial^{\gamma-e_i}F
\]
vanishes at the point and has
\(\partial_i G_{s+1}\ne0\).
Branch among the boundedly many
possibilities for \((\gamma,i)\),
adjoin this equation, and remove
\(i\) from the current coordinates.
It depends only on the preceding
current coordinates and maintains
the weighted degree and coefficient
height bounds, including after
taking the derivative. This gives
exactly the stated triangular
conditions for the next stage,
even on singular points of the
undifferentiated auxiliary locus.

## Terminal reconstruction and cost

When \(j=1\) the corrected good-prime
reference sets consist of isolated
points, with boundedly many choices
by Bézout. The lift choices and
bad-prime indicators cost \(m^{o(1)}\).
Two original tuples in the same
record then agree at a combined
modulus exceeding \(m^{1/2}\), since
the retained log weight is at least
\(9/10\) and the correction precision
is at least \(4e/5\).
Each rational cross difference
has size at most \(2m^{2\max_i C_i/k}\),
which is less than \(m^{1/2}\) by
the fullness threshold. Each record
therefore contains at most one
chosen representative. Summing the
chart costs over \(j=d,\ldots,2\)
proves (ATC2).

If product amplification by a fixed \(\ell\)
is desired, choose prime sign colors
of the exact rational outputs before
multiplying. This retains a subpower
fraction and ensures subpower presentation
multiplicity, precisely as in the
two-basis argument. The resulting
coordinate heights are raised to
the \(\ell\)-th power. For prime-resolved
weighted valuation digits on two free letters,
of width \(R_i\) in coordinate \(i\), the coordinate
height is at most \(X^{O(L R_i)}\) with an absolute
constant: bound each prime's digit by \(R_i\) and
sum the total retained radical mass, which is \(O(L)\)
independently of \(d\). Use a prescribed full modulus of
numerical exponent \(w\ge w_0>0\),
with \(k=\lfloor\kappa/L\rfloor\).
One may take \(\ell\asymp_{d,\kappa}1/\eta\)
for normalized chart cost at most
\(\eta\). The condition becomes
\[
             L\prod_{i=1}^d\max(1,R_i)
                    \ll_{d,w_0,\kappa} \eta^d .
\]
Thus the direct release at this chart
margin needs a small anisotropic
digit-width product, independently
of any compression to two bases.


## An exponential dimension constant

The fullness threshold can be taken \(B_d\le A^d\)
for an absolute \(A\). At a stage of dimension \(j\),
choose \(A_j=a^{j-1}\) for a sufficiently large fixed
\(a\), in place of one \(A_d\) in the rectangular
support construction. Then
\(r_0^{1/(j-1)}\ge c A_j^{j/(j-1)}P_I\)
and the Taylor mean order is
\(\gg (j-1)r_0^{1/(j-1)}\) with an absolute constant.
Indeed
\(\binom{u+j-1}{j-1}
 \le(e(u+j-1)/(j-1))^{j-1}\)
shows that fewer than half the required Taylor monomials
have order below \(c'(j-1)r_0^{1/(j-1)}\), for
a sufficiently small absolute \(c'\) and
\(r_0^{1/(j-1)}\) sufficiently large.
The ratio \(A_j^{1/(j-1)}=a\) therefore supplies the
needed fixed margin against
\(H_I\le j(A_j+1)P_I\).
Conversely the last Taylor order is at most
\(O(j)(A_j+2)^{j/(j-1)}P_I\), by the opposite
binomial bound. This is at most \(A^{d}P\) after
increasing the absolute base \(A\).
The row coefficient height exponent and weighted
degree are \(O(j)(A_j+1)P_I/k\) and
\(O(j)(A_j+1)P_I\), respectively.
Summing the denominator and coefficient budgets
over at most \(d\) rows of each Jacobian minor
costs \(O(d^2)a^dP/k\), also bounded by
\(A^dP/k\) for a larger absolute base. These are all
the requirements on \(k/P\); the coefficient norms,
chart complexity and degree counts enter only fixed
multiplicative constants at each budget.

Consequently the count also handles dimensions
varying in the outer small-budget limit. In the actual
full-modulus application with fixed
\(w_0,\kappa,\eta>0\), take
\(\ell\asymp d^2/(\kappa\eta)\) for the chart allowance.
If the retained digit box satisfies
\(\prod_i\max(1,R_i)\le L^{-1+\delta}\) for fixed
\(\delta>0\), one has \(A^d\prod_i C_i\le k\)
for all sufficiently small budgets whenever
\[
             d\le c_\delta
                    \frac{\log(1/L)}{\log\log(1/L)}
\]
and \(c_\delta>0\) is sufficiently small.
Indeed \(k\asymp\kappa/L\), and the left side is at
most \((C d^2)^d L^{-1+\delta}\), where \(C\) depends
on fixed parameters. Fix the displayed \(d,\ell\),
all bases, and all diagram choices at each individual
budget before taking the height limit. This permits
a larger family of sparse mixed grids than a theorem
requiring a single rank bound throughout the small-budget
limit, while the chart exponent remains \(\eta L\).
