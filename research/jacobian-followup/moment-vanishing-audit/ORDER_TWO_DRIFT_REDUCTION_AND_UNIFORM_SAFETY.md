# Order-two drift reduction and uniform safety

15 November 2026. This extends the [homogeneous quadratic theorem](QUADRATIC_OPERATORS_AFFINE_TERMS_AND_CLOSED_SAFETY.md)
to every constant-coefficient operator of order at most two. It is an internal
written argument, with independent review outstanding and no novelty claim.
The [registry](../../MATH_STATUS.json) is authoritative. Earlier proof notes and
their narrower hypotheses are preserved unchanged.

## Scope and previous evidence

Over a characteristic-zero field, write
\[
 \Lambda=\Delta_A+D_b+c,\quad
 \Delta_A=\sum_{i,j}A_{ij}\partial_i\partial_j,\quad
 D_b=\sum_i b_i\partial_i,
\]
where A is symmetric and all coefficients are constant in the polynomial
variables. The coefficients A,b,c may vary in a family. Phases P and
multipliers Q are arbitrary polynomials. Purity means polynomial identities
\(\Lambda^mP^m=0\), not values at the origin. Safety means full eventual
mixed vanishing, subject to purity.

The [earlier binary drift theorem](../../archive/non-elliptic/extended-geometry/SEPARABLE_GVC_ESCAPE_OBSTRUCTIONS.md)
already closes the two-variable case, including much higher operator orders;
the extension here is in every dimension. [Zhao0704.1691v2](https://arxiv.org/html/0704.1691v2)
explicitly treats homogeneous second-order operators; its rank reductions do
not silently remove a first-order part. The retained arXiv version is
28 October 2007, despite an August 2026 date in the HTML body. The
[commutator lemma](https://arxiv.org/html/1206.2836v2) and previously imported
Gabber/Wang theorems remain the classical inputs. No new external theorem is
needed for the drift reduction below.

## 1. Complete pure criterion and a degree-free finite prefix

If c is nonzero, \(\Lambda\) is invertible on every bounded-degree polynomial
space: its derivative part lowers degree, so a finite geometric inverse
exists. Thus the first pure equation forces P=0. This branch is trivially safe.

Now suppose c=0 and put \(J=A\operatorname{Hess}P\). We prove
\[
 \boxed{\ \Lambda^mP^m=0\ (1\le m\le n+1)
 \iff J^n=0\ \text{and}\ D_bP=0
 \iff\Lambda^mP^m=0\ (m\ge1).\ }                              \tag{1}
\]
The prefix n+1 is sufficient; no optimality is asserted.

Use the scalar Hamilton--Jacobi potential from the preceding theorem:
\[
 F_t=x-tA\nabla P,\quad G_t=F_t^{-1},\quad
 U_t=P(G_t)-\frac t2\nabla P(G_t)^TA\nabla P(G_t).
\]
The t-adic identities \(\nabla U_t=\nabla P(G_t)\) and
\(\partial_tU_t=\nabla U_t^TA\nabla U_t/2\) hold without any purity assumption.
For
\[
 H=\exp\!\left(\frac{t\Lambda}{2s}\right)e^{sP},\qquad L=\log H,
\]
the normalized t-adic Laurent interpretation is the same as before. Its PDE is
\[
 \partial_tL=\frac{\Delta_A L+D_bL+\nabla L^TA\nabla L}{2s}.       \tag{2}
\]
Each coefficient of L still has s degree at most one; the new drift term
divided by s has degree at most zero. The coefficient of s is still U_t.
The first n+1 pure equations imply
\[
 L=sU_t\pmod{t^{n+2}},\qquad \Lambda U_t=0\pmod{t^{n+1}}.       \tag{3}
\]
Here \(\Lambda(1)=0\) is essential; the c!=0 branch was handled separately.
The exact trace identity becomes
\[
 (\Lambda U_t)(F_t(x))
  =\operatorname{Tr}\big(J(I-tJ)^{-1}\big)+D_bP(x).              \tag{4}
\]
Thus (3) forces
\[
 \operatorname{Tr}(J^k)=0\quad(2\le k\le n+1),\qquad
 \operatorname{Tr}J+D_bP=0.                                    \tag{5}
\]

The missing first trace cannot be assumed zero. Instead work over the
algebraic closure of K(x). If J had distinct nonzero eigenvalues
\(\lambda_1,\ldots,\lambda_r\), with positive integer multiplicities
\(m_1,\ldots,m_r\), the equations for k=2,...,r+1 would give
\[
 \sum_i m_i\lambda_i^{j+2}=0\quad(0\le j<r).
\]
The matrix is an invertible Vandermonde matrix times
\(\operatorname{diag}(\lambda_i^2)\). It would force every m_i=0, impossible
in characteristic zero. Hence J is nilpotent. Now its first trace is zero,
and (5) gives \(D_bP=0\). This proves the difficult implication of (1).

Conversely, when J is nilpotent, P is pure for \(\Delta_A\) by the preceding
theorem. If \(D_bP=0\), that derivative kills every P power and commutes with
\(\Delta_A\), so \(\Lambda^mP^k=\Delta_A^mP^k\) for every k,m. This proves
the converse. For arbitrary c the full criterion is the closed union
\[
 \{P=0\}\ \cup\ \{c=0,\ J^n=0,\ D_bP=0\}.                    \tag{6}
\]
It is still detected by the first n+1 pure polynomial equations.

## 2. The drift preserves qualitative safety, not multiplier values

On the nontrivial pure branch, equation (4) gives
\(\Lambda U_t=0\). The heat identity and positive s coefficients are
\[
 H=e^{sU_t},\qquad
 U_t^k=k!\sum_{j\ge0}
 \frac{t^j\Lambda^jP^{j+k}}{2^j j!(j+k)!}.                       \tag{7}
\]
Therefore the following are equivalent:

1. Full mixed safety for Lambda.
2. Eventual vanishing for the single multiplier Q=P.
3. U_t is polynomial in x,t.
4. F_t has a polynomial inverse over K[t].
5. Full mixed safety for the homogeneous quadratic part \(\Delta_A\).

Equivalence with the last item uses the same scalar U_t; for P powers the two
operators already agree. The commutator lemma gives the explicit bound
\[
 \Lambda^m(QP^m)=0\quad\text{if }m>(1+\deg_tU_t)\deg Q.          \tag{8}
\]
Constant Q follows from purity. The P=0,c!=0 branch is handled directly and
does not use the heat identity (7).

It would be false to remove D_b from marked contractions or their optimal
cutoffs. For example,
\[
 \Lambda=\partial_x^2+\partial_y,\quad P=x,\quad Q=y^d
\]
has all pure moments zero, but
\[
 \Lambda^{2d}(y^dx^{2d})=\binom{2d}{d}d!(2d)!
   =\frac{((2d)!)^2}{d!}\ne0.                                 \tag{9}
\]
Only the term with d copies of each derivative survives. Giving x weight1
and y weight2 shows that every multiplier of total degree d is killed for
m>2d. Thus slope2 is sharp here, whereas \(\partial_x^2\) alone has sharp
slope1 for the same phase. In particular the previous D<=1 homogeneous
quadratic bound cannot simply be copied to an operator with drift.

## 3. An exact quotient and a dimension-threshold relation

Suppose b!=0. Choose linear coordinates (z,y) with \(D_b=\partial_y\).
By (1), P=P(z) has n-1 variables. Write A in these coordinates in blocks,
and let A_0 be its z,z block. The Hamilton--Jacobi equation of P stays
independent of y, with exactly the potential for \(\Delta_{A_0}\) in n-1
variables. Cross derivatives involving y and the y,y block do not enter this
scalar equation. Consequently
\[
 (\Lambda,P)\text{ is safe}
 \iff (\Delta_{A_0},P)\text{ is safe in }n-1\text{ variables}.    \tag{10}
\]
This is a qualitative quotient for every polynomial multiplier; it does not
assert equality of arbitrary marked contractions. If the pair fails, Q=P
already detects failure and belongs to the quotient variables.

Let d_0 be the least failing dimension for homogeneous quadratic operators,
and d_1 the least failing dimension for order-at-most-two operators with a
nonzero first-order part. Allow infinity if no failure exists. Then
\[
 \boxed{d_1=d_0+1.}                                           \tag{11}
\]
One direction is (10). For the other, append an unused variable y to a
quadratic failure and add \(\partial_y\) to the operator; its original phase
and detecting multiplier remain independent of y. The same statement holds
under any phase degree bound, or for homogeneous phases of a prescribed
degree, since linear changes and unused-variable padding preserve these
conditions. This relates thresholds; it does not supply their unknown numbers
or certify the previously reported quartic witnesses.

## 4. Closed safety and uniformity for every operator of order at most two

Let S(A,P) be the closed pure-safe field locus from the preceding theorem.
Then the full pure-safe field locus for all A,b,c,P is
\[
 \{P=0\}\ \cup\ \{c=0,\ D_bP=0,\ (A,P)\in S\}.               \tag{12}
\]
Both branches are closed, including where b becomes zero, A changes rank or
c becomes nonzero. The explicit finite equations for S are retained inverse
composition equations; no general integration or period algorithm is needed.

For fixed n>=1,D>=2 every safe field pair has the uniform cutoff
\[
 m>C(n,D)\deg Q,\qquad
 C(n,D)=2n(D-1)^{n-1}-n+1.                                    \tag{13}
\]
For active quadratic rank at least two, the preceding proof actually bounded
\(1+\deg_tU_t\) by this C. At rank one, a pure phase is affine in the sole
active coordinate, and the scalar potential has t degree at most one; thus
\(1+\deg_tU_t\le2\le C(n,D)\). At rank zero, U=P has t degree zero.
Equation (8) therefore proves (13) also with drift. For degree D=1, slope2
suffices, and (9) proves that it is needed in general. For D=0, slope1
suffices. These are nonoptimal bounds in many special cases.

Thus all polynomial coefficient-family pullbacks have closed safe field loci,
and the increasing closed field cutoff chain stabilizes at the stated bound.
This closes the field uniformity question for **all** constant-coefficient
operators of order<=2. It is not a full ordinary-ideal result over nonreduced
coefficient algebras, nor a Gaussian origin-moment theorem. The universal
inverse-equation producer remains unimplemented.

## 5. Cubic operator-order minimum and remaining boundaries

The preceding Wang-based cubic theorem, combined with (1)--(8), shows that
every pure phase of degree<=3 is safe for every constant-coefficient operator
of order<=2, in every dimension, including nonhomogeneous operators.
The independently retained DVG1 example has a cubic phase and an order-three
operator. Hence the least operator order permitting a cubic GVC failure,
over all dimensions, is **exactly three**. No new Dvorsky replay is needed.

Any obstruction-preserving lift of that cubic to an operator of order<=2
must raise phase degree above three. Allowing a first-order or constant term
does not evade the obstruction. Higher-degree lifts, genuinely higher-order
operators, general Gaussian uniformity and ordinary nilpotent geometry remain
open; the protected DVG4 one-pair obstruction and ternary Delta-power frontier
are unchanged. In particular a zero-constant-term assumption for higher-order
operators does not imply their drift annihilates P: the s-degree-one heat-log
argument here depends on order at most two.

The [checker](verify_order_two_drift.py), [receipt](order_two_drift_verification.json)
and [source audit](order_two_drift_source_audit.json) retain finite controls
for the trace identity, a misleading first-moment cancellation, the quotient
with cross terms, the constant-symbol branch and the sharp drift cutoff.
These controls do not establish the universal trace/Vandermonde proof or
replace independent mathematical review.
