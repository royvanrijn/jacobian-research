# The mixed cubic Jacobian pencil and the homogenizing chain

This note identifies an exact compiler objective and an all-affine
obstruction for the retained eighteen-variable index17 construction.
The obstruction concerns this same mixed compiler, not every possible
cubic reduction of a Keller map.

## Companion resolvent

Let K(X)=X+Q(X)+C(X) be a tangent-normalized cubic Keller map in n variables.
Factor its homogeneous cubic part as C=iota C_*, where iota is a constant
n-by-m matrix of rankm. The economical companion is
\[
 H(X,W,\tau)=(\tau Q(X)-\tau^2\iota W,\ C_*(X),\ 0).
\]
At tau=1, its old-variable Jacobian block is
\[
 A=\begin{pmatrix}DQ&-\iota\\ DC_*&0\end{pmatrix}.
\]
Put M(lambda)=I+lambda DQ+lambda^2 iota DC_*=DK(lambda X).
Since detDK=1, S(lambda)=M(lambda)^(-1) is a polynomial matrix.
**Proposition CompanionResolventCriterion.** Block elimination gives the exact identity
\[
 (I+\lambda A)^{-1}=
 \begin{pmatrix}
 S&\lambda S\iota\\
 -\lambda DC_*S&I-\lambda^2DC_*S\iota
 \end{pmatrix}.                                                     \tag{P1}
\]
Its maximum lambda degree is exactly one less than the nilpotency index
of A. This is an objective on the actual differential circuit, rather
than the number of its auxiliary variables.

For general tau, over K(tau),
\[
 A_\tau=
 \begin{pmatrix}\tau DQ&-\tau^2\iota\\DC_*&0\end{pmatrix}
 =\tau\,D^{-1}AD,\qquad D=\operatorname{diag}(I,\tau I).
\]
With b_tau=(Q-2tau iota W,0), the full Jacobian has powers
\[
 JH^k=
 \begin{pmatrix}A_\tau^k&A_\tau^{k-1}b_\tau\\0&0\end{pmatrix}.         \tag{P2}
\]
The identities extend polynomially to tau=0. To reduce the full index,
one must shorten both the resolvent chain and the extra tau-column
chain; reducing a dimension count alone need not do either.

## Exact retained profile

For the retained twelve-variable map at scale11=2, n=12 and m=5.
Literal multiplication gives
\[
 \deg_\lambda S=15,\quad
 \deg_\lambda(S\iota)=14,\quad
 \deg_\lambda(DC_*S)=14,\quad
 \deg_\lambda(DC_*S\iota)=13.
\]
Thus A has exact nilpotency index16. In the full homogeneous map JH^16
has five nonzero entries, all in the tau column, and JH^17=0.

Write N(X)=DK(X)^(-1) and N_15 for its degree15 homogeneous part.
The only nonzero column of N_15 is columnx_8 (zero-based7); also
N_15 iota=0. In rowx_4,
\[
 (N_{15}Q)_{x_4}=-288x_1^{10}x_2^4x_3^3\ne0 .                       \tag{P3}
\]
Here Q_(x8)=x_2^2. Thus the longest homogeneous chain is seeded by that
quadratic carrier. The companion's seventeenth nilpotency step is the
extra tau derivative, not a nonzero sixteenth power of A.

The entire A^15 has ten nonzero entries, only in columnsx_8 andw_5.
The two carrier rows themselves are
dx_8 A=2x_2 dx_2 and dw_5 A=3x_2^2 dx_2. Their exact cancellation
annihilates A^16. This provides a concrete shared-carrier target for
a genuinely different source/target compiler.

## An all-affine normalization obstruction

The preceding obstruction survives every affine source recentering and
tangent normalization of this same mixed cubic compiler.

In general let N=DK^(-1) have maximum degree d and nonzero highest part N_d.
The highest homogeneous degree in N DK=I gives
\[
 N_d DC=0.                                                         \tag{P4}
\]
At any source point a, define
\[
 K_a(X)=DK(a)^{-1}\bigl(K(a+X)-K(a)\bigr).
\]
Its quadratic and cubic parts are
\[
 Q_a=DK(a)^{-1}(Q+DC\,a),\qquad C_a=DK(a)^{-1}C .
\]
Since DK_a(X)^(-1)=N(a+X)DK(a), its top homogeneous inverse-Jacobian part is
N_(a,d)=N_d DK(a). Consequently
\[
 N_{a,d}Q_a=N_dQ,\qquad
 N_{a,d}\iota_a=N_d\iota,\quad \iota_a=DK(a)^{-1}\iota .             \tag{P5}
\]
The equality for Q uses(P4). Additional invertible linear source changes
transport the first expression by
N_dQ -> U^(-1)(N_dQ)(UX), so its nonvanishing persists. An independent
affine target change disappears upon tangent normalization. A change of
basis for the cubic output space is a constant fibre conjugacy.

**Theorem AffineCompilerIndexObstruction.** Suppose N_dQ is nonzero. In the homogeneous companion, the top portion
of A^d b at W=0 is (-1)^d N_dQ. Therefore its full Jacobian nilpotency
index is at least d+2. Formula(P5) proves the same lower bound for every
affinely recentered/tangent-normalized version of the mixed compiler.

**Corollary RetainedCompilerIndexLowerBoundSeventeen.** For the actual data d=15 and(P3) is nonzero. Hence no affine choice of
source basepoint, independent affine target coordinates, or constant
cubic-output basis can improve the homogeneous index below17.

This remains true after every constant cubic target shear fixing tau,
H -> H+(c tau^3,0), with c in K^(n+m). It leaves A unchanged and adds
3c tau^2 to b. At tau=1, the added top A^d c terms have X-degree at mostd,
and the bottom-to-top resolvent block has degree at mostd-1. They cannot
cancel the nonzero degree d+2 term N_dQ. In the actual rowx4 one can
equivalently distinguish tau-degrees15 versus17or18 before specialization.

This is a precise obstruction to affine normalization and that specified
target-shear family. It does not exclude nonlinear stable equivalences,
a different factor-aware cubic reduction, or a shorter failing index.

## A legal product-splitting move

**Lemma LegalProductSplitting.** For a product fg in outputi, append independent variablesU,V. The source
shear(U,V)->(U-f,V-g), followed by the target shear
outputi -> outputi-outputU*outputV, gives
\[
 K'_i=K_i-fg-UV+Ug+Vf,\quad K'_U=U-f,\quad K'_V=V-g.                 \tag{P6}
\]
Both are genuine polynomial shears of determinant1. A collision lifts
by U=f(p),V=g(p). Replacing a monomial product merely on a constraint
graph would not establish the same conclusion.

All displayed block identities are proved directly. The construction
context is the retained [POWER89 companion argument](../cubic-reduction/POWER89.md)
and de Bondt's
economic stabilization arXiv2610.24889v1 Sections4.6--4.7. No priority
claim about companion linearization, restriction of scalars or degree
reduction is made.


## Exact reproduction

Run from the repository root:

    python3 research/jacobian/subtree-shuffling/verify_compiler_jacobian_pencil.py \
      --repo-root . --output /tmp/compiler_jacobian_pencil_review.json

The accompanying input selects the pinned canonical45-term seed, the
exact degree profile and three affine-center controls. The standard-library
checker reconstructs the entire resolvent, verifies both inverse orders,
the full tau-block identities and top-chain support, the invariant
N_d DC=0 and N_d Q witness, and complete affine-normalization identities.
It also checks the universal product-splitting identity and negative
controls for its matrix order and resolvent signs. The universal affine
and target-shear conclusions follow from the written degree proof, not
from sampling the three centers.

The mandatory ledger pins this proof, its checker/input, and the exact
retained seed, arithmetic module and companion/index17 proofs. No new
external, formal-verification or priority assurance is claimed.

Updated 20 November 2026.
