# A 34-cube lower bound for the fixed padded companion

Let K=X+Q+C be the scale11=2 twelve-variable base map used in
[KCUB5](POWER89.md), with C supported in its first five outputs. In
Z=(X,W,tau), put

\[
P(Z)=(\tau Q-\tau^2\iota W,\ C_1,\ldots,C_5,\ \tau^3).
\]

Every representation of this particular vector of eighteen cubics as

\[
P(Z)=\sum_{a=1}^r b_a\ell_a(Z)^3
\]

requires **r at least 34**, even over a characteristic-zero extension
field. Together with KCUB5 this gives **34 <= r <= 89** for this fixed
tensor. This is not a minimum-dimension bound for all Keller maps.

For every output i and variable j,

\[
\partial_jP_i=3\sum_a (b_a)_i(\partial_j\ell_a)\ell_a^2.
\]

All first derivatives therefore lie in the span of r squares, so their
linear span has dimension at most r. The exact derivative matrix for P
has rank 34. It has 324 rows, 77 nonzero rows and 51 occurring quadratic
monomial columns. Rational rank is preserved by characteristic-zero
field extension. Invertible linear source or output changes also
preserve this derivative-space bound.

[The primary computation](compute_fixed_tensor_bound.py) reconstructs P
from the audited source with SymPy and computes the rank using FLINT.
[The independent checker](verify_fixed_tensor_bound.py) reconstructs it
from the separately retained scale11=1 coefficient artifact by diagonal
conjugation, uses only standard-library Fraction arithmetic, and records
34 independent derivative labels and pivot monomials. Both paths agree
on a canonical hash of all tensor coefficients. See
[primary receipt](fixed-tensor-bound-primary.json) and
[independent receipt](fixed-tensor-bound-validation.json).

The initial scratch diagnostic counted one zero-polynomial monomial as
an extra ambient column; its rank 34 was correct. The reproducible
primary and independent receipts exclude that zero column.

The one-step Koszul flattening of the derivative system has modular
rank 562 at prime 1000003, giving ceil(562/17)=34. This primary-only
[diagnostic](fixed-tensor-koszul.json) does not improve the checked bound.

The warm start in the 178-form pool made 16,980 exact exchanges in
60 seconds without improving 89. Enlarging to 858 primitive weighted
factor sums (weights from -2 through 2) left cube-span rank 121 and made
1,966 exchanges in 90 seconds without improvement. Both runs are
bounded heuristic searches, not proofs that 89 is optimal. Their
[178-form checkpoint](power-basis-from89/checkpoint.json) and
[858-form checkpoint](power-basis-weighted/checkpoint.json) pin the
sources, seed and inputs.

A next construction test should mix factors across the existing
term groups or change the companion tensor, rather than repeat either
bounded run. Any smaller candidate still requires a complete integer
matrix, determinant argument and independently verified collision.


The [36-term base](PROOF12.md) does not define a new companion orbit.
If T is the diagonal change from the scale11=1 source to scale11=2,
and U is the source matrix displayed in that proof, then R=T^{-1}U
has only these nonidentity source coordinates:

\[
(RX)_9=(x_9+3x_8)/2,\qquad (RX)_{11}=x_{11}+3x_{10}.
\]

It fixes the first five coordinates. With V=diag(R,I_5,1), the associated
padded tensor is V^{-1}P(VZ). Thus its minimum shared-cube count is the
same. A factorization P=B(LZ)^{circ3} transports to factors
B'=V^{-1}B and L'=LV, with L'B'=LB; the same identity holds after
zeroing the pure tau row. This explains why transporting the old
factorization does not itself give a smaller cubic-linear map.
[The exact bridge](power36_bridge.py) checks every tensor coefficient
and the unchanged matrix product; its [receipt](power36-bridge.json)
is primary evidence, not a new independent theorem.
A separate50-restart termwise test gave91 forms on this base, as
recorded in its [candidate](power36-termwise/candidate.json). It does not
supersede the transported89-form upper bound. Fewer displayed monomials
need not produce fewer cubes under a particular polarization algorithm.
