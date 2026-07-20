# Triangular-chart attack: exact scope and current results

## Problem encoded

For

\[
\Phi_p(x,y)=(u,v)=(y+p(1/x),1/x)
\]

the implementation verifies `det(D Phi_p)=v^2`. It searches Laurent
coordinates with bracket `v^-2` and imposes exact polynomialization after
pullback. The planted pair

\[
A=1/v,\qquad B=u-p(v)
\]

always pulls back to `(x,y)` and has the required bracket. Consequently, a
blanket statement that one exceptional divisor cannot polynomialize two
algebraically independent reciprocal-Jacobian coordinates is false. Any valid
obstruction must use noninjectivity, generic degree, or stronger boundary data.

All `p` charts are related by the triangular target automorphism
`(u,v)->(u+p(v),v)`. They share one canonical geometric signature, although a
finite Laurent support window changes under this automorphism and is therefore
retained as separate computational data.

## Complete pole phase in the requested box

The box `0<=i<=6`, `-8<=j<=8` has 119 Laurent monomials. Exact cancellation was
exhausted for normalized sparse `p` strata through degree four. Kernel
dimensions are:

| `deg p` | kernel dimension |
|---:|---:|
| 0 or 1 | 63 |
| 2 | 59 |
| 3 | 49 |
| 4 | 39 |

For every stratum the first forbidden valuation before solving is the `x^-1`
pole. Ranks over `Q` agree with ranks modulo `1000003`, `1000033`, and
`1000037`; no bad prime occurred.

Sparse staircase runs and independent bound increases (`i<=7`, `j>=-9`, and
`j<=9`) are stored as generated JSON. Increasing the negative Laurent range
adds polynomializable directions; increasing the positive range mostly adds
new pole equations. Higher-degree `p` sharply reduces the kernel inside a
fixed box because translating `u` spreads support toward higher positive
`v`-powers.

## Nonlinear Stage A milestone

After exact polynomialization, the collision-normalized Keller systems of
ordinary polynomial degrees 1, 2, and 3 were solved completely with msolve F4.
For each degree, over `Q` and over all three large primes, the reduced Gröbner
basis contains `1`. These are machine-checkable characteristic-zero emptiness
certificates for the named degree/support families.

This is not a new low-degree theorem—the corresponding cases are classically
known. Its purpose is validation: normalization, F4 export, modular comparison,
and exact emptiness certification agree end to end.

## Next finite milestone

1. Enumerate sparse kernel subspaces containing 10–30 Laurent monomials, modulo
   the triangular chart equivalence and source/target swaps.
2. Saturate by required leading coefficients before F4, removing coordinate
   and degree-drop components.
3. Classify each empty ideal by its earliest pole/Newton row and cluster its
   Gröbner certificate signature.
4. Run 30–100-term one-divisor supports only after stable obstruction classes
   appear.
5. If collision-normalized survivors remain, compute generic degree,
   collisions, line restrictions, and Puiseux/dicritical data before homotopy.

## Continuation results

Coefficient-rank stratification was exhaustively sampled over every normalized
lower-coefficient tuple in `F_3` and `F_5`, and randomly at the three large
primes. For every degree 1 through 4, rank was constant across the entire small
finite field and every large-prime sample. This strongly supports—but does not
yet prove over `Q`—that the rectangular pole rank depends only on `deg p`.

Three two-divisor unimodular toric charts were then searched with 10- and
15-term cones. All six collision-normalized systems have exact reduced basis
`[1]` over `Q`, with the same outcome over the three primes. Thus the small
two-divisor toric extension does not yet remove the obstruction.

For the degree-three triangle, a greedy modular deletion pass reduced the 23
equations to a 19-equation core that is still exactly empty over `Q`. Four
equations are redundant: two high-edge Jacobian equations, the constant
Jacobian equation (already fixed by `DF(0)=I`), and `P_y(0)=0`. The stored core
is not claimed inclusion-minimal because its single-deletion ideals become
positive-dimensional and expensive.

The subsequent translated two-divisor experiment is documented separately in
`TRANSLATED_TWO_DIVISOR.md`. Its central conclusion is that symmetric Laurent
boxes reproduce ordinary dense degree bounds exactly, so further progress
requires Newton-directed thin edge supports rather than larger boxes.
