# The extra degree-24 Long direction cannot solve the second quantum jet

Use the canonical full global polynomial Long tuple
\(f=(R,T,D,S)\), standard source brackets
\(\{P,X\}=\{Z,Q\}=1\), and Moyal ordering.
Assign source and parameter weights
\[
 \operatorname{wt}(X,Q,P,Z,h)=(-1,1,4,2,3),\qquad
 \operatorname{wt}(R,T,D,S)=(-1,2,4,1).
\]

**Theorem.** No even equivariant full scalar Long Weyl jet
\[
 F_i=f_i+h^2g_i+h^4k_i+O(h^6),\qquad
                      \deg g_i,\deg k_i\leq24
\]
satisfies all six canonical relations modulo \(h^7\).
Thus an even equivariant full lift with a common spatial bound
must have bound at least 25. Non-equivariant and odd-jet lifts,
and degree 25, are not excluded.

The complete coefficient enumeration, Moyal normalization and
target-gauge calculation are as in the
[cap-23 proof](GRADED_LONG_CAP23_JET_OBSTRUCTION.md).
This result retains the additional degree-24 direction rather
than removing it by an unjustified cap-preserving source gauge.

## The complete first family has one additional parameter

For the coefficient at \(h^{2m}\) of image weight \(w_i\), enumerate
exactly
\[
 X^aQ^bP^cZ^d,\quad
 a=b+4c+2d+6m-w_i\geq0,\qquad a+b+c+d\leq24.
\]
At \(m=1\), the counts are \((62,91,114,80)\), totaling 347.
The six equations are
\[
 A_3g=-E_3,\qquad E_{3,ij}=\Pi^3(f_i,f_j)/24.
\]
The matrix has 2129 rows. A retained \(344\)-square minor has
determinant \(381548\) modulo \(1000003\), so its rational rank
is at least 344.

Three exact independent rational kernels are
\[
 \begin{split}
 K_1&=(0,0,-3R^2,0)=\{R^3,f\},\\
 K_2&=(0,R^4,-4R^3S,0)=\{R^4S,f\},\\
 K_3&=\{X^3,f\}.
 \end{split}                                             \tag{1}
\]
All supports satisfy the cap and weights, and exact Fraction
substitution verifies each kernel. The first two are independent
as before. The \(S\)-component of \(K_3\) is nonzero, while those
of \(K_1,K_2\) vanish, so all three are independent. Their
existence gives the opposite rank bound 344.

The cap-23 particular \(g^0\), with its retained 59 coefficients,
still solves every first equation. Consequently the ENTIRE
first family over every characteristic-zero scalar extension is
\[
                  g=g^0+\alpha K_1+\beta K_2+\gamma K_3.  \tag{2}
\]
Each first-image component is zero, so \(r_2=0\) throughout.

The extra generator \(X^3\) is a source Hamiltonian. Its exact
source Weyl shear fixes \(X,Q,Z\) and sends
\(P\mapsto P-3\gamma h^2X^2\). This explains its first variation;
it does not justify removing \(\gamma\) while retaining a
degree-24 two-jet. We keep that parameter in every second equation.

## Only the two target gauges are removed

The exact target momentum shears from
\(G=\alpha R^3+\beta R^4S\) remove the first two parameters.
As \(r_2=0\), their only additional \(h^4\) term, up to its scalar
factor, is in \(D\):
\[
 R^3\bigl(g^0_S+\gamma(K_3)_S\bigr)
                            +\Pi^2(R^3,S)/8.             \tag{3}
\]
The constant part has degree at most 16. Since
\((K_3)_S=-3X^2\partial_P S\) has degree at most 12, the
\(\gamma\)-part has degree at most 21. Hence these target shears
preserve the required cap 24 through \(h^4\) and preserve the
grading. The ordering term in (3) is included explicitly.
No other \(h^4\) target-shear component changes.

Thus every hypothetical jet can be tested with
\(g=g^0+\gamma K_3\), without a source-gauge assumption.

## An exact rational functional rules out every parameter value

For \(g=g^0+\gamma K_3\), the second equation is
\[
                      A_5 k=b_0+\gamma b_1+\gamma^2b_2,   \tag{4}
\]
where
\[
 \begin{split}
 b_{0,ij}=-\bigl(&\{g^0_i,g^0_j\}
  +(\Pi^3(f_i,g^0_j)+\Pi^3(g^0_i,f_j))/24
  +\Pi^5(f_i,f_j)/1920\bigr),\\
 b_{1,ij}=-\bigl(&\{(K_3)_i,g^0_j\}+\{g^0_i,(K_3)_j\}\\
 &+(\Pi^3(f_i,(K_3)_j)+\Pi^3((K_3)_i,f_j))/24\bigr),\\
 b_{2,ij}=-&\{(K_3)_i,(K_3)_j\}.
 \end{split}
\]
The coefficient matrix allows every \(h^4\) monomial at cap 24.
Its counts are \((24,40,54,34)\), totaling 152, including all
24 first-image coefficients. The complete row set for
\([A_5\mid b_0,b_1,b_2]\) has 1289 rows.

The retained witness gives a rational linear functional \(\ell\)
on this row space, supported on 118 specified rows, such that
\[
             \ell A_5=0,\qquad
             \ell b_0=1,\qquad \ell b_1=\ell b_2=0.        \tag{5}
\]
The replay evaluates every one of these identities exactly
over Fraction arithmetic, on every full column of the generated
matrix, not merely on a selected minor or a numerical sample.
Applying (5) to (4) gives \(0=1\), for EVERY \(\gamma\).
This includes algebraic or transcendental choices of the
parameter after any characteristic-zero scalar extension.

The two target gauges were treated cap-preservingly, and the
extra source parameter was retained. Thus (5) covers (2)
completely. Higher corrections cannot change the coefficient
at \(h^5\), proving the theorem.

## Data, replay and scope

[The witness](graded_long_cap24_witness.json) retains the first
minor indices, prime and determinant, and every rational
coefficient and row key of \(\ell\). It reuses the frozen
[59-coefficient first particular](graded_long_cap23_witness.json).
[The checker](check_graded_long_cap24.py) reconstructs both
complete coefficient families, verifies the three kernels and
particular exactly, checks the nonzero first minor with bare
modular elimination, and evaluates (5) over the rationals.
The replay also rejects deliberately corrupted particular and functional data.
FLINT was used to discover the data, but is unnecessary for replay.
[The receipt](graded_long_cap24_controls.json) pins those inputs.

Evenness refers to the Moyal symbols. PBW reordering preserves a spatial
cap but need not preserve parameter parity.

This is a complete finite even/graded two-jet obstruction with
all four images allowed to deform; it is not an exact-\(R\)
hypothesis or a scalar finite-order mate exclusion. It gives
no non-equivariant or odd-jet conclusion, no bounded
normalization theorem, no cap-25 exclusion and no DC1/DC2
resolution. Independent whole-proof replay, formal verification,
external review and novelty priority remain unclaimed.
