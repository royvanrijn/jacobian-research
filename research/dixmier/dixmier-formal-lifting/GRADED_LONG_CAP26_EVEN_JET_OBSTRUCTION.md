# Complete even graded Long obstruction at spatial cap 26

Use the canonical global polynomial Long tuple \(f=(R,T,D,S)\) in
source variables \((X,Q,P,Z)\), with \(\{P,X\}=\{Z,Q\}=1\), and the
standard Moyal product. The input is
[long_images](../dixmier-audit/verify_long_pbw.py). Assign
\[
 \operatorname{wt}(X,Q,P,Z,h)=(-1,1,4,2,3),\qquad
 \operatorname{wt}(R,T,D,S)=(-1,2,4,1).                    \tag{1}
\]

**Theorem.** Over any characteristic-zero field there is no even
equivariant scalar full Weyl jet
\[
             F=f+h^2F_2+h^4F_4+O(h^6),\qquad
             \deg F_{2,i},\deg F_{4,i}\leq26,             \tag{2}
\]
satisfying all six canonical commutators through \(h^5\).
Every coefficient in all four images is allowed. In particular,
the first image is allowed to deform.

The proof uses the complete first-coefficient kernel, which grows
from four directions at cap25 to six at cap26. All six parameters
are retained. No source or target gauge is removed.

Consequently an even full Long lift in (1), with a common spatial
degree bound, needs bound at least 27. This makes no claim about
odd first jets, normalization into the grading, nongraded maps,
or DC1/DC2. A finite polynomial jet of common cap39 through the
same commutator order is constructed below as a positive control.

## Entire coefficient spaces and equations

For image weight \(w_i\), the full coefficient basis at \(h^{2k}\) is
\[
 X^\alpha Q^\beta P^\gamma Z^\delta,\qquad
 \alpha=\beta+4\gamma+2\delta+6k-w_i\geq0,\qquad
 \alpha+\beta+\gamma+\delta\leq26.                        \tag{3}
\]
Every such monomial is used. Let \(A_{2k}\) denote the six-pair
linear map
\[
       (A_{2k}U)_{ij}=\{U_i,f_j\}+\{f_i,U_j\},\quad i<j.
\]
The two equations needed for (2) are
\[
 \begin{split}
 A_2F_2&=B_3,\qquad B_{3,ij}=-\Pi^3(f_i,f_j)/24,          \tag{4}\\
 A_4F_4&=B_5(F_2),                                      \tag{5}\\
 -B_{5,ij}(F_2)&=\{F_{2,i},F_{2,j}\}
       +\frac{\Pi^3(f_i,F_{2,j})+\Pi^3(F_{2,i},f_j)}{24}
       +\frac{\Pi^5(f_i,f_j)}{1920}.
 \end{split}
\]
These follow directly from the Moyal commutator. Here the coefficient
is even; equation (5) is not used for an odd first-jet family.

## Complete first family, with six exact kernel directions

The component counts in (3) at \(k=1\) are
\((80,114,141,102)\), totaling 437 columns. The retained particular
\(g_0\), with all 59 rational coefficients, is the same exact tuple
used in [the cap25 theorem](GRADED_LONG_CAP25_JET_OBSTRUCTION.md).
Direct substitution verifies \(A_2g_0=B_3\).

The witness specifies a 431-square minor of \(A_2\) with determinant
\[
                        181171\pmod {1000003}.           \tag{6}
\]
All denominators are invertible at this prime. It also retains six
exact independent rational kernel tuples \(K_0,\ldots,K_5\).
For clarity, their normalization is fixed by the identity coefficient
matrix on the following six coordinates, in the displayed order:
\[
 \begin{array}{c|c}
 0&D:\ X^4Q^2\\
 1&D:\ X^{14}Q^4P^2\\
 2&S:\ X^9Q^4\\
 3&S:\ X^{10}QP\\
 4&S:\ X^{11}Q^2P\\
 5&S:\ X^{15}Q^6P .
 \end{array}                                           \tag{7}
\]
The replay checks each entire tuple, not just these coordinates.
The rank lower bound (6) and the six independent exact kernels give
\(\operatorname{rank}_{\mathbb Q}A_2=431\). Thus the whole first family is
\[
                       F_2=g_0+\sum_{j=0}^5k_jK_j.       \tag{8}
\]
The scalar parameters \(k_j\) are unrestricted. All \(R\)-components
in (8) vanish as a computed consequence of the full equation (4);
this was not an exact-\(R\) hypothesis.

The same kernel has a compact polynomial Hamiltonian basis. Put
\[
 \begin{split}
 t&=XQ,\quad B=3X^2P+2(1-3XQ)Z,\\
 A_0&=Q+3XQ^2-XB/3,\quad Y=RA_0-3t,\\
 J_4&=9Y^4+66Y^3+36tY^2+150tY+171t^2-380t .
 \end{split}
\]
Then the six Hamiltonians may be taken as
\[
              R^3,\quad R^4S,\quad X^3,\quad
              R^3(Y+2),\quad X^4Q,\quad R^3J_4.          \tag{9}
\]
The expression \(R^3(Y+2)=R^4a\) is an actual polynomial in the
original variables. No rational chart is required to define it.
The checker verifies the exact invertible transition between the
Hamiltonian variations in (9) and the normalized tuples in (7).
It also verifies the polynomial Hamiltonian generating every \(K_j\).

The two new directions in (9) are included in (8). There is no
assumption that earlier target-gauge arguments preserve cap26.

## Complete next equation and two rational functionals

The \(F_4\) component counts are \((34,54,71,47)\), totaling 206
columns. In particular all 34 first-image coefficients remain
unknown. Substituting the entire family (8) into (5) yields 1590
support rows and twenty nonzero parameter-coefficient vectors.

Two retained rational row functionals \(\ell_0,\ell_1\), supported
on 168 and 167 specified rows, satisfy
\[
 \begin{split}
 \ell_0A_4&=0,&\quad \ell_0B_5&=p_0
                  =1-\frac{32}{3969}k_3k_4,\\
 \ell_1A_4&=0,&\quad \ell_1B_5&=p_1=k_4^2.              \tag{10}
 \end{split}
\]
Both identities are checked over \(\mathbb Q\) on every one of the
206 coefficient columns and all twenty parameter coefficients.
They are not consequences inferred from modular inconsistency or
from sampled parameter values.

Let \(a=32/3969\). The exact polynomial identity
\[
                 1=(1+ak_3k_4)p_0+a^2k_3^2p_1          \tag{11}
\]
shows that these two necessary compatibility polynomials generate
the unit ideal. Therefore (5) has no solution with (8), including
over every complex or algebraic parameter extension. This proves
the theorem. No root branch is removed and no parameter is assumed
generic.

For completeness, \(A_4\) has rank \(205/206\). Its retained
205-square minor has determinant
\[
                         732211\pmod {1000003},
\]
and its exact nonzero kernel is
\(K_{R^9}=(0,0,-9R^8,0)\). This tail rank is not needed for the
contradiction (11).

## Positive polynomial control by weighted de Rham homotopy

The bounded-family obstruction does not obstruct the same linear
equation without cap26. Fix \(F_2=g_0\), and set
\(E_{ij}=B_{5,ij}(g_0)\), \(E_{ji}=-E_{ij}\).
The replay checks its closedness:
\[
       \{f_i,E_{jk}\}+\{f_j,E_{ki}\}+\{f_k,E_{ij}\}=0 .
\]
Define the four polynomial multipliers
\[
                  t=(-4D,-S,-R,2T)
\]
and the actual finite polynomial tuple
\[
                        C_i=-\frac19\sum_jt_jE_{ji}.    \tag{12}
\]
Direct exact substitution, separately for all six pairs, verifies
\[
                             A_4C=B_5(g_0).
\]
Its component degrees \((R,T,D,S)\) are
\[
                              (19,31,39,27).            \tag{13}
\]
Thus \(f+h^2g_0+h^4C\) is a concrete even polynomial jet satisfying
all canonical commutators through \(h^5\), of common cap39.
In particular its first image has an actual nonzero \(h^4\)
correction. Formula (12) is regenerated by the checker from the
input, rather than read from a precomputed correction list.
This is the weighted de Rham homotopy applied to this closed defect.
It is a finite-order positive control; the truncation has a nonzero
commutator defect at order \(h^7\), with 434 source-monomial terms
and maximum degree 48. It supplies no all-order finite lift or
termination theorem.

## Retained replay and scope

The [witness](graded_long_cap26_even_witness.json) retains the first
particular, six full kernel tuples and their Hamiltonians, both
rank minors, and both functionals with precise row locators.
The [checker](check_graded_long_cap26_even.py) uses ordinary Python
Fraction arithmetic and reconstructs all columns and right sides
from the canonical Long input. It checks (4), completeness of (8),
the compact-basis transition, every identity (10), the polynomial
Bezout certificate (11), and positive control (12). Corrupted
particular and functional controls must fail.

Run from the repository root:

    python3 research/dixmier/dixmier-formal-lifting/check_graded_long_cap26_even.py

The [receipt](graded_long_cap26_even_controls.json) pins the proof,
checker, witness, input and reused exact arithmetic engines.
Modular nonzero minors certify rational rank lower bounds.
The obstruction itself is checked by exact rational identities.

This is an internal written proof with retained finite certificates.
It has no formal verification, independent whole-proof replay or
external review claim. The only excluded family is the complete
even, graded, scalar, cap26 family (2).
