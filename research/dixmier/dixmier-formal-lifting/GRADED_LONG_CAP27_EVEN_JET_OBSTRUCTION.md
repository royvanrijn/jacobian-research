# Entire even graded Long family at cap27: positive at h5, obstructed at h7

Use the canonical global polynomial Long tuple \(f=(R,T,D,S)\) in
source variables \((X,Q,P,Z)\), with \(\{P,X\}=\{Z,Q\}=1\), and the
standard Moyal product. Its input is
[long_images](../dixmier-audit/verify_long_pbw.py). Assign
\[
 \operatorname{wt}(X,Q,P,Z,h)=(-1,1,4,2,3),\qquad
 \operatorname{wt}(R,T,D,S)=(-1,2,4,1).                    \tag{1}
\]

**Theorem.** Over any characteristic-zero field there is no even
equivariant scalar full Weyl jet
\[
      F=f+h^2F_2+h^4F_4+h^6F_6+O(h^8),\qquad
      \deg F_{2,i},\deg F_{4,i},\deg F_{6,i}\leq27,        \tag{2}
\]
satisfying all six canonical commutators through \(h^7\).
Every coefficient in all four images is allowed, including
momentum-dependent deformations of the first image.

At this same cap the complete equations through \(h^5\) DO have
points. The proof classifies their field-valued points into two
components, retains every fourth-coefficient solution on both,
and proves a contradiction at \(h^7\) on each entire component.
This is stronger than testing only a fixed first particular.

Consequently an even full Long lift in (1), with a common spatial
bound, requires bound at least 28. There is no assertion that a
general lift can be normalized to (1), no odd-jet or nongraded
obstruction, no cap28 conclusion and no DC1/DC2 resolution.

## Complete coefficient maps

The coefficient at \(h^{2k}\), in image weight \(w_i\), has the entire
basis
\[
 X^\alpha Q^\beta P^\gamma Z^\delta,\qquad
 \alpha=\beta+4\gamma+2\delta+6k-w_i\geq0,\qquad
 \alpha+\beta+\gamma+\delta\leq27.                        \tag{3}
\]
No monomial from this space is omitted. Let
\[
        (A_{2k}U)_{ij}=\{U_i,f_j\}+\{f_i,U_j\},\quad i<j.
\]
For scalar Moyal symbols the coefficient at commutator order \(h^N\)
is
\[
 \sum_{\substack{a+b+\ell=N\\ \ell\geq1\ {\rm odd}}}
 \frac{\Pi^\ell(F_{a,i},F_{b,j})}{2^{\ell-1}\ell!},        \tag{4}
\]
with \(F_0=f\) and odd-indexed coefficients zero here. Moving known
terms to the negative right side gives
\[
                   A_2F_2=B_3,\quad A_4F_4=B_5,\quad
                   A_6F_6=B_7.                         \tag{5}
\]
In particular \(B_3=-\Pi^3(f_i,f_j)/24\), and
\[
 -B_{5,ij}=\{F_{2,i},F_{2,j}\}
 +\frac{\Pi^3(f_i,F_{2,j})+\Pi^3(F_{2,i},f_j)}{24}
 +\frac{\Pi^5(f_i,f_j)}{1920}.                           \tag{6}
\]
For the next step the complete known defect is
\[
 \begin{split}
 -B_{7,ij}={}&\{F_{2,i},F_{4,j}\}+\{F_{4,i},F_{2,j}\}\\
 &+\frac{\Pi^3(f_i,F_{4,j})+\Pi^3(F_{2,i},F_{2,j})
                           +\Pi^3(F_{4,i},f_j)}{24}\\
 &+\frac{\Pi^5(f_i,F_{2,j})+\Pi^5(F_{2,i},f_j)}{1920}
                           +\frac{\Pi^7(f_i,f_j)}{322560}.
 \end{split}                                           \tag{7}
\]

## The entire first correction has eight parameters

The first component counts are \((91,127,156,114)\), totaling 488.
The retained 59-coefficient particular \(g_0\) satisfies
\(A_2g_0=B_3\) exactly. An explicit 480-square minor of \(A_2\) has
determinant
\[
                         386834\pmod {1000003}.          \tag{8}
\]
Eight exact kernel tuples \(K_0,\ldots,K_7\) are retained. Their
coefficient matrix is the identity on these eight coordinates:
\[
 \begin{array}{c|c}
 0&D:X^4Q^2\\
 1&D:X^{14}Q^4P^2\\
 2&S:X^9Q^4\\
 3&S:X^{10}QP\\
 4&S:X^{10}Q^5\\
 5&S:X^{11}Q^2P\\
 6&S:X^{13}P^2\\
 7&S:X^{15}Q^6P .
 \end{array}                                           \tag{9}
\]
Thus the kernels are independent. Together with (8), they prove
rational rank 480 and give the COMPLETE family
\[
                        F_2=g_0+\sum_{j=0}^7k_jK_j.      \tag{10}
\]
No parameter is normalized away.

A compact alternative Hamiltonian basis is
\[
 R^3,\ R^4S,\ X^3,\ R^3(Y+2),\
                   X^4Q,\ R^3J_4,\ R^5T,\ X^5Z,         \tag{11}
\]
where
\[
 \begin{split}
 t&=XQ,\quad B=3X^2P+2(1-3XQ)Z,\\
 A_0&=Q+3XQ^2-XB/3,\quad Y=RA_0-3t,\\
 J_4&=9Y^4+66Y^3+36tY^2+150tY+171t^2-380t .
 \end{split}
\]
The replay verifies the invertible rational transition between
their Hamiltonian variations and (9).
In particular
\[
                         \{X^5Z,R\}=-3X^7
\]
is admitted at this cap. The first image is not held fixed.
No claim that a finite source shear preserves cap27 is used.

## Complete field-valued classification at h5

The \(F_4\) component counts are \((40,62,80,54)\), totaling 236.
All forty first-image coefficients remain unknown. Substituting
(10) into (6) produces 1722 support rows and 36 nonzero parameter
coefficient vectors.

Four retained rational row functionals annihilate EVERY column
of \(A_4\) and give the exact polynomials
\[
 \begin{split}
 k_5^2,\qquad k_6^2,\qquad k_7^2,\qquad
 L=k_2+\frac{10}{3}k_4+3k_5+9k_6-\frac{1340}{27}k_7 .
 \end{split}                                           \tag{12}
\]
Their supports have sizes \(196,43,195,195\), respectively.
These are individual squares, not a sum of squares. Over every
field, including arbitrary characteristic-zero complex extensions,
compatibility therefore forces
\[
                k_5=k_6=k_7=0,\qquad k_2=-10k_4/3.       \tag{13}
\]

After precisely this forced substitution, a fifth rational row
functional, supported on 180 rows, annihilates all columns of
\(A_4\) and evaluates the right side to
\[
                                  k_3k_4.               \tag{14}
\]
Thus the field-valued first-correction points lie in the two components
\[
 \begin{array}{c|cccc}
       &k_0&k_1&k_3&k_4\\ \hline
 A&u&v&t&0\\
 B&u&v&0&t .
 \end{array}                                           \tag{15}
\]
Their intersection is included. Equations (12)--(14) are exact rational
functional identities, not results inferred from a parameter sample
or an unretained Groebner computation. We classify field-valued
points; no equality of nonreduced ideals is asserted.

Explicitly, put
\[
         L_A=K_3,\qquad L_B=K_4-\frac{10}{3}K_2 .
\]
The two entire first families are
\[
                       F_2^{A,B}=g_0+uK_0+vK_1+tL_{A,B}. \tag{16}
\]
Both have zero first-image component as a consequence of the full
compatibility equations, including elimination of the admitted
\(X^5Z\) direction. This is not an imposed exact-\(R\) condition.

## Both entire fourth-coefficient families are retained

An explicit 235-square minor of \(A_4\) has determinant
\(673174\) modulo \(1000003\). Its exact nonzero kernel is
\[
                            K_{R^9}=(0,0,-9R^8,0).       \tag{17}
\]
Therefore its rank is exactly 235 and its kernel is precisely
the one-dimensional space in (17).

The witness specifies all rational coefficients of the following
two families:
\[
 \begin{split}
 F_4^A&=U_0+tU_1+t^2U_2+vU_3+vtU_4+wK_{R^9},\\
 F_4^B&=U_0+tV_1+vU_3+wK_{R^9}.                         \tag{18}
 \end{split}
\]
The respective coefficient counts of the displayed tuples are
\[
 \begin{array}{c|rrrrrr}
 A&74&25&60&5&18&9\\
 B&74&21&5&9&& .
 \end{array}
\]
All \(u,v,t,w\) remain parameters. The \(u\) terms cancel in the
particular solution without any gauge normalization.
Every tuple identity \(A_4F_4^{A,B}=B_5(F_2^{A,B})\) is verified
coefficientwise. In combination with rank \(235/236\) and (17),
these are the ENTIRE \(F_4\) families on both components. Hence
(12)--(18) give necessity and sufficiency at order \(h^5\).

The shared base tuple \(U_0\) has component degrees
\((-,19,27,15)\), with zero \(R\)-component and 74 coefficients.
It is an actual compatible polynomial jet through \(h^5\).
Together with [the complete cap26 obstruction](GRADED_LONG_CAP26_EVEN_JET_OBSTRUCTION.md),
this establishes spatial cap27 as the minimum for an even graded
jet through \(h^5\). This positive finite-order statement does not
give a full lift.

## The entire two-component family fails at h7

The \(F_6\) counts are \((13,24,34,20)\), totaling 91 columns.
In particular every one of the thirteen first-image coefficients
remains unknown. A retained 91-square minor has determinant
\(604791\) modulo \(1000003\), so \(A_6\) is injective on this full
space. The obstruction below does not need this injectivity.

Use (7) with each whole family (16)--(18). Both row sets have
974 entries. The nonzero parameter monomials in \(B_7\) are
\[
 \begin{array}{c|l}
 A&1,\ t,\ v,\ t^2,\ t^3,\ vt^2\\
 B&1,\ t,\ v .
 \end{array}                                           \tag{19}
\]
The \(u,w\) terms cancel identically; no transformation removed them.

For EACH component the witness retains a rational functional on
73 specified rows satisfying
\[
          \ell_AA_6=\ell_BA_6=0,\qquad
          \ell_AB_7^A=1,\quad \ell_BB_7^B=1.             \tag{20}
\]
The identities hold on every one of the 91 columns and on every
parameter-coefficient vector in (19). Therefore either full
component would require \(0=1\). This excludes their intersection
and every other field-valued point, proving the theorem.

A redundant 74-row common functional is also retained and replayed
on the union of both right-side coefficient sets.

## Replay boundary

The [witness](graded_long_cap27_even_witness.json) retains the first
particular, all eight normalized kernels and their polynomial
Hamiltonians, three rank minors, the five staged \(h^5\) functionals,
both complete fourth-coefficient families, and the two \(h^7\)
functionals. Row locators are explicit.
The [checker](check_graded_long_cap27_even.py) reconstructs all full
column spaces and parameter coefficients with Python Fraction
arithmetic. It verifies every identity above, the compact-basis
transition, and corrupted-input rejection controls.

At one nontrivial rational parameter point on EACH component, a
separate normal-PBW multiplication control verifies all six relations
through \(h^5\). These finite controls check conventions and the
retained positive families; the universal exclusion rests on the
exact coefficient identities (12), (14) and (20).

Run from the repository root:

    python3 research/dixmier/dixmier-formal-lifting/check_graded_long_cap27_even.py

The [receipt](graded_long_cap27_even_controls.json) pins this note,
checker, witness, canonical input and reused arithmetic engines.
Modular nonzero minors certify rational rank lower bounds; no
modular inconsistency is promoted alone to a rational obstruction.

This is an internal written proof with retained finite certificates.
It asserts no formal verification, independent whole-proof replay,
or external review. Its scope is exactly the even, graded, scalar,
cap27 family (2).
