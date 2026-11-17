# No degree-24 graded full Long jet, including the complex branches

Use the canonical global polynomial Long tuple \(f=(R,T,D,S)\) in
source variables \((X,Q,P,Z)\), with \(\{P,X\}=\{Z,Q\}=1\), and the
standard Moyal product. Its input is
[long_images](../dixmier-audit/verify_long_pbw.py). Assign
\[
 \operatorname{wt}(X,Q,P,Z,h)=(-1,1,4,2,3),\qquad
 \operatorname{wt}(R,T,D,S)=(-1,2,4,1).                    \tag{1}
\]

**Theorem.** Over any characteristic-zero field containing the input
coefficients, there is no equivariant scalar full Weyl jet
\[
 F_i=f_i+hF_{1,i}+h^2F_{2,i}+h^3F_{3,i}+h^4F_{4,i}+O(h^5),
 \qquad \deg F_{m,i}\leq24,                              \tag{2}
\]
satisfying the six canonical commutators through \(h^5\).
No parity restriction is imposed. All four images, including every
admitted first-image coefficient, may deform.

Thus a full Long lift in the grading (1), with a common spatial-degree
bound, requires bound at least 25. This asserts no normalization of an
arbitrary bounded lift into (1), no nongraded obstruction, no
full-parity cap-25 exclusion, and no DC1/DC2 resolution.

## The entire lower family, including both field branches

The [cap24 complex first-jet theorem](GRADED_LONG_CAP24_COMPLEX_FIRST_JET_FAMILY.md)
classifies the full first-coefficient space and all second coefficients.
Only its uniformly cap-preserving source \(XQ\) and target \(RS\)
shears are used to normalize. Before the nonlinear conditions the
remaining first Hamiltonian is
\[
 G=aY+bR^2T+cH_4+dH_6+eX^2Z+fH_5.                       \tag{3}
\]
For notation, put
\[
 \begin{split}
 t&=XQ,\quad B=3X^2P+2(1-3XQ)Z,\\
 A_0&=Q+3XQ^2-XB/3,\quad Y=RA_0-3t,\\
 H_4&=Y^4+4tY^2-\frac{16}{3}tY+10t^2,\\
 H_6&=R^4D+
 \frac{9Y^5+45tY^3+45t^2Y-60tY^2+110tY-720t^2}{25},\\
 H_5&=9Y^5+45tY^3-60tY^2+135t^2Y-10tY-555t^2.
 \end{split}                                           \tag{4}
\]
The complete \(h^3\) field-valued compatibility conditions are
\[
                d^2-250df+35625f^2=0,\qquad e^2=0,\qquad ef=0.         \tag{5}
\]
Over fields \(e=0\). Adjoin \(\omega\), \(\omega^2=-2\), if necessary.
Every point of (5) is covered by
\[
 f=z,\qquad d=(125+100\omega)z
 \quad\hbox{or}\quad d=(125-100\omega)z.                 \tag{6}
\]
The branch \(z=0\) is included. Both genuine complex branches are
retained; there is no restriction to rational or real parameter values.

Write
\[
 F_1=\{aY+bR^2T+cH_4+
                z(H_5+(125+100\omega)H_6),\,f\}.         \tag{7}
\]
The complete \(F_2\) family on this branch has eleven particular
coefficient tuples, with monomials
\[
        1,a^2,ab,ac,b^2,bc,c^2,az,bz,cz,z^2,              \tag{8}
\]
plus
\[
 uK_{R^3}+vK_{R^4S}+\gamma K_{X^3},\qquad K_G=\{G,f\}.  \tag{9}
\]
All coefficients are retained in the prior complex witness as pairs
in \(\mathbb Q+\omega\mathbb Q\). The new replay verifies these
identities coefficientwise before using them. Conjugating (7)--(9)
gives the second branch. No further source or target shear removes
\(a,b,c,z,u,v,\gamma\).

## Complete spaces and coefficient equations

The coefficient at \(h^m\) of image weight \(w_i\) has exactly the basis
\[
 X^\alpha Q^\beta P^\gamma Z^\delta,\qquad
 \alpha=\beta+4\gamma+2\delta+3m-w_i\geq0,\qquad
 \alpha+\beta+\gamma+\delta\leq24.                        \tag{10}
\]
Here \(\gamma\) in a monomial exponent is unrelated to the parameter
in (9). The linear coefficient map \(A_m\) sends a tuple \(U\) to
the six polynomials
\(\{U_i,f_j\}+\{f_i,U_j\}\).
Every column from (10) is used.

For scalar Moyal symbols, the coefficient at commutator order \(h^N\)
is
\[
 \sum_{\substack{a+b+\ell=N\\ \ell\geq1\ {\rm odd}}}
 \frac{\Pi^\ell(F_{a,i},F_{b,j})}{2^{\ell-1}\ell!},
 \qquad F_0=f.                                         \tag{11}
\]
The terms involving \(F_{N-1}\) give \(A_{N-1}\);
put every known term on the negative right side \(B_N\).
In particular,
\[
 \begin{split}
 -B_{4,ij}={}&\{F_{1,i},F_{2,j}\}+\{F_{2,i},F_{1,j}\}\\
 &+\frac{\Pi^3(f_i,F_{1,j})+\Pi^3(F_{1,i},f_j)}{24},     \tag{12}\\
 -B_{5,ij}={}&\{F_{1,i},F_{3,j}\}+\{F_{2,i},F_{2,j}\}
                          +\{F_{3,i},F_{1,j}\}\\
 &+\frac{\Pi^3(f_i,F_{2,j})+\Pi^3(F_{1,i},F_{1,j})
                          +\Pi^3(F_{2,i},f_j)}{24}
                           +\frac{\Pi^5(f_i,f_j)}{1920}. \tag{13}
 \end{split}
\]
The \(\Pi^3(F_1,F_1)\) term is essential; this is not the even-jet
equation.

## Order \(h^4\) excludes both complex branches

The \(F_3\) space has component counts \((40,62,80,54)\), hence
236 columns. Using the full family (7)--(9), its row set has
1673 entries. A retained rational functional \(\ell_{\rm field}\),
supported on 187 specified rows, satisfies exactly
\[
             \ell_{\rm field}A_3=0,\qquad
             \ell_{\rm field}B_4=z^3+\omega\,0.           \tag{14}
\]
These are polynomial identities in all seven parameters. The
functional applies to the combined field-valued equation. We do
**not** impose separate real and imaginary equations on complex
parameter values.

The entire right side is also replayed under \(\omega\mapsto-\omega\).
The same rational functional gives \(z^3-\omega\,0\). Consequently
both branches in (6) force \(z=0\) over any field. This conclusion
does not rely on sampling or on a reality assumption.

After \(z=0\), the remaining coefficients have rational polynomial
expressions in \(a,b,c,u,v,\gamma\). A second retained rational
functional, supported on 180 rows, gives
\[
 \ell_LA_3=0,\qquad
 \ell_LB_4=\left(c-\frac{27a}{544}\right)^3.              \tag{15}
\]
Thus compatibility over fields requires
\[
                             a=\frac{544}{27}c.         \tag{16}
\]
Equations (14)--(16) classify parameter points; no assertion about
nilpotent parameters or equality of nonreduced schemes is needed.

## The whole surviving third coefficient

After (16), the retained seven particular tuples have monomials
\[
                   c,\ b,\ c\gamma,\ cv,\ b\gamma,\
                                c^3,\ bc^2.             \tag{17}
\]
All 253 nonzero rational coefficients are specified in the new witness.
Direct substitution verifies the complete polynomial identity
\(A_3F_3=B_4\) for each coefficient in (17). The component degrees are
\[
 \begin{array}{c|ccc}
     &T&D&S\\ \hline
 c&13&21&-\\
 b&-&14&-\\
 c\gamma&18&24&14\\
 cv&-&22&-\\
 b\gamma&-&19&-\\
 c^3&21&23&-\\
 bc^2&-&22&-
 \end{array}
\]
Every tuple has zero \(R\)-component.

An explicit \(235\)-square minor is nonzero modulo \(1000003\),
with determinant \(43500\), so
\(\operatorname{rank}_{\mathbb Q}A_3\geq235\).
The exact nonzero kernel
\[
                         K_{R^6}=(0,0,-6R^5,0)         \tag{18}
\]
gives the opposite bound. Thus the rank is exactly 235 and the
ENTIRE surviving family is
\[
                 F_3=P_3(b,c,v,\gamma)+wK_{R^6}.        \tag{19}
\]
All \(b,c,u,v,\gamma,w\) remain. In particular no cap-preserving
normalization of the \(v\) or \(\gamma\) direction is assumed.
The disappearance of \(u\) from the particular coefficient is an
identity, not a removal of that parameter.

## Order \(h^5\) gives a constant contradiction

Allow all 152 coefficients of \(F_4\) from (10), with component
counts \((24,40,54,34)\). In particular all 24 first-image
coefficients remain unknown. The row set has 1289 entries.

The full right side of (13), using (19), has exactly the twelve
parameter monomials
\[
 1,\gamma,v,\gamma^2,\gamma v,c^2,bc,c^2\gamma,c^2v,
                       bc\gamma,c^4,bc^3.              \tag{20}
\]
The \(u,w\) terms cancel identically without a gauge transformation.

The retained rational functional \(\ell_5\), supported on 119
specified rows, has the exact identities
\[
 \ell_5A_4=0\quad\hbox{on all 152 columns},\qquad
 \ell_5B_{5,1}=1,\qquad
 \ell_5B_{5,M}=0\quad(M\ne1).                            \tag{21}
\]
Therefore the coefficient equation would require \(0=1\) at every
remaining parameter point. This proves the theorem.

For completeness, \(A_4\) has rank 151, not full column rank:
a retained \(151\)-square minor has determinant \(765461\) modulo
\(1000003\), and
\(K_{R^9}=(0,0,-9R^8,0)\) is its exact nonzero kernel.
This rank statement is not needed for contradiction (21).

## Retained evidence and replay boundary

The [witness](graded_long_cap24_full_jet_witness.json) stores all
three rational functionals with row locators, both rank minors,
and the full seven-tuple family (17). The
[checker](check_graded_long_cap24_full_jet.py) rebuilds every complete
column and every parameter coefficient from the global Long input,
using exact Fraction arithmetic. It checks (14) for both field
conjugates, (15), all of (19), and every identity (21). It also rejects
deliberately corrupted functionals and a corrupted particular tuple.

A separate normal-PBW multiplication control checks all six
commutators through \(h^4\) at a nontrivial surviving rational
parameter point. This is a finite convention/family control, not
the evidence for the universal parameter identities.

Run from the repository root:
\[
 \texttt{python3 research/dixmier/dixmier-formal-lifting/}
 \texttt{check\_graded\_long\_cap24\_full\_jet.py}.
\]
The [receipt](graded_long_cap24_full_jet_controls.json) pins this note,
the witness, checker, canonical input, lower-family source, and reused
polynomial engines. Modular nonzero minors certify rational rank lower
bounds; all obstruction functionals themselves are verified over
\(\mathbb Q\). No modular inconsistency is promoted alone to a
characteristic-zero obstruction.

This is an internal written proof with retained exact finite
certificates. It has not received formal verification, independent
whole-proof replay, or external review. The theorem concerns precisely
the original polynomial, graded, cap24 family (2).
