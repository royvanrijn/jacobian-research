# An admitted Airy-Long cap25 matrix jet and its fixed next-step obstruction

The canonical Long symbols admit a concrete rank-two polynomial matrix
jet, with common Bernstein degree 25, satisfying every Weyl relation
through \(h^4\). This is an actual finite jet in the original global
source algebra. Its next equation has no degree-25 correction when
these lower coefficients are fixed, even allowing every matrix entry
and every first-image coefficient.

The positive jet is compositional. The fixed extension obstruction
does not exclude other Airy matrix jets, fourth-order gauge choices,
or changes in the earlier scalar parameters.

## Exact input and formulas

Use source order \((X,Q,P,Z)\), brackets
\(\{P,X\}=\{Z,Q\}=1\), and the standard matrix Moyal product.
Let \(f=(R,T,D,S)\) be
[the canonical full Long polynomial tuple](../dixmier-audit/verify_long_pbw.py).
The target brackets are \(\{D,R\}=\{S,T\}=1\).
Assign source weights \((-1,1,4,2)\), parameter weight 3, and image
weights \((-1,2,4,1)\).

The scalar polynomial tuple \(g_0\), with all 59 rational coefficients
retained in [the earlier witness](graded_long_cap23_witness.json),
satisfies
\[
 \{g_{0,i},f_j\}+\{f_i,g_{0,j}\}
                          =-\Pi^3(f_i,f_j)/24.            \tag{1}
\]
Set
\[
 \begin{split}
 t&=XQ,\quad M(t)=\begin{pmatrix}0&1\\t&0\end{pmatrix},
          \quad V=Q\partial_P+X\partial_Z,\\
 B&=3X^2P+2(1-3XQ)Z,\quad
 A_0=Q+3XQ^2-XB/3,\quad Y=RA_0-3t,\\
 G_4&=R^3(Y+2)=R^4a,\qquad K_i=\{G_4,f_i\},\\
 \delta&=-2/27,\qquad g_i=g_{0,i}+\delta K_i .
 \end{split}                                           \tag{2}
\]
Here \(G_4\) is an actual original polynomial; the phase-chart
notation \(a\) is not needed for its definition. The Hamiltonian
variation \(K\) solves the homogeneous equation associated with (1).

Define
\[
 F_{0,i}=f_iI_2,\qquad
 F_{1,i}=M(t)Vf_i,\qquad
 F_{2,i}=\left(g_i+\frac t2V^2f_i\right)I_2.              \tag{3}
\]
Thus the first correction is genuinely matrix valued, while the
second correction in this chosen jet is scalar.

## The actual third coefficient

Let \(\iota\) be the admitted source representation
\[
 X\mapsto XI_2,\quad Q\mapsto QI_2,\quad
 P\mapsto PI_2+hQM(t),\quad Z\mapsto ZI_2+hXM(t).
\]
Its exact source relations and boundedness scope are proved in
[the polynomial connection theorem](POLYNOMIAL_CONNECTION_BOUNDEDNESS_REFLECTION.md).
Write \(J_n(f)\) for the \(h^n\) Weyl-symbol coefficient of
\(\iota(\operatorname{Weyl}(f))\). Then
\[
 J_0(f)=fI_2,\quad J_1(f)=MVf,\quad
 J_2(f)=\frac t2V^2f\,I_2 .
\]
The retained third coefficient is
\[
                              F_{3,i}=J_3(f_i)+J_1(g_i). \tag{4}
\]

For an explicit finite recipe for \(J_3\), put
\(M_X=QM(t)\), \(M_Q=XM(t)\).
With multi-indices \(\alpha\in\mathbb N^2\), initialize
\(U_0=W_0=I_2\) and compute
\[
 U_{\alpha+e_j}=\partial_jU_\alpha-U_\alpha M_j,\qquad
 W_{\alpha+e_j}=M_jW_\alpha+\partial_jW_\alpha.
\]
Define
\[
 T_\gamma=\sum_{\alpha\leq\gamma}
 (-1)^{|\alpha|}\binom\gamma\alpha U_\alpha W_{\gamma-\alpha},
 \qquad
 J_n(f)=\sum_{|\gamma|=n}
       \frac{T_\gamma\,\partial_{P,Z}^\gamma f}{2^n\gamma!}. \tag{5}
\]
Only indices of size at most three are used. One way to derive (5)
is to choose a formal position solution \(\partial_jB=BM_j\);
then \(U_\alpha=(\partial^\alpha B^{-1})B\),
\(W_\alpha=B^{-1}\partial^\alpha B\), and expand
\(B^{-1}\star f\star B\). The resulting coefficients in (5) are
finite polynomials and independent of any analytic or global gauge
existence. The checker regenerates them and verifies the resulting
Weyl identities directly.

For the unmodified \(g_0\), the compositional \(F_3\) has degree 27.
The choice \(\delta=-2/27\) in (2) cancels its degree-27 and
degree-26 terms. The actual coefficient degrees are
\[
 \begin{array}{c|rrrr}
       &R&T&D&S\\ \hline
 F_1&-&17&25&13\\
 F_2&-&17&25&13\\
 F_3&-&19&25&-
 \end{array}                                           \tag{6}
\]
A dash denotes a zero matrix. All entries are original global
polynomials. Their entire rational coefficient lists are retained
in the [witness](airy_long_cap25_matrix_jet_witness.json).

Exact full matrix Moyal multiplication gives, for
\(F_i=F_{0,i}+hF_{1,i}+h^2F_{2,i}+h^3F_{3,i}\),
\[
       [F_D,F_R]_\star=[F_S,F_T]_\star=hI_2+O(h^5),
       \qquad\hbox{all cross commutators }O(h^5).         \tag{7}
\]
This verifies all six relations at every coefficient through \(h^4\).
It does not assert that the truncation is an exact Weyl tuple.

## Complete next-coefficient space for this fixed point

Hold every coefficient (2)--(4) fixed and seek \(F_4\) at spatial
cap25 in the same grading. For image weight \(w_i\), the entire scalar
basis is
\[
 X^\alpha Q^\beta P^\gamma Z^\nu,\qquad
 \alpha=\beta+4\gamma+2\nu+12-w_i\geq0,\quad
 \alpha+\beta+\gamma+\nu\leq25 .
\]
Its component counts are \((29,47,62,40)\), totaling 178.
Allowing all four matrix entries gives 712 columns. All first-image
entries are included.

Let \(\delta_i=\{f_i,-\}\) and \(A_i=F_{1,i}\).
The unknown coefficient contributes at \(h^5\) by the covariant map
\[
 (\mathcal L C)_{ij}
  =\delta_iC_j-\delta_jC_i+[A_i,C_j]-[A_j,C_i].            \tag{8}
\]
This uses ordinary matrix commutators in the last two terms.
The full known \(h^5\) right side is calculated with the matrix Moyal
formula, including even bidifferential orders when applicable; the
scalar odd-only commutator formula is not imposed on matrix entries.

The complete map has 7100 support rows and 239744 nonzero column
entries. The right side has 402 nonzero entries. The witness specifies
an exact rational functional \(\ell\), supported on 274 matrix rows,
such that
\[
                  \ell\mathcal L=0\quad
                  \hbox{on all 712 columns},\qquad
                  \ell(\text{right side})=1.             \tag{9}
\]
Thus no \(F_4\) in this whole space extends the fixed lower point
through commutator order \(h^5\).

The functional was obtained by the classical scalar trace projection:
ordinary matrix commutators in (8) have zero trace, and the projected
scalar map has rank \(177/178\), augmented rank 178. The retained proof
does not rely on a numerical rank discrepancy. The checker verifies
(9) coefficientwise on EVERY full matrix column, including off-diagonal
and first-image columns. No cyclic quantum trace is assumed.

## Exact scope and replay

The [checker](check_airy_long_cap25_matrix_jet.py) regenerates the point,
checks its grading and degrees, verifies all six relations (7), builds
all 712 columns of (8), checks (9), and rejects a corrupted functional.
The [receipt](airy_long_cap25_matrix_jet_controls.json) pins the complete
input and replay sources. Run from the repository root:

    python3 research/dixmier/dixmier-formal-lifting/check_airy_long_cap25_matrix_jet.py

The positive construction remains compositional; no boundedness gain
over its scalar formal input is asserted. The negative conclusion fixes
\(F_1,F_2,F_3\) exactly. Other homogeneous \(F_3\) directions, different
\(F_2\) scalar or traceless parameters, noncompositional matrix jets and
larger caps are not excluded. Neither formal covariant gauge rigidity
nor this fixed obstruction implies a bounded gauge theorem or a
classification of all matrix carriers.

No bounded full Long lift, scalar realization, DC1/DC2 resolution,
formal verification, independent whole-proof replay or external review
is claimed.
