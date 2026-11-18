# Complete cap25 matrix obstruction for the prescribed Airy first coefficient

Fix the canonical global Long symbols \(f=(R,T,D,S)\), source
\((X,Q,P,Z)\), brackets \(\{P,X\}=\{Z,Q\}=1\), and the matrix Moyal
product. Put
\[
 t=XQ,\quad M(t)=\begin{pmatrix}0&1\\t&0\end{pmatrix},
 \quad V=Q\partial_P+X\partial_Z,\qquad A_i=M(t)Vf_i.       \tag{1}
\]
Use source weights \((-1,1,4,2)\), parameter weight 3, and image
weights \((-1,2,4,1)\). Matrix indices have weight zero.

**Theorem.** There is no graded \(2\)-by-\(2\) matrix jet
\[
 F_i=f_iI_2+hA_i+h^2B_i+h^3C_i+h^4E_i+O(h^5),            \tag{2}
\]
with common spatial degree at most 25, satisfying all six canonical
Weyl commutators through \(h^5\).
The first matrix coefficient \(A\) is prescribed by (1).
EVERY matrix entry in \(B,C,E\), including every first-image entry,
is allowed. No parity or diagonal ansatz is imposed on these later
coefficients, and no source or target gauge is removed.

This closes one complete bounded family beyond its positive finite
jet. It is not a theorem about other first coefficients, other
connections, larger caps, nongraded matrices, or arbitrary finite graph
carriers. It gives no bounded formal-gauge theorem or DC1/DC2 result.

## The known positive point

[The admitted matrix-jet source](AIRY_LONG_CAP25_MATRIX_JET_AND_FIXED_EXTENSION_OBSTRUCTION.md)
retains exact global polynomial matrices \(B^*,C^*\) of common cap25
satisfying every relation through \(h^4\). They arise from the
polynomial-connection source embedding with \(M(t)\). In its notation,
\[
 \begin{split}
 g&=g_0-\frac2{27}K_{R^4a},\\
 B_i^*&=\left(g_i+\frac t2V^2f_i\right)I_2,\\
 C_i^*&=J_3(f_i)+M(t)Vg_i .
 \end{split}                                           \tag{3}
\]
Here \(g_0\) is the exact 59-coefficient scalar first particular,
and \(K_G=(\{G,f_i\})_i\). The expression
\(R^4a=R^3(Y+2)\) is an original global polynomial:
\[
 \begin{split}
 B_0&=3X^2P+2(1-3XQ)Z,\\
 Y&=R(Q+3XQ^2-XB_0/3)-3XQ .
 \end{split}
\]
The letter \(B_0\) here is only the displayed source polynomial,
not the matrix coefficient \(B_i\).

The coefficient \(B^*\) is scalar and \(C^*\) is traceless.
The positive point alone does not exclude other later coefficients;
the following full kernel and compatibility calculation does.

## Complete covariant coefficient spaces

Let
\[
 \delta_i=\{f_i,-\},\qquad \mathcal D_i=\delta_i+\operatorname{ad}A_i .
\]
The unknown coefficient at \(h^m\), \(m\geq2\), contributes at
commutator order \(h^{m+1}\) by
\[
                 (\mathcal LC)_{ij}
                        =\mathcal D_iC_j-\mathcal D_jC_i. \tag{4}
\]
The scalar and traceless parts split under this map.

All scalar monomials for image weight \(w_i\), order \(m\), are
\[
 X^\alpha Q^\beta P^\gamma Z^\nu,\qquad
 \alpha=\beta+4\gamma+2\nu+3m-w_i\geq0,\quad
 \alpha+\beta+\gamma+\nu\leq25.                           \tag{5}
\]
No support is omitted.

For a traceless coefficient \((a,b,c)\), representing
\(\left(\begin{smallmatrix}a&b\\c&-a\end{smallmatrix}\right)\), the
covariant derivative is exactly
\[
 \mathcal D_i(a,b,c)=
 (\delta_i a+\alpha_i(c-tb),\
  \delta_i b-2\alpha_i a,\
  \delta_i c+2t\alpha_i a),\qquad \alpha_i=Vf_i.           \tag{6}
\]
This gives the entire sparse column map used below.

At \(m=2\), (5) has component counts \((71,102,127,91)\), totaling
391. The scalar map has rank 387 and exactly the four-dimensional
kernel
\[
                  K_{R^3},\quad K_{R^4S},\quad K_{X^3},
                                      \quad K_{R^4a}.    \tag{7}
\]
The prior exact scalar witness supplies a nonzero rank-387 minor,
and the replay checks all four independent polynomial kernel tuples.

The traceless map has 1173 columns and 8364 rows. The retained
1173-square minor is nonzero modulo \(1000003\), with determinant
\(936633\). Thus it is injective over \(\mathbb Q\) and every
characteristic-zero field extension. In particular there is NO
traceless homogeneous second coefficient in this bounded space.

The known \(h^3\) defect from \(fI+hA\) is scalar. Using (3), (7),
and injectivity, the ENTIRE second coefficient is therefore
\[
 B=B^*+uK_{R^3}I_2+vK_{R^4S}I_2+
                    \gamma K_{X^3}I_2+eK_{R^4a}I_2.     \tag{8}
\]
Here \(e=\delta+2/27\) is centered at the positive point. Every
parameter in (8) is retained before the next equation.

## The complete h4 equation forces two parameters

At \(m=3\), (5) has component counts \((47,71,91,62)\), totaling
271. Its scalar map has rank 270, certified by a 270-square minor
with determinant \(167431\) modulo \(1000003\); the exact nonzero
kernel is
\[
                             K_{R^6}=(0,0,-6R^5,0).       \tag{9}
\]
The traceless map has 813 columns and 6677 rows. Its retained
813-square minor has determinant \(131985\) modulo \(1000003\),
so it too is injective.

The known \(h^4\) right side after (8) is traceless and linear in
\(u,v,\gamma,e\). Two retained rational row functionals, each on
629 specified rows, annihilate ALL 813 traceless columns. Their
values on this right side are exactly
\[
                                  \gamma,\qquad e.       \tag{10}
\]
These identities are verified on the constant coefficient and every
parameter coefficient. They imply
\(\gamma=e=0\) over any characteristic-zero field. There is no
reality assumption, parameter sampling, or nonreduced scheme
assertion.

A particular for the \(v\)-coefficient is \(MV(K_{R^4S})\), while
the \(u\)-coefficient is identically zero. Direct substitution verifies
these identities. The complete scalar and traceless kernel results
therefore make the ENTIRE family through \(h^4\)
\[
 \begin{split}
 B&=B^*+uK_{R^3}I_2+vK_{R^4S}I_2,\\
 C&=C^*+vMV(K_{R^4S})+wK_{R^6}I_2.                      \tag{11}
 \end{split}
\]
All \(u,v,w\) remain. No gauge normalization is used to obtain (11).
Its compatibility is verified coefficientwise, not merely at
selected parameter values.

## A constant h5 obstruction on the entire family

At \(m=4\), (5) has scalar component counts \((29,47,62,40)\).
Thus the unknown \(E\) in (2) has \(178\cdot4=712\) columns when
all matrix entries are allowed.

Using the full matrix Moyal commutator with (11), the known \(h^5\)
right side is exactly
\[
                                  R_0+vR_v,             \tag{12}
\]
with 402 and 192 nonzero matrix-polynomial terms, respectively.
All \(u,w\) and quadratic parameter coefficients cancel identically.
They were not removed from the family.

The exact 274-row rational functional already retained for the
positive point satisfies
\[
                 \ell\mathcal L=0
                 \quad\hbox{on ALL 712 columns},\qquad
                 \ell R_0=1,\quad \ell R_v=0.             \tag{13}
\]
The expanded replay rebuilds every full matrix column and every
parameter coefficient in (12). Consequently every member of the
complete family (11) would require \(0=1\). This proves the theorem.

Ordinary matrix trace appears only in the derivation of this
functional after scalar coefficient projection; no cyclic quantum
trace is assumed. The checker verifies (13) directly in the full
matrix coefficient space.

## Corollary for the displayed scaled lower family

There is a further conclusion for a specified lower family at every
scalar \(\lambda\), including zero. It does not assert that this
family exhausts the lower kernels when \(\lambda\ne1\).

Put \(M'=\partial_tM=E_{21}\), \(H=[M,M']=\operatorname{diag}(1,-1)\),
and
\[
 \begin{split}
 g_\lambda&=g_0-\frac{2\lambda^2}{27}K_{R^4a}
                         +uK_{R^3}+vK_{R^4S},\\
 J_{3,\lambda}(f)
  &=\left(\frac{\lambda^3tM}{6}
                     +\frac{\lambda^2H}{12}\right)V^3f
                      +\frac{\lambda M'}4V(\partial_P\partial_Zf).
 \end{split}
\]
The displayed scaled lower coefficients are
\[
 \begin{split}
 F_1&=\lambda MVf,\\
 F_2&=\left(g_\lambda+\frac{\lambda^2t}{2}V^2f\right)I_2,\\
 F_3&=J_{3,\lambda}(f)+\lambda MVg_\lambda+wK_{R^6}I_2.   \tag{14}
 \end{split}
\]
The coefficient at \(K_{R^4a}\) is \(-2\lambda^2/27\), not a fixed
\(-2/27\) independent of the connection scale. All \(u,v,w\) are retained.

Every coefficient in (14) has spatial degree at most 25.
The replay verifies all parameter coefficients of every canonical
commutator through \(h^4\). Nevertheless no \(F_4\) of cap25 extends
this displayed family through \(h^5\), for any characteristic-zero
\(\lambda,u,v,w\), including all complex values.

Indeed the known \(h^5\) right side has exactly the parameter
monomials
\[
                         1,\quad v,\quad\lambda^2,\quad
                                  \lambda^2v,\quad\lambda^4,   \tag{15}
\]
with respective full matrix support sizes \(192,138,300,192,390\).
The witness retains a 139-row rational functional \(\ell_{\rm sc}\)
on the scalar coefficient space, annihilating every one of its
178 derivative columns. Extend it by ordinary coefficient trace:
\[
                    \widetilde\ell(C)
                          =\ell_{\rm sc}(\operatorname{tr}C/2).
\]
This lift has 278 matrix rows. Its values on the five coefficient
vectors in (15) are exactly
\[
                              (1,0,0,0,0).              \tag{16}
\]

The unknown map is now
\[
 \mathcal L_\lambda E
   =\delta_iE_j-\delta_jE_i+
                  \lambda([A_i,E_j]-[A_j,E_i]).
\]
At this coefficient order the matrix entries are commuting
polynomials, so the ordinary trace of each displayed commutator is
zero. Consequently \(\widetilde\ell\mathcal L_\lambda=0\) for ALL
712 columns and every \(\lambda\). The checker also verifies both
affine coefficients of this operator directly. Equation (16) would
therefore require \(0=1\), proving the corollary.

This full comparison matters: a zero of a single nonconstant
quadratic trace test in \(\lambda^2\) is only a necessary condition,
not an actual compatible point. The complete trace conditions
above contain the constant \(1\); no algebraic root survives within
the displayed family (14).

For \(\lambda\ne1\), additional homogeneous matrix directions or
different second and third coefficients have not been classified.
The corollary supplies no universal statement about all jets with
first coefficient \(\lambda MVf\).

## Retained certificates and scope

The [witness](airy_first_coefficient_cap25_witness.json) retains the
two large full-column minors, the scalar third-coefficient minor,
both exact compatibility functionals (10), and the scaled-family trace
functional (16), with precise coefficient locators.
It identifies the frozen positive-point and scalar-kernel inputs.
The [checker](check_airy_first_coefficient_cap25.py) reconstructs all
column maps and parameter coefficients. No stored rank flag is
accepted as a proof.

The two large rank controls require the installed **python-flint**
package, used only for exact modular determinants. All operator,
kernel, parameter and functional identities use Fraction arithmetic.
Absence of the determinant dependency is an error; those checks
are never skipped. Run from the repository root:

    .venv/bin/python research/dixmier/dixmier-formal-lifting/check_airy_first_coefficient_cap25.py

The [receipt](airy_first_coefficient_cap25_controls.json) records
the exact environment, source hashes and coverage. Modular nonzero
minors give rational rank lower bounds; exact independent kernels
give the scalar upper bounds. No modular inconsistency is promoted
alone to an unrestricted characteristic-zero nonexistence theorem.

This strengthens the fixed-point extension calculation to ALL
later coefficients with the prescribed first coefficient (1).
It does not classify other \(h\)-linear matrix data or other
source-dependent connections. The main complete-kernel classification fixes \(\lambda=1\);
the scaled corollary covers only the explicitly displayed later
coefficients (14).
No formal verification, independent whole-proof replay, external
review, bounded formal gauge, scalar realization or DC1/DC2
conclusion is asserted.
