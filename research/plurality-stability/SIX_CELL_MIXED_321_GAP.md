# The six-cell \((3,2,1)\) interface block

For the balanced six-cell Gaussian simplex and each fixed \(0<\rho<1\),
the full \(V_{(3,2,1)}\) interface isotypic component has a positive
\(L^2\) spectral gap in every ambient dimension \(d\ge5\), with a
constant independent of \(d\). This controls every multiplicity of
this shape. The [cubic](SIX_CELL_222_GAP.md) and
[quadratic](SIX_CELL_2211_QUADRATIC_GAP.md) proofs give the
remaining nonstandard six-cell estimates at every positive
correlation. The balanced standard sector is the remaining
all-correlation obligation in six cells.

Use the facet conventions of (S67)–(S75) in
[the interface reduction](SIMPLEX_INTERFACE_SYMMETRY.md). In the
positive off-diagonal seminormal Young basis for \((3,2,1)\), retain
the eight coordinates in the \((12)\)-minus space. Order tableaux by
recursive bottom-up induction on the number of boxes, appending the
largest label in each removable corner in order from the top row to
the bottom. The rational basis is specified as in (S74), with diagonal
entry \(1/r\) at axial distance \(r\), off-diagonal entry \(1\) towards
the greater inversion count of the row-reading word, and opposite
entry \(1-r^{-2}\). Write \(R\) for these matrices and \(P\) for the
minus projection. All representations below are in this basis.

The positive equivariant field on the ordered facet \(W\) is \(Lx\).
Its first two columns vanish. The last four columns are

\[
 \begin{pmatrix}
  0&1&-1&0\\
  1&-1/3&-2/3&0\\
  1/2&-1/2&0&0\\
  0&1/2&1/2&-1\\
  1/2&-1/6&1/3&-2/3\\
  3/8&3/8&-1/4&-1/2\\
  3/8&-3/8&0&0\\
  3/16&3/16&-3/8&0
 \end{pmatrix}.
\]
One way to obtain it is to sum the \(\{4,5,6\}\)-permutation orbit of
the all-ones vector in the minus space, then send label \(3\) to each
\(l\in\{3,4,5,6\}\) to form column \(l\) and multiply by \(16/169\).
The checker reconstructs \(L\), checks its remaining-label
equivariance, and checks positivity on \(W^\circ\) by its values on
the four rays \(r_m\), \(m=2,3,4,5\). Here \(r_m\) has first \(m\)
entries \(6-m\) and all other entries \(-m\). Every ray value is
coordinatewise nonnegative, and its sum over the four rays is
coordinatewise strictly positive.

Put \(T_l=2PR((2l))P\) on the minus space. The compressed interaction
has the raw form
\(\mathsf K(x,y)=2M_0(x,y)I+\sum_{l=3}^6M_l(x,y)T_l\).
The folded kernel is entrywise nonnegative by (S69).
The normalization is self-adjoint after a positive diagonal rescaling,
which preserves positivity and the inequalities below.

Let \(C_l=I-2w_l w_l^T\) be the rotation-hook transport of (S71) and
\(m_l(x)=\int_{F_{12}}M_l(x,y)y\,d\nu_{12}(y)\).
The rotation identity gives
\[
   a_\rho Lx-\mathsf K(L\,\cdot)(x)
       =\sum_{l=3}^6 (LC_l-T_lL)m_l(x).                            \tag{SM1}
\]
For \(g=(2l)\pi\in S_{\{2,3,4,5,6\}}\) with
\(\pi\in S_{\{3,4,5,6\}}\), put
\[
 c_{g,m}=(LC_l-T_lL)\pi r_m,\qquad
 \Phi_m(A)=\sum_g c_{g,m}\prod_{i=2}^6 A_{i,g^{-1}(i)}.
\]
Assign \(c_{g,m}=0\) if \(g\) fixes \(2\). Every other permutation
decomposes uniquely in the stated way. Since
\(y=\sum_{m=2}^5(y_m-y_{m+1})r_m/6\) on the ordered facet, folding
(SM1) over the 24 remaining-label chambers writes its integrand as
a positive common Gaussian factor times a positive combination
of \(\Phi_m([e^{\tau x_i y_j}]_{i,j=2}^6)\). The fixed coordinate 1
contributes only to this common factor and \(\tau=\rho/(1-\rho^2)\).
These five-by-five matrices are strictly totally positive in the
interior by (S47).

### Finite positive-coefficient calculation

The two Jacobi words are \(1213214321\) for the lower factor and its
reverse for the upper factor, with one independent positive variable
at every position. Five positive diagonal variables give a
25-parameter product \(A=L_{\rm J}DU_{\rm J}\). The checker computes
every coefficient of every component of \(\Phi_m\) under this
substitution. It finds 744735 nonzero coefficients across the 32
polynomials, all strictly positive, with every polynomial nonzero.
The receipt records counts, minimum coefficients and hashes of the
coefficient arrays.

Here is a direct covering proof for that substitution. For a strictly
totally positive matrix \(A\), Gaussian elimination yields \(A=BDU\)
with unit triangular \(B,U\) and positive diagonal \(D\), since the
leading principal minors are positive. For \(B\), every flag minor
\(\det B[I,1:r]\), \(|I|=r\), is positive, since it equals the
corresponding flag minor of \(A\) divided by the first \(r\) diagonal
pivots. Let \(n\) denote the size of \(B\), put \(b_j=B_{nj}\), and
set \(a_j=b_j/b_{j+1}>0\), with \(b_n=1\). For
\(x_i(t)=I+tE_{i+1,i}\) put
\[
  Q=x_{n-1}(a_{n-1})\cdots x_1(a_1),
  \qquad Q^{-1}=I-\sum_{j=1}^{n-1}a_jE_{j+1,j}.
\]
The last row of \(BQ^{-1}\) is \(e_n^T\), since column \(j<n\)
is replaced by column \(j\) minus \(a_j\) times column \(j+1\).
Its leading block \(B'\) remains unit lower triangular. For any flag
minor of \(B'\), apply the same column differences simultaneously
to the matrix with rows \(I\cup\{n\}\) and columns \(1,\ldots,r+1\).
The determinant is unchanged and the last row is zero on its first
\(r\) entries. Expanding on that row gives
\[
  \det B'[I,1:r]=\frac{\det B[I\cup\{n\},1:r+1]}{b_{r+1}}>0.
\]
Induction yields the concatenated lower word \(1,21,321,4321\).
Applying the same reduction to \(A^T\) gives the reverse word for the
upper factor. This proves coverage. It is also the specific lower
and upper word case of Fomin and Zelevinsky's general factorization
theorem. [Fomin–Zelevinsky, *Total positivity: tests and
parametrizations*, Definition 22 and Theorem 23](https://www.math.kth.se/math/GRU/2011.2012/SF2741/FZ.pdf)

Thus each coordinate of (SM1) is strictly positive on the facet
interior. For \(T=a_\rho^{-1/2}\mathsf K a_\rho^{-1/2}\), the folded
operator, the positive field \(\psi=\sqrt{a_\rho}\,Lx\) satisfies
\(T\psi<\psi\) coordinatewise. If \(Th=\lambda h\) with \(h\ne0\),
then \(|\lambda||h|\le T|h|\). Self-adjointness gives
\[
 |\lambda|\langle\psi,|h|\rangle
       \le\langle T\psi,|h|\rangle
       <\langle\psi,|h|\rangle.
\]
All nonzero eigenvalues have magnitude less than one. Each
\(M_l\) on the reference facet is dominated by \(M_0\), a scalar
multiple of the restricted Mehler kernel. The unweighted kernels
are Hilbert–Schmidt, and the lower bound
\(a_\rho/\rho\ge\ell/\sqrt{1-\rho^2}\) preserves this after
normalization. Finite folding also preserves compactness. Thus
\(\|T\|<1\), giving the claimed positive active-space gap.
In the Hermite decomposition underlying (S42), unused-coordinate
degree \(m\) multiplies \(T\) by \(\rho^m\) and preserves the bound.

The checker uses exact rational Young matrices and integer monomial
accumulation after clearing denominators. Gaussian positivity,
the matrix-factor coverage, and the operator implication are
written arguments. Run the [checker](check_six_cell_321_gap.py)
from the repository root with standard-library Python:

~~~sh
python3 research/plurality-stability/check_six_cell_321_gap.py
~~~

The [verification record](six_cell_321_gap_verification.json)
contains the complete coefficient counts and digests. The proof
and exact replay have received internal review.
