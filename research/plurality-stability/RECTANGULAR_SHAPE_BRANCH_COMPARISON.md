# Comparing a rectangular Young shape with its neighbor

In six cells, let $N_\lambda$ denote the normalized,
ordered six-cell kernel on the shape $\lambda$ minus space,
with diagonal normalization $a_\rho^{-1/2}$ as in (S67).
The two shapes are related by an isometric embedding $J$ of the
$(2,2,2)$ ordered multiplicity space into the $(2,2,1,1)$
space and a positive semidefinite operator $B$:
\[
 J^*N_{(2,2,1,1)}J=\frac23N_{(2,2,2)}+\frac13B.
 \tag{QB1}
\]
Consequently $\|N_{(2,2,2)}\|\le(3/2)\|N_{(2,2,1,1)}\|$.
With the [quadratic gap proof](SIX_CELL_2211_QUADRATIC_GAP.md),
this gives a strict gap in $(2,2,2)$ for $0<\rho\le2/3$.
For $\rho<2/3$ the gap constant is at least
$\ell(1-3\rho/2)/\sqrt{1-\rho^2}$, uniformly in ambient
dimension $d\ge5$. At the endpoint the gap is qualitative.

## The common branch and the coefficient

Use reverse label order $6,5,4,3,2,1$ in the orthonormal Young basis.
The subgroup fixing the original label $1$ is the first five labels
in this order. Shape $(2,2,2)$ restricts to the single branch
$\mu=(2,2,1)$ by removing the last box $(3,2)$.
Shape $(2,2,1,1)$ restricts to $\mu\oplus\nu$, where
$\nu=(2,1,1,1)$, with $\mu$ obtained by removing $(4,1)$.
The Young matrices of the subgroup on the common branch are
identical in the two reverse bases.

In the $\mu$-tableaux the original label $2$ is the last of
the first five labels. It lies in the corner $(2,2)$ or $(3,1)$.
For shape $(2,2,2)$ the last two labels are in the same column
in the first case and in the same row in the second. Thus its
$(12)$-minus basis is the set of three $\mu$-tableaux with
original label $2$ at $(2,2)$. In shape $(2,2,1,1)$ the same
tableau in the $\mu$ branch has original label $1$ at $(4,1)$.
The axial distance for the transposition of the last two reverse
labels is $-3$. In its two-dimensional Young matrix the diagonal
entry at the common $\mu$-tableau is $-1/3$ and the other
tableau belongs to $\nu$. The orthogonal minus projection $P$
has diagonal $2/3$ on this $\mu$-tableau.

Map a selected unit $\mu$-tableau $e$ in $(2,2,2)$ to
$Pe/\sqrt{2/3}$ in the $(2,2,1,1)$ module.
Distinct selected tableaux give disjoint pairs under the last
adjacent transposition, so this map is an isometry. Its two
components are $\sqrt{2/3}\,e$ in $\mu$ and $\sqrt{1/3}$ times
an orthonormal companion in $\nu$. This isometry extends
pointwise to the ordered-chamber $L^2$ spaces. In
[S67](SIMPLEX_INTERFACE_SYMMETRY.md) the operator before the minus
compression is a sum over the subgroup fixing label $1$, hence
block diagonal in $\mu,\nu$. Applying the isometry on both sides
gives (QB1), with $B$ the restriction of the $\nu$ kernel by the
companion isometry. All scalar Gaussian factors and the local
normalization are common.

## Positivity of the subgroup kernel

For clarity, the positivity here is positivity of the integral
operator. Before the minus compression, the group-dependent
part of (S67) is $\sum_{g\in S_5}e^{\tau x\cdot gy}R(g)$,
where the extra factor for the fixed original first label has
been included in the dot product. The other Gaussian factors
are a positive constant times separate scalar factors in $x$
and $y$. Let
\[
 \Phi_\tau(x)=\bigoplus_{j\ge0}\sqrt{\tau^j/j!}\,x^{\otimes j},
 \qquad U(g)=\bigoplus_{j\ge0}g^{\otimes j}.
\]
Then $\langle\Phi_\tau(x),U(g)\Phi_\tau(y)\rangle
=e^{\tau x\cdot gy}$. The average of $U(g)\otimes R(g)$ over
the finite subgroup is its orthogonal invariant projector.
It follows that the sum kernel is a positive semidefinite
Gram kernel on vector data. For compactly supported bounded
fields, its integrated positivity follows by integrating this
Gram representation. Scalar normalization and compression
preserve the inequality, which extends to the full $L^2$
space by boundedness of the Mehler kernels. In particular,
both $B$ and $N_{(2,2,2)}$ are positive semidefinite.

Therefore (QB1) gives the claimed operator norm comparison.
The integrated strict supersolution in (QH7), the positive
folded kernel, and compactness give
$\|N_{(2,2,1,1)}\|<\rho$. Hence
$\|N_{(2,2,2)}\|<3\rho/2$, strictly below one even at
$\rho=2/3$. The local lower bound and the unused-coordinate
Hermite factors give the stated gap for the full cylindrical
block. The [direct cubic estimate](SIX_CELL_222_GAP.md) covers every
positive correlation. The compression gives a separate consequence
of the quadratic estimate.

## All rectangular shapes

**Theorem (Rectangular branch comparison).** Let
\(\lambda=(b^r)\), \(b,r\ge2\), and
\(\eta=(b^{r-1},b-1,1)\), both of size \(k=br\). At each
\(0<\rho<1\) the corresponding normalized ordered operators
satisfy an isometric compression comparison
\[
            c\,N_\lambda\ \preceq\ J^*N_\eta J,
             \qquad c=\frac{b+2}{2(b+1)} .
 \tag{QB2}
\]
In particular \(\|N_\lambda\|\le c^{-1}\|N_\eta\|\). Whenever
\(\|N_\eta\|<c\), this supplies a positive rectangular gap off
the level one, with lower bound
\((\ell/\sqrt{1-\rho^2})(1-c^{-1}\|N_\eta\|)\).

The subgroup fixing the last reverse label acts on the sole
branch \(\mu=(b^{r-1},b-1)\) of \(\lambda\).
In \(\lambda\) the last reverse label is at \((r,b)\).
The penultimate reverse label is at either \((r-1,b)\) or
\((r,b-1)\), with the usual same-column and same-row
signs respectively. The minus tableaux are exactly those
of the first kind. In \(\eta\), the common branch \(\mu\)
appends its last reverse label at \((r+1,1)\). For a selected
tableau the axial distance from \((r-1,b)\) to \((r+1,1)\)
is \(-(b+1)\). The \((12)\)-minus projection on this basis
entry has diagonal \(c=(1+1/(b+1))/2\).
Swapping the last two reverse labels puts the last label
at \((r-1,b)\), hence in the distinct branch with diagram
\((b^{r-2},b-1,b-1,1)\).
This diagram remains valid when \(r=2\) or \(b=2\).
On distinct selected tableaux the two-label pairs are
disjoint, so normalized projection gives the same isometry
as above, with common-branch norm \(\sqrt c\) and
companion-branch norm \(\sqrt{1-c}\).
Other branches of \(\eta\) have no contribution to this
compression. The Gram positivity of each subgroup branch
then gives (QB2).

This comparison uses the complete compressed interaction
in each shape, including every multiplicity. The transfer
condition concerns the normalized interaction norm at the
same \(\rho\), so its scope depends quantitatively on the
neighboring-shape estimate.
