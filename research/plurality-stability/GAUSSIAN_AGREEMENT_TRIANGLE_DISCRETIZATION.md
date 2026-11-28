# Finite approximation of the Gaussian agreement triangle relaxation

Fix \(0<\rho<1\) and let \(\gamma=\gamma_2\) be standard planar
Gaussian measure. Write \(M_\rho\) for the relative Mehler kernel,
the density of a correlated standard Gaussian pair against
\(\gamma\otimes\gamma\). The objective is
\(\int M_\rho F\,d\gamma^2\). The continuum relaxation asks for a
measurable \(F\) with \(0\le F\le1\), symmetry, row mean \(1/3\),
and
\[
                 F(x,y)+F(y,z)-F(x,z)\le1 .
\]
All conditions hold almost everywhere for the relevant product
measure. One can require simultaneous rotation invariance. Indeed,
Haar averaging over rotations preserves every constraint and the
objective, by invariance of the Gaussian pair law.
[The agreement triangle counterexample](GAUSSIAN_AGREEMENT_TRIANGLE_RELAXATION.md)

## Integrated finite values converge to the continuum value

Let \(\mathcal P\) be a finite measurable partition into cells
\(B_a\) of Gaussian masses \(w_a>0\). The integrated finite LP is
\[
 \begin{aligned}
 &\text{maximize}\quad \sum_{a,b} W_{ab}c_{ab},
  \qquad W_{ab}=\int_{B_a\times B_b}M_\rho\,d\gamma^2,\\
 &0\le c_{ab}\le1,\qquad c_{ab}=c_{ba},\qquad
      \sum_b w_bc_{ab}=\tfrac13,\qquad
      c_{ab}+c_{bc}-c_{ac}\le1 .
 \end{aligned}                                               \tag{AD1}
\]
The triangle rows include repeated indices. In particular the
\(a,b,a\) row is \(2c_{ab}-c_{aa}\le1\). The diagonal entry
represents independent points in the same cell, so the raw LP
uses the box \(0\le c_{aa}\le1\). If the grid admits a finite
group of rotations, one may impose invariance under its action
on cells. Averaging a feasible matrix under the group preserves
the objective and gives that invariance.

Every feasible matrix lifts to the corresponding block-constant
continuum kernel with the same objective. It can subsequently be
averaged under all rotations at no cost. Conversely, for a feasible
continuum \(F\), put
\[
 c_{ab}=\frac{1}{w_aw_b}\int_{B_a\times B_b}F\,d\gamma^2 .
\]
The row condition follows by averaging. For each triangle row take
three independent points drawn from the respective normalized cell
measures. Their expected triangle inequality gives the finite row
even if two or three of the cells coincide. A rotationally invariant
\(F\) gives a matrix invariant under any rotations of the grid.

Put \(M_{\mathcal P}=\mathbb E[M_\rho\mid
\mathcal P\otimes\mathcal P]\). The objective of the projected
matrix is \(\int M_{\mathcal P}F\,d\gamma^2\), so its difference
from the continuum objective has absolute value at most
\(\|M_\rho-M_{\mathcal P}\|_{L^1(\gamma^2)}\). In particular, if
\(V_{\mathcal P}\) and \(V\) are the raw finite and continuum suprema,
\[
          0\le V-V_{\mathcal P}
                  \le \|M_\rho-M_{\mathcal P}\|_{L^1(\gamma^2)}.
                                                               \tag{AD2}
\]
For any refining finite partitions generating the measurable sets
modulo null sets, the right side tends to zero by conditional
expectation convergence. The finite values are nondecreasing under
refinement, by representing the coarse kernel on the finer cells.
The Hilbert–Schmidt error gives the same conclusion, since
\(\|M_\rho\|_2^2=(1-\rho^2)^{-2}\) by the planar Hermite
expansion and \(\|F\|_2\le1\).

Bare operator norm of the kernel error is insufficient without
a trace norm bound on \(F\). For example, on \(4n\) equal-mass
cells take a symmetric Sylvester Hadamard matrix \(W\) of order
\(n\), with \(n\) a power of two, and set
\[
     G=\begin{pmatrix}0&W\\W&0\end{pmatrix},
        \qquad H=\begin{pmatrix}G&-G\\-G&G\end{pmatrix}.
\]
The kernel with entries \(H\) has zero rows and diagonal, entries
of absolute value at most one, and operator norm \(1/(2\sqrt n)\),
since the matrix is divided by \(4n\) as an integral operator.
The kernel \(F=1/3+H/6\) meets every box, row, and triangle
constraint, including diagonal lower bound \(1/3\). Yet
\(\int HF\,d\gamma^2=1/12\), because half the entries of \(H\)
have absolute value one. This example can use radial cells only.
It rules out a uniform bound on the pairing with general errors
from these constraints and operator norm alone. For ordinary
agreement kernels, including their rotation averages, the operator norm bound is available with coefficient \(2/3\).
Subtract the constant kernel \(1/3\) to express the remainder
as the sum of the three rank-one kernels of the centered
indicators. This is positive semidefinite with trace \(2/3\),
and \(T_\rho-P_\mathcal P T_\rho P_\mathcal P\) annihilates constants,
where \(P_\mathcal P\) is the conditional-expectation projection.
Averaging the ordinary agreement kernel over rotations
preserves these properties.

## The tightened diagonal gives the same limiting supremum

Ordinary agreement projected to a cell has diagonal entry
\(\sum_i \gamma(A_i\mid B_a)^2\ge1/3\). This remains true after
rotation averaging. The lower bound is independent of the raw
LP constraints. For example, with two cells of masses \(1/4,3/4\),
the entries \(c_{11}=0,\ c_{12}=4/9,\ c_{22}=8/27\) satisfy
all raw rows and triangles. Nevertheless imposing
\(1/3\le c_{aa}\le1\) has vanishing cost under refinement.

Here is an explicit repair of any raw feasible matrix. Let
\(m=\max_a w_a\le1/24\) and first mix the matrix with the
constant \(1/3\) by weight \(\epsilon=12m\), obtaining
\(\widetilde c\). Each box endpoint then has margin at least
\(\epsilon/3\), and each triangle has slack at least
\(2\epsilon/3\). Raise each diagonal to
\(\max(\widetilde c_{aa},1/3)\). Its extra row contribution is
\(\delta_a=w_a(1/3-\widetilde c_{aa})_+\), between zero and \(m/3\).
Subtract \(u_a+u_b\) on every off-diagonal position, with
\[
 \overline u=
 \frac{\sum_a w_a\delta_a/(1-2w_a)}
      {1+\sum_a w_a/(1-2w_a)},\qquad
 u_a=\frac{\delta_a-\overline u}{1-2w_a}.
\]
These equations give \(\overline u=\sum_a w_au_a\). The row
subtraction is \((1-2w_a)u_a+\overline u=\delta_a\), as
required. Also \(0\le\overline u\le m/3\) and
\(|u_a|\le2m/3\), so the correction of an off-diagonal
entry is at most \(4m/3\) in absolute value. The box and
all-distinct triangle slack cover these changes. In a repeated
row of the form \(2c_{ab}-c_{aa}\le1\) the raised diagonal
helps and the off-diagonal correction changes the left side
by at most \(8m/3\). Other repeated rows just bound a diagonal
by one. The repair preserves every grid symmetry.

The off-diagonal change from the original matrix is at most
\(10m\), while all diagonal blocks together have product
Gaussian measure at most \(m\). Since \(\int M_\rho\,d\gamma^2=1\),
the objective changes by at most
\(10m+\|M_\rho\|_{L^2(\gamma^2)}\sqrt m\). Thus (AD2) gives
the same limiting supremum for the tightened LP along a
refinement with maximum cell mass tending to zero. Tightened
finite values are also nondecreasing, since each finer diagonal
inside a coarse cell inherits that cell's diagonal value.

In particular, a strict continuum counterexample to the simplex
gives a strict finite counterexample against the **true Gaussian
simplex value** at a sufficiently fine level, for both diagonal
choices. Project its feasible kernel and make the two errors just
bounded smaller than its strict excess over the simplex. For the
radial-cut witness the diagonal repair is unnecessary on a fine
polar grid. Two points in the same angular bin have circular
distance at most the bin width \(\Delta\), so the witness is at
least \(1-\Delta/a_\epsilon-\epsilon>1/3\) there for sufficiently
small \(\Delta\). An optimal simplex projection on selected
coarse grids does not change this refining-limit conclusion.

## A convergent hierarchy using finite label consistency

A stronger condition can retain all pairwise probabilities of a
three-label assignment on a fixed number of sampled points.
For \(m\) points, let \(\mathcal C_m\) be the convex hull of the
matrices \((1_{\ell_i=\ell_j})_{i,j=1}^m\) over assignments
\(\ell\in\{1,2,3\}^m\). The level-\(m\) constraint on \(F\) is that
for almost every \(m\)-tuple \((x_1,\ldots,x_m)\) its off-diagonal
entries \(F(x_i,x_j)\), with diagonal filled by ones, form a
matrix in \(\mathcal C_m\). Keep the range, symmetry and exact
row mean in (AT1). This constraint includes the corresponding
constraints on smaller tuples, by taking principal submatrices.
It is preserved under simultaneous Gaussian rotations and
averaging. For ordinary partitions the matrices are vertices of
\(\mathcal C_m\) or, after averaging over transformations of
the partition, convex combinations of vertices.

Let \(V_m(\rho)\) be the supremum under this additional
consistency constraint and \(V_{\mathrm{bal}}(\rho)\) the
ordinary balanced optimum. For any finite Gaussian partition
\(\mathcal P\) with \(m\) cells of masses \(w_a\) one has
\[
  0\le V_m(\rho)-V_{\mathrm{bal}}(\rho)
     \le \|M_\rho-M_{\mathcal P}\|_{L^1(\gamma^2)}
                                      +\sqrt{3\sum_a w_a^2}.
                                                               \tag{AD3}
\]
In particular along successively finer generating partitions
with maximum cell mass tending to zero, the finite consistency
levels converge to the genuine optimum.

For the proof draw one independent point in each cell from
its normalized Gaussian law. The average sampled matrix
\(C\) belongs to \(\mathcal C_m\) by convexity. Its off-diagonal
entries are the block averages of \(F\). Its diagonal entries
are one. A finite convex decomposition of \(C\) corresponds
to a random assignment of all the cells to the three labels.
If the resulting mass vector of a realization is \(p\),
then
\[
  \mathbb E\|p-\mathbf1/3\|_2^2
    =\sum_{a,b}w_aw_b C_{ab}-\tfrac13
    \le \sum_a w_a^2 ,
\]
using the exact row mean of \(F\). The ordinary but possibly
unbalanced cell partition of a realization can be balanced
by changing labels on a set of Gaussian mass
\(\frac12\|p-\mathbf1/3\|_1\). This is done by reassigning
surplus mass from labels with more than \(1/3\) to the
deficient labels, using nonatomicity. At the true Mehler
objective the absolute change is at most twice this mass,
since each coordinate has marginal \(\gamma\).
The average change is at most \(\sqrt{3\sum_a w_a^2}\)
by Cauchy–Schwarz. The integrated objective of \(C\) is at
least that of the block averages of \(F\), by nonnegativity
of the diagonal weights. The expected agreement value of
the cell assignments is exactly the integrated objective
of \(C\). Replacing the original \(F\) by its block averages
costs at most \(\|M_\rho-M_{\mathcal P}\|_1\). These estimates
give (AD3).

This gives consistency levels whose limit encodes the
balanced optimization problem. In contrast to (AD1), their
constraints retain the joint assignment of all labels in a
finite sample. This convergence does not give an upper
certificate at any fixed level.

## Recovering a partition from its averaged agreement

If a rotationally averaged ordinary balanced agreement kernel
equals the simplex kernel \(F_T\) almost everywhere, then the
partition is a simplex up to rotation, relabeling and null sets.
To see this, let \(p(r)\) be the circular label mass vector and
let \(\nu\) be Rayleigh measure. Integration in the angular
difference gives
\(p(r)\cdot p(s)=1/3\) for \(\nu^2\)-almost every \((r,s)\).
By balance, \(\int p\,d\nu=(1/3,1/3,1/3)\). Integrating the
centered Gram identity on \(E\times E\) for every measurable
\(E\) shows \(p(r)=(1/3,1/3,1/3)\) almost everywhere.

Write \(a_i(r)=\int 1_{A_i}(r,\theta)e^{\mathrm i\theta}
\,d\theta/(2\pi)\) for the first circular moment. The first
Fourier coefficient of the averaged agreement is the Gram
kernel \(a(r)^*a(s)\). For \(F_T\) it is the constant
\(9/(4\pi^2)\). Integrate first on the whole square and then
in each variable, and subtract the resulting mean vector.
The centered vector has Gram kernel zero almost everywhere.
Integrating on \(E\times E\) again shows that \(a(r)\) is a
common vector, with squared norm \(9/(4\pi^2)\).
On every slice each individual moment has norm at most
\(\sin(\pi/3)/\pi\), by the one-set cosine rearrangement.
The total squared norm saturates all three caps. Equality
forces each cell on the circle to be the \(120^\circ\) arc
centered at its moment direction, since a cosine superlevel
arc uniquely maximizes its integral at its prescribed mass.
These three directions are shared at almost all radii, and
the three arcs tile the circle. This gives the asserted simplex.
[Circular moment caps and equality](THREE_CELL_SIMULTANEOUS_FIRST_MOMENT_ENVELOPE.md#feasible-moments-on-a-single-circle)

This implication requires equality of the kernels. A finite
objective comparison alone yields neither that equality nor
uniqueness. The strict relaxation counterexample has no
implication for the optimality of ordinary Gaussian partitions
or for the balanced low-influence voting conjecture.
