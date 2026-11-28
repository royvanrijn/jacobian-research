# Mixed wall conditions in simplex interface spectra near one

The triple-sign projections governing the fixed-field wall penalty
also govern moving fields in every shape with at most three rows and
in further shapes with a bounded sum of triple-sign projections. The
coercivity argument uses matrix-valued quadratic forms, so it applies
without a positive polynomial field.

Let \(k\ge5\), let \(H=\mathbf1^\perp\cap\{x_1=x_2\}\), and write
\(g_l(x)=w_l\cdot x\) for \(w_l=(e_1+e_2)/2-e_l\),
\(3\le l\le k\). The open full reference facet is
\(F=\{g_l>0:\ 3\le l\le k\}\). Use the standard Gaussian
\(\gamma=\gamma_H\) and let \(G_0\) permute the untied labels.
For an irreducible shape \(\lambda\), write \(V=V_\lambda^-\)
for its \((12)\)-minus space. The fields in the multiplicity
space are the \(G_0\)-equivariant \(V\)-valued functions on \(F\).
Let \(A_l\) be the orthogonal triple-sign projection of (S85)
for labels \(1,2,l\). For \(F_l=\overline F\cap\{g_l=0\}\) the
following is a closed Sobolev form domain:
\[
 {\cal D}_\lambda=
 \left\{h\in W^{1,2}(F,\gamma;V):\
    h\text{ is }G_0\text{-equivariant},\
    A_l\operatorname{Tr}_{F_l}h=0\text{ for each }l
 \right\}.
\]
This trace condition prescribes Dirichlet data in the minus
component \(A_lV\). The unprescribed component \((I-A_l)V\)
has the natural Neumann condition of the gradient form.

**Theorem (Moving mixed-wall limit).** Assume \(\lambda\) is
neither the standard shape nor the rotation hook. Suppose either
that \(\lambda\) has at most three rows or that
\[
                  \left\|\sum_{l=3}^k A_l\right\|<4 .
 \tag{MW1}
\]
On this entire multiplicity space the fixed-index min–max
levels of the interface form, normalized by
\(\rho\sqrt{1-\rho^2}/(2\phi_0(1-\rho))\) with
\(\phi_0=(2\pi)^{-1/2}\), converge as \(\rho\uparrow1\) to
those of
\[
           \int_F (|\nabla h|^2-|h|^2)\,d\gamma,
                   \qquad h\in{\cal D}_\lambda .
 \tag{MW2}
\]
Fixed finite lists of levels are discrete eigenvalues near one.
The lowest limiting level is strictly positive. At fixed
ambient dimension the corresponding cylindrical limit adds each
unused-coordinate Hermite degree, counted with its multiplicity.
Positivity of the active gap near one gives a gap at the same
correlation independent of the ambient dimension.

In particular (MW1) holds for the six-cell \((2,2,1,1)\) block,
with projector-sum norm \(8/3\), while \((2,2,2)\) meets both
conditions with norm \(4/3\). The theorem and the existing
standard and hook estimates give the full balanced interface
gap modulo rotations sufficiently near one in both six and
seven cells.

## Coercivity of the matrix quadratic form

Use the planar Mehler kernel \(m=m_\rho\), the adjacent kernels
\(m_l=m_{\rho,l}=m\exp[-\rho g_l(x)g_l(y)/(1-\rho^2)]\), and the
normalized scalar row
\(b_\rho=m^F+\frac12\sum_l m_l^F\).
The superscript \(F\) means the row integral with respect to
\(\gamma|_F\). These are the quantities in
[(SN2)](STANDARD_SECTOR_NEAR_ONE.md#same-facet-reduction-and-the-two-sided-test).
In particular \(p_F=\gamma(F)\le b_\rho\le1\) and
\(b_\rho(x)\to1\) in the interior. By (S85) of the
[interface reduction](SIMPLEX_INTERFACE_SYMMETRY.md), define
\[
 B(x,y)=I+\tfrac12\sum_l (m_l/m)(x,y)(I-3A_l),\qquad
 D_l(h)=\int_F m_l^F(x)|A_l h(x)|^2\,d\gamma(x).
\]
For any symmetric matrix-valued kernel \(K\), write \(J_K(h)\)
for one half of the integral over \(F^2\) of
\(\langle h(x)-h(y),K(x,y)(h(x)-h(y))\rangle\).
When \(K\) is scalar use the identity on \(V\).
The rescaled form is
\[
 \frac{Q_\rho(h)}{1-\rho}-\int_F b_\rho|h|^2\,d\gamma,\qquad
 Q_\rho(h)=J_{mB}(h)+\tfrac32\sum_l D_l(h).
 \tag{MW3}
\]
This follows by expanding the difference form. The matrix row
of \(mB\) is \(b_\rho I-\tfrac32\sum_l m_l^F A_l\).

### Shapes with at most three rows

For these shapes \(B(x,y)\) is positive semidefinite.
Since \(m_l/m\in[0,1]\) and \(B\) is affine in these numbers,
it suffices to consider the vertex with indicators of a subset
\(L\) of the untied labels, of size \(s\). The class-sum identity
underlying (S89) gives on the minus space
\[
 3\sum_{l\in L}A_l=(s-1)I-Z_{\{1,2\}\cup L}+Z_L ,
\]
where \(Z\) denotes the sum of transpositions on the indicated
labels. Restricting to the subgroup on \(\{1,2\}\cup L\)
gives diagrams \(\beta\) with at most three rows, by the Young
branching rule. On the \((12)\)-minus space the spectator type
\(\mu\) removes a vertical two-strip from \(\beta\). The two class
sums act by the sums \(c_\beta,c_\mu\) of box contents. If the
removed boxes are in rows \(i<j\), then
\(c_\beta-c_\mu=\beta_i-i+\beta_j-j\ge-3\). It follows that the
vertex eigenvalues of \(B\) are
\((3+c_\beta-c_\mu)/2\ge0\). The empty subset gives \(B=I\).

Since \(I\preceq B+\tfrac32\sum_l(m_l/m)A_l\) and
\(J_{m_l}(A_l h)\le2D_l(h)\), (MW3) gives
\(J_m(h)\le2Q_\rho(h)\) and
\(\sum_l D_l(h)\le\tfrac23Q_\rho(h)\).
On each compact interior pair set \(B\to I\) uniformly.

### The projector-sum condition

Even if \(B\) is indefinite, let \(M=\|\sum A_l\|<4\) and
choose \(1/4<\theta<1\) with \(\theta M<1\). Expanding (MW3) in
the orthogonal parts of each \(A_l\) and using
\(J_{m_l}(A_lh)\le2D_l(h)\) gives
\[
 Q_\rho(h)\ \ge\
 J_{mR}(h)+(2\theta-\tfrac12)\sum_l D_l(h),\qquad
 R=I-\theta\sum_l(m_l/m)A_l .
 \tag{MW4}
\]
In detail the \(l\)-term in \(J_{mB}\) is
\(\tfrac12J_{m_l}((I-A_l)h)-J_{m_l}(A_l h)\).
The bound keeps the part with coefficient \(-\theta\) and
absorbs the remaining negative term into \(\tfrac32D_l\).
Here \(R\succeq(1-\theta M)I\) on all pairs, and \(R\to I\)
uniformly on compact interior pair sets.

In both cases \(Q_\rho\) controls \(J_m+\sum_lD_l\) with
a positive constant independent of \(\rho\). For the sharp
interior lower limit, the same estimates give a lower bound
by a nonnegative pair form with fiber matrix tending to \(I\)
on compact subsets of \(F^2\).

## Strong limits and boundary traces

For \(t_i=1-\rho_i\downarrow0\), bounded norm and bounded above
rescaled energy imply \(Q_{\rho_i}(h_i)=O(t_i)\), since \(b_\rho
\le1\). The confinement and compactness arguments
(SN4)–(SN8) apply componentwise using \(J_{m_{\rho_i}}(h_i)/t_i
=O(1)\). Thus a subsequence converges strongly in \(L^2\)
to an equivariant \(h\in W^{1,2}(F,\gamma;V)\). Inserting the
compact interior lower bound from (MW3) or (MW4) into the
finite-difference proof (SN7)–(SN8) gives
\[
              \liminf_i Q_{\rho_i}(h_i)/t_i
                    \ \ge\ \int_F|\nabla h|^2\,d\gamma .
\]
The potential term converges to \(-\|h\|^2\) by boundedness of
\(b_\rho\), interior convergence and strong \(L^2\) convergence.

Fix a point of the relative interior of \(F_l\) and a smooth
compactly supported scalar cutoff \(\chi\) near that point.
The support can be chosen so that \(g_l\) is the smallest
gap on its intersection with \(F\), and \(\chi=1\) in a
smaller neighborhood of the point. The adjacent row estimate
proved in [(SG5)](NONSTANDARD_SECTOR_NEAR_ONE.md#pointwise-comparison-and-loss-of-row-mass)
gives \(m^{H\setminus F}\le C_k m_l^F\) there. In that proof
the probability of crossing any of the walls is bounded by
the sum of their one-dimensional tails. The adjacent row for
the minimal gap dominates the largest tail by conditioning on
the normal of this wall.

Extend \(\chi A_l h_i\) by zero from \(F\) to \(H\).
Its full Mehler difference energy is bounded by a constant
times
\[
 J_m(h_i)+t_i\|h_i\|^2+D_l(h_i).
\]
Indeed, in \(F^2\) split the difference of the product into
the difference of \(h_i\) and the difference of \(\chi\).
The Mehler row of \(|\chi(x)-\chi(y)|^2\) is \(O(t_i)\)
uniformly in \(x\). For \(x\) in a bounded set this follows
from the Lipschitz bound and the Gaussian second moment of
\(y-x\). For \(x\) outside a larger ball it follows from
exponentially small row mass on the cutoff support. The
pairs across the boundary are controlled by the adjacent
row comparison and \(D_l\). The full-space Hermite argument
of (SG7) shows that the strongly convergent zero extensions
have a \(W^{1,2}(H,\gamma)\) limit. This limit is zero outside
\(F\), so its trace vanishes on the flat face. Therefore
\(A_l\operatorname{Tr}_{F_l}h=0\) on the smaller patch.
The relative interior covers the wall up to a surface-null
set. Repeating for all \(l\) gives \(h\in{\cal D}_\lambda\).

## Recovery and min–max values

Smooth compactly supported restrictions satisfying the
projected trace constraints are dense in \({\cal D}_\lambda\).
First truncate the support radially and truncate the vector
norm, preserving trace constraints in the linear subspaces
\(\ker A_l\). For a bounded compactly supported field, remove
neighborhoods of each pairwise wall intersection by a
logarithmic cutoff. In distance \(r\) to the codimension-two
linear intersection take the cutoff zero for \(r\le\epsilon^2\),
one for \(r\ge\epsilon\), and logarithmically linear between
those radii. Its squared gradient has integral \(O(1/|\log
\epsilon|)\) on bounded sets by integrating in the two
orthogonal coordinates. Hence these cutoffs converge in
the Gaussian Sobolev form norm after multiplication by the
bounded field. Away from the pairwise intersections, use a
partition of unity into interior sets and sets meeting a
single wall. Extend the \(A_l\) component oddly through that
wall, and the other component evenly, then mollify with an
even smooth kernel. The odd extension is Sobolev by the zero
trace, and the mollification preserves this trace. Choose the
mollification radius small enough that the support stays away
from other walls. Interior pieces use ordinary mollification.
Finally average over the joint \(G_0\) action on facet and
fiber, which permutes the projected constraints.

The fixed-field calculation (JP2) in the
[nonstandard sector note](NONSTANDARD_SECTOR_NEAR_ONE.md#the-leading-junction-penalty-for-fixed-interface-fields)
gives recovery on this dense space. The embedding into
\(L^2(F,\gamma;V)\) is compact. For example, on bounded
truncations use ordinary Sobolev compactness and for the
tails use the Gaussian energy of the radial cutoff.
The latter can be seen by \(u=\sqrt{\phi_H}h\), for which
\(\int_F|\nabla h|^2\,d\gamma
=\int_F(|\nabla u|^2+(|x|^2/4-n/2)|u|^2)\,dx\)
on smooth compactly supported restrictions, where
\(n=\dim H\). There is no flux from \(x\cdot\nu\) on the
conical faces. Approximation and lower semicontinuity
give the tail bound on the closed form domain.

Use a fixed smooth approximating plane for the min–max upper
limit, and along a lower-limit sequence take orthonormal
bases of almost minimizing planes. The compactness and lower
bound apply to their finitely many elements and all fixed
linear combinations. This proves convergence of the min–max
values. At fixed \(\rho\) the matrix interaction is
Hilbert–Schmidt, dominated in norm by a constant times the
full Mehler kernel. The rescaled diagonal is
\(\rho b_\rho/(1-\rho)\ge\rho p_F/(1-\rho)\). The bounded
fixed-index levels therefore lie below the essential spectrum
near one.

In extra Hermite degree \(r\), the normalized form is
\[
 \frac{\rho^r Q_\rho(h)}{1-\rho}
     +\left(\frac{1-\rho^r}{1-\rho}-1\right)\int_F b_\rho|h|^2\,d\gamma.
\]
For fixed \(r\) the limit shifts by \(r\). The nonnegativity of
\(Q_\rho\) and \(b_\rho\ge p_F\) bound all degrees \(r\ge L\ge1\)
below by \(((1-\rho^L)/(1-\rho)-1)p_F\), which tends to
\((L-1)p_F\). Thus for each fixed index only finitely many
inactive degrees matter. Positivity in every ambient dimension
also follows directly from tensorization of the active gap
as in (S42).

## The positive first level and the six- and seven-cell consequences

A trivial \(G_0\) fiber type can occur only in the standard
and rotation shapes, by adding a vertical two-strip to the
one-row spectator diagram. Thus for all shapes of this theorem
the Gaussian mean of an equivariant field vanishes in every
component. The unrestricted Gaussian Neumann Rayleigh value on
the full \(F\) is bounded below by \(\min(2,\beta_F)>1\), where
\(\beta_F(\beta_F+n-2)\) is its first nonconstant spherical
Neumann value. The strict angular inequality above \(n-1\)
follows from the pointed-cone case of Andrews–Clutterbuck–Hauer,
Theorem 9.1, as checked in the
[standard-sector proof](STANDARD_SECTOR_NEAR_ONE.md#the-invariant-gaussian-neumann-gap).
This applies componentwise on the smaller mixed form domain,
so the first level of (MW2) is positive.

For a shape with \(r\) rows, the content of a vertical two-strip
is at least \(3-2r\). Taking all \(k-2\) untied labels in the
class-sum identity yields
\(\|\sum_l A_l\|\le (k+2r-6)/3\).
For \((2,2,2)\) the only spectator type is \((2,1,1)\),
giving \(4/3\). For \((2,2,1,1)\) the types are
\((2,2),(2,1,1),(1^4)\), giving \(8/3,2,2/3\).
A nonhook shape has at most \(k-2\) rows, so for \(k=6,7\)
its projector sum has norm at most \(8/3,11/3\), respectively.
Taking the minimum over the finite collection of active shapes
and using the already proved standard and hook bounds gives
the asserted full gaps near one. The cutoff remains qualitative.

The coercivity, mixed trace proof and min–max argument have
received internal mathematical review.
