# The adaptive Gaussian surface problem with unrestricted direction rank

Allowing Hilbert-valued directions gives an adaptive surface
problem with a unique scalar diagonal certificate at each
finite-energy partition. The certificate also controls
how the weights of near-attaining dual tests can drift.
For actual three-cell noise-stability maximizers the
unrestricted surface problem still gives exactly the
correlation derivative. Thus a global simplex minimum
for this auxiliary problem would still exclude distinct
right-contact ties.

The Gram-operator formulation and its scalar diagonal have
the finite vector-program analogue in the Max-Cut
semidefinite relaxation. Here compactness of the trace
variables and the fresh-coordinate variation supply
attainment and dual equality on continuous interface measures.
[Goemans and Williamson, Sections 4.1–4.2](https://math.mit.edu/~goemans/PAPERS/maxcut-jacm.pdf)

## Hilbert directions and attainment

Fix \(0<r<1\) and an ordinary planar partition \(A\) of
finite Gaussian perimeter. Write \(\nu_e\) for the
Gaussian measure on the tagged interface \(e=(i,j)\)
and \(d_e=e_i-e_j\) for its label vector. Use the tagged
measure space \(\mathcal E\) with measure \(\nu=\sum_e\nu_e\);
the copies of distinct tags are disjoint in \(\mathcal E\).
Put \(s=\sqrt r\). For a measurable field \(u:\mathcal E\to
\ell^2\) bounded in norm by one, set
\[
 Z(u)=T_s\sum_e d_e\otimes u_e\,\nu_e,\qquad
 F_\infty(A)=\sup_{|u|\le1}\|Z(u)\|_2 .
\]
Tensor norms sum the three label rows and the Hilbert
direction coordinates. The smoothing on measures uses
the relative Mehler kernel. The scalar tagged kernel is
\[
 K((e,x),(f,y))=(d_e\cdot d_f)M_r(x,y).
\]
For bounded fields of any rank, its quadratic pairing
is exactly the squared norm of their smoothed tensor.

Assume the original two-direction functional is finite.
Taking a direction supported on one tag shows that
\(T_s\nu_e\in L^2\) for every \(e\). The smoothed norm
pointwise is at most \(\sqrt2\sum_e T_s\nu_e\).
Consequently \(F_\infty(A)<\infty\), and the kernel \(K\)
is absolutely integrable on \(\mathcal E^2\). Truncating
the direction coordinates shows that \(F_\infty\)
is the supremum of the finite-rank values.
[The original adaptive surface functional and dual tests](THREE_CELL_ADAPTIVE_SURFACE_VARIATIONAL_PROBLEM.md#a-dual-formula-with-metric-surface-weights)

The supremum is attained. To see this, on the real
space \(L^2(\nu)\) represent a direction field by the
positive trace-class operator
\[
             X_u=\sum_n |u_n\rangle\langle u_n| .
\]
Its trace is at most \(\nu(\mathcal E)\), and for every
bounded \(\phi\ge0\),
\[
       \operatorname{Tr}(X_u M_\phi)
                       \le\int \phi\,d\nu ,
\]
where \(M_\phi\) multiplies by \(\phi\). Conversely a
positive trace-class operator with these inequalities
gives a direction field by its spectral decomposition:
take coordinates \(u_n=\sqrt{a_n}\,h_n\) for eigenvalues
\(a_n\) and orthonormal eigenfunctions \(h_n\). The
inequalities give \(\sum_n|u_n|^2\le1\) almost everywhere.
These operators form a compact set for the topology of
testing against compact operators. Indeed the trace ball
is weak-star compact as the dual ball of the compact
operators. Positivity is closed. For a fixed \(\phi\ge0\),
the left side of the inequality is the supremum of
the finite-rank tests of \(M_\phi^{1/2}X M_\phi^{1/2}\)
against finite-rank orthogonal projections, so the
inequality is closed as well.

For the kernel cut off to a compact source region in
both variables, the objective is pairing \(X_u\) with a
Hilbert–Schmidt operator and is continuous on this compact
set. The uncut objective is its uniform limit because
the absolute error is bounded by the integral of
\(|K|\) on pairs with at least one source outside the
cutoff. This is uniform since
\(|\langle u_e(x),u_f(y)\rangle|\le1\). Thus the objective
reaches its supremum, with the maximizing operator
realized by countably many direction coordinates.

For a unit-ball Hilbert test \(G\in L^2(\gamma_2;
\mathbb R^3\otimes\ell^2)\) define
\[
 w^G_e=|T_sG_i-T_sG_j|,\qquad
 W_G(A)=\sum_e\int w^G_e\,d\nu_e .
\]
Pointwise alignment of directions gives
\(F_\infty(A)=\sup_{\|G\|_2\le1}W_G(A)\).
The smoothing is locally bounded with all fixed-order
local derivative bounds uniform on the unit ball, by
the point-source Mehler kernel in Gaussian \(L^2\).
For an attaining direction field and positive value
\(F=F_\infty(A)\), the test \(G=Z(u)/F\) attains the
dual value.

The minimum over balanced ordinary partitions is also
attained. The finite-rank functionals are lower
semicontinuous in label distance by the weighted
perimeter dual formula, which uses contractions from
the direction space to the physical plane. Their
supremum is therefore lower semicontinuous. It bounds
the original two-direction functional and hence
Gaussian perimeter. Compactness of balanced
bounded-perimeter partitions applies, and the simplex
has finite unrestricted energy.
[Lower semicontinuity and existence in the two-direction problem](THREE_CELL_ADAPTIVE_SURFACE_VARIATIONAL_PROBLEM.md)

## The unique scalar diagonal and quantitative weight rigidity

Use a positive value \(F=F_\infty(A)\) in this section. At zero
value all tagged measures vanish.

For an attaining direction field \(u\), first variation
in the unit balls gives
\[
          Ku=\lambda u,\qquad
          \lambda_e=|(Ku)_e|=F w^G_e ,
          \qquad G=Z(u)/F .
\]
Indeed a one-sided variation toward every admissible
field gives pointwise maximization against \(Ku\).
Where \(\lambda_e>0\) the field has norm one. Here
\(Ku\) is the Hilbert-valued integral of the scalar
tagged kernel against \(u\), and the formula with
\(G\) follows from the Mehler semigroup. Its product
integrals converge by the absolute kernel bound.

Vary toward a fresh coordinate orthogonal to the
entire direction Hilbert space, using
\((\sqrt{1-t^2f^2}\,u,tf)\) for any bounded real
\(f\) and small \(t\). The field remains in the unit
balls, and the extended countable Hilbert space has
the same unrestricted value. The quadratic term at
the maximum gives the scalar inequality
\[
         \langle f,Kf\rangle
                   \ \le\ \int \lambda f^2\,d\nu .           \tag{UR1}
\]
The expansion follows from the uniform Taylor formula
and absolute product integrability. A positive-measure
zero set of \(\lambda\) on one tag would contradict
(UR1) by using its indicator and the strictly positive
same-tag kernel. Thus \(\lambda>0\) almost everywhere,
\(|u|=1\) there, and \(\int\lambda\,d\nu=F^2\).

In particular,
\[
 F_\infty(A)^2
 = \min\left\{\int l\,d\nu:\ l\ge0,\
        \langle f,Kf\rangle\le\int l f^2\,d\nu
                 \ \text{for all bounded real }f\right\}.    \tag{UR2}
\]
The inequality for any feasible \(l\) follows by
summing the scalar inequalities for the direction
coordinates, and \(\lambda\) attains equality.
This minimizing diagonal is unique up to \(\nu\)-null
sets. If \(l\) attains equality, every coordinate
of one fixed attaining \(u\) has zero quadratic
defect for \(l-K\), since the sum of their
nonnegative defects vanishes. Polarizing against
arbitrary bounded scalar functions shows
\(l u_n=K u_n\) almost everywhere for each coordinate.
Since \(K u_n=\lambda u_n\) and \(\sum_n u_n^2=1\)
almost everywhere, this forces \(l=\lambda\).

The same scalar inequality gives a quantitative
comparison of the weights of every unit-ball test
to the unique weights on the interfaces:
\[
 \sum_e \int \frac{(w^G_e-\lambda_e/F)^2}{\lambda_e}\,d\nu_e
       \ \le\ 2\bigl(1-W_G(A)/F\bigr).                       \tag{UR3}
\]
For any bounded nonnegative scalar field \(f\), use
directions of length \(f_e\) aligned with
\(T_sG_i-T_sG_j\). Hilbert Cauchy–Schwarz followed by
(UR1) applied to their coordinates gives
\(\int fw^G\,d\nu\le(\int\lambda f^2\,d\nu)^{1/2}\).
Truncating \(w^G/\lambda\) in this inequality and passing
to the limit yields \(\int (w^G)^2/\lambda\,d\nu\le1\).
Expanding the square and using \(\int\lambda=F^2\)
proves (UR3). In particular every attaining test
has the same weights on the tagged supports, by
continuity and equality almost everywhere. The
directions or weights away from those supports may
still depend on the test.

At an actual planar noise-stability maximizer, the
first extra-coordinate stability inequality with
normal multiplier bounds every scalar field as
in (UR1). The normal field attains equality, and
summing the scalar inequality proves
\(F_\infty(A)^2=S'_r(A)\). The same conclusion
holds at the regular simplex for every positive
correlation, by the all-correlation interface gap.
Thus simplex minimality for \(F_\infty\), with a
strict inequality for nonsimplices, would exclude
a nonsimplex right tie, whose correlation derivative
is at most that of the simplex.
[Normal saturation and the scalar extra-coordinate inequality](THREE_CELL_ADAPTIVE_NORMAL_SATURATION.md#normal-saturation-and-the-scalar-interface-form)
[The simplex all-correlation interface gap](SIMPLEX_INTERFACE_SYMMETRY.md#reflected-kernels-close-the-three-cell-spectral-problem)

## Finite optimizing rank for compact interfaces

On any compact tagged source region \(D\), testing (UR1)
on one tag with \(f=1_D/(\lambda+\epsilon)\) gives
\(\int_D \lambda^{-1}\,d\nu_e<\infty\).
In fact \(M_r\) has a positive minimum on \(D\times D\)
and the right side is at most \(\int f\,d\nu_e\).
The left side is at least
\(2\min M_r\,(\int f\,d\nu_e)^2\). Passing to the limit
proves the bound. In particular \(\lambda\) is positive
as a continuous weight at every point of positive
lower one-dimensional density of that tag: if it
vanished there, local Lipschitz continuity and the
density would make the reciprocal integral infinite.

If all tagged interfaces are supported in one compact
set, the scalar operator on \(L^2(\nu)\) with kernel
\(K(x,y)/\sqrt{\lambda(x)\lambda(y)}\) is Hilbert–Schmidt.
The positivity and the upper bound by the identity follow
from the Mehler representation and (UR1), first testing
away from small \(\lambda\) and extending by density.
It follows that the eigenspace for eigenvalue one is
finite-dimensional. For any attaining field each
\(\sqrt\lambda\,u_n\) lies in this eigenspace, since
\(Ku_n=\lambda u_n\). Expanding all these coordinates in
one orthonormal basis shows that the Hilbert-valued
field takes values in a fixed finite-dimensional
direction subspace. Its dimension is at most that
eigenspace dimension. Thus compact interfaces need only
a finite optimizing rank, although the argument does not
give a universal rank bound of two.

## A quantitative patch comparison where all tags have density

Let \(A\) be a balanced unrestricted-direction minimizer,
and suppose each tagged interface at \(p\) has positive
lower one-dimensional density. Choose any attaining
test to define its frozen positive metric weights
\(c_e=\phi_2(p)\lambda_e(p)/F\). Here each value
\(\lambda_e(p)\) is the common continuous value provided
by the attaining tests on the support of that tag.
For spherical patches at \(p\), the minimizer is an
almost-minimizer of this fixed metric perimeter. More
precisely, for ordinary competitors \(B\) differing
compactly inside \(B(p,t)\), one has for small \(t\)
\[
   \sum_e c_e P_e(A;B(p,t))
      \le \sum_e c_e P_e(B;B(p,t)) + C_p t^{3/2}.             \tag{UR4}
\]
Here \(P_e\) is ordinary Euclidean pair length. The same
estimate holds for balls \(B(q,t)\) with
\(|q-p|\le C_0t\) and any fixed \(C_0\), with constants
also allowed to depend on \(C_0\).

Indeed, write \(\delta=F-W_G(A)\ge0\) for a unit-ball
test. On a fixed compact neighborhood, all weights
are uniformly Lipschitz, including the weights of a
fixed attaining test. If their difference on tag \(e\)
at \(p\) has magnitude \(d\), it has magnitude at
least \(d/2\) on a ball of radius proportional to \(d\)
for small \(d\). The lower density and the local
upper bound for \(\lambda\) make the left side of (UR3)
at least a constant times \(d^3\). For larger \(d\),
one can use a fixed small ball and adjust the constant.
Thus
\[
            |w^G_e(p)-\lambda_e(p)/F|\le C_p\delta^{1/3}.
                                                               \tag{UR5}
\]

An exterior compact flow restores the masses for any patch
with cost \(O(t^2)\) in \(W_G\), uniformly for all rank and all
unit-ball tests. The mass directions span because each
cell still has positive mass outside a fixed small
disk about \(p\). The continuous local derivative bounds
of the weights and the Gaussian tangential Jacobian give
the uniform cost. To obtain a local linear perimeter
upper bound, fill \(B(p,t)\) with one label at a good
slicing radius and use an attaining \(G\) of the corrected
fill. Its only new local perimeter has length at most
\(2\pi t\). Minimality and the dual formula give
\[
 0\le -\sum_e\int_{B(p,t)}w^G_e\,d\nu_e^A+Ct-\delta .
\]
Consequently \(\delta\le Ct\), and (UR5) makes all
local weights uniformly positive for small \(t\).
This gives \(P(A;B(p,t))\le Ct\), extended to every
small radius by monotonicity of perimeter in the ball.
The same upper bound holds on the stated comparable
nearby balls by containment.
[The exterior mass correction mechanism](THREE_CELL_ADAPTIVE_SURFACE_TWO_PHASE_REGULARITY.md#nearby-dual-tests-and-an-exterior-mass-correction)

For an arbitrary competitor, require consideration
only when its frozen weighted perimeter is no larger
than that of \(A\), since otherwise (UR4) holds.
Positivity of \(c_e\) makes both local unweighted
perimeters \(O(t)\). Use an attaining \(G\) of the
balanced correction of \(B\). Its weighted perimeter
increment over \(W_G(A)\), minus \(\delta\), is
nonnegative by minimality. Freezing the physical
weights at \(c_e\) makes an error at most
\(C_p t(\delta^{1/3}+t)\) by (UR5), uniform Lipschitz
bounds, and the volume repair. Finally
\(\sup_{\delta\ge0}(C_p t\delta^{1/3}-\delta)=O_p(t^{3/2})\).
This proves (UR4). Positive lower density for all
tags is used to compare with a fixed metric even
on newly inserted interfaces. At a missing tag
the same argument does not supply the pointwise
control (UR5).

## A three-tag partition can have only collinear optimizing tests

Unrestricted directions do not by themselves give strict
triangle inequalities. There are balanced ordinary
partitions with compact interfaces, all three tags of
positive measure, and every attaining test of rank one
for all sufficiently small positive \(r\).

Write the three tags in cyclic orientation, so that their
distinct label vectors have inner product \(-1\).
Suppose the three Gaussian interface masses are
\(a,b,c>0\), with \(a>2\max(b,c)\). Put a sign \(\eta=1\)
on the tag of mass \(a\) and \(\eta=-1\) on the other
two tags. Let \(m_r,M_r^*\) be the minimum and maximum
of the relative Mehler kernel on the union of the
compact interfaces squared. For sufficiently small
positive \(r\), \(m_r,M_r^*\to1\) and
\(m_r a>2 M_r^*\max(b,c)\).

Use the candidate direction field \(\eta\) times a
constant unit vector, and let \(l=\eta K\eta\) denote
the scalar product of \(\eta\) and \(K\eta\) pointwise.
This is positive. On the large tag every term in its
row sum after sign conjugation is positive. On each
small tag the positive contribution from the large
tag outweighs the possibly negative row contribution
from the other small tag by the assumed inequality.
For a scalar \(f=\eta v\), the defect of this diagonal is
\[
 \int l f^2\,d\nu-\langle f,Kf\rangle
  =\tfrac12\iint\eta(x)K(x,y)\eta(y)
                               (v(x)-v(y))^2\,d\nu(x)d\nu(y).
                                                               \tag{UR6}
\]
In the conjugated kernel the only negative block is
between the two small tags. If \(I_{ab},I_{ac},I_{bc}\)
denote the unweighted product integrals of the
squared difference between the indicated tags, then
\[
 I_{bc}\le \frac{2c}{a} I_{ab}+\frac{2b}{a} I_{ac}.
\]
This follows by bounding \((v(y)-v(z))^2\) by
\(2(v(y)-v(x))^2+2(v(z)-v(x))^2\) and averaging \(x\)
on the large tag. Its positive couplings dominate
that negative block strictly unless
\(I_{ab}=I_{ac}=0\). The same-tag terms in (UR6)
are also nonnegative, so (UR6) is a feasible diagonal
whose kernel consists of \(f\) proportional to
\(\eta\). The constant signed candidate attains
equality; applying the scalar equality coordinate by
coordinate makes every optimizing direction a
constant vector times \(\eta\). Every attaining
unit test is its normalized smoothed tensor, by
alignment and Hilbert equality. Thus its three
smoothed score rows are collinear everywhere.

Such interface masses can occur with exactly balanced
cells. Take label three outside a centered disk of
Gaussian mass \(2/3\). To form label one inside the
disk, use the radial graph
\(R(\theta)=b_0+\epsilon\sin(n\theta)\), and add a
fixed small disjoint cap cut from label two that
touches a short arc of the outer circle. The cap is
separated from the radial graph by a fixed gap.
Choose the offset \(b_0\) by continuity so that the
cap and the radial interior have total mass \(1/3\).
For a sufficiently small cap and amplitude this is
possible with a uniform gap to the outer circle.
The radial interior mass is independent of the
positive integer \(n\) for fixed \(b_0,\epsilon\),
by angular integration, so the same offset works
for all \(n\). All three tags are then present and
compact. Only the interface between labels one and
two has length tending to infinity, since the
oscillating radial boundary has length at least
\(\int_0^{2\pi}|R'(\theta)|\,d\theta=4\epsilon n\)
and lies in a fixed compact annulus with positive
Gaussian density. This gives the required interface
mass inequality. This construction concerns ordinary
competitors. It leaves strict triangle inequalities
at special classes of minimizers open.
