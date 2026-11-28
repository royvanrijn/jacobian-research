# A thin third phase lowers the pair surface bound

For every \(0<r<1\), there are balanced ordinary planar
partitions whose unoriented pair surface energy is strictly
below that of the centered regular simplex. This holds for
the equilateral pair directions and for the optimal
directions that are constant on each pair. Both surface
functionals fail lower semicontinuity under strong Gaussian
label convergence. The construction inserts a thin layer
of the third label along a short part of a simplex ray.

The two new interfaces converge to the same segment.
Their combined pair measure changes the unoriented
interaction, although the indicator partition converges
back to the simplex. Allowing the directions to vary along
each pair recovers the original tensor measure in this
example.
[The right-contact surface bounds (OA2) and (OA12)](THREE_CELL_ORIGIN_ARM_SURFACE_COMPARISON.md#optimal-pair-directions-and-strictness-at-a-maximizer)

## The change in the surface pairing

Label the three interface types of the centered regular
simplex \(T\) by \(e,f,g\). Write \(\nu_i\) for their
Gaussian arclength measures and
\(\langle\alpha,T_r\beta\rangle=\iint M_r\,d\alpha\,d\beta\)
for the Mehler pairing of surface measures. Set
\(H_i(x)=\int M_r(x,y)\,d\nu_i(y)\).
The equilateral energy of the tagged measures is
\[
 \mathcal J_r(\nu)=
 2\sum_i\langle\nu_i,T_r\nu_i\rangle
        +\sum_{i<j}\langle\nu_i,T_r\nu_j\rangle .
\]
Let \(I_\ell\) be a segment of length \(\ell>0\) starting
at \(a v_e\) on ray \(e\), where \(v_e\) is the ray's outward
unit tangent and \(a>0\). Let \(\mu_\ell=\nu_e|_{I_\ell}\).
Consider first the limiting surface measures
\[
 \widetilde\nu_e=\nu_e-\mu_\ell,\qquad
 \widetilde\nu_f=\nu_f+\mu_\ell,\qquad
 \widetilde\nu_g=\nu_g+\mu_\ell .
\]
Expanding the quadratic form gives exactly
\[
 \mathcal J_r(\widetilde\nu)-\mathcal J_r(\nu)
 =\int_{I_\ell}\big[4(H_f+H_g)-2H_e\big]\,d\mu_\ell
                 +5\langle\mu_\ell,T_r\mu_\ell\rangle .    \tag{PW1}
\]

Let \(\sigma=\sqrt{1-r^2}\). Conditional Gaussian integration
over a ray with unit tangent \(v_i\) and unit normal \(n_i\)
gives
\[
 H_i(x)=\frac1{\sigma}
       \phi_1\!\left(\frac{r x\cdot n_i}{\sigma}\right)
       \Phi\!\left(\frac{r x\cdot v_i}{\sigma}\right).
\]
Here \(\phi_1,\Phi\) are the standard normal density and
distribution function. Along \(x=s v_e\), \(H_e\) tends
to \(\phi_1(0)/\sigma>0\), while \(H_f,H_g\) tend to zero
because their supporting lines are nonparallel to ray \(e\).
Fix \(a\) large enough that the linear integrand in (PW1)
is strictly negative, bounded away from zero, on
\([a,a+1]v_e\). For \(0<\ell\le1\) its integral is bounded
above by a negative constant times \(\ell\), since the
Gaussian arclength density has a positive minimum on
this fixed compact segment. The quadratic term is
\(O(\ell^2)\), because \(M_r\) is bounded on the segment
product. Thus (PW1) is strictly negative for sufficiently
small \(\ell>0\).

The optimal constant pair directions have the same first
variation. Indeed, put \(m_{ij}=\langle\nu_i,T_r\nu_j\rangle\).
The three off-diagonal values for \(T\) coincide and are
positive. Their changes under \(\nu\mapsto\widetilde\nu\)
are \(O(\ell)\). In the formula (OA12), the three numbers
\(c_i=\sqrt{m_{ij}m_{ik}/m_{jk}}\) therefore differ from
their common positive reference value by \(O(\ell)\).
They satisfy strict triangle inequalities for small
\(\ell\). The gain of the optimized energy over
\(\mathcal J_r\) is then
\(\tfrac12\sum_{i<j}(c_i-c_j)^2=O(\ell^2)\), and at \(T\)
it is zero. Consequently the optimized energy of
\(\widetilde\nu\) is strictly below its value at \(T\) too.
[Exact optimization over constant pair directions (OA12)](THREE_CELL_ORIGIN_ARM_SURFACE_COMPARISON.md#optimal-pair-directions-and-strictness-at-a-maximizer)

## Realizing the measures by balanced partitions

Fix such a segment \(I_\ell\). In its tangent and normal
coordinates assign the rectangle
\([a,a+\ell]\times[-\delta,\delta]\) to the label absent
from the original \(e\) interface. For sufficiently small
\(\delta>0\) the rectangle meets only the two original
incident cells. It removes the original interface on
\(I_\ell\) and creates one long wall of each of types
\(f,g\). The two short caps have total length \(O(\delta)\).
These are ordinary piecewise smooth interfaces.

The volume errors are \(O(\delta)\). Use two fixed compact
normal bumps supported on different simplex rays away
from the rectangle and the origin. Their two first-order
volume-change vectors are independent in the plane of
zero total volume change, because each is a nonzero
multiple of its pair's label-transfer vector. The inverse
function theorem therefore gives flow times \(O(\delta)\)
that repair the masses exactly. Their supports leave all
four sides of the rectangle untouched. Denote the repaired
partition by \(A_\delta\).

For each pair label the surface measures of \(A_\delta\)
converge locally weakly to \(\widetilde\nu_i\). The long
walls converge with arclength to \(I_\ell\), the caps
vanish in total arclength, and the two repaired wall
pieces converge smoothly to the originals. All changes
lie in a fixed compact set. The remaining ray tails have
uniform Gaussian decay in every Mehler surface pairing
by the joint Gaussian density in (OA3). Compact weak
convergence and this tail bound imply convergence of
every \(m_{ij}(A_\delta)\) to the corresponding pairing
for \(\widetilde\nu\). Both surface energies are continuous
in these pairings near the indicated limit. They are
therefore strictly below their simplex values for
sufficiently small \(\delta>0\).

At the same time \(A_\delta\to T\) in Gaussian label
measure. This proves the lower semicontinuity failure at
the simplex with fixed \(\ell\). One may instead choose
\(\ell\to0\) and sufficiently small positive \(\delta\)
at each step. Then the total Gaussian perimeter also
converges to that of \(T\), while the individual
inequalities remain strict.

## Varying the pair directions along the walls

Orient the pair labels cyclically and write \(d_i\) for
their label vectors, so \(d_e+d_f+d_g=0\). A choice
of a unit direction field on each pair gives the tensor
surface measure with density \(d_i\otimes u_i\).
Keep the simplex equilateral directions on all unchanged
walls and assign direction \(-u_e\) to both newly inserted
long walls. In the thin layer limit their combined tensor
density on \(I_\ell\) is
\[
          -d_f\otimes u_e-d_g\otimes u_e=d_e\otimes u_e .
\]
It recovers the tensor density removed from ray \(e\).
The caps vanish and the balancing flows have the same
vanishing effect as before. Hence this particular
direction field has energy tending to the original
simplex energy. The counterexample exploits the use of
a single direction on each entire pair rather than a
restriction on all measurable direction fields.
[The direction-field supremum for actual maximizers](THREE_CELL_WEIGHTED_PERIMETER_CONTACT.md#a-direction-supremum-for-the-correlation-derivative)

In this construction the adaptive surface energy is
in fact strictly above the simplex value for sufficiently
small layer width. Use the original simplex tensor as a
comparison in the smoothed Hilbert space. The pair row
differences for that comparison, before division by its
norm, are
\[
           F_i=2H_i u_i-\sum_{j\ne i}H_j u_j .
\]
On ray \(e\), symmetry gives \(H_f=H_g=h>0\);
write \(a=H_e\). Then \(F_e=(2a+h)u_e\). The other
two vectors have common \(u_e\) component
\(-a-h/2\) and opposite perpendicular components of
magnitude \(3\sqrt3h/2\). Therefore
\[
                        |F_f|+|F_g|>2a+h .
\]
Keep all unchanged directions and align the two new
walls with \(F_f,F_g\), evaluated at their projections
onto the original segment. In the limiting tensor
the cross term with the simplex strictly increases,
by the displayed inequality integrated on the
segment. Its squared norm is the simplex squared
norm plus twice that positive difference plus the
nonnegative squared norm of the tensor change.
The same compact and tail convergence used above
passes this strict inequality to sufficiently
thin actual strips.
