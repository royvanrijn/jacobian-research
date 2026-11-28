# Conditional masses do not give a sharp surface bound for the Gaussian simplex

Fix \(0<r<1\). The unoriented surface functional used in
the three-cell right-contact test is an average of a
quadratic perimeter functional in conditional Gaussian
partitions. The conditional partitions of the centered
simplex are translated \(120^\circ\) tripods. At every
nonzero translation, such a tripod fails to minimize the
conditional quadratic perimeter at its own masses.
Consequently an argument that minimizes the conditional
surface functional independently at each posterior mass
vector loses a strictly positive amount at the simplex.
The loss is at least order \(r\) as \(r\) tends to zero.
This concerns the lower bound on the unoriented surface
functional, which in turn is used to bound the correlation
derivative of an actual maximizer.

## The conditional surface identity

For a planar partition \(A\), use its three tagged pair
interfaces, their Gaussian arclength measures \(\nu_i\)
and the unoriented energy \(\mathcal J_r(A)\) from (OA2).
Here the index \(i\) denotes a pair. Assume the interfaces
are locally rectifiable and the energy is finite, as for
actual maximizers. For \(z\in\mathbb R^2\), let
\(A^z=(A-\sqrt r\,z)/\sqrt{1-r}\) denote the conditional
partition in standard Gaussian coordinates. Set
\[
 \mathscr Q(\tau)=2\sum_i\tau_i^2+\sum_{i<j}\tau_i\tau_j,
 \qquad \tau_i(A^z)=\nu_i(A^z)(\mathbb R^2).
\]
Changing variables in arclength gives
\[
 \tau_i(A^z)=\sqrt{1-r}
      \int M_{\sqrt r}(z,x)\,d\nu_i(x).
\]
Indeed the conditional density of \(X\) given the shared
standard Gaussian \(Z=z\) is
\(\phi_2(x)M_{\sqrt r}(z,x)\), and the transformation of
arclength combines with the planar density to give the
displayed factor. The Mehler semigroup and the nonnegative
entries of \(\mathscr Q\) then give
\[
           \mathcal J_r(A)=\frac1{1-r}
                         \mathbb E\mathscr Q(\tau(A^Z)).
                                                               \tag{CS1}
\]
The ordered off-diagonal sum in (OA2) has coefficient
\(1/2\), giving the unordered cross coefficient one here.
[Unoriented surface energy normalization](THREE_CELL_ORIGIN_ARM_SURFACE_COMPARISON.md#an-unoriented-surface-bound-from-stability)

## The translated tripod is not a conditional minimizer

Let \(T_p\) be a translated regular tripod with vertex \(p\)
and \(120^\circ\) unit outgoing directions \(v_1,v_2,v_3\).
Its interface areas are
\[
 \tau_i(p)=\phi_2(p)\int_0^\infty
                 \exp(-s^2/2-s p\cdot v_i)\,ds .
\]
The integral strictly decreases with \(p\cdot v_i\).
Since the three \(v_i\) span the plane and sum to zero,
the areas can all be equal only at \(p=0\).

Set \(w_i=\partial_{\tau_i}\mathscr Q=3\tau_i+\sum_j\tau_j\)
at the translated tripod. The only linear relation among
the \(v_i\) is their sum, so
\[
                  b_p=\sum_iw_i v_i=3\sum_i\tau_i(p)v_i
\]
is nonzero when \(p\ne0\). A smooth compact vector field
\(X\) supported in \(B_\epsilon(p)\), constant near \(p\)
with value \(b_p/|b_p|\) and bounded by one, gives
\[
 \dot\tau_i=-\phi_2(p)\,v_i\cdot X(p)
             -(p\cdot N_i)\int_{\Gamma_i} X\cdot N_i\,d\nu_i .
                                                               \tag{CS2}
\]
Here \(N_i\) is either unit normal to the straight ray
\(\Gamma_i\). To prove (CS2), differentiate the Gaussian
line integral under the flow: its integrand is
\(\partial_s(X\cdot v_i)-x\cdot X\). Integration by parts
in arclength gives the boundary term and cancels the
tangential position contribution. The remaining normal
position equals \(p\cdot N_i\) on the whole ray.

The cell mass derivatives of this flow are \(O(\epsilon)\),
by boundary flux on the three short pieces. Choose two
fixed compact normal bumps on separate regular interfaces
away from \(p\), whose mass derivative vectors are
independent. For small \(\epsilon\) they are disjoint
from the junction variation. Composing the flow with
their two flows and applying the implicit function theorem
preserves the three masses exactly. The two initial
correction velocities are \(O(\epsilon)\), and their
perimeter derivative on each ray is \(O(|p|\epsilon)\):
they have no junction contribution and the Gaussian
normal drift in (CS2) is constant on each ray.
The first variation of \(\mathscr Q\) is therefore
\[
              -\phi_2(p)|b_p| + O(|p|\epsilon).
\]
For each \(p\ne0\), fixing sufficiently small
\(\epsilon>0\) makes this negative. The deformation is an
ordinary partition with the exact original masses.

There is a uniform estimate near the centered tripod. Taylor
expansion of the three areas and
\(\sum_i v_iv_i^{\mathsf T}=\tfrac32 I\) give
\[
             b_p=-\tfrac92\phi_2(0)\,p+O(|p|^2).
\]
The error bound \(O(|p|\epsilon)\) just obtained is uniform
for \(|p|\) sufficiently small. Indeed use the same two
bumps translated with \(p\). The inverse mass derivative
matrix stays bounded near zero, and the boundary lengths
and Gaussian densities on their compact supports are
uniformly bounded. Choose \(\epsilon>0\) fixed small.
The initial derivative of the exactly corrected
\(\mathscr Q\) is at most \(-c|p|\). With this fixed
\(\epsilon\), the second derivative for small flow time
is uniformly bounded. This follows by differentiating the
transported line integrals, and the first two derivatives
of the balancing corrections are uniformly bounded by the
same invertible flux matrix and smooth Gaussian densities.
All the vector fields involved range over a bounded
smooth family with compact supports in translated
coordinates, including all unit junction directions.
A positive flow time proportional to \(|p|\), with a
small fixed proportionality constant, therefore decreases
\(\mathscr Q\) by at least \(c'|p|^2\), uniformly for small \(p\).

## The gap in a posterior mass comparison

Suppose \(F\) is a measurable lower bound depending only
on the three masses \(m\) of an ordinary conditional
Gaussian partition:
\[
                 F(m)\le\mathscr Q(\tau(C))
                     \quad\text{whenever }\gamma(C)=m .
\]
Assume evaluations in (CS1) are integrable. For the
centered simplex \(T\), the conditional partition \(T^z\)
has vertex \(p_z=-\sqrt{r/(1-r)}\,z\) and the same three
outgoing directions. Thus for almost every \(z\),
\[
 F\bigl(\gamma(T^z)\bigr)<\mathscr Q(\tau(T^z)).
\]
In particular the lower estimate obtained by integrating
\(F\) loses a strictly positive amount from
\(\mathcal J_r(T)=S'_r(T)\). The squared quantitative
improvement above gives
\[
 \mathcal J_r(T)-\frac1{1-r}\mathbb E F(\gamma(T^Z))
  \ \ge\ \frac{c}{1-r}\,
       \mathbb E\bigl[|p_Z|^2\,1_{\{|p_Z|\le c_0\}}\bigr]
  \ \ge\ c' r
\]
for all sufficiently small \(r>0\). For the last inequality
retain the fixed event \(0<|Z|\le1\).

For the optimized constant pair-vector test, the centered
simplex uses the equilateral directions, so its value there
is again the left side of (CS1). The same obstruction applies
to a conditional minimization of that equilateral functional.
It leaves open a comparison which uses a single partition
across the different conditional laws.

## Optimizing the constant pair directions after conditioning

The loss persists if the lower estimate first minimizes at
each posterior mass and then chooses the best common pair
directions. For unit \(u_i\) in the calibration plane, put
\[
 \mathscr Q_u(\tau)=2\sum_i\tau_i^2
                    -2\sum_{i<j}(u_i\cdot u_j)\tau_i\tau_j .
\]
For any measurable lower bound \(F_u(m)\le\mathscr Q_u(\tau(C))\)
at each mass vector \(m=\gamma(C)\), there is a positive
gap between the best resulting bound at \(T\) and its
derivative. Uniformly for small positive \(r\),
\[
 S'_r(T)-\sup_u\frac{1}{1-r}
                  \mathbb E F_u(\gamma(T^Z))\ \ge\ c r .        \tag{CS3}
\]
The same argument gives a strictly positive gap for each
fixed \(r\in(0,1)\).

To see this, write \(b_u(p)=\sum_i
\partial_{\tau_i}\mathscr Q_u(\tau(p))v_i\) and
\(q=|u_1+u_2+u_3|\). The identity
\[
 u_i\cdot u_j+\tfrac12
     =\tfrac12 q^2-(u_1+u_2+u_3)\cdot u_k
     \qquad (\{i,j,k\}=\{1,2,3\})
\]
shows that \(b_u(p)-b_p=O(q)\) uniformly for small \(p\)
and \(q\). Thus if \(q\le c_1|p|\) and \(p\) is
sufficiently small, then \(|b_u(p)|\ge c_2|p|\), with
positive constants independent of \(u,p\).
The exactly balanced junction motion from (CS2) applies
with this \(b_u(p)\) in place of \(b_p\).
All derivative weights are bounded for \(p\) in the
fixed small disk and all \(u\), so the \(O(|p|\epsilon)\)
first derivative error and the quadratic flow error stay
uniform. Consequently \(\mathscr Q_u\) also decreases
by at least a fixed multiple of \(|p|^2\) in this range.

On the other hand, before any conditional minimization
the integrated loss at \(T\) from using \(u\) instead
of equilateral directions is exactly \(m_o(r)q^2\):
the three cross energies \(m_{ij}(T)\) are the same
positive number \(m_o(r)\), and the expansion in (OA12)
gives this difference. The number \(m_o(r)\) has a
positive limit as \(r\downarrow0\), by Gaussian domination
in the ray integral. If \(q\ge c_1\sqrt r/2\), this
integrated loss alone is at least a fixed multiple of \(r\).
If \(q<c_1\sqrt r/2\), retain the event
\(1/2\le |Z|\le1\); for small \(r\) it has
\(\sqrt r/2\le|p_Z|\le c_0\).
The conditional \(\mathscr Q_u\) improvement on this event
has expectation at least a fixed multiple of \(r\).
In both cases the conditional lower bounds can only
increase the indicated loss. For a fixed positive \(r\)
use a sufficiently small fixed annulus in \(p_Z\) of
positive Gaussian probability in place of this event,
and split the two cases by \(c_1\) times its inner radius.
This proves the claims.
[Optimized pair directions, OA12](THREE_CELL_ORIGIN_ARM_SURFACE_COMPARISON.md#optimal-pair-directions-and-strictness-at-a-maximizer)
