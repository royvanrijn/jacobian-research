# Stationary correlations of a fixed Gaussian maximizer

Let \(A\) be a balanced three-cell Gaussian partition that maximizes
noise stability at one positive correlation. Work in its active plane.
If \(A\) is nonsimplex, some compact smooth exactly balanced deformation
has a nonzero first shape derivative as a holomorphic function of the
correlation. This derivative vanishes at every positive correlation
where the same \(A\) is stationary. Its zeros on compact subsets of
the open unit disk can be counted uniformly over nonsimplex maximizers
whose reference maximizing correlations lie in a compact positive
band. The same witness tests higher-order contact with the simplex.
For any cell count the degree-one witness exists whenever the moments
have full rank.

## A shape derivative with a nonzero degree-one term

For a compactly supported smooth vector field \(X\), use the linear
mass fluxes
\[
 M_i(A,X)=\int_{A_i}\operatorname{div}_\gamma X\,d\gamma_2,
 \qquad \operatorname{div}_\gamma X=\operatorname{div}X-x\cdot X.
\]
Their sum is zero. The vector of fluxes ranges over all of
\(\mathbf1^\perp\). In fact, if \(c\) annihilates its image, the
weighted divergence identity for all compact smooth fields implies
that the distributional gradient of \(\sum_i c_i1_{A_i}\) is zero.
To pass to the ordinary divergence test, absorb the positive
Gaussian density into the field. This function is constant almost
everywhere. Positive masses force all \(c_i\) to be equal. Select
two fields realizing independent fluxes at \(A\).
For a field with zero mass flux, use its flow followed by the flows of
two such fields. The implicit function theorem sets the two independent
volumes exactly. The two correction parameters vanish to second order,
so this gives a smooth balanced deformation \(A_t\) with initial field
\(X\).
[Regular interface arcs and isolated triple junctions](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#regular-interface-arcs-and-isolated-triple-junctions)

Write \(T_z\) for the operator with Hermite multipliers \(z^n\) for
\(|z|<1\), and define
\[
 H_{A,X}(z)=\sum_i\int_{A_i}
   \bigl[X\cdot\nabla T_z1_{A_i}+
             (\operatorname{div}_\gamma X)T_z1_{A_i}\bigr]d\gamma_2.
                                                               \tag{SCW1}
\]
For real \(0<z<1\), this is half the first derivative of
\(S_z(A_t)\) at \(t=0\). To compute it, transport both integration
variables in the Mehler integral to the fixed cells, differentiate,
and use symmetry to combine the two terms. Differentiating the
Gaussian Jacobian gives the weighted divergence term. The deformation
is accepted for both signs of \(t\), so stationarity of \(A\) at \(z\)
implies \(H_{A,X}(z)=0\).

The function in (SCW1) is holomorphic on the open disk. For instance
the Hermite multipliers and the gradient series converge in operator
norm locally uniformly, using
\(\|\nabla\Pi_n\|_{L^2\to L^2}=\sqrt n\) in degree \(n\) and the
summability of \(\sum_n\sqrt n R^n\) for \(R<1\).
Since the cell masses are equal, \(H_{A,X}(0)=0\). If
\(m_i=\int_{A_i}x\,d\gamma_2\), then \(\Pi_1 1_{A_i}=m_i\cdot x\),
and Gaussian divergence gives
\[
 H_{A,X}'(0)=\sum_{i<j}\int_{\Sigma_{ij}}
       (m_i-m_j)\cdot x\, (X\cdot N_{ij})\,d\nu_{ij}.       \tag{SCW2}
\]
Here \(N_{ij}\) points out of \(A_i\), and \(\nu_{ij}\) is Gaussian
surface measure. The integrals use compact support and the usual
locally finite perimeter divergence formula.

Suppose (SCW2) vanished for every field with zero mass flux. By the
surjectivity of the flux map there would be constants \(c_i\) such
that \(H_{A,X}'(0)=\sum_i c_i M_i(A,X)\) for every compact smooth
field. Test on each open regular pair facet, with arbitrary compact
normal speed there. It follows that its containing line has equation
\((m_i-m_j)\cdot x=c_i-c_j\). The three moments sum to zero and span
the active plane. The map \(x\mapsto (m_i\cdot x)_i\) is therefore
an isomorphism from the active plane to \(\mathbf1^\perp\).
Subtract the average of the \(c_i\) and choose \(p\) with
\(m_i\cdot p=c_i\) for all \(i\). Every regular facet then lies
on a line through \(p\). Thus \((x-p)\cdot D1_{A_i}=0\) for each
cell, since the regular facets carry the perimeter measure. In polar
coordinates about \(p\) this says that the radial weak derivative
of every indicator is zero away from \(p\). The indicators are
invariant under positive dilations there up to null sets. The
classification of conical maximizers about an arbitrary center
gives the centered regular simplex. This proves that a nonsimplex
maximizer admits a field \(X\) with zero flux and \(H_{A,X}'(0)\ne0\).
[Full rank of the cell moments](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#one-active-dimension-is-impossible)
[Conical maximizers about arbitrary centers](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#cones-centered-at-an-arbitrary-point)

## Uniform finite stationary correlations and finite shape jets

Let \(J\subset(0,1)\) be a compact interval of reference maximizing
correlations, and fix \(0<b<R<1\). There exist \(M\ge\delta>0\)
depending on \(J,R\) such that every nonsimplex balanced maximizer
\(A\) at a reference correlation in \(J\) has a witness \(X\) with
\[
 F(z)=H_{A,X}(z)/z,\qquad |F(0)|\ge\delta,
               \qquad \sup_{|z|\le R}|F(z)|\le M.           \tag{SCW3}
\]
At the origin use the removable value \(H_{A,X}'(0)\).

To prove the uniform assertion, consider one such \(A\) and the two
fixed flux fields from the volume repair at \(A\). For any sufficiently
nearby partition \(B\), correct \(X\) by a linear combination of those
same fields to give zero flux at \(B\). The fluxes are continuous for
strong Gaussian label convergence by their volume-integral formulas.
Thus the flux matrix stays invertible and the corrected field varies
continuously in its finite-dimensional family. The series bound in
(SCW1) gives a uniform bound for \(H\) on \(|z|\le R\), and the maximum
principle gives the asserted bound for \(H/z\) using its boundary
values. The derivative at zero stays separated from zero by the
continuity of (SCW1) and its degree-one term. In the fixed active plane,
the nonsimplex maximizers over \(J\) form a strongly compact family.
This uses strong positive-correlation maximizer compactness and uniform
measurable local isolation of the simplex on \(J\), which keeps every
subsequential limit of nonsimplex maximizers off the simplex orbit.
Take a finite subcover of the neighborhoods just constructed.
[Strong maximizer compactness](SMALL_POSITIVE_CORRELATION_INTERVAL.md#a-rigorous-continuation-boundary-in-correlation)
[Uniform local comparison in Gaussian label distance](GAUSSIAN_SIMPLEX_DIMENSION_FREE_LOCAL_GROWTH.md#uniform-repair-and-comparison)

The number of zeros of this \(F\) in \(|z|\le b\), counted with
multiplicity, is at most
\[
                 N=\left\lfloor
                      \frac{\log(M/\delta)}{\log(R/b)}
                    \right\rfloor.                         \tag{SCW4}
\]
In fact, scale the disk of radius \(R\) to the unit disk and divide
out the Blaschke factor at each zero in the smaller disk, with
multiplicity. The quotient has boundary modulus at most \(M\)
and its value at zero has modulus at least \(\delta(R/b)^n\)
when \(n\) factors are used. The maximum principle gives (SCW4).
At a positive stationary correlation \(r\le b\), the zero of \(F\)
has the same order as that of \(H\). In particular the same fixed
shape cannot be stationary at more than \(N\) distinct such
correlations, and this witness has at least one nonzero correlation
derivative of order at most \(N\) at each of them. The fixed-shape
bound does not count correlations with a different optimizing shape.

## The shape test at a higher-order contact

Suppose \(A\) ties the simplex \(T\) at \(r\), and distinct
simplex-optimal correlations tend to \(r\). Use the fixed-shape
deficit \(D(s)=S_s(T)-S_s(A)\). If it has order \(m\) at \(r\),
write \(D(r+h)=d_m h^m+o(|h|^m)\) with \(d_m\ne0\).
For any compact smooth mass-tangent field \(X\), all derivatives
of \(H_{A,X}\) of orders strictly less than \(m/2\) vanish at \(r\).
If \(m=2j\), the derivative at order \(j\) satisfies
\[
 \left|\frac{H_{A,X}^{(j)}(r)}{j!}\right|^2
             \le d_m q_X,\qquad
 q_X=-\tfrac12\partial_t^2 S_r(A_t)|_{t=0}\ge0.             \tag{SCW5}
\]
The second derivative is taken on the exact balanced path constructed
from \(X\). For speeds supported in regular facet interiors, \(q_X\)
is \(Q_0(X\cdot N)\). The argument also works when the field meets a
junction, by differentiating the transported bulk integral.
[The ordinary balanced shape form and the second-order contact test](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#a-first-crossing-obstruction-at-the-level-of-hermite-weights)

For a fixed \(X\) and its exact balanced deformation, Taylor expansion
in the shape parameter gives
\[
 S_{r+h}(A_t)-S_{r+h}(A)
    =2tH_{A,X}(r+h)-t^2 q_X+o(t^2)
                      \qquad ((h,t)\to(0,0)).
\]
The second shape derivative is jointly continuous here: transport both
variables to fixed cells. The smooth compactly supported fields and a
correlation band inside \((0,1)\) give integrable Gaussian domination
of the kernel and its derivatives, including the polynomial factors.
At each simplex-optimal \(r+h\) the increment is at most \(D(r+h)\).
If the first nonzero Taylor coefficient of \(H_{A,X}(r+h)\) is
\(a_jh^j\) with \(j<m/2\), use \(t=c h^j\) with a real constant
\(c\). Here \(j\ge1\) by stationarity at \(r\). Dividing by
\(h^{2j}\) and passing along the optimal correlations gives
\(2c a_j-c^2 q_X\le0\) for every \(c\). A sufficiently small
\(c\) with the sign of \(a_j\) contradicts this. If \(m=2j\), the
same choice gives \(d_m-2c a_j+c^2q_X\ge0\) for every \(c\),
which yields (SCW5), including if \(q_X=0\).

For a nonsimplex fixed shape the witness (SCW3) is nonzero, so at an
endpoint of simplex optimality its least nonzero mixed order \(j\)
is at most \(N\) if \(r\le b\). The contact order obeys \(m\le2j\).
When \(m=2j\), (SCW5) also requires strictly positive shape curvature
in this detecting direction. These constraints apply to an arbitrary
ordinary maximizer, without an assumption about the ends or the
reference pair lines. They leave the existence of a nonsimplex tied
maximizer open.

## Approximate stationarity and motion with correlation

For a balanced partition \(A\), let its linear assignment
regret be
\[
 \mathcal R(A,r)=\sup_{B\ {\rm balanced}}
       \sum_i\langle 1_{B_i}-1_{A_i}, T_r1_{A_i}\rangle .
                                                               \tag{SR1}
\]
A stationary score winner has zero regret. This follows by
integrating its winner inequalities and canceling the label
prices by balance. More generally,
\[
                         V(r)-S_r(A)\ \ge\ 2\mathcal R(A,r),
                                                               \tag{SR2}
\]
where \(V(r)\) is the global balanced maximum. Indeed the
difference between the full stability gain and twice the
linear gain for any \(B\) is the nonnegative Ornstein–Uhlenbeck
quadratic form of \(1_B-1_A\).

Fix a compact positive correlation band \(J\). For each
fixed nonsimplex actual maximizer \(A\) at some correlation
in \(J\), choose the stationary witness \(H_A\) from the
finite family above, corrected to be tangent to the two
volume constraints at \(A\). There is a uniform \(c_J>0\)
such that
\[
                    \mathcal R(A,r)\ge c_J |H_A(r)|^2
                         \qquad(r\in J).                   \tag{SR3}
\]
Neither this inequality nor its proof needs regular
interfaces of \(A\) at the tested correlation \(r\).

To see it, let \(X_A\) be the corrected compactly supported
smooth field and \(Y_1,Y_2\) the volume repair fields in
its witness neighborhood. Compose the time-\(t\) map of
\(X_A\) with the time-\(u_i\) maps of \(Y_i\).
The two independent transported Gaussian masses are
smooth functions of \((t,u)\), by pulling their densities
back to the fixed cells. Their first \(t\)-derivative at
the origin is zero and their \(u\)-Jacobian is invertible.
Over the witness neighborhood the inverse stays bounded.
The parameter derivatives through second order are bounded
uniformly over all measurable cells because the fields
are supported on a fixed compact set and their coefficients
stay bounded. Consequently the implicit function theorem
gives a balanced deformation with \(u=O(t^2)\) on a common
interval \(|t|\le t_0\). The constants can be chosen for
the whole finite family of witnesses.

The gain in the linear assignment functional along this
deformation is \(tH_A(r)+O(t^2)\) uniformly for \(r\in J\).
For this estimate, use the same pullback with the score
\(T_r1_{A_i}\) held fixed and differentiate in space, whose
first two derivatives are uniformly bounded on the compact
support by Gaussian convolution on \(J\). The \(O(t^2)\)
volume repairs leave the first derivative unchanged.
Taking \(t=cH_A(r)\) for a sufficiently small uniform
\(c>0\) makes the gain at least \(c|H_A(r)|^2/2\),
proving (SR3).

In particular, for some \(N_J\ge1\), \(C_J<\infty\),
and every sufficiently small \(\epsilon>0\), the set
\[
                    \{r\in J:V(r)-S_r(A)\le\epsilon\}
\]
has length at most \(C_J\epsilon^{1/(2N_J)}\), uniformly
over the fixed nonsimplex maximizing shapes \(A\) selected
anywhere on \(J\). Here is one version of the uniform
analytic estimate. Rescale the complex parameter so
\(H_A(z)/z\) is bounded on the closed unit disk and the
tested band lies in \(|z|\le b_0<1\). Its modulus at zero
has a uniformly positive lower bound from the stationary
witness. Fix \(b_0<r_0<1\). The finite Blaschke argument
bounds the number of zeros in \(|z|\le r_0\) by a uniform
integer \(N_J\). After dividing out their Blaschke factors,
Harnack's inequality bounds the nonvanishing product below
on \(|z|\le b_0\). Hence on the tested real band the modulus
of \(H_A\) is bounded below by a constant times the
\(N_J\)-th power of the distance to this finite set of
complex zeros, with the distance capped at one. The
bands here stay away from zero, so division by the
parameter costs only a uniform factor. Apply (SR2) and
(SR3), and cover the real sublevel set by the intervals
near those zeros. This is the same bounded-holomorphic
argument as for the fixed value contact, applied now to
a shape derivative.
[Uniform holomorphic sublevel bound](GAUSSIAN_FIXED_PARTITION_CONTACT_BOUNDS.md#a-uniform-bound-on-contacts-and-near-contacts)

There is also a lower bound on shape motion. Fix a
nonsimplex maximizer \(A_0\) at \(r_*\in J\) and use
its witness neighborhood. If \(H_{A_0}\) has a zero of
order \(m\) at \(r_*\), any exact stationary winner
\(A'\) at a nearby \(r\) satisfies
\[
       \inf_{O,\pi}\,\sum_i\gamma_2(A'_i\triangle O A_{0,\pi(i)})
                         \ \ge\ c_{A_0,r_*}|r-r_*|^m .     \tag{SR4}
\]
Here \(O\) ranges over orthogonal transformations and
\(\pi\) over label permutations. The flux correction
depends Lipschitz continuously on the label distance,
since it uses fixed bounded Gaussian divergences and
a uniformly invertible finite matrix. The difference
of the two shape derivatives on the same compact
observation support is also Lipschitz in that distance.
Indeed the Mehler kernel relative to Gaussian measure
and its first spatial derivatives are uniformly bounded
there in the source coordinate, since \(J\) stays away
from zero and one. Stationarity makes the shape
derivative of the aligned \(A'\) zero at \(r\), so its
distance to \(A_0\) bounds \(|H_{A_0}(r)|\) from above
up to a constant. The first nonzero Taylor term bounds
it from below by a constant times \(|r-r_*|^m\).
Partitions outside the witness neighborhood already
have a fixed positive orbit distance to \(A_0\).
This proves (SR4). The estimate allows the maximizing
shape to change with \(r\). It quantifies a necessary
change rather than excluding a moving maximizing branch.

## Full-rank partitions with any number of labels

Let \(A\) be a balanced actual maximizer with \(k\ge3\) cells and
centroids of rank \(k-1\). Work in a fixed dimension containing
its active space. Unless \(A\) is the regular simplex with its
cylindrical extension, it has a compact smooth mass-tangent
witness with \(H_{A,X}'(0)\ne0\). For example, every fixed
nonsimplex maximizer of this rank can be stationary only at
finitely many correlations in \((0,b]\) for each \(b<1\).
The higher-contact test (SCW5) and pointwise orbit motion
estimate (SR4) apply in this rank class. For a strongly compact
family of nonsimplex actual maximizers of full rank in fixed
dimension, the stationary-root and linear-regret estimates
(SCW3)–(SCW4) and (SR1)–(SR3) are uniform on compact positive
bands. Such compactness of the full-rank nonsimplex subfamily
is an additional hypothesis in the general \(k\) case.

The bulk derivative and mass correction work with the same
formulas in any fixed dimension, with \(k-1\) independent
flux fields. The regular pair facets represent the perimeter
measure. One can use the regularity of positive-correlation
maximizers here. Alternatively, the strict adjusted-score
winners form a locally semianalytic partition. On a
codimension-one stratum there are two sides, and strata
of smaller dimension carry no perimeter measure.
[Heilman–Tarter, regularity of maximizers](https://arxiv.org/pdf/2011.05583)

If \(H_{A,X}'(0)=0\) for all zero-flux fields, the annihilator
argument in (SCW2) puts every regular label-changing facet
in one of the planes
\[
                (m_i-m_j)\cdot x=c_i-c_j.
\]
The full-rank hypothesis gives a common point \(p\) in the
centroid span. On each connected chamber of this arrangement
all indicators have zero distributional gradient, because
all perimeter measure is carried by the displayed planes.
Thus the chamber labels represent the partition, making it
a finite polyhedral cone partition about \(p\). Centering of
balanced stationary polyhedral cones makes it conical about
the origin at the one reference correlation. The full-rank
classification of origin-conical adjusted-score winners then
identifies the standard simplex. The holomorphic zero count,
mixed shape test, assignment gain and motion estimates use
only the smooth transported densities, the witness and the
finite-dimensional flux correction, so they apply as claimed.
[Centering of balanced polyhedral cones](AFFINE_SIMPLEX_FAN_RIGIDITY.md#centering-balanced-polyhedral-cones-from-one-correlation)
[Full-rank conical classification](AFFINE_SIMPLEX_FAN_RIGIDITY.md#harmonic-taylor-terms-for-general-conical-winners)
