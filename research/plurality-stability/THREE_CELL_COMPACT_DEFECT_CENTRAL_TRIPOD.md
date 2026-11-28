# A centered tripod germ inside every compact simplex defect

Let \(A\) be a balanced planar Gaussian noise-stability maximizer at
\(0<r<1\) agreeing with the centered regular simplex \(T\) outside a
bounded set. Then a disk about the origin has exactly the simplex or
its antipode as its strict-winner partition, up to null sets. The
radius of the disk has a positive lower bound uniform on compact
positive-correlation bands. The location and size of the exterior
exceptional set do not enter this bound.

An antipodal germ forces a genuine triple junction on each of the six
open reference half-rays, hence at least seven including the origin.
With a simplex germ, a strict spectator interval on a forward reference
ray forces at least two additional genuine triple junctions on that
ray. The same holds on its reverse half-ray if the incident pair
strictly beats the spectator there. An isolated dormant tie does not
supply such an interval.

## Exact germ and uniform radius

Write \(U_i=T_r1_{A_i}\), and let \(L_{ij}\) be the supporting line
of the reference pair \(i,j\). By exterior agreement, the adjusted
pair contrast vanishes on its exterior ray. Analytic continuation
makes the unadjusted contrast constant on all of \(L_{ij}\). On the
forward tail the conditional distribution sees the two incident
half-planes with equal probability up to a vanishing error, so the
constant is zero. Thus all three pair pressures vanish and
\(U_i=U_j\) on \(L_{ij}\).
[Zero pressures and the three full score lines](THREE_CELL_COMPACT_SIMPLEX_DISCREPANCY.md#every-label-must-receive-a-bounded-discrepancy)

Let \(w_i\) be the unit centers of the three reference sectors, with
the cyclic convention \(\ell_{ij}(x)=(w_i-w_j)\cdot x\).
The gradient at zero of each pair contrast has the form
\(\nabla(U_i-U_j)(0)=c_{ij}(w_i-w_j)\), since the contrast vanishes
on \(L_{ij}\). The cyclic sum of the three gradients is zero, and
two edge vectors of the regular triangle are independent. Hence all
three coefficients equal a single scalar \(c\). All three scores
tie at the origin, so central maximal-tie gradient separation gives
\(c\ne0\).
[Uniform central gradient separation](THREE_CELL_CENTRAL_INTERFACE_COMPACTNESS.md#uniform-separation-at-every-maximal-tie)

For \(x\) near zero, let \(z\) be its orthogonal projection to
\(L_{ij}\), and put \(n=(w_i-w_j)/|w_i-w_j|\). Then
\[
 (U_i-U_j)(x)=(x\cdot n)
     \int_0^1\partial_n(U_i-U_j)(z+u(x-z))\,du .             \tag{CG1}
\]
On a compact positive-correlation band the magnitude of the normal
derivative at zero has a uniform positive lower bound by central
gradient separation. The second derivatives have a uniform bound
because the scores convolve bounded indicators with Gaussian
densities whose variances stay bounded away from zero. Thus a
uniform disk keeps all three factors in the integral of (CG1) of
the same sign as \(c\). It gives the exact reference simplex germ
when \(c>0\), and its antipode when \(c<0\).

## Genuine triple points on the reference lines

Parametrize \(L_{ij}\) by \(s\tau_{ij}\) where \(\tau_{ij}\)
is the unit vector along the forward reference ray, and let \(k\)
be the remaining label. On this line define
\[
 h(s)=U_i(s\tau_{ij})-U_k(s\tau_{ij})
     =U_j(s\tau_{ij})-U_k(s\tau_{ij}).                       \tag{CG2}
\]
The exterior simplex makes \(h(s)>0\) for all sufficiently large
positive \(s\) and \(h(s)<0\) for all sufficiently large negative
\(s\). The limiting scores in these two cases are respectively
\((U_i,U_j,U_k)=(1/2,1/2,0)\) and \((0,0,1)\).

Every interval on which \(h>0\) is an actual pair interface.
Indeed the normal derivative of \(U_i-U_j\) on the line is an
analytic function of \(s\) with nonzero limit on the forward tail
by the binary half-space formula. Its zeros are isolated. Where
it is nonzero and \(h>0\), labels \(i,j\) are strict winners on
the two sides of the line. The reciprocal-gradient bound (G20)
excludes a zero of this normal derivative in the closure of any
such interval. Otherwise local Lipschitz continuity bounds the
gradient magnitude above by a constant times the distance to
that zero, and its reciprocal is not integrable along the
adjacent true pair wall. At a sign-changing zero of \(h\),
the strict winners of both \(i,j\) approach along the positive
side of the interval. Label \(k\) strictly wins on the nearby
part of the line with \(h<0\). The zero is therefore a genuine
triple junction.
[Reciprocal-gradient control and genuine junctions (G20)](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#regular-interface-arcs-and-isolated-triple-junctions)

For the antipodal germ, (CG1) gives \(h<0\) just forward of
zero and \(h>0\) just backward. Comparing to the exterior
signs yields a sign-changing zero on each reference half-ray.
For the simplex germ the signs near zero and at the respective
infinities agree. An excursion to the opposite strict sign
forces at least two sign-changing zeros on that half-ray. An
isolated zero with no sign change can be a dormant equality and
the argument makes no junction claim there.

The same argument yields a finite-arm consequence for a
nonsimplex tie approached from the right by simplex-optimal
correlations. On a compact positive-correlation band there is
a uniform finite length \(L\) such that at least one of the
three active germ half-rays encounters a strict spectator win
before radius \(L\). Otherwise \(h\ge0\) from the origin
through length \(L\) on all three active rays. Each trace is
positive except isolated zeros and therefore contains the
corresponding actual pair arm up to surface-null sets. These
three centered arms of length \(L\) contradict the finite-arm
comparison at a right contact. Before the strict spectator win,
there is a genuine triple junction where the sign changes.
[The finite-arm surface comparison at a right contact](THREE_CELL_ORIGIN_ARM_SURFACE_COMPARISON.md#arcs-close-to-the-origin-and-endpoint-implications)

## Other straight active segments lie on their reference lines

Suppose a nontrivial straight active segment for \(i,j\) lies
on an affine line \(L\). The analytic contrast \(U_i-U_j\)
vanishes on all of \(L\). If the direction of \(L\) is not a
reference line direction, its two ends asymptote to antipodal
directions in the interiors of sectors. They cannot both lie
in the spectator sector, whose opening is \(120^\circ\). At
one end the pair contrast tends to \(1\) or \(-1\), a
contradiction. If the direction is parallel to either of the
other two reference lines, the reverse direction of that
reference ray lies in the interior of sector \(i\) or \(j\),
giving the same contradiction. Thus \(L\) is parallel to
\(L_{ij}\). With signed normal offset \(b\), its forward binary
limit is \(2\Phi(rb/\sqrt{1-r^2})-1\), so \(b=0\).
