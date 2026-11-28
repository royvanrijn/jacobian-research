# Thin folds on all three arms of a compact theta

Let \(A\) be a balanced planar Gaussian noise-stability maximizer at
\(0<r<1\). Suppose its entire actual interface is a compact theta
consisting of three regular pair arcs with common junctions \(p,q\).
Two strict-winner cells are bounded and the third is unbounded.
Choose orthonormal coordinates \((x,y)\) through the Gaussian origin,
with \(x(p)=\alpha<\beta=x(q)\). Suppose all three arms lie in the
closed strip \(\alpha\le x\le\beta\). On almost every vertical fiber,
suppose the highest and lowest crossings of each individual arm
are at most \(\sqrt{1-r^2}/10\) apart. Then this theta admits a
mass-preserving planar speed with negative second variation.

The hypothesis allows backtracking on every arm. The width concerns
each arm separately, leaving the distances between arms unrestricted.
When both exterior arms are graphs, the middle arm can have a larger
width.
[Thin folds with graph exterior arms](THREE_CELL_THETA_THIN_FOLD_STABILITY.md)

## Signed crossings of each arm

Orient every arm from \(p\) to \(q\), with unit tangent \(\tau\) and
normal \(N=J\tau\), where \(J\) is counterclockwise rotation by
\(\pi/2\). Its vertical translation speed \(s=e_y\cdot N\) is
\(dx/d\ell\) in arclength \(\ell\). Close each arm separately by a
simple path from \(q\) to \(p\) that meets the arm only at its
endpoints. The path uses the two open exterior half-planes and a
right-to-left horizontal bridge below the theta. This construction
still works if the arm touches a strip boundary away from its junctions.
On almost every interior vertical fiber the original arm has
finitely many transverse crossings by the arclength area formula
and regularity. Along the closed Jordan curve the crossing signs
alternate, since each transverse crossing switches between the
interior and exterior. The bridge furnishes the bottom crossing
with negative sign, leaving \(+,-,\ldots,+\) on the original arm.
The signed fiber measure given by the area formula for
\(s\,d\mathcal H^1\) therefore has these same signs.

Assign the incidence vector \(d_e=e_i-e_j\) on an arm \(e\) when \(N\)
points from cell \(i\) to cell \(j\). Each bounded cell borders two
arms, and the Euclidean vertical integral of \(s\) on each arm
is \(\beta-\alpha\). The outward vertical flux of a bounded cell
is zero, so the two arms have opposite incidence signs at that
label. The third coordinate of \(\sum_e d_e\) is also zero because
its coordinate sum is zero. As each vector has squared norm two,
\(d_e\cdot d_f=-1\) for distinct arms.

## Pairing two signed fibers

Scale the vertical coordinates by \(\sigma=\sqrt{1-r^2}\). The
joint one-dimensional Gaussian density is, up to a positive factor,
\[
             G(y,z)=\exp\!\bigl(-(y^2+z^2-2r yz)/2\bigr).
\]
Consider two alternating signed crossing measures, each starting
and ending with positive atoms. Let their intervals have scaled
lengths \(\delta,\eta\le 1/10\). Choose a corner of the product rectangle
where \(G\) has least value among its corners. Use inward coordinates
\(0\le u\le \delta\), \(0\le v\le \eta\) from that corner. Reversal preserves
the alternating sign order. Expanding the exponent gives
\[
  G=G(0,0)U(u)V(v)e^{\varepsilon r uv},\qquad
  U(u)=e^{p_0u-u^2/2},\quad V(v)=e^{q_0v-v^2/2},\quad
  \varepsilon\in\{-1,1\}.
\]
Comparison with the two corners along the axes gives
\(p_0\ge \delta/2\), \(q_0\ge \eta/2\) when the respective lengths are
positive.

For the alternating measure on \([0,\delta]\), write \(L\) for its
functional. Its total mass is one and its cumulative mass lies
in \([0,1]\). Integration by parts gives
\[
  L(F)=F(0)+\int_0^\delta H(u)F'(u)\,du,\qquad 0\le H\le1 .
\]
Then \(L(U)>0\) and for each integer \(n\ge1\)
\[
              |L(Uu^n)|\le 4(n+1)\delta^n L(U).            \tag{TA1}
\]
If \(\delta=0\), evaluation at zero proves both statements.
If \(p_0\ge \delta+1/\delta\), then \(U'\ge U/\delta\), so \(L(U)\ge1\) and
\(0\le(Uu^n)'\le(n+1)\delta^n U'\).
If \(p_0<\delta+1/\delta\), then
\[
 U\le e^{1+\delta^2},\qquad
 \int_0^\delta |U'|\le(1+\delta^2)e^{1+\delta^2}.
\]
For \(p_0\ge \delta\) the negative variation vanishes. For \(p_0<\delta\) it
is at most \(\delta^2 e^{\delta^2}\) because \(U\le e^{\delta^2}\) and
\(U'\ge-\delta U\). Hence \(L(U)\ge1-\delta^2e^{\delta^2}>0.98\) in both
cases. The product derivative bounds
\[
       |L(Uu^n)|\le(2+\delta^2)e^{1+\delta^2}\delta^n<6.1\delta^n,
\]
which implies (TA1). The same bounds apply to \(V\) and the second
functional \(L'\).

The product pairing of \(G/G(0,0)\) has leading term \(L(U)L'(V)>0\)
from the exponential series in \(uv\). With \(t=r\delta\eta\le1/100\),
(TA1) bounds the ratio of the sum of the absolute values of the
remaining terms to the leading term by
\[
  16\sum_{n\ge1}\frac{(n+1)^2t^n}{n!}
     =16\bigl(e^t(t^2+3t+1)-1\bigr)
     \le\frac{16(t^2+4t)}{1-t}<1.
\]
Here \(e^t\le(1-t)^{-1}\). The pairing of the two signed fibers
is therefore strictly positive. The estimate is uniform in the
locations of the intervals and the number of crossings.

## Negative planar speed

For distinct arms \(e,f\), integrate the signed fiber pairings
against the positive joint Gaussian density in the horizontal
coordinates, using the arclength area formula in both factors.
The surface pair measure is the joint Gaussian density times
\(d\mathcal H^1(p)d\mathcal H^1(q)\), and the joint density
factors over the horizontal and vertical coordinates. The scalar
pairing of \(s1_e\) and \(s1_f\) is consequently strictly positive.
Their incidence factor is \(-1\), so
\(\langle s1_e,K(s1_f)\rangle<0\) for the interface operator \(K\).

Take \(g=s1_e\) on any one arm. With \(a\) denoting the surface
gradient weight, the translation identity \(a s=rKs\) gives
\[
   Q^{(1)}(g)=\int a g^2\,d\nu-r\langle g,Kg\rangle
             =r\langle g,K(s-g)\rangle<0.
\]
Let \(P\) be the \(L^2(a\,d\nu)\) projection onto the translations
and put \(h=g-Pg\). Write \(B\) for the volume-flux map and
\(C:v\mapsto B(v\cdot N)\), which is invertible by full active
moment rank. Then \(t_0=(C^{-1}Bh)\cdot N\) has \(Bt_0=Bh\).
The translation identity and \(h\perp t_0\) in \(L^2(a\,d\nu)\) give
\[
   Q_0(h-t_0)=r^{-1}Q^{(1)}(g)
                   -(r^{-1}-1)(\|h\|_a^2+\|t_0\|_a^2)<0 .
\]
The speed \(h-t_0\) is balanced. Cutting off near the junctions
and correcting the resulting small flux error on regular arc
patches preserves the strict sign.
[Translation identity and stability](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#translation-equations-on-an-unconstrained-space-of-speeds)
[Planar flux repair](THREE_CELL_THETA_THIN_FOLD_STABILITY.md#a-strictly-negative-mass-preserving-second-variation)

## A positive measure of thick fibers

If the network lies in a fixed disk of radius \(R\) and
\(r\) belongs to a compact interval \(I\subset(0,1)\), at least
one arm has crossing spread greater than \(\sigma/10\) on a set
of base abscissae with Gaussian measure at least \(c(R,I)>0\).
Indeed the two bounded cells have combined mass \(2/3\) and
lie in this strip, so its projection interval has Gaussian
measure at least \(2/3\).

On two fibers satisfying the width condition, the lower bound
from (TA1) and the series estimate is uniformly positive.
Use \(L(U),L'(V)>0.98\) and the positive minimum of the joint
Gaussian density on the fixed disk product. On any two fibers,
the absolute scalar pairing is uniformly bounded independently
of the crossing count. Apply
\(L(F)=F(0)+\int H F'\) in each coordinate, now without a width
restriction. The constant, first-derivative and mixed-derivative
terms of the Gaussian are bounded on the compact coordinate
product, and the integration lengths are at most \(2R\).
These upper and lower bounds persist after including the positive
horizontal joint density relative to product Gaussian measure
on the fixed base interval.

Let \(D_e\) be the set of fibers of arm \(e\) violating the
width bound. The scalar pairing of arms \(e,f\) is bounded below
by
\[
 m\bigl(2/3-\gamma_1(D_e)\bigr)\bigl(2/3-\gamma_1(D_f)\bigr)
       -M\bigl(\gamma_1(D_e)+\gamma_1(D_f)\bigr)
\]
if both exceptional masses are less than \(2/3\), with
\(m,M>0\) depending only on \(R,I\). In particular all scalar
pairings are strictly positive if every exceptional mass is
sufficiently small. The preceding negative variation gives
the assertion.
