# Long regular pair arcs cannot pass close to the origin

Let \(I\subset(0,1)\) be a compact correlation band. There are
\(\delta>0\) and \(L<\infty\), depending only on \(I\), with the
following property. No balanced three-cell planar maximizer at
\(r\in I\) has an injective regular arc of a single pair type,
parametrized by arclength on \([-L,L]\), whose midpoint is within
\(\delta\) of the Gaussian origin. Other interfaces may occur
anywhere off this arc. In particular no proper two-ended
pair arc passes through the origin. The same conclusion holds
for positive prescribed masses if the constants also depend on
a positive lower bound for the third mass.

A whole centered straight line of actual pair interface is thus
excluded. More precisely, suppose a whole straight line at signed
normal displacement \(b\) is a dormant interface. Write \(c_j\) for
the constant unadjusted score of its dormant third label along
the line. Then
\[
 \exp\!\left[-\frac{1-r}{2(1+r)}b^2\right]
 < \frac r3 \exp\!\left[-\frac12
      \bigl(\Phi^{-1}((1+c_j)/2)\bigr)^2\right].
 \tag{PA1}
\]
The variations use the surface measure of an actual pair interface.
A line of score equality alone supplies no such support.

## Comparison by arclength

Set \(\sigma=\sqrt{1-r^2}\), \(\kappa=(1-r)/(2(1+r))\), and use
\(\phi_d\) for standard Gaussian density in dimension \(d\). Orient
the arc from phase \(i\) toward phase \(k\) and write \(a\) for the
negative outward normal derivative of \(U_i-U_k\). If \(Z\) is
standard Gaussian in the plane and \(n\) is any unit vector, then
\[
 \partial_n(U_i-U_k)(x)=\frac r\sigma\mathbb E\bigl[(Z\cdot n)
       (1_{A_i}-1_{A_k})(rx+\sigma Z)\bigr].
 \tag{PA2}
\]
Thus \(a\le 2r\phi_1(0)/\sigma\), and the inequality is strict at
each point because the remaining phase has positive mass. This
uses only the score derivative and the positive mass of the third
label.

There is a uniform improvement of this bound for \(|x|\le2\)
when \(r\in I\) and \(\gamma_2(A_j)\ge\eta>0\). Fix a source
disk whose standard Gaussian complement has measure below
\(\eta/2\). On its product with \(|x|\le2\), the conditional
density relative to standard Gaussian has a uniform positive
lower bound. Hence \(\mathbb P(rx+\sigma Z\in A_j)\ge c>0\)
uniformly. Of events of mass at least \(c\), the centered band
\(|Z\cdot n|\le q_c\), \(q_c=\Phi^{-1}((1+c)/2)\), minimizes
the integral of \(|Z\cdot n|\). This follows by moving mass
from larger to smaller values of \(|Z\cdot n|\). In particular
\(a\le 2r\phi_1(q_c)/\sigma\) on the central disk.

The first extra-coordinate form on a single pair, after a change
from weighted area to Euclidean arclength, uses the kernel
\[
 k_r(x,y)=M_r(x,y)\sqrt{\phi_2(x)\phi_2(y)}
 = \frac{1}{2\pi\sigma^2}
 \exp\!\left[-\frac{\kappa}{2}(|x|^2+|y|^2)
              -\frac{r}{2\sigma^2}|x-y|^2\right].
 \tag{PA3}
\]
For an arclength parametrization \(x(s)\) of the given arc
with \(x(0)=p\), use the speed
\[
 f(x(s))=\sqrt{\frac{\phi_2(se_1)}{\phi_2(x(s))}},\qquad
 -L\le s\le L,
\]
and zero speed elsewhere. This is a bounded compactly supported
measurable speed on the given pair. The local term in (G1) is
\(\int_{-L}^L a(x(s))\phi_2(se_1)\,ds\). When \(L\ge1\) and
\(|p|\le1\), the uniform improved score derivative on \(|s|\le1\)
therefore gives
\[
 \int a f^2\,d\nu\le \Lambda_r-d_I,\qquad
 \Lambda_r=2r\phi_1(0)^2/\sigma,\qquad d_I>0.
 \tag{PA4}
\]
Here \(d_I\) can also depend on the lower mass bound.

The arclength conditions imply \(|x(s)-x(t)|\le|s-t|\) and
\(|x(s)|\le|s|+|p|\). Consequently
\[
 k_r(x(s),x(t))\ge k_r(se_1,te_1)
                  e^{-\kappa |p|^2-\kappa |p|(|s|+|t|)}.
\]
At \(p=0\), on the whole reference line, the \(rK\) interaction
for this test equals \(\Lambda_r\). There \(f=1\), surface measure
is \(\phi_1(0)\,d\gamma_1\), and the Mehler kernel is
\(\sigma^{-1}M_r^{(1)}\) with the same-pair factor two. On truncating
to \([-L,L]^2\) and moving to \(p\), the displayed pointwise lower
bound gives
\[
 r\langle f,Kf\rangle\ge
              \Lambda_r-C_I|p|-C_Ie^{-c_I L^2}.
 \tag{PA5}
\]
Indeed all reference integrands are bounded by an integrable
Gaussian on the compact correlation band, and the displacement
loss is controlled by \(|p|(1+|s|+|t|)\) times this Gaussian.
Choosing \(\delta\) small and then \(L\) sufficiently large
contradicts (G1) by (PA4)–(PA5).
Thus a regular point in the smaller central disk reaches a
genuine junction within arclength \(L\) in at least one direction
along its pair edge, or lies on a junction-free closed component
of length at most \(2L\). Here one uses the local continuation
of regular pair edges up to a junction, a closed curve, or
infinity.
[The first extra-coordinate form](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#translation-equations-on-an-unconstrained-space-of-speeds)
[Continuation of regular pair edges](THREE_CELL_CONNECTED_NETWORK_TOPOLOGY.md#the-local-graph-and-its-three-faces)

## A one-dimensional comparison with missing interface portions

Let \(x:\mathbb R\to\mathbb R^2\) be a \(1\)-Lipschitz path
with \(x(0)=p\), and let \(E\subset\mathbb R\) be measurable.
Suppose that on \(E\) the path is injective and has unit
arclength density almost everywhere, and that \(x(E)\) is part
of a single actual regular pair interface up to arclength zero.
No same-pair condition is imposed on the omitted portions of
the path. Define
\[
 g(s)=\frac{\sigma a(x(s))}{2r\phi_1(0)},\qquad
 w_p(s)=e^{-\kappa |p||s|} .
\]
For every real \(h\in L^2(\gamma_1)\) supported on \(E\),
\[
 \int g(s)h(s)^2\,d\gamma_1(s)
 \ \ge\ e^{-\kappa |p|^2}
        \langle w_p h,T_r^{(1)}(w_p h)\rangle_{\gamma_1} .
                                                               \tag{PA6}
\]
Indeed use the speed \(|h(s)|\sqrt{\phi_2(se_1)/\phi_2(x(s))}\)
on the retained part of the path. First truncate its size and
support. The density conversion in (PA3) and the Lipschitz
inequalities on the entire path give the kernel comparison
used in (PA5). The area formula identifies the retained
arclength integrals with the integrals on \(E\). Divide (G1)
by \(\Lambda_r\), and pass to the limit by monotone convergence
for the nonnegative integrands. Positivity of the reference
kernel permits replacing \(|h|\) by \(h\) on the right. This
proves (PA6).

At \(p=0\), write \(S_r^{(1)}(E)=\langle1_E,T_r^{(1)}1_E\rangle\).
Testing (PA6) with \(h=1_E\) gives the exact restriction
\[
 \gamma_1(E)-S_r^{(1)}(E)
           \ \ge\int_E (1-g(s))\,d\gamma_1(s).            \tag{PA7}
\]
The right side measures the loss from the largest possible
normal derivative on a binary wall. On the fixed central
interval it is bounded below by a positive constant times
\(\gamma_1(E\cap[-1,1])\), uniformly on compact correlation
bands with a positive third mass. The same assertion holds
with an additive error \(C_I|p|\) for \(|p|\le1\), since
the loss from the factor \(e^{-\kappa |p|^2}w_p(s)w_p(t)\)
in the reference expectation is bounded by \(C_I|p|\).
On the other hand, \(S_r^{(1)}(E)\ge\gamma_1(E)^2\) by
positive definiteness on mean-zero functions. In particular
there is a uniformly positive omitted Gaussian mass whenever
\(|p|\) is sufficiently small. This excludes paths that
nearly fill the pair interface on the Gaussian arclength
scale even if the omissions occur in many separate places.

## The offset bound for a dormant straight wall

If the whole straight line is a dormant wall, its third-label
conditional mass is \(c_j\) at every point on the line. Equation
(PA2) and the centered-band optimization above therefore yield
\(a\le 2r\phi_1(q_{c_j})/\sigma\). Along the line at offset \(b\),
\(d\nu=\phi_1(b)d\gamma_1\) and the Mehler normal factor on
the product of the line with itself is
\(\sigma^{-1}\exp(rb^2/(1+r))\). Apply (NT2) with speed one,
using bounded compact-support cutoffs and their limit. Its strict
version persists in the limit because the dormant-label mass flux
of the strip insertion tends to the positive measure \(\phi_1(b)\),
and the positive translation term in (NT1) is a continuous
quadratic form of this flux. The result is
\[
 (2r/\sigma)\phi_1(q_{c_j})\phi_1(b)
 > (6/\sigma)\phi_1(b)^2 e^{rb^2/(1+r)},
\]
which gives (PA1).
[Dormant-wall insertion inequality (NT1)–(NT2)](THREE_CELL_CONNECTED_NETWORK_TOPOLOGY.md#inserting-the-dormant-label-along-a-wall)

## Two separate arms of one pair and the right-contact count

The two halves in (PA4)–(PA5) need not be joined. More precisely,
on each compact interior correlation band there are \(\delta>0\)
and \(L<\infty\) such that no balanced maximizer has two regular
arclength-parametrized arms of a single pair type, each of
length \(L\), starting within \(B_\delta\), with disjoint
arclength interiors. Parametrize them as \(x(s)\) on
\([-L,0)\cup(0,L]\) with starting points \(p_-,p_+\). The bound
\(|x(s)|\le |s|+\delta\) holds on both halves. For two parameters
of the same sign their chord is at most \(|s-t|\), and for
opposite signs
\[
 |x(s)-x(t)|\le |s|+|t|+|p_+-p_-|
                    \le |s-t|+2\delta .
 \tag{PA8}
\]
Use the same surface speed on the two arms as in (PA4). They
can use a single orientation for that pair irrespective of
their tangent directions. The local estimate (PA4) is
unchanged, since \(|x(s)|\le2\) for \(|s|\le1\) and
\(\delta\le1\). In (PA3), the extra chord loss from (PA8) is
at most \(C_I[\delta|s-t|+\delta^2]\) in the negative exponent.
Together with the radius loss this gives (PA5) with \(C_I\delta\)
in place of \(C_I|p|\), using the same reference Gaussian
integrability. Choosing \(\delta\) and \(L\) as before gives
the contradiction.

It follows that at a right-accumulating tie a nonsimplex balanced
maximizer has at most two disjoint sufficiently long regular pair
arms starting sufficiently near the origin, uniformly on compact
interior correlation bands. Indeed if three arms include a
repeated pair, apply (PA8). If their types are distinct, apply
the positive derivative comparison (OA8) and its finite-arm
version.
[Long arms with three distinct pair types and the right-contact derivative sign](THREE_CELL_ORIGIN_ARM_SURFACE_COMPARISON.md#arcs-close-to-the-origin-and-endpoint-implications)

## An annulus crossed by all three strict winners

Let \(\delta,L\) satisfy the right-contact count. In a
nonsimplex right-accumulating tie, suppose that \(0<a<\delta\),
\(b-a>L\), and that each of the three strict-winner sets has a
properly embedded crosscut from the inner to the outer circle
of the annulus \(a\le |x|\le b\). Here the interiors of the
crosscuts lie in the open annulus. Then the closed annulus
contains a genuine triple junction.

Suppose otherwise, and first take the two circles transverse to
the jump graph \(G\). All of \(G\) in the annulus is regular, with
locally finitely many analytic branches. By compactness and
transversality the intersections with the circles are finite,
and the components in the closed annulus are finitely many
disjoint regular arcs and loops. A crosscut in a strict-winner
set is disjoint from \(G\), including at its endpoints.

The three disjoint winner crosscuts divide the annulus into
three topological rectangles. Each has two opposite sides
on the inner and outer circles and its other two sides on
crosscuts of distinct labels. In any such rectangle at least
one of its jump arcs joins the inner and outer sides. To see
this, if there were none, every arc would cut off a Jordan cap
on either the inner or the outer side. The remaining jump
components would be Jordan loops in the interior. The cap
closures and loop interiors with their boundaries taken
together are either disjoint or nested by Jordan separation.
The finitely many maximal such regions can be avoided by a
path between the two winner sides: detour around each closed
region along its arc or loop on the outside. These detours
can be taken in disjoint neighborhoods and away from the
inner and outer sides of the rectangle. Since the jump graph
is disjoint from the resulting path, its endpoints must have
the same label. This contradicts the choice of the winner
sides.

The three spanning jump arcs are disjoint and have constant
pair type on each arc, since there is no genuine junction.
Each starts in the smaller central disk and has arclength at
least \(b-a>L\). The right-contact count excludes these arcs.

The transversality assumption entails no additional restriction.
With no genuine junction in the closed annulus there is a
slightly larger junction-free annulus. Sard's theorem on the
regular edges gives transverse choices of both nearby radii.
The strict-winner crosscuts can be extended a short distance
at their endpoints within their open winner sets. Taking
subarcs at the new circles preserves the crosscuts and the
strict inequalities \(a<\delta\), \(b-a>L\).
[Regularity, local finiteness, and face labels of the jump graph](THREE_CELL_CONNECTED_NETWORK_TOPOLOGY.md#the-local-graph-and-its-three-faces)
