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

## The transverse source in the pair-wall inequality

Suppose the scores of labels \(i,j\) agree on their full centered
reference line. Use coordinates \((t,v)\) along and normal to the
line. Let \(q(t)\) be the spectator's transverse mass with respect to
\(V\sim N(0,\sigma^2)\), and let \(P_s\) denote expectation in
\(t\sim N(rs,\sigma^2)\). The tilted slice identity (CD1) makes the
other two transverse masses both \((1-q(t))/2\) for almost every \(t\).
Set
\[
 \Gamma(u)=\frac{\phi_1(\Phi^{-1}((1+u)/2))}{\phi_1(0)},\qquad
 m(t)=\frac{\mathbb E[V(1_{A_i}-1_{A_j})(t,V)]}
                    {2\sigma\phi_1(0)} .
\]
Put \(\Gamma(1)=0\) by continuity. For a portion \(E\) of the line
that is an actual regular \(i,j\) interface, the normalized slope
in (PA6) is
\[
                  g(s)=|P_sm|\le P_s\Gamma(q).
                                                            \tag{PA9}
\]
The equality is Gaussian differentiation in the observer's normal
coordinate. The inequality follows by placing the two equal-mass
transverse sets in the opposite Gaussian tails to maximize the
absolute signed first moment.
[The exact tilted slice identity (CD1)](THREE_CELL_COMPACT_SIMPLEX_DISCREPANCY.md#every-label-must-receive-a-bounded-discrepancy)

There is a quantitative version. For each \(s\) choose a sign
\(\epsilon_s\) with \(\epsilon_s P_sm=|P_sm|\). Let
\(\Delta_i(t;s)\) and \(\Delta_j(t;s)\) be the \(N(0,\sigma^2)\)
masses of the symmetric differences from the two tails that
maximize \(\epsilon_s m(t)\). The upper tail is assigned to \(i\)
when \(\epsilon_s=1\), and to \(j\) when \(\epsilon_s=-1\).
Then
\[
 \Gamma(P_sq)-g(s)\ \ge\
 \frac{\pi}{4}\left(\operatorname{Var}_{P_s}q+
            P_s[\Delta_i(\,\cdot\,;s)^2+\Delta_j(\,\cdot\,;s)^2]\right).
                                                            \tag{PA10}
\]
Indeed, writing \(z=\Phi^{-1}((1+u)/2)\) gives
\(\Gamma''(u)=-1/(4\phi_1(0)\phi_1(z))\le-\pi/2\).
Thus the loss from Jensen's inequality is at least
\((\pi/4)\operatorname{Var}_{P_s}q\), including endpoint values by
continuity. For the tail estimate, suppose a set differs from a
tail of the same mass by symmetric difference \(\Delta\). On each
side of the tail threshold there is a displaced mass \(\Delta/2\).
Since the Gaussian density is at most \(D=\phi_1(0)/\sigma\),
its distance integral from the threshold is at least
\((\Delta/2)^2/(2D)\) on each side. This is the integral obtained by
placing the mass at density \(D\) immediately next to the threshold.
The first-moment loss is therefore at least \(\Delta^2/(4D)\).
Sum over the two labels, divide by \(2\sigma\phi_1(0)\), and
average with \(P_s\). This proves (PA10).

Substituting (PA9) into (PA6) yields the source-only condition
\[
 \int_E P_s\Gamma(q)\,h(s)^2\,d\gamma_1(s)
       \ge \langle h,T_r^{(1)}h\rangle_{\gamma_1},
       \qquad h=1_Eh.
                                                            \tag{PA11}
\]
For \(h=1_E\), the right side of (PA7) is at least
\(\int_E[1-\Gamma(P_s q)]\,d\gamma_1(s)\) plus the integral of the
two terms in (PA10). Thus the one-pair test measures variance of
the unsmoothed spectator profile and departure from the same
transverse ordering on nearby slices.

## Source oscillations with a strict one-pair gap

For each fixed \(0<r<1\), ordinary balanced partitions can differ
from the simplex on arbitrarily small compact sets, satisfy all
three tilted slice identities, and have arbitrarily many source
sign changes for all three spectators. Their scores choose the
reference simplex and satisfy a strict (PA6) on each reference
winner ray using their own score slopes. The changed cells fail
the winner equation, so these tests are on score interfaces.

For spectator \(k\), orient \(t\) opposite its reference center and
choose \(v\) so that near the upper spectator boundary the cone is
\(v<-\sqrt3t\). The reference profile for \(t<0\) is
\(q_T(t)=2\Phi(-\sqrt3 t/\sigma)-1\). Its threshold \(q_T=1/3\)
occurs at \(t_*=-\sigma\Phi^{-1}(2/3)/\sqrt3<0\). For small
\(d>0\), choose any finite sequence of disjoint ordered intervals
\(I_\ell\) in \(|t-t_*|<d^2\). At each interval use the box
\[
 t\in I_\ell,\qquad v+\sqrt3t\in[d,11d/10]
                    \quad\text{or}\quad[-11d/10,-d],
\]
alternating the choice. Exchange the labels of the two incident
cones on each box and on all its images under the six symmetries
of the reference simplex, acting on points and labels together.
The boxes lie strictly on one side of the upper ray between \(k\)
and a neighbor. Reflection in that ray takes
\(t\) to \(t-(\sqrt3/2)(v+\sqrt3t)\). Thus its image near the same
ray misses every \(I_\ell\). Reflection in the \(t\)-axis gives
the lower copy and preserves \(k\) and \(t\). The remaining two
images are near the pair ray of the two other labels and do not
change the spectator. These copies are disjoint up to null sets.
On each of the two spectator rays, the unreflected boxes are
disjoint by their \(t\)-intervals, and the reflected boxes are
their isometric images with disjoint \(t\)-range from the originals.
On the remaining ray the same disjointness follows by rotation.

The resulting assignment is equivariant, so its Gaussian masses
are equal and reflection gives all three slice identities.
For \(t\in I_\ell\), it changes \(q_T(t)\) by a quantity of
alternating sign and magnitude at least \(c_r d\), from the two
copies at the same \(t\). The reference error
\(|q_T(t)-1/3|\) is \(O_r(d^2)\), so the new profile alternates
sides of the threshold on the successive intervals. The three
sources are rotations of this example.

The total changed area is \(O(d^3)\) since the intervals have
summed length at most \(2d^2\). Each score difference and its
first derivatives change uniformly by \(O_r(d^3)\) by convolution
with the Gaussian kernel. The score difference perturbation
vanishes on the corresponding full pair line. Integrating the
normal derivative from the line sharpens the bound to
\(O_r(d^3)\min\{\operatorname{dist}(x,L_{ij}),1\}\).
The saturated reference margin on each winning sector then
preserves its winner.
[Saturated reference score margin](POSITIVE_CORRELATION_CONTINUATION.md#relative-exterior-control-below-every-fixed-correlation-bound)

At the reference simplex, the difference between the left and
right sides of (PA6) for a single ray is bounded below by
\(c'_r\|h\|_{L^2(\gamma_1)}^2\) with \(c'_r>0\).
Apply (S37) in one extra coordinate \(z\) to the balanced speed
\(zh(s)\) supported on that ray. The rotations in its degree-one
block have the form \(z\,b\cdot N\) for a fixed planar vector \(b\).
The three rays have equal Gaussian surface mass, and their unit
normals satisfy \(\sum_e N_e N_e^{\mathsf T}=(3/2)I\).
At most two thirds of the squared norm of such a rotation lies
on one ray. Cauchy–Schwarz then bounds the squared projection of
the test speed onto rotations by two thirds of its own squared
norm. The normalization in (PA6) divides the corresponding
first extra-coordinate form by the positive constant \(\Lambda_r\).
This proves the strict one-pair bound. Since the perturbed slopes
\(g\) differ uniformly by \(O_r(d^3)\), the same bound with a
smaller constant holds on their score rays for small \(d\).
Thus the exact slice identities, reference score winners, and
the entire family of one-ray tests permit arbitrarily many
transverse source oscillations.
[Three-cell simplex gap modulo rotations (S37)](SIMPLEX_INTERFACE_SYMMETRY.md#reflected-kernels-close-the-three-cell-spectral-problem)

## A strict one-pair gap even with spectator excursions

For each fixed \(0<r<1\), the tests on the individual score
lines also permit genuine excursions of the spectator score
on all three forward rays. Take the partition \(A^R\) equal
to the antipodal simplex on \(R<|x|<2R\) and to the reference
simplex elsewhere. Its mass is equal for each label and it
has all three full pair-score equalities, by the rotational
and reflection symmetries of the reference simplex. Along each
forward ray at \(s=3R/(2r)\), the conditional Gaussian lies
with probability tending to one in the annulus and in the
central angular part of the antipodal spectator sector.
Thus the score difference \(h(s)=U_i(s,0)-U_k(s,0)\) tends
to \(-1\) there. It is positive near zero on this ray for
large \(R\) and eventually positive again because the defect
is compact.

Let \(E_R=\{s:h(s)>0\}\) on one full line, and on \(E_R\)
use its normalized pair slope \(g_R\) from (PA9). There is
a fixed \(c>0\) with \(g_R\ge c\) on \(E_R\) for all
sufficiently large \(|s|\), independently of \(R\). To see
this, write \((t,V)\) for the conditional Gaussian at
\((s,0)\). The event \(G_s=\{st>0,\ |V|<\sqrt3|t|\}\)
has probability tending to one as \(|s|\to\infty\). On
this event, the partition is either in the spectator phase
or uses the pair split with signed contrast
\(\operatorname{sign}(V)\) for \(s>0\), after orienting
the transverse coordinate. For \(s<0\) the sign is
reversed. Indeed for positive \(t\) the two alternatives
inside \(G_s\) are the antipodal spectator sector and the
reference split. For negative \(t\) they are the reference
spectator sector and the antipodal split.

Since \(U_k(s,0)<1/3\) on \(E_R\), the nonspectator part
of \(G_s\) has probability at least \(1/2\) when \(|s|\)
is large. Its integral of \(|V|\) has a positive uniform
lower bound, by discarding \(|V|\le\sigma/4\), an event
of probability less than \(1/4\). The integral of
\(|V|\) on \(G_s^c\) tends to zero by Cauchy–Schwarz.
After subtracting this error and applying (PA9), the
uniform bound on \(g_R\) follows.

On a fixed observation window \([-M,M]\), the scores
and their derivatives converge uniformly to the reference
scores as \(R\to\infty\), since modifications are outside
the growing source disk. The zero at \(s=0\) remains exact
by the three pair equalities. The reference has a simple
zero there and strict signs on both open half-rays, so for
large \(R\) the set \(E_R\cap[-M,M]\) is \((0,M]\).
The slopes converge uniformly there.

The reference strict single-ray bound above applies to
speeds supported in \((0,M]\). For a speed \(h_0\)
supported on \(E_R\), decompose it into its restrictions
to \([-M,M]\) and its complement. The Mehler pairings
involving the complement have total absolute value at most
\(2\|T_r^{(1)}1_{\{|s|>M\}}\|\|h_0\|_2^2\), where the
indicator denotes multiplication on the operator's right.
This operator norm tends to zero as \(M\to\infty\) by
compactness of \(T_r^{(1)}\) on \(L^2(\gamma_1)\), whose
Hermite eigenvalues are \(r^j\). Choose \(M\) first so
this error is smaller than half the smaller of \(c\)
and the strict reference gap divided by two. Then take
\(R\) large enough to retain half the reference gap in the
window. The positive tail slope and the retained reference
gap prove the strict (PA6) on all of \(E_R\), for every
line.

These tests are still on score interfaces. For example,
choose \(1<\lambda<\min\{2,1/r\}\). Near the point
\(\lambda R\tau_k\), the assigned label is the antipodal
spectator \(k\), while the conditional Gaussian is
concentrated near the forward reference ray at radius
\(r\lambda R<R\). Its spectator score tends to zero
there. Hence this assignment violates its own winner
equation on an open set despite the strict score
excursions and the single-ray bounds. The excursions in
this family escape every fixed observation window as
\(R\to\infty\).

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
