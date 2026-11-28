# Angular compatibility of two isolated ends

The two ends and connectedness of their common strict winner already
force the centered simplex by a conditional-score comparison. The
argument below also proves a bounded-slope estimate for a single
backward graph using the translation and stability operators.
[Rigidity without a third graph](THREE_CELL_TWO_END_CONDITIONAL_RIGIDITY.md)

Fix \(0<r<1\) and a balanced planar Gaussian noise-stability
maximizer. Suppose two distinct pairs have full-half-plane bounded
isolated graph ends in strictly obtuse directions. Suppose the strict
winner of their common label is connected. Let \(\theta\) be the
angle between the two outward directions. Then
\[
               2\pi/3\ \le\ \theta\ \le\ 3\pi/4 .             \tag{TA1}
\]
If the third interface in the backward bisector half-plane is
eventually one continuous graph, the partition is a centered simplex.
No height or slope bound on that third graph is required. The proof
also gives two directional constraints on the noncommon labels
without assuming any third graph.

## The common phase and the two slice constraints

The two ends are exactly straight by shifted-contour straightening,
and their supporting lines pass through the origin by the
two-end pressure argument. All pair pressures vanish. Choose the
positive horizontal bisector and put
\[
 a=\theta/2,\quad c=\cos a,\quad d=\sin a,\quad
 u_+=(c,d),\quad u_-=(c,-d),\quad n_+=(-d,c).
\]
Denote the common label by \(0\), and the noncommon labels of the
upper and lower ends by \(+\) and \(-\), respectively.
[Exact straightening of bounded isolated ends](THREE_CELL_TAIL_FLATNESS.md#reaching-the-straightness-threshold)
[Two isolated ends and the zero pressures](THREE_CELL_ONE_STRAIGHT_END_REFLECTION.md#two-isolated-ends-and-a-second-order-tie)

The connected strict winner of label \(0\) is confined to the
open wedge between \(u_+,u_-\). Indeed it avoids both supporting
lines, and it wins near infinity in the overlap of the two forward
half-planes, where it is the only label present in both. It agrees
with this wedge outside a disk up to null sets. Every point of the
wedge sufficiently far out is in at least one of the two defining
forward half-planes, since its projection onto one outward
direction is at least \(c\) times its norm. In that half-plane
the wedge is on the common side of the exact interface. The
wedge's Gaussian mass is \(\theta/(2\pi)\), so balance gives the
lower bound in (TA1). At equality the common cell fills a regular
simplex sector and the single-cell rigidity result gives the
whole simplex.
[One simplex cell determines the whole partition](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#one-simplex-cell-determines-the-whole-partition)

Write \(x=s u_+ + v n_+\) and \(\sigma=\sqrt{1-r^2}\).
Set \(\kappa=\sin\theta/(-\cos\theta)>0\).
For all sufficiently negative \(s\), the common cell on the
slice is precisely \(v<-\kappa|s|\): this is the wedge on the
slice and it lies far outside the exceptional disk. The forward
half-plane of the lower end excludes label \(+\) whenever
\[
 v < \frac{s\cos\theta-S_-}{\sin\theta}
                =\frac{|s|}{\kappa}-\frac{S_-}{\sin\theta},
                                                               \tag{TA2}
\]
where \(S_-\) is its longitudinal cutoff. Here we used
\(u_-\cdot x=s\cos\theta-v\sin\theta\).
The smoothed \(0,+\) contrast vanishes on the whole line \(v=0\)
by analytic continuation from the straight tail. Injectivity of
the longitudinal Gaussian convolution, as in (PG10), gives
\[
       \int_{\mathbb R}1_{A_+}(s,v)p_\sigma(v)\,dv
             =\Phi(-\kappa |s|/\sigma)
                           \qquad\text{for a.e. late }s .   \tag{TA3}
\]
The same holds with indices exchanged and coordinates reflected
across the horizontal bisector for the other end. Comparing (TA3)
with the largest possible mass allowed by (TA2) forces
\(\kappa\ge1/\kappa\) by the Gaussian tail rate, hence
\(\theta\le3\pi/4\).
[The Gaussian slice identity on a straight line, PG10](THREE_CELL_PRESSURE_GRADIENT_TAIL.md#exterior-geometry-when-all-three-ends-straighten)

## Directional concentration and an exclusion for the third graph

For each \(0<\varepsilon<\kappa\), label \(+\) occurs on
positive transverse measure in the band
\[
    (\kappa-\varepsilon)|s|\ \le v\le\
    (\kappa+\varepsilon)|s|
                   \qquad\text{for a.e. sufficiently negative }s .
                                                               \tag{TA4}
\]
In fact it cannot be a strict winner for
\(|v|\le(\kappa-\varepsilon)|s|\) at large negative \(s\).
To see this, condition the source given the observation \((s,v)\).
Its mean is \((rs,rv)\) and its coordinate variance is
\(\sigma^2\). On the source box
\(\{|s'-rs|\le\delta|s|,\ |v'-rv|\le\delta|s|\}\), choose
fixed \(\delta>0\) small enough that \(s'<0\) and
\(|v'|\le(\kappa-\varepsilon/2)|s'|\).
At almost every \(s'\) there, (TA3) bounds the Lebesgue
length of the \(+\) slice in this box by
\[
 \frac{\Phi(-\kappa |s'|/\sigma)}
      {p_\sigma((\kappa-\varepsilon/2)|s'|)}
                         \le C e^{-C'|s|^2},
\]
because the denominator is a lower bound on its Gaussian
density. The conditional transverse density is uniformly
bounded, and the conditional probability of leaving the
box tends to zero at a Gaussian rate. Thus \(U_+(s,v)\to0\)
uniformly in the stated region. With all pressures equal,
\(\max_i U_i\ge1/3\) excludes \(+\) there. At negative
\(v\) outside this central region it is also excluded
by (TA2). Finally the slice \(v>(\kappa+\varepsilon)|s|\)
has Gaussian mass strictly smaller than the positive mass
in (TA3), which proves (TA4). The same reasoning localizes
label \(-\) in the reflected band along the lower end.

In the backward horizontal half-plane the common label is
absent. The \(+\) and \(-\) labels are respectively the
upper and lower regions at large positive and negative
vertical coordinates by the two exact ends. Suppose their
interface there is a single graph \(y=h(x_1)\) for \(x_1\)
sufficiently negative. If \(\theta>2\pi/3\), then
\(\kappa<d/c\), so \(D=d-\kappa c>0\). Apply (TA4)
to the upper and lower slices with the same longitudinal
parameter \(s=-L\), using the intersection of their two
full-measure sets of \(L\). In ordinary horizontal and
vertical coordinates this gives a \(+\) point and a
\(-\) point both within \(d\varepsilon L\) horizontally of
\(-(c+\kappa d)L\). Their vertical coordinates are at most
\(-DL+c\varepsilon L\) for the \(+\) point and at least
\(DL-c\varepsilon L\) for the \(-\) point. Thus the values of
\(h\) at those two abscissae differ by at least \(DL\) for
sufficiently small fixed \(\varepsilon\). Their abscissae
differ by at most \(2d\varepsilon L\).

This is impossible with a fixed Lipschitz bound, since
\(\varepsilon\) can be made arbitrarily small before taking
\(L\) large. More generally, it is impossible if the
rescaled heights are equicontinuous at relative longitudinal
scales on the tail, meaning
\[
 \lim_{\delta\downarrow0}\ \limsup_{x_1\to-\infty}
   \ \sup_{|x'_1-x_1|\le\delta|x_1|}
       \frac{|h(x'_1)-h(x_1)|}{|x_1|}=0 .
\]
Under either hypothesis (TA1) forces \(\theta=2\pi/3\),
which gives the simplex as above. Without such a hypothesis
the slice proof imposes the directional concentration. The
translation equation supplies the slope bound whenever there is
a single graph.

## A third graph has bounded slope

Assume that the entire third interface for sufficiently negative
horizontal coordinate is one continuous graph. Every point of this
graph is a regular pair point, since the common label has no strict
winners in the backward half-plane. Orient the graph from the upper
label \(+\) to the lower label \(-\), and write its unit normal as
\(N=(f,-w)\), so \(w\ge0\). Vertical tangencies are allowed at this
stage. On every open portion of its parameter interval there is a
point with \(w>0\), since the horizontal coordinate of a graph cannot
be constant on an open arc.

The graph cannot enter either defining forward half-plane, since
one of its two labels is absent there. In particular, for any fixed
\(C\) with \(c/d<C<1\) it is contained sufficiently far back in
\[
                  x_1<0,\qquad |x_2|\le C|x_1| .          \tag{TA5}
\]
The two forward half-planes together with a backward horizontal
half-plane cover the complement of a compact set. Indeed, failure
of the two forward inequalities bounds both vertical sides whenever
the horizontal coordinate stays bounded below, and their sum bounds
that coordinate above. Thus outside this compact set the only
interfaces are the graph and the two exact straight rays.

Let \(\Gamma\) be the part of the graph to the left of a sufficiently
late cut, to be chosen below. The vector translation identity on
\(\Gamma\) gives
\[
 a(p)N(p)=2r\int_\Gamma M_r(p,z)N(z)\,d\nu(z)+F(p),       \tag{TA6}
\]
where \(a\) is the positive normal score derivative, \(d\nu\) is
Gaussian surface measure and \(F\) is the contribution from the
complementary interfaces. Orient the upper straight ray from \(+\)
to the common label and the lower straight ray from the common label
to \(-\). Their normals are \((d,-c)\) and \((-d,-c)\), respectively.
Both incidence dot products with the graph incidence are one, so
their contribution is
\[
           F^{\mathrm{ray}}=r(d(J_+-J_-),-c(J_++J_-)),
           \qquad J_\pm(p)=\int_{\mathrm{ray}_\pm}M_r(p,z)d\nu(z).
\]
It satisfies \(|F_x^{\mathrm{ray}}|\le(d/c)(-F_y^{\mathrm{ray}})\).
[Translation identity (G2)](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#translation-equations-on-an-unconstrained-space-of-speeds)

Fix an initial backward cut at a point with \(w>0\), beyond which the
graph description and (TA5) hold. The complementary interfaces other
than the two ray tails form a fixed compact portion \(E\). For any
later cut at a point with \(w>0\), the remaining contribution to
\(F\) comes from \(E\) and the compact graph segment \(I\) between
the two cuts. On \(I\) there is no sign change in the downward normal
component. We prove that its total forcing obeys a cone bound
\[
                    |F_x(p)|\le C_1(-F_y(p))              \tag{TA7}
\]
on the whole later tail, for a finite constant \(C_1\).

For a compact surface piece, write the conditional boundary
integral as
\[
  M_r(p,z)\,d\nu(z)=
   \frac{e^{-r^2|p|^2/(2\sigma^2)}}{2\pi\sigma^2}
   e^{(r/\sigma^2)p\cdot z}e^{-|z|^2/(2\sigma^2)},ds(z),
   \qquad \sigma^2=1-r^2 .                               \tag{TA8}
\]
Let \(\mathcal C\) be the compact set of unit directions in the
closed cone (TA5). A fixed compact graph arc at sufficiently
negative horizontal coordinate with \(w\) bounded below can be
chosen with projection \(q\cdot z\) strictly larger than every
projection of \(E\), uniformly in \(q\in\mathcal C\). Indeed a
graph point \(z\) in (TA5) has
\(q\cdot z\ge |q_1||z_1|(1-C^2)\), and \(|q_1|\) is bounded below.
Use a compact nonvertical patch about such a point. Formula (TA8)
makes the norm of the forcing from \(E\) exponentially smaller,
uniformly in \(q=p/|p|\), than the positive downward contribution
of this patch for large \(|p|\). In particular, a fixed threshold
ensures it is at most half that downward contribution. Choose
the later cut beyond the patch and threshold and at a nonvertical
point. Then the patch lies in \(I\), supplies a positive downward
contribution to all observation points on the tail, and half
of it suffices to dominate \(E\).

For the now fixed compact segment \(I\), the ratio
\(\int_I M_r(p,z)|f(z)|d\nu(z)/\int_I M_r(p,z)w(z)d\nu(z)\)
is bounded uniformly for \(p\) in (TA5). In fact for any
\(q\in\mathcal C\), a maximum of \(q\cdot z\) on \(I\) cannot
occur at a vertical tangent. The two endpoints are nonvertical.
At an interior vertical tangent with \(q_2\ne0\) the tangent
derivative of the projection is nonzero, and for \(q_2=0\) the
maximum is at the left endpoint. By compactness there is a
uniform positive projection gap between the maximum and points
where \(w\) is below a sufficiently small fixed positive number.
A finite cover of \(\mathcal C\) gives near-maximizing patches
with \(w\,ds\) mass uniformly bounded below. Formula (TA8)
therefore bounds the contribution of the small-\(w\) portion
relative to the vertical integral at large \(|p|\). On the
remaining portion \(|f|/w\) is bounded pointwise. At bounded
\(|p|\) the same ratio of integrals is bounded by continuity
and strict positivity of the denominator. This proves the
graph segment's cone bound. Combining it with the ray bound
and the domination of \(E\) proves (TA7).

In the weighted space \(L^2(a\,d\nu|_\Gamma)\) let \(B\) denote
the operator with kernel \(2r M_r(p,z)/a(p)\) against \(d\nu(z)\).
It is positive semidefinite and preserves nonnegativity. The
first extra-coordinate form satisfies
\(Q^{(1)}(g)\ge(1-r)\operatorname{dist}_a(g,V)^2\), with
\(V\) the two-dimensional translation space. On speeds
supported on \(\Gamma\) it is exactly the form of \(I-B\).
The distance from the unit sphere of \(V\) to that supported
subspace is positive, because a nonzero translation cannot
have zero normal speed on both exterior straight rays. It
follows that \(\|B\|<1\). The Neumann series for \((I-B)^{-1}\)
thus converges in weighted operator norm and preserves
nonnegativity.
[Weighted gap off translations](THREE_CELL_GLOBAL_INTERFACE_BOUNDS.md#a-global-matrix-bound-on-volume-flux)

By (TA6)–(TA7), \((I-B)(C_1w\pm f)\ge0\) on the tail.
Apply the positive inverse to get \(|f|\le C_1w\). The
translation speeds used here belong to the weighted space
because \(\int a\,d\nu<\infty\) by (G3). Normal continuity
extends the almost-everywhere inequality to the regular
graph. Since \(f^2+w^2=1\), it makes \(w\) uniformly positive
and bounds the slope by \(C_1\). The angular argument above
now forces \(\theta=2\pi/3\) and simplex rigidity.
