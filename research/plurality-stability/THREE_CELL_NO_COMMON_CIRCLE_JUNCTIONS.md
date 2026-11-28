# Junctions when no centered circle meets all three winners

Let \(A\) be an actual planar Gaussian noise-stability maximizer
at \(0<r<1\) for three prescribed positive masses. Write \(W_i\)
for the open strict-winner sets of the adjusted scores
\(V_i=T_r1_{A_i}-\beta_i\). Suppose no circle centered at the
Gaussian origin meets all three \(W_i\). A regular pair arc
at radii met by its spectator must lie in a centered circle.
There is no genuine triple junction at the origin.

Every noncollinear genuine triple junction lies at a positive
radius on a centered circular pair arc. For that pair, one of
the labels has a separate positive-mass group on the wrong side
of the circle, across a full annulus that it does not meet.
A collinear genuine triple has a different restriction: the
two walls of its thin middle-label region must bracket the
centered circle at the junction. In particular their radial
Taylor coefficients vanish through the order preceding the
cusp opening.

## Regular pair arcs in radii met by the spectator

Put
\[
 \mathcal R_i=\{s>0:W_i\cap\{|x|=s\}\ne\varnothing\}.
\]
These sets are open. The strict-winner sets coincide with
the ordinary cells almost everywhere, so almost every radial
slice has positive angular mass of label \(i\) exactly when
its radius is in \(\mathcal R_i\). A circle meeting all three
strict winners would imply the same property on an interval
of radii, by openness.

Let \(x\ne0\) be a regular point on an \(i/j\) jump arc with
\(|x|\in\mathcal R_k\), where \(k\) is different from \(i,j\).
If the arc is transverse to the centered circle through
\(x\), that circle meets both strict winners near \(x\).
This contradicts the hypothesis since it also meets \(W_k\).
The spectator-radius set is open, so on a sufficiently small
connected regular arc around \(x\), the derivative of radius
along the arc is identically zero. The arc lies in the
centered circle through \(x\), and real analyticity of the
angular restriction gives
\[
          V_i(s\omega)=V_j(s\omega)
               \quad\text{for every }\omega\in S^1
          \quad\text{at }s=|x|.                              \tag{RJ1}
\]
For each pair the radii with this exact full-circle equality
are locally finite, including no accumulation at zero.
Otherwise for every fixed \(\omega\) the one-variable
analytic contrast \(V_i(s\omega)-V_j(s\omega)\) would
vanish identically by accumulation at a finite radius.
This would make the planar contrast identically zero,
contrary to Gaussian convolution injectivity.

## Noncollinear junctions produce radial separation

The local jump classification gives a noncollinear genuine
triple three distinct tangent rays, and each of its strict
linear winner sectors has angle less than \(\pi\). If such
a junction is at zero, every sufficiently small centered
circle meets all three strict winners.
[Local analytic jump classification, G20](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#regular-interface-arcs-and-isolated-triple-junctions)

Let the junction \(p\ne0\) have radius \(R\). If neither
direction of the circle tangent at \(p\) is a wall direction,
the two opposite tangent directions enter distinct strict
winner sectors. The circle of radius \(R\) meets both of
those strict winners arbitrarily close to \(p\), by the
tangent sector description. Both labels occur on all
nearby circles, by openness of their radial images.
The third strict winner approaches \(p\), so it occurs on
some circle in this open radial interval, a contradiction.

Consequently a tangent half-ray is the direction of an
incident \(i/j\) wall. The opposite tangent direction lies
in the remaining strict sector \(k\), since both sectors
incident to the wall have angle less than \(\pi\).
Thus \(R\in\mathcal R_k\). The preceding regular-arc
argument makes the \(i/j\) wall circular on an open
segment ending at \(p\). After exchanging \(i,j\) if
necessary, \(i\) is on the inner side of that segment and
\(j\) is on the outer side. A small two-label neighborhood
of an ordinary point on the segment gives \(\varepsilon>0\)
such that
\[
 \begin{array}{ll}
 (R-\varepsilon,R)\subset\mathcal R_i\cap\mathcal R_k,
     & (R-\varepsilon,R)\cap\mathcal R_j=\varnothing,\\
 (R,R+\varepsilon)\subset\mathcal R_j\cap\mathcal R_k,
     & (R,R+\varepsilon)\cap\mathcal R_i=\varnothing .
 \end{array}                                                  \tag{RJ2}
\]
The exclusions apply to the whole annuli, because a strict
winner \(j\) in the first annulus or \(i\) in the second
would create a common circle.

For an actual maximizer the radial overlap of \(i,j\)
\[
 \delta_{ij}(A)=\int_0^\infty
                \min(p_i(s),p_j(s))\,s e^{-s^2/2}\,ds
\]
is positive, where \(p_i(s)\) is the normalized angular
mass of the \(i\)-slice. It is bounded below uniformly on
compact positive correlation bands and compact interior
mass sets.
[Radial overlap and its uniform lower bound](THREE_CELL_RADIAL_FIBER_EXCLUSIVITY.md#a-stability-bound-from-the-radial-overlap)
The portion of that integral at \(s<R\) is at most the mass
of label \(j\) at \(s\le R-\varepsilon\). The portion at \(s>R\)
is at most the mass of label \(i\) at \(s\ge R+\varepsilon\).
In particular one of those two wrong-side masses is at least
\(\delta_{ij}(A)/2\). Each is separated from its corresponding
component at the circular arc by a full annulus in (RJ2).
Thus at least one of \(i,j\) is disconnected. If the inner
wrong-side mass is positive, it is a bounded group separated
from all other mass of its own label by an annulus. In
particular, with two connected labels, a circular arc
supplied by a noncollinear triple in this situation cannot
be between those two labels.

## Collinear junctions bracket a centered circle

At a collinear genuine triple, the two extreme linear
winners have opposite tangent half-planes and the third
label occupies a thin one-sided or two-sided cusp.
Near the origin the extreme labels meet every small
positive circle. In the middle-label open cusp there is a
continuous path tending to the junction, for instance the
analytic midgraph of its two cusp walls on an active side.
If the vertex is the origin the radii of this path include
all sufficiently small positive values, giving a common
circle. This completes the origin exclusion in the
collinear case.

Now let a collinear junction be at \(p\ne0\) of radius \(R\).
If the circle tangent is transverse to the cusp direction,
its two tangent directions enter the two extreme strict
winner regions. Both extreme labels are then present on
an open radial interval about \(R\). The middle label
approaches \(p\), again giving a common circle. Thus the
cusp must be tangent to the centered circle.

On a small active side of the cusp, use local polar
coordinates \((\theta,s)\), with \(\theta=0\) at \(p\).
The two cusp walls are analytic radial graphs
\(s_-(\theta)<s_+(\theta)\), both equal to \(R\) at zero.
Here the strict middle winner lies between them, and
the two extremes are below and above, respectively.
Graph representation follows from the common tangent
being angular and the nonzero gradients of the pair
contrasts at a genuine triple.
The absence of a common circle forces
\[
                        s_-(\theta)\le R\le s_+(\theta)
                         \quad\text{on that side}.            \tag{RJ3}
\]
Indeed if \(s_-(\theta_0)>R\), choose \(u\) strictly between
\(R\) and \(s_-(\theta_0)\), sufficiently close to \(R\).
At \(\theta_0\) the circle \(s=u\) meets the lower extreme,
and near \(\theta=0\) it meets the upper extreme. At an
intermediate crossing of \(s_-(\theta)=u\), the upper wall
is strictly above \(u\). Moving slightly to the preceding
side of the first such crossing gives a point strictly
between the walls. Thus the middle label also occurs on
the same circle, a contradiction. If \(s_+(\theta_0)<R\),
the reversed argument gives the same contradiction.

Suppose the opening of the cusp has order \(q\ge2\),
so \(s_+(\theta)-s_-(\theta)=O(|\theta|^q)\).
This is the same order as the difference of the two
analytic wall heights in the local cusp coordinates
of (G20). An analytic change to polar coordinates is
invertible at \(p\), and both pairs of graphs have the
same tangent. Substituting and solving for the new
tangential coordinate preserves agreement through order
\(q-1\), and the nonzero normal derivative of the
coordinate change preserves the first nonzero normal
difference at order \(q\). By (RJ3),
\[
              |s_\pm(\theta)-R|\le s_+(\theta)-s_-(\theta)
                               =O(|\theta|^q).
\]
Thus at a cusp of order at least three the two walls
have the curvature of the corresponding centered circle
at the junction. This restriction and (RJ1)–(RJ2) follow
from the absence of a common circle. They do not exclude
all partitions of that kind.
