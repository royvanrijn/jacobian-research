# Reflection from a straight isolated end

Let \(A\) be a balanced planar three-cell Gaussian noise-stability
maximizer at \(0<r<1\). Suppose there is a bounded isolated graph
end for each pair. Thus, in an orthonormal coordinate system
\((t,z)\) about the Gaussian origin for that pair, its entire
interface in a full forward half-plane \(t>S\) is one bounded
continuous graph \(z=h(t)\). The two incident labels fill
opposite sides of the graph, up to null sets.

The [shifted-contour theorem](THREE_CELL_TAIL_FLATNESS.md)
now proves straightness of every such bounded isolated end.
The compact-source uniqueness statement below also applies
to curved Gaussian-flat graphs without a full isolated
source half-plane.

If the three outward directions contain no antipodal pair and
one end is straight, its supporting line passes through the
origin. Reflection in this line exchanges the incident labels
and preserves the third label outside a compact set. The other
two ends are reflections of one another. If they are curved,
the reflection holds on the whole plane. If either is straight,
both are straight and the exterior is a centered regular simplex.
In the antipodal case neither member of the antipodal pair
can be straight. In particular, at every positive correlation,
two straight ends force three straight ends and a centered
regular-simplex exterior.

The passage from exterior to global reflection for curved ends
uses an analytic uniqueness statement for compactly supported
contrast changes. A Gaussian convolution of a compactly
supported function cannot vanish on a curved graph which
approaches an affine line at a Gaussian rate.

## Slices through a straight end

Write \(\sigma=\sqrt{1-r^2}\), let \(p_\sigma\) denote the
centered Gaussian density of standard deviation \(\sigma\),
and let \(\Phi\) be the standard normal cumulative distribution
function. For a pair contrast \(f\), its smoothed score
\(T_r f\) is real analytic and has constant value \(\beta\)
on the interface. If the end is straight at offset \(b\),
analytic continuation gives \(T_r f(t,b)=\beta\) for all \(t\).
Since the longitudinal convolution has a nowhere-zero Fourier
multiplier, this implies
\[
     \int f(s,v)p_\sigma(v-rb)\,dv=\beta
                       \quad\text{for almost every }s.    \tag{ER1}
\]
The convolution is with variance \(\sigma^2\) and evaluated
at \(rt\). Its injectivity on bounded functions follows from
the Fourier multiplier identity on tempered distributions,
testing against compactly supported smooth functions in frequency.
[Pressure offsets and isolated-end regularity](THREE_CELL_ISOLATED_END_RIGIDITY.md#conditional-equation)

The three outward directions have pairwise nonpositive dot
products. A positive dot product would put one end eventually
inside another's forward half-plane, where the label belonging
only to the first pair is excluded. Boundedness of the
transverse heights makes this impossible on a positive-area
subset of the two-label side. Suppose no two directions are
antipodal. Select coordinates \((s,v)\) for the straight
end and call its upper and lower labels \(i,j\), with
spectator label \(k\). The other directions have the form
\[
   (-c_+,d_+),\quad (-c_-,-d_-),\qquad
             c_\pm\ge0,\quad d_\pm>0 .
\]
Indeed they lie in the closed left semicircle. If they lay on
the same side of the horizontal axis, their mutual angular
separation would be less than \(\pi/2\), except when one was
antipodal to the straight direction. Since they are also
non-antipodal to each other, at least one \(c_\pm>0\).
By overlap with the right forward half-plane at very large
positive or negative \(v\), the upper arm has outer label
\(i\), and the lower arm has outer label \(j\).
Their shared label is \(k\).

A small cone about the reverse of the straight end belongs
to label \(k\). If both \(c_\pm>0\), this cone lies in the
overlap of the other arms' forward half-planes. If one
\(c_\pm=0\), the negative horizontal direction is strictly
inside the other arm's half-plane and has nonzero transverse
projection there. Boundedness of that arm's height gives a
pure-label cone about the direction. Tilt into the remaining
forward half-plane to see that this pure label is \(k\).
Consequently along the straight supporting line the pair
score \(T_r(1_{A_i}-1_{A_j})(t,b)\) tends to zero as
\(t\to-\infty\). This uses the bounded transverse spread of
the conditional Gaussian and its mean \(r(t,b)\) in the
interior of the reverse cone for large \(-t\). The pressure
is therefore zero. The flat pressure equation
\(\beta=2\Phi((r-1)b/\sigma)-1\) gives \(b=0\).

For sufficiently negative \(s\), the other two forward
half-planes cover the entire transverse source slice.
In fact their intersections with that slice are of the
form \(v>(S_++c_+s)/d_+\) and \(v<-(S_-+c_-s)/d_-\),
which overlap for large \(-s\). Neither \(c_\pm\) can
vanish. For example, if \(c_+=0\), the upper arm's
supporting longitudinal direction is vertical and its
transverse coordinate is bounded. The upper forward
half-plane on this slice then contains only \(k\) for
large negative \(s\), while the lower forward half-plane
excludes \(i\) and contains \(j\) on a transverse lower
tail. This contradicts (ER1) with \(b=\beta=0\).
The case \(c_-=0\) is symmetric.

When \(c_\pm>0\), the \(s\) coordinate along either of the
other ends has derivative \(-c_\pm+o(1)\), by the convergence
of the transverse derivatives to zero. It tends to
\(-\infty\), and thus each arm has exactly one intersection
with each sufficiently negative \(s\) slice. Increase the
forward cutoffs before taking these slices, so that the
derivative is strictly negative throughout each parameter
interval. The intersection heights \(v_+(s)\) and \(v_-(s)\)
grow respectively to \(+\infty\) and \(-\infty\), with
\(v_\pm(s)=\mp (d_\pm/c_\pm)s+O(1)\).
The two interfaces separate the pure upper and lower tail
labels from \(k\) in the middle. To see this, each arm's
line of intersection lies outside the other's open forward
half-plane, since that other half-plane excludes its outer
label. The two half-planes overlap between the intersections,
and their overlap has the common label \(k\).
Applying (ER1) gives
\[
     1-\Phi(v_+(s)/\sigma)=\Phi(v_-(s)/\sigma),
     \qquad v_-(s)=-v_+(s) .                              \tag{ER2}
\]
The equality extends from almost every slice by continuity
of the two intersection heights. This proves exact
reflection on the far left. On the right, the straight
graph \(v=0\) gives exact reflection. On the remaining
bounded interval in \(s\), sufficiently large \(v\) lies
in the upper arm's forward half-plane on the outer-label
side of its bounded graph, uniformly for all such \(s\).
Sufficiently negative \(v\) similarly lies in label \(j\).
Thus the reflection fails, if at all, only in a compact set.

If one of the two other arms is straight, (ER2) makes both
straight. In this case every pressure vanishes by the
reverse-cone argument, so the three supporting lines are
centered. Repeating (ER2) with each of the three ends as
the reference makes the two neighboring angular gaps equal
at every ray. Hence all angles are \(2\pi/3\). Their forward
half-planes cover the complement of a large disk and give
exactly the labeled simplex there.
[Related exterior geometry and the compact-core moment conditions](THREE_CELL_PRESSURE_GRADIENT_TAIL.md#exterior-geometry-when-all-three-ends-straighten)

## A compact contrast cannot vanish on a curved flat end

Let \(q\) be bounded and compactly supported in the plane.
Suppose \(T_r q(t,h(t))=0\) for every sufficiently large
\(t\), where \(h\) is continuous, is analytic on a tail,
is not eventually constant, and satisfies
\(|h(t)-b|\le C e^{-c t^2}\) there for some \(c>0\).
Then \(q=0\) almost everywhere.

With \(\alpha=r/\sigma^2\), remove the nonvanishing Gaussian
factor in the observation from the convolution. The remaining
function, up to a constant prefactor, is
\[
 H(t,z)=\int e^{\alpha(ts+zv)}
               q(s,v)e^{-(s^2+v^2)/(2\sigma^2)}\,ds\,dv.
                                                               \tag{ER3}
\]
It is entire in both variables. For \(t\ge0\) and
\(|z-b|\le1\), each fixed \(z\) derivative is bounded by
\(C_m e^{Lt}\), with \(L\) fixed by the support.

For a finite signed or complex measure \(\mu\) supported
on \([-B,B]\), its real Laplace transform
\(L_\mu(t)=\int e^{ts}\,d\mu(s)\) cannot decay faster than
every exponential on \(t\to+\infty\) unless \(\mu=0\).
In the faster-decay case, the function
\[
    z\longmapsto \int_0^\infty e^{-zt}L_\mu(t)\,dt
\]
is entire by compact-uniform domination. In \(\Re z>B\),
Fubini identifies it with the Cauchy transform
\(\int (z-s)^{-1}\,d\mu(s)\). They agree on the
connected complement of \([-B,B]\) by analyticity.
The Cauchy transform tends uniformly to zero at infinity
there, so the entire function vanishes by Liouville's
theorem. Its power expansion at infinity shows that every
moment of \(\mu\) vanishes. Polynomial density on the
support interval gives \(\mu=0\).

Now \(H(t,h(t))=0\) implies
\(|H(t,b)|\le C_1e^{Lt}|h(t)-b|\), by the derivative bound.
The function \(H(t,b)\) is a one-dimensional Laplace
transform of a compactly supported measure in \(\alpha s\),
so it is identically zero. Inductively suppose the
\(z\)-derivatives of orders below \(m\) vanish at \(b\)
for every \(t\). Taylor's formula at each \(t\) where
\(h(t)\ne b\) gives
\[
      |\partial_z^m H(t,b)|\le C'_m e^{Lt}|h(t)-b|.
\]
The zeros of \(h-b\) on the analytic curved tail are
isolated, so the bound extends to them by continuity.
This derivative is again a compact-measure Laplace
transform, and hence vanishes identically. The induction
makes every transverse jet zero on the supporting line.
Entire analyticity in \(z\) now gives \(H=0\) on the plane,
and injectivity of the Gaussian convolution proves \(q=0\).

More generally, if two stationary partitions at the same
correlation agree outside a compact set and share such an
end, the pair contrasts agree everywhere. Indeed, the
difference of the pair pressures on the shared end is the
smoothed difference of two source contrasts. That source
difference has compact support, so its smoothed value
tends to zero there and the pressures agree.
[Gaussian flattening for a bounded isolated end](THREE_CELL_TAIL_FLATNESS.md)

Apply this to \(A\) and its reflection with \(i,j\) exchanged.
The reflected partition is also a maximizer, and (ER2)
and the slab argument give agreement outside a compact set.
If the two other arms are curved, apply the compact
contrast result on each. The differences of the two
partitions' \(i,k\) and \(j,k\) indicator contrasts then
vanish everywhere. Since the three indicator discrepancies
sum to zero, each vanishes. This proves global reflection.
In this case the known single-reflection theorem gives
the global monotone graph and connected strict-winner cells.
[Single-reflection reduction](THREE_CELL_REFLECTION_REDUCTION.md)

## The antipodal case

If two directions are antipodal, the third is perpendicular
to both. Rotate and relabel so that the \(1/2\) direction
is rightward, \(2/3\) is upward, and \(1/3\) is downward.
In the upper half-plane label 2 is to the right of its
bounded graph and label 3 to its left. In the lower half-plane
label 1 is to the right and label 3 to its left. The
orientations follow by overlap on the far right with the
\(1/2\) half-plane, where only the shared label is possible.

Suppose the upper arm is straight at horizontal offset \(p\)
and use contrast \(1_{A_2}-1_{A_3}\). On every sufficiently
late lower source slice its transverse expectation with
mean \(rp\) is \(-\Phi((h_{\rm low}-rp)/\sigma)\), because
label 2 is absent. Formula (ER1) makes this independent of
the lower longitudinal coordinate. Strict monotonicity of
\(\Phi\) makes the lower arm straight at an offset \(q\).
The upper pressure equation gives
\[
  2\Phi((r-1)p/\sigma)-1=-\Phi((q-rp)/\sigma).
\]
Its right side is negative, so \(p>0\). Writing
\(u=(1-r)p/\sigma>0\), the identity gives
\(\Phi((q-rp)/\sigma)=2\Phi(u)-1<\Phi(u)\), hence \(q<p\).
The same argument with the lower arm gives \(p<q\).
Thus neither antipodal member can be straight.

With three bounded isolated ends, the shifted-contour
theorem applies to each and precludes an antipodal pair.
The compact core of the resulting simplex exterior is a
separate rigidity problem.

## Two isolated ends and a second-order tie

The existence of a graph end for the third pair is unnecessary
for two pressure conclusions. Suppose two distinct pairs have
bounded isolated graph ends. Their outward directions have
nonpositive dot product by the absent-label argument. Both
ends are straight by shifted-contour straightening.

They cannot be antipodal. Otherwise rotate the two ends to
right and left half-planes and write their height offsets as
\(p,q\), respectively. Write \(a\) for the common label, and
choose the height orientation so that \(a\) is above the right
line. The pair contrast oriented positively for \(a\) is
constant along each entire supporting line. At its opposite
end the other, noncommon label in that contrast is absent.
Thus each pressure there equals the positive Gaussian
probability of \(a\) in the opposite end.

For example, on the right supporting line its own end gives
the pressure \(2\Phi((r-1)p/\sigma)-1>0\), so \(p<0\).
If \(a\) is above the left line too, its probability on the
reverse tail of the right line is \(\Phi((rp-q)/\sigma)\).
Since \(2\Phi(u)-1<\Phi(u)\) for finite \(u\), matching
pressures yields \(q>p\). Reversing the two ends gives
\(p>q\). If \(a\) is below the left line, its pressure on
that line is \(2\Phi((1-r)q/\sigma)-1>0\), so \(q>0\).
Matching the opposite-tail probabilities on the two lines
now gives \(q<(1-2r)|p|\) and \(|p|<(1-2r)q\), which
are incompatible for \(0<r<1\).

If the dot product of the two outward directions is strictly
negative, both lines pass through the origin. Indeed the
reverse tail on the first line lies increasingly deep in
the forward half-plane of the second. The directions are
neither parallel nor antiparallel, so it is also increasingly
far from the second line. Hence on conditioning at that
tail, the first contrast, oriented positively for \(a\),
has limit either \(0\) or \(1\). Since a pair pressure has
absolute value strictly below one, it must be zero. The
binary pressure formula then puts the first line through
the origin. Interchanging the two proves the same for the
second. The two zero pair pressures make all three prices
equal.

If the strict-winner set of the common label is connected,
the angle between these outward directions is at least
\(2\pi/3\). In fact the overlap of the two positive
longitudinal half-planes contains a sector at infinity
between the directions and has the common label. The
strict-winner set avoids both supporting lines, so by
connectedness it lies in the sector of angle \(\theta\)
between the two outward rays. Its Gaussian mass \(1/3\)
is at most \(\theta/(2\pi)\). At equality the common
cell fills that sector and single-simplex-cell rigidity
identifies the entire partition. Thus a nonsimplex with
the common winner connected needs \(\theta>2\pi/3\).
[Single-simplex-cell rigidity](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#one-simplex-cell-determines-the-whole-partition)

Suppose additionally that \(A\) ties the centered simplex
\(T\) at \(r\) and that simplex-optimal correlations
accumulate at \(r\). Put \(D(s)=S_s(T)-S_s(A)\). With two
strictly obtuse bounded isolated ends, the conditions
\(D'(r)=D''(r)=0\) force \(A=T\) up to rotation and
relabeling. In fact (G10) makes the fixed-partition
correlation-derivative functional annihilate all compact
balanced interface speeds. The flux map is onto the
centered label space by complete adjacency and regular
spanning-tree bumps. The derivative on each pair is
therefore a difference of constant label numbers, even
if the interface has several components. On each of the
two straight rays through the origin this derivative
tends to zero. The centered binary halfspace contrast
vanishes at every correlation on its supporting line,
and differentiation of the source error outside the
isolated half-plane preserves Gaussian tail decay.
Two zero differences force the third difference to vanish.
The dilation inequality (G6) then forces zero dilation
speed, hence centered conicality and the regular simplex.

This excludes higher-order contacts with these two ends,
including when the third pair has no controlled graph.
[Mixed inequality (G10)](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#a-first-crossing-obstruction-at-the-level-of-hermite-weights)
[Dilation inequality (G6)](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#translation-equations-on-an-unconstrained-space-of-speeds)
[Conical rigidity](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#convex-and-conical-maximizers)
