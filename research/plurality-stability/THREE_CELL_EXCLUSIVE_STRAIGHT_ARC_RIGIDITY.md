# Straight or algebraic arcs in an exclusive-fiber maximizer

Let \(A\) be an actual balanced planar three-cell Gaussian maximizer
at \(0<r<1\). Suppose labels 1 and 2 have exclusive vertical fibers
and that the resulting third-cell epigraph has a finite continuous
height \(h\). Use the
[exclusive-fiber representation and height laws](THREE_CELL_FIBER_EXCLUSIVITY_RIGIDITY.md#sorting-fibers-and-identifying-their-base-labels).
In particular there is a locally finite nonempty set \(B\) of
base switches, and the two signed height laws agree under
\(p_\sigma(s-rb)\,ds\) at each \(b\in B\), where
\(\sigma=\sqrt{1-r^2}\).

Each of the following conditions forces the centered regular simplex.

- Some open \(3/i\) interface arc, \(i\in\{1,2\}\), is straight.
- One of the other two label pairs also has exclusive fibers in
  some direction.
- Some open \(3/i\) interface arc is contained in a real algebraic
  curve.
- On a full connected component of the existence domain of one
  pair-score root \(H_i\), that root is definable with parameters
  in a polynomially bounded o-minimal expansion of the real field.

For \(r\ge 2/3\), the height law also gives a convexity criterion
without assuming exclusive fibers in advance. If two cells have
convex representatives, the maximizer is the centered simplex.
More generally, an affine-separated pair and any one convex cell
suffice.

## A straight arc propagates across the base label

Suppose an open \(3/i\) arc lies on a line. The graph over the
horizontal coordinate makes this a nonvertical line \(y=L(x)\),
with \(L\) affine. Its adjusted pair contrast \(V_3-V_i\) vanishes
along this whole line by real analyticity. The contrast is
strictly increasing in the vertical coordinate. In fact its
unsmoothed contrast is vertically nondecreasing and varies on
a positive-measure family of fibers, and Gaussian convolution
gives the derivative its strict sign. Therefore on every open
base interval of label \(i\), the unique graph root satisfies
\(h(x)=L(x)\).

The affine function \(L\) is nonconstant. Otherwise every
height in label \(i\) would have the same value. At any base
switch the signed height law (FX3), tested against all other
heights, would force the other label's heights to have that
value too. Continuity would then make \(A_3\) a half-plane,
contrary to the
[individual-cell dimension restriction](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#each-cell-uses-both-active-directions).

If the other base label has a bounded interval \((a,b)\), continuity
from its neighboring intervals of label \(i\) gives \(h(a)=L(a)\)
and \(h(b)=L(b)\). Let \(J\) be the open interval between these two
distinct heights. A positive-length subset of \((a,b)\) has heights
in \(J\), by continuity. No interval of label \(i\) can have a point
with height in \(J\), since on that label \(h=L\) and the preimage
\(L^{-1}(J)\) is precisely \((a,b)\). This contradicts the signed
height law (FX3).

Thus the other base label has at most two components, each
unbounded, and there are at most two switches. The
[finite-switch rigidity theorem](THREE_CELL_PROPER_SWITCH_HEIGHT_RIGIDITY.md#critical-height-accumulation-excludes-multiple-finite-switches)
makes \(B\) a singleton. At its point \(b\), the height law for
the entire height range gives equal masses of the two sides
of \(b\) under \(N(rb,\sigma^2)\). This forces \(b=0\), for
every \(0<r<1\).

Orient the horizontal coordinate so that label \(i\) is on
the right. Its half-profile \(h(s)=L(s)\) is monotone on
\(s>0\). The [common height law and likelihood-ratio
comparison](THREE_CELL_AFFINE_GRAPH_REFLECTION_TESTS.md#the-common-height-distribution)
give \(h(-s)=h(s)\). The resulting height is affine in
\(|x|\), and the two graph rays and the downward base
ray form a partition conical about \((0,h(0))\). The
[translated-cone rigidity theorem](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#cones-centered-at-an-arbitrary-point)
identifies the maximizer as the centered regular simplex.

## A second exclusive direction

Suppose labels \(i,3\) also have exclusive fibers, with \(i=1\)
or \(2\). Use the open representatives of (FX1), ignoring the
discrete switch lines for the base labels. Projections of open
sets are open, and any nonempty fiber of an open set contains
an open interval in that fiber. Thus overlap of the two
projections would give a positive-length set of projection
coordinates where both slices have positive measure.
Fubini transfers this statement to the original indicators
because the representatives agree almost everywhere. Hence
the projections of the open representatives are disjoint.

Let \(n=(n_x,n_y)\) be a unit normal to these second fibers.
The case \(n_y=0\) gives the original vertical fibers, which
contain both labels at each interior abscissa of a base-\(i\)
interval. After orienting \(n\), we may take \(n_y>0\).
The projection of \(A_3\) is then the upper interval with
endpoint
\[
 L_*=\inf_{x\in\mathbb R} [n_xx+n_y h(x)],
\]
and that of \(A_i\) is the lower interval with endpoint
\[
 U_*=\sup_{x\ \mathrm{in\ base}\ i}[n_xx+n_y h(x)].
\]
Disjointness and nonemptiness of these intervals give
finite \(U_*\le L_*\). On every base-\(i\) interval the
expression inside the supremum is also at least \(L_*\),
by the definition of the infimum. Therefore
\(U_*=L_*\) and \(n_xx+n_y h(x)=L_*\) at each such
abscissa. This supplies an open straight \(3/i\) arc,
so the preceding argument applies.

## One algebraic arc extends to a global pair-score root

Suppose a nonzero real polynomial vanishes on an open
\(3/i\) graph arc. Put \(F_i=V_3-V_i\) and
\(\pi_i=T_r^{[1]}1_{\{\ell=i\}}\). The common price is
\(\beta=\beta_3-\beta_1=\beta_3-\beta_2\) by the signed
height law. At horizontal observation \(x\), the lower
vertical limit of \(F_i\) is \(-\pi_i(x)-\beta\), and
the upper limit is \(1-\beta>0\). These follow by
dominated convergence since all heights are finite.
The strict inequality \(\beta<1\) holds at the active
pair interface. Therefore \(F_i\) has a unique
real-analytic zero \(H_i(x)\) on the open set
\(J_i=\{x:\pi_i(x)+\beta>0\}\). This zero equals \(h\) on
the base-\(i\) intervals.

On the connected component of \(J_i\) containing the
algebraic arc, choose a real irreducible factor \(P\)
of its carrier polynomial such that \(P(x,H_i(x))=0\).
One factor must vanish identically along the arc,
since a finite union of discrete analytic zero sets
cannot cover an interval. The relation then extends
along the component by real analyticity.

Suppose this component has a finite endpoint \(a\).
Here \(\pi_i(a)+\beta=0\). The function \(\pi_i+\beta\)
is real analytic and does not vanish identically,
by injectivity of Gaussian smoothing and the positive
masses of both base labels. Hence it is bounded below
along the component near \(a\) by \(c|x-a|^n\) for
some fixed integer \(n\ge1\) and \(c>0\). At \(x=a\)
the contrast \(F_i(a,y)\) is strictly positive for
every finite \(y\), since its lower limit is zero
and it is strictly increasing. This implies
\(H_i(x)\to-\infty\) as \(x\to a\) within the component.

At the root, the exact conditional formula is
\[
 \pi_i(x)+\beta
   =\mathbb E\big[(1+1_{\{\ell(X)=i\}})
                         1_{\{Y>h(X)\}}\big],
 \qquad (X,Y)\sim N\big((rx,rH_i(x)),\sigma^2 I\big).
\]
By the [linear height bound](THREE_CELL_FIBER_EXCLUSIVITY_RIGIDITY.md#a-growth-bound-even-with-infinitely-many-switches),
\(h(X)\ge-C(1+|X|)\). Write \(H=|H_i(x)|\). For
large \(H\), the event in the expectation lies in
\(\{|X|>rH/(4C)\}\cup\{Y>-rH/2\}\). Since \(x\) is
bounded near \(a\), both Gaussian probabilities are
\(O(e^{-c_1H^2})\) for a fixed \(c_1>0\).
The power lower bound for \(\pi_i+\beta\) consequently
gives \(H_i(x)=O(\sqrt{\log(1/|x-a|)})\) in absolute
value.

Such slow divergence is impossible on \(P=0\).
Let \(d\) be the degree in the height variable and put
\(Q(x,z)=z^dP(x,1/z)\). The polynomial \(Q(a,z)\)
is not identically zero, since otherwise \(x-a\)
would divide the irreducible polynomial \(P\), which
depends on height. At \(z=1/H_i(x)\to0\), the identity
\(Q(x,z)=0\) gives \(Q(a,0)=0\) and
\(|Q(a,z)|\le C_2|x-a|\). If \(k\ge1\) is the
vanishing order of \(Q(a,z)\) at zero, this forces
\(|H_i(x)|\ge c_2|x-a|^{-1/k}\), a contradiction.
The component thus has no finite endpoint and is all
of \(\mathbb R\).

## The algebraic height determines both tails

The global algebraic graph is nonconstant,
as a constant \(H_i\) would force all heights constant
by the signed height law. The polynomial \(P\) therefore
depends on both coordinates. In the height variable its
leading coefficient and discriminant are nonzero polynomials
in \(x\), by irreducibility over \(\mathbb R[x,y]\) and
separability over \(\mathbb R(x)\) in characteristic zero.
Off their finitely many zeros, the real roots are distinct
and can be ordered. Continuity makes \(H_i\) choose one
fixed root index on each of the resulting open intervals.
The condition that a root is of a prescribed index is
semialgebraic by real quantifier elimination, by counting
the roots strictly below it. Adding the finitely many
endpoint graph points makes \(H_i\) semialgebraic on
the whole line.
[Coste, Tarski–Seidenberg, Theorem 2.6](https://arxiv.org/pdf/2611.31550)

The semialgebraic monotonicity theorem gives extended-real
limits of \(H_i\) in both directions. Its derivative is
nonzero eventually on each nonconstant tail. Indeed a
critical point where the height root is simple would be
a common zero of \(P\) and \(P_x\). They are coprime
since \(P_x\ne0\) and has lower degree in \(x\), and
there are only finitely many common plane points.
[Coste, Monotonicity Theorem 2.1](https://arxiv.org/pdf/2611.32426)

Choose a bounded open height band with positive
base-\(i\) preimage measure and with closure avoiding the
two extended-real limits. This is possible on a
nonconstant base-\(i\) interval. Its base-\(i\)
preimage is bounded, since there \(h=H_i\).
If \(B\) were infinite, it would escape in at least
one direction, where every height band of positive
preimage measure must recur in both base labels by
(FX3). Consequently \(B\) is finite and thus consists
of one point by the finite-switch theorem. As above,
that point is zero.
[Height-band recurrence in both base labels](THREE_CELL_FIBER_EXCLUSIVITY_RIGIDITY.md#counting-switches-with-a-height-band)

The half-profile on the side of label \(i\) is the
restriction of \(H_i\), oriented in the outward
coordinate \(s>0\). This gives the hypothesis of the
following height-law criterion.

## One regular monotone tail forces reflection

Consider a centered affine-separated maximizer with finite
continuous height, and suppose one outward half-profile
has nonzero derivative on a tail. Neither half-profile is
constant. Constancy of one would propagate along its whole
analytic branch, and the common height law would make the
other constant as well. A monotone tail cannot have a finite
limit. Such a limit makes it a bounded isolated end approaching
its pressure offset from one side, so exact one-sided straightening
would make its profile constant. The monotone tail therefore
escapes to \(+\infty\) or \(-\infty\).
[Bounded isolated ends and one-sided straightening](THREE_CELL_ISOLATED_END_RIGIDITY.md)

The opposite half-profile has the same height law
under \(p_\sigma(s)\,ds\). For an upward escape,
above a sufficiently large height the known
profile has only its single monotone
tail branch and its height law has a locally bounded
density. A critical point of the opposite profile
at any such height would give nonnegative local mass
of order \(\varepsilon^{1/m}\), \(m\ge2\), in a
symmetric band of width \(2\varepsilon\). This is
incompatible with the \(O(\varepsilon)\) upper
bound from the density of the known profile.
There are therefore no critical points of the
opposite profile at those high levels. Its
high superlevel set is nonempty and has an
unbounded component, since any bounded component
away from zero would have an interior maximum,
and the levels can be taken above the value at
zero. It follows that its tail is also eventually
monotone with upward escape. For a downward
escape, the same argument uses sublevel sets and
critical minima.

On the corresponding extreme height range,
each inverse tail threshold is determined by the
same \(p_\sigma\)-height law, as the preimages
of the extreme superlevel or sublevel sets are
upper half-lines in \(s\). Thus the profiles
agree on their tails. Analytic continuation
gives exact reflection on the whole half-line. Thus an asymmetric
maximizer in this class has critical points arbitrarily far out
on both half-profiles.

In the algebraic case the full graph is now semialgebraic.
The strict winner regions are connected. In particular the
equal scores \(U_1(0,y)=U_2(0,y)\), followed by
horizontal Gaussian tilting and the opposite-side
supports, give \(V_1>V_2\) on the right and the
reverse on the left. The other two comparisons
are the strictly increasing contrasts against
the finite continuous height. The
[connected semialgebraic rigidity theorem](THREE_CELL_CIRCULAR_INTERFACE_OBSTRUCTION.md#connected-semialgebraic-maximizers)
therefore gives the centered regular simplex.

## Polynomially bounded definable roots

For a pair root \(H_i\) definable on a full component of
\(J_i\), the same finite-endpoint estimate gives
\(|H_i(x)|=O(\sqrt{\log(1/|x-a|)})\) at a hypothetical finite
endpoint \(a\), where \(H_i\) diverges to \(-\infty\). By the
definable monotonicity theorem this branch is strictly
monotone near \(a\), and its inverse \(x(-t)\) is definable
for large \(t>0\). Polynomial boundedness applied to
\(1/|x(-t)-a|\) would force \(|H_i(x)|\ge c|x-a|^{-1/N}\)
for some \(c,N>0\). This contradicts the logarithmic estimate.
Thus \(H_i\) is definable on \(\mathbb R\).
[Coste, Theorem 2.1](https://arxiv.org/pdf/2611.32426)

The global root has extended-real limits at both ends by
definable monotonicity. Its derivative is definable by the
first-order formula for a derivative, and its zero set contains
no interval by analyticity and nonconstancy. O-minimality makes
that set finite. As for the algebraic root, a height band away
from the two limits has bounded preimage on base label \(i\)
and positive preimage measure. This yields one centered switch
by (FX3) and finite-switch rigidity, and the regular monotone
tail criterion gives reflection. Hence \(h\) itself is definable
in the same structure.

Now \(h''\) is definable, so it has a fixed weak sign on a
tail. An eventually convex reflected maximizer is the simplex.
In the eventually concave case the limiting slope \(m\) exists
by the global slope bounds for reflected maximizers, and \(m\le1\)
also gives the simplex. For \(m>1\), the nontrivial concave case
has a strictly positive affine height defect \(q(t)\) with
\(\lim_{t\to\infty}t^{-2}\log q(t)<0\) by (CS1). The defect is
definable, since it is obtained from the graph by an orthogonal
linear change of coordinates and subtraction of a constant.
Polynomial boundedness of \(1/q\) contradicts this Gaussian
decay. This proves the definable-root criterion.
[Reflected eventual convexity and small-slope concavity](THREE_CELL_REFLECTED_CONSTANT_CORRELATION_RIGIDITY.md#a-convex-tail-forces-simplex-rigidity).
[Concave Gaussian defect, CS1](THREE_CELL_CONCAVE_STEEP_TAIL_ASYMPTOTICS.md#the-height-rate-and-the-source-locations).

## Convex cells

Suppose \(2/3\le r<1\). This covers every possible first
transition correlation, since the explicit global comparison
reaches \(339/500\). In this range affine-separation rigidity gives
a centered separating line and a finite continuous graph. If the
third cell is convex, its epigraph height is convex and the
convex-height classification gives the simplex. If one of the
two base cells is convex, its
half-profile is concave, since its hypograph on that half-plane is
convex.

These convexity conclusions apply to representatives agreeing
almost everywhere. For a positive-measure convex representative,
its interior is canonical up to closure, and the closure of its
interior is its essential closure. The epigraph or hypograph
of a finite continuous height, on its open horizontal domain,
is the interior of its essential closure on that domain. The
convex representative therefore makes the open epigraph or
hypograph convex there.

A nonconstant analytic concave half-profile has nonzero
derivative eventually, since the derivative is nonincreasing
and cannot be zero on an interval. The regular monotone tail
criterion makes both half-profiles equal, so the reflected
profile is concave on its entire half-line. Global concavity
rigidity gives the simplex. Finally, any two convex cells have
disjoint nonempty interiors and so admit a separating affine
line. This proves the two-convex-cell criterion.
[Explicit global interval](THREE_CELL_CENTRAL_SHIFT_BOOTSTRAP.md).
[Affine separation](THREE_CELL_AFFINE_SEPARATION_RIGIDITY.md).
[Convex exclusive-fiber height, FX5](THREE_CELL_FIBER_EXCLUSIVITY_RIGIDITY.md#classification-for-a-convex-height).
[Global reflected concavity](THREE_CELL_REFLECTED_SHAPE_RIGIDITY.md#global-convexity-or-concavity-forces-the-simplex).
