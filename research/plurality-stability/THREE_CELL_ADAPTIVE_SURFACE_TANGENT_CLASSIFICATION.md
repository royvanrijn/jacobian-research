# The possible tangent cones of adaptive surface minimizers

Fix \(0<r<1\) and let \(A\) be a balanced planar minimizer
of the adaptive surface functional \(\mathcal F_r\).
Every nonconstant tangent partition is either a line
separating two labels or a tripod with all three labels
and positive tensions in force balance. At a tripod
the scaled tagged interface measures converge to the
actual tripod measures. At a line any excess interface
is precisely the measure of a collapsed layer of the
absent label. Such excess forces equality of the
triangle inequality at every attaining test that is
sharp on replacement by the line.

The local upper perimeter bound and the stationarity
equation ensure homogeneous tangent partitions. The
step here is to compare patches while retaining the
supremum over attaining tests. A patch selects its own
test, so it is important to account for all such tests
and for any interface carried by the limiting measure
that the limiting cells do not show.
[Local perimeter bounds and homogeneous tangents](THREE_CELL_ADAPTIVE_SURFACE_TANGENTS.md)
[Attaining tests, reciprocal weights and pressure stationarity](THREE_CELL_ADAPTIVE_SURFACE_STATIONARITY.md)

## The comparison at a tangent scale

Write \(E_n=(A-p)/t_n\) for a sequence of rescaled
partitions at a fixed observation \(p\), with \(t_n\downarrow0\).
Use Euclidean pair interface measures \(\sigma_{n,e}\)
for \(E_n\), counting each pair once. After passing to a
subsequence the cells converge in \(L^1_{\mathrm{loc}}\)
to a homogeneous partition \(C\) and the measures converge
locally weakly to \(\eta_e\). These conclusions use the
local perimeter bound at \(p\) and the averaged
stationarity equation. A homogeneous planar partition
of locally finite perimeter has finitely many rays:
on a polar annulus almost every angular slice has finite
one-dimensional perimeter, and homogeneity makes the
angular partition independent of radius.

Let \(\mathcal M(A)\) be the compact set of attaining
normalized Gaussian Hilbert tests of \(A\). For
\(G\in\mathcal M(A)\) use the frozen tensions
\[
 w_e(G)=|T_{\sqrt r}G_i(p)-T_{\sqrt r}G_j(p)|,
 \qquad e=(i,j).
\]
If \(w_e(G)=0\) for even one attaining test \(G\), then
\(\eta_e=0\). Indeed its weight \(w_e^G(x)\) is locally
Lipschitz, so, for \(t\downarrow0\),
\[
 \nu_e(B(p,t))\le
     C t\int_{B(p,t)} \frac{1}{w_e^G}\,d\nu_e=o(t)
\]
by the reciprocal estimate (ASV4). Here \(\nu_e\) is
the Gaussian tagged measure of \(A\). The Euclidean
and Gaussian measures differ by a positive smooth
factor near \(p\), giving the claim at every fixed
rescaled radius.

Fix a radius \(R>0\) where every \(\eta_e\) gives zero
mass to \(\partial B_R\), and let \(D\) be an ordinary
finite-perimeter competitor for \(C\) differing from it
in a compact subset of \(B_R\). Denote the pair perimeter
of \(D\) in the open ball by \(P_e(D;B_R)\). Then
\[
 \max_{G\in\mathcal M(A)}
     \sum_e w_e(G)\bigl[P_e(D;B_R)-\eta_e(B_R)\bigr]
       \ \ge 0.                                             \tag{TC1}
\]
More precisely, the left side multiplied by the
Gaussian density \(\phi_2(p)\) is the first-order
increment in \(\mathcal F_r\) for the following patches.

By \(L^1_{\mathrm{loc}}\) convergence and polar coarea,
one can choose radii \(R_n\to R\) where the traces of
\(E_n\) differ from those of \(C\) in circle length
\(o(1)\). To do this, choose annular widths tending
to zero sufficiently slowly that the volume mismatch
divided by the width tends to zero, and take a
good slice in each annulus. The BV slicing and gluing
formula applies at almost every radius. Replace \(E_n\)
inside \(B_{R_n}\) by \(D\), which equals \(C\) in those
annuli. The new tagged interface on the circle has
total length at most the length of the trace mismatch.
After returning to physical scale the seam thus costs
\(o(t_n)\). Each cell mass changes by \(O(t_n^2)\).
Two exterior flow parameters correct the masses exactly,
at cost \(O(t_n^2)\) for every unit-ball dual test.
They may be chosen outside a fixed disk about \(p\)
where all three cells still have positive mass. The
surjectivity and uniform cost are the same as for
the two-phase patch. Call the corrected partition
\(B_n^*\).
[Exterior mass repair and compact perturbations, AR1](THREE_CELL_ADAPTIVE_SURFACE_TWO_PHASE_REGULARITY.md#nearby-dual-tests-and-an-exterior-mass-correction)

Both the deleted and inserted local perimeters are
\(O(t_n)\). The compact transport estimate in (AR1)
also applies here with three labels and shows that
any attaining test \(G_n\) of \(B_n^*\) approaches
\(\mathcal M(A)\) in norm. Specifically, use a
direction field of \(B_n^*\) whose smoothed tensor
attains its value. Copy its directions on unchanged
interfaces to the original partition, pull back its
directions on the exterior flow supports, and use
zero directions on the deleted interfaces of \(A\).
The tensor norm difference is \(O(t_n)\) because
the changes are compactly supported with the stated
perimeter and flow estimates. Compactness of the
reference tensors and minimality of \(A\) then identify
every limit of \(G_n\) as attaining at \(A\).

Uniform local Lipschitz bounds for the weights of
unit tests justify freezing the Gaussian density and
weights at \(p\), with \(o(t_n)\) error. Weak convergence
and the zero boundary masses give
\(\sigma_{n,e}(B_{R_n})\to\eta_e(B_R)\). The inserted
interior lengths converge to \(P_e(D;B_R)\), since
near the circle \(D\) is a finite cone. For attaining
\(G_n\) at \(B_n^*\), the increment of \(\mathcal F_r\)
is at most \(W_{G_n}(B_n^*)-W_{G_n}(A)\).
After dividing by \(t_n\) its upper limit is at most the
maximum in (TC1) times \(\phi_2(p)\). For every fixed
\(G\in\mathcal M(A)\) the same difference is a lower
bound for the increment. This proves the asserted
first-order formula. Minimality gives (TC1).

## Equality tests and missing interface in the limit

The tensions of a test are a three-point semimetric.
For each \(G\), write \(w_{ij}=\lambda_i+\lambda_j\),
where \(\lambda_i=(w_{ij}+w_{ik}-w_{jk})/2\ge0\).
Put \(Q_i=\sum_{e\ni i}\eta_e\), which dominates
the Euclidean perimeter measure of \(C_i\) by BV
lower semicontinuity. It follows that
\[
 d(G):=\sum_e w_e(G)\eta_e(B_R)-
                \sum_e w_e(G)P_e(C;B_R)
       =\sum_i\lambda_i(G)\bigl[Q_i(B_R)-P(C_i;B_R)\bigr]
       \ge0.                                                \tag{TC2}
\]
Applying (TC1) to \(D=C\) shows that the compact
set \(\mathcal M_0=\{G\in\mathcal M(A):d(G)=0\}\)
is nonempty.

If \(C\) is constant, every \(\eta_e\) vanishes.
For an equality test (TC2) gives \(\sum_e w_e\eta_e(B_R)=0\),
and any term with zero weight already has \(\eta_e=0\).
Taking continuity radii tending to infinity gives the
claim.

If \(C\) is nonconstant, every equality test has
strictly positive weight on every pair used by \(C\).
If all its weights vanished then all \(\eta_e\) would
vanish, which is impossible by lower semicontinuity.
If \(w_{ij}=0\) and the other two weights are \(w>0\),
then \(\eta_{ij}=0\). Hence
\[
 \begin{split}
   \sum_e w_e\eta_e(B_R)
      &=w\bigl[Q_i(B_R)+Q_j(B_R)\bigr]\\
      &\ge w\bigl[P(C_i;B_R)+P(C_j;B_R)\bigr]\\
      &= \sum_e w_eP_e(C;B_R)+2wP_{ij}(C;B_R).
 \end{split}
\]
For an equality test this forces \(P_{ij}(C;B_R)=0\).
Thus no zero-weight tag of an equality test can be
used by the limiting cone. The argument detects
the extra cost of a collapsed pair even when the
ordinary label limit alone loses that information.

## Shortening rays without increasing another pair length

For a cone with two nonopposite rays, shorten the
two segments on the same tag by connecting their
endpoints on a smaller circle with a chord through
the convex sector. Keep the outer cells unchanged.
Both sectors have well-defined labels, so the
chord carries the original tag and reduces its length.

Suppose there are \(m\ge4\) rays, and write the
successive sector labels cyclically. Orient the label
cycle as \(1\to2\to3\to1\). A sector whose two
neighbors have the same label corresponds exactly to
a sign change between successive oriented transitions.
If there is any such sector, there are at least two,
since sign changes have even total count in a cyclic
sequence. At most one angular gap can be at least
\(\pi\), so one of these sectors has aperture below
\(\pi\). Its two rays have the same tag. Join their
endpoints on a smaller circle by a chord within
that sector and assign the separated central portion
to their common neighboring label. No other tagged
length changes, and this pair shortens strictly.

In the remaining case all transitions have the same
sign. The labels repeat the three-cycle, so \(m\)
is a multiple of three and \(m\ge6\). The average
span of three consecutive rays, measured as the sum
of their two intervening gaps, is \(4\pi/m<\pi\).
Choose a block of three within this span. The two
sector labels just outside the block agree. Work
inside a slightly larger angular wedge of aperture
less than \(\pi\), with no other rays. On a smaller
circle join the three endpoints to a common
interior point displaced toward the middle of the
wedge. The displacement has positive dot product
with each of the three unit tangents, so a
sufficiently small displacement shortens all three
segments. Convexity keeps the resulting star inside
the wedge. Its labels are given by the circle trace,
and the matching outside labels introduce no
interface on the wedge sides. Thus all three pair
lengths strictly decrease.

Each of these competitors makes all terms
\(P_e(D;B_R)-P_e(C;B_R)\) nonpositive and makes at
least one tag used by \(C\) strictly shorter. For a
test outside \(\mathcal M_0\), the expression in
(TC1) is at most \(-d(G)<0\). For a test in
\(\mathcal M_0\), its positive weights on the
shortened tags make the same expression strictly
negative. Compactness and continuity then give
a strictly negative maximum, contradicting (TC1).
The surviving cones are a constant partition,
a straight two-phase line or a three-label tripod.
The chord and star replacements are related to
the planar Steiner argument for fixed surface
energies. The exact comparison (TC1) addresses
the variable choice of attaining test and the
possible limiting interface excess here.
[Franceschi, Pratelli and Stefani, Theorem A and the cone replacements in Section 2.3](https://www.numdam.org/item/10.5802/jep.238.pdf)

## Tripod tensions and tagged convergence

For a tripod, let \(v_e\) be the outward unit
tangent of each ray. Perturb its common vertex
from zero to \(h z\), keeping the endpoints on
a strictly smaller circle fixed. These are
ordinary competitors; \(z\) is any planar
vector and \(h>0\) is small. The change in
tagged length is \(-h z\cdot v_e+O(h^2)\),
uniform for the finitely many tags for fixed \(z\).
In (TC1) any maximizing test for this patch has
\(d(G)=O(h)\). Taking limits of such tests as
\(h\downarrow0\) yields
\(\max_{G\in\mathcal M_0}(-z\cdot\sum_e w_e(G)v_e)\ge0\).
Finite-dimensional separation supplies a convex
combination of equality tests whose averaged
tensions \(\widehat w_e>0\) satisfy
\[
                       \sum_e \widehat w_e v_e=0.            \tag{TC3}
\]
All three triangle inequalities are strict for
these averaged tensions. Indeed an equality
would force the two other positively weighted
unit tangents to coincide, by the equality case
of the Euclidean triangle inequality. The three
rays are distinct. All three averaged coefficients
\(\widehat\lambda_i\) in (TC2) are therefore positive.
Equality for the convex combination now forces
\(Q_i(B_R)=P(C_i;B_R)\) for each \(i\).
Because \(Q_i\) dominates the perimeter measure
of \(C_i\), this is equality of measures inside
\(B_R\). The system \(Q_i=\eta_{ij}+\eta_{ik}\)
is invertible, so each \(\eta_e\) is the actual
tagged measure of the cone. Taking larger
continuity radii gives local tagged convergence
throughout the tangent plane.

The fixed pressure mixture of the original
minimizer also has force balance at a tripod
tangent. Its rescaled varifold is stationary
in the limit. Tagged convergence identifies the
limit with the three rays carrying their fixed
averaged weights at \(p\), times the common
Gaussian density there. These weights are
positive by the reciprocal bound, since each
tag has positive rescaled limit. The first
variation of each limiting ray at the origin
then gives (TC3) with the pressure-mixture
weights. In particular those weights also
obey the strict triangle inequalities.

## The only possible hidden interface at a line

Let \(C\) be the line separating \(i,j\), with
the third label \(k\) absent, and let \(\ell\)
be length measure on the line. Suppose there
is a genuine excess on a continuity ball in
the sense that the tagged limits are different
from the cone's tagged measures there. For
any \(G\in\mathcal M_0\), its line weight
\(w_{ij}>0\). If \(\lambda_k>0\) then (TC2)
forces \(Q_k=0\) on the ball. At least one of
\(\lambda_i,\lambda_j\) is positive and then
forces \(\eta_{ij}=\ell\) there. This would
give no excess. Thus \(\lambda_k=0\).
If either \(\lambda_i\) or \(\lambda_j\) also
vanishes, its pair with \(k\) has zero weight
and therefore zero limiting measure by (ASV4).
The other label has exact perimeter \(\ell\)
by (TC2), while lower semicontinuity for
the first label requires at least \(\ell\)
on the common pair. Again no excess results.

Consequently \(\lambda_i,\lambda_j>0\) and
\(\lambda_k=0\) for every equality test if
there is excess. In particular \(w_{ij}=w_{ik}+w_{jk}\)
and \(T_{\sqrt r}G_k(p)\) lies strictly between
the other two score vectors on a line. Equality
in (TC2) imposes \(Q_i=Q_j=\ell\) on the ball,
so the only possible tagged measures are
\[
            \eta_{ij}=\ell-\xi,\qquad
            \eta_{ik}=\eta_{jk}=\xi,\qquad 0\le\xi\le\ell .
                                                               \tag{TC4}
\]
For larger continuity balls the same argument
applies if excess is present on any smaller
ball. The measure \(\xi\) represents a collapsed
layer of the missing phase. Eliminating such
layers is an additional obligation for a global
regularity theorem for the adaptive problem.

## The tangent type at a point and isolation of tripods

At every point of the interface support the tangent
partitions are nonconstant, by the uniform raw density
bound and the tangent density alternative. At a fixed
point the line and tripod alternatives cannot occur
along different sequences. In fact the entire tangent
set of ordinary partitions is connected in
\(L^1_{\mathrm{loc}}\). The rescaling map
\(t\mapsto A(p+t\,\cdot)\) is continuous for \(t>0\)
in this topology, as follows by approximating the
indicators on compact sets by continuous functions.
The local perimeter bound gives relative compactness
of all sufficiently small rescalings. The tangent
set is the intersection of the nested compact
connected closures of the images of \((0,t)\), so
it is connected.

For a fixed \(p\) admitting a tripod tangent, the
pressure-mixture tensions there are positive and
satisfy strict triangle inequalities. Their force
equation determines the angles between its tagged
rays by the Euclidean cosine rule. Every tripod
tangent at \(p\) uses these same weights and angles.
The set of these possible tripod cones is a union
of compact rotation orbits, separated in
\(L^1_{\mathrm{loc}}\) from the compact set of
two-phase lines. This proves that every tangent
at this point is a tripod. If there is no tripod
tangent, all tangents are lines. In that case the
line uses the same pair of labels along every
sequence, since the line orbits for different
pairs are disjoint compact sets. Similarly, at
a tripod point the cyclic order of the labels is
the same along every sequence. These conclusions
do not require uniqueness of the orientation
angle.

Tripod points are isolated among tripod points.
To prove this, suppose distinct tripod points
\(p_n\) converge to a tripod point \(p\), and put
\(t_n=|p_n-p|\). At a tripod point the three weights
of the pressure mixture are positive, hence they
are bounded below near \(p\). The pressure
formula gives the first variation bound
\(|\delta V(X)|\le C\int |X|\,d\mu\) there for
the weighted varifold, where
\(\mu=\sum_e\bar w_e\nu_e\). The radial first
variation identity gives monotonicity of
\(e^{Ct}\mu(B(q,t))/t\) for all centers \(q\)
in a smaller neighborhood. Thus for any such
\(q\) and all small \(t\),
\[
 \mu(B(q,t))\ge
          e^{-Ct}\,t\,\Theta(q),\qquad
 \Theta(q)=\lim_{u\downarrow0}\mu(B(q,u))/u .
\]
For a tripod point, tagged convergence identifies
\(\Theta(q)=\phi_2(q)\sum_e\bar w_e(q)\).
Meanwhile a subsequence of \(t_n\) at \(p\) gives
a tripod tangent, and the scaled centers
\(u_n=(p_n-p)/t_n\) converge to a unit vector \(u\).
Choose \(a>0\) so small that \(B(u,a)\) intersects
at most one ray of this tangent. For \(a'<a\),
the balls \(B(u_n,a')\) eventually lie in \(B(u,a)\).
The density lower bound at \(p_n\) passes to the
limit and shows that the tangent weighted
measure of \(B(u,a)\) is at least
\(a'\phi_2(p)\sum_e\bar w_e(p)\). Let \(a'\uparrow a\).
The same measure is at most
\(2a\phi_2(p)\max_e\bar w_e(p)\), by the single-ray
intersection. Strict triangle inequalities at
\(p\) make this impossible. This isolation result
allows possible accumulation at a different
type of point. In particular a local topology
theorem still needs to address collapsed layers
and the possibility of third-phase incursions
at a line point.

The support itself converges locally with every
subsequential tangent as a closed set. The possible
hidden measure in (TC4) stays on the line, so the support
of the limiting total tagged measure is exactly the
line or the three tripod rays. Any point of that
support is approached by rescaled support points
by weak measure convergence. Conversely, let support
points of the rescaled partitions converge to a
point \(x\) in a fixed compact set. For each small
fixed \(a>0\), the compact-uniform raw lower bound
(AT7) applied at the original points gives a positive
lower bound proportional to \(a\) for the rescaled
mass within distance \(a/2\) of these support
points. Weak convergence on a slightly larger
closed ball shows that \(x\) is in the support of
the limiting measure. Thus no additional branch
disappears from the closed-set limit.
