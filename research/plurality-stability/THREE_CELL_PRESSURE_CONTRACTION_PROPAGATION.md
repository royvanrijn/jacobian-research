# Three-cell price bounds and contraction of occupied phases

Fix \(0<\rho<1\) and an actual planar Gaussian noise-stability
maximizer with three positive prescribed masses. Write
\(U_i=T_\rho1_{A_i}\) and use centered prices satisfying the
pointwise winner equations. The prices give a fixed positive
lower bound on the conditional probability of a winning label.
This yields quantitative phase recurrence at contractively
related centers and excludes every superlinear isolated graph
end in a full half-plane.
[Fixed-mass assignment and stationary winner equations](MAXIMIZER_PERTURBATIVE_RIGIDITY.md#assignment-and-the-exact-information-supplied-by-global-maximality)

## A positive own-label score on its winning set

Let \(U_i=T_\rho1_{A_i}\) be the unadjusted scores and center
the prices by \(\sum_i\beta_i=0\). Every actual three-cell
maximizer with positive masses has complete adjacency.
Take a regular interface of the two labels \(j,k\) distinct
from \(i\). Its adjusted scores satisfy
\(U_j-\beta_j=U_k-\beta_k\ge U_i-\beta_i\).
Since the three adjusted scores sum to one, this implies
\(U_i-\beta_i\le 1/3\). At every finite observation \(U_i>0\)
by positive mass and the strictly positive Gaussian kernel.
Therefore
\[
                \theta_i:=\tfrac13+\beta_i>0                 \tag{PP1}
\]
for every label. Conversely \(U_i-\beta_i\ge1/3\) at any
point where label \(i\) is an adjusted-score winner, by the
same sum over scores. Hence \(U_i\ge\theta_i\) on its full
closed winner set. Both inequalities include actual
interface points by continuity. The constraint (PP1) uses
only positivity of the three specified masses, in addition
to actual maximality.
[Complete adjacency of actual maximizers](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#complete-adjacency)

The three \(\theta_i\) are uniformly bounded below on
compact positive correlation bands and compact sets of
positive mass triples for actual maximizers in the active
plane. Otherwise use strong maximizer compactness at
limiting parameters. The centered prices are bounded
because any point in cell \(i\) gives
\(\beta_i-\beta_j\le U_i-U_j<1\), with the reverse
bound obtained from cell \(j\). Pass to convergent prices.
Local uniform convergence of the smoothed scores and
strong label convergence pass the winner rule to the limit.
The limit is an actual maximizing partition with positive
masses, so the limiting prices must obey (PP1).
[Near-maximizer compactness at varying correlations and masses](GAUSSIAN_NEAR_MAXIMIZER_COMPACTNESS.md)

## Contractive recentering and isolated ends

For every winning point \(x\) of label \(i\), the Gaussian
source distribution with mean \(\rho x\) assigns \(A_i\)
probability at least \(\theta_i\). Choose a fixed \(R\) with
\(\mathbb P(|\sqrt{1-\rho^2}Z|>R)<\min_i\theta_i/2\).
The source density is bounded above, so there is
uniformly positive Lebesgue area of \(A_i\) in
\(B_R(\rho x)\). This is (PP2). The same holds for each incident label at
any actual interface or higher junction point. The constants
are uniform on the compact parameter bands of (PP1).
Since Gaussian and Lebesgue null sets agree, these sets
of positive area contain strict winner points.

This recurrence can be followed through any number of
successive contractions. Starting at a winning point \(x\),
choose strict winner points \(x_\ell\) of the same label, with
\(x_0=x\) and \(|x_{\ell+1}-\rho x_\ell|\le R\). Then
\[
        |x_\ell-\rho^\ell x|\le R/(1-\rho).
\]
The same estimate shows that the ball of radius \(R/(1-\rho)\)
about \(\rho^\ell x\) contains a positive area of the label
bounded uniformly from below for every \(\ell\ge1\):
use the area in \(B_R(\rho x_{\ell-1})\), which is contained
in that ball.

The recurrence holds for **each** incident label at an active
two-label arc, starting at every point of the arc or nearby
on either of its strict sides. At a triple point it holds
for all three labels. The resulting source ball centered at
\(\rho^\ell x\) contains the specified positive-area pieces
simultaneously. Two distinct recurring labels force a boundary
of the source partition in the ball, though a third label can
separate the two pieces there. No ordering of source normals
is asserted.

In particular suppose an isolated binary end in a full
half-plane \(s>S\) is given by a continuous graph
\(y=h(s)\), with one label filling its lower side and the
other its upper side, up to null sets. Then
\[
                            |h(s)|=O(1+s)                    \tag{PP3}
\]
as \(s\to+\infty\). For \(s\) large, apply the fixed-radius
source recurrence to both incident labels at \((s,h(s))\).
Their source points lie in the same ball in the half-plane.
The continuous graph crosses the segment joining two such
points on opposite sides, giving a wall point
\((s',h(s'))\) within \(R\) of \((\rho s,\rho h(s))\).

Iterate until the first horizontal coordinate \(s_\ell\)
at or below a fixed cutoff \(L\) on the positive tail.
Choose \(L\) so that \(\rho L-R\) exceeds both \(S\) and
\(R/(1-\rho)\), and that \(\rho s+R<s\) for \(s>L\).
The iteration terminates since \(s_\ell\le\rho^\ell s+
R/(1-\rho)\), and its terminal horizontal coordinate
belongs to the fixed compact interval \([\rho L-R,L]\).
The recurrent estimates imply
\[
  |\rho^\ell s-s_\ell|\le R/(1-\rho),\qquad
  |\rho^\ell h(s)-h(s_\ell)|\le R/(1-\rho).
\]
The first bounds \(\rho^\ell s\) below by a positive
constant, and the second bounds \(\rho^\ell |h(s)|\)
above by a constant by continuity on the terminal interval.
This proves (PP3). It rules out superlinear isolated ends,
without proving that an at-most-linear end has a limiting
slope.
[Boundedness and straightening for sublinear isolated graphs](THREE_CELL_ISOLATED_END_RIGIDITY.md)

## The outer cell at a distant circular junction

Suppose a noncollinear triple \(p\) has the circular-wall
geometry in (RJ2) in a partition without a common centered
circle for all three winners. Write \(r_p=|p|\) and let \(j\)
be the outer label of the circular wall. That geometry
excludes \(j\) from the full inner annulus
\(r_p-\varepsilon<|x|<r_p\), and supplies strict \(j\)
points just outside the wall. When \((1-\rho)r_p>R\)
for the one-step source radius in (PP2), the source ball
\(B_R(\rho p)\) lies strictly inside radius \(r_p\). Its
positive area of label \(j\) must lie at radii at most
\(r_p-\varepsilon\), by the excluded annulus. Thus \(j\)
has a bounded positive-area group separated from its
outer strict winner region by the entire inner annulus.
This identifies the disconnected label in this distant
circular-junction regime.
[Inner and outer annuli in (RJ2)](THREE_CELL_NO_COMMON_CIRCLE_JUNCTIONS.md#noncollinear-junctions-produce-radial-separation)
