# No dormant wall in the central disk of a three-cell maximizer

Let \(A\) be an ordinary balanced planar three-cell Gaussian
noise-stability maximizer at \(339/500\le\rho<1\).
A true two-label interface along which the third adjusted
score also ties cannot meet the closed disk of radius three
about the Gaussian origin. In this statement “along” means
on an open interface arc. Isolated three-score ties can occur.

Consequently, on each compact correlation interval in
\([339/500,1)\), the number of all three-score tie points
in this disk is uniformly bounded over balanced maximizers.
On each strictly smaller disk the numbers of intersecting
strict-winner components and local graph faces are also
uniformly bounded. These conclusions require no assumption
on connectedness of the full cells or on the exterior graph.
In particular they apply at any nonsimplex right endpoint
in the remaining intermediate range.

## Long dormant arcs

Write \(\sigma=\sqrt{1-\rho^2}\),
\(\kappa=(1-\rho)/(2(1+\rho))\), and
\(P(h)=2\Phi(h)-1\), with \(\Phi\) the standard Gaussian
cumulative distribution. Let \(\Gamma\) be the interface
component of an open arc where all three adjusted scores tie
and just two labels win on the two sides. The continuation
argument for common-score walls makes \(\Gamma\) a proper
regular line or a simple closed regular curve, with the
third label dormant throughout. The insertion variation on
this component gives, for each nonzero nonnegative compact
wall speed \(f\),
\[
      \int_\Gamma a f^2\,d\nu
         >6\iint_{\Gamma^2} M_\rho(x,y) f(x)f(y)
                                  \,d\nu(x)d\nu(y).           \tag{CD1}
\]
Here \(a\) is the magnitude of the normal pair-score gradient,
\(d\nu=\phi_2\,ds\), and \(ds\) is arclength.
[Continuation and dormant-label insertion, NT1–NT2](THREE_CELL_CONNECTED_NETWORK_TOPOLOGY.md#a-dormant-score-along-an-entire-interface-component)

Suppose \(\Gamma\) contains \(p\) and has an arclength
parametrization \(x(s)\), \(-h\le s\le h\), with \(x(0)=p\)
and no repeated points except possibly the two endpoints
of a closed curve. A bounded speed supported on this segment
is allowed in (CD1) by approximation. Indeed the local term,
interaction term and strip insertion's volume vector extend
continuously along bounded compact approximations. In (NT1)
the volume of the inserted third phase to first order is
\(\int f\,d\nu>0\), so the strictly positive translation
term persists for this nonzero nonnegative limit.

Test with \(f(x(s))=\sqrt{\phi_2(se_1)/\phi_2(x(s))}\) on
the segment and zero elsewhere. The universal derivative
bound \(a\le2\rho\phi_1(0)/\sigma\) gives
\[
            \int_\Gamma af^2\,d\nu
                \le \frac{2\rho\phi_1(0)^2}{\sigma}P(h).
                                                               \tag{CD2}
\]
For the interaction write
\(k_\rho(x,y)=M_\rho(x,y)\sqrt{\phi_2(x)\phi_2(y)}\).
Its exponential factor is
\[
 \exp\!\left[-\frac{\kappa}{2}(|x|^2+|y|^2)
                         -\frac{\rho|x-y|^2}{2\sigma^2}\right].
\]
The chord and norm bounds for the parametrization give
\[
 k_\rho(x(s),x(t))
 \ge k_\rho(se_1,te_1)
                \exp[-\kappa(|p|^2+|p|(|s|+|t|))].
\]
On the centered reference line the Mehler factor is the
one-dimensional kernel divided by \(\sigma\). Thus the
interaction in (CD1), before its factor six, is at least
\[
 \frac{\phi_1(0)^2}{\sigma}\,
 \mathbb E\!\left[
    1_{\{|S|,|T|\le h\}}\,
     e^{-\kappa(|p|^2+|p|(|S|+|T|))}
 \right],                                                     \tag{CD3}
\]
where \(S,T\) are standard one-dimensional Gaussians of
correlation \(\rho\).
[Arclength kernel comparison, PA3–PA5](THREE_CELL_PAIR_ARC_CENTRAL_EXCLUSION.md#comparison-by-arclength)

The probability of the truncation in (CD3) is at least
\(P(h)^2\), by positivity of the one-dimensional Mehler
operator on the centered interval indicator. Also
\[
 \mathbb E(|S|\mid |S|,|T|\le h)\le\mathbb E|S|
       =\sqrt{2/\pi}.
\]
Indeed the conditional probability of \(|T|\le h\)
given \(S=s\) is even and decreases with \(|s|\).
After first conditioning on \(|S|\le h\), reweighting
by this decreasing function can only decrease the mean
of \(|S|\). The same argument applies to \(T\). Jensen's
inequality in (CD3) and comparison with (CD2) imply
\[
        e^{-\kappa(|p|^2+2\sqrt{2/\pi}|p|)}
                        <\frac{\rho}{3P(h)} .                \tag{CD4}
\]
A proper regular line supplies the segment for every
\(h>0\), and a closed curve of length at least \(2h\)
supplies it for that value of \(h\).

Take \(h=3/2\), \(R=3\), \(\rho_0=339/500\), and
\(\kappa_0=161/1678\). The left side of (CD4) for
\(|p|\le R\) is at least \(e^{-\kappa(69/5)}\), since
\(\sqrt{2/\pi}<4/5\). The function
\(e^{-\kappa(69/5)}/\rho\) increases on
\([\rho_0,1)\). Its logarithmic derivative equals
\((69/5)/(1+\rho)^2-1/\rho>0\) on that range.
At \(\rho_0\),
\[
  \kappa_0(69/5)=11109/8390
       <\log(430/113)<\log\bigl(3P(3/2)/\rho_0\bigr).
\]
For the second inequality, \(P(3/2)>43/50\):
integrating the fifth odd Taylor truncation of
\(e^{-s^2/2}\) on \([0,3/2]\) and using
\(\sqrt{2/\pi}>797/1000\) gives the rational lower bound
\(34928203809/40370176000>43/50\). For the first
inequality, the first four positive terms of the series
\(2\sum_{n\ge0}z^{2n+1}/(2n+1)\) at \(z=317/543\)
already exceed \(11109/8390\). These bounds contradict
(CD4) for every \(\rho\ge\rho_0\).

## Short closed dormant walls

If the whole dormant component has length \(2h<3\),
the arclength mean-square radius \(m_2\) satisfies
\[
 m_2\le |p|^2+|p|h+h^2/3<57/4 .
\]
This integrates the pointwise bound
\(|x(s)|^2\le(|p|+|s|)^2\) over the whole curve, with
\(|p|\le3\). The closed-wall bound (NT4) requires
\[
 \kappa m_2>\log\!\left(
          3P(e\sqrt{\pi/(2\rho)})/\rho^{3/2}\right)
       \ge \log(8/3)+6\kappa.
\]
Indeed \(e\sqrt{\pi/(2\rho)}>2\) and \(P(2)>8/9\),
the latter by the one-dimensional Gaussian tail bound
\(1-P(2)\le\phi_1(2)<1/9\). Also
\(-\log\rho\ge 2(1-\rho)/(1+\rho)=4\kappa\).
But \(\kappa(m_2-6)<\kappa_0(33/4)=5313/6712
<\log(8/3)\), a contradiction.
[Closed dormant-wall perimeter bound, NT4](THREE_CELL_CONNECTED_NETWORK_TOPOLOGY.md#a-radius-bound-for-compact-dormant-walls)

## Uniform counts of ties and local components

The disk \(\overline B_3\) lies in the uniform
maximal-tie gradient regime of the central compactness
theorem. In its notation,
\(\kappa R_*(\rho)^2=1-\log 2-2\log\rho
\ge1-\log 2+8\kappa>9\kappa\) for
\(\rho\ge\rho_0\). If the number of three-score ties in
\(\overline B_3\) were unbounded on a compact correlation
interval, take a limiting maximizer, limiting prices and
an accumulation point of these ties. Maximal pair gradients
are distinct there by the same central estimate. A nonzero
transverse derivative of one pair score gives a holomorphic
equality graph in a fixed small neighborhood, with
convergence of the nearby graphs for the sequence.
Restrict a second pair score to those graphs. The
unbounded number of zeros accumulating at the point
forces the limiting restriction to vanish identically,
by the holomorphic zero count on a smaller disk.
This is the central junction-count argument, which
uses only the two equality equations and therefore
counts all three-score ties here.
[Uniform central gradient separation and local holomorphic zero count](THREE_CELL_CENTRAL_INTERFACE_COMPACTNESS.md#collapse-of-components-and-accumulation-of-junctions)

The three limiting gradients on a small equality arc
are distinct and differ only in the transverse direction.
The two extreme transverse derivatives give the strict
winners on the respective sides. This is an actual
pair arc with a dormant third score, whose continuation
is excluded by the estimates above. Hence the total
central tie count is uniformly bounded. The central
gradient theorem also gives a uniform bound on the
arclength of the jump graph in \(\overline B_3\).
[Uniform bound on central graph length](THREE_CELL_CENTRAL_INTERFACE_COMPACTNESS.md#uniform-separation-at-every-maximal-tie)

For completeness fix \(R<3\). The bounded-component
diameter criterion (IC5) gives a uniform positive
diameter for pockets surrounded by one label and
contained in \(B_3\) on any compact band in
\([\rho_0,1)\). To verify its strict threshold, the
function \(e^{-9\kappa}/\rho^2\) increases on
\([\rho_0,1)\), since
\(9\rho-2(1+\rho)^2=(2\rho-1)(2-\rho)>0\).
At \(\rho_0\) it exceeds \(4/(3e)\), because
\(9\kappa_0<1\) and \(4\rho_0^2/3<1\).
Continuity gives (IC5) on a neighborhood of each
correlation, and a finite cover of the compact band
gives one positive diameter lower bound.
[Minimum diameter of a surrounded component, IC5](THREE_CELL_ISOLATED_SUBCLUSTERS.md#minimum-diameter-in-a-central-disk)

By coarea and the uniform bound on central length, one
can choose a radius \(s\in(R,3)\) at which the circle
is transverse to all regular jumps, avoids the genuine
junctions and meets the graph in a uniformly bounded
number of points. Indeed the integral of the count
over \(s\in(R,3)\) is at most the central length.
Critical radii of the analytic edges and radii of
junctions form a measure-zero set. The graph restricted
to the closed disk of radius \(s\), together with the
circle boundary, is a finite embedded graph. Vertices
come from junctions, circle intersections, and auxiliary
vertices on any remaining closed loops.

A closed loop with no junction lies on a single pair
interface. Its outer side has one fixed label locally
along the loop. The component of that label's open
complement abutting the loop on the inner side is
contained in the Jordan interior and is surrounded by
the outer label, even if it has holes. The Jordan
interior is inside the convex hull of the loop and
its diameter is at most half the loop length. Thus
each such loop uses at least twice the uniform
minimum pocket diameter in central arclength. This
bounds the number of auxiliary vertices. Bounded
vertex valence then bounds the number of edges, and
Euler's formula bounds the number of local faces.

Inside each face the partition has a single almost
everywhere constant label by the local definition of
the jump graph. Points of the face where that label is
not a strict winner are isolated dormant ties and can
be avoided by path detours. Thus a local face cannot
meet two distinct strict-winner components. Every
component meeting \(B_R\) hits some local face,
which proves the asserted uniform bound.
[Analytic graph topology and face labels](THREE_CELL_CONNECTED_NETWORK_TOPOLOGY.md#finitely-many-strict-winner-components-suffice-for-a-finite-graph)
[Isolation of dormant interior ties, FC4](THREE_CELL_FIRST_CONTACT_VARIATIONS.md#point-transfers-at-a-dormant-interior-score-tie)
