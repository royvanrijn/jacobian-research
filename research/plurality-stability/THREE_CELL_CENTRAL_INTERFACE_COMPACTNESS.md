# Uniform central interface regularity for three-cell maximizers

Fix \(0<r<1\) and a balanced planar Gaussian noise-stability
maximizer \(A\). Use the adjusted scores \(V_i=T_r1_{A_i}-\beta_i\)
with \(\sum_i\beta_i=0\), and let \(O_i=\{V_i>\max_{j\ne i}V_j\}\).
Define
\[
 R_*(r)^2=\frac{2(1+r)}{1-r}(1-\log 2-2\log r).
\]
At any point of \(\overline B_{R_*(r)}\), all gradients of the
maximal tied scores are distinct. This includes a tied label
with no strict wins nearby. The gradient separation is uniform on
closed central disks and compact positive correlation intervals whose
parameters satisfy the radius bound.

Consequently, arbitrarily small strict-winner components in these
disks can occur in a sequence of maximizers only by collapsing
to a triple equality with collinear gradients. An unbounded number
of genuine triple junctions in a fixed central disk would force a
limit maximizer with three scores identically tied along an open
analytic arc on which only two labels meet.

## Inserting a missing label next to a regular wall

Suppose \(V_i(p)=V_j(p)=V_k(p)\),
\(\nabla V_i(p)=\nabla V_j(p)\ne\nabla V_k(p)\), and \(O_j\) misses
a neighborhood of \(p\). Set
\[
 H=\nabla^2(U_i-U_j)(p),\quad U_i=T_r1_{A_i},
 \qquad n=\frac{\nabla(V_i-V_k)(p)}{|\nabla(V_i-V_k)(p)|}.
\]
For a vector \(z\) with \(n\cdot z>0\), \(p+hz\) lies above the
\(i,k\) score wall for all sufficiently small \(h>0\).
Since \(j\) is absent and the tie sets have measure zero,
continuity gives \(V_i-V_j\ge0\) throughout this open side.
Its value and gradient vanish at \(p\), so \(z^{\mathsf T}Hz\ge0\).
The quadratic form is even, hence \(H\succeq0\).

Let \(P\) be an invertible real matrix with determinant one
and let \(E_R=p+R\{z\in P B_1:n\cdot z>\epsilon\}\)
for a fixed small \(\epsilon>0\). For small \(R>0\), \(E_R\) lies in
\(A_i\) up to null sets. Put \(t=\gamma_2(E_R)\).
Translate the whole partition by \(w_R\) and transfer
\(E_R+w_R\) to \(j\). The exact balance equations have a solution
\(w_R=t v+o(t)\). Indeed their derivative in \(w\) at \(R=0,w=0\)
is the invertible translation volume map onto the centered label
plane, by full moment rank. For a fixed cap the extra volume is a
smooth integral in \(w,R\) up to \(R=0\), using Gaussian density in
rescaled coordinates. The vector \(v\ne0\) is the translation
velocity compensating a unit transfer \(e_j-e_i\).

Apply the priced objective \(S_r-2\sum\beta_i\gamma_2(A_i)\).
The translation has zero linear term and contributes
\((r^{-1}-1)\|s_v\|_a^2t^2+o(t^2)\), by the translation identity.
The transfer contributes a quadratic self-interaction
\(2M_r(p,p)t^2+o(t^2)\). Its linear cost divided by \(t^2\)
has the same limit as
\(2t^{-2}\int_{E_R}(V_i-V_j)\,d\gamma_2\).
In fact translating a source by \(w_R\) shifts the score argument
by \(-w_R/r\), while the points of \(E_R\) are shifted by \(w_R\).
Because \(\nabla(V_i-V_j)(p)=0\), the resulting score error is
\(O(R|w_R|+|w_R|^2)=o(t)\) on the cap.

Let \(R\) tend to zero and then let \(\epsilon\) tend to zero.
The half ellipse has area \(\pi/2\), and its even second moments
are exactly half those of \(P B_1\). The limiting normalized
linear cost is \(\operatorname{tr}(P^{\mathsf T}HP)/(2\pi\phi_2(p))\).
Maximality and minimization over \(\det P=1\) give
\[
 \frac{\sqrt{\det H}}{\pi\phi_2(p)}
 \ \ge\ 2M_r(p,p)+(r^{-1}-1)\|s_v\|_a^2 .                  \tag{CC1}
\]
The infimum of the trace is \(2\sqrt{\det H}\), including
zero for singular \(H\). In particular \(H\) is positive definite.
This is the half-ellipse version of the interior point transfer.
[Interior point transfer, FC4](THREE_CELL_FIRST_CONTACT_VARIATIONS.md#point-transfers-at-a-dormant-interior-score-tie)
[Translation identity and full rank](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#translation-equations-on-an-unconstrained-space-of-speeds)

Differentiating the conditional convolution yields
\(\operatorname{tr}H\le4r^2/[e(1-r^2)]\), as in (FC5).
Use \(\sqrt{\det H}\le \operatorname{tr}H/2\) and
\[
 2\pi\phi_2(p)M_r(p,p)
       =(1-r^2)^{-1}\exp\!\left(-\frac{1-r}{2(1+r)}|p|^2\right).
\]
The positive translation term in (CC1) implies
\[
 \exp\!\left(-\frac{1-r}{2(1+r)}|p|^2\right)<2r^2/e.
                                                               \tag{CC2}
\]
Hence the assumed endpoint coincidence is excluded for
\(|p|\le R_*(r)\).

## Uniform separation at every maximal tie

Suppose two maximal scores at \(p\) have coincident gradients.
If their two strict winners both accumulate at \(p\), either
there is a genuine triple junction or some regular interface
between those two labels accumulates there. The first case
violates distinctness at a triple junction. In the other case
the reciprocal-gradient bound (G20) forbids the vanishing
pair gradient. If only one strict winner accumulates at \(p\)
among all three labels, the ball about \(p\) is filled by that
label up to null sets. Apply (FC5) to it and a distinct tied
label. Its excluded radius is larger than \(R_*(r)\).

The remaining case has exactly two strict winners accumulating
at \(p\). One is among the tied pair with coincident gradients,
and the other is the third label. All scores tie at \(p\).
The two winners have distinct gradients by (G20) on their
accumulating pair interfaces, and the missing label is absent
in a neighborhood. Apply (CC2).
[Regular arcs and genuine triple junctions, G20](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#regular-interface-arcs-and-isolated-triple-junctions)

For uniformity take a compact \(J\subset(0,1)\) and
\(R\le\min_{r\in J} R_*(r)\). Balanced planar maximizers are
strongly compact along convergent correlations in \(J\).
Normalized prices lie in a compact set since every price
difference is less than one in magnitude. The smoothed scores
and every fixed spatial derivative converge locally under
this compactness, after taking a price subsequence. The
Mehler integral also gives locally uniform convergence of
their complex extensions. For complex observation points
in a compact set, the absolute conditional kernels and their
derivatives have an integrable Gaussian majorant in the source
point, uniform on \(J\). Strong convergence of the indicators
and dominated convergence give the stated limits.

The limiting adjusted scores are prices for the limit maximizer.
Indeed the winning inequalities pass to the limit almost
everywhere, using a further almost everywhere convergent subsequence
of the indicators. The limiting pair contrasts are nonconstant
by Gaussian injectivity, so the tie sets are null. If the gradient
separation failed, a convergent sequence of positions in
\(\overline B_R\) and of their maximal tied scores would give
a coincident pair at a limiting maximal tie. The preceding
exclusions make this impossible. Thus some constant
\(c(J,R)>0\) bounds all maximal pair gradient differences below.
[Strong compactness](POSITIVE_CORRELATION_CONTINUATION.md)

For each fixed pair, balls centered on its active graph in the
disk can be taken with a common positive radius on which its
pair score is strictly monotone in a fixed directional coordinate.
This follows from the uniform gradient lower bound and uniform
score Hessian upper bound. The zero set in such a ball lies on
a single graph over the perpendicular coordinate, with a
uniform slope bound. A bounded number of these balls covers
the pair's active graph by a radius-separated center selection
in the disk. In particular the active graph has uniformly bounded
length in the disk. The graph bound allows inactive portions of the
same pair equality curve.

## Collapse of components and accumulation of junctions

Suppose a sequence of connected strict-winner components
of label \(i\) meets \(\overline B_R\) and has diameter tending
to zero. Pass to a limiting maximizer, limiting prices, and
a limiting point \(p\). Label \(i\) has maximal score there.
If a vector \(v\) has strictly positive dot product with
all \(\nabla(V_i-V_j)(p)\) for \(j\ne i\) tied at \(p\), the
same directional derivatives remain positive in a fixed
ball for all sufficiently late members of the sequence.
The originally untied scores are separated there by
continuity. Starting inside the shrinking component,
the segment in direction \(v\) to the edge of a fixed
smaller ball remains in the same strict-winner component.
This contradicts the shrinking diameter. The same argument
works when \(i\) is the sole maximal score, since an
entire fixed neighborhood then consists of strict wins.
By the central separation result, every tied pair difference
has nonzero gradient. Therefore two competitor gradients
relative to \(i\) must occur and be negatively proportional.
Equivalently the three adjusted gradients are collinear,
with that of \(i\) strictly between the other two.

Suppose instead that a sequence of maximizers in this band
has an unbounded number of genuine triple junctions in
\(\overline B_R\). Pass to a converging subsequence and a
point \(p\) in the disk such that every neighborhood contains
arbitrarily many of the triple junctions in individual members
of the subsequence. Such a point exists by compactness of the
disk and a finite covering argument. All limit scores tie at
\(p\), and the central result separates their gradients.
Choose coordinates in which one limiting pair difference
has nonzero transverse derivative at \(p\). The holomorphic
implicit function theorem, or Rouché's theorem in the transverse
coordinate on a small product of disks, gives a common
longitudinal disk on which the corresponding near-\(p\)
equality curves are holomorphic graphs converging to the
limit graph. Restrict a second pair contrast to these graphs.
If the limiting restriction were nonzero as a holomorphic
function, Rouché's theorem on a smaller disk would bound the
number of its nearby zeros in the sequence. Thus the three
limiting adjusted scores agree identically on the local curve.
All three gradient differences are normal to this curve and
distinct on a small arc. On the two sides near that arc the
strict winners are the two extreme gradients, respectively.
The middle label is tied along the arc and absent in the
adjacent open regions.

These conclusions isolate two necessary central degenerations.
They allow genuine common-tangent junctions and a middle score
identically tied along a two-label wall. No connectivity or
uniform junction-count assertion beyond the stated implications
is used.
