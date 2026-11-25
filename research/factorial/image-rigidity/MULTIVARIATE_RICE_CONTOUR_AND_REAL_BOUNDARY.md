# A single Rice contour for the multivariate moments

For distinct nonnegative integer exponent columns \(a,b\), delete
equal rows and put \(d_i=b_i-a_i\), \(x_i(q)=a_i+d_i q\),
\(m_i=(a_i+b_i)/2\), and \(\Delta=\sum_i d_i\).
The ordinary moment polynomial at order \(n\) is
\[
 M_n(t)=\sum_{k=0}^n\binom nk t^k
                  \prod_i\Gamma(1+n a_i+d_i k).
\]
There is an exact one-variable contour for every such support,
including rows with a zero exponent. On a compact part of the
open upper half-plane its Stirling phase has stationary equation
\(X(q)=-\operatorname{Log}\zeta\) at the scaled parameter
\(\zeta=n^\Delta t\). The corresponding lower phase has no
stationary point when \(\Im\zeta>0\). For fully shared supports,
we also identify all real boundary intervals through which an
upper descending contour can pass.

## Exact contour and its phases

Let \(\Im\zeta>0\) and \(\arg(-\zeta)=\arg\zeta-\pi\).
With \(q=z/n\), write
\[
 \mathcal I_n(q;\zeta)=
 \frac{\Gamma(n(q-1))}{\Gamma(nq+1)}
 \prod_i\frac{\Gamma(1+n x_i(q))}{\Gamma(1+n m_i)}
 \exp\!\bigl(nq[\operatorname{Log}(-\zeta)-\Delta\log n]\bigr).
\]
For any positively oriented contour \(\Gamma\) surrounding
\(0,1/n,\ldots,1\) and no other poles,
\[
 M_n(n^{-\Delta}\zeta)
  =\frac{(-1)^n n!\,n\prod_i\Gamma(1+n m_i)}{2\pi i}
              \int_\Gamma\mathcal I_n(q;\zeta)\,dq.          \tag{1}
\]
Indeed \(\Gamma(z-n)/\Gamma(z+1)=1/[z(z-1)\cdots(z-n)]\).
The residue of this ratio at \(z=k\) is
\((-1)^{n-k}/[k!(n-k)!]\), which gives (1).
This is the usual finite-difference contour applied to the
factorial sequence. [Flajolet and Sedgewick, equation 2](https://doi.org/10.1016/0304-3975(94)00281-M)

For example, the contour is the difference of two vertical
integrals on \(\alpha_n<0\) and \(\beta_n>1\), both taken
upward, with the integral on the right minus the one on the
left. Choose \(\alpha_n,\beta_n\) between the terminal
coefficient poles and the nearest poles of the Gamma factors.
Such choices exist for private rows as well, since the first
Gamma poles lie strictly beyond \(0\) or \(1\) at distances
of order \(1/n\). The horizontal sides can be sent to infinity
for every fixed \(n\) if \(D=\sum_i|d_i|\ge2\).
For \(q=x+iy\) in the bounded real range of the rectangle
and fixed \(n\), the Gamma factors have exponential rate
\(-n\pi D|y|/2\) in both directions by vertical Stirling
estimates. The factor with \(\operatorname{Log}(-\zeta)\)
adds \(n(\pi-\arg\zeta)y\), strictly smaller in absolute
value than \(n\pi D|y|/2\) for nonzero \(y\).
The rational Gamma ratio and powers of the radii contribute
only polynomial factors in \(|y|\) for fixed \(n\).
The one-row unit-gap case \(D=1\) still has the finite
contour (1).

Write the principal logarithms for \(q,1-q,x_i(q)\) in
each open half-plane. Define there
\[
 \phi(q;\zeta)=
 \sum_i\bigl(x_i(q)\log x_i(q)-m_i\log m_i-x_i(q)+m_i\bigr)
 +(q-1)\log(1-q)-q\log q+q\operatorname{Log}\zeta .
\]
On a compact subset of the upper half-plane, after
multiplication by a nonzero scalar, the integrand is
\[
 A(q)e^{n\phi(q;\zeta)}\bigl(1+O((1+W)/n)\bigr),\qquad
 A(q)=\frac{\exp(\tfrac12\sum_i\log(x_i(q)/m_i))}
                 {\sqrt{q}\sqrt{q-1}},\quad
 W=\sum_i d_i^2/m_i .                                      \tag{2}
\]
The square roots in \(A\) take the principal branches.
On the lower half-plane the analogous formula uses
\(\phi_-(q;\zeta)=\phi(q;\zeta)-2\pi i q\).
All estimates are relative to the displayed amplitude. The
constants in the error are uniform in the exponent profile and
dimension on compact \(q\)-sets and bounded \(W\)-bands.
The scalar can be chosen the same in the two half-planes
at integer \(n\).

To see this, the Stirling formula for
\(\Gamma(n(q-1))/\Gamma(nq+1)\) gives, up to a constant,
the amplitude \(1/\sqrt{q(q-1)}\) with the specified two
branches, and phase
\((q-1)\log(q-1)-q\log q\).
In the upper half-plane
\(\log(q-1)=\log(1-q)+i\pi\). Combining this with
\(\operatorname{Log}(-\zeta)=\operatorname{Log}\zeta-i\pi\)
adds the constant \(-i\pi\) to the phase.
In the lower half-plane the first sign reverses, adding
\(-2\pi i q+i\pi\) instead. Both constants exponentiate to
\((-1)^n\).
For the row factors use
\(\log\Gamma(1+z)=(z+1/2)\log z-z+\tfrac12\log(2\pi)+R(z)\),
where \(R(z)=O(|z|^{-1})\) on closed sectors away from the
negative real axis. Cauchy's estimate on a slightly wider
sector gives \(R'(z)=O(|z|^{-2})\) away from zero.
For \(q\) in a compact set of either open half-plane, the
points \(x_i(q)/m_i=1+(2q-1)d_i/(2m_i)\), and their segments
to \(1\), stay in a compact part of the slit plane because
\(d_i/(2m_i)\in[-1,1]\). The segment is uniformly separated
from zero. One therefore has
\[
 R(n x_i(q))-R(n m_i)=
         O\left(\frac{|d_i|}{n m_i^2}\right).
\]
Each \(m_i\ge1/2\) and each \(|d_i|\ge1\), so these errors sum
to \(O(W/n)\); the two denominator factors add \(O(1/n)\).
The remaining row terms in powers of \(n\) combine with the
external scaling, since
\(\sum_i(x_i(q)-m_i)-\Delta q=-\Delta/2\).
This leaves only a nonzero scalar.
[DLMF Stirling expansion](https://dlmf.nist.gov/5.11.E1)

Differentiating gives
\[
 \phi'(q;\zeta)=X(q)+\operatorname{Log}\zeta,\qquad
 X(q)=\sum_i d_i\log x_i(q)+\log(1-q)-\log q.
\]
Consequently the [upper-sheet count](MULTIVARIATE_SADDLE_MAP_TWO_SHEETS.md)
applies directly to the full upper half-plane Rice phase,
without the strip restriction of the conditional integral.
In the lower half-plane \(\Im X(q)<\pi\), since
\(d_i\Im\log x_i(q)<0\) and
\(\arg(1-q)-\arg q\in(0,\pi)\). It follows that
\(\Im\phi_-'(q;\zeta)<\arg\zeta-\pi<0\).
In particular the lower phase has no stationary point
and its real part strictly decreases as one moves
vertically down. At an upper stationary point the phase,
up to a \(q\)-independent constant, is
\(F_a(q)-\Delta q\). To check this, substitute the stationary
equation in \(q\operatorname{Log}\zeta\) and use
\(x_i(q)-q d_i=a_i\). Thus the phase comparison is the
same as for the conditional representation.

## Where the upper descending contours can meet the real axis

The following boundary analysis assumes all \(a_i,b_i>0\).
Put \(p_L=\max_{d_i>0}(-a_i/d_i)\) and
\(p_R=\min_{d_i<0}(-a_i/d_i)\), with \(-\infty,+\infty\)
for empty indexing sets. On the real intervals
\((p_L,0)\) and \((1,p_R)\), the upper limit of
\(\Im X\) is \(-\pi\). On the other real intervals,
excluding the branch points, it is nonnegative. For the
latter intervals the descent flow for \(\Re\phi\) points
strictly into the upper half-plane at the real boundary:
its vertical component is
\(-\partial_y\Re\phi=\Im\phi'=\Im X+\arg\zeta>0\).
A descending trajectory therefore cannot exit there.

A descending trajectory also cannot land at a real branch
point. To prove this, fix \(b\in\{0,1,p_i\}\) and write
\(q=b+z\), \(z=x+iy\) in the upper half-plane. Near \(b\),
\[
 \phi'(b+z)=c\log z+A+iB+O(z)
\]
for a real nonzero logarithmic coefficient \(c\) and real
\(A,B\). The imaginary derivative on the real axis to the
right has limit \(B\), and on the left it has limit \(B+c\pi\).
Integrating gives, as \(z\to0\),
\[
 \Im(\phi(b+z)-\phi(b))
   =c y\log|z|+(B+c\arg z)x+O(y+|z|^2).                     \tag{3}
\]
The value \(\phi(b)\) is a finite limit. The imaginary
part is constant along the descent flow
\(\dot q=-\overline{\phi'(q)}\), so a trajectory landing at
\(b\) would have the left-hand side of (3) equal to zero.

At \(0\), \(c=-1\), the right boundary derivative is
\(\arg\zeta>0\), and the left derivative is
\(\arg\zeta-\pi<0\). At \(1\), \(c=1\) and these signs
are reversed. At a finite \(p_L\), \(c>0\), the right
boundary derivative is \(\arg\zeta-\pi\) and the left
one is at least \(\arg\zeta\). At a finite \(p_R\),
\(c<0\), the signs are reversed. In each case the
term \(c y\log|z|\) and the term involving \(x\) have
the same strict sign near either real ray. In angular
directions bounded away from those rays, the \(y\log|z|\)
term dominates. Therefore the imaginary difference in
(3) is nonzero throughout a sufficiently small upper
half-disk.

At any farther real pole \(p_i\), both boundary imaginary
derivatives are strictly positive. In this case (3) can
vanish near \(b\) only on the right side when \(c>0\),
or on the left side when \(c<0\). Indeed away from
the real directions the first term dominates, and on
the other side the two leading terms agree in sign.
But the real component of the descent vector has sign
\(c\) near \(b\), since it is \(-c\log|z|-A+O(|z|)\).
Thus it points away from \(b\) on the side where the
imaginary difference can vanish. This excludes landing
at the remaining poles.

It follows that all finite real exits of descending upper
trajectories lie in the two intervals \((p_L,0)\) and
\((1,p_R)\). The phase derivative analytically continues
across them to the lower phase derivative, by the opposite
logarithmic determinations, and the lower vertical descent
inequality remains available. The exact contour and this
boundary exclusion do not yet identify the contributing
upper saddle cycles or control their transitions.

## Actual two-saddle expansion on the equal-degree graph

**Theorem.** For fixed fully shared supports of equal total
degree with \(\max\rho\ge1\), both upper saddles contribute to
the Rice contour in a neighborhood of every interior point
of their equal-modulus graph. In particular the consecutive
ordinary moment polynomials have no common root at all
sufficiently large orders near each compact subarc of that graph.

Use the equal-modulus graph and the two upper saddles proved
in the [phase separation theorem](MULTIVARIATE_REGULAR_SADDLE_PHASE_SEPARATION.md#equal-total-degrees-at-and-above-the-real-cusp).
Here \(\Delta=0\), both signs of gaps occur, and \(p_L<0<1<p_R\)
are finite. Join the two open half-planes through the gaps
\((p_L,0)\) and \((1,p_R)\). The real parts of the two phases
agree across each gap. Their derivatives agree there because
the jumps of the logarithms of \(q,1-q\) cancel the derivative
of \(-2\pi i q\). The Stirling expansions continue across
the gaps. For the left gap use reflection,
\[
 \frac{\Gamma(n(q-1))}{\Gamma(nq+1)}
 =(-1)^{n+1}\frac{\Gamma(-nq)}{\Gamma(1+n(1-q))},
\]
whose two arguments are positive on the real gap. For the
right gap the original ratio already has positive arguments.
All \(x_i\) are positive on the real gaps.
[Gamma reflection formula](https://dlmf.nist.gov/5.5.E3)

Consider the two curves descending from each simple saddle for
\(\dot q=-\overline{\phi'(q)}\), using the lower phase after
crossing a gap. They start in the two downward directions of
a local analytic quadratic coordinate. The real phase decreases
strictly away from the saddle, with derivative \(-|\phi'(q)|^2\),
and its imaginary part is constant along the upper flow.
At an equal-modulus point no descending curve can approach
the other saddle, since their real phases agree. There are
no other stationary points in the domain.

Each descending curve escapes to infinity. At a non-gap
real interval away from a branch point the upper flow
points into the upper half-plane. Such a point cannot be a
forward accumulation point, since the locally extended
backward flow from that alleged limit reaches the opposite
side while all the approaching trajectories stay in the
upper half-plane. At a real branch point (3) rules out
accumulation. In the farther-pole case where (3) can vanish,
the phase level has \(y/|x|\to0\) near the allowed real ray.
A curve entering a sufficiently small half-disk
\(|z|\le\epsilon\) on this level has \(|x|\ge\epsilon/2\)
on entry, and its real coordinate moves away from the pole
as long as it stays inside. It therefore cannot accumulate
there. Once a curve crosses a gap into the lower half-plane,
its vertical velocity is bounded above by
\(\arg\zeta-\pi<0\), so it cannot return. A bounded curve
staying away from these boundaries would have a stationary
limit point: the real phase would be bounded below and
\(\int |\phi'|^2\,dt\) finite, giving a subsequence on which
\(\phi'\to0\). Strict real-phase loss excludes this, since
both saddles have the same real phase.

With \(\theta=\arg\zeta\), the derivatives of the two phases
tend uniformly at infinity to constants \(C_\pm\) in the
corresponding half-planes. This follows from the expansion
of \(X\) at infinity with \(\Delta=0\). They satisfy
\[
 \Im C_+=\theta+\pi(D_--1)>0,\qquad
 \Im C_-=\theta-\pi(D_-+1)<0,\qquad
 D_-=\sum_{d_i<0}|d_i|\ge1 .
\]
Thus when \(|\Im q|\) exceeds a fixed bound the vertical
velocity points outward, uniformly on a compact parameter
neighborhood. Every unbounded trajectory eventually
reaches one of these regions. If it goes arbitrarily far
horizontally with bounded imaginary part, the velocity
there has a nonzero outward vertical component and bounded
horizontal speed. Starting at a sufficiently large real
coordinate, it reaches a large height before returning
to a fixed bounded disk. The end type of a descending
curve is locally constant along the parameter graph.
Indeed start on the appropriate local quadratic direction,
use continuous dependence for the finite time needed to
reach an outward region, and then use its strict vertical
velocity bound.

At the negative-axis end \(\theta=\pi\), both stationary
points lie in the real gaps, with \(\phi''=X'>0\).
Each has one descending direction into the upper half-plane
and one into the lower. Above the real line \(\Im X>-\pi\),
so the upper curve moves strictly upwards everywhere.
The lower curve similarly moves strictly downwards.
They escape to infinity by the same boundary and stationary
point arguments. The upper phase has no stationary point
there at \(\theta=\pi\), since \(X\) cannot have height
\(-\pi\) in the open upper half-plane. In this geometry
each curve is homotopic to the corresponding vertical
contour by horizontal interpolation at the same imaginary
part. It crosses the real line in the same gap. The
Gamma representation holds at this endpoint by taking
\(\operatorname{Log}(-\zeta)\) real, and continuity into
the graph follows from the regular saddle in each gap.

Connectedness of the graph now keeps the end types and the
contour homotopies fixed. To see the homotopy for exact
integrands, truncate the curves at heights \(+Y,-Y\)
beyond the outward regions and join them horizontally to
the original vertical lines. Along any compact parameter
path starting at the negative-axis end one can choose a
common sufficiently large \(Y\), by finite time continuous
dependence and the strictly outward vertical velocities.
The truncated contours vary continuously in the domain
with the three real cuts removed. At infinity,
\(\phi_\pm(q)=C_\pm q+O(\log(1+|q|))\). The descending
ends have direction \(-\overline{C_\pm}+o(1)\), since their
velocity approaches that constant. At those ends and at
the fixed vertical lines, the real part of the linear
phase has strict linear loss in \(|Y|\), uniformly on a
compact parameter path. It has the same loss along the
horizontal connectors by linearity, up to the logarithmic
error. The angles of \(q\) there are separated from real
directions, since the velocities are asymptotically constant
with nonzero vertical components. The exact Stirling
bounds give decay of the connector integrals as \(Y\to\infty\).
The two actual vertical periods therefore equal the
corresponding descending integrals, with the orientations
fixed at the negative-axis end.

At each interior graph point the two upper stationary
values have imaginary difference strictly between zero
and \(\pi\). Thus in a small parameter neighborhood no
upper descending trajectory can connect the saddles.
The same finite time argument keeps both end types and
contour classes there. The local analytic quadratic
coordinate at each saddle and the relative expansion (2)
give a Gaussian amplitude, analytic and nonzero because
\(A(q_j)\ne0\) and the saddle is simple. The rest of the
contour has a uniform strict loss after a fixed local
saddle neighborhood, and at infinity it has the linear
loss just established. These give uniform relative
\(o(1)\) errors for both descending integrals on a smaller
parameter disk.

The consecutive comparison applies to these actual periods.
Because \(\Delta=0\), the external scale \(\zeta=t\) is
unchanged between orders. The ratio of the two stationary
exponential multipliers is \(e^{D(w)}\ne1\) at the graph.
By continuity this ratio stays separated from one on a
smaller disk, excluding common zeros there. A finite cover
gives the uniform exclusion near a compact subarc.
