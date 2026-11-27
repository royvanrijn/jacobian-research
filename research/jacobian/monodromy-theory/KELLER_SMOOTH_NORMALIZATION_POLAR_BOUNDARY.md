# Polar boundary of a smooth finite plane Keller normalization

Let \(F:X=\mathbb A^2_{\mathbb C}\to\mathbb A^2_{\mathbb C}\) be a
polynomial map with nonzero constant Jacobian and generic sheet degree
\(N>1\). Let \(\pi:W\to\mathbb A^2_{\mathbb C}\) be the full finite
normalization in the source field. The open immersion \(X\hookrightarrow W\)
over the target is understood. Fix a polynomial target coordinate system
\(p,q\). Write \(g_p,r_p,R_p\) for the generic compactified slice genus,
number of poles and weighted directional ramification cost of
[KellerDirectionalRamification](KELLER_DIRECTIONAL_RAMIFICATION_AND_POLAR_ENDS.md).

**Theorem KellerSmoothNormalizationPolarBoundary.** If \(W\) is smooth,
then \((g_p,r_p)\) is neither \((0,2)\) nor \((1,1)\). In particular
\(R_p\ge N+1\), and equality forces \((g_p,r_p)=(0,3)\). If in addition
some ramified omitted prime has normal index at least three and its
analytic map germ at every point of the prime is a cyclic power map
along an immersed image branch, then \(g_p\ge1\) and \(R_p\ge N+2\).
These conclusions hold in every polynomial target coordinate.

## Compactification over the \(p\)-line

Normalize \(\mathbb A^1_p\times\mathbb P^1_q\) in the source field to
obtain \(S\). The normalization is finite. Its projection \(f:S\to
\mathbb A^1_p\) is projective and its open inverse image over the affine
target is \(W\). Set \(H=(S\setminus W)_{\mathrm{red}}\).

The nonproperness curve \(A_F\) is
[finite over the \(p\)-line](KELLER_DIRECTIONAL_RAMIFICATION_AND_POLAR_ENDS.md#meridians-and-polar-ends).
Near any point \((p_0,\infty)\) on the partially compactified target,
choose a product of small \(p\)- and \(1/q\)-disks whose intersection with
the affine target avoids \(A_F\). The cover of the punctured product is
finite étale. Each connected component is the cyclic cover of a
punctured disk, pulled back along the second projection, since the
fundamental group is \(\mathbb Z\). Its normal extension to the product
is the smooth power map \((p,u)\mapsto(p,u^j)\) in \(p,1/q\) coordinates,
by uniqueness of finite normal extension. Thus \(S\) is smooth along
\(H\), and \(H\to\mathbb A^1_p\) is finite étale. Since the complex affine
line is simply connected, \(H\) consists of \(r_p\) disjoint sections.
They meet each fibre at smooth points.

Stein factorization and generic connectedness give connected fibres of
\(f\). Every irreducible component of a fibre maps onto its target
\(\mathbb P^1_q\) because the map is finite and its fibres have dimension
zero. It therefore contains a pole on \(H\) and also meets \(W\).
Each omitted prime of \(W\) is horizontal over \(p\), since its image
is a component of \(A_F\) and \(A_F\) is finite over \(p\). So every
fibre component meets \(X\). There \(p\circ F\) is a submersion, so each
component has multiplicity one. The fibres are Cartier divisors in the
smooth surface \(S\) and hence reduced. Each of the \(r_p\) marked points
belongs to exactly one component, giving at most \(r_p\) components
for every fibre.

## Canonical class on the finite affine normalization

The class group of \(W\) is free on the classes of the omitted primes
\(E\). Indeed, divisor localization gives generation because
\(\operatorname{Cl}(X)=0\). A relation would give a rational function
on \(W\) whose divisor restricts to zero on \(X=\mathbb A^2\). This
function and its inverse are regular on \(X\) by normality. Its being a
unit forces it to be constant, so its divisor on \(W\) is zero.

On smooth \(W\) the pullback of \(dp\wedge dq\) gives the canonical class
\[
 [K_W]=\sum_E(e_E-1)[E].
\]
Here \(e_E\) is the normal ramification index. Ramification can occur
only outside \(X\). At least one \(e_E>1\), since otherwise purity and
simple connectedness of the affine target would force the finite map
to have degree one. Freeness therefore gives \([K_W]\ne 0\).

## Low genus obstructions

Suppose \(g_p=0\). A smooth projective completion of \(S\) preserving
\(f\) over \(\mathbb A^1_p\) is obtained by taking its closure in a
projective space over \(\mathbb P^1_p\) and resolving over infinity.
In its resulting rational fibration the support of each fibre is a
simple normal crossings tree of projective lines.
This is [Gurjar–Miyanishi, Lemma 7(i)](sources/gurjar-miyanishi-contractible-curves-1996.pdf).
Over the finite base the fibres are already reduced.

If \(r_p=2\), write \(H=H_1+H_2\). Every reducible finite fibre has two
projective lines meeting transversally once, with one marked point
on each line. Its components have self-intersection \(-1\), by
intersecting each with the full fibre. Contract one component in each
such fibre. The contractions can be made in the smooth projective
completion and then restricted over the finite base. They give
\(\rho:S\to S_0\) over \(\mathbb A^1_p\), with smooth projective-line
fibres and disjoint image sections \(H_{1,0},H_{2,0}\). Each contraction
collapses a curve meeting exactly one section transversally once.
Thus each reverse blowup has its center on exactly one section.
The canonical divisor and strict-transform formulas give
\[
 \mathcal O_S(K_S+H_1+H_2)
 \simeq \rho^*\mathcal O_{S_0}(K_{S_0}+H_{1,0}+H_{2,0}).
\]
The bundle on \(S_0\) has degree zero on every projective-line fibre,
so it restricts to the trivial bundle on every fibre. Cohomology and
base change and the evaluation map give a pullback from
\(\mathbb A^1_p\). This pullback is trivial because
\(\operatorname{Pic}(\mathbb A^1)=0\). Restricting to \(W=S\setminus H\)
then gives \([K_W]=0\), a contradiction.

Suppose \(g_p=1\) and all finite fibres are integral. The map \(f\)
is flat since \(S\) is integral and the base is a smooth curve. Its
fibres are Gorenstein because they are Cartier in smooth \(S\), and
flatness gives arithmetic genus one on every fibre. On an integral
Gorenstein curve of arithmetic genus one, including a singular
curve, the dualizing bundle has one independent section by duality
and degree zero by Riemann–Roch. A nonzero section of a line bundle
on an integral curve has an effective zero divisor of degree equal
to the degree of the bundle. The section therefore trivializes the
bundle. Consequently, cohomology and base change and the evaluation
map give
\(\omega_{S/\mathbb A^1_p}\simeq f^* f_*\omega_{S/\mathbb A^1_p}\)
with \(f_*\omega_{S/\mathbb A^1_p}\) a line bundle. The base has
trivial canonical bundle and trivial Picard group. Again
\([K_W]=0\), a contradiction. In particular \(r_p=1\) is impossible
when \(g_p=1\), since the single marked point makes every fibre integral.

Finally let an omitted prime have the asserted power-map germ at
every point and normal index \(e\ge3\). Its image branch has
normalization \(\mathbb A^1\), and the restriction of \(p\) to that
normalization has degree at least two. Degree zero contradicts
finiteness over \(p\) and degree one would make the branch a polynomial
graph, contrary to Chau's exclusion of an affine-line component.
At a finite critical point the restriction has order \(m\ge2\) after
subtracting its critical value. Choose a target normal coordinate
\(z\) for its immersed analytic branch. For the power-map germ
\((t,v)\mapsto(t,v^e)\) the fibre equation is
\(t^m+v^e b(t,v^e)=0\) with \(b(0,0)\ne0\), after normalizing the
analytic parameter along the branch. The unit property follows from
\(dp\ne0\) and its vanishing along the tangent at the critical point.
Taking an analytic \(e\)-th root of \(b\) gives the fibre germ
\(t^m+v^e=0\). For \(m\ge3\) the multiplicity is at least three.
For \(m=2\) the tangent cone is a double line. Neither can occur in
the simple normal crossings fibre of a smooth rational fibration.
This excludes \(g_p=0\). The numerical conclusions now follow from
\(R_p=N+2g_p+r_p-2\), \(R_p\ge N\), and \(r_p\ge1\).

## Six-sheet consequences

**Corollary SixSheetSmoothProjectionBounds.** Every hypothetical
six-sheet plane Keller map with group \(S_6\) satisfies the following
bounds in every polynomial target coordinate \(p\).

| Ramified meridians | Directional bound | Weighted sum of total branch degrees in any frame |
|---|---|---:|
| One transposition | \(w_1(p)\ge7\) | \(\ge8\) |
| Two transpositions | \(w_1(p)+w_2(p)\ge7\) | \(\ge8\) |
| Transposition and 3-cycle | \(w_T(p)+2w_C(p)\ge8\) | \(\ge9\) |

By [SixSheetFiniteCuspProfiles](KELLER_BOUNDARY_CUSP_PROFILES_AND_OREVKOV_INVARIANT.md#consequences-in-generic-sheet-degree-six),
the \(S_6\) normalization is smooth and its mixed index-three prime
has zero excess, giving precisely the everywhere-immersed cyclic
power-map hypothesis. This hypothesis concerns a single finite
source germ. If another image branch passes through the critical
value, its sheets outside the cyclic germ do not change that germ
or its fibre equation. Apply the theorem with \(N=6\).
For the total-degree bounds in a fixed frame, choose a linear
combination of its two coordinates canceling the leading power of
the normalization parameter on one branch. Complete this to a
linear coordinate change. The resulting width on that branch is
strictly below its prior total degree, and the widths for the other
branches are at most their prior total degrees. All weights are
positive integers.

**Corollary SixSheetSmoothExceptionalFoldDegree.** In the exceptional-fold setting of
[SixSheetExceptionalFoldDegree](KELLER_DIRECTIONAL_RAMIFICATION_AND_POLAR_ENDS.md#six-sheet-consequences),
the sole-transposition case requires \(\max(\deg_a P,\deg_a Q)\ge48\).
In the original frame at least one coordinate restricts to degree at
least eight on the branch. Apply the weighted fold bound in that
original coordinate, using the coordinatewise estimate that
\(\deg A=j\) forces \(\deg_a P\ge6j\), or the corresponding estimate
for \(Q\).

**Corollary SixSheetImmersedAlternatingWidth.** In the six-sheet
\(A_6\) case with two unramified omitted primes, the immersed sole
3-cycle branch has width at least five in every polynomial target
coordinate.

The [affine cusp classification](KELLER_BOUNDARY_CUSP_PROFILES_AND_OREVKOV_INVARIANT.md#consequences-in-generic-sheet-degree-six)
gives smooth \(W\) and an everywhere-immersed index-three power map.
By the theorem \(g_p\ge1\), and \(R_p=2w_1(p)\ge8\).
Width four would force \(R_p=8\) and \((g_p,r_p)=(1,2)\).
Each finite fibre is then a connected reduced curve of arithmetic
genus one with at most two components. Its singular points are
among the critical points of \(p\) on the immersed branch. Over
an unramified sheet, or at a point where the \(p\)-projection of the
branch is regular, the fibre is smooth by the local power map.
At the remaining points the germ is \(t^m+v^3=0\) for \(m\ge2\).
For \(m\ge3\) this germ has delta invariant at least three, since a
reduced plane curve germ of multiplicity three has
\(\delta\ge\binom{3}{2}\). The normalization genus formula bounds
the total delta invariant of a connected arithmetic-genus-one
fibre with at most two components by two. So every singularity is
an ordinary unibranch cusp, with \(m=2\). Such points cannot join
distinct components, and connectedness makes every fibre integral.
The preceding dualizing-sheaf argument for integral genus-one
fibres again forces \([K_W]=0\), excluding width four.

The [immersed-branch argument](KELLER_IMMERSED_SOLE_BRANCH_PROJECTION.md)
strengthens this bound to composite width at least six and genus at least two.
