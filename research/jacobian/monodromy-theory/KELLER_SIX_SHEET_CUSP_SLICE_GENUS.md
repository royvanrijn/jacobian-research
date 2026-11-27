# Slice genus in the six-sheet cusp profiles

Let \(F\) be a hypothetical complex plane Keller map of generic sheet
degree six, with full finite normalization \(\pi:W\to\mathbb A^2_\mathbb C\).
For a polynomial target coordinate frame \(p,q\), let \(g_p\) and \(r_p\)
be the compactified generic \(p\)-slice genus and number of poles.
Write \(w(p)\) for the directional width of the sole branch in the
alternating case.

**Theorem SixSheetAlternatingCuspProjectionGenus.** Let the monodromy
be \(A_6\). If there are no unramified omitted primes, then
\(g_p\ge2\) and \(w(p)\ge5\) for every \(p\). If there is one unramified
omitted prime, then \(g_p\ge1\). In the latter case width four forces
\(g_p=1\), \(r_p=2\), and the projection of the ordinary affine cusp
has local degree two on the branch normalization.

The case with two unramified omitted primes is the immersed case of
[SixSheetImmersedAlternatingCompositeWidth](KELLER_IMMERSED_SOLE_BRANCH_PROJECTION.md#six-sheet-consequence).

## Fibre genus and canonical bundle

Use the [relative compactification over the \(p\)-line](KELLER_SMOOTH_NORMALIZATION_POLAR_BOUNDARY.md#compactification-over-the-p-line)
by normalizing \(\mathbb A^1_p\times\mathbb P^1_q\). Its polar sections
are smooth and disjoint, including when \(W\) has the isolated \(A_1\)
singularity of the one-unramified-prime case. If necessary, resolve
that point by the minimal resolution to obtain a smooth relative
surface. The new exceptional curve lies in the fibre over its image
under \(p\). The fibres are connected. Away from the exceptional
curve each fibre component has multiplicity one by the original
source submersion, and each such component meets a pole. Reducedness
along the exceptional curve is checked below.

For a connected, reduced fibre with \(k\) irreducible components, let
\(\delta_x\) and \(b_x\) be the delta invariant and number of analytic
branches of a singular point \(x\). If the fibre has arithmetic
genus \(g_p\), the normalization sequence and connectedness give
\[
 g_p=\sum_{i=1}^k g_i+\sum_x\delta_x-k+1,\qquad
 \sum_x(\delta_x-b_x+1)\le g_p ,
\]
where \(g_i\) are the normalization genera. The inequality uses
\(\sum_x(b_x-1)\ge k-1\).

At a smooth immersed index-three germ, the only possible singularities
of the fibre are \(t^m+v^3=0\), with \(m\ge2\) the local degree of \(p\)
on the image branch. As in the
[immersed genus argument](KELLER_IMMERSED_SOLE_BRANCH_PROJECTION.md#genus-of-slices-with-immersed-index-three-ramification),
these each contribute at least one to \(\sum_x(\delta_x-b_x+1)\).
The contribution is one exactly for an ordinary cusp or an ordinary
triple point. An arithmetic-genus-one fibre whose singularities are
of these two types is integral or is the union of three smooth
rational curves at their common ordinary triple point. In the latter
case each component has self-intersection \(-2\), by intersecting it
with the full fibre. Both types have trivial dualizing bundle, by
adjunction and the restriction of a nonzero dualizing section to
the components, as in that argument. The same reasoning works for
two smooth rational components meeting transversely at two distinct
points, since each component again has self-intersection \(-2\).

If \(g_p=1\) and every finite fibre has trivial dualizing bundle,
relative duality and base change make the canonical bundle on the
whole relative smooth surface a pullback from \(\mathbb A^1_p\),
hence trivial. Its restriction to the smooth locus of \(W\) is
nontrivial. Indeed the Weil class group of \(W\) is free on the
omitted primes by localization and the absence of nonconstant
units on the affine plane. The canonical Weil divisor from
\(\pi^*(dp\wedge dq)\) is twice the ramified prime, and restriction
across a codimension-two set induces an isomorphism of class groups.
Thus the trivial-dualizing alternative is impossible.

## Smooth normalization with one \((2,5)\) affine cusp

By the [six-sheet cusp classification](KELLER_BOUNDARY_CUSP_PROFILES_AND_OREVKOV_INVARIANT.md#consequences-in-generic-sheet-degree-six),
when there are no unramified omitted primes, \(W\) is smooth and the
only nonimmersed branch point has local finite degree five. Use the
rank-one local model from that proof, with \(\pi=(a,f(a,v))\),
omitted prime \(a=v^2\), and
\[
 f(a,0)=0,\qquad
 f_v(a,v)=(a-v^2)^2\epsilon(a,v),\qquad \epsilon(0,0)\ne0 .
\]
The model with \(f(a,0)=0\) uses a local analytic target change, so
\(p\) is evaluated as an analytic target coordinate with nonzero
differential in this model. Its fibre at the cusp is smooth if
\(dp\) does not annihilate the rank-one image of \(d\pi\).
Then \(p\) has local degree two on the branch normalization.
If \(dp\) annihilates that rank-one image, the implicit function
theorem writes the fibre equation as \(f(a,v)-h(a)=0\), with
\(\operatorname{ord}(h)\ge2\).

Put \(\lambda=\epsilon(0,0)\). If \(h\) has order two and leading
coefficient \(c\), the \((5,2)\)-weighted initial part in \(a,v\)
is \(-c a^2+\lambda v^5/5\). There is one plane branch of
characteristic orders \(2,5\), with delta invariant two. In detail,
the substitution \(v=s^2,\ a=s^5 A\) yields two distinct simple
roots for \(A\) at \(s=0\), exchanged by changing the sign of \(s\).
They give one branch of degree two over \(v\), exhausting the local
degree from Weierstrass preparation in \(a\).

If \(\operatorname{ord}(h)\ge3\), the substitution \(a=v^2 A\), after
division by \(v^5\), gives at \(v=0\) the quadratic
\(\lambda(A^2-2A/3+1/5)\). Its two simple roots give smooth branches
\(a=v^2 A_i(v)\) of mutual intersection multiplicity two. A third
smooth branch follows by setting \(v=a V\) and dividing by \(a^3\).
At \(a=0\) the equation is \(\lambda V-h_3=0\), with \(h_3\) the
cubic coefficient of \(h\). The third branch meets the first two
transversely. These exhaust the multiplicity three of the original
equation. Thus \(\delta_x=4\) and \(b_x=3\).

In both tangent cases \(\delta_x-b_x+1=2\), forcing \(g_p\ge2\).
In the transverse case the cusp contributes no singularity to the
source fibre. Width is already at least four by
[SixSheetProjectionWidths](KELLER_DIRECTIONAL_RAMIFICATION_AND_POLAR_ENDS.md#six-sheet-consequences).
The projection has local degree two at the cusp, so at least one
other finite critical point lies on the immersed locus. Its fibre
germ excludes \(g_p=0\) by the normalization genus inequality.
If \(g_p=1\), every singularity in a finite fibre is an ordinary
cusp or an ordinary triple point. The trivial-dualizing contradiction
above excludes this genus too. Now \(2w(p)=6+2g_p+r_p-2\)
and \(r_p\ge1\) give \(w(p)\ge5\).

At width five, the same formula gives \((g_p,r_p)=(2,2)\).
The local degree of \(p\) on the cusp normalization can be two,
four or five. These are, respectively, the transverse case, the
tangent case with \(\operatorname{ord}(h)=2\), and the tangent case
with \(\operatorname{ord}(h)\ge3\), since the branch orders in
\(a,f\) are two and five. Local degree five is impossible at width
five. In that case all five punctures in a nearby generic target
line coalesce at the cusp, so their meridians, with a common stem,
belong to the local complement group of the cusp. The image of that
group preserves the sheet block of the degree-five germ and fixes the
sixth sheet. This contradicts the [punctured-line generation theorem](KELLER_DIRECTIONAL_RAMIFICATION_AND_POLAR_ENDS.md#meridians-and-polar-ends).

## Singular normalization with one ordinary affine cusp

Here the
[local degree-four classification](KELLER_BOUNDARY_CUSP_PROFILES_AND_OREVKOV_INVARIANT.md#small-ramification-indices-at-a-contracted-endpoint)
gives \(W\) one \(A_1\) point over the sole ordinary cusp. The cusp has
local projection degree at most three in every polynomial target
coordinate, since its orders are two and three. Width at least four
therefore supplies an immersed critical branch point. Its source fibre
germ is unchanged on resolving the \(A_1\) singularity and cannot
lie in a simple normal crossings rational fibre. This excludes
\(g_p=0\) by the [rational-fibration lemma](sources/gurjar-miyanishi-contractible-curves-1996.pdf),
Lemma 7(i). Width four then gives \((g_p,r_p)=(1,2)\).

Use the explicit analytic degree-four model from the local
classification on the double cover of the \(A_1\) point:
\[
 (u,v)=(x^2+4xy^3+\tfrac83 y^6,\ xy+y^4)
\]
with the involution \((x,y)\mapsto(-x,-y)\). In the resolution chart
\(x=sy,\ z=y^2\), the pullbacks are
\[
 u=z(s^2+4sz+\tfrac83 z^2),\qquad v=z(s+z).
\]
The exceptional curve is \(z=0\), and a target coordinate vanishing at
the cusp has linear part \(A u+Bv\). Its order along the exceptional
curve is one, since its pullback divided by \(z\) restricts in this
chart to \(As^2+Bs\) and this polynomial is nonzero. The resolved
fibres are therefore reduced.

The projection has local degree three on the image cusp if and only
if \(B=0\). Then \(A\ne0\), and after dividing \(p\) by \(z\), the
leading homogeneous part at \(s=z=0\) is
\(A(s^2+4sz+\tfrac83z^2)\). All contributions from the higher terms
in the target coordinates have order at least three there. The
quadratic has distinct factors, neither proportional to \(z\).
Thus the resolved fibre has an ordinary triple point, consisting
locally of the exceptional curve and two smooth transverse strict
branches. If \(g_p=1\), the triple point uses the full genus allowance
and its fibre has trivial dualizing bundle. All other finite fibres
also have trivial dualizing bundle by the immersed argument. This
contradicts the canonical class calculation, so width four can
occur only in a transverse coordinate, with \(B\ne0\).

**Corollary SixSheetSingularAlternatingWidthFourFibre.** In the
one-unramified-prime \(A_6\) case at width four, the resolved fibre
over the affine cusp value consists of the exceptional \((-2)\)-curve
between a smooth rational \((-1)\)-curve and an integral
arithmetic-genus-one component. Contracting the rational end and
then the exceptional curve leaves a relatively minimal
genus-one fibration with precisely two finite cuspidal fibres.

For \(B\ne0\), the two distinct tangent lines of the lifted \(p=0\)
germ are the roots of \(x(Ax+By)\). On the resolution their strict
transforms are transverse to the exceptional curve and meet it at
distinct points. These are ordinary nodes of the resolved fibre.
There are no other fibre components through the exceptional curve.
The \(p\)-projection of the branch normalization has one simple
critical point at the affine cusp. The other critical points must
be simple, since a triple point away from the exceptional curve
gives total delta at least three and cannot occur in a genus-one
fibre with at most two components. A critical point of local degree
three on the immersed locus cannot lie over the cusp's \(p\)-value,
because the two local projection degrees would sum to more than
four. Thus there are exactly two simple critical points on the
immersed locus. They give distinct \(p\)-values by the genus inequality
and can give at most one cusp in the resolved fibre through the
exceptional curve.

If the two strict branches there belong to one global component,
it and the exceptional curve meet transversely at two points.
The genus formula allows no other singularity and makes both
components smooth rational. This gives trivial fibre dualizing
bundle. The other fibres are integral or of the immersed triple
type, also with trivial dualizing bundle. Hence the branches belong
to distinct global components. Each meets the exceptional curve
once, and the normalization genus formula forces one smooth rational
component and one integral arithmetic-genus-one component, smooth or
cuspidal. They cannot meet at any other point: away from the
exceptional curve all additional singularities are ordinary
unibranch cusps. The smooth rational end has self-intersection
\(-1\) by intersecting it with the full fibre. Its contraction
makes the exceptional curve a \((-1)\)-curve. Contracting this
next curve leaves the integral genus-one component. These
contractions are at smooth points of the genus-one component in
the reverse direction and preserve its cusp if present. All other
finite fibres are integral. The two cusps from the two simple
immersed critical points are exactly the finite singular fibres of
the resulting relatively minimal fibration.

## Mixed symmetric profile at directional cost eight

**Proposition SixSheetMixedEqualityEllipticFibre.** In the mixed
six-sheet \(S_6\) profile, suppose that \(w_T(p)+2w_C(p)=8\).
Then \((g_p,r_p)=(1,2)\) and \((w_T(p),w_C(p))\) is either \((2,3)\)
or \((4,2)\). In the \((4,2)\) case the transposition projection
has at least one simple critical point. For \((2,3)\), and also
for \((4,2)\) with one simple and one double transposition critical
point, the finite-base compactification has exactly one reducible
fibre. It consists of a smooth rational \((-1)\)-curve and an
integral arithmetic-genus-one component, joined at an ordinary
node. Contracting the rational component gives a relatively
minimal genus-one fibration with precisely two finite cuspidal fibres.

The [smooth-normalization bound](KELLER_SMOOTH_NORMALIZATION_POLAR_BOUNDARY.md#six-sheet-consequences)
excludes genus zero in this mixed profile, and the directional
formula gives the displayed genus and pole count. Each width is at
least two by the no-line result. The saturated mixed profile has
only immersed power-map germs, so at an index-two critical point of
local projection degree \(m\) the fibre is \(t^m+v^2=0\), by the
same local calculation as at index three. At most two components
occur in any fibre, giving total delta invariant at most two.

At width three the index-three projection has two simple critical
points, since local degree three gives an ordinary triple point of
delta invariant three. Their values are distinct by the arithmetic
genus inequality. At width two there is one simple critical point.
The index-three singularities of the fibres are therefore ordinary
cusps. The index-two projection at width two gives one ordinary
node. At width four it gives three ordinary nodes, an ordinary
cusp and an ordinary node, or one ordinary tangency of two smooth
branches, according to its critical multiplicities.

In the tangency case the only way a fibre is reducible is for two
components to meet there. Those two components are smooth rational
curves with intersection multiplicity two and self-intersection
\(-2\), by the genus formula and intersection with the full fibre.
No cusp fits on the tangency fibre because total delta is at most two.
The dualizing bundle on the tangency fibre is trivial, using degree
zero on each component and a global nonzero section as for a triple
fibre. All other fibres are integral, so their dualizing bundles are
also trivial. The canonical class contradiction rules out the
tangency case. This proves existence of a simple critical point.

In the cases with just one ordinary node in all finite fibres,
the other singularities are exactly two ordinary cusps, with
distinct values by the genus inequality. Only the node can connect
distinct components. If it joins one component to itself, every
finite fibre is integral and the canonical class contradiction
applies again. Thus the unique reducible fibre has two components
joined at this node. The genus formula gives a smooth rational
component and an integral arithmetic-genus-one component, possibly
cuspidal. The rational component has self-intersection \(-1\).
Contracting it preserves any cusp on the other component, since
it meets that component at a smooth point. All other fibres
remain integral, giving the stated relatively minimal model.

## The two Weierstrass models of the elliptic residuals

**Corollary SixSheetLowWidthWeierstrassResidual.** For the singular
\(A_6\) width-four case and the mixed \(S_6\) cases with the unique
reducible fibre just proved, the relatively minimal genus-one
fibration has a section. After an affine change in \(p\) and
constant Weierstrass rescaling, its projective relatively minimal
completion has one of the two Weierstrass models
\[
 \begin{split}
 Y^2 &= X^3+(p^2-1),\\
 Y^2 &= X^3-3(p^2-1)X+2p(p^2-1).
 \end{split}
\]
They have two finite type-\(II\) fibres and an infinity fibre of
type \(IV^*\) and \(I_2^*\), respectively. The claim concerns the
contracted elliptic fibration. The omitted primes, both pole
sections and the undoing of the finite contractions still have
to satisfy the original polynomial-source conditions.

The [width-five lower bound](KELLER_ALTERNATING_ELLIPTIC_TRIPLE_POLES.md#directional-and-total-degree-bounds)
excludes the singular alternating width-four residual.
The mixed symmetric cases continue to require marked admission.

Write \(S_0\) for the finite-base model obtained by the indicated
contractions. The original pole sections descend to sections. A
smooth projective completion over \(\mathbb P^1_p\) can be made
relatively minimal by contracting fibre \((-1)\)-curves only at
infinity, preserving \(S_0\). Denote it by \(\bar S_0\). It is rational
because its function field is the original plane function field.
The section extends across infinity.

In the mixed case the Picard group before contraction has rank
four, freely generated by the two omitted primes and two pole
sections by divisor localization from \(\mathbb A^2\). One
contraction leaves rank three. In the singular alternating case
the corresponding group on the resolution has rank five, because
the exceptional curve gives one more boundary component. Two
contractions also leave rank three. The group of units on \(S_0\)
is constant in both cases, because a unit pulls back and restricts
to a unit on the original affine plane.

By
[Miranda, Section 1, equation (1.2) and Table (1.1)](sources/miranda-persson-rational-elliptic-fibres-1990.pdf)
([author copy](https://www.math.colostate.edu/~miranda/preprints/Perssonslist.pdf)),
the Picard group of \(\bar S_0\) has rank ten and its singular
fibres have total Euler number twelve. Divisor localization and
the constancy of units on \(S_0\) show that the fibre at infinity
has \(10-3=7\) irreducible components. Its Euler number is eight,
since the two finite type-\(II\) fibres contribute two each.
The table gives \(IV^*\) or \(I_2^*\) for that pair of numbers.

Normalize the finite singular values to \(p=\pm1\). In Miranda's
Weierstrass form \(Y^2=X^3+A(p)X+B(p)\), the coefficient degrees are
at most four and six. Put \(D=p^2-1\). Both coefficients vanish at
the two cusps. The finite Weierstrass surface is the smooth model
there, since those fibres are integral. Smoothness of its total
space at the cusp gives a simple zero of \(B\). Hence the
discriminant has a double zero at each finite singular value and
no others on the affine line. Write \(A=D\alpha\), \(B=D\beta\).
Then
\[
 4D\alpha^3+27\beta^2=c,\qquad c\ne0 .
\]
If the fibre at infinity is \(IV^*\), the \(j\)-function has no pole
and equals zero, by the type-\(II\) finite values and the table.
This gives \(\alpha=0\) and nonzero constant \(\beta\). For \(I_2^*\)
the \(j\)-function has a pole of order two at infinity by the same
table. The formula for \(j\) makes it a nonzero constant multiple
of \(D\alpha^3\), so \(\alpha\) is a nonzero constant. The displayed
identity then gives \(\beta=bp\) with \(27b^2=-4\alpha^3\).
Constant Weierstrass rescaling gives the two stated forms over
\(\mathbb C\).
