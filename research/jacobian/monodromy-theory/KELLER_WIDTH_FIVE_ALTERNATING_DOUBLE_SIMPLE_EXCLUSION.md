# Exclusion of the double and simple critical pattern in the width-five elliptic slice

**Theorem SixSheetWidthFiveAlternatingDoubleSimpleExclusion.**
In the six-sheet \(A_6\) profile with one unramified omitted prime,
a width-five coordinate of slice genus one cannot have one double
and one simple critical point on the immersed locus of the sole
ramified branch. This includes coincidence of the double critical
value with the projection of the affine cusp.

Write \(p_0\) for the projection of the affine cusp and \(p_d,p_s\)
for the double and simple critical values on the immersed locus.
By the [width-five critical-pattern reduction](KELLER_WIDTH_FIVE_ALTERNATING_THREE_SIMPLE_EXCLUSION.md#the-remaining-critical-point-pattern),
\(p_d\ne p_s\), the cusp projection has local degree two, and the
pole count is four. Use the
[singular relative compactification](KELLER_SIX_SHEET_CUSP_SLICE_GENUS.md#singular-normalization-with-one-ordinary-affine-cusp).
Resolve the \(A_1\) point, with exceptional curve \(D_{\rm exc}\).
On the resolved finite-base surface the fibres are reduced and
connected. Each component other than \(D_{\rm exc}\) meets at
least one of the four pole sections at a smooth point of the fibre.

## The special fibres

At \(p_d\) the immersed index-three power map gives an ordinary
triple point of the fibre, from \(t^3+v^3=0\).
At \(p_s\) it gives an ordinary cusp, from \(t^2+v^3=0\).
At \(p_0\) the transverse cusp model gives two nodes on
\(D_{\rm exc}\). The two strict branches through that curve
meet it at distinct points, and no other branches meet it.

If \(p_d\ne p_0\), the fibre at \(p_d\) has three smooth
rational components at their common triple point. In fact
this singularity has delta invariant three, so the genus
formula requires at least three components, while
connectivity at the triple point requires at most three.
Away from \(p_0,p_d\) the fibres are integral. If the two
branches on \(D_{\rm exc}\) belong to the same component,
the fibre at \(p_0\) consists of two smooth rational curves
meeting twice, also by the genus formula. The dualizing
bundle is then trivial on every finite fibre. For the
two- and three-component fibres this follows by adjunction:
each component's self-intersection is \(-2\), so the
dualizing bundle has degree zero on it. A global nonzero
section exists by arithmetic genus one and is nonvanishing
on all components by connectivity. Integral arithmetic-genus-one
fibres have trivial dualizing bundle by duality and degree
zero. Evaluation and base change would then make the
relative dualizing bundle trivial, contradicting the
[free horizontal ramification class](KELLER_SMOOTH_NORMALIZATION_POLAR_BOUNDARY.md#canonical-class-on-the-finite-affine-normalization)
on the smooth locus of the original normalization. The
class calculation holds on this locus because its complement
in the normal surface has codimension two.

Thus the \(p_0\)-fibre in this case has two strict components,
one rational end \(A\) and one integral arithmetic-genus-one
component, linked by \(D_{\rm exc}\) at ordinary nodes.
Contract \(A\), a \((-1)\)-curve, then the image of
\(D_{\rm exc}\), now a \((-1)\)-curve. Any cusp on the
surviving component is preserved. The finite singular
fibres of the resulting relatively minimal model \(S_0\)
are of types \(IV\) and \(II\) at \(p_d,p_s\), respectively.

If \(p_d=p_0\), that fibre has the two nodes on \(D_{\rm exc}\)
and an ordinary triple point away from it. No other
singularity is possible, since \(p_s\) is different.
The total delta invariant is five. Arithmetic genus
one and at most five components therefore force five
smooth rational components. Connectivity requires the
triple point to join three different components and the
two nodes to link the fourth strict component \(A\) to
one of those three through \(D_{\rm exc}\). Every strict
component contains exactly one pole. Again \(A\) is a
\((-1)\)-curve, and contracting it and then the image
of \(D_{\rm exc}\) leaves type \(IV\) at \(p_0\).
The only other finite singular fibre is integral of
type \(II\) at \(p_s\).

In either case, the Picard group of the resolved finite-base
surface before contraction has rank seven by localization
from the affine plane: there are the two omitted primes,
the four poles and the exceptional curve. The two contractions
leave rank five. Units on \(S_0\) are constant since they
pull back to units on the original affine plane. A relatively
minimal completion is a rational elliptic surface with a
section. Localization from its rank-ten Picard group gives
five components at infinity. Its Euler number there is
six because types \(IV\) and \(II\) contribute six at the
finite values. By
[Miranda's Table 1.1](sources/miranda-persson-rational-elliptic-fibres-1990.pdf),
the infinity fibre has type \(I_0^*\). The \(j\)-function
has no poles and has value zero at the type-\(II\) fibre.
Thus, with any marked section as origin, the Weierstrass
equation is
\[
                         Y^2=X^3+D(p),
\]
where \(D\) is a cubic with a double zero at \(p_d\) and
a simple zero at \(p_s\).

## The poles and their distribution

The canonical bundle on \(S_0\) is trivial by fibrewise
triviality and evaluation, as above. The pullback of a
nowhere-zero regular two-form from this model to the
resolved surface has divisor \(2A+D_{\rm exc}\).
Consequently the
[marked factor argument](KELLER_TWO_POLE_ELLIPTIC_SLICE_MARKED_FACTORS.md#the-singular-alternating-profile-at-width-four)
on the original affine plane gives \(dq=\varphi^2\omega\)
on the generic fibre, with \(\varphi\) a nonzero constant
multiple of the reciprocal of the polynomial factor cutting
out \(A\cap X\), and \(\omega=dX/(2Y)\). Every pole of \(q\)
has odd order, so the four orders are \(3,1,1,1\).

This square equation gives odd ramification for \(q\) on
each strict component of the type-\(IV\) fibre, where the
dualizing differential on the normalization of a component
has a double pole at the triple point and otherwise is
nonvanishing. It gives the same parity on \(A\):
in the twice blown relative normal coordinate the pullback
of \(\omega\) vanishes to order two along \(A\), and
\(\varphi\) has pole order one there. If the strict
component opposite \(A\) through \(D_{\rm exc}\) is
integral of arithmetic genus one, its normalization
also has the odd ramification property. For a cuspidal
component its dualizing differential has a double pole
at the cusp and no other zeros or poles.

Suppose first that \(p_d\ne p_0\). At type \(IV\) the
three components each contain a pole, so one carries
two marks. That component cannot have two simple
poles: a degree-two rational map has some even local
degree by Riemann–Hurwitz. Thus the triple pole and
one simple pole lie on the same component. The other
two simple marks belong to its two companion components.
At the chain over \(p_0\), the only possible groups on
\(A\) are one simple pole, the triple pole alone, or
all three simple poles. A degree-two map on \(A\) or
on the opposite component is impossible by the odd
ramification property and Riemann–Hurwitz. A degree-one
map on the opposite integral genus-one component is
impossible because the finite birational map to the
normal projective line would be an isomorphism.

Suppose instead that \(p_d=p_0\). There is exactly one
mark on \(A\), and the triple pole lies on \(A\) or on
the strict component adjacent to \(D_{\rm exc}\).
Indeed the projection degree five is exhausted at this
value by local degrees two and three at the affine cusp
and the double critical point. There is no other finite
ramification of \(q\) on the strict components away from
their intersections with \(D_{\rm exc}\) and their common
triple point. At the triple point \(q\) is unramified on
each of the three strict branches: in the local cubic
power model \(q\) is the immersed tangent parameter, and
has nonzero derivative on each branch of \(t^3+v^3=0\).
A strict component not attached to \(D_{\rm exc}\) would
therefore have no finite ramification away from its pole.
A pole of order three on such a rational component would
contradict Riemann–Hurwitz. After the two contractions the
triple pole and one simple pole thus lie on the component
that was attached to \(D_{\rm exc}\), while the other two
simple poles lie on the other components of the type-\(IV\)
fibre.

## The marked-section equations

In both cases choose the triple pole section as origin
\(O\) and write the three simple sections as \(T_i=(u_i,v_i)\).
The two simple sections on the two nonidentity components
of the type-\(IV\) fibre are labeled \(T_2,T_3\); \(T_1\)
lies on the identity component at that value. None meets
\(O\) over the finite base, so \(u_i,v_i\in\mathbb C[p]\).
Near the type-\(IV\) value set \(k=p-p_d\), and write
\(D=k^2 d(k)\), \(d(0)\ne0\). At both nonidentity
sections \(u_i,v_i\) are divisible by \(k\). In the
blowup chart \(X=kU,Y=kV\) the equation is
\[
                              V^2=kU^3+d(k).
\]
The two components of its exceptional locus have
\(V=\pm\sqrt{d(0)}\). Each section on a nonidentity
component lies in this chart, since \(k\) is a local
parameter on the section and \(X,Y\) both vanish there.
Thus \(v_2/k,v_3/k\) have opposite nonzero residues,
and \(u_2,u_3\) vanish at \(p_d\). The section \(T_1\)
has nonzero abscissa at \(p_d\), since it meets the
identity component away from the origin and the
contracted Weierstrass cusp.

As in the [three-simple marked calculation](KELLER_WIDTH_FIVE_ALTERNATING_ELLIPTIC_MARKS.md#the-square-differential-at-three-simple-sections),
write \(\varphi=U_0X+\sum_i C_i f_i+L\), where
\(f_i=(Y+v_i)/(X-u_i)\), and \(\operatorname{Res}_{T_i}
(f_i\omega)=1\) on the generic fibre. The corresponding
leading polar coefficient of \(\varphi\) is unit on
each unblown section over the finite base. It has
order \(s\) at \(p_0\) on a twice blown section where
\(q\) has order \(s\). This follows from the same
\(\epsilon^2\) scaling of the relative normal parameter
at the blowups, for \(\epsilon=p-p_0\). The argument
applies at type \(IV\) as well, since every marked
section there meets a smooth point of its fibre, and
the invariant differential is nonzero in a relative
parameter at that point.

If \(A\) carries the triple mark alone, \(U_0\) has
order three at \(p_0\) and the \(C_i\) are nonzero
constants. Pole comparison in \(dq/\omega=\varphi^2\)
gives
\[
                 q=U_0^2Y/3 + U_0(\textstyle\sum_i C_i)X
                                     -\sum_i C_i^2 f_i+\text{constant}
\]
on the generic fibre. The order-two polar coefficient
at \(O\) must vanish to order at least four on the
minimal model, by the two blowups of the order-three
pole section. The displayed coefficient
\(U_0\sum_i C_i\) therefore gives \(\sum_i C_i=0\).
The simple-pole residue conditions then imply
\(\sum_i C_i u_i=0\): multiply each condition
\(L+U_0u_i+\sum_{j\ne i} C_j f_j(T_i)=0\) by \(C_i\)
and add, using the skew evaluations. At \(p_d\) the
terms for \(T_2,T_3\) vanish but \(C_1u_1(p_d)\ne0\).
This contradiction excludes the triple mark on \(A\)
in both configurations.

If \(A\) carries all three simple marks, necessarily
\(p_0\ne p_d\), and these sections have pairwise contact
exactly two at \(p_0\). In this model an integral
section cannot have abscissa of degree four: the
[integral-section calculation](KELLER_WIDTH_FIVE_ALTERNATING_ELLIPTIC_MARKS.md#integral-sections-of-the-i_0-model)
forces a centered cubic \(D=a t^3+c\), which has no
double root unless it has a triple root. Hence all
three abscissas have degree at most two. The
[shared-jet comparison](KELLER_WIDTH_FIVE_ALTERNATING_ELLIPTIC_MARKS.md#all-three-simple-marks-on-the-rational-end)
shows \(u_i=\alpha_i\epsilon^2\) with distinct
\(\alpha_i\) and a common ordinate \(v\), with
\(v(p_0)\ne0\). In fact comparing degrees six, five
and four in \(v^2-u_i^3=D\) gives
\(v=v(p_0)+\gamma\epsilon^3\), where
\(\gamma^2=\alpha_i^3\ne0\). Thus \(D\) is a
constant plus a nonzero cube term in \(\epsilon\),
and is squarefree, contrary to its double root.

Finally let \(A\) carry one simple mark. All sections
are then pairwise disjoint over the finite base,
the origin is unblown, and every evaluation \(f_j(T_i)\)
for \(i\ne j\) is polynomial in \(p\). Indeed this is
the usual regular evaluation away from type \(IV\).
At type \(IV\), an identity and a nonidentity mark
have distinct abscissas. For the two opposite
nonidentity marks the identity
\[
                    f_3(T_2)=
            \frac{u_2^2+u_2u_3+u_3^2}{v_2-v_3}
\]
shows regularity: the numerator vanishes to order at
least two and the denominator exactly to order one.
This identity also holds for opposite generic
sections with identical abscissas, by the regular
evaluation formula.

After scaling the invariant equation and reindexing
the simple marks so that the blown one has index
three, the residue and origin equations are exactly
(1)–(3) of the
[three-simple exclusion](KELLER_WIDTH_FIVE_ALTERNATING_THREE_SIMPLE_EXCLUSION.md#equations-for-the-marked-sections).
Its algebraic contradiction uses only the cubic
degree of \(D\), polynomial abscissas and ordinates,
regular finite evaluations, disjoint sections and the
stated leading polar coefficients. It does not require
that the roots of \(D\) be simple. Here degree-four
abscissas have already been excluded by the
integral-section calculation. Thus that contradiction
applies, ruling out the last group of marks on \(A\).

## Consequences for width and genus

**Corollary SixSheetSingularAlternatingEllipticWidthAtLeastSix.**
In the six-sheet \(A_6\) profile with one unramified omitted
prime, every coordinate with a genus-one generic slice has
directional width at least six. Any width-five coordinate in
this profile has slice genus and pole count \((2,2)\).

The [width bound](KELLER_ALTERNATING_ELLIPTIC_TRIPLE_POLES.md#directional-and-total-degree-bounds)
gives width at least five in every coordinate and leaves
\((1,4)\) and \((2,2)\) as the width-five possibilities.
The critical-pattern reduction and the theorem exclude the
first possibility.
