# Exclusion of two triple poles in the singular alternating width-four slice

**Theorem SixSheetSingularAlternatingTriplePoleExclusion.**
In the singular six-sheet \(A_6\) profile with one unramified omitted
prime and directional width four, the two poles of the complementary
coordinate on the generic elliptic slice cannot both have order three.

Let \(\epsilon=p-p_0\), where \(p_0\) is the projection of the affine
branch cusp. The
[resolved fibre and its contractions](KELLER_SIX_SHEET_CUSP_SLICE_GENUS.md#singular-normalization-with-one-ordinary-affine-cusp)
give a relatively minimal elliptic fibration after contracting the
rational end and then the exceptional curve. One pole section met the
rational end, and the other met the arithmetic-genus-one component.
These sections stay disjoint and meet smooth points of the contracted
fibres. Indeed the first contraction moves the end's section to a
transverse intersection with the exceptional curve. The second moves
that section to a transverse intersection with the surviving component
at the former attaching point, away from the other section.

Choose the unblown section as origin \(O\) and write the blown section
as \(T=(u,v)\). The two Weierstrass models are
\[
                 Y^2=X^3+a X+b ,
   \qquad (a,b)=(0,p^2-1)
   \quad\text{or}\quad
   (a,b)=(-3(p^2-1),\,2p(p^2-1)).
                                                        \tag{1}
\]
In particular \(u,v\) are regular at \(p_0\), and \(T(p_0)\)
does not meet the cusp if the fibre there is singular.

## Principal parts at the two sections

Write \(\omega=dX/(2Y)\) for the invariant differential. By the
[marked factor theorem](KELLER_TWO_POLE_ELLIPTIC_SLICE_MARKED_FACTORS.md#the-singular-alternating-profile-at-width-four)
there is a rational function \(f\), a nonzero scalar multiple of
the reciprocal of the marked polynomial factor, with
\[
                            dq/\omega=f^2 .
                                                        \tag{2}
\]
The choice of \(\omega\) is valid in (2) since its product with
\(dp\) trivializes the canonical bundle of the smooth finite
Weierstrass surface, as does the volume form in that theorem.
Their quotient is a unit on the finite-base elliptic surface,
hence a nonzero constant.

Use a local fibre parameter \(t\) at either section with
\(\omega=dt\) on the fibre. The complementary coordinate \(q\)
has order-three poles on the original sections.
Near the unblown section its leading polar coefficient remains
a unit at \(\epsilon=0\).
Near the blown section the two inverse blowups, both at
the section in the special fibre, replace \(t\) by
\(\epsilon^2 t_2\) in the chart of its strict transform.
On that chart \(q\) has a pole of order three with leading
coefficient a unit, by the
[finite cover at infinity](KELLER_SMOOTH_NORMALIZATION_POLAR_BOUNDARY.md#compactification-over-the-p-line).
Thus in the parameter \(t\) at \(T\) the coefficient of \(t^{-3}\)
has order exactly six in \(\epsilon\), and the coefficient of
\(t^{-2}\) has order at least four. By (2), the principal part
of \(f\) at \(O\) has a unit coefficient at \(t^{-2}\). At \(T\)
its \(t^{-2}\) coefficient has order three in \(\epsilon\),
and its \(t^{-1}\) coefficient has order at least one, by comparing
the \(t^{-3}\) coefficients of \(dq/\omega\) and \(f^2\).

## The derivative-square equations

On the generic fibre put
\[
 s=\frac{Y+v}{X-u},\qquad X_T=s^2-X-u,\qquad
 Y_T=-Y+s(X-X_T).
\]
Thus \(X_T,Y_T\) are the coordinates after translation by \(-T\).
The function \(s\) has simple poles at \(O,T\), with principal
parts \(-1/t,1/t\) respectively. Write \(\delta=d/\omega\)
for the invariant derivative on the fibre. Translation gives
\(\delta X_T=2Y_T,\ \delta Y_T=3X_T^2+a\), and
\(\delta s=X-X_T\). The last equality also follows by comparing
principal and constant terms at the two poles.

The functions \(X,X_T,s,1\) form a basis for \(L(2O+2T)\)
by Riemann–Roch and their principal parts. Similarly
\(Y,X,Y_T,X_T,s,1\) form a basis for \(L(3O+3T)\).
Write
\[
                      f=L X+M X_T+N s+D .
\]
The orders proved above say that \(L\) is a unit at \(p_0\),
\(M\) has order exactly three there, and \(N\) has order at
least one. These conclusions use the invariant local
parameters and therefore also hold when a section passes
through a smooth point on a singular fibre.

The local expansions at \(O\) needed below are
\[
\begin{aligned}
 X&=t^{-2}-(a/5)t^2+O(t^4),&
 Y&=-t^{-3}-(a/5)t+O(t^3),\\
 X_T&=u-2v t+(3u^2+a)t^2+O(t^3),&
 s&=-t^{-1}-u t+v t^2+O(t^3).
\end{aligned}
\]
They follow from the Weierstrass equation, the invariant
derivative, and the formula for \(s\). At \(T\) the expansions
are obtained by translation, exchanging \(X,Y\) with \(X_T,Y_T\)
and replacing \(v,s\) by \(-v,-s\).

The absence of a \(t^{-1}\) term in a derivative gives the
first equation below. The \(t^{-4},t^{-3}\) terms at both
poles force the \(Y,X,Y_T,X_T\) coefficients in \(q\) to be
\(L^2/3,LN,M^2/3,-MN\), respectively. Comparing the
\(t^{-2}\) terms gives the second equation, because the
coefficient of \(s\) in \(q\) is both
\(N^2+2L(Mu+D)\) and \(-N^2-2M(Lu+D)\).
Finally comparison of the constant term at \(O\), eliminating
\(N^2\) with the second equation, gives the third:
\[
\begin{aligned}
 2LMv+N[(L+M)u+D]&=0,\\
 N^2+2LMu+(L+M)D&=0,\\
 D^2-(L+M)uD+2LMu^2+2(L+M)Nv
                  +(2LM-(L^2+M^2)/3)a&=0 .
\end{aligned}                                                \tag{3}
\]
In the last comparison, before eliminating \(N^2\), the sum of
the terms involving \(uD,u^2,N^2\) in \(f^2-\delta q\) is
\(2(L+M)uD+8LMu^2+3N^2u\).

## Contradiction at the centre of the blowups

The second equation makes \(D\) regular and zero at \(p_0\).
The third then requires \(a(p_0)=0\). In the second model in (1)
this forces \(p_0=\pm1\). Since \(T(p_0)\) avoids the cusp,
\(u\) is a unit in that model.

In either model, if \(u\) is a unit, the first two equations
imply
\[
                 N=-2M v/u+O(\epsilon^6),\qquad
                 D=-2M u+O(\epsilon^6).
\]
Indeed the coefficient of \(N\) in the first equation is a unit,
so first \(N=O(\epsilon^3)\), and then the two leading terms follow.
In the third equation these substitutions, together with
\(v^2=u^3+a u+b\), give
\[
                 -\frac{L^2}{3}a -2LM a-4LM b/u
                                      +O(\epsilon^6)=0.
\]
In the second model \(a\) has a simple zero at \(p_0\), which
is impossible. In the first model \(a=0\), while \(b\)
has order either zero or one, so the leading order is three or
four, again impossible.

It remains to consider the first model with \(u(p_0)=0\).
Then \(v\) is a unit since \(T\) avoids a cusp.
The function \(u\) is not identically zero, since \(p^2-1\) is
not a rational square. Let \(t_0\ge1\) be its order at \(p_0\).
If \(N\) vanished identically the first equation would fail.
Put \(n_0\ge1\) for its order. The second equation gives
\(\operatorname{ord} D\ge\min(2n_0,3+t_0)\). In the third
equation, if \(n_0<3+2t_0\), each term other than \(2(L+M)Nv\)
has order strictly greater than \(n_0\). Thus
\(n_0\ge3+2t_0\). But now the first equation has its term
\(2LMv\) of order three and its other term of strictly greater
order, the final contradiction.

## Directional and total degree bounds

**Corollary SixSheetSingularAlternatingWidthAtLeastFive.**
In the six-sheet \(A_6\) profile with one unramified omitted prime,
every polynomial target coordinate has directional width at least five.
The sole ramified target curve has total degree at least six in every
polynomial target frame. If the directional width is five, the
generic slice genus and pole count are \((1,4)\) or \((2,2)\).

The previous [directional bound and genus theorem](KELLER_SIX_SHEET_CUSP_SLICE_GENUS.md#singular-normalization-with-one-ordinary-affine-cusp)
give width at least four and genus at least one. At width four the
[marked factor](KELLER_TWO_POLE_ELLIPTIC_SLICE_MARKED_FACTORS.md#the-singular-alternating-profile-at-width-four)
forces the pole pair \((1,5)\) or \((3,3)\). The
[simple-pole exclusion](KELLER_SINGULAR_ALTERNATING_SIMPLE_POLE_EXCLUSION.md)
rules out the first pair, and the theorem above rules out the second.

In any target frame, the total degree of the ramified curve is the
maximum of the two coordinate degrees on its affine-line normalization.
This follows by taking the preimage of a generic line under the
birational polynomial parametrization. If the two degrees are equal,
cancel their leading coefficient by an affine target change.
The smaller resulting degree is at least five by the directional
bound, so their maximum is at least six.
Finally, the [polar degree formula](KELLER_DIRECTIONAL_RAMIFICATION_AND_POLAR_ENDS.md#six-sheet-consequences)
gives \(2w(p)=4+2g_p+r_p\). At width five this leaves
\(2g_p+r_p=6\), and \(g_p\ge1,\ r_p\ge1\) yield the stated pairs.

**Corollary SixSheetSingularAlternatingTangentSliceGenus.**
In the same singular alternating profile, a coordinate tangent
to the affine cusp has generic compactified slice genus at least
two, at any directional width. Thus at width five the slice
with genus and pole count \((1,4)\) has transverse cusp projection.

Indeed, the
[local resolution calculation](KELLER_SIX_SHEET_CUSP_SLICE_GENUS.md#singular-normalization-with-one-ordinary-affine-cusp)
makes the only singular point on the exceptional curve an
ordinary triple point for a tangent coordinate. All other
possible singularities of finite fibres are immersed index-three
germs \(t^m+v^3=0\). If the generic genus were one, the
[normalization genus inequality](KELLER_SIX_SHEET_CUSP_SLICE_GENUS.md#fibre-genus-and-canonical-bundle)
would allow only ordinary cusps and ordinary triple points.
Each contributes one to the inequality, so there is at most
one in a fibre. An ordinary cusp is unibranch, hence its
connected fibre is integral. At a triple point a genus-one
fibre must have three smooth rational components, since its
delta invariant is three and three is the maximum number of
components meeting there. Each component has self-intersection
\(-2\), so adjunction makes the dualizing bundle have degree
zero on each. A nonzero dualizing section trivializes the bundle
because the restrictions match at the common point.
The integral genus-one fibres likewise have trivial dualizing
bundles. The canonical class contradiction on the smooth locus
of the affine normalization then applies as in the cited
genus argument. This proves the additional transversality
assertion without a width hypothesis.

The [width-five marked elliptic exclusion](KELLER_WIDTH_FIVE_ALTERNATING_DOUBLE_SIMPLE_EXCLUSION.md#consequences-for-width-and-genus)
now excludes the width-five pair \((1,4)\). A genus-one
coordinate in the singular alternating profile has width at
least six.
