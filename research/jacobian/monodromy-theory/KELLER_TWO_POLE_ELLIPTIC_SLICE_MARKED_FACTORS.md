# Polynomial factors in two-pole elliptic Keller slices

Let \(F:X=\mathbb A^2_\mathbb C\to\mathbb A^2_\mathbb C\) be a plane
Keller map of generic sheet degree \(N>1\), with full finite
normalization \(W\). In a polynomial target frame \(p,q\), let
\(g_p,r_p\) denote the genus and pole count of its compactified
generic \(p\)-slice. Where a normal model is singular, divisors
below are taken on its smooth locus.

**Theorem KellerSmoothTwoPoleEllipticSliceFactors.** Suppose \(W\)
is smooth and \((g_p,r_p)=(1,2)\). The finite-base elliptic surface
has at least one reducible fibre consisting of a rational
\((-1)\)-component attached transversely at one point to an
integral arithmetic-genus-one component. Let \(A_i\) be all the
rational components of this type. There is a polynomial \(G\)
on the original plane with simple divisor \(\sum_i(A_i\cap X)\).
It has order \(-(e_E-1)\) at each omitted prime \(E\), where
\(e_E\) is the normal ramification index. In particular the
marked rational components come from irreducible factors of
the corresponding values of \(p\circ F\), and the map
\((p\circ F,G)\) has generic sheet degree \(N+2\).

## Smooth normalizations

Use the
[finite-base surface](KELLER_SMOOTH_NORMALIZATION_POLAR_BOUNDARY.md#compactification-over-the-p-line)
\(S\to\mathbb A^1_p\), finite over \(\mathbb A^1_p\times\mathbb
P^1_q\), with its two pole sections \(H_1,H_2\). The surface
is smooth, its fibres are reduced and connected, and every
component of every fibre meets at least one section at a smooth
point of the fibre. Thus a reducible fibre has exactly two
components with one pole on each. Flatness gives arithmetic
genus one on all fibres.

Write a reducible fibre as \(A+B\). Its components are integral,
and
\[
                    1=p_a(A)+p_a(B)+A\cdot B-1 .
\]
The intersection is positive by connectivity. If \(A\cdot B=2\),
both components are smooth rational curves. Intersecting the
full fibre with each gives self-intersection \(-2\), so the
fibre dualizing bundle has degree zero on each projective line.
A nonzero global section exists by duality. Its restriction
to either component is zero or nowhere zero, and the restrictions
agree at every common point. Thus it trivializes the fibre
dualizing bundle, including when the intersection is tangential.

If \(A\cdot B=1\), one component is smooth rational and the
other has arithmetic genus one. Both are smooth and transverse
at their common point, since the intersection multiplicity there
is one. The rational component has self-intersection \(-1\),
by intersection with the full fibre. Contract it in a smooth
projective completion and restrict the resulting contraction
over the finite \(p\)-line, as in the
[polar boundary argument](KELLER_SMOOTH_NORMALIZATION_POLAR_BOUNDARY.md#low-genus-obstructions).
The surviving component remains integral and meets the
contraction center in a smooth point, by the transverse
intersection with the contracted curve.

These contractions are disjoint and finite in number,
by properness and generic smoothness of the fibration.
Write \(\rho:S\to S_0\) for their composition over the finite
\(p\)-line. The fibres of \(S_0\to\mathbb A^1_p\) are reduced,
connected and have arithmetic genus one, as any changed
fibre is now integral and the contraction preserves the
remaining component's curve germ. Every integral
arithmetic-genus-one fibre has trivial dualizing bundle,
by duality and degree zero. The other fibres have
trivial dualizing bundle by the degree-zero argument above.
The relative dualizing bundle is therefore trivial by
cohomology and base change and evaluation, since the base
has trivial Picard group. This is the same step as in the
[integral-fibre argument](KELLER_SMOOTH_NORMALIZATION_POLAR_BOUNDARY.md#low-genus-obstructions).
Because the base has trivial canonical bundle, \(S_0\) has a
nowhere-zero regular two-form \(\Omega_0\). Each inverse blowup
is at a smooth point, so its pullback satisfies
\(\operatorname{div}(\rho^*\Omega_0)=\sum_i A_i\).

The rational function
\[
                 G=\frac{\rho^*\Omega_0}{\pi^*(dp\wedge dq)}
\tag{1}
\]
has divisor \(\sum_i(A_i\cap X)\) on \(X\), as the denominator
is a nonzero constant times the original affine volume form.
By normality of \(X\) it is polynomial there. Every \(A_i\)
meets \(X\) by the finite-base surface argument, and its
intersection with \(X\) is a dense open part of an integral
fibre component. Its irreducible polynomial equation is a
factor of \(p\circ F-p(A_i)\). On \(W\) the denominator
in (1) has divisor \(\sum_E(e_E-1)E\), and the numerator
has only vertical zeros. This gives the stated orders at
the omitted primes. If there were no \(A_i\), the polynomial
would have no zeros on the plane, hence would be constant.
This contradicts a nonzero pole order at a ramified omitted
prime.

On the generic compact fibre choose the nowhere-zero
differential \(\omega_0\) induced by \(\Omega_0\) and \(dp\).
Then \(G=\omega_0/dq\). Away from the poles of \(q\) the
finite critical points of \(q\) lie over the omitted
ramified primes. A generic \(p\)-slice meets such primes
transversely, and the orders in (1) give a pole divisor of
degree \(R_p\) for \(G\). At a pole of \(q\) of order \(s\),
\(dq\) has pole order \(s+1\), so \(G\) has no additional
poles there. Thus the function field degree over
\(\mathbb C(p,G)\) is \(R_p=N+2g_p+r_p-2=N+2\).

## The singular alternating profile at width four

**Corollary SixSheetSingularAlternatingMarkedEllipticFactor.**
In the six-sheet \(A_6\) profile with one unramified omitted
prime, suppose \(w(p)=4\).
Let \(p_0\) be the projection of the affine cusp. The
polynomial \(p\circ F-p_0\) has an irreducible factor \(g\)
with order \(-1\) at the ramified omitted prime and order
zero at the unramified omitted prime. The restriction of
\(g\) to the generic compact elliptic fibre has degree four.
Both poles of \(q\) on that fibre have odd orders, so their
unordered orders are \((1,5)\) or \((3,3)\).

Use the
[resolved special fibre and its contractions](KELLER_SIX_SHEET_CUSP_SLICE_GENUS.md#singular-normalization-with-one-ordinary-affine-cusp).
On the resolved surface \(\widetilde S\) the special fibre
contains a rational \((-1)\)-end \(A\) and the exceptional
\((-2)\)-curve \(D\) of the \(A_1\) point, with \(A\) meeting
\(D\) transversely once. First contract \(A\), then the image
of \(D\), which becomes a \((-1)\)-curve. Denote the composition
by \(\rho\). The result \(S_0\) has integral arithmetic-genus-one
fibres over the finite \(p\)-line, so it again has a nowhere-zero
regular two-form \(\Omega_0\). The rule for the canonical
divisor under a blowup now gives
\[
                     \operatorname{div}(\rho^*\Omega_0)=2A+D .
\]
The coefficient two on \(A\) results from pulling back the
image of \(D\), which meets the center of the first blowup
with multiplicity one.

The ratio (1) on the original affine plane has divisor
\(2(A\cap X)\); it is therefore a nonzero constant times
\(g^2\) for an irreducible polynomial \(g\) defining \(A\cap X\).
Since the two rational functions agree on \(X\), this equality
also holds in the full source field. The curve \(D\) disappears
on the smooth locus of \(W\), and the denominator of (1)
has order two at the index-three prime and order zero at the
unramified prime. Hence \(g\) has the asserted orders and
divides \(p\circ F-p_0\). Its degree on the generic fibre
equals the degree four of the ramified prime's intersection
with that fibre, where its poles are simple.

On that generic compact fibre a nonzero constant times \(g^2\)
is \(\omega_0/dq\), with \(\omega_0\) nowhere zero. At a pole
of order \(s\) of \(q\), the order \(s+1\) of this ratio
must be even. Thus both \(q\)-pole orders are odd. They
sum to the sheet degree six, giving the two stated
possibilities. These conclusions are necessary conditions
on the original polynomial source in the low-width elliptic
residual. The [width-five lower bound](KELLER_ALTERNATING_ELLIPTIC_TRIPLE_POLES.md#directional-and-total-degree-bounds)
uses them to exclude the singular alternating width-four case.
