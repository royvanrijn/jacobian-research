# The pole markings in a width-five alternating elliptic slice

**Theorem SixSheetWidthFiveAlternatingEllipticMarks.** Consider the
six-sheet \(A_6\) profile with one unramified omitted prime. Suppose
the target coordinate \(p\) has width five, the compactified generic
slice has genus one, and the three critical points of the projection
on the immersed locus of the ramified branch are simple. Then the
relatively minimal finite-base elliptic model has equation
\[
                           Y^2=X^3+D(p)
\]
with \(D\) a squarefree cubic. It has three finite fibres of type
\(II\) and a fibre of type \(I_0^*\) at infinity. The four poles
of the companion coordinate have orders \(3,1,1,1\). The rational
end in the resolved fibre above the projection of the affine cusp
contains exactly one pole, of order one. In particular the triple
pole remains on the surviving arithmetic-genus-one component,
together with two simple poles.

## The relative model and its markings

Use the [width-five bound and transverse projection](KELLER_ALTERNATING_ELLIPTIC_TRIPLE_POLES.md#directional-and-total-degree-bounds)
and the [singular relative compactification](KELLER_SIX_SHEET_CUSP_SLICE_GENUS.md#singular-normalization-with-one-ordinary-affine-cusp).
The polar count is four, since \(2w(p)=6+2g_p+r_p-2\).
The projection of the affine cusp has local degree two.
On the resolution the exceptional \((-2)\)-curve \(D_{\rm exc}\)
meets the strict fibre in two distinct ordinary nodes, with no
other branches through it. Away from this curve all singularities of
finite fibres are unibranch cusps, by the immersed index-three
local power map and the simple-critical-point hypothesis.
Two such cusps cannot lie in the same genus-one fibre: the
normalization genus formula bounds the sum of \(\delta_x-b_x+1\)
by one, with \(b_x\) the number of analytic branches at \(x\).
Thus the three immersed cusps have distinct \(p\)-values.

All fibres except the fibre over the affine cusp projection
\(p_0\) are integral. At \(p_0\), the two strict branches
attaching to \(D_{\rm exc}\) are either on the same component
or on distinct components. In the first case the fibre has two
smooth rational components meeting twice, by the arithmetic
genus formula. Its dualizing bundle is trivial. Every other
fibre also has trivial dualizing bundle by integrality and
arithmetic genus one. Cohomology and base change would make the
relative dualizing bundle trivial. Restricting to the smooth
locus of the affine normalization contradicts the
[free horizontal ramification class](KELLER_SMOOTH_NORMALIZATION_POLAR_BOUNDARY.md#canonical-class-on-the-finite-affine-normalization).
This uses the same class calculation on the smooth locus of
the normal surface, since the removed singular point has
codimension two.

In the other case the fibre is a chain with \(D_{\rm exc}\)
between a smooth rational \((-1)\)-end \(A\) and an integral
arithmetic-genus-one component \(B\). In fact its only
intercomponent intersections are the two transverse nodes.
The arithmetic genera of the strict components sum to one,
and the \((-1)\) assertion follows by intersection with the
full fibre. Contract \(A\), then the image of \(D_{\rm exc}\).
Each contraction leaves the attaching point on the surviving
component smooth. Thus it preserves the finite cusps, including
a cusp on \(B\) if an immersed critical value is \(p_0\).
There are precisely three finite fibres of type \(II\) on
the resulting relatively minimal model \(S_0\).

The smooth relative surface before contraction has free Picard
group on the two omitted prime transforms, the exceptional
curve, and the four pole sections, by localization from the
original affine plane and constancy of its units. After the
two contractions this rank is five. The completion of the
relatively minimal model is a rational elliptic surface with
a section and has Picard rank ten. Units on \(S_0\) are constant,
as they pull back to units on the original affine plane.
Divisor localization thus gives five components at infinity.
The infinity fibre has Euler number six since the three
finite type-\(II\) fibres contribute six in total. By
[Miranda's Table 1.1](sources/miranda-persson-rational-elliptic-fibres-1990.pdf)
this is type \(I_0^*\). There are no poles of the \(j\)-function,
and a type-\(II\) value makes it zero. Choosing any image pole
section as origin and using the usual rational elliptic
Weierstrass degree bounds gives \(Y^2=X^3+D(p)\);
\(D\) has one simple zero at each of the three type-\(II\)
values and no other zeros.

The poles on \(A\) become sections through the same smooth
point \(T_0\) on the contracted fibre \(B\). The other sections
stay disjoint and avoid that point. Two former \(A\)-sections
have contact exactly two: in a local parameter \(s\) vanishing
on one of them and base coordinate \(\epsilon=p-p_0\), the
reverse blowups have coordinates \(s=\epsilon s_1=\epsilon^2
s_2\). The sections meet \(A\) at distinct finite values of \(s_2\).
This local description applies to each former \(A\)-section
separately by taking \(s=0\) on that section.

The
[marked factor calculation](KELLER_TWO_POLE_ELLIPTIC_SLICE_MARKED_FACTORS.md#the-singular-alternating-profile-at-width-four)
applies to this chain as well: the pulled-back regular
two-form from \(S_0\) has divisor \(2A+D_{\rm exc}\), and
its ratio to the pulled-back target form is proportional
to \(g^2\) for the irreducible polynomial defining \(A\cap X\).
On the generic fibre, for the invariant differential \(\omega\)
of this Weierstrass model, write
\[
                       dq=\varphi^2\omega ,
\]
where \(\varphi\) is a nonzero constant multiple of \(1/g\).
Both invariant volume forms on \(S_0\) differ by a constant
unit. The pole orders of \(q\) are odd, since \(dq/\omega\)
is a square and \(\omega\) is nonvanishing on that fibre.
They sum to six, hence are \(3,1,1,1\).

Both restrictions \(q|_A,q|_B\) are nonconstant since these
components dominate the target line of the finite normalization.
The map on \(A\) has only odd ramification: the pullback
of \(\omega\) along the relative blowup coordinate \(s_2\)
has leading term a nonzero constant times
\(\epsilon^2\,ds_2\), while \(\varphi\) has a simple pole
along \(A\). On \(A\), \(dq/ds_2\) is therefore a constant
times a rational square. The restriction to \(B\) has the
same parity property, since the dualizing differential of
a smooth elliptic fibre is nonvanishing, and on the
normalization of a cuspidal fibre has a double pole at the
cusp and no other zeros or poles. These parity statements
rule out degree two on either component by Riemann–Hurwitz.
Degree one is impossible on \(B\), since a degree-one finite
morphism from an integral arithmetic-genus-one curve to the
normal projective line would be an isomorphism. The possible
groups on \(A\) are now one simple mark, the sole triple mark,
or all three simple marks.

## Integral sections of the \(I_0^*\) model

We need the following elementary restrictions on an integral
section \(T=(u(p),v(p))\) of \(Y^2=X^3+D(p)\), where \(D\)
is a squarefree cubic. Here integral means disjoint from the
chosen origin section over the finite line, so \(u,v\) are
polynomials. A constant \(u\) is impossible since \(u^3+D\)
would have odd degree three. For \(\deg u\ge2\) we have
\(\deg u=2m,\deg v=3m\). The identity
\[
                 3u'D-uD'=v(3u'v-2uv')
\tag{1}
\]
bounds \(\deg u\le4\): for \(\deg u>4\) the degree on the
left is less than \(\deg v\), forcing \(u^3/D\) constant,
impossible by degree. Thus the possible degrees of \(u\)
are one, two, and four.

For degree four, center \(D=a t^3+b t+c\) by a base
translation, and write \(u=d t^4+u_3t^3+u_2t^2+u_1t+u_0\).
In (1) the multiplier \(3u'v-2uv'\) is a nonzero constant
since \(\deg v=6\). The leading coefficients give
\[
                   d(3u'D-uD')^2=81a^2(u^3+D).
\]
Its degree-eleven through degree-eight coefficients give
\[
         u_3=0,\quad u_2=22bd/(21a),\quad
         u_1=8cd/(9a),\quad u_0=41b^2d/(441a^2).
\]
The degree-six coefficient is then
\(-640 b^3d^3/(343a)\), so \(b=0\). The
[coefficient checker](verify_keller_width_five_alternating_elliptic_marks.py)
replays these comparisons. Here \(c\ne0\) since \(D\)
is squarefree. Thus the degree-four sections all have
\(u=d(t^4+8ct/(9a))\), and \(v\) is proportional to
\(27a^2t^6+36act^3+8c^2\) by (1).
Two with different \(u\) can meet only at a simple zero
of their \(u\)-difference, so cannot have contact two.
Two with the same \(u\) and opposite \(v\) can meet only
at a simple zero of \(v\): the quadratic in \(t^3\) has
distinct nonzero roots. Hence any contact-two group
contains at most one degree-four section.

## The square differential at three simple sections

In either of the two cases still to be excluded, choose
the triple pole as origin \(O\), and denote the three
simple sections by \(T_i=(u_i,v_i)\).
The sections are integral in this Weierstrass frame
because none meets \(O\) over the finite base.
On the generic elliptic curve put
\(f_i=(Y+v_i)/(X-u_i)\). Its only poles are simple at
\(O,T_i\) and \(\operatorname{Res}_{T_i}(f_i\omega)=1\)
for \(\omega=dX/(2Y)\). The constant term of \(f_i\)
in the invariant local parameter at \(T_i\) is zero:
if \(v_i\ne0\), expand \(X-u_i=2v_i s+(3u_i^2)s^2+
O(s^3)\) and \(Y+v_i=2v_i+3u_i^2s+O(s^2)\).
If \(v_i=0\), translation by the two-torsion point
gives the same result by oddness in the invariant
parameter. The functions \(1,X,f_1,f_2,f_3\) are a
basis for \(L(2O+T_1+T_2+T_3)\) by their principal
parts and Riemann–Roch. Write
\[
                     \varphi=U X+\sum_i C_i f_i+L.
\]
A mark not on \(A\) gives a nonzero constant leading
coefficient of \(\varphi\) there. A mark on \(A\) where
\(q\) has order \(s\) gives a leading coefficient of
\(\varphi\) vanishing exactly to order \(s\) at \(p_0\).
Indeed on the original pole section \(q\) has a
nowhere-zero leading polar coefficient by the finite
power-map chart. In the minimal model the two reverse
blowups scale the local normal parameter by
\(\epsilon^2\), so this coefficient of \(q\) on the
minimal section vanishes to order \(2s\), and the
square-differential equation halves that order.
These polar coefficients are regular along the whole
finite base and have no other zeros. The invariant
normalization of \(X\) at \(O\) and of the residues
of \(f_i\) at \(T_i\) makes them precisely the
coefficients \(U,C_i\), up to constants.

Since \(\varphi^2\omega\) is exact, its residue at
each \(T_i\) is zero. The preceding constant-term
calculation gives
\[
             L+Uu_i+\sum_{j\ne i} C_j f_j(T_i)=0 .
\tag{2}
\]
When \(u_i\ne u_j\) the evaluation is
\((v_i+v_j)/(u_i-u_j)\). When \(T_i=-T_j\) on the
generic elliptic curve it is the regular value
\(-3u_j^2/(2v_j)\). These formulas suffice to compare
growth at infinity even if the sections have poles
there in the affine Weierstrass chart.

### All three simple marks on the rational end

Now \(U\) is a nonzero constant, \(C_i=h_i\epsilon\)
with \(h_i\ne0\), and the three sections have pairwise
contact exactly two at \(p_0\). There is at most one
degree-four section. If there is one, evaluation
between it and another section grows at infinity
at most as \(p^2\), since the denominator has degree
four and the numerator degree at most six. Between
the two lower-degree sections, distinct \(u_i,u_j\)
differ by a nonzero multiple of \(\epsilon^2\).
The evaluation then grows at most as \(p\).
In the opposite-section case with equal \(u_i\),
contact two forces \(v_i\) to have a multiple zero,
so \(\deg u_i=2,\deg v_i=3\), and the regular-value
formula has the same bound. In (2) the degree-four
section would make \(L\) have growth of degree four,
while either other section bounds it by degree three.

Otherwise all three sections have degree at most
two. Sharing their first jets at \(p_0\), write
\[
 u_i=u_0+u_1\epsilon+\alpha_i\epsilon^2,\qquad
 v_i=v_0+v_1\epsilon+\beta_i\epsilon^2+\gamma_i\epsilon^3.
\]
If \(v_0=0\), then \(u_0\ne0\) because the section
avoids a cusp. The \(\epsilon^2\) curve equation
forces all \(\alpha_i\) equal, giving at most two
distinct sections. Therefore \(v_0\ne0\), and the
\(\alpha_i\) are distinct because identical \(u_i\)
would give opposite ordinates on distinct sections.
The coefficients of orders two and three make
\(\beta_i=\lambda\alpha_i+\beta\),
\(\gamma_i=\mu\alpha_i+\gamma\), with
\(\lambda=3u_0^2/(2v_0)\). Order five gives
\(2\beta_i\gamma_i=3u_1\alpha_i^2\);
this quadratic identity at three distinct values
gives \(\lambda\gamma+\mu\beta=0,\ \beta\gamma=0\).
If \(u_0\ne0\), then \(\lambda\ne0\).
If \(\mu=0\) this gives \(\gamma=0\), and order six
would give all \(\alpha_i=0\). If \(\mu\ne0\),
then \(\beta=\gamma=0\), and order six gives only
\(\alpha_i=0\) or \(\mu^2\). Hence \(u_0=0\).
Now \(\lambda=0\), the leading coefficient at order
five forces \(u_1=0\), and then the order-three
comparison gives \(\mu=0\). Thus \(u_i=\alpha_i\epsilon^2\)
and all \(v_i\) are the same polynomial \(v\), with
\(v(p_0)\ne0\).

In (2) the sum is now
\[
                   \frac{2v}{\epsilon}
                        \sum_{j\ne i}\frac{h_j}{\alpha_i-\alpha_j}.
\]
The function \(L\) is regular at \(p_0\), because it is
the constant term of \(\varphi\) in the invariant
parameter at the unblown origin, after subtracting
\(UX+\sum C_i f_i\). The functions \(X,f_i\) have zero
constant term there. Thus each displayed residue
sum is zero, and (2) would give \(L=-Uu_i\) for all
three distinct \(u_i\). This is impossible.

### The lone triple mark on the rational end

Now \(U=h\epsilon^3\), \(h\ne0\), and the three \(C_i\)
are nonzero constants. The sections \(T_i\) stay
pairwise disjoint at finite points. If one has degree
four, (2) forces \(L\) to have growth degree seven
with leading coefficient \(-h\) times that section's
leading \(u\)-coefficient. Indeed every evaluation
grows by degree at most six, or at most two for
opposite sections of degree four. Therefore all three
must have degree four and the same leading \(u\)-coefficient.
The explicit degree-four form then makes all their
\(u_i\) identical, which permits at most two distinct
sections.

If all have degree at most two and one has degree
two, then all must have degree two, since evaluations
grow by degree at most three while \(Uu_i\) has degree
five in degree two and four in degree one. Subtracting
(2) for two indices now forces each \(u_i-u_j\) to be
constant. Identical \(u_i\) for distinct sections would
make their ordinates opposite; the cubic ordinate has
a finite zero, giving a forbidden meeting. Thus the
three constant terms are distinct. Center the common
quadratic part and write \(u_i=a t^2+c_i\),
\(v_i=d_i t^3+e_i t^2+f_i t+g_i\). The degree-six and
degree-five equations give \(d_i^2=a^3,\ e_i=0\), and
degree four gives \(f_i=3a^2c_i/(2d_i)\). Degree two
then makes \(c_i^2\) independent of \(i\), allowing at
most two distinct constant terms. This too is impossible.
If all sections have degree one, the evaluations
grow by degree at most two, including at opposite
sections. Subtracting (2) would force all \(u_i\)
identical and again permit at most two distinct sections.

The only remaining group on \(A\) is one simple mark.
