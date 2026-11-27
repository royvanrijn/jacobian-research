# The quartic sequence \((14,4,11)\) in conductor degree twenty-two

**Theorem ConductorTwentyTwoQuarticTetradecicExclusion.** Let \(F=(P,Q)\)
be a complex polynomial plane Keller map with \(P\) monic of positive
degree in the longitudinal coordinate \(x\). Its generic directional
coordinate order cannot have reduced characteristic sequence
\((14,4,11)\).

## The square-root tail

Put \(K=\mathbb C(y)\) and \(J=\operatorname{Jac}_{x,y}(P,Q)\ne0\).
The [exact coordinate-order theorem](GENERIC_COORDINATE_ORDER_CONDUCTOR_AND_CLOSED_COLLISIONS.md)
gives generic colength eleven and a nonconstant specialized conductor
on every finite line. At each finite place there is a matching pair
with distinct finite physical limits. A generic conductor-root branch
approaches a finite root of the specialized conductor. Monicity of
\(P\) keeps the other endpoint finite, and equal limits would force
\(P_x=Q_x=0\) by the divided differences. In particular, some physical
gap has valuation zero at each finite place.

As in the [quartic–decic reduction](CONDUCTOR_TWENTY_QUARTIC_DECIC_EXCLUSION.md#1-the-cubic-correction-and-the-parameter-invariants),
the [one-place approximate-root theorem](sources/assi-garcia-sanchez-1407.0490v1.pdf)
gives a degree-eleven approximate root, monic quadratic in the
degree-fourteen generator over the polynomial ring in the quartic
generator. Complete the square by a triangular target change. Translate
the quartic to remove the sextic term in the monic degree-seven
polynomial in the resulting approximate root. Terms of degree at most
two in that polynomial contribute at most eight to the longitudinal
degree and can be discarded. Depress the monic quartic in \(x\) by a
translation \(h(y)\in K\). If \(\lambda\in K^\times\) is the determinant
of the target changes, put \(\rho^{18}=\lambda\), \(x=h+\rho z\), and
multiply the two generators by \(\rho^{-4},\rho^{-14}\), respectively.
In these normalized generators write
\[
 U=z^4+az^2+bz+c,\qquad \Delta=a^2-4c,\qquad
 W=V^2-R(U),\quad \deg_z W=11,
\]
where
\[
 R(U)=U^7+\alpha U^5+\beta U^4+(\gamma+\alpha^2/4)U^3 .
\]
The monic Laurent square root \(S_0=\sqrt{R(U)}\) has expansion
\[
 S_0=U^{7/2}+\frac{\alpha}{2}U^{3/2}
             +\frac{\beta}{2}U^{1/2}+\frac{\gamma}{2}U^{-1/2}
             +O(z^{-6}).
\]
The highest exponent of \(V-S_0\) is \(-3\). Thus \(V\) is the
polynomial part of \(S_0\), its coefficients at \(z^{-1},z^{-2}\) vanish,
and \(E=[z^{-3}]S_0\ne0\).

The normalized target change has determinant one. The
[parameter identity](CONDUCTOR_TWELVE_CABLE_EXCLUSIONS.md#1-common-source-and-target-identities)
gives
\[
 U_zV_y-U_yV_z=J\rho+\partial_z H(U,V).
\]
In \(U_z(S_0)_y-U_y(S_0)_z\), the terms from differentiating \(U\)
cancel. The possible leading powers from \(\alpha',\beta',\gamma'\)
are \(9,5,1\), respectively, each with nonzero scalar coefficient.
Replacing \(S_0\) by \(V\) changes only the constant and negative
powers. The only positive nongaps of \(\langle4,14,11\rangle\) at most
ten are four and eight. Thus the coefficients at \(9,5,1\) successively
force \(\alpha',\beta',\gamma'=0\). The degree-zero term, from the
\(y\)-derivative of \(V-S_0=-Ez^{-3}+O(z^{-4})\), gives
\[
                     -4 E'=J\rho. \tag{1}
\]

Put \(X=b^2\), \(Y=ab^2\), and \(L=56Y-64\alpha-35\Delta^2\). The two
vanishing Laurent coefficients and the expression for \(E\) give
\[
\begin{aligned}
 35\Delta^3-280Y\Delta+192\alpha\Delta+56X^2-512\beta&=0,\\
 2240Y^2+3920Y\Delta^2-3072\alpha Y-3456\alpha\Delta^2
                  -665\Delta^4+8192\beta\Delta+16384\gamma&=0,\\
                 2048E &= b^3 L . \tag{2}
\end{aligned}
\]
These identities follow by binomial expansion and are verified by the
[exact checker](verify_conductor_twenty_two_quartic_tetradecic_exclusion.py).
In particular \(bL\ne0\). A constant \(\Delta\) forces \(Y\), \(X\), and
then \(E\) constant by (2), contrary to (1).

## Leading shapes and finite gaps

At any pole of \(\Delta\), the second equation of (2) forces \(Y\)
to have twice the pole order of \(\Delta\). Its leading ratio
\(\kappa=Y/\Delta^2\) is a root of \(64\kappa^2+112\kappa-19=0\).
The first equation gives \(b^4/\Delta^3\to5\kappa-5/8\ne0\).
Since \(\kappa\ne0\), the pole orders of \(a,b,\Delta\) are in the
ratio \(2:3:4\), and \(L/\Delta^2\to56\kappa-35\ne0\).
For \(\Delta\) of pole order \(4q\), \(E\) has pole order \(17q\).

Rescale \(z\) by a factor of pole order \(q\), and \(U,V,W\) by its
inverse fourth, fourteenth and twenty-eighth powers. The limiting
quartic and degree-fourteen generator are obtained by setting
\(\alpha=\beta=\gamma=0\) after rescaling. Set the nonzero limiting
\(a\) equal to one by a constant scale. The leading relations give
\[
 \Delta=\frac{1600\kappa+2040}{361},\qquad
 b^2=\frac{70400\kappa+32000}{6859}.
\]
The checker computes the resultant of the quartic derivative and the
remainder of \(V'\) by that derivative. After division by \(b\) and
reduction modulo \(64\kappa^2+112\kappa-19\), it is a nonzero scalar
multiple of \(20161770792104\kappa-6158358133417\), coprime to the
quadratic. Thus the limiting plane parametrization is immersed.

The limiting \(W\) has degree eleven, with leading coefficient
\(-2E\) in the scaled variables. Matching endpoints stay bounded in
the rescaled coordinate. Indeed, escaping endpoints would produce a
nonzero common pair for the divided differences of the leading powers
of \(U,W\), of coprime degrees four and eleven. Their limits on the
nonzero diagonal are nonzero as well. Limiting immersion prevents
coalescence. Consequently every \(z\)-gap has pole order \(q\). If
the place is finite, \(\rho\) has a pole by (1) and the pole of \(E\).
Every physical gap then diverges, a contradiction. Hence \(\Delta\)
has no finite pole, and \(Y,X,L\) are regular at all finite places by
(2).

For any matching pair let \(s\) be the sum and \(g\) the difference
of its endpoints in \(z\). Since \(b\ne0\), the quartic divided
difference gives \(s\ne0\), the product
\((s^3+as+b)/(2s)\), and \(g^2=-s^2-2a-2b/s\).
Set \(k_0=s^2+2a\). Substitute this product in the divided difference
of \(V\), and use the first equation of (2). The checker verifies the
resulting sum equation
\[
\begin{split}
 0={}&8s^5 k_0^3+28\Delta s^3 k_0^2-8bs^4 k_0^2
        -28b\Delta s^2 k_0-16b^2s^3 k_0-Ls k_0\\
     &+8b^3s^2-28b^2\Delta s+bL. \tag{3}
\end{split}
\]
At a finite place of \(K(\rho)\) use integer valuation \(v\),
and put \(k=v(b)\ge0\), \(A=\max(0,-v(a))\le2k\) and
\(\ell=v(L)\ge0\). Equation (3) is monic up to a nonzero scalar in \(s\).
If \(A=0\), all coefficients are regular and \(s\) is bounded. For
\(A>0\), replace \(s\) by \(\sqrt a\) times a new variable, divide (3)
by \(a^{11/2}\), and take the limit. Its leading part is
\(8s_{\rm new}^{\,5}(s_{\rm new}^{\,2}+2)^3\), with regular scaled
coefficients. In both cases \(v(s)\ge-A/2\).

If \(v(s)\le k+A\), this already gives \(v(g^2)\ge-A\). If
\(v(s)>k+A\), divide (3) by \(s\). The term \(Lk_0\)
has valuation at least \(\ell-A\), and all other terms apart from
\(bL/s\) have valuation at least \(2k\). This proves
\[
                   v(g^2)\ \ge\ \min(-A,2k-\ell). \tag{4}
\]
For gap extensions the same valuation is extended with its original
normalization. If \(3k+\ell>0\) and the ramification index over
\(\mathbb C(y)\) is \(e\), then \(v(\rho)=3k+\ell-e\) by (1).
The squared physical gap therefore has valuation at least
\[
               \min(6k+2\ell-2e-A,\ 8k+\ell-2e). \tag{5}
\]

## The inhomogeneous cases

The weights of \(\alpha,\beta,\gamma\) under longitudinal scaling
are \(8,12,16\). Their unscaled values lie in \(K\). Hence the cyclic
extension \(K(\rho)/K\), whose degree divides eighteen, has degree
dividing two if \(\alpha\ne0\) or \(\gamma\ne0\), and degree dividing
six if \(\beta\ne0\).

First suppose \(\alpha\ne0\) or \(\gamma\ne0\). The quantities
\(\Delta,Y,X,L\) have even weights and belong to \(K\). Since they
have no finite poles, they are polynomials in \(y\). Differentiating
the third equation of (2) makes \(\rho\) a nonzero constant multiple
of \(bT\), where
\[
                        T=\tfrac32 X'L+XL' .
\]
A root of \(T\) cannot be a root of \(X\). At a place over such a root
one would have \(k\ge1\) and \(v(\rho)>k\). With \(e\le2\) and
\(A\le2k\), the lower bound (5) is positive except possibly for
\(e=2,k=1,A=2,\ell=0\). In that exception \(v(\rho)=k=1\),
also a contradiction.

At a root of \(T\) where \(X\ne0\), one has \(k=A=0\) and \(e=1\).
Indeed, if the place ramifies, the inertia involution acts by sign
on \(b\) of weight three, which would force its finite residue to be
zero. If \(L\) is a unit, (4) bounds every gap, while \(\rho\) has a
zero. If \(L\) is zero, then \(v(\rho)=\ell-1>0\), and the finite unit
gap together with (5) gives \(\ell=2\). Thus every root of \(T\)
is simple and is a double root of \(L\). In particular
\(\deg T\le\deg L/2\). The nonconstant polynomial \(\Delta\) has a
pole at infinity, so its leading ratios give \(\deg X\ge3\) and
\(\deg L\ge4\). The leading terms in \(T\) give
\(\deg T=\deg X+\deg L-1\), a contradiction.

Next suppose \(\alpha=\gamma=0\) and \(\beta\ne0\). Here both
\(X=b^2\) and \(r=Y/\Delta^2\) are invariant under the group of degree
dividing six. Equations (2) give
\[
 \Delta^3=\frac{-8192\beta}{35(64r^2+112r-19)},\qquad
 X^2=\frac{64\beta}{7}\frac{(8r+1)(8r-3)}{64r^2+112r-19}.
\]
The double cover of the \(r\)-line has four distinct simple branch
points, at the two numerator roots and the two denominator roots.
Its smooth completion has genus one by Riemann–Hurwitz. Rational
\(r,X\) in \(y\) therefore give a constant map from
\(\mathbb P^1_y\) to this curve. A constant \(r\) forces \(\Delta\)
constant by the first equation, already excluded.

## The homogeneous case

Finally let \(\alpha=\beta=\gamma=0\). The homogeneous leading
relations now hold exactly. In particular \(\Delta\ne0\), and after
scaling by \(\tau\) with \(\tau^4=\Delta\), the coefficients of both
generators are complex constants. All matching gaps in \(z/\tau\)
are fixed and nonzero. The physical scale \(w=\rho\tau\) has
neither a zero nor a pole at any finite place by the finite unit
gap. Since \(w^4=\rho^4\Delta\in K\), it is a nonzero constant.
Write \(E=E_*\tau^{17}\) with \(E_*\ne0\). Equation (1) now makes
\(\tau^{18}\) affine nonconstant in \(y\). Consequently
\(\lambda^{-1}=\tau^{18}/w^{18}\) has a finite zero.

The unscaled generators in \(u=x-h(y)\) are fixed. Their matched
endpoints are fixed, so the finite endpoint condition forces \(h\)
to be regular at all finite places. In \(u\), the actual polynomials
\(P,Q\) belong coefficientwise to the fixed order. Indeed, membership
in each bounded-degree part of the fixed order is defined by
constant linear equations and holds generically. Thus \(P,Q\)
have polynomial ambient lifts in the fixed target generators. At
a finite zero of \(\lambda^{-1}\), their specialized order coincides
with the fixed order by containment and the fact that the colength
is at most the generic value eleven.

At the image of a fixed matched pair the fixed plane equation and its
first derivatives vanish. The regular ambient lifts and the generic
inverse target automorphism agree on the fixed curve, so their
determinants agree at this singular point. The determinant of the
inverse target automorphism is \(\lambda^{-1}\). By equality of orders,
the specialized lifts induce an isomorphism of image curves. The
[singular tangent determinant argument](HYPERELLIPTIC_CONDUCTORS_AND_SINGULAR_TANGENT_DETERMINANT.md#the-determinant-at-a-finite-singular-point)
then makes their derivative invertible on the full plane tangent
spaces, contradicting the finite zero. This completes the exclusion.
