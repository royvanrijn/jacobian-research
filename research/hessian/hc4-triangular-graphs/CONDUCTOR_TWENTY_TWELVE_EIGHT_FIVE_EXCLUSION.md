# The octic–degree-twelve exclusion at conductor twenty

**Theorem ConductorTwentyTwelveEightFiveExclusion.** Let \(F=(P,Q)\)
be a complex polynomial plane Keller map with \(P\) monic of positive
degree in the longitudinal coordinate \(x\). Its reduced generic
directional order cannot have characteristic sequence \((12,8,5)\).

The proof determines the degree-five approximate root from a quartic
square root of the octic generator. A coefficient ratio of that
square-root decomposition would have to omit two values at every
finite parameter.

## The approximate-root equation

Write \(K=\mathbb C(y)\) and \(J=\operatorname{Jac}_{x,y}(P,Q)\ne0\).
The [exact coordinate-order theorem](GENERIC_COORDINATE_ORDER_CONDUCTOR_AND_CLOSED_COLLISIONS.md)
applies to every finite parameter. At each finite place at least one
generic matching pair has distinct finite limits in the original
coordinate \(x\), by the finite-pair argument in the
[degree-nine cable proof](CONDUCTOR_FOURTEEN_DEGREE_NINE_CABLE_EXCLUSION.md).
Indeed, a generic conductor root approaches a finite root of the
nonconstant specialized conductor. Its mate stays finite by monicity
of \(P\), and if the two limits agreed, the original divided differences
would force \(P_x=Q_x=0\). Thus some physical gap is a unit at every
finite place.

After triangular reduction and monic scaling, call the generators
\(W,V\) of degrees eight and twelve. The degree-five approximate
root is monic cubic in \(W\), by the one-place theorem as in the
[three-to-two Laurent obstruction](CONDUCTOR_FOURTEEN_CUBIC_ORDER_ALTERNATE_EMBEDDINGS.md#the-general-three-to-two-laurent-obstruction).
For weights \((8,12)\) the exponents zero, one and two in \(W\)
give distinct residues modulo twelve. Hence terms above weight
twenty-four cannot cancel in the approximate root, and the top
terms are \(W^3-V^2\). Complete the square in \(V\) by a linear
polynomial in \(W\) and translate \(W\) to remove the \(W^2\)
term on the right. Absorbing a constant in the approximate root
gives
\[
                         T=V^2-W^3-kW,\qquad \deg T=5 .    \tag{1}
\]

Use the source monic formal square root \(S=(W^3+kW)^{1/2}\).
Equation (1) gives \(S-V\) of degree exactly \(-7\) at infinity,
where degree of a Laurent series means its highest exponent.
Write \(W=L^2+M\) with \(L\) monic quartic and \(\deg M\le3\).
Expansion at infinity gives
\[
 \begin{split}
 V&=L^3+\tfrac32 LM+\tfrac38\operatorname{quot}(M^2,L),\\
 S-V&=\frac{3R/8+k/2}{L}
                 -\frac{M^3}{16L^3}+O(z^{-8}),\qquad
 R=\operatorname{rem}(M^2,L).
 \end{split}                                             \tag{2}
\]
In particular \((6R+8k)L^2-M^3\) has degree exactly five.
If \(\deg M\le2\), the \(L^2\) term would have to vanish
and \(M^3\) could not have degree five. Thus \(M=r m\) for
\(r\ne0\) and a monic cubic \(m\), and \(6R+8k\) is linear
with leading coefficient \(r^3\). Center the source at its
root \(h(y)\in K\), writing \(u=x-h(y)\). Then
\(6R+8k=r^3u\).

Let \(\lambda\in K^\times\) be the determinant of the target
change to this ordered pair \((W,V)\). Now take
\(u=\rho z,\ \rho^{20}=\lambda\), and scale the generators
by \(\rho^{-8},\rho^{-12}\). The identities (1)–(2) persist
with the scaled coefficients and \(6R+8k=r^3z\). In particular
\[
                           \deg(m^3-zL^2)=5 .             \tag{3}
\]

## Coefficients and the parameter Jacobian

Write \(m=z^3+az^2+bz+c\) and \(\Delta=a^2-4b\).
The coefficients of degrees eight through six in (3) give
\[
 L=z^4+\tfrac32 az^3+
 \left(\tfrac32 b+\tfrac38 a^2\right)z^2+
 \left(\tfrac32 c+\tfrac34ab-\tfrac1{16}a^3\right)z+\ell .
\]
The coefficients of degrees three and two in
\(R/r^2=\operatorname{rem}(m^2,L)\) vanish. They give
\[
 c=-a\Delta/8,\qquad
 \ell=(5a^4-24a^2b+80ac+16b^2)/64 .
\]
With these substitutions the degree-one and degree-zero
coefficients of \(R/r^2\) are \(a\Delta^2/128\) and
\(\Delta^2(2a^2+\Delta)/512\), and the degree-five coefficient
in (3) is \(\Delta^2/64\). Comparing \(6R+8k=r^3 z\) yields
\[
         r=\frac{3a\Delta^2}{64},\qquad
         k=-\frac{3r^2\Delta^2(2a^2+\Delta)}{2048}.         \tag{4}
\]
Thus \(a,\Delta\ne0\). Put \(G=r^3\Delta^2\); equation (2)
gives \([z^{-7}](S-V)=-G/1024\) and therefore
\([z^5]T=G/512\).

The determinant-one target normalization gives the
[parameter identity](CONDUCTOR_TWELVE_CABLE_EXCLUSIONS.md#1-common-source-and-target-identities)
\[
                  W_zV_y-W_yV_z=J\rho+\partial_z H(W,V)
\]
with \(H\) polynomial in the target variables. The contribution
from \(V-S\) to the polynomial part of the left side is
\(8([z^{-7}](V-S))'\), since \(W_z\) has degree seven.
The chain rule gives the other contribution as
\(k'(\sqrt W)_z\) in nonnegative Laurent degrees. Thus the
left side has degree at most three. The pole semigroup
\(\langle5,8,12\rangle\) has no positive degrees through four,
so the evaluated \(H\) is constant in \(z\). Degree three
forces \(k'=0\); the remaining constant term gives
\[
                           G'=128J\rho .                 \tag{5}
\]
The [exact checker](verify_conductor_twenty_twelve_eight_five_exclusion.py)
verifies the coefficient identities in (3)–(4).

## Excluding two finite values

Let \(X=\Delta/a^2\). This is a rational function of \(y\):
the unscaled centered coefficients \(a_0,\Delta_0\) lie in
\(K\), and \(a=a_0/\rho,\ \Delta=\Delta_0/\rho^2\).
If \(k\ne0\), equation (4) gives a nonzero constant multiple of
\(a^{16}X^6(2+X)\) as \(k\), and \(G\) is a nonzero scalar
multiple of \(a^{19}X^8\). Differentiating these identities gives
\[
                G'/G=\frac{(28-5X)X'}{16X(2+X)} .         \tag{6}
\]
At a finite parameter where \(X=28/5\), the coefficient \(a\)
is a unit by (4), and \(G\) is a unit. All normalized coefficients
are regular, and the polynomial functions \(W,T\) have unit
leading coefficients of coprime degrees eight and five.
Every normalized matched pair is therefore bounded. Otherwise
scaling an escaping pair by its largest coordinate would give
a nonzero simultaneous zero of the two divided differences
of the monic powers. Coprimality excludes this off the diagonal,
and the derivatives exclude it on the nonzero diagonal.
Since \(X'\) is regular at the finite parameter, (6) gives a
zero of \(\rho\). This shrinks every physical gap to zero,
contrary to the finite-pair condition.

At a finite parameter where \(X=-2\), equation (4) instead
makes \(a\) polar. Scale \(z=aw\) and scale \(W,V,T\) by
\(a^{-8},a^{-12},a^{-24}\). Their limits are the polynomials
in (2)–(4) at \(a=1,\Delta=-2\). Here \(r=3/16,\ k=0\).
The limiting \(W,T\) have degrees eight and five and no common
critical point, by the exact derivative gcd test in the checker.
Their degree-eight and degree-five leading coefficients are
units through the limit. The coprime-degree argument bounds
all scaled matching pairs, and the derivative test makes each
scaled gap a unit. Thus their normalized gaps in \(z\) all
have poles. Moreover \(G\) is a nonzero scalar multiple of
\(a^{19}X^8\), so \(G\) has a pole, and (5) gives a pole of
\(\rho\). Every physical gap diverges, also impossible.

A nonconstant rational function has finite preimages of every
value except possibly its value at infinity. Thus \(X\) cannot
omit both \(28/5\) and \(-2\) at finite parameters unless it is
constant. If \(X\) is constant and \(k\ne0\), equation (4) makes
\(a\) constant, and then \(G'=0\). This contradicts (5).

## The fixed shape when \(k=0\)

Now (4) gives \(X=-2\) identically. Put \(w_0=\rho a\).
After scaling \(z=a w\), the normalized generators are the
fixed shape just computed, and all their matching gaps are
fixed nonzero constants. The original gaps are these gaps
multiplied by \(w_0\). The finite-pair condition forces \(w_0\)
to have no finite zero or pole. Since \(w_0=a_0\in K\),
it is a nonzero constant. The equation \(G=g_*a^{19}\)
with \(g_*\in\mathbb C^\times\), together with (5), implies
that \(a^{20}\) is affine nonconstant in \(y\). Hence
\(\lambda^{-1}=a^{20}/w_0^{20}\) has a finite zero \(y_0\).

The unscaled pair in \(u=x-h(y)\) is now fixed. A finite
pole of \(h\) would take every original matching pair to
infinity, so \(h\) is polynomial. Let \(R_*\subset\mathbb C[u]\)
be the fixed order and \(C_*\) its fixed target image curve.
The original polynomials \(P,Q\) in \(u,y\) belong
coefficientwise to \(R_*\), since generic membership in
each finite-degree part of the fixed order is given by
constant linear equations. They therefore have polynomial
ambient lifts \(\Psi_y\) in the two fixed target generators.
The actual specialized order is contained in \(R_*\), of
colength ten, and has colength at most ten by exact conductor
specialization. Thus the two orders coincide at \(y_0\).

Choose a singular point of \(C_*\) from a fixed matching pair,
which exists by the finite-pair condition. The generic inverse
target change from the fixed pair to the original pair has
determinant \(\lambda^{-1}\). Its difference from \(\Psi_y\)
is divisible componentwise by the fixed curve equation, whose
first derivatives vanish at the singular point. Therefore
\(\det D\Psi_y\) at that point equals \(\lambda^{-1}\).
At \(y_0\) the equality of the orders makes the map of curves
an isomorphism. The induced map on their full two-dimensional
singular tangent spaces is invertible. This contradicts
\(\lambda^{-1}(y_0)=0\) by the
[singular tangent determinant argument](HYPERELLIPTIC_CONDUCTORS_AND_SINGULAR_TANGENT_DETERMINANT.md),
and proves the theorem.
