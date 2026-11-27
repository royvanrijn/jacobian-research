# Exclusion of a simple pole in the singular alternating elliptic slice

**Theorem SixSheetSingularAlternatingEqualPolarOrders.** In the
six-sheet \(A_6\) profile with one unramified omitted prime, suppose
a target coordinate \(p\) has width four. Then both poles of the
companion target coordinate \(q\) on the generic compactified slice
have order three.

## The markings and the differential equation

Use the
[singular alternating elliptic factor](KELLER_TWO_POLE_ELLIPTIC_SLICE_MARKED_FACTORS.md#the-singular-alternating-profile-at-width-four)
and the relatively minimal model \(S_0\to\mathbb A^1_p\). On the
generic compact fibre \(dq=\varphi^2\omega\), where \(\omega\)
is a nowhere-zero relative differential on \(S_0\) and \(\varphi\)
is a constant multiple of the reciprocal of that polynomial
factor. The constant has a square root over \(\mathbb C\).
The pole orders of \(q\) are \((1,5)\) or \((3,3)\).

Assume the orders are \((1,5)\). On the resolved special fibre,
one pole is on its rational end and the other is on the
surviving integral arithmetic-genus-one component. The latter
is nonconstant under \(q\) by finiteness on the original
normalization. It cannot have only a simple pole, since a
degree-one finite map from an integral curve to the normal
projective line would be an isomorphism. Thus the pole section
through the rational end is the simple one.

Write \(O\) for the unblown high pole section and \(T\) for the
simple pole section on \(S_0\). The
[two contractions](KELLER_SIX_SHEET_CUSP_SLICE_GENUS.md#singular-normalization-with-one-ordinary-affine-cusp)
move the mark of the rational end first to the exceptional
curve and then to the surviving component, away from \(O\).
The reverse blowups are successively centered on \(T\) and its
strict transform over the cusp projection value \(p_0\).
Apply the
[Weierstrass classification](KELLER_SIX_SHEET_CUSP_SLICE_GENUS.md#the-two-weierstrass-models-of-the-elliptic-residuals)
with zero section \(O\). Its two models are
\[
 Y^2=X^3+aX+b,\qquad (a,b)=(0,p^2-1)
                  \quad\hbox{or}\quad (-3(p^2-1),2p(p^2-1)).
\tag{1}
\]
Here \(T=(u,v)\) with \(u,v\in\mathbb C[p]\), because the
sections remain disjoint over the finite base and the finite
Weierstrass surface is smooth. The differential \(dX/(2Y)\)
and \(\omega\) differ by a nonzero constant since both
trivialize the relative dualizing bundle over the finite base,
where the only units are constants.

On the original surface the leading polar coefficient of \(q\)
nowhere vanishes along each pole section, by the local
power-map description. Near \(T(p_0)\) on \(S_0\) use
\(k=p-p_0\) and a normal coordinate \(s\) to \(T\), and
write \(q=N(k,s)/s\), with \(N\) regular. The two reverse
blowups have coordinates \(s=k s_1=k^2s_2\) along the
successive strict transforms of \(T\). The leading polar
coefficient on the original section is \(N(k,0)/k^2\), a
unit. Thus the polar coefficient of \(q\) on \(T\) has a
double zero at \(p_0\) and no other finite zero. On \(O\)
the corresponding coefficient is nowhere zero.

Put \(f_T=(Y+v)/(X-u)\). It has simple poles at \(O,T\),
and \(\operatorname{Res}_T(f_T\,dX/(2Y))=1\). The spaces
with pole divisors \(3O+T,\ 5O+T\) have bases
\(1,X,Y,f_T\) and \(1,X,Y,X^2,XY,f_T\), respectively,
by Riemann–Roch and these pole orders. The rational
function \(\varphi\) has pole divisor \(3O+T\) on the
generic fibre and no vertical pole on \(S_0\).
Indeed on the original normalized affine chart its only
possible vertical pole is on the contracted rational end.
Its leading coefficient at \(O\), in the formal parameter
\(-X/Y\), is a unit over the finite base by the differential
equation. Thus its coefficient of \(Y\) is a nonzero
constant. Normalize it to one, scaling \(q\) by a constant
to keep the equation:
\[
         \varphi=Y+BX+C f_T+L,\qquad dq=\varphi^2 dX/(2Y).
\tag{2}
\]
The coefficient \(B\) is polynomial in \(p\) by extracting
the order-two coefficient at \(O\) after subtracting \(Y\).
The coefficient \(C\) is the residue of \(\varphi\,dX/(2Y)\)
at \(T\). It is regular over the finite base, and (2) and the
double zero of the polar coefficient of \(q\) give
\(C=h(p-p_0)\), \(h\ne0\).

Set \(D_\omega=2Y\partial_X+(3X^2+a)\partial_Y\).
The coefficient of \(f_T\) in \(q\) is \(-C^2\) by the
pole at \(T\). Matching the highest orders at \(O\) and
the terms linear in \(Y\) gives, up to a function of \(p\),
\[
 \begin{split}
 q&=XY/5+(B/2)X^2+(B^2+2C)Y/3+(BC+L)X-C^2 f_T,\\
 L&=-Bu-v .
 \end{split}
\]
The remaining equations from \(D_\omega q=\varphi^2\) are
\[
\begin{split}
 uB^2+vB&=C^2+Cu+a/5, \tag{3}\\
 0&=5aB^2-15B^2u^2-30BCv-30Buv-15C^2u
                  -20Ca-30Cu^2+9au+9u^3-24v^2 . \tag{4}
\end{split}
\]
The [exact checker](verify_keller_singular_alternating_simple_pole.py)
replays these identities, with \(v^2=u^3+au+b\).

## The \(IV^*\) model

Here \(a=0,\ b=p^2-1\). If \(u\) is nonconstant,
\(\deg u=2m,\deg v=3m\). For \(D=p^2-1\), differentiation
of \(v^2-u^3=D\) gives \(3u'D-uD'=v(3u'v-2uv')\).
If \(\deg u>2\), the left degree is below \(\deg v\),
so the left side vanishes. This would make \(u^3/D\)
constant, impossible by degree. Thus \(\deg u\le2\).

If \(u\) is constant, the curve equation gives \(u^3=1\),
\(v=\pm p\). Scale \(X\) by a cube root of unity and
change the sign of \(Y\) to take \(u=1,v=p\).
Then (3) bounds \(\deg B\le1\). The degree-three term in
(4) comes only from \(-30BCv\), forcing \(B\) constant,
and (3) contradicts \(\deg C=1\).

If \(\deg u=2\), center the quadratic by translating the
parameter temporarily. The coefficients of degrees five,
four and three in \(v^2-u^3\) force \(v\) to be odd in
the centered parameter. The right side \(p^2-1\) is even
only if the center is zero. Comparing the quadratic and
constant terms and scaling \(X\) by a cube root of unity gives
\[
                 u=1+A^2p^2,\quad
                 v=A p(3/2+A^2p^2),\quad A^2=-4/3 .
\]
By (3), \(B\) has degree at most one. For \(B\) constant
the degree-six term in (4) comes only from \(9u^3-24v^2\)
and is nonzero. Thus (3) in degree four gives \(B=-Ap+b_0\).
The coefficients of degrees three, two and one in (3) give
\[
     A b_0+h=0,\qquad 4b_0^2+3h^2+4hp_0-2=0,\qquad
     -A b_0+4h^2p_0-2h=0 .
\]
The first two give \(h p_0=1/2\), and the third then gives
\(h=0\). This excludes the \(IV^*\) model.

## The \(I_2^*\) model

Set \(K=v+uB\), \(D_0=B^2-u-2C\).
Using (3), (4) and the equation of \(T\) gives
\[
\begin{aligned}
 BK&=C^2+Cu+a/5,\\
 u^2D_0&=2uC^2+\tfrac75 au+b-K^2,\\
 0&=5aD_0+15u^2D_0+30CuD_0
                    -(15C^2+16a)u-30C^3-16aC-24b .
\end{aligned}
\tag{5}
\]
If \(\deg u\ge2\), the curve equation gives
\(\deg u=2m,\deg v=3m\). By (3), \(\deg B=m\), except
possibly if \(m=1,\deg B=0\). That exception makes the
degree-six coefficient of (4) nonzero from \(9u^3-24v^2\).
The leading coefficient of \(B\) is \(-v_{3m}/u_{2m}\),
whose square is \(u_{2m}\). Also \(\deg K=m+1\) by (5).
For \(m\ge2\) the second equation forces \(D_0=0\) by
degree. The third would then give \(15C^2+16a=0\), since
the other terms have degree at most three. This is
impossible because \(a=-3(p^2-1)\) is squarefree.
For \(m=1\), the second equation makes \(D_0\) constant.
The coefficients of \(p^4\) in the second and third
equations give, respectively,
\[
             u_2D_0=h^2+\tfrac75a_2,\qquad
             u_2D_0=h^2+\tfrac{16}{15}a_2 .
\]
Here the first comparison uses \(K_2^2=u_2h^2\), since
\(B_1K_2=u_2h\) and \(B_1^2=u_2\). The two equalities
contradict \(a_2=-3\).

It remains that \(\deg u=1\), since constant \(u\) makes
\(u^3+au+b\) cubic and not a square. For slope \(u_1\) the
cubic term must cancel, giving \((u_1-1)^2(u_1+2)=0\)
and \(\deg v\le1\). Equation (3) now makes \(B\) constant.
Its degree-two coefficient is \(h^2+u_1h-3/5=0\).
For slope \(1\), the cubic coefficient in (4) gives
\(-15h^2+30h-18=0\), which forces \(h=3/5\) after using
(3). This value does not satisfy (3).

For slope \(-2\), put \(u=-2p+c,\ v=d p+e\).
The curve equation gives \(c^2=16/3,\ d^2=9c,\ de=-20\):
indeed \(v^2=9cp^2-(6c^2+8)p+c^3+3c\), which has
degree one if \(c=0\) and must have zero discriminant
otherwise. Set \(B=d z\), \(C=h p+c w\).
The coefficients of degrees two through zero of (3)
and the degree-two coefficient of (4) give, respectively,
\[
\begin{split}
 5h^2-10h-3&=0,\\
 -2hw-h+2w-18z^2+9z&=0,\\
 -80w^2-80w+720z^2-300z-9&=0,\\
 h^2-4hw+18hz-8h+4w+45z^2-36z+9&=0 .
\end{split}
\]
These generate the unit ideal over \(\mathbb Q\), as
the exact checker verifies. Thus there is no simple
pole in the \(I_2^*\) case either.
