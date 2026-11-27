# The quintic–sextic conductor-twenty exclusion

**Theorem ConductorTwentyQuinticSexticExclusion.** Let \(F=(P,Q)\)
be a complex polynomial plane Keller map with \(P\) monic of
positive degree in the longitudinal coordinate \(x\). Its generic
directional order cannot have reduced characteristic sequence
\((6,5)\).

Put \(K=\mathbb C(y)\) and \(J=\operatorname{Jac}_{x,y}(P,Q)\ne0\).
The [exact coordinate-order theorem](GENERIC_COORDINATE_ORDER_CONDUCTOR_AND_CLOSED_COLLISIONS.md)
identifies the generic colength as ten and the restricted conductor
on every finite line. In particular, at each finite parameter some
generic matching pair has distinct finite limiting physical endpoints
in \(x\). A generic conductor-root branch approaches a finite root
of the nonconstant specialized conductor, and monicity of \(P\)
keeps its mate finite. Equal limits would force \(P_x=Q_x=0\)
by their divided differences. At each finite place some physical
gap is therefore a unit.

## The parameter equations

Reduce the actual generators over \(K\) by triangular target changes.
Take monic degree-five and degree-six generators, depress the
quintic by \(u=x-h(y)\), and remove its constant. Subtract a
multiple of the quintic and a constant from the sextic to remove
its coefficients in degrees five and zero. If \(\lambda\in K^\times\)
is the determinant of this target change, take \(\rho^{11}=\lambda\),
set \(u=\rho z\), and multiply the two generators by \(\rho^{-5}\)
and \(\rho^{-6}\). The normalized pair has determinant one in the
target change and takes the form
\[
 U=z^5+az^3+bz^2+cz,\qquad
 V=z^6+dz^4+ez^3+fz^2+gz .
\]
The [parameter identity](CONDUCTOR_EIGHTEEN_QUARTIC_SEPTIC_EXCLUSION.md#the-parameter-identities)
gives
\[
              U_zV_y-U_yV_z=J\rho+\partial_z H(U,V)       \tag{1}
\]
for a target polynomial \(H\). The left side has degree at most
eight, and the only pole degrees through nine in \(\langle5,6\rangle\)
are \(0,5,6\). Thus \(H\) lies in the span of \(1,U,V\).

The coefficients in degrees eight, seven and six of (1) give
\[
 d=6a/5+D,\qquad e=6b/5+E,\qquad f=4Da/5+3a^2/25+6c/5+F,
 \qquad D,E,F\in\mathbb C .
 \tag{2}
\]
In degrees five and four one obtains the coefficients of \(V,U\)
in \(H\):
\[
\begin{aligned}
 \theta_V&=(25g'-20Db'-15Ea'-6ab'-6ba')/30,\\
 \theta_U&=(20Daa'-100Dc'-75Eb'-50Fa'
                   +12a^2a'-30ac'-30ca'-30bb')/125 .
\end{aligned}
\]
To state the remaining first integrals, use the polynomial differential
\(\mathrm d\) in \(a,b,c,g\), treating \(D,E,F\) as constants.
Replace the primes in the two \(\theta\)'s by this differential.
Define \(\omega_j\) by
\[
 U_z\,\mathrm dV-V_z\,\mathrm dU-\theta_U U_z-\theta_V V_z
              =\omega_3 z^3+\omega_2 z^2+\omega_1 z+\omega_0 .
\]
For a polynomial one-form \(\omega\), let
\[
 {\cal I}(\omega)=\int_0^1\omega_{(ta,tb,tc,tg)}(a,b,c,g)\,\mathrm dt .
\]
The integral of each power of \(t\) defines a polynomial over the
coefficient field. Put
\[
 \begin{aligned}
 I_7&={\cal I}(\omega_3),&
 I_8&={\cal I}(\omega_2),\\
 I_9&={\cal I}(\omega_1+2I_7\,\mathrm da/5),&
 G_{10}&={\cal I}(\omega_0+(I_7\,\mathrm db+I_8\,\mathrm da)/5).
 \end{aligned}
\]
All four forms being integrated are exact with the displayed
primitives, as the
[exact checker](verify_conductor_twenty_quintic_sextic_exclusion.py)
verifies. Thus (1) gives the constants
\[
  R_7=I_7,\qquad R_8=I_8,\qquad R_9=I_9-2R_7a/5
\]
and the parameter equation
\[
           \Gamma=G_{10}-(R_7 b+R_8 a)/5,\qquad
           \Gamma'=J\rho .                                \tag{3}
\]
The longitudinal weights of \(a,b,c,g\) are \(2,3,4,5\).
The weights of \(D,E,F,R_7,R_8,R_9\) are \(2,3,4,7,8,9\),
and \(\Gamma\) has weight ten.

## The leading shapes

Setting the constants to zero in the three first integrals
gives the following leading equations without their nonzero denominators:
\[
\begin{aligned}
 18a^2b-25ag-30bc&=0,\\
 9a^4-60a^2c-90ab^2+125bg+75c^2&=0,\\
 6a^3b-25a^2g+60abc+20b^3-125cg&=0 .
\end{aligned}                                             \tag{4}
\]
The corresponding polynomial for the leading part of \(\Gamma\) is
\[
               \Gamma_*=
             -(108a^5-2700a^2b^2-4500ac^2+15625g^2)/37500 .
\]
If \(a=0\), equations (4) force \(b=c=0\). In a nonzero shape
one then has \(g\ne0\). The pair \((z^5,z^6+gz)\) has no
common longitudinal critical point, and
\(\Gamma_*=-5g^2/12\ne0\).

If \(a\ne0\), make \(a=1\) by a constant scale. The first
equation gives \(g=b(18-30c)/25\). The remaining equations give
either \(b=0,\ 75c^2-60c+9=0\), or
\[
               b^2=3/5-15c^2/2,\qquad
               375c^3+25c^2-50c+3=0 .
\]
In both cases the limiting pair is
\[
 U_*=z^5+z^3+bz^2+cz,\qquad
 V_*=z^6+\tfrac65 z^4+\tfrac65 b z^3+
                       (\tfrac3{25}+\tfrac65 c)z^2+gz .
\]
The checker verifies the reductions of \(\Gamma_*\) and the
longitudinal derivative resultant:
\[
\begin{array}{c|c|c}
 &\Gamma_*& \operatorname{Res}_z(U_*',V_*')\\ \hline
 b=0 &6(50c-9)/3125&1296(365c-72)/3125\\
 b\ne0 &91(125c^2-25c-3)/6250&
                  -2548(24200c^2-7555c+132)/3125 .
\end{array}
\]
Each reduced expression is coprime to the corresponding
quadratic or cubic in \(c\). Thus every nonzero leading
shape is immersed and has nonzero \(\Gamma_*\).

## Finite places and infinity

Suppose \(a,b,c\), or \(g\) has a pole at a finite place.
Let \(q>0\) be the largest pole order divided by its weight
\(2,3,4,5\), respectively. Rescale the source coordinate by
a factor of pole order \(q\), and the two monic generators
by its inverse fifth and sixth powers. After this scaling the
coefficients are regular and the nonzero leading shape obeys
(4). The constants of positive weight vanish under the
scaling. The nonzero leading \(\Gamma_*\) gives a pole of order
\(10q\) for \(\Gamma\), so \(\rho\) has a pole by (3).

All matching pairs stay bounded in the rescaled coordinate,
since the degrees five and six are coprime and both leading
coefficients remain units. If a pair escaped, scaling by the
largest endpoint would give a nonzero common zero of the divided
differences of the two monic powers. Coprimality rules this out
off the diagonal, and the derivatives rule it out on the nonzero
diagonal. Limiting immersion prevents bounded endpoints from
coalescing. Thus every normalized \(z\)-gap has pole order \(q\),
and multiplying by \(\rho\) makes every physical gap diverge.
This contradicts the finite-pair condition.

So all coefficients and \(\Gamma\) are regular at finite places.
The normalized pairs are bounded there by the same coprime-degree
test. A finite zero of \(\rho\) would shrink every gap, which is
impossible. If any of \(D,E,F,R_7,R_8,R_9\) is nonzero, its weight
is coprime to eleven. Its unscaled counterpart lies in \(K\),
so its nonzero constant scaled value gives the corresponding
power of \(\rho\) in \(K\). Since \(\rho^{11}\in K\), this gives
\(\rho\in K\). Regularity of \(\Gamma\) then excludes finite
poles of \(\rho\) by (3). Thus \(\rho\) is constant and
\(a,b,c,g\) are polynomials in \(y\). If one has a pole at
infinity, the same leading equations give \(\deg_y\Gamma=10q>1\),
since \(q\ge1/5\). If all are constant, \(\Gamma'=0\).
Both contradict \(\Gamma'=J\rho\in\mathbb C^\times\).

## The fixed order in the homogeneous case

Suppose instead all six constants vanish. Equations (4) now
hold exactly. Not all \(a,b,c,g\) can vanish by (3). If \(a=0\),
scale by \(\tau\) with \(\tau^5=g\). If \(a\ne0\), take
\(\tau^2=a\). The rescaled coefficients are fixed complex
constants by the finite shape list above. The matching pairs
in \(z/\tau\) are fixed, and finite in number since the field
index divides both five and six. All physical gaps are their
nonzero gaps multiplied by \(w=\rho\tau\). Hence \(w\)
has neither a finite zero nor a finite pole. In the first
case \(w^5\in K\) from the unscaled coefficient \(g\), and in
the second case \(w^2\in K\) from the unscaled coefficient \(a\).
Thus \(w\) is a nonzero constant. Also
\(\Gamma=\gamma_0\tau^{10}\) with \(\gamma_0\ne0\), and
\(\Gamma'=Jw/\tau\) makes \(\tau^{11}\) affine nonconstant in
\(y\). Consequently \(\lambda^{-1}=\tau^{11}/w^{11}\) has a
finite zero \(y_0\).

The unscaled normalized pair in \(u=x-h(y)\) is fixed. Since
its matching points are fixed, a finite pole of \(h\) would
take every physical matching pair to infinity. Thus \(h\) is
polynomial. Write \(R_*\subset\mathbb C[u]\) for the fixed
order generated by the two unscaled polynomials. Membership
in each finite-degree part of \(R_*\) is given by constant
linear equations. Generic membership of \(P,Q\) therefore
gives membership of all their parameter coefficients, so the
original pair admits polynomial ambient lifts \(\Psi_y\) in
the fixed target generators. The specialized order is
contained in \(R_*\), of colength ten, and has colength at most
ten by the exact specialized conductor. It equals \(R_*\) at
\(y_0\).

Choose a singular point of the fixed image curve, for example
the image of a fixed matched pair. The generic inverse target
change and the regular lifts differ by multiples of the
fixed plane equation, whose first derivatives vanish at that
point. Their target determinants agree there, and thus the
determinant of \(\Psi_y\) is \(\lambda^{-1}\). At \(y_0\) the
specialized lift induces an isomorphism of image curves by
the equality of orders. On their full two-dimensional
singular tangent spaces its derivative is invertible.
The [singular tangent determinant argument](HYPERELLIPTIC_CONDUCTORS_AND_SINGULAR_TANGENT_DETERMINANT.md#the-determinant-at-a-finite-singular-point)
therefore contradicts the zero of \(\lambda^{-1}\). This
excludes the homogeneous case and proves the theorem.
