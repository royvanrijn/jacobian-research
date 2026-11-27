# The cubic–octic conductor-fourteen exclusion

**Theorem ConductorFourteenCubicOcticExclusion.** Let \(F=(P,Q)\) be a
complex polynomial plane Keller map in a frame \((x,y)\) with \(P\) monic in
\(x\) of positive degree. Its exact directional coordinate conductor cannot
have global degree fourteen and reduced characteristic sequence \((8,3)\).

We use the [exact coordinate-order
specialization](GENERIC_COORDINATE_ORDER_CONDUCTOR_AND_CLOSED_COLLISIONS.md)
and the finite-gap property from the [cubic–septic
argument](CONDUCTOR_TWELVE_HYPERELLIPTIC_AND_CUBIC_EXCLUSIONS.md#the-cubicseptic-branch).
The [characteristic-sequence
classification](CONDUCTOR_FOURTEEN_CHARACTERISTIC_SEQUENCES.md) identifies
the semigroup of this branch as \(\langle3,8\rangle\). At every finite place
of every relevant finite parameter extension at least one generic matching
pair has distinct finite original limiting coordinates. A conductor root
stays finite at such a place, its mate stays finite by monicity, and a
collapsed gap would force both longitudinal derivatives to vanish.

## Normal form and coefficient equations

Let \(K=\mathbb C(y)\). Perform triangular target reduction, centering \(u=x-h(y)\) with \(h\in K\), and monic scalings. Subtract a multiple of the square of the cubic generator, a multiple of the cubic generator, and a constant from the octic generator to remove source degrees six, three and zero. If \(\lambda\in K^\times\) is the resulting target determinant, take \(u=\rho z\), \(\rho^{11}=\lambda\), and scale the generators by \(\rho^{-3},\rho^{-8}\). Over \(K(\rho)\) one obtains a determinant-one target change with
\[
 U=z^3+az,\qquad V=z^8+bz^7+cz^5+dz^4+ez^2+fz .
\]
The target parameter identity is \(U_z V_y-U_y V_z=J\rho+\partial_z H(U,V)\) as in the cubic–septic proof, where \(J\) is the original nonzero constant Jacobian. Its left side has degree at most nine. The semigroup \(\langle3,8\rangle\) has basis degrees \(0,3,6,8,9\) through degree ten, so \(H\) is a linear combination over the parameter field of \(1,U,U^2,V,U^3\). Comparing \(z^9\) first gives \(b'=0\). The coefficients of \(U^3,V,U^2,U\) are then
\[
 -8a'/9,\quad (3c'-7ba')/8,\quad
 (ac'-5ca')/6,\quad 8a^3a'/9+(ae'-2ea')/3 .
\]
The equations at \(z^6,z^4,z^3\) integrate, respectively, to
\[
\begin{aligned}
d&=D-\frac{28a^2}{9}-\frac{49ab^2}{24}+\frac{7bc}{8},\\
e&=E+\frac{4Da}{3}-\frac{176a^3}{81}-\frac{49a^2b^2}{48}
       -\frac{7abc}{24}+\frac{5c^2}{16},\\
f&=F-\frac{7Dab}{6}+\frac{Dc}{2}+\frac{98a^3b}{81}
       +\frac{343a^2b^3}{288}-\frac{10a^2c}{9}
       -\frac{49ab^2c}{48}+\frac{7bc^2}{32},
\end{aligned}
\]
with \(D,E,F\in\mathbb C\). The \(z^1\) equation gives the constant
\[
R=\frac76Da^2b-\frac12Dac+\frac74Eab-\frac34Ec-Fa
 -\frac{28}{81}a^4b-\frac{343}{1728}a^3b^3
 +\frac5{27}a^3c-\frac{49}{192}a^2b^2c+\frac{21}{64}abc^2-\frac5{64}c^3.
\]
The degree-zero equation is \(G'=J\rho\) with
\[
\begin{aligned}
G={}&\frac4{27}Da^3-\frac{49}{96}Da^2b^2+\frac7{16}Dabc-\frac3{32}Dc^2
+\frac13Ea^2+\frac78Fab-\frac38Fc\\
&-\frac8{243}a^5+\frac{343}{1296}a^4b^2+\frac{2401}{6912}a^3b^4
-\frac{77}{216}a^3bc-\frac{343}{768}a^2b^3c
+\frac5{48}a^2c^2+\frac{49}{256}ab^2c^2-\frac7{256}bc^3 .
\end{aligned}
\]
The constants \(b,D,E,F,R\) have longitudinal weights \(1,4,6,7,9\) under scaling of \(z\) with the fixed monic degrees. \(G\) has weight ten.

## Coefficient poles and finite gaps

For a matching pair write \(S\) for its sum and \(p\) for its product in \(z\). The divided difference of \(U\) gives \(p=S^2+a\), so the squared gap is \(-3S^2-4a\). The divided difference of \(V\) gives
\[
 f_S(S)=S^7+bS^6+2aS^5+(4ab-c)S^4+(-2a^2-d)S^3
 +(3a^2b-ac)S^2+(-4a^3-2ad+e)S-a^3b+a^2c+f .
\]
At any pole of \(a\) or \(c\), take the positive scale \(q\) equal to the maximum of their pole orders divided by \(2,3\), respectively, and write \(A,C\) for their residues after multiplying by a local parameter to powers \(2q,3q\). A finite extension may make \(q\) integral. At least one of \(A,C\) is nonzero. The top term of the constant \(R\) is \(5A^3C/27-5C^3/64\), so \(A\ne0\) and either \(C=0\) or \(C^2=64A^3/27\). The top term of \(G\) is \(-8A^5/243+5A^2C^2/48\), equal to \(-8A^5/243\) or \(52A^5/243\), respectively. Thus \(G\) has pole order \(10q\).

After scaling \(S\) by the same \(q\), the monic initial polynomial is
\[
 P_{A,C}(S)=S^7+2AS^5-CS^4+\tfrac{10}9A^2S^3-ACS^2+
 (\tfrac4{81}A^3+\tfrac5{16}C^2)S-\tfrac19 A^2C .
\]
Every scaled sum is integral by this monic polynomial. In the balanced case the resultant of \(P_{A,C}\) and \(3S^2+4A\), reduced using \(C^2=64A^3/27\), is \(20800A^7/9\). In the \(C=0\) case the remainder upon division by \(3S^2+4A\) is \(-20A^3 S/81\), with no common root because \(A\ne0\). Every normalized gap has pole order \(q\). At a finite place, the derivative of a pole increases the pole order by the ramification index. Hence \(\rho=G'/J\) also has a pole and *all* physical gaps diverge. This contradicts the finite-gap property, so all normalized coefficients and \(G\) are regular at every finite place. The monic \(f_S\) then bounds all normalized gaps at finite places, making any finite zero of \(\rho\) impossible.

A pole of \(\rho\) at an unramified finite place is impossible by regularity
of \(G\). The Kummer extension \(K(\rho)/K\) has group contained in the
eleventh roots of unity, since the ground field contains them and
\(\rho^{11}=\lambda\). Its degree divides the prime eleven. The only
ramification index other than one is therefore eleven. At such a place
inertia sends \(\rho\) to \(\zeta\rho\) for every eleventh root of unity
\(\zeta\). The normalized coefficients transform by the corresponding
inverse longitudinal weights. None of the weights \(1,4,6,7,9\) is divisible
by eleven, so the corresponding constants \(b,D,E,F,R\) must all vanish.
Therefore for nonhomogeneous constants, \(\rho^{11}\) has no finite poles or
zeros and \(\rho\in\mathbb C^\times\). Then \(a,c\) and \(G\) are
polynomials in \(y\). If \(a,c\) are both constant, \(G'=0\). If one has a
pole at infinity, the same leading-weight computation gives pole order
\(10q>1\) for \(G\) (and forces \(a\) to have a pole). Neither case permits
\(G'=J\rho\in\mathbb C^\times\).

## Homogeneous constants and the fixed singular order

Now suppose \(b=D=E=F=R=0\). Then \(c(64a^3-27c^2)=0\) exactly and \(a\ne0\)
by \(G'=J\rho\). If \(c=0\), take \(\tau^2=a\). If \(c\ne0\), take
\(\tau=c/a\); then \(a=(27/64)\tau^2\) and \(c=(27/64)\tau^3\). After
scaling by \(\tau\), \(U,V\) are fixed monic polynomials, and the
corresponding matching-pair polynomial is precisely \(P_{A,C}\) above, with
fixed \(A,C\). These pairs are finite in number since the rational-function
field index divides both degrees three and eight. Their gaps are nonzero by
the resultant and remainder calculation.

All original gaps are therefore fixed nonzero multiples of \(w=\rho\tau\).
The finite-gap property shows that \(w\) is a unit over every finite
parameter. In the first case \(w^2\in K\) is the linear coefficient of the
unscaled cubic generator. In the second \(w\in K\) is the ratio of the
degree-five coefficient of the unscaled octic generator to that same linear
coefficient. Thus \(w\in\mathbb C^\times\). The formula for \(G\) is
\(G=g_0\tau^{10}\) for \(g_0\in\mathbb C^\times\), by the leading forms
already checked. Since \(G'=Jw/\tau\), \(\tau^{11}\) is nonconstant affine
in \(y\). Hence \(\lambda^{-1}=(\tau/w)^{11}\) has a finite zero.

The unscaled reduced generators \(U_0,V_0\) are constant polynomials in the
centered coordinate \(u=x-h(y)\). The center \(h\) is polynomial, because
all original pairs are \(h(y)\) plus fixed offsets and at least one pair
stays finite at every finite parameter. In the centered coordinates,
\(P,Q\in\mathbb C[y,u]\) have all their parameter coefficients in the fixed
order \(R_0=\mathbb C[U_0,V_0]\): membership in each finite-degree part of
\(R_0\) is a system of constant-coefficient linear equations that holds
identically over \(K\). Thus polynomial-in-\(y\) target lifts
\(\Psi_y(T,W)\) of \(P,Q\) exist.

The exact conductor polynomial \(\kappa(u+h(y),y)\) equals
\(A_0(y)\kappa_0(u)\), where \(\kappa_0\) is the monic generator of the
conductor of \(R_0\). This equality holds generically because the generic
order is fixed, and comparison of leading coefficients shows \(A_0\in\mathbb
C[y]\). No finite zero of \(A_0\) is possible since the exact specialized
conductor is nowhere the zero polynomial. In particular the colength of the
actual specialized order is constantly seven. By the inclusion in \(R_0\)
and its colength seven, \(\Psi_y\) generates \(R_0\) at every finite
parameter.

Apply the [finite singular tangent determinant
argument](HYPERELLIPTIC_CONDUCTORS_AND_SINGULAR_TANGENT_DETERMINANT.md#the-determinant-at-a-finite-singular-point)
to the fixed order. Take a fixed singular point of the plane curve given by
\(R_0\). Such a point exists because its colength relative to the
normalization is positive. The generic ambient target inverse of the
reduction from \(P,Q\) to \(U_0,V_0\) has determinant \(\lambda^{-1}\). Its
difference from \(\Psi_y\) lies coordinatewise in the principal ideal of the
fixed reduced relation. At the singular point the relation and its first
derivatives vanish, so the Jacobian determinants of these lifts agree. The
Jacobian determinant of \(\Psi_y\) at this point is polynomial in \(y\). At
every specialization the two specialized components generate \(R_0\), so
their differences from their values generate the maximal ideal at the fixed
point. They span its two-dimensional cotangent space, and the specialized
determinant is nonzero. This contradicts the finite zero of
\(\lambda^{-1}\). Thus the homogeneous case is also impossible.

## Exact controls

The [checker](verify_conductor_fourteen_cubic_octic_exclusion.py)
replays the coefficient integrations, the collision equation, both
highest-weight gap tests and the homogeneous scaling. The exclusion
at finite places and the specialized tangent argument are written
proofs. Replay from the repository root.

~~~sh
.venv/bin/python research/hessian/hc4-triangular-graphs/verify_conductor_fourteen_cubic_octic_exclusion.py
~~~
