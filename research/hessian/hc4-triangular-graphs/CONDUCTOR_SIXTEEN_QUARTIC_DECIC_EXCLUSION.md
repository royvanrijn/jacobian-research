# The quartic–decic sequence \((10,4,9)\) in conductor degree sixteen

**Theorem (The \((10,4,9)\) conductor-sixteen exclusion).** Let
\((P,Q)\) be a complex polynomial plane Keller map in a source frame
\((x,y)\) with \(P\) monic of positive degree in \(x\). If its exact
coordinate-order conductor has global \(x\)-degree sixteen, its reduced
generic characteristic sequence cannot be \((10,4,9)\). The proof
allows degree drops on individual source lines.

## 1. Laurent equations for the actual generators

Write \(K=\mathbb C(y)\) and \(J=\operatorname{Jac}_{x,y}(P,Q)\ne0\).
The
[coordinate-order theorem](GENERIC_COORDINATE_ORDER_CONDUCTOR_AND_CLOSED_COLLISIONS.md)
gives the finite birational orders on the generic and all specialized
lines. The generic order has colength eight. Apply the
[Assi–García-Sánchez approximate-root theorem](sources/assi-garcia-sanchez-1407.0490v1.pdf),
Propositions 2 and 4, to the reduced pair. Its one-place hypothesis and
monic relation conditions are those verified in the
[cubic–quintic proof](CONDUCTOR_EIGHT_CUBIC_QUINTIC_EXCLUSION.md).

Normalize the degree-four and degree-ten generators \(U,V\) to be
monic. The degree-nine approximate root is monic of degree five in
the quartic generator. For lower exponents \(i\) of \(U\) the weights
\(4i\), \(0\le i<5\), are distinct modulo ten. Thus no terms above
weight twenty can cancel. Completing the square in \(V\) and
translating \(U\) gives an approximate root of the form
\[
 W=V^2-(U^5+\alpha U^3+\beta U^2+\gamma U+\delta),\qquad \deg W=9.
\tag{1}
\]
This uses only polynomial target changes over \(K\), as for
[sequence \((10,4,7)\)](CONDUCTOR_FOURTEEN_QUARTIC_DECIC_EXCLUSION.md).
After depressing the quartic by a source translation \(h(y)\in K\),
let \(\lambda\in K^*\) be the determinant of the target change.
Choose \(\rho^{14}=\lambda\), write \(x=h(y)+\rho z\), and rescale
the generators by \(\rho^{-4},\rho^{-10}\). The new target change
has determinant one. In these scaled coordinates continue to denote
the generators and coefficients by \(U,V,\alpha,\beta,\gamma,\delta\).
Then
\[
 U=z^4+az^2+bz+c,\qquad \Delta=a^2-4c.
\]
Let \(S\) be the monic formal square root of the polynomial in \(U\)
in (1) at \(z=\infty\). Since \(V+S\) has degree ten, the leading
term of \(V-S\) is a nonzero multiple of \(z^{-1}\), with
\[
 V=[S]_{\ge0},\qquad E=V-S=\sum_{i\ge1}e_i z^{-i},\qquad e_1\ne0.
\]
The parameter-dependent target automorphism, now of determinant one,
gives the
[Hamiltonian identity](CONDUCTOR_TWELVE_CABLE_EXCLUSIONS.md#1-common-source-and-target-identities)
\[
 U_z V_y-U_y V_z=J\rho+\partial_z H(U,V).
\tag{2}
\]
The polynomial part of the contribution from \(S\) on the
left is \([U_z(\alpha' U^{1/2}+\beta' U^{-1/2})/2]_{\ge0}\);
the contribution from \(E\) is
\[
 4e_1'z^2+4e_2'z+4e_3'+2ae_1'+a'e_1 .
\]
The degree is at most five. Since the semigroup is
\(\langle4,9,10\rangle\), the evaluated \(H\) through degree six
belongs to the span of \(1,U\).
The degree-five coefficient gives \(\alpha'=0\), and then the
degree-three coefficient makes \(H\) constant.
Since \([U_z U^{-1/2}]_{\ge0}=4z\), comparison of degrees two,
one and zero gives
\[
 e_1'=0,\qquad 2e_2'+\beta'=0,\qquad (4e_3+ae_1)'=J\rho .
\]

Set \(\kappa=-e_1\in\mathbb C^*\). The expansion of
\(S=U^{5/2}+\alpha U^{1/2}/2+\beta U^{-1/2}/2+O(z^{-6})\)
gives
\[
 \begin{split}
 \kappa&=b(15\Delta^2-40ab^2+64\alpha)/256,\\
 -e_2-\beta/2&=5(-\Delta^3+24\Delta ab^2-8b^4)/1024-\alpha\Delta/16,\\
 e_3&=a\kappa/2-5\Delta b^3/128 .
 \end{split}
\]
In particular \(b\ne0\), and the two constant parameters and the
remaining derivative equation can be written
\[
 \begin{split}
 a&=(15\Delta^2+64\alpha-256\kappa/b)/(40b^2),\\
 C&=\Delta^3+\tfrac{16}5\alpha\Delta-\tfrac{96}5\kappa\Delta/b-b^4
       \ \in\mathbb C,\\
 G&=a\kappa-5\Delta b^3/32,\qquad G'=J\rho.
 \end{split}
\tag{4}
\]
The nonzero constant \(\kappa\) has scaling weight eleven.
Its unscaled value is in \(K\), so \(\rho^{11}\in K\). Together
with \(\rho^{14}\in K\) this gives \(\rho\in K\), and hence the
coefficients in (4) are all rational in \(y\).

## 2. Matched-pair equations

For distinct paired points on a generic line in the \(z\)-coordinate,
write \(s\) for their sum and \(g\) for their gap. The quartic divided
difference gives \(s\ne0\) since \(b\ne0\), and their product is
\((s^3+as+b)/(2s)\). The decic divided difference, after inserting
that product and (4), equals \(s/32\) times the polynomial
\[
 \begin{split}
 F(s)={}&2s^4(s^2+2a)^2+5\Delta s^2(s^2+2a)-2bs^3(s^2+2a)\\
        &\quad -2b^2s^2-5\Delta b s+64\kappa/b .
 \end{split}
\tag{5}
\]
For a paired sum \(F(s)=0\), and \(g^2=-s^2-2a-2b/s\).
Equivalently, with \(Y=-g^2=s^2+2a+2b/s\),
\[
 2s^4Y^2+(5\Delta s^2-10bs^3)Y
           +10b^2s^2-15\Delta bs+64\kappa/b=0 .
\tag{6}
\]
The proof needs only necessity of these equations. The
[exact checker](verify_conductor_sixteen_quartic_decic_exclusion.py)
verifies the Laurent coefficients, the differentiated parameter identities,
and (5)–(6).

At each finite parameter there is a generic matched pair with a finite
distinct pair of limits in the original \(x\)-coordinate. This is
the finite-unit-gap property used in the
[lower cable exclusions](CONDUCTOR_TWELVE_CABLE_EXCLUSIONS.md#1-common-source-and-target-identities).
The exact specialized conductor is nonzero and nonconstant, so some
generic conductor root reaches a finite specialized root. Its mate
stays finite by monicity of the original \(P\). If their limits
coincided, the two original divided differences would give
\(P_x=Q_x=0\) there. Thus some physical gap \(\rho g\) is a unit at each
finite parameter, after passage to the endpoint field.

## 3. A finite pole of \(b\) is impossible

Suppose \(b\) has a finite pole. Equation (4) then forces
\(\Delta^3\sim b^4\), where \(\sim\) denotes equality of nonzero
leading terms. Hence the pole orders of \(b,\Delta,a\) are
\(3k,4k,2k\) for \(k>0\), and
\(a\sim3\Delta^2/(8b^2)\). The pole order of \(G\) is \(13k\),
so \(\rho=G'/J\) has pole order \(13k+1\).

Put \(\sigma=\Delta/b\) and \(z=\sigma w\). The scaled polynomials
\(\sigma^{-4}U(\sigma w)\) and \(\sigma^{-10}V(\sigma w)\)
have bounded coefficients and leading limits
\[
 U_*=w^4+\tfrac38 w^2+w-\tfrac{55}{256},\qquad V_*=[U_*^{5/2}]_{\ge0}.
\]
Indeed \(b/\sigma^3,\Delta/\sigma^4\to1\), and the constant \(\alpha\)
term of \(V\) vanishes in the scale. The limiting derivative resultant
\(\operatorname{Res}_w(U_*',V_*')=208000\) is checked in the exact replay.

If an original gap is a unit, its \(w\)-gap tends to zero since
\(\rho\sigma\) has a pole. Its two \(w\)-endpoints must be bounded:
if an endpoint had a pole, the other would have the same pole and
leading term, since their difference tends to zero. The leading quartic
divided difference would then be four times the cube of that leading
term, which cannot vanish. The two bounded \(w\)-endpoints therefore
have a common limit, and their divided differences give a common
critical point of the limiting \(U_*,V_*\), contrary to the resultant.
This excludes a finite pole of \(b\).

## 4. A finite zero of \(b\) is impossible

Suppose \(v(b)=l>0\) at a finite parameter. If \(\Delta\) is
bounded, (4) gives \(\Delta=O(b)\) and \(a\sim-32\kappa/(5b^3)\).
If \(\Delta\) has a pole, the same equation gives
\(\Delta^2\sim96\kappa/(5b)\) and \(a\sim4\kappa/(5b^3)\).
In either case \(G\sim a\kappa\), so \(v(\rho)=-3l-1\).
We show that every physical gap diverges, including the pairs whose
normalized gaps get small.

For the sum equation (5), the nonzero endpoint coefficients of the
lower Newton polygon have valuations \(-l,-6l,0\) at degrees zero,
four and eight. The degree-two coefficient has valuation at least
\(-7l/2\), and the degree-six one has valuation \(-3l\);
all other coefficients lie above the two lower segments.
Thus four sums have valuation \(5l/4\), and four have valuation
\(-3l/2\), counted with multiplicity after extension to the endpoint
field. For the first group \(g^2\sim-2a\), so the physical gap diverges.

In the second group the face is \(2s^4(s^2+2a)^2\), so
\(s^2\sim-2a\). Use (6). The coefficients of \(Y^2,Y,1\)
have valuations respectively \(-6l\), at least \(-7l/2\), and
at least \(-l\). The constant term has valuation exactly \(-l\):
if \(\Delta\) is bounded, it has leading coefficient \(192\kappa\)
relative to \(b^{-1}\), using \(a\sim-32\kappa/(5b^3)\).
If \(\Delta\) has a pole, write the leading terms as
\(\Delta\sim d b^{-1/2}\) and \(s\sim q b^{-3/2}\) over a local
extension. Then \(d^2=96\kappa/5\), \(q^2=-8\kappa/5\), and
the constant term has leading coefficient \(48\kappa-15dq\)
relative to \(b^{-1}\). Its vanishing would force
\(dq=16\kappa/5\), impossible on squaring because
\(d^2q^2=-768\kappa^2/25\).
The quadratic Newton polygon of (6) therefore gives
\(v(Y)=5l/2\), even if its leading roots are repeated.
Thus \(v(g)=5l/4\), while
\[
 v(\rho g)=5l/4-3l-1<0.
\]
This excludes a unit physical gap also in the second group,
and so excludes every finite zero of \(b\).

Since \(b\in K\) is nonzero and has neither a finite zero nor a finite
pole, it is constant. The cubic relation (4) then makes \(\Delta\)
constant, and its first equation makes \(a\) constant. Now \(G'=0\)
contradicts \(J\rho\ne0\). This excludes the characteristic
sequence \((10,4,9)\).
