# The cubic–undecic conductor-twenty exclusion

**Theorem ConductorTwentyCubicUndecicExclusion.** Let \(F=(P,Q)\)
be a complex polynomial plane Keller map with \(P\) monic of positive
degree in its longitudinal coordinate \(x\). Its generic exact
directional order cannot have reduced characteristic sequence \((11,3)\).

Use the [exact specialized conductors](GENERIC_COORDINATE_ORDER_CONDUCTOR_AND_CLOSED_COLLISIONS.md)
and the [finite-gap argument](CONDUCTOR_TWELVE_HYPERELLIPTIC_AND_CUBIC_EXCLUSIONS.md#the-cubicseptic-branch).
At every finite parameter place, a generic conductor root reaches a
finite specialized root, its mate stays finite by monicity, and a
collapsed gap would give \(P_x=Q_x=0\). Thus some generic pair has a
unit original gap at every such place.

## The parameter equations

Let \(K=\mathbb C(y)\). Use monic reduced generators of degrees three and
eleven. Center the cubic at \(u=x-h_0(y)\) and subtract powers of the
cubic from the degree-eleven generator to remove degrees nine, six,
three and zero. Let \(\lambda\in K^\times\) be the determinant of these
target changes. Take \(u=\rho z,\ \rho^{14}=\lambda\) and multiply the
two generators by \(\rho^{-3},\rho^{-11}\), respectively, making the
normalized target determinant one. In this frame write
\[
 U=z^3+az,\qquad V=z^{11}+bz^{10}+cz^8+dz^7+ez^5+fz^4+gz^2+hz .
\]
The [parameter identity](CONDUCTOR_EIGHTEEN_CUBIC_DECIC_EXCLUSION.md#normal-form-and-first-integrals)
gives \(U_zV_y-U_yV_z=J\rho+\partial_z H(U,V)\).
The left side has degree at most twelve. The degrees up to thirteen
of \(\langle3,11\rangle\) are \(0,3,6,9,11,12\), so the evaluated
\(H\) lies in the span of \(1,U,U^2,U^3,V,U^4\). Degree twelve
gives \(3b'=0\). Put \(b=B\in\mathbb C\). The coefficients of
\(U^4,V,U^3,U^2,U\) in \(H\), obtained in that order, are
\[
 -\tfrac{11}{12}a',\quad \tfrac{3c'-10Ba'}{11},\quad
 \tfrac{ac'-8ca'}9,\quad \tfrac{22a^3a'+ae'-5ea'}6,\quad
 \tfrac{-a^4c'+8a^3ca'+3ag'-6a'g}{9}.
\]

Subtract these Hamiltonian terms. The remaining equations in
degrees nine, seven, six, four and three successively solve for
\(d,e,f,g,h\). The [exact algebra](verify_conductor_twenty_cubic_undecic_exclusion.py)
gives these polynomials by recurrence. At the step for coefficient \(v\),
the current residual coefficient is \(3v'+M\), where \(M\) is linear
in the remaining derivatives and contains no \(v'\). Replace \(v\) by
its constant of integration plus the polynomial primitive of \(-M/3\).
The primitive is obtained by integrating the coefficient of each
derivative in turn, in the order \(a,c\) followed by the remaining
coefficients, after subtracting the derivative of the primitive already
obtained. The recurrence checks the full derivative identity at each
step. The residual is then \(zR'+L'\), with \(R'=0,\ L'=J\rho\).
The integration constants \(B,D,E,F,G,H_0\) and the primitives \(R,L\)
have weights \(1,4,6,7,9,10,12,13\), respectively. The weights of
\(a,c\) are two and three.

For a matching pair let \(S\) be its normalized sum. The cubic divided
difference makes its product \(S^2+a\) and its squared gap
\(-3S^2-4a\). The negative divided difference of \(V\) evaluated at
this product is monic of degree ten in \(S\).

## The leading forms

At a pole of \(a\) or \(c\), let \(q>0\) be the larger of their pole
orders divided by two and three respectively. Let \(A,C\) be their
residues after scaling by a uniformizer to powers \(2q,3q\),
passing to an extension to take these powers if necessary. The
highest forms of \(R,L\) are
\[
 R_*=-\frac{2}{2910897}(29282A^6-98010A^3C^2+10935C^4),\qquad
 L_*=-\frac{4}{29403}A^2C(242A^3-135C^2).
\]
The constant \(R\) gives \(R_*=0\), hence \(A,C\) are both nonzero.
In particular, for an integer-valued normalized valuation on
\(K(\rho)\) the pole orders are \(2q,3q\) with \(q\) an integer.
Set \(A=1\) by a common scaling when testing nonvanishing. Then
\(C^2\) is a root of \(f(v)=29282-98010v+10935v^2\).
Since \(f(242/135)=-556358/5\), the leading term \(L_*\)
is nonzero, so \(L\) has pole order \(13q\).

The pair sum scaled by the same \(q\) has a monic integral equation.
At \(A=1\) the resultant of its initial polynomial with \(3S^2+4\),
reduced modulo \(f(C^2)\), is \(4576(45C^2-11)/9\). This is
nonzero since \(f(11/45)=29887/5\). The scaled sum is integral and
the scaled squared gap is a unit. The exact checker replays the
leading form identities and both nonvanishing tests.

## Finite places and the Kummer degree

If a coefficient pole occurs over finite \(y\), every normalized
gap has pole order \(q\). Since \(L\) has a pole, \(L'=J\rho\)
makes \(\rho\) polar there as well. Every original gap would
diverge, contrary to the finite-gap property. Thus \(a,c\) and
all normalized coefficients are regular above every finite \(y\).
The monic sum polynomial bounds the normalized gaps there, so
\(\rho\) cannot have a finite zero. Consequently \(\lambda=\rho^{14}\)
is rational with no finite zeros and has no pole at infinity.

If any of \(B,D,E,F,G,H_0,R\) is nonzero, the Kummer degree
\([K(\rho):K]\) divides both fourteen and its longitudinal weight.
Indeed the Galois action multiplies \(\rho\) by a fourteenth root of
unity, and the weight is the inverse scaling action on that constant.
Every listed greatest common divisor is at most seven, so the
ramification index \(e_\infty\) at any place above infinity is at
most seven. If either \(a\) or \(c\) is nonconstant, it has a pole
on the projective parameter curve, necessarily above infinity.
Then \(L\) has pole order \(13q\) with integer \(q\ge1\). The
valuation of \(L'\) is \(e_\infty-13q<0\), as \(y\) has pole order
\(e_\infty\) and the characteristic is zero. Thus \(\rho\) and
hence \(\lambda\) are polar above infinity, a contradiction.
If \(a,c\) are constant then \(L'=0\), also a contradiction.

## The fixed order

Suppose all seven constants vanish. The homogeneous equation
\(29282a^6-98010a^3c^2+10935c^4=0\) holds, with \(a,c\) both
nonzero since \(L'=J\rho\). Thus \(c^2/a^3=\nu\) is a nonzero
constant. Put \(\tau=c/a\) and \(A=C=1/\nu\), giving
\(a=A\tau^2,\ c=C\tau^3\). After substituting \(z=\tau v\) and
rescaling by their monic degrees, the generators are fixed
polynomials. Their matching gaps are fixed nonzero constants by the
resultant test. All original gaps are these constants times
\(w=\rho\tau\), so the finite-gap property forces \(w\) to be a
unit over every finite parameter. The unscaled cubic and undecic
coefficients \(Aw^2,Cw^3\) belong to \(K\), so \(w\in K\) and
\(w\in\mathbb C^\times\). Now \(L=L_*\tau^{13}\) with
\(L_*\ne0\), and \(L'=Jw/\tau\) gives \(\tau^{14}\) affine
nonconstant in \(y\). Hence \(\lambda^{-1}=\tau^{14}/w^{14}\)
has a finite zero.

The unscaled reduced generators in \(u=x-h_0(y)\) are a fixed
pair. Since a matching pair stays finite at every finite parameter,
\(h_0\) is polynomial. The original generators belong coefficientwise
to the fixed order \(R_0\) in \(u\) because membership in each
finite-degree part of the fixed order is defined by constant linear
equations. This gives polynomial-in-\(y\) ambient lifts to the fixed
plane curve. Each specialized order is contained in \(R_0\), of
colength ten, and has colength at most ten by exact specialization.
The orders therefore coincide.

At a fixed singular point, which exists by the positive colength,
apply the [singular tangent determinant](HYPERELLIPTIC_CONDUCTORS_AND_SINGULAR_TANGENT_DETERMINANT.md#the-determinant-at-a-finite-singular-point).
The difference between the generic target inverse and the polynomial
lifts is divisible by the fixed relation, so their derivatives
agree there. The determinant is \(\lambda^{-1}\). At every finite
parameter the specialized lifts generate \(R_0\), so their two
components after subtracting their values span the full
two-dimensional cotangent space. This forces the determinant to
be nonzero at the finite zero of \(\lambda^{-1}\), a contradiction.

Replay the algebra from the repository root.

~~~sh
.venv/bin/python research/hessian/hc4-triangular-graphs/verify_conductor_twenty_cubic_undecic_exclusion.py
~~~
