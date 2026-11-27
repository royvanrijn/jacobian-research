# The cubic–decic conductor-eighteen exclusion

**Theorem ConductorEighteenCubicDecicExclusion.** Let \(F=(P,Q)\) be a
complex polynomial plane Keller map in a frame \((x,y)\) with \(P\) monic in
\(x\) of positive degree. Its exact directional coordinate conductor cannot
have global degree eighteen and generic reduced characteristic sequence
\((10,3)\).

Use the [exact specialization theorem](GENERIC_COORDINATE_ORDER_CONDUCTOR_AND_CLOSED_COLLISIONS.md)
and the [finite-gap argument](CONDUCTOR_TWELVE_HYPERELLIPTIC_AND_CUBIC_EXCLUSIONS.md#the-cubicseptic-branch).
At every finite parameter place, some generic matching pair has both
original coordinates finite and its original gap a unit. A conductor root
stays finite, its matching mate stays finite by monicity, and a collapsed
gap would make both longitudinal derivatives zero. These facts hold after
any finite extension used to split the pair equations.

## Normal form and first integrals

Let \(K=\mathbb C(y)\). The reduced sequence gives monic cubic and decic
generators over \(K\). Center the cubic at \(u=x-h(y)\) and subtract
powers of the cubic from the decic to remove degrees nine, six, three and
zero. Write \(\lambda\in K^\times\) for the target determinant of the
resulting reduction. Take \(u=\rho z\) with \(\rho^{13}=\lambda\) and
multiply the generators by \(\rho^{-3},\rho^{-10}\), respectively. The new
target determinant is one and the two polynomials are
\[
 U=z^3+az,\qquad
 V=z^{10}+bz^8+cz^7+dz^5+ez^4+fz^2+gz .
\]
As in the [cubic–octic reduction](CONDUCTOR_FOURTEEN_CUBIC_OCTIC_EXCLUSION.md#normal-form-and-coefficient-equations),
the parameter Jacobian satisfies
\[
 U_zV_y-U_yV_z=J\rho+\partial_z H(U,V),
\]
where \(J\in\mathbb C^\times\) is the original Jacobian and \(H\) is a
polynomial in the target variables. The left side has degree at most ten.
The semigroup \(\langle 3,10\rangle\) has degrees \(0,3,6,9,10\) up to
eleven, so \(H\) evaluated in the order is a linear combination of
\(1,U,U^2,U^3,V\). Comparing degree ten gives \(b=B+10a/3\) with \(B\)
constant. The coefficients of \(V,U^3,U^2,U\) in \(H\) then equal
\[
 \frac{3c'}{10},\qquad
 -\frac{2a'(12B+35a)}{27},\qquad
 \frac{ad'-5da'}6,\qquad
 \frac{24Ba^3a'+70a^4a'+9af'-18a'f}{27}.
\]
The residual equations in degrees seven, six, four and three give
\[
\begin{aligned}
d&=D+\tfrac45Bc+\tfrac73ac,\\
e&=E-\tfrac{28}9Ba^2-\tfrac{490}{81}a^3+\tfrac7{20}c^2,\\
f&=F-\tfrac{176}{81}Ba^3+\tfrac15Bc^2+\tfrac12Dc+\tfrac43Ea
       -\tfrac{910}{243}a^4+\tfrac7{15}ac^2,\\
g&=M-\tfrac89Ba^2c-\tfrac{10}9Da^2+\tfrac25Ec
        -\tfrac{112}{81}a^3c+\tfrac7{150}c^3 .
\end{aligned}
\]
Here \(D,E,F,M\in\mathbb C\) are constants. The residual is \(zR'+L'\), with
\[
\begin{aligned}
8100R={}&1200Ba^3c-324Bc^3+1500Da^3-1215Dc^2\\
       &-3240Eac-4860Fc-8100Ma+700a^4c-378ac^3,\\
-4374000L={}&144000Ba^5-291600Ba^2c^2-729000Da^2c-648000Ea^3\\
       &+262440Ec^2-1458000Fa^2+1312200Mc
        +70000a^6-226800a^3c^2+15309c^4.
\end{aligned}
\]
Thus \(R\in\mathbb C\) and \(L'=J\rho\). The longitudinal weights of
\(B,D,E,F,M,R,L\) are \(2,5,6,8,9,11,12\), respectively. The weights of
\(a,c\) are two and three.

## Finite places and infinity

For a matching pair let \(S,p\) be its sum and product in \(z\). The cubic
divided difference gives \(p=S^2+a\), and the squared gap is
\(-3S^2-4a\). Substitute this product into the negative decic divided
difference to obtain a monic polynomial of degree nine in \(S\).

At a pole of \(a\) or \(c\) let \(q>0\) be the maximum of their two
pole orders divided by two and three, respectively. After extending the
local field to make \(q\) integral, let \(A,C\) be the residues after
multiplying \(a,c\) by powers \(2q,3q\) of a uniformizer. At least one
residue is nonzero. The leading form of the constant \(R\) is
\(7AC(50A^3-27C^2)/4050\). Thus \(C=0\), \(A=0\), or
\(C^2=50A^3/27\) with \(AC\ne0\). In the same order the leading form of
\(L\) is \(-35A^6/2187\), \(-7C^4/2000\), or \(595A^6/8748\),
so \(L\) has pole order \(12q\).

After multiplying \(S\) by the uniformizer to power \(q\), the scaled
sum satisfies a monic integral equation with initial polynomial
\[
\begin{aligned}
P_{A,C}(S)={}& S^9+\tfrac83 AS^7-CS^6+\tfrac73 A^2S^5
  -\tfrac53 ACS^4+(\tfrac{50}{81}A^3+\tfrac7{20}C^2)S^3\\
 &-\tfrac23 A^2CS^2+(-\tfrac5{243}A^4+\tfrac7{30}AC^2)S
  +\tfrac4{81}A^3C-\tfrac7{150}C^3 .
\end{aligned}
\]
The resultant of this polynomial with \(3S^2+4A\) is
\(4900A^9/9\), \(107163C^6/2500\), or \(41650A^9/9\) in the three
regimes. Every scaled sum is integral and every scaled squared gap is a
unit. Thus every normalized gap has pole order \(q\). At a finite place
\(L'=J\rho\) gives a pole of \(\rho\) since differentiation in \(y\)
increases the pole order of \(L\) by the ramification index. Then all
original gaps diverge. The finite-gap property rules this out. It
follows that \(a,c,L\) and all normalized coefficients are regular at
finite places. The monic sum polynomial bounds the normalized gaps at
every finite place, so a finite zero of \(\rho\) is also impossible.

If \(\rho\) has a finite pole, the place is ramified over \(K\) since
\(L\) is regular. The Kummer extension \(K(\rho)/K\) has degree dividing
thirteen, so inertia at that place acts by all thirteenth roots of unity
on \(\rho\). This action sends each normalized coefficient to its
multiple by the root of unity raised to its inverse longitudinal weight.
The constants \(B,D,E,F,M,R\) must be fixed, and their weights are nonzero
modulo thirteen. Therefore a finite pole forces all six constants to
vanish. If any constant is nonzero, \(\rho^{13}\) has no finite poles or
zeros and is a nonzero constant. Then \(a,c,L\) are polynomials in \(y\).
If either \(a\) or \(c\) is nonconstant, the same leading forms at infinity
give \(\deg L=12q>1\). If both are constant, \(L'=0\). This contradicts
\(L'=J\rho\in\mathbb C^\times\) in either case.

## Homogeneous case

Suppose \(B=D=E=F=M=R=0\). Then \(ac(50a^3-27c^2)=0\)
identically. At least one of \(a,c\) is nonzero by \(L'=J\rho\).
Write \(a=A\tau^2,\ c=C\tau^3\), choosing \(\tau^2=a\) if \(c=0\),
\(\tau^3=c\) if \(a=0\), and \(\tau=c/a\) in the third case. The fixed
pair \((A,C)\) is \((1,0)\), \((0,1)\) or \((27/50,27/50)\),
respectively. Substituting \(z=\tau v\) gives fixed monic polynomials
after dividing by \(\tau^3,\tau^{10}\). Their matching pairs form a
finite set by the degree-nine sum equation and the product formula.
Their gaps are nonzero by the resultant tests. Consequently all original
gaps are fixed nonzero multiples of \(w=\rho\tau\). By the finite-gap
property \(w\) has neither a zero nor a pole over any finite parameter.

If \(A\ne0\), the linear coefficient of the unscaled reduced cubic is
\(Aw^2\in K\). If \(C\ne0\), the degree-seven coefficient of the
unscaled reduced decic is \(Cw^3\in K\). Either nonzero coefficient
shows that the corresponding power of \(w\) is a rational function
without finite poles or zeros. Hence \(w\in\mathbb C^\times\).
The leading form for \(L\) gives \(L=L_0\tau^{12}\) with
\(L_0\in\mathbb C^\times\), and \(L'=Jw/\tau\) makes
\(\tau^{13}\) nonconstant affine in \(y\). In particular
\(\lambda^{-1}=\tau^{13}/w^{13}\) has a finite zero.

The unscaled reduced generators \(U_0,V_0\) are fixed monic polynomials
in \(u=x-h(y)\). Since at least one pair is finite at every finite
parameter and all original coordinates of matching pairs are \(h(y)\)
plus fixed offsets, \(h\in\mathbb C[y]\). Thus \(P,Q\), expressed in
\((u,y)\), have polynomial parameter coefficients. Membership in each
finite-degree part of the fixed order \(R_0=\mathbb C[U_0,V_0]\) is
defined by constant-coefficient linear equations. Since \(P,Q\) belong
to the generic order, all their parameter coefficients lie in \(R_0\).
They admit ambient target lifts \(\Psi_y(T,W)\) polynomial in \(y\).
The specialized order is contained in \(R_0\), of colength nine, and has
colength at most nine by the exact specialized conductor. It therefore
equals \(R_0\) at every finite parameter.

Take an affine singular point of \(R_0\), which exists because \(R_0\)
has positive colength in its normalization. The [singular tangent
argument](HYPERELLIPTIC_CONDUCTORS_AND_SINGULAR_TANGENT_DETERMINANT.md#the-determinant-at-a-finite-singular-point)
identifies the ambient Jacobian of \(\Psi_y\) at the singular point
with \(\lambda^{-1}\). Indeed the difference from the generic ambient
target inverse is divisible by the fixed plane relation in each
component, so the first derivatives agree at the point. At every finite
parameter the specialized lift generates \(R_0\), so its two components
after subtraction of their values span the two-dimensional cotangent
space at the singular point. Its determinant is therefore nonzero there.
This contradicts the finite zero of \(\lambda^{-1}\).

## Exact algebra

The [checker](verify_conductor_eighteen_cubic_decic_exclusion.py) replays
the parameter integrations, weights, collision equation and the three
leading-form tests. The valuation and singular tangent steps are given
above. Replay from the repository root.

~~~sh
.venv/bin/python research/hessian/hc4-triangular-graphs/verify_conductor_eighteen_cubic_decic_exclusion.py
~~~
