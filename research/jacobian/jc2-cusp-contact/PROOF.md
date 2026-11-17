# F2 cusp contact: Wronskian bounds and the adapted-quintic Euler obstruction

14 November 2026. New scoped deductions from the retained F2 support
envelope, whole-source Keller condition, and k=1 implicit quintic.
The ambient (75,125) problem, JC2 and HC4 remain unresolved.

This note separates three objects: actual Keller restrictions to the
source cusp, commuting polynomial controls which are NOT Keller maps,
and an auxiliary filtration defined by the actual target quintic.
None is silently substituted for another. Claims JC2CFW1, JC2CFQ1 and
JC2CFC1 and JC2CFN1 are scoped below. No novelty priority, formal verification or
external review is claimed.

## 1. Replace a large power difference by a small differential

Keep the actual global source coordinates x,y of the F2 envelope and put

    x=z^5, y=z^2,
    f(z)=P(z^5,z^2), g(z)=-(5/9)*Q(z^5,z^2),
    D(z)=f(z)^5-g(z)^3,
    B(z)=5*f'(z)*g(z)-3*f(z)*g'(z).

The [retained cusp-contact theorem](../jc2-cusp-filtration-rigidity/PROOF.md)
gives exact source degrees delta_(2,5)(P)=195 and
delta_(2,5)(-(5/9)Q)=325. Its unique weighted supporting vertices prove
more: deg_z(f)=195 and deg_z(g)=325. Their leading terms are pure powers
of z; restriction to x^2-y^5=0 cannot remove those unique maxima.

The exact identity

    g*D'-3*g'*D=f^4*B                                     (1.1)

therefore gives, whenever B is nonzero,

    deg(D)>=4*195-325+deg(B)+1=456+deg(B).                (1.2)

The coefficient polynomial B has much lower computational complexity than
the full fifth/cubic power difference. An exact upper bound on deg(D)
can be tested through B, retaining all cancellations in (1.1).

For a polynomial source expression R, its cusp-normal form has x exponent
0 or 1. On restricting to x=z^5,y=z^2, the surviving w^0 monomials have
distinct z exponents. Thus

    delta_(2,5)(R)>=deg_z R(z^5,z^2)                      (1.3)

when that restriction is nonzero. This is an inequality, not identification
of the auxiliary filtration with a pre-existing dicritical valuation.

## 2. The global Keller hypothesis forces a precise critical ideal

Let J be the nonzero constant Jacobian of (P,-5Q/9).
On the source cusp the derivative column is

    (f',g')^T = DF(z^5,z^2)*(5*z^4,2*z)^T.

The polynomial matrix is unimodular over C[z], since its determinant is J.
Its polynomial inverse gives

    (f',g')=(5*z^4,2*z)=(z)                              (2.1)

as ideals of C[z]. This uses the WHOLE polynomial Keller map, not just
invertibility of its derivative at one point.

In particular B is divisible by z. Also B cannot be zero.
If it were zero, f^5/g^3 would be a nonzero constant; coprimality of 3,5
then gives f=alpha*h^3,g=beta*h^5 with deg h=65. Their derivative gcd has
degree 3*65-1=194, contrary to (2.1).

Equations (1.2)--(1.3) now improve the previous 456 endpoint:

    delta_(2,5)(P^5+(125/729)*Q^3)>=457.                 (2.2)

Thus the entire declared F2 contact-degree-456 stratum is excluded.
This is an additional necessary condition, not exclusion of all F2 maps.

More generally, if source cusp exponents are coprime a,b>=2 and the
restrictions of a Keller pair have degrees qe,pe, then their derivative
ideal is (z^(min(a,b)-1)). If min(p,q)*e>min(a,b), the same common-power
argument makes B=p*f'*g-q*f*g' nonzero and gives

    deg(f^p-g^q)>=(pq-p-q)*e+min(a,b).

Here the last term is min(a,b): ord B>=min(a,b)-1 and the differential
identity contributes one more. The F2 case has min(a,b)=2.
No general restriction-degree equality is asserted without proving its
leading source vertices or equivalent input.

### 2.1 The auxiliary source cusp has a birational degree325 image

There is also a global normalization constraint. The polynomial
parametrization (f,g) factors through the affine-line normalization of its
image as (f,g)=(p_0(H),q_0(H)); this is the same polynomial-normalization
fact used in the retained [target atlas](../../foundations/plane-jc/F2_AFFINE_TARGET_CURVE_ATLAS.md).
The covering degree d=deg H divides gcd(195,325)=65. Since H' divides
both derivatives, (2.1) gives d-1<=1, hence d<=2. Therefore d=1.

The auxiliary cusp image is birationally parametrized, has projective
degree325, and cannot coincide with any nonproper-value component in the
retained F2 atlas, whose degree is at most120. Consequently, for any such
component equation G, the polynomial G(P,Q) is NOT divisible by the source
cusp equation x^2-y^5.

This does not prohibit intersection with those components or factorization
along some other source curve. It gives an exact whole-curve exclusion for
one deliberately chosen auxiliary cusp, without claiming that it is an
omitted boundary divisor.

## 3. Complete lowest-jet classification, including target centering

Write the restricted source jets as

    f=f0+A*z^2+B2*z^4+C*z^5+O(z^6),
    g=g0+E*z^2+F2*z^4+G*z^5+O(z^6).

The source Jacobian at the origin is J=C*E-A*G !=0.
The first possible coefficients of B(z)=5*f'*g-3*f*g' are

    b1=10*A*g0-6*f0*E,
    b2=0,
    b3=20*B2*g0-12*f0*F2+4*A*E,
    b4=25*C*g0-15*f0*G.                                 (3.1)

If (f0,g0)!=(0,0), the literal identities

    50*g0*J=2*E*b4-5*G*b1,
    30*f0*J=2*A*b4-5*C*b1                               (3.2)

show that b1,b4 cannot both vanish. Hence the exact order of B is
one of 1,3,4, with the first nonzero coefficient selecting the branch.
The resulting necessary source-contact bounds are respectively

    457, 459, 460.                                      (3.3)

Now CENTER THE TARGET at the actual image F(0). These translated
polynomials are different from the uncentered ones, but retain their
degrees and constant Jacobian. Set f0=g0=0.
Then

    b3=4*A*E,
    b5=14*B2*E-2*A*F2,
    b6=19*C*E-5*A*G.                                    (3.4)

The vector (A,E) cannot be zero. Thus there are exactly these branches:

| source tangent/jet condition | exact ord B | necessary centered contact degree |
| --- | ---: | ---: |
| A*E!=0 | 3 | >=459 |
| A=0, B2!=0 (so C*E=J!=0) | 5 | >=461 |
| A=B2=0 (so b6=19J!=0) | 6 | >=462 |
| E=0, F2!=0 (so -A*G=J!=0) | 5 | >=461 |
| E=F2=0 (so b6=5J!=0) | 6 | >=462 |

The centered contact polynomial here is exactly

    (P-P(0))^5+(125/729)*(Q-Q(0))^3.

All orders in both tables occur for simple polynomial automorphisms,
verified as controls in the checker. Those controls do not have F2
degrees and do not claim attainment of the F2 global degree bounds.
Equations (3.2)--(3.4) are exact local certificates exhausting these
lowest-jet branches, independent of a huge coefficient elimination.

## 4. What the retained commuting freedoms can do

The archived top-layer calculation has a genuine Q-only kernel generated
by C_0^4,C_0^3,C_0^2,C_0,1 at resonant descents 5,10,15,20,25.
It does NOT prove those freedoms survive the full nonlinear circuit.

A precise polynomial control is

    P=s^3,
    R=s^5+c4*s^4+c3*s^3+c2*s^2+c1*s+c0.

The first five coefficients of P^5-R^3, at s degrees 14 through 10,
are triangular with successive leading terms -3*c4,-3*c3,-3*c2,-3*c1,-3*c0.
Their exact coefficient ideal is (c4,c3,c2,c1,c0), with no multiplicity
or hidden branch. If s is replaced by any genuine source carrier with
delta_(2,5)(s)=65, the successive first nonzero contact degrees are

    910,845,780,715,650.

In particular, any nonzero retained commuting coefficient on this declared
slice defeats a contact upper bound <=455. A polynomial target shear
can remove c3*P and c0, but it does not remove the other displayed
coefficients by that operation.

These pairs are algebraically dependent and have Jacobian zero.
They are exact controls for early commuting freedom, NOT F2 solutions,
full-circuit points or evidence of characteristic-zero Keller examples.
A full-circuit exclusion still requires the nonlinear equations.

## 4A. Two retained faces give the complete leading-degree spectrum

There is a stronger actual-Keller conclusion, independent of the small jet
bounds. Use ordinary positive weight omega=(5,2), the original carrier
weight eta=(5,-1), and the source cusp semidegree delta=delta_(2,5).
Let L=deg_eta(W) and D=deg_omega(W) for the ACTUAL polynomial

    W=P^5-g^3, g=-(5/9)*Q.

Write the retained carrier as C_R=y^(-5)*c(v), v=x*y^5.
Then c=v*(v-1)^2*R(v), deg c=5, with a simple zero at0.
The two classified carrier types have root multiplicities (1,2,1,1)
or (1,2,2), hence respectively four or three distinct roots.

The exact identity {g,W}=-5*J*P^4 is used on both Newton faces.
For eta, its degree shift is -(5-1)=-4. Thus

    L>=4*15-25+4=39.

The eta-leading band75 of P^5-g^3 cancels identically, so L<75.

If L>39, write W_L=y^(-L)*H(v), with H a polynomial.
The homogeneous leading bracket gives

    5*c*H'-L*c'*H=0.

Its valuation at the simple root v=0 forces 5|L; equivalently
H^5 is a scalar times c^L. Hence

    L=5m, 8<=m<=14, H=kappa*c^m, kappa!=0, deg H=L.      (4A.1)

If L=39, the leading nonzero bracket is instead

    5*c*H'-39*c'*H=-J*c^8.                              (4A.2)

At a root of c of multiplicity mu=1 or2, its valuation gives
ord H=7mu+1, since 5*ord H-39mu cannot vanish for either mu.
Therefore deg H>=35+r>=38, where r is the number of distinct roots.
At infinity a degree d!=39 would force d+4=40, hence d=36.
Thus necessarily deg H=39. In either branch we have proved

    D>=15*deg H-2L=13L>=507.                            (4A.3)

Now use omega. The same bracket identity has degree shift -7, with
deg_omega(g)=325 and deg_omega(P^4)=780.
Since D>=507>462, the top omega bracket vanishes.
The omega-leading part of g is a nonzero scalar times x^25*y^100.
For a monomial x^i*y^j in W_D, its bracket coefficient is 25j-100i.
Distinct monomials give distinct output monomials, so each nonzero term
requires j=4i. At a fixed omega degree this leaves exactly ONE monomial:

    W_D=w_D*x^n*y^(4n), D=13n, w_D!=0.

Its eta degree n is at most L, so D<=13L. Combining with (4A.3)
gives D=13L and n=L. The top omega face is a single monomial, so its
top term cannot cancel in cusp normal form: delta(W)=D.

Consequently every actual retained F2 Keller pair has exactly one of
the following necessary leading strata:

    L=39,40,45,50,55,60,65,70,
    delta(W)=507,520,585,650,715,780,845,910.             (4A.4)

These eight cases exhaust the highest-degree alternatives. Their
nonemptiness after all lower equations is NOT asserted.
In the resonant cases W_L=kappa*C_R^m with m=L/5, and kappa is the
nonzero coefficient defining that stratum. Successive vanishing of the
higher eta bands selects the lower cases. This is a classification from
the actual leading Keller equations, not from the Jacobian-zero controls
of Section4.

The lowest branch has an explicit forced coefficient. Write H=c^7*d(v).
Then (4A.2) is equivalent to

    5*d'-4*(c'/c)*d=-J.                                 (4A.5)

For the squarefree carrier, let

    R=(v^2-3v+3)/25,
    d=(25J/3)*v*(v-1)*R.

For a double carrier, let rho^2-3rho+1=0, alpha=3(rho+1)/5, and

    R=(v-rho)^2/[25*(1-rho)^2],
    d=J*v*(v-1)*(v-rho)*(v-alpha)/(rho*alpha).

Both satisfy (4A.5) exactly. Uniqueness follows because the homogeneous
solution would have order39/5 at the simple root0 and cannot be a
polynomial. The denominators are nonzero on the retained
carrier locus. Thus on L=39 the coefficient of x^39*y^156 in W is

    c_lead^7*d_lead,

with d_lead=J/3 in the squarefree case and J/(rho*alpha) in the double
case. In particular it is nonzero. The checker verifies these formulas
directly over Q and modulo rho^2-3rho+1.

Restriction to the source cusp preserves this single top monomial,
so deg_z D(z)=13L. Since 13L<975, the leading term in (1.1) does not
cancel. Hence the EXACT Wronskian-degree spectrum is

    deg_z B=13L-456
           =51,64,129,194,259,324,389,454.               (4A.6)

On the L=39 stratum its leading coefficient is -468*d_lead, since the
leading coefficients of f and g are c_lead^3 and c_lead^5.
This gives a much smaller differential object for the next coefficient
elimination.

Target translations preserve both retained leading faces, so the same
spectrum applies separately to the centered difference of Section3.
The small jet bounds remain correct but are strictly weaker here.

The original zero-defect upper target455 is therefore incompatible with
the actual leading Keller equations: they already force at least507.
This is a verified limitation of that proposed proof route, not a
contradiction in a hypothetical F2 map. A useful exclusion must now
eliminate one of the eight remaining leading strata, or introduce
additional actual-map information.

## 5. Test an adapted key polynomial before assuming an Euler derivation

The natural refinement is to replace the binomial u^5-v^3 by the actual
implicit target curve. First treat the complete lower-weight family

    K=u^5-v^3+t*u^3*v+B*u*v^2+A*u^4+U*u^2*v
        +C*v^2+V*u^3+Z*u*v+L2*u^2+L3*v+L1*u+L0.        (5.1)

The seven explicitly named coefficients have weights 14,13,12,11,10,9,8
for wt(u,v)=(3,5); the L terms have weight at most7.

Adjoin k=K and give normal monomials u^i*v^j*k^l, i<5, weights 3i+5j+7l.
Since (5.1) is monic in u, division gives a unique normal basis.
Its leading relation is X^5=Y^3, so this defines a positive source-ring
semidegree delta_K, regardless of the lower coefficients.

The degree-zero Poisson bracket is the same as for the (5,3) binomial
cusp; its graded fraction field has scalar Poisson centre. Therefore the
divergence-injectivity lemma of JC2CFR1 still applies. Existence of an
actual preserving polynomial derivation, however, must be checked.

The complete result is:

    delta_K has a nonzero preserving polynomial derivation
    if and only if all four following coefficients vanish:

    I1=t,
    I2=U+B^2/3,
    I3=V-2*A^2/5-B^3/27,
    I4=Z+2*B*C/3.                                      (5.2)

When they vanish, its derivations form one line, generated by

    D(u)=3*(u+A/5),
    D(v)=5*v-(2B/3)*u-5C/3+A*B/5.                      (5.3)

This is an all-coefficient theorem, not a sampled parameter search.

Proof: because delta_K(u)=3 and delta_K(v)=5, every preserving derivation
has the form

    L(u)=r0+r1*u, L(v)=r2+r3*u+r4*v.

The top condition gives 5*r1=3*r4. A nonzero preserving derivation cannot
have zero divergence by the scalar-centre lemma. Normalize divergence
to8, hence r1=3,r4=5. Successive weights14,13,12,10 force

    t=0, r3=-2B/3, r0=3A/5, r2=A*B/5-5C/3.

The entire remaining high-weight part of L(K)-15K is then

    -4*I2*u^2*v
    +(-6*I3-(2B/3)*I2)*u^3
    +(-7*I4+(6A/5)*I2)*u*v.                            (5.4)

Distinct normal monomials force exactly (5.2), and conversely those
conditions make (5.3) preserve all three generators and the whole ring.

Equivalently set alpha=A/5,beta=-B/3,gamma=-C/3. Conditions (5.2)
say K differs from

    (u+alpha)^5-(v+beta*u+gamma)^3

by terms of filtration degree at most7. Thus exactly the eligible
filtrations are affine filtered copies of the binomial filtration.
The affine change has determinant1. Their Euler derivation (5.3) is the
actual transported polynomial Euler field.

If any I_i is nonzero, NO polynomial coordinate change can turn delta_K
into any binomial-cusp filtration: a filtered polynomial isomorphism
would transport a nonzero preserving Euler derivation, which we have
proved does not exist. This conclusion covers nonlinear coordinate
changes too, without classifying those changes individually.

## 6. Apply the Euler test to the entire retained k=1 quintic family

Use the exact retained parametrization

    p(t)=t^3+a*t, q(t)=t^5+b*t^4+c*t^2+d*t,

and its twelve-support implicit quintic from
[PF2K1I1](../../foundations/plane-jc/F2_AFFINE_TARGET_K1_IMPLICIT_CONDUCTOR.md).
The [checker](verify.py) transcribes that polynomial, checks its
parametrization identity and independently recomputes its resultant.

Its seven high-weight coefficients in (5.1) are

    t=3b, B=-5a, A=ab+b^3+3c,
    U=-5a^2-4ab^2+3bc+3d,
    C=2a^2b-2ac,
    V=-abc+4ad+3b^2d+3c^2,
    Z=a^3b-3a^2c-5abd+3cd.

The exact Euler-eligibility ideal (I1,I2,I3,I4) is

    (b, d+10*a^2/9, c^2-25*a^3/81, a^2*c).             (6.1)

For a compact two-sided certificate, first reduce each I_i modulo b;
the discarded terms are visibly multiples of b. Set
g2=d+10a^2/9,g3=c^2-25a^3/81,g4=a^2c. The remaining three identities are

    I2|_(b=0)=3*g2,
    I3|_(b=0)=4*a*g2-(3/5)*g3,
    I4|_(b=0)=3*c*g2+g4/3.                             (6.2)

They have invertible constant triangular pivots, proving both inclusions.
Also

    a^5=(81/25)*(c*g4-a^2*g3).                         (6.3)

The lexicographic Groebner basis in b,d,c,a is therefore

    b, d+10a^2/9, c^2-25a^3/81, a^2*c, a^5,

with standard monomials 1,a,a^2,a^3,a^4,c,a*c.
The eligibility scheme has length7 and just one geometric point:
a=b=c=d=0. Nonreduced multiplicity is preserved; it is not counted as
seven different geometric cases.

Consequently every nonmonomial k=1 implicit-quintic filtration has NO
preserving polynomial derivation, and admits no polynomial filtered
conjugacy to a binomial-cusp filtration. The sole geometric exception is
the retained E8 endpoint u^5-v^3=0, whose monodromy and global attachment
obligations were already open. This does not exclude the associated
Keller map: it excludes this proposed Euler-preserving adaptation
mechanism on every other k=1 parameter row.

## 7. What remains for the ambitious F2 goal

The two-face spectrum (4A.4) excludes every degree below507 and all gaps
between its eight surviving degrees. The first branch has a forced nonzero
coefficient and a degree51 Wronskian. The small jet and commuting-slice
certificates remain useful controls. The full adapted-quintic Euler test
is also closed in its stated scope.

Still unproved are: a bound on the full actual F2 contact from the entire
954-variable retained circuit; whether any of the early commuting
freedoms survives all equations; a key-polynomial mechanism beyond the
Euler-compatible class; higher k target charts; and global F2 closure.
The missing historical circuit artifact was NOT regenerated.

A useful next calculation is to compile the tangent differential B
directly from the original source-band variables, compare its first
potential coefficients to the retained nonlinear equations, and obtain
a small exact coefficient-ideal certificate. The derivative-ideal and
jet certificates above are already available before that large step.
A filtration without a polynomial Euler derivation needs a genuinely
different transport argument; its central Poisson symbol alone does
not restore the existing one.

The all-degree bounds, ideal exhaustion and filtration statements have
written proofs. Finite checks replay identities, small coefficient ideals,
resultants and explicit controls. They are not formal verification or an
independent end-to-end review of the geometry.
