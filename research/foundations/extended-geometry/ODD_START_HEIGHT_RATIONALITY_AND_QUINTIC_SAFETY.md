# Odd-start height rationality, the zero center and a quintic safe box

Work over any characteristic-zero field K on the sphere
S=K[x,y,t]/(t²+2xy-1), with normalized angular functional

    Lambda(f)=CT_x (1/2)integral_(-1)^1
                        f(x,(1-t²)/(2x),t) dt.

The integral means its exact rational polynomial formula. Consider
regular phases

    p=t^a(C(x)+t²D(x))/x,      a>=1 odd, C,D inK[x].

Regularity is C(0)+D(0)=0. Purity means vanishing of every positive
pure moment; safety means eventual vanishing for every fixed polynomial
marker. Least degree is the sphere polynomial filtration, equivalently
the radial-primitive balanced degree of its homogeneous two-pair lift.
These statements concern this two-height-term form. They do not classify
arbitrary height-five phases or arbitrary quintics.

## 1. Rational descent for every odd starting height

**Theorem OddStartRationalCenter.** If p is pure and C(0)!=0,
normalize by a nonzero scalar so C(0)=-D(0)=1. Put

    X=sqrt(-D/C), X(0)=1,
    A=C/X^a, L=1/X, u=x/A.

Then u,X andA all belong to K(x). Rationality is a consequence of
purity, not an additional ansatz.

For even m, substitution s=Xt gives the endpoint integrand
s^(am)(1-s²)^m. Its primitive has orderm+1 contact at both endpoints,
independently ofa. SinceX=1+O(x), replacing the actual endpoints
plus/minusX by plus/minus1 changes only terms of order at leastm+1.
The nonzero rational factor is

    kappa_(a,m)=integral_0^1 s^(am)(1-s²)^m ds
       =2^m m!/product_(j=0)^m(am+2j+1).

It is nonzero in characteristic zero; one may verify positivity over
the rationals without assuming that K is ordered. Thus purity is exactly

    [x^m]A^m L=0 for every positive evenm.                (A1)

Odd moments vanish by t parity. Formal residue substitution with
v(u) the inverse branch gives

    F(u)=L(v(u))*u*v'(u)/v(u)=1+uH(u²)
         =v'(u)*W(v(u)),       W=X^(a-1)/C inK(x).       (A2)

The last rationality uses that a-1 is even. Let V be the normalized
formal primitive of W. Integration gives V(v(u))=u+E(u²), with
E(z) in zK[[z]]. Therefore

    R=u²=x²(-D)^a/C^(a+2) inK(x),
    S=W/R' inK(x),
    S(x)=1/(2u(x))+E'(R(x)).                            (A3)

R has order two at zero. Its inverse branches v(u),v(-u) have
S-values differing by1/u. Apply the retained
[QuadraticInverseBranchDescent](ODD_CUBIC_HEIGHT_LOW_DEGREE_CLASSIFICATION.md),
proved there by CharacterSpanDescent. It forces sqrt(R) into K(x).
The specified formal square root has leading termx, hence u belongs
to K(x), without an extension of the coefficient field. Neither V nor
F was assumed rational.

Now X^a=uC/x is rational, and X²=-D/C is rational. Becausea is odd,

    X=(uC/x)/(X²)^((a-1)/2) inK(x),       A=x/u inK(x).

This proves the theorem.

Over the algebraic closure, parity of factor multiplicities consequently
gives polynomials

    C=HJ², D=-HK², H squarefree, H(0)=J(0)=K(0)=1.

Here the polynomial namedK is distinct from the coefficient field.
The chosen square root is X=K/J. Thus the exact necessary coefficients
are

    f_m=[x^m]H^m J^((a+2)m+1)K^(-am-1).                (A4)

The least-degree formula for a unit-center phase is

    max(degC+a-1,degD+a+1).                            (A5)

Indeed its unique t-degree-at-most-one normal form is

    t(1-2xy)^((a-1)/2)*((C+D)/x-2yD).

The two leading contributions have different highest y-degrees.
Reduction by t²=1-2xy never increases total degree. If leading
coefficients of C+D cancel, the D contribution dominates anyway.
This proves(A5) as the least degree, not just the displayed degree.

## 2. Complete zero-center classification in all degrees

**Theorem OddStartZeroCenterPurity.** Suppose C(0)=D(0)=0.
With no degree bound on C,D, the following are equivalent:

1. Every positive pure moment vanishes.
2. Lambda(p²)=Lambda(p⁴)=0.
3. C and D are both divisible byx².

Every such pure phase is safe for every fixed sphere marker and every
fixed original polynomial marker.

Write C=x(c+xC2),D=x(d+xD2). All terms of p except the core
c t^a+d t^(a+2) have strictly positive rotation weight, under
wt(x)=1,wt(y)=-1,wt(t)=0. All weights of p are nonnegative.
Therefore its pure moments equal those of the core, whose second
and fourth moments are

    c²/(2a+1)+2cd/(2a+3)+d²/(2a+5),
    c⁴/(4a+1)+4c³d/(4a+3)+6c²d²/(4a+5)
                    +4cd³/(4a+7)+d⁴/(4a+9).

Ifd=0 the second moment forcesc=0. Otherwise setz=c/d, and call
the two normalized polynomials P2(z),P4(z). Put

    D2=(2a+1)(2a+3)(2a+5),
    D4=(4a+1)(4a+3)(4a+5)(4a+7)(4a+9),
    A2=D2P2, A4=D4P4.

The quadratic is

    A2=(2a+3)(2a+5)z²+2(2a+1)(2a+5)z+(2a+1)(2a+3).

The coefficient ofz^(4-j) in A4 is binomial(4,j)D4/(4a+2j+1).
The six-by-six Sylvester matrix, with four shifted quadratic rows
and two shifted quartic rows, has the exact identity over Z[a]

    det Sylvester(A2,A4)=4096P(a),                     (A6)
    P(a)=13312a^8+132096a^7+544768a^6+1209600a^5
           +1565824a^4+1195296a^3+518064a²+115200a+10125.

All coefficients are positive. For every positive integera this is a
nonzero integer, hence is nonzero in K. A common rootz ofA2,A4
would make the Sylvester matrix annihilate(z^5,z^4,z³,z²,z,1)^T,
a nonzero vector. This contradicts(A6). Thusc=d=0.

This proves2 implies3, while1 implies2 is immediate. For3 implies1,
the phase isx t^a(C2+t²D2), of strictly positive weights unless zero.
Its powers have no weight-zero term, and the same holds after a fixed
marker for all sufficiently large powers. This proves sphere safety.

For original safety, every homogeneous lift contains Z in each
monomial and otherwise uses R,T,Z. Under xi1->s xi1,z1->s^(-1)z1,
withxi2,z2 fixed, R,T have weight zero andZ has weight one.
The two-pair contraction E2 is equivariant and takes values in
K[z1,z2], whose weights are nonpositive. Thus it kills every
STRICTLY POSITIVE input weight. The lift has positive minimum
weight and a fixed original marker has a finite weight range.
Their products therefore have only positive weights for sufficiently
large powers, proving original safety. This argument does not claim
that E2 kills arbitrary negative weights; for example E2(z1)=z1.
The zero phase is safe.

## 3. The degree-five safe box with starting height three

**Theorem OddStartThreeDegreeFiveSafety.** Every pure phase
t³(C+t²D)/x of least degree at most five is safe. At a unit center,
up to nonzero scalar normalization it is2t³y. At a zero center it is
x t³(alpha+beta*x), including zero.

Formula(A5) gives degC<=3,degD<=1. At a unit center the HJK
boxes are exactly degH=0 or1,degJ<=1,K=1. IfH=1, putJ=1+bx.
Equation(A4) at m=2 gives55b²=0, so b=0. Thus C=1,D=-1.

For linearH, rescale x and multiply the phase by a nonzero scalar
to normalizeH=1+x, and again putJ=1+bx. The necessary second
and fourth coefficients are

    A(b)=1+22b+55b²,
    B(b)=1+84b+1260b²+5320b³+5985b⁴.

They have the exact Bezout identity

    (20153693+671446930b+3272197005b²+3832387020b³)A(b)
      +(-12852499-35218260b)B(b)=7301194.                (A7)

Thus no linear-H branch is pure in characteristic zero. The sole
unit-center phase2t³y has strictly negative weight, so every fixed
sphere marker is eventually annihilated. For ORIGINAL safety, its
balanced lift has weight-1 and bidegree(b,b). Split a fixed marker
into its finitely many bidegrees(u,v). The output of E2 after
multiplication by the mth power has degreev-u in z, independent
ofm; it is zero ifu>v. Otherwise its possible weights are confined
to[-(v-u),0]. Equivariance preserves input weights, which are
w_marker-m and eventually lie below that fixed interval. Hence
the mixed contraction is eventually zero, for every marker component.
The zero-center theorem gives C=x²(alpha+beta*x),D=0 in this degree
box and supplies its all-marker safety. This proves the statement.

## 4. A prime-degree restriction at finite nonzero infinity

**Corollary OddStartFiniteNonzeroInfinity.** Suppose p is pure at
a unit center and its rational coordinate u has a finite NONZERO
value at infinity. Let r be the integer growth order of the rational
X there. Thenr>=1, and its exact least degree is

    d=(a+2)(r+1),                                    (A8)

which is composite. Equivalently, a pure unit-center odd-start phase
of prime least degree must have u vanishing or having a pole at
infinity. This is a necessary search restriction; neither infinity
condition is sufficient to construct a witness.

For proof, u=xX^a/C with a finite nonzero limit gives
c=degC=1+ar. The identity D=-CX² givese=degD=1+(a+2)r.
Since D is a nonzero polynomial,r<0 is impossible. Ifr=0 then
both C,D are linear. Rationality of X forces their roots to coincide;
normalization givesD=-C,X=1. WriteC=1+h x withh!=0. Then
p=2(1+h x)t^a y has only weights-1 and0, and its zero-weight
core ish t^a(1-t²). Its second moment ish² kappa_(a,2)!=0,
contradicting purity. Hencer>=1. Formula(A5) now gives(A8).
Both factors are integers at least two, and every nonzero unit-center
phase has degree at leasta+1>=2. Thus the stated prime-degree
restriction is exactly the contrapositive in this class.

## 5. Sources and exact verification

The only field-descent input is the retained
[CharacterSpanDescent and quadratic inverse-branch corollary](ODD_CUBIC_HEIGHT_LOW_DEGREE_CLASSIFICATION.md).
The endpoint argument above extends their input to every odd starting
height and proves rational X as well as rational u. The zero-center
classification is independent of rational descent and uses only the
displayed exact moment expressions and polynomial determinant.

The [stdlib checker](verify_odd_start_height_rationality_and_quintic_safety.py)
and [receipt](odd_start_height_rationality_and_quintic_safety_verification.json)
multiply the literal Bezout identity, construct the Sylvester matrix
over Z[a], sum all720 permutation terms, and verify its coefficient
identity(33 terms are nonzero). They also compare direct angular
moments with the odd-start endpoint coefficients in exact rational
arithmetic. These finite controls support the written all-degree
and all-order proofs; no bounded prefix is treated as a proof of
unrestricted purity or safety.

The work has internal mathematical review and exact computational
support. No external human review, formal verification, minimum-prime
claim or historical priority is asserted.
