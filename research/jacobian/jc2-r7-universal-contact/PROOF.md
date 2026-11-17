# Universal gap-21 contact obstruction and the three infinity escapes

14 November 2026. PF2R7PAR2 is an exact algebraic theorem for the entire
birational polynomial degree-(6,10), gap-21 family. PF2R7INF1 is its
geometric consequence for a degree-six Keller finite normalization with
exactly one ramified affine prime of generic (e,f)=(2,1). Omitted
unramified primes are retained. Neither statement assumes source finality.
The final corollary explicitly adds that external input. All geometric
arguments are written proofs; no formalization, external review or
novelty priority is claimed.

## 1. Normalization and retained gap equations

Use the affine parameter and affine target normalizations of PF2D6O1:

    p=t^6+a*t^4+b*t^3+c*t^2+d*t,
    q=[p^(5/3)+B*p^(4/3)+D*p^(2/3)+E*p^(1/3)]_+.

Here the bracket is polynomial part at infinity. The free target terms
C*p+F have been removed by an affine shear and translation. These
operations preserve quadratic graphs with the fixed target infinity
point [0:1:0] and the Hermite identity below. The gap-21 equations are
rho1=rho3=rho5=rho7=rho9=0 and rho11!=0. The literal coefficient
construction is retained in the archived gap checker, not recomputed
from missing generated artifacts.

The retained identity

    (rho5+(2*a/3)*rho3+(c/3)*rho1)|_(b=0)=-5*d^3/81

shows that b=0 forces d=0. In that case both p and q are even, so
C(p,q) is contained in C(t^2) and the parametrization is not birational.
Thus b!=0 throughout the stated family.

## 2. No fourth-order zero of p'

Suppose p' has a zero alpha of multiplicity at least four. Since p is
centered, its quintic derivative must be

    p'=6*(t-alpha)^4*(t+4*alpha).

This includes the multiplicity-five case alpha=0; that case would give
b=0 and is already impossible. For alpha!=0 scale t=alpha*u, p by
alpha^(-6), and q by alpha^(-10). Vanishing of the gap coefficients is
preserved. The normalized coefficients become

    (a,b,c,d)=(-15,40,-45,24).

Let r_j be rho_j at this specialization. Exact coefficient arithmetic gives

    r7 + (1882241/41553)*r1 - (8818601/41553)*r3
       + (2107/285)*r5 = 1661075456/6561.

The nonzero constant contradicts the first four gap equations. This is
a literal affine-linear polynomial identity in B,D,E, not a numerical
sample or a genericity assumption. Hence every finite zero of p' has
multiplicity at most three, including in special immersed strata.

## 3. No quadratic-contact involution

Put s=q'/p'. We prove there is no nontrivial Mobius involution iota of
P1 satisfying the rational identity

    (p(t)-p(v))*(s(t)+s(v))=2*(q(t)-q(v)),  v=iota(t).   (H)

If iota moves infinity to a finite alpha, then v-alpha has order t^(-1)
at infinity. Section 2 says s(v)=O(t^3). Meanwhile p(v),q(v) are bounded,
p(t)~t^6, q(t)~t^10 and s(t)~(5/3)*t^4. The leading coefficient at
t^10 of left minus right in (H) is 5/3-2=-1/3. This is impossible.
No squarefreeness of p' is required.

An involution fixing infinity has v=2*h-t for a complex constant h.
Clear denominators, differentiating before evaluation, and define

    T_h=(p(t)-p(v))*(q'(t)*p'(v)+q'(v)*p'(t))
        -2*(q(t)-q(v))*p'(t)*p'(v).

For the whole normalized polynomial-part family, exact identities are

    [t^17] T_h = 3840*h^3,
    [t^11] T_0 = (160/9)*b^3.

If T_h=0, the first forces h=0; the second then forces b=0,
contradicting Section 1. These two identities hold before specializing
the gap equations. The checker reconstructs the full symbolic contact
polynomial; it does not check a finite list of parameter values.

## 4. No degree-six constant-leading graph family

This part holds for any polynomials of degrees (6,10). Let g=gcd(p',q'),
k=deg(g)<=5, d0=p'/g and n0=q'/g. For a constant A consider, over C(x),

    y=A*(x-p(t))^2+q(t)+(n0(t)/d0(t))*(x-p(t)).         (L)

The reduced denominator is d0: modulo d0 its numerator is n0*(x-p),
coprime to d0 over C(x). The two relevant numerator degrees are

    A!=0: 17-k,     A=0: 15-k.

For A!=0 the highest term is from A*d0*p^2. For A=0 the coefficient is
nonzero because the leading terms in d0*q and n0*p have ratio 6:10.
Both degrees exceed six for every 0<=k<=5. In particular (L) never has
function-field degree six. This replaces the separate degree14/16
calculations for the two cusp examples by a uniform degree obstruction.

There is a second uniform contact obstruction useful below. For every
constant r put P=p-r and Q=(p-r)*q, of degrees (6,16). No nontrivial
Mobius involution satisfies (H) for this pair either. Its finite slope
poles still have order at most three, and the non-affine leading
obstruction is now 16/6-2=2/3. For the cleared contact polynomial the
two exact reflection coefficients are

    [t^23] T_h(P,Q)=-30720*h^3,
    [t^17] T_0(P,Q)=-(1280/9)*b^3.

They hold symbolically for every r, before the gap specialization, and
again force h=0,b=0. This polynomial target substitution is used only
to transport a contact identity; no Keller property is asserted for it.

Together Sections 1--4 prove PF2R7PAR2. They use no Keller compactification
or Borisov finality assertion.

## 5. The rank-eight pencil for every gap-21 image

Now assume the stated degree-six Keller map and entire single ramified
affine prime. Its branch curve has projective degree10. On the minimal
smooth resolution Y of its finite projective normalization W set
H=g^*L_infinity, K=K_Y and Q=2K+3H. JC2POS1--2 apply to this actual
model and whole source A2. Generic-line Riemann--Hurwitz gives

    sum f_i = 10-6+2-2*gamma = 6-2*gamma.

Since gamma>=2 and the sum is positive, gamma=2, sum f_i=2 and K.H=-4.
Hodge and Q^2>0 give 8<=rho<=11. This calculation does not require a
cusp. The proof in Section 2 of the earlier adjoint note now applies
verbatim: |K+H| has h0=2; its movable square is zero by integral Hodge;
a degree-one rational fibre violates Q nefness; a degree-two fibre
makes the fixed divisor exceptional; K on that divisor is nonnegative.
Thus rho=8 and the fixed divisor vanishes. We obtain

    F=K+H, F^2=0, K.F=-2, H.F=2,
    |F| a basepoint-free pencil of connected rational curves.

For an immersed image, PF2OMIT2 proves the affine normalization smooth
and counts exactly three unramified omitted primes. There are then
exactly two exceptional curves, both over the common infinity point A.
For a nonimmersed image, PF2R7CUSP1 is the complete two-class reduction.
Exclude its double-cusp alternative without the singular-at-infinity
assertion: the local normal double cover is z^2=f(x,y), with f defining
the ordinary cusp, and is the cyclic quotient A2 singularity. Its local
pullback divisor (f) is exactly 2*(z=0), supported only on E. Borisov's
Theorem 3.6 requires a coexceptional curve in g^*C-2E through that cyclic
quotient point. There is none in this local cover, a contradiction.
The other local factors lie at different normalization points and cannot
supply it. This excludes the double-cusp alternative directly. The
remaining cubic-cusp cover is smooth in the affine part and there are
three exceptional curves, all over A. These alternatives exhaust the
given polynomial family. The source input here is Theorem 3.6, whose
proof precedes and does not use Theorem 3.8; it has not been reproved here.

Import Borisov's type1 residue theorem and common-point boundary theorem,
with the target still P2. There are exactly two generic-infinity curves
C1,C2, of index pairs (1,5),(2,4),(3,3). Each C_i is smooth P1:
g|C_i is finite of degree one onto the normal target line. In particular

    d_i=F.C_i=-1-C_i^2 >=0.                            (I)

We have not asserted that d_i=0 or that C_i^2=-1.
Neither singularity of W at A nor any other conclusion of Theorem 3.8
is needed for Sections 5--8. In particular the old rank-count proof of
the double-cusp exclusion is not being used to hide that dependency.

## 6. Exclusion when both infinity curves are vertical

Suppose d1=d2=0. The actual divisor H consists of e1*C1+e2*C2 and
positive multiples m_j of exceptional curves B_j over A. Thus

    2=H.F=sum_j m_j*(F.B_j).

There is a horizontal exceptional curve mapping to the single target
infinity point a. A general F maps either two-to-one onto a line through
a or isomorphically onto a smooth conic. The line alternative is excluded
as in the earlier proof: the affine linear target pencil component would
factor as a cubic polynomial in a nonconstant source polynomial, and
would have a critical point, contradicting the Keller condition.

In the conic alternative, injectivity on F leaves a unique horizontal
exceptional section, of pole multiplicity two. Every general conic meets
target infinity only at a with multiplicity two. The branch curve also
passes through a, hence a=[0:1:0] in the retained coordinates. Consequently

    C(Y)=C(z,x),  y=a(z)*x^2+b(z)*x+c(z), a(z)!=0.

The affine normalization is smooth by Section 5. After shrinking the
pencil base, its resolved surface is a P1-bundle with the exceptional
infinity section removed at finite x. The generic coordinate x extends
over an open base as its affine fibre coordinate. Hence its only finite-x
ramification divisor is the specified simple prime E. No omitted e1
prime has been removed from this chart. The derivative equation
a'(z)*x^2+b'(z)*x+c'(z)=0 shows deg(E->P1_z) is one or two; E cannot be
vertical because H.E=10>H.F=2.

Degree one forces a'=0: otherwise the quadratic critical equation has
another root or a double root, violating respectively uniqueness or the
simple ramification coefficient. Tangency gives exactly (L), contradicted
by Section 4. Degree two gives a quadratic extension C(t)/C(z), whose
nontrivial automorphism is a Mobius involution. The two contacts with
the same parabola satisfy (H), contradicted by Section 3.

Thus the vertical-infinity case is impossible throughout the family.
This argument never uses the source-finality assertion.

## 7. Three remaining numerical escapes without finality

Intersect the actual divisor H with F and retain every exceptional prime:

    2=e1*d1+e2*d2+sum_j m_j*(F.B_j).

All terms are nonnegative integers. Section 6 excludes d1=d2=0. Up to
permuting the two infinity curves, the exhaustive remaining possibilities
are therefore:

| (e1,e2) | (d1,d2) | (C1^2,C2^2) | horizontal exceptional pole budget |
|---|---|---|---|
| (1,5) | (1,0) | (-2,-1) | 1: exactly one (m,F.B)=(1,1) |
| (1,5) | (2,0) | (-3,-1) | 0 |
| (2,4) | (1,0) | (-2,-1) | 0 |

In particular (3,3) is excluded independently of finality. This table is
an exhaustive list of necessary fibre-degree budgets, not a construction,
full boundary-graph classification or exclusion of the three escapes.

The actual determinant divisor sharpens the contact degrees in these
escapes. Since R_tot.F=(K+3H).F=4, its coefficient one on E, coefficients
e_i-1 on C_i, and zero coefficients on the unramified affine primes give

    4=F.E+sum_i(e_i-1)*d_i+sum_j r_j*(F.B_j).

Every contracted curve has r_j>=1: a morphism contracting a divisor has
zero differential determinant along it. Consequently, in the table's
order, the degrees F.E are respectively in {1,2,3}, exactly 4, and
exactly 3. In the first row its horizontal exceptional section has
r_B=4-F.E; in the other rows there is no horizontal exceptional curve.
These degree-three and degree-four alternatives explain why the
degree-two involution argument alone cannot close the remaining shapes.

If the source-finality assertion used inside Borisov's Theorem 3.8 is
added, the previously proved Q-intersection transport forces both d_i=0.
Section 6 then excludes the entire immersed family. Combined with
PF2R7EX1, this excludes every gap-21 image under the declared entire
single-ramified-prime hypothesis, relative to that additional source input.
The assertion was re-read in its primary-source proof; it has not been
independently proved or avoided for these last three escapes. Thus this
corollary must retain that dependency explicitly.

## 8. Explicit conic families for two escapes

In the first escape the general fibre meets infinity at the fixed point
a=[0:1:0] and at a moving point on C1. Since both pi|C1 and g|C1 have
degree one, choose the base coordinate z so that the moving point is
[1:z:0]. The image cannot be a line: its two distinct infinity points
would lie on that line. Writing the smooth conic equation and dividing
by its nonzero XY coefficient gives the exact rational graph form

    y=z*x+alpha(z)+beta(z)/(x-r(z)), beta!=0.            (G)

Here x is the generic fibre coordinate. Differentiating at fixed x gives
the numerator

    G_z=(x+alpha')*(x-r)^2+beta'*(x-r)+beta*r'.

Its reduced degree in x is three if r'!=0, two if r'=0,beta'!=0, and
one if r'=beta'=0. Indeed its value at x=r is beta*r'; in the latter
cases cancellation leaves respectively a denominator x-r or no
denominator. Each remaining critical polynomial must be irreducible
over C(z) with multiplicity one, since the complete affine ramification
divisor is the one simple prime E; its roots must avoid the pole x=r.

In fact the degree-one alternative is impossible. Here r,beta are
constants, z is a rational coordinate on E, p=-alpha'(z), and
q=z*p+alpha+beta/(p-r). At the unique infinity point of E, p has pole
order six. If its z-coordinate is finite, alpha has pole order five
and q has pole order at most six. If that coordinate is infinity,
alpha has pole order seven and q has pole order at most seven.
Both contradict q's pole order ten.

The degree-two alternative is also impossible. Here r is constant.
Apply the polynomial substitution X=x-r,Y=(x-r)*y. The family becomes

    Y=z*X^2+(alpha+r*z)*X+beta.

The map on E has polynomial parametrization P=p-r,Q=(p-r)*q. At its
two generic contacts over z, X!=0 and the substitution is locally
invertible, so the common quadratic graph gives the Hermite identity
for P,Q. The degree-two extension on E supplies its Mobius involution,
contradicted by Section 4's shifted-product obstruction. Hence in the
first infinity escape r'!=0 and F.E=3, rather than merely {1,2,3}.

In the third escape the sole generic infinity point is moving and has
intersection multiplicity two. Again write its direction [1:z:0].
First exclude the possible double-line image. If a general F mapped
two-to-one onto a line of direction z, that line has equation
y=z*x+c(z). Its unique infinity point is totally ramified under the
degree-two map from F, so there is an affine fibre coordinate w with
x=A(z)*w^2+B(z)*w+C(z), A!=0. The determinant is

    det(d(x,y)/d(z,w))=-x_w*(x+c').

The first factor has a simple finite root, and the second is quadratic.
If all finite roots agree, their product has multiplicity three there;
otherwise there is more than one ramification prime on a generic fibre.
Both contradict the unique simple E. Hence the images are smooth conics.

For such a conic the quadratic homogeneous part is (y-z*x)^2. Put
w=y-z*x. Its equation then solves as

    x=A(z)*w^2+B(z)*w+C(z), y=z*x+w, A!=0.              (V)

The determinant is now x_z-x*x_w, a cubic in w with leading coefficient
-2*A^2. It must be irreducible and simple over C(z); its normalization E
has projection degree three. These are concrete rational-function
families still to be excluded or realized, not Keller-source constructions.
The second escape has two moving infinity points on the same C1, with
pi|C1 of degree two and ramified-curve projection degree four; no involution
on E follows from that degree-four map.

## 9. Verification boundary and next gate

[verify.py](verify.py) checks the symbolic coefficient identities, literal
fourth-order-root certificate, gap parity identity, degree arithmetic,
the three complete integer budgets, a genuine two-contact positive control
and mutations of each principal certificate. [receipt.json](receipt.json)
binds local inputs. Surface geometry, local smoothness and the cited
boundary theorems are written inputs, not outputs of the checker.

The next independent geometric gate is to exclude or realize the three
horizontal-infinity cases without finality. The entire single-ramified-row
hypothesis must also be derived in any proposed application to the full
normal-r7 terminal passport. Normal r5, repaired even-r10 omitted-ledger
questions, nonnormal slices, (75,125) and JC2 remain open.
