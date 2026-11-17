# Complete ruling-six finite-end coverage and the quadratic boundary residue algebra

PF2RULED162. 17 November 2026. Actual complex polynomial PLANE Keller maps.

The [ruling-five theorem](POLYNOMIAL_RESIDUE_RAMIFIED_LINE_AND_QUADRATIC_RULING_FIVE_CLOSURE.md)
excludes actual quadratic ruling-five coordinates. At ruling degree six,
every M=0 end is holomorphic. A central restriction of degree two leaves an
exact rank-two residual algebra and forces an original ramified prime with
a quadratic residue extension and explicit conductor. A central restriction
of degree three has no extra horizontal central prime.

The complete-curve argument confines every survivor to three off-M
weighted-face classes, with the necessary affine-response condition.
The ruling-six class remains open. All primes of the original full
normalization, including P-vertical primes invisible to the generic
P fibre, remain in scope. Neither generic derivative integrality nor
normal-tail termination is proved.

## 1. Actual image and notation

Use PF157 complete quadratic source-image coverage and PF161 notation:

    V: w^2=N+M*z, N=u*Sigma^2, M=u*A,
    u=t^2,w=t*(Sigma+A*x),z=2Sigma*x+A*x^2,
    ord_0 A=k>=1, m0=k+1>=2, Sigma(0)!=0.                   (1)

The original target B=C[P,Q] is in O(V), with full normalization S in
L=C(x,t), intermediate E=C(u,w), and full intermediate normalization S_E.
Normalize dP wedge dQ=Omega=du wedge dw/M. Its source pullback is
-2 dx wedge dt. Let an ACTUAL noninvertible target coordinate have ruling6:

    P=b3*z^3+w*(a2*z^2+a1*z+a0)+b2*z^2+b1*z+b0, b3!=0.     (2)

All coefficients are arbitrary polynomials in u, without degree bounds.
On every reached M-root line the restriction must have degree>=2 by PF52.
Consequently, if b3(c)=0 and M(c)=0,

    b2(c)+w_c*a2(c)!=0 for EACH reached w_c^2=N(c).           (3)

At the central zero this says b2(0)!=0 if b3(0)=0. At nonzero M roots
both signs occur. For the central restriction d=deg_x P(x,0),

    d=3 if b3(0)!=0; d=2 if b3(0)=0.                      (4)

The entire normalization S, Lambda and all primes remain unchanged by
any legitimate target changes below. Ruling6 is NOT generic Keller degree6.

## 2. Every finite M-root end is holomorphic

Work on every normalized branch of every geometric generic component P=p,
with coframe omega=-partial_p(w)*du/M. All affine points are regular.
An end over an M root with b3(c)!=0 is impossible by even z-leading
dominance (PF149). Suppose b3(c)=0; write v=u-c, e=ord_c M>=2,
b=ord_c b3>=1, and take v=s^alpha on a normalized branch.

If w tends to a finite value w0 with w0^2!=N(c), then b>=e; otherwise
b3*(w^2-N)^3 dominates all terms in the cleared equation. Dividing
M^3(P*-p) by M gives a polynomial in w over C[[v]], with p term -p*M^2.
Its reduction is

    (w^2-N(c))^2 * [c0*(w^2-N(c))+a2(c)*w+b2(c)],            (5)

where c0=(b3/M)(c). The bracket is not identically zero by (3), and
has degree<=2. Its root w0 has multiplicity ell<=2. The first p-dependent
coefficient argument of PF159, now with parameter term of order2e*alpha,
gives r>=ceil(2e*alpha/ell)>=e*alpha. Therefore ord omega>=alpha-1.

If w is polar, again b>=e, and reciprocating y=1/w after the same
division gives a polynomial with p term -p*M^2*y^6. At y=0 the first
unit coefficient has index ell<=2: either a2(c) is a unit at index1,
or b2(c) is a unit at index2, by (3). The constant term must vanish
for a small root. If beta=ord y>0, the first-coefficient argument gives

    ord(y_p)>=ceil((2e*alpha+6beta)/ell)>=e*alpha+3beta,
    ord omega>=beta+alpha-1>=0.                             (6)

No double reciprocal root or ramification index is deleted.

It remains that w->w0, w0^2=N(c), and z is polar. Put y=1/z. Its
exact equation from (2) is

    b3+y*(b2+w*a2)+y^2*(b1+w*a1)+y^3*(b0+w*a0-p)=0.         (7)

Its y derivative at (v,w-w0,y)=(0,0,0) is the UNIT in (3). Formal
implicit lifting over C[p][[v,w-w0]] gives y=b3*U, U a unit with
p-dependence divisible by v^(2b): differentiate (7) to obtain
 y_p=y^3/(unit). Thus every p-dependent coefficient of U is divisible
by v^(2b). Since z is polar, beta=b*alpha, and the relation becomes

    w^2=N+v^(e-b)*C(v,w,p), C a unit,
    C(v,w,p)-C(v,w,0) divisible by v^(2b).                   (8)

Here b<e is forced by w0^2=N(c). Its p-dependent part has order at
least(e+b)*alpha.

If N(c)!=0, the w derivative of (8) is a unit, so
ord(w_p)>=(e+b)*alpha and omega is holomorphic. If N(c)=0, then c=0
and N has a SIMPLE zero. When e-b>=2, w has exact order alpha/2;
2w dominates the other w derivative, giving

    ord omega >=(b+1/2)*alpha-1>=0.                         (9)

When e-b=1, including every leading contact, the coefficient of w^2
in the p-independent equation is a unit. Applying the same first-
coefficient argument at that index gives

    ord(w_p)>=ceil((e+b)*alpha/2),
    ord omega>=ceil(alpha/2)-1>=0.                         (10)

Equations(5)--(10) exhaust all finite M-root ends, including the central
simple N zero, both noncentral reached signs, every high contact, and
all original ramification. The calculation uses the entire equation,
not a finite prefix or a generic simple-root assumption.

## 3. Exact residual algebra after the two whole charts

Suppose first d=2 in (4), so b3(0)=0,b2(0)!=0. In ORIGINAL source
coordinates let f(t,X)=P(t,X)-p, degree6 in X. At t=0 it has degree2
with constant nonzero leading coefficient c2. There is a unique monic
quadratic factor W(t,X,p) over C[p][[t]], reducing to f(0,X)/c2.
Here is the all-order proof with the coefficient ring retained.

Write f=W*U, W monicdegree2 and deg U<=4, initially U0=c2. At each
new t order, divide the residual polynomial by the MONIC W0. The quotient
is the correction to U; the remainder divided by c2 is the degree<2
correction to W. No p-dependent denominator occurs. This recursively
proves the full factorization, not only a bounded truncation.

The actual deck action is (t,X)->(-t,-X-d0(t)),

    d0(t)=2Sigma(t^2)/A(t^2), ord_t d0=-2k.

It produces a second monic factor W2=W(-t,-X-d0,p). The resultant of
W and W2 has leading term d0^4, a nonzero constant times t^(-8k).
It is a UNIT in C[p]((t)). One may see this in the actual rank-two
multiplication algebra C[p][[t]][X]/W: multiplication by W2 has
leading matrix d0^2*identity and unit determinant. Therefore

    f=a6(t)*W*W2*R,
    R=X^2+c1(t,p)X+c0(t,p), c_i in C[p]((t)),               (11)

where a6(t)=b3(t^2)*A(t^2)^3 is a p-independent unit Laurent series.
This is an EXACT finite-free quadratic residual algebra, with ACTUAL
multiplication matrix

    M_X=[[0,-c0],[1,-c1]], M_X^2+c1*M_X+c0*I=0.             (12)

Its trace and determinant are -c1,c0. It is not original S or an HC4
normal-tail coefficient quotient. It retains BOTH extra generic sheets.
If W=X^2+W1(t,p)X+W0(t,p), comparing X^5 coefficients also gives

    c1=d0+t*a2(t^2)/b3(t^2)-W1(t,p)+W1(-t,p).              (13)

Thus the p-dependent part of the residual centre after shifting by d0/2
starts at order t>=1. All its polar centre terms are p-independent.
If d=3 instead, the same construction uses monic degree3 factors W,W2;
they exhaust all six sheets, and there is no residual central factor.

These counts concern P-HORIZONTAL points over t=0. P-vertical original
boundary primes are retained separately; absence from a generic fibre
never means deletion from S.

## 4. The central degree-two case forces a quadratic residue extension

In the d=2 case, the two extra sheets from R are outside BOTH whole
charts; their source X and translated source coordinate are polar.
By Section2 their coframe is holomorphic, so the ACTUAL companion is
finite on every such end.

For a first p-dependent coefficient lambda after a finite p-independent
Puiseux truncation, PF159 gives H(lambda)=p with H polynomial. Its degree
is at most TWO: each root of H(Y)-p lifts to a distinct extra Puiseux root
with the same prior polar chart data, and there are only two extra sheets.
The leading coefficient of every polynomial source function is polynomial
in lambda. In particular the finite companion residue is G(lambda) with
G in C[lambda]. Later coefficients lie in C(lambda) by ordinary Hensel.

An affine H gives polynomial residue in p and is excluded by PF161's
closed-graph/full-normalization argument. The same is true if the closed
prime has residue degree1 over C(p): with quadratic H,
C[lambda] intersect C(p)=C[H(lambda)], so G is again polynomial in p.
The graph is a full closed image in the finite normalization, not a
punctured parameter curve. The retained Chau COMPONENT Remark excludes it.

The fundamental local degree sum for R is two. Consequently the only
remaining possibility is ONE original horizontal prime R0, with

    ord_R0(t)=1, k(R0)=C(lambda), [C(lambda):C(P)]=2.          (14)

After an affine target P change and a constant lambda change take
P|R0=lambda^2. The actual Q residue is G(lambda). If G were in C[lambda^2],
its image would again be a forbidden graph. Hence C(P,G)=C(lambda),
and the residue degree of R0 over its target image is ONE.

The deck involution stabilizes the unique R0 and fixes P,Q. It therefore
fixes its entire residue field. Its inertia is TWO. In the intermediate
D0 below R0, u has order1 and w has an integral Laurent expansion in u.

w cannot tend to zero at D0. If it did, w^2=N+Mz and (7) force
ord_0 b3=k=m0-1 and b3*N/M=b2+O(u). At a p-independent positive-order
centre for w, the w^2 Taylor coefficient of P* has EXACT order
1-2m0, with nonzero leading term N*b2/M^2. Other lower and shifted
higher terms have larger order. A quadratic first response would require
its order to be -2r with r an integer, impossible for the odd number
1-2m0. Thus w tends to a nonzero constant or is polar.

Write ord_D0(w)=-a, a>=0. A simple initial root would lift over C[p]
and give a forbidden graph. Therefore the TWO residual roots have a
DOUBLE nonzero leading root in

    (b3/M)*w^2+a2*w+b2-(b3/M)*N.

Balancing its three terms against the unit b2(0) gives the NECESSARY data

    ord_0 a2=a, ord_0 b3=m0+2a,
    M_m*(a2_a)^2=4*b3_(m0+2a)*b2(0).                       (15)

Here M_m,a2_a,b3_(m0+2a) are their indicated nonzero leading constants.
A finite p-independent Laurent centre w_c(u) has

    w=w_c(u)+u^(m0+a)*lambda+...,
    P_initial=H(lambda), deg H=2.                           (16)

To check the exponent, the quadratic Taylor coefficient of P* has order
-2(m0+a): it is the nonzero second coefficient of the double quadratic
factor, times (w^2-N)^2/M^2. Lower z terms cannot change that leading
order. Making its lambda^2 term have order zero proves (16).

## 5. Exact original ramification, image and conductor

Put m=m0. The ACTUAL original prime in (14) has

    e_R=2a+2, f_(R/image)=1,
    -ord_R(x)=q=2(m+a)-1,
    ord_R(partial_p x)=2a+1=e_R-1.                          (17)

Indeed (16) gives the intermediate coframe order a and normal index a+1;
t^2=u doubles it. The source formula in (1) gives q. Every other original
prime, all singular finite points and all P-vertical primes remain retained.
There is no assumption that this is the only ramified prime of S. The
normal index e_R is NOT the finite fibre length at a cusp or node; no
equality at singular finite images is asserted.

The two visible central source/chart primes are unramified and each has
P-degree2. R0 has P-degree2 as well. These account for all six sheets
above t=0. The image normalization of R0 is C[lambda]: it is the integral
closure of C[P] in C(lambda), with P=lambda^2, so no finite lambda point
is missing. The actual ring S/R0 may be an intermediate nonnormal ring;
it is NOT silently identified with its normalization.

Subtract the even part of G by a legitimate target shear. In this image
frame, with T the new target coordinate,

    P=lambda^2, T=lambda*H(P), H in C[P] NONCONSTANT,
    Gamma: T^2=P*H(P)^2.                                  (18)

H=0 gives a graph, and nonzero constant H gives an embedded A1; both are
excluded by the retained Chau component theorem. Nonconstant H is the
explicit remaining image class. It can have a cusp over P=0 or multiple
image branches over nonzero H roots. Such components are NOT excluded by
that theorem. The target shear in (18) need not preserve a reduced companion
ruling degree; no simultaneous-frame claim is made.

Let A_Gamma=C[P,T] and B_tilde=C[lambda]. The ACTUAL residue multiplication
and image conductor are

    B_tilde=C[P] direct_sum lambda*C[P], lambda^2=P,
    M_lambda=[[0,P],[1,0]],
    (A_Gamma:B_tilde)=H(P)*B_tilde.                         (19)

For f(P)+lambda*g(P), membership and membership after multiplication by
lambda require BOTH f and g divisible by H. This proves the conductor.
The trace discriminants in bases(1,lambda) and(1,T) are respectively
4P and4P*H(P)^2. These are residue/image algebras over C[P], NOT S/B
or the original full codifferent. At P=0, B_tilde/(P) is dual numbers;
at P!=0 it splits geometrically. Globally lambda is not nilpotent since
lambda^(2j)=P^j. Even actual nilpotence in a finite conductor quotient
supplies no normal-index power containment for HC4 graph coefficients.

For the forced ramified prime, PF141/PF143 give the codifferent threshold
1-e_R and generic first-velocity pole q+e_R. Its contributions to FULL
first and marked-quadratic target clearing exponents are respectively

    2+ceil((2m-3)/(2a+2)) >=3,
    3+ceil((4m-5)/(2a+2)) >=4.                              (20)

Every other prime sharing an image curve remains in the full maximum.
This is pole detection, not a supply of marked regularity. Since the
original target degree is even and some source prime also lies above the
image (a nonconstant defining polynomial cannot be a source unit),
[L:F]>=e_R+2, so [E:F]>=a+2. No equality is inferred.

## 6. Complete global pole localization and the three surviving faces

PF149 excludes constant b3 on every complete quadratic model. Thus b3
is nonconstant. Its infinity bound, with deg M>=2 and deg b3>0, makes
EVERY u-infinity coframe end holomorphic. Section2 covers every M-root
end. At a finite zero c of b3 off M, the specialization P*(c,w) has
degree>=2 by PF52; degree>=3 gives holomorphic ends by PF158.

An actual nonzero rational exact coframe on the complete generic curve
must have a pole of order at least TWO. Therefore EVERY actual ruling-six
coordinate has a bad point c off M with

    b3(c)=a2(c)=b2(c)=a1(c)=0, b1(c)!=0.                    (21)

These are exactly the degree-two specialization conditions, not independent
assumptions on four unrelated divisors. On some polar branch at this c,
write u-c=s^alpha,w=c0*s^-h+..., c0!=0, and weighted initial I(Y).
PF159's multiplicity bound says r>=ceil(H0/ell)-h, H0>=2h.
If ell<=2 the end is holomorphic. Thus ell=3 or4. Its slope
L=h/alpha is an INTEGER: a denominator q>1 would give q conjugate
nonzero initial roots of multiplicity ell, whereas there are at most
four escaping roots. This yields exactly THREE face classes:

    A: I=C*Y^2*(Y-c0)^4, H0=2h, L>=3.
    B: I=C*Y^2*(Y-c0)^3*K(Y), H0=2h, L>=4,
       deg K<=1, K(0)*K(c0)!=0.
    C: I=C*Y^3*(Y-c0)^3,
       ord_c p3=B0, 4<=B0<L,
       (ord p3,ord p4,ord p5,ord p6)
           =(B0,B0+L,B0+2L,B0+3L).                         (22)

All displayed C,c0 are nonzero constants. In A every face coefficient
has its binomial value and ord p_i=(i-2)L for i=2,...,6. In B only
nonzero face coefficients obey this order; an absent highest term is
above the face. Here p_i are coefficients of w^i in P*. The lower
bounds follow from r+alpha<0 at a pole of order>=2. In C, p2 is the
unit specialization term above the face, giving B0<L. These classes
include all cancellations and multiple roots; no simple-root subset is
substituted for the classification.

The LATER first parameter response at every actual pole must be AFFINE
by PF159. A triple or quartic face alone is not that response and is
not sufficient for a mate. With P fixed, legitimate companion reductions
produce s>6,6 not dividing s,d=gcd(6,s) in{1,2,3}; the actual leading
Jacobian identity also requires (6/d) divides ord_c b3 at every c off M.
This companion-degree frame is separate from the optional image shear(18).

## 7. Boundaries and the full objective

For d=3, the two central whole charts exhaust all P-horizontal sheets;
its central omitted image has normalization A1 with P-degree3, not a
forced graph. Indeed the field index of C(P(0,x),Q(0,x)) in C(x) divides
THREE; index3 would put Q(0,x) polynomially in P(0,x) and give a forbidden
graph. Thus the index is1 and its finite normalization is C[x]. For d=2, (14)--(20) are the complete necessary quadratic
residue alternative. Neither alternative is declared realizable or excluded.
Together with (21),(22), they supply coverage before another coordinate
architecture expansion.

The original normalizations S,S_E, codifferent Lambda and EVERY ramified
prime remain the ones attached to the actual map. The rank-two formal
factor and conductor algebra do not replace them. The scalar coframe,
whole source charts, source-line criterion, actual deck inertia and retained
Chau plane theorem are identifiable two-dimensional inputs. Known nonproper
three-dimensional Keller controls continue to rule out dimension-free
integrality conclusions.

PF142's mixed minimal, PF156's fixed-v2 quintic graph, PF158's ruling-four
and PF161's ruling-five actual-target branches remain closed. Ruling-six
faces(22), their finite residue image class(18), higher cyclic images,
trivial deck, full marked first-derivative regularity and normal-tail power
growth remain OPEN. HC4/JC2 is not proved and no counterexample is constructed.
The uniform arguments, not finite coefficient prefixes, justify the scopes.
The checker passes 88 exact arithmetic and symbolic controls. No new
external classification or predecessor checker replay is used.
Independent whole-proof replay, formal verification and external review
remain unclaimed.
