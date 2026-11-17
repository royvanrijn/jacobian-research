# Polynomial mate responses and complete quadratic ruling-five escape coverage

PF2RULED159. 17 November 2026. Actual complex polynomial PLANE Keller maps.

This note proves the polynomial-mate response lemma and necessary
local escape conditions for the complete quadratic ruling-five source-image
class. The lemma is an input to the
[full ruling-five closure](POLYNOMIAL_RESIDUE_RAMIFIED_LINE_AND_QUADRATIC_RULING_FIVE_CLOSURE.md).
The local test alone does not exclude every rational mate; its affine-response
boundary is retained in the proof. All original ramified primes and the
full finite normalization remain in scope.

## 1. A polynomial mate forces an affine first parameter response at a pole

Let k be an algebraic closure of C(p), with the extension of partial_p.
Its constants are C, as retained in PF149. Let u=a+s^alpha, alpha>=1,
a in C, and let w in k((s)) satisfy P(u,w)=p, with P polynomial in w
and p-independent coefficients in C((s)). Assume the actual local form is

    Omega=k0(s)*s^(alpha-1)*ds wedge dw,
    k0 in C[[s]]*, k0(0)!=0.                                  (1)

For the quadratic image this is du wedge dw/M(u), with M(a)!=0;
then k0(0)=alpha/M(a). No finite-normalization stability is assumed.

Let r be the index of the FIRST coefficient of w outside C, and write

    w=w0+s^r*(lambda+O(s)), w0 in C[s,s^-1], lambda not in C.     (2)

The truncation is finite. Evaluating any polynomial T in w over C((s))
by its Taylor polynomial at w0 shows that its leading s coefficient is

    G_T(lambda), G_T in C[Y] nonzero.                          (3)

Indeed take the minimum order of the coefficients of T(w0+s^r*Y).
Its nonzero initial polynomial cannot vanish at lambda, because lambda
is transcendental over the algebraically closed C. Thus later Laurent
coefficients of w, even if rational in lambda, cannot alter this leading
coefficient. No bound on deg_w T is used.

In particular, for P(w)=p the minimum order is zero, and

    H(lambda)=p, H in C[Y], deg H>=1,
    partial_p lambda=1/H'(lambda).                            (4)

Suppose a mate Q polynomial in w over C((s)) satisfies dP wedge dQ=Omega.
The relative form along P=p is

    omega=-k0(s)*s^(alpha-1)*partial_p w*ds,
    ord omega=r+alpha-1.                                     (5)

If r+alpha=0, this has a nonzero simple pole and cannot be dQ. If
q=r+alpha<0, Q must have leading term G(lambda)*s^q, with G in C[Y]
nonzero by (3). Comparing dQ with (5) gives the POLYNOMIAL identity

    q*G(Y)*H'(Y)=-k0(0).                                    (6)

Consequently H is AFFINE and G is a nonzero CONSTANT. Therefore:

    at every pole of an actual polynomial-mate coframe,
    pole order>=2 and the first parameter response is affine.  (7)

A nonlinear H at a pole excludes ALL polynomial-in-w mates, of arbitrary
degree. This uses the actual mate, not only constant determinant of a
one-coordinate recurrence. It does NOT exclude rational-in-w mates.
For example P=s^6*(w-s^-8)^3 and
Q=s^-5*(w-s^-8)^(-2)/3 satisfy Jac_(s,w)(P,Q)=1. Their initial
response is cubic and the coframe has a double pole; Q has a denominator
in w, exactly outside (3). This is only a local Laurent-field control.

## 2. Multiplicity of an escaping initial root sharpens PF158

Let P(v,w) in C[[v]][w] have degree D, specialize to degree d>=1,
and take a polar branch v=s^alpha, w=c*s^-h+..., h>0. The leading
coefficient c is in C: the negative-order initial equation is independent
of p. Define

    -H0=min_i{alpha*ord_v p_i-i*h},
    I(Y)=the nonzero weighted initial polynomial of P,
    ell=multiplicity of the nonzero root c of I.               (8)

The coefficient p_d is a unit, so H0>=d*h. Let r be as in (2).
In P(s^alpha,w0+s^r*Y), the coefficient of Y^ell has EXACT order
ell*(r+h)-H0. This follows by differentiating I at its root c:
its ell-th derivative is nonzero, and all later terms of w0 have
higher order than c*s^-h. By (3),(4) the minimum order is zero, so

    r>=ceil(H0/ell)-h.                                       (9)

This retains every normalized branch and all earlier p-independent
polar terms. It is a refinement using actual root multiplicity, not a
simple-root approximation or a claimed uniform Laurent truncation.

Apply it to D=5,d=2. Every term of I has degree at least two, since
p_i are regular and i<2 cannot attain a weight at most -2h.
Thus ell<=3. If ell<=2, (9) gives r>=0, so the coframe (5) is
holomorphic. A polar end can therefore occur only when

    ell=3, H0=2h,
    I(Y)=C*Y^2*(1-Y/c)^3, C*c!=0.                            (10)

Here degree two MUST occur in I: otherwise a degree-five polynomial
with smallest degree at least three cannot have a triple nonzero root.
Since p_2 is a unit its weight is -2h, proving H0=2h.
All coefficients of degrees2,3,4,5 in (10) are nonzero. Hence

    L=h/alpha is a POSITIVE INTEGER,
    ord_v p_3=L, ord_v p_4=2L, ord_v p_5=3L.                  (11)

With gamma=p_2(0), and eta_i the indicated leading coefficients,

    eta_3^2=3*gamma*eta_4,
    eta_3^3=27*gamma^2*eta_5.                                (12)

The integer assertion follows already from ord_v p_3=h/alpha.
For a pole of order at least two, (9) also forces L>=4: when L<=3,
r+alpha>=0. At every actual polynomial-mate pole, (7) further forces
the LATER first parameter response H to be affine. A triple weighted
initial root alone does not determine that later response.

## 3. Complete audit of the degree-five quadratic-image ends

Use the COMPLETE [PF157 quadratic coverage](CANONICAL_DECK_JET_FIELDS_AND_RATIONAL_POLE_ATLAS_COVERAGE.md):

    V_R0: w^2=N(u)+M(u)*z,
    N=u*Sigma^2, M=u*A, A(0)=0, gcd(A,Sigma)=1,
    u=t^2, w=t*(Sigma+A*x), z=2*Sigma*x+A*x^2.                 (13)

Every M root has order at least two; N has a simple zero at zero,
and N is nonzero at every other M root. The original target coordinate
of ruling degree five has its UNIQUE expression

    P=w*(a2*z^2+a1*z+a0)+(b2*z^2+b1*z+b0), a2!=0.             (14)

All coefficients are arbitrary polynomials in u. There is no common
leading-coefficient admissibility hypothesis and no companion-degree
bound. We assume an ACTUAL noninvertible plane Keller target.

PF52's whole-source-line criterion forces, at EVERY M root c,

    b2(c)+w_c*a2(c)!=0 for EACH reached value w_c^2=N(c).       (15)

At c=0 this says b2(0)!=0. At other roots both signs occur on whole
fixed-t lines in (13), with z affine-linear. Their degree would fall
below two if (15) failed. In particular a2(c)=0 implies b2(c)!=0;
if a2(c)!=0, Delta=b2^2-N*a2^2 is a unit at c.

At any M-root infinity end with a2(c)!=0, a polar w is impossible
by highest-degree dominance. If w is finite and z polar, use PF150's
exact reciprocal coframe, with n=2 and unit a2. Since Delta(c)!=0,
its full noncontact calculation gives

    ord omega=(ord_c(M)+1)*alpha-1>=0.                       (16)

This applies also at the central N zero, because (15) makes Delta a
unit there. No contact or ramified branch has been suppressed.

If a2(c)=0, b2(c) is a unit. Finite w and polar z are impossible
by the unit leading z^2 term. For polar w put y=1/w and set

    F(u,y)=M^2*y^5*P*(u,1/y)
      =a2*(1-N*y^2)^2+b2*y*(1-N*y^2)^2
       +M*a1*y^2*(1-N*y^2)+M*b1*y^3*(1-N*y^2)
       +M^2*a0*y^4+M^2*b0*y^5.                              (17)

The equation is F=p*M^2*y^5, and its y derivative at (c,0) is the
UNIT b2(c). Thus, on every such branch,

    partial_p y=M^2*y^5/(F_y-5pM^2*y^4),
    omega=M*y^3*du/(F_y-5pM^2*y^4),                          (18)

which is holomorphic. This treats leading roots ON M without any
coprimality assumption. These are generic-curve computations, not an
identification of (17) with the original full finite normalization.

At u=infinity PF150's monic parameter argument applies without its
finite-root admissibility assumptions. With m=deg M>=2 and b=deg a2,

    ord omega>=ceil(((3m-5+b)/5)*alpha)-1>=0.                 (19)

All components, all normalized branches and all ramification indices
are included. At finite M*a2!=0 there are no ends with w,z polar.
At a finite root of a2 OFF M, P*(u,w) is polynomial in w with regular
coefficients. Its specialization degree is at least TWO by the actual
whole-source-line criterion; w is an affine-linear coordinate on that
fixed-t line. If that degree is at least THREE, PF158 makes every such
end holomorphic. Only degree-two specializations remain, and Section2
covers ALL of their polar branches.

The whole affine generic curve has regular nonzero omega. Since an
actual mate makes it rationally exact, some complete-curve point must
have a pole of order at least TWO. By (16)--(19) and Section2, that
point MUST lie off M at a triple initial root (10)--(12), with L>=4,
and its first parameter response must be AFFINE by Section1.

## 4. The coupled leading powers force orders(15j,10j,5j)

Keep P fixed. By PF158 all actual quadratic-image coordinates have
ruling degree at least five. Repeatedly subtract constant powers of P
from Q whenever 5 divides its ruling degree. This is a legitimate
polynomial TARGET automorphism, preserving B,F,S,Lambda and ALL primes.
It terminates with a companion of degree s>5, 5 not dividing s: a
resulting degree below five would already contradict PF158.

The actual scalar Jacobian leading identity retained in PF151 now gives

    p_5=alpha0*H^5, leading_s(Q*)=beta0*H^s,
    H in C(u)*, alpha0*beta0!=0.                             (20)

At the necessary bad point c off M, (11) gives ord_c p_5=3L.
Equation (20) makes 5 divide 3L. Thus L=5j for an integer j>=1.
For EVERY actual ruling-five target coordinate there consequently exists
c off M with the following canonical data, v=u-c:

    ord_c a2=15j, ord_c b2=10j, ord_c a1=5j,
    b1(c)!=0, j>=1.                                         (21)

Writing A1,A2,B2 for the leading coefficients of a1,a2,b2 at those
orders, the exact constants in (12) become

    M(c)*A1^2=3*b1(c)*B2,
    M(c)*A1^3=27*b1(c)^2*A2.                                (22)

Indeed p_5=a2/M^2, p_4=b2/M^2,
p_3=a1/M-2N*a2/M^2, p_2=b1/M-2N*b2/M^2. These identities prove
(21),(22) without discarding lower terms.

This is COMPLETE NECESSARY escape coverage, not an existence theorem.
Every ruling-five coordinate lacking a point (21),(22) is excluded,
with arbitrary lower coefficients and unbounded companion retained.
At every actual polar end the first parameter response must also pass
(7). The full deck order is TWO by PF157's ruling bound r>=2|D|.
No theorem here eliminates all affine-response survivors.

## 5. An unbounded excluded family passing the reached-line tests

On S2 take M=u^2,N=u, Sigma=-1,A=u, put v=u-1, and let h=5j,
j>=1. Consider the REGULAR surface polynomial

    P_h=w^2*(1-v^h*w)^3+v^(4h+1)*z^2+v^(h+1)*w+v.            (23)

Its ruling-five leading coefficient is -v^(3h), a FIFTH POWER up to
constant because h=5j. Its triple data at v=0 have L=h, so it passes
the coupled necessary leading conditions (20)--(22). Its surface a2
also vanishes on M, outside PF150's coprime-leading rational theorem.

In the whole source marking, the restriction to t=0 and t=+/-1 has
degree TWO, and every other fixed-t line has degree FIVE. The added
v^(h+1)*w makes the source t derivative nonzero on the ENTIRE t=0
line; the added v makes dP nonzero on the ENTIRE u=1 fibre, where
P=w^2 and P_u=1. These named differential and line-degree tests thus
do not themselves eliminate (23). No global submersion is asserted.

Put g=gcd(3,h), alpha=3/g, and make the EXACT rational substitution

    u=1+s^alpha,
    w=s^(-3h/g)+s^(-h/g)*Y,
    z=(w^2-u)/u^2.                                          (24)

The initial polynomial of P_h is -Y^3, and all remaining substituted
terms have STRICTLY positive s order. The last three repair terms
have order at least alpha. Hensel at -lambda^3=p gives a branch with
first parameter index r=-h/g. Its coframe pole order is

    ord omega=(3-h)/g-1<=-2.                                (25)

Section1 excludes EVERY mate polynomial in w over C(u), of any degree,
because the response is cubic. In particular no whole-plane Keller
target whose coordinate is (23) exists. ANY finite regular perturbation
whose WHOLE substitution in (24) has positive order leaves this local
obstruction intact. This last condition refers to the full substituted
polynomial, not a bounded coefficient prefix or a normal-tail claim.

## 6. A residue-free affine-response boundary control

The distinction in (7) is real. Still on S2, define

    v=u-1, S=1-v^10*w, g0(u)=3-2u=1-2v,
    P_aff=w^2*S^3+v^2*g0(u)*w*S+v^41*z^2+v^11*w+v.            (26)

It has fifth-power ruling-five leading coefficient -v^30, triple data
L=10, the same reached-line degrees2/5, and nonzero differential along
all leading-drop source lines. At u=1 it has TWO finite generic roots
and exactly THREE escaping normalized branches, all unramified over v.

For C=0 use w=v^-10+v^-2*Y. Its initial polynomial is -Y.
For C=+i,-i use

    w=v^-10+C*v^-4-C*v^-3+v^-2*Y.

Its initial polynomial is 2Y+C. Thus the three first responses are
AFFINE, and their coframes have double poles. Direct expansion gives

    P_w/v^2=-g0+O(v^2) for C=0,
    P_w/v^2=2*g0+O(v^2) for C=+/-i.                           (27)

These first two coefficients are independent of Y. Since
u^2*g0=(1+v)^2*(1-2v)=1+O(v^2), omega=-du/(u^2*P_w) has ZERO
residue at all three poles. Section3 makes every other end holomorphic.
Hence all local residues vanish, while the present polynomial-mate
leading obstruction is also absent. This does NOT prove exactness,
a global Keller pair, global submersion or a counterexample. Global
periods/polynomial-mate compatibility remain to be tested. It guards
against replacing the missing argument by pole existence or residue
cancellation alone.

## 7. Original full normalization and the remaining goal

Retain the original B=C[P,Q], R=C[x,t], L/F, FULL finite normalization S
and codifferent Lambda. No original ramified or unramified boundary
prime is deleted. Sections1--6 use polynomial source/actual image
functions, the real scalar Jacobian and COMPLETE geometric generic
curve branches. Parameter differentiation is not asserted to preserve
S or Lambda, and no field trace is promoted to marked target regularity.

The plane-specific steps are the one-dimensional complete generic
curve/scalar coframe, the whole-line injectivity theorem, and the actual
2x2 leading Jacobian identity. The retained three-dimensional nonproper
Keller control does not have this curve geometry and still forbids
any dimension-free integrality claim from constant Jacobian or derivative
generation. No new external classification is imported.

PF142's mixed minimal case, PF156's fixed-v2 quintic graph branch and
PF158's quadratic ruling-four branch stay closed. The ruling-five
residual is now localized by (21),(22) and the affine-response condition,
but remains OPEN. Trivial deck, higher quadratic/cyclic targets, full
marked quadratic regularity, generic first-derivative integrality and
normal-tail increasing ideal powers remain OPEN. Full HC4/JC2 is not
proved and no counterexample is constructed.

The universal proofs are (3)--(12) and the full-end audit (15)--(22).
The checker passes 82 exact arithmetic and symbolic controls. Predecessor
proofs/receipts remain pinned without checker replay. Independent whole-
proof replay, formal verification and external review remain unclaimed.
