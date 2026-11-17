# Cubic-height exclusion by a formal Airy--Riccati obstruction

14 November 2026. This extends the [affine-W classification](CUBIC_AFFINE_W_CLASSIFICATION.md)
using the retained [formal Gaussian--Lagrange identity](../../foundations/extended-geometry/FORMAL_GAUSSIAN_LAGRANGE_LEMMA.md).
It excludes a genuinely additional cubic monomial, T^3, in the same fixed
circular frame. The proof stays in formal series: it does not shift a Gaussian
integral by a singular saddle or use convergence of a cubic moment generating
function. The [Airy equation](https://dlmf.nist.gov/9.2.E1) supplies terminology;
all required differential and algebraic-function arguments are given below.
This is an internal proof without formalization, external review or a priority claim.

Put E(ZW)=E(T^2)=1, all other Wick pairings zero, and consider

    P=W*h(Z)+s(Z)*T^3+v(Z)*T^2+u(Z)*T+r(Z),          (1)

where all five coefficient functions are arbitrary complex polynomials.

## Complete classification and the cubic consequence

If every positive pure moment of (1) vanishes, then:

- On h(0)!=0, s is identically zero. The three-way classification in the
  preceding note applies unchanged: the constant-h Gaussian shear is safe,
  and the nonconstant linear-h family gives persistent GMC failures of
  degree at least four.
- On h(0)=0, h is divisible by Z^2 and s,v,u,r are divisible by Z. Every
  monomial then has strictly positive Z-minus-W weight, so the polynomial
  is GMC-safe with mixed vanishing for m>deg Q.

Both converse assertions hold. Thus the entire thirteen-coefficient cubic
family

    deg h<=2,  s constant,  deg v<=1,
    deg u<=2,  deg r<=3                            (2)

contains no GMC failure: its moment-null locus has s=0 and is exactly the
previous twelve-coefficient locus. In degree at most four, the extended
family has eighteen coefficients. Its failure locus is still exactly the
previous linear-h, linear-repair family with s=0. Only the safe zero branch
gains the possible term Z*T^3.

The proof below also gives a sufficient finite pure-moment gate, uniform
in the coefficients and their specializations:

    first 19D-22 pure moments zero => all pure moments zero
       within (1), for total degree D>=3.          (3)

This is 35 for cubics and 54 for quartics. These are sufficient bounds,
not optimality claims. When s is already zero, use the smaller 5D-3 gate.
The remaining cubic monomials outside (1) are
W^2,W^3,W^2Z,W^2T,WT,WZT,WT^2. No frame reduction is supplied for them,
so the global ternary cubic question remains open.

## Formal cubic Gaussian system

For series A(z),B(z),C(z) with A(0)=C(0)=0 and B(0)=-1/2, define

    I=E_T exp(A*T^3+(B+1/2)*T^2+C*T),
    J=partial_C I.

This is coefficientwise finite when A,B+1/2,C have positive z-order.
Before substitution treat A,B+1/2,C as independent formal source variables.
Gaussian integration by parts and source differentiation give

    3A*I_CC+2B*I_C+C*I=0,
    I_A=I_CCC,       I_B=I_CC.                     (4)

No integration contour is involved. Work in the Laurent field with A!=0.
For n=2,3,4 write I_(C^n)=p_n*J+q_n*I. Repeated differentiation of (4) gives

    p2=-2B/(3A),                   q2=-C/(3A),
    p3=(4B^2-3AC)/(9A^2),          q3=(2BC-3A)/(9A^2),
    p4=(-8B^3+12ABC-18A^2)/(27A^3),
    q4=(6AB-4B^2*C+3A*C^2)/(27A^3).                (5)

A prime denotes d/dz. Along any rational parameter curve A,B,C,

    I'=U*I+V*J,             J'=U2*I+V2*J,
    U=A'*q3+B'*q2,
    V=A'*p3+B'*p2+C',
    U2=A'*q4+B'*q3+C'*q2,
    V2=A'*p4+B'*p3+C'*p2.                         (6)

These identities hold in the formal Laurent field after substitution.

## No algebraic Airy Riccati solution

The equation f_x=x-f^2 has no solution algebraic over C(x).
Here is the elementary argument used in this proof.

At a finite x0, a finite branch value f(x0) has a unique analytic solution
of the regular first-order equation. A pole is handled by g=1/f:

    g_x=1-x*g^2,       g(x0)=0.

Again the solution is analytic, with g_x(x0)=1, so f has an ordinary
meromorphic simple pole and no branching. Equivalently, a lowest fractional
power in a local algebraic Puiseux expansion contradicts the regular ODE.
The normalized algebraic curve C(x,f) therefore has no ramification above
any finite x. If its degree over C(x) is d and it has r points above infinity,
Riemann--Hurwitz gives

    2g-2=-2d+(d-r)=-d-r.

Since g>=0 and r>=1, this forces d=1. Thus f would be rational. At infinity,
a pole of integer order k>=1 would make f^2 dominate both x and f_x;
a bounded f would instead leave the x term uncanceled. Neither is possible.
So there is no rational solution, completing the contradiction.

## Cubic Gaussian series cannot have rational logarithmic derivative

Assume A,B,C are rational and A has positive order q at zero, while B(0)=-1/2
and C(0)=0. In an algebraic Puiseux extension choose a=(3A)^(1/3) and set

    X=(B^2-3AC)/a^4,            f=-a*F-B/a^2.       (7)

The numerator of X is a unit at zero. Hence X has a pole of order 4q/3,
is nonconstant and has X'!=0. Direct algebra in (5)--(6) gives

    V=-a*X'.

For any solution of the system with I!=0, F=J/I satisfies

    F'=U2+V2*F-F*(U+V*F).

Substitution into (7) gives the exact identity

    f'=X'*(X-f^2).                                 (8)

The checker verifies both universal rational identities after substituting
A=a^3/3; this uses indeterminates, not numerical parameter samples.

If I'/I were rational in z, (6) and V!=0 would make
F=(I'/I-U)/V rational. Then f would be algebraic in z. Because X is a
nonconstant algebraic function, the one-variable function field C(z,a)
is a finite algebraic extension of C(X). Thus f would be algebraic over
C(X) and (8) would contradict the preceding Riccati lemma. We conclude
that I'/I is not rational.

For later use the same argument applies to *any* nonzero solution of (6),
not just the particular Gaussian series. In particular no solution of its
eliminated scalar equation can have rational logarithmic derivative.

## Applying the obstruction to the polynomial architecture

On h(0)!=0 let g=t*h(g) and use z=g, t=z/h(z). GL1 gives

    E(exp(tP))=exp(t*r(g))*I/(1-t*h'(g)),
    A=z*s/h,       B=-1/2+z*v/h,       C=z*u/h.

If s is not identically zero, A has positive finite order, B(0)=-1/2 and
C(0)=0. Vanishing of all positive pure moments would give

    I=H*exp(-R),
    H=(h-z*h')/h,             R=z*r/h.

This has rational logarithmic derivative H'/H-R', contradicting (8).
Therefore s=0, after which the preceding affine-W theorem applies.

On h(0)=0 the formal fixed point is g=0. The same formula has

    A=t*s(0), B=-1/2+t*v(0), C=t*u(0),
    H=1-t*h'(0), R=t*r(0).

If s(0)!=0 it gives the same contradiction. If s(0)=0 the scalar quadratic
formula and rational-exponential rigidity force h'(0)=v(0)=u(0)=r(0)=0,
as proved in the preceding note. Thus all five coefficient divisibilities
stated in the zero branch are necessary. Their converse and mixed safety
follow termwise from positive isotropic weight.

## A finite pure-moment bound

Eliminating J from (6) gives the scalar equation

    V*I''-(U*V+V'+V*V2)*I'
       +(-V*U'-V^2*U2+U*V'+U*V*V2)*I=0.           (9)

Since V!=0, any nonzero solution of (9) lifts back to (6).
The preceding argument therefore proves that substituting
I0=H*exp(-R) into (9) has a nonzero rational residual after division by I0.

Here are explicit degree bounds, retaining all vanishing-coefficient strata.
In total degree D>=3 put p=z*s, q=-h/2+z*v, c=z*u, j=h-z*h'. Then

    A=p/h, B=q/h, C=c/h, H=j/h,
    deg(p,q,c,h,j)<=(D-2,D-1,D,D-1,D-1).

Use the indicated common denominators without assuming they are reduced:

| Quantity | Common denominator | Numerator degree at most |
|---|---|---:|
| U | h^2*p^2 | 4D-5 |
| V | h^2*p^2 | 4D-6 |
| U2 | h^2*p^3 | 5D-6 |
| V2 | h^2*p^3 | 5D-7 |

Multiplying (9) by h^6*p^7 gives polynomial coefficients Q2,Q1,Q0 with

    deg Q2<=13D-20, deg Q1<=13D-19, deg Q0<=13D-18.

The target logarithmic derivative L=H'/H-R' has denominator h^2*j and
numerator degree at most 3D-2. Therefore

    Q2*(L'+L^2)+Q1*L+Q0

has, on denominator (h^2*j)^2, a numerator of degree at most 19D-24.
This numerator is nonzero by the Riccati obstruction. Its denominator is
a unit at zero, since h(0)=j(0)!=0.

If the first K pure moments vanish, GL1 and the invertible formal change
of parameter give I-I0=O(z^(K+1)). Applying the polynomial differential
operator Q2*d^2+Q1*d+Q0 makes the residual vanish to order at least K-1.
Taking K=19D-22 makes K-1=19D-23 exceed 19D-24, a contradiction whenever
s!=0. Thus that finite prefix already forces s=0, and the smaller 5D-3
quadratic-height gate completes the pure vanishing conclusion.

For h(0)=0 and s(0)!=0, the rational curve has A=t*s(0), B=-1/2+t*v(0),
C=t*u(0), H=1-t*h'(0), R=t*r(0). The same estimates apply using denominator
h=1 and numerator j=H, all of degree at most one; the degree table remains
valid for every D>=3. If s(0)=0, four moments suffice for the remaining
quadratic scalar expression. This completes both charts of (3).

The [checker](verify_cubic_height_obstruction.py) validates the universal
Riccati algebra, compares direct Wick moments to the GL1 formula including
nonzero cubic height, and checks polynomial differential residual degrees
on explicit generic cubic and quartic inputs. It does not replace the
algebraic-function or all-coefficient degree proofs.

## An additional closed WT face

The same obstruction closes a second thirteen-coefficient cubic family:

    P=W*(a+bZ+cT)+G(Z,T),       deg G<=3.           (10)

Every pure-moment-zero member is GMC-safe. If c!=0 it is actually quadratic:
the complete homogeneous cubic part G3 must vanish.

At fixed T, the circular Lagrange fixed point is now explicitly

    g(t,T)=alpha(t)+beta(t)*T,
    alpha=a*t/(1-bt),       beta=c*t/(1-bt),

with determinant 1-bt. Expand

    t*G(alpha+beta*T,T)=A(t)*T^3+V(t)*T^2+C(t)*T+R(t).

Pure vanishing would make the scalar cubic Gaussian series equal
(1-bt)*exp(-R). If A is nonzero, the Riccati obstruction applies and gives
a contradiction. But A=t*G3(beta(t),1). For c!=0, beta is a nonconstant
rational function, so A=0 forces the polynomial G3(w,1) to vanish identically.
Thus G3=0 and P is quadratic; the completed quadratic Gaussian theorem
gives mixed safety. If c=0, (10) belongs to the preceding univariate-h
cubic architecture, already closed including a=0. No additional coefficient
or frame-genericity assumption is used.

A coarse finite pure gate of 92 moments also suffices on (10). Put
q=1-bt. Each coefficient A,V,C,R has common denominator q^3 and a numerator
of degree at most four, since a monomial Z^i*T^j contributes
 t^(i+1)*q^(3-i) after clearing that denominator. Represent the target
H=q by h=q^3, j=q^4. These rational functions satisfy every degree-table
hypothesis above with the bound parameter D=6 (even though (10) is cubic).
Thus the same nonzero residual has numerator degree at most 90, whereas
92 zero moments force order at least 91. This forces A=0 when c!=0; the
remaining quadratic is closed by its first seven moments. When c=0, the
previous 35-moment gate is already stronger. Ninety-two is sufficient and
is not asserted optimal.

In the fourteen-coefficient family W*(a+bZ+dZ^2+cT)+G(Z,T), the faces c=0
and d=0 are now both closed. A failure in that declared family therefore
requires c*d!=0. This is a conditional local-family statement, not a frame
normal form for arbitrary ternary cubics. WZT, WT^2, higher W powers and
the simultaneous WT/WZ^2 interior remain distinct gates.

The subsequent [mixed WT/WZ^2 analysis](MIXED_WT_WZ2_CHART.md) excludes its
entire homogeneous-G3 spine with an optimal eight-moment gate. Its full
parent was then closed, including every exceptional parameter stratum, by
the separate [rescaled-period proof](MIXED_PERIOD_AT_INFINITY.md). The earlier
standard-basis run remains an UNKNOWN computation.
