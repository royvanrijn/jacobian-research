# Complete positive-unit/opposite-unit family through degree seven (KCUB34)

Let K have characteristic zero, and let r>=4 be an integer. Every polynomial
Keller map in three variables, with source and output weights (1,-1,-r)
and ordinary degree at most seven, is an automorphism. After diagonal output
normalization it is exactly

    Z=c+eta*b^r+theta*a*b^(r+1),
    F(a,b,c)=(a, b+lambda*a^(r-1)*Z, Z),

with the following complete restrictions:

| r | coefficient restrictions |
|---|---|
| 4 | lambda*theta=0 |
| 5 | lambda*eta=lambda*theta=0 |
| 6,7 | theta=0, lambda*eta=0 |
| >=8 | lambda=eta=theta=0 |

The inverse is explicit: a=X, b=Y-lambda*X^(r-1)*Z, and
c=Z-eta*b^r-theta*X*b^(r+1), using the recovered b. These are compositions
of two elementary shears, including every coefficient-zero boundary.

There is no support or invariant-profile assumption. This closes the
proposed ratios r=4..7 and all larger r; it supplies no smaller counterexample
or unrestricted profile exclusion. The separate r=2 class contains the
known raw13 example by [KCUB15](WEIGHTED_COMPLETE_DEGREE7_CLASSIFICATION.md);
r=3 is excluded by [KCUB18](WEIGHTED_113_DEGREE7_EXCLUSION.md). Those theorems
and their external imports are not premises of this new r>=4 proof.

## Complete modules and degree budgets

Distinct nonzero weights force a diagonal linear part. A nonzero constant
Jacobian makes it invertible, so scale the outputs to normalize the linear
part to identity and the determinant to one. Set u=ab,w=a^r*c. On a!=0
write F=(a*f,a^-1*g,a^-r*h). Exponent comparison determines the a exponent
uniquely from each b,c monomial. Complete source enumeration, with no support
deletions, gives the following models.

For r=4,5:

    f=phi(u)+n*w,                  phi=1+mu*u+p*u^2+z*u^3,
    g=B(u)+A(u)*w,                 B=u+b2*u^2+b3*u^3+b4*u^4,
    h=E(u)+D(u)*w+C*w^2,           D=1+d1*u+d2*u^2+d3*u^3,
    deg A<=1, C constant,          E=e0*u^r+e1*u^(r+1).

Each model has 15 nonlinear source atoms: 4,5,6 by output.
For r=6,7 the same notation applies with n=C=e1=0 and A constant.
These models have 11 atoms: 3,4,4 by output. For r>=8 also A=E=0,
so f=phi, g=B, h=D*w. These are written degree bounds covering infinitely
many r, not conclusions from testing a finite list of ratios.

Write {P,Q}=P_u*Q_w-P_w*Q_u. The full source determinant is the polynomial
obtained by substituting u=ab,w=a^r*c in

    J=g*{f,h}+r*h*{g,f}+f*{g,h}=1.                         (1)

For example this follows by differentiating the quotient invariants
F1*F2=f*g and F1^r*F3=f^r*h. The identity initially holds on a!=0;
both sides are source polynomials, so it holds globally and no a=0
boundary is discarded. The independent checker compares it directly to
all four full source determinants, rather than trusting the formula.

## Nonconstant univariate first quotient is impossible

Suppose n=0 and f=phi is nonconstant. Work over an algebraic closure.
At a multiple root of phi both phi and phi' vanish, contradicting (1),
so phi is squarefree. With g=A*w+B and h=C*w^2+D*w+E, the three
coefficients of (1) in w are

    w^2: (2-r)*phi'*A*C+2*phi*A'*C,
    w^1: phi'*[(1-r)*A*D+2*B*C]+phi*(A'*D+2*B'*C-A*D'),
    w^0: phi'*(B*D-r*A*E)+phi*(B'*D-A*E').                (2)

If C!=0, at any root of phi the w^2 coefficient forces A=0 there,
since r>=4. The w coefficient then forces B=0 there. The constant
coefficient becomes zero, contradicting J=1. Thus C=0 globally.

If A=0, equation (1) becomes D*(phi*B)'=1. Polynomial units and D(0)=1
give D=1 and (phi*B)'=1. Since B(0)=0, integration gives phi*B=u.
But phi(0)=1, so a nonconstant phi cannot divide u.

If A!=0, the w coefficient is

    (1-r)*phi'*A*D+phi*(A'*D-A*D')=0.

Here D(0)=1. The derivative of A/(phi^(r-1)*D) is zero, so in
characteristic zero

    phi^(r-1)*D=k*A

for a nonzero field constant k. Its left degree is at least r-1>=3;
the right degree is at most one. This is impossible. Consequently every
univariate first quotient is constant, hence equal to one. This also
handles every r>=6; for r>=8 it immediately yields the identity map.
The rational-derivative argument has no generic-parameter denominator:
the case A=0 was handled separately, D is nonzero, and phi is assumed
nonconstant. Characteristic-zero coefficient descent causes no issue
because the final normal forms have coefficients in the original field.

## Graph first quotient when n!=0

Only r=4,5 permit n!=0. A diagonal conjugacy in c makes n=-1, so
f=N(u)-w with N=1+mu*u+p*u^2+z*u^3. Put v=f and define
G(u,v)=g(u,N(u)-v), H(u,v)=h(u,N(u)-v). Their complete shapes are

    G=G0(u)+v*G1(u),       deg G1<=1,
    H=H0(u)+v*H1(u)+T*v^2, deg H1<=3, T constant.

Initially deg G0<=4, deg H0<=6. Chain rule transforms (1) into

    (G+v*Gv)*Hu-Gu*(r*H+v*Hv)=1.                         (3)

At v=0, G0*H0'-r*G0'*H0=1. If m=deg G0,n=deg H0 and n!=r*m,
the leading coefficient has degree m+n-1, forcing m+n=1.
If n=r*m, the bound n<=6 and r=4,5 forces m<=1. Neither polynomial
can be zero. Thus in all cases

    G0=A*u+B, deg H0<=r, G1=alpha*u+beta.

The source point u=w=0 is u=0,v=1, where G=H=0. Thus B=-beta
and Hv=-1 there. The v^3 and v^2 coefficients of (3) are

    -(r+2)*alpha*T=0,
    2*G1*H1'-(r+1)*alpha*H1-(r+2)*A*T=0.                (4)

### Affine G is impossible

If alpha=0 then G=A*u+beta*(v-1) is affine, for arbitrary H and T.
If A=0, (3) is beta*(2v-1)*Hu=1, impossible in a polynomial ring.
If A!=0, use q=G as a coordinate. The vector field
(G+v*Gv)*partial_u-Gu*v*partial_v sends q to A*q and v to -A*v.
Equation (3) is therefore

    A*(q*Hq-v*Hv-r*H)=1.

Its polynomial solutions have constant term -1/(r*A) and all other
monomials q^(r+j)*v^j, j>=0. At the source q=0,v=1 they have value
-1/(r*A), contradicting H=0. This covers every affine G and any degree
of H, including all T!=0 cases forced affine by (4).

### The only nonaffine resonance is r=5

We may now assume alpha!=0, so (4) forces T=0. Since Hv=-1 at the
source, H1 is nonzero. Leading degrees in the remaining equation give

    2*deg H1=r+1.

For r=4 this is impossible. For r=5 the full polynomial identity gives
H1=gamma*G1^3. Its value at zero forces beta!=0 and gamma=-beta^-3.
A diagonal conjugacy in b normalizes beta=1; it leaves alpha unchanged
and covers every nonzero beta without taking roots. Now

    G0=A*u-1, G1=1+alpha*u, H1=-(1+alpha*u)^3.

At w=0 the original h has no u^j coefficients for j<5. Therefore,
using deg H0<=5,

    H0=[N*(1+alpha*u)^3]_(degree<5)+T5*u^5.

The normalized g has u coefficient one. This and the constant coefficient
of G0*H0'-5*A*H0=1 give respectively

    A+alpha+mu=1, mu+3*alpha+5*A=-1.

Thus A=-(1+alpha)/2 and mu=(3-alpha)/2. No parameter denominator is
inverted; the only denominators here are nonzero rational constants.

Put

    E0=G0*H0'-5*A*H0-1,
    E1=2*G1*H0'+G0*H1'-5*alpha*H0-6*A*H1.

The u^1,u^2,u^4 coefficients of E0 successively and uniquely give

    p=(2*alpha^2-alpha+3)/2,
    z=-(5*alpha^3-alpha^2+7*alpha-3)/4,
    T5=-alpha*(alpha+1)*(5*alpha^3-3*alpha^2+3*alpha-9)/40.

The corresponding pivot coefficients are -2,-3,-5. After these complete
eliminations, the u^1,u^2 coefficients of E1 become

    Q=alpha^2+4*alpha+6,
    S=(4*alpha^3+19*alpha^2+30*alpha+9)/2.

The five-monomial exact identity

    (-8*alpha^2/27-26*alpha/27-1/9)*Q
       +(4*alpha/27+10/27)*S = 1

excludes the remaining branch over every characteristic-zero extension.
This is a multiplied rational identity, not a primary Groebner verdict.
All graph branches are now excluded, so f=1 for every r>=4.

## Constant first quotient: exact two-shear normal forms

Now J={g,h}=1. Its w^2,w,constant coefficients are

    2*A'*C=0,
    A'*D+2*B'*C-A*D'=0,
    B'*D-A*E'=1.                                         (5)

We first exclude C!=0, possible only for r=4,5. Equation (5) gives
A'=0, and A=0 would force B'=0 and then contradict the last equation.
Thus A is a nonzero constant. Set t=C/A!=0. Integration and the normalized
values B(0)=E(0)=0,D(0)=1 give

    D=1+2*t*B, B+t*B^2-A*E=u.                            (6)

If d=deg B, leading degrees imply 2d<=deg E<=r+1. For r=4, d<=2,
so B=u+p*u^2. The missing u^2,u^3 terms in E give
p+t=0, 2tp=0. But 2t^2=2t(p+t)-2tp, contradiction.
For r=5, d<=3, so B=u+p*u^2+q*u^3. The missing u^2,u^3,u^4 terms give

    e2=p+t=0, e3=q+2tp=0, e4=t*(2q+p^2)=0.

The identity 5t^3=e4-2t*e3-t*(p-5t)*e2 is a contradiction.
Both identities are independently checked. Thus C=0.

With C=0, (5) implies either A=0, B=u,D=1 directly, or D=k*A
for a nonzero constant k. In the latter case
A*(k*B'-E')=1 makes A constant, and D(0)=1 gives D=1.
In all cases A=lambda is constant and integration gives B=u+lambda*E.
Substitution back into the complete source models gives precisely

    F3=c+eta*b^r+theta*a*b^(r+1),
    F2=b+lambda*a^(r-1)*F3.

The degree budgets give the restrictions in the opening table. Conversely
every listed parameter choice gives degree<=7, determinant one and the
explicit polynomial inverse there. This proves both necessity and
sufficiency, including the full zero-coefficient branches.

## Exact replay and limits

Run from this directory:

    python verify_exceptional_weight1.py

The [independent checker](verify_exceptional_weight1.py) uses stdlib Fraction
arithmetic and the existing [polynomial helper](verify.py). It enumerates
all 52 nonlinear source atoms and reconstructs all 72 native determinant
coefficient generators in the four complete finite models. It checks the
three universal univariate coefficients, four universal graph coefficients,
affine vector-field identities and 28 monomial eigenvalues, every pivot
and residual in the r=5 reduction, the five-term unit, both quadratic-tail
exclusions and all eight complete generic inverse charts.

The [receipt](exceptional-weight1/validation.json) pins the inputs, producers,
checker, arithmetic helper and this note. The [primary model producer](exceptional_weight1.py)
and [structural producer](exceptional_weight1_structure.py) use SymPy only
for discovery/export. Their primary nine-equation Groebner result is unused
as a proof premise. See the [artifact record](exceptional-weight1/README.md)
for the corrected initial serialization mismatch and retained discovery data.

The complete weight modules, all-r degree bound, root/field arguments,
Wronskian degree split, rational derivative consequences and diagonal
normalizations are written universal mathematics. They are not formalized
or externally reviewed; exact replay checks the finite algebraic premises.
No plane Keller classification or other external theorem is imported here.

The whole-profile counts remain 62 excluded / 12 remaining. The raw3D13,
normalized3D[10,48], identity-linear-part cubic12/36, free-linear-part
cubic12/31+17 and cubic-linear89 bounds do not change. Global support,
three-variable degree and special-form dimension minima remain open.
