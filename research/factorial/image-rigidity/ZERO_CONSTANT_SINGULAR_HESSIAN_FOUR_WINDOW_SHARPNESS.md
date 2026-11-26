# Sharpness of the singular-Hessian four-moment window

The uniform four-window for zero-constant singular-Hessian binary
quadratics cannot be shortened to three. This note gives an exact
five-monomial example with its first three factorial moments zero.

The only imported theorem is the uniform four-window of
ZERO_CONSTANT_SINGULAR_HESSIAN_BINARY_QUADRATIC_STRONG_FACTORIAL.md,
SHA256 2de47281b84b6a05d64b9e31849e72e3b2aa9a47adf1b23ffc891d4b56f02e4a.
It proves the fourth moment of the example is nonzero; the first three
zeros are verified independently below.

Define

 P(T)=9T^12+12T^11+19T^10+16T^9+T^8-8T^7-17T^6
           -8T^5+T^4+16T^3+19T^2+12T+9,
 J(T)=(T^2-1)(3T^2+T+3),
 D(T)=-2(6T^6+4T^5-4T^2+3),
 E(T)=2(3T^6-4T^4+4T+6).

Choose any root t of P in an algebraic closure of Q, and put

 d=D(t)/J(t), e=E(t)/J(t),
 f(X,Y)=(tX+Y)^2+dX+eY.                                 (1)

**Theorem 1 (optimal uniform window).** Formula(1) is well defined and
has exactly five nonzero monomials. Its quadratic Hessian is singular,
and

 Lambda_2(f)=Lambda_2(f^2)=Lambda_2(f^3)=0,
 Lambda_2(f^4)!=0.

Therefore the uniform four-window on the entire zero-constant
singular-Hessian family is optimal. This is sharpness of that stronger
uniform window, not a claim that the actual five-term SFC bound is sharp.

## 1. The denominator is nonzero

P is a nonconstant polynomial, so a root t exists in Qbar. Its constant
coefficient9 shows t!=0. Direct evaluation gives P(1)=81 and P(-1)=1.
To exclude the other factor of J, set

 R(U)=9U^6+12U^5-35U^4-44U^3+6U^2+4U+1.

The exact reciprocal identity is

 P(T)=T^6 R(T+T^(-1)),  T!=0.

If 3t^2+t+3=0, then t+t^(-1)=-1/3, but R(-1/3)=121/81!=0.
Thus P(t)!=0 in that case, a contradiction. Hence J(t)!=0.

No irreducibility or numerical root choice is a premise.

## 2. Literal verification of the three zeros

For a symbolic phase (tX+Y)^2+dX+eY, write M_r=Lambda_2(f^r).
Native expansion gives

 M_1=d+e+2(t^2+t+1).

After setting e=-2(t^2+t+1)-d, the next two native identities are

 M_2=2[d^2+(6t^2+2t-2)d+12t^4+8t^3+4],                 (2)
 M_3=12(t^2+1)M_2+24[J(t)d+B(t)],                       (3)
 B(T)=12T^6+8T^5-8T^2+6=-D(T).

The polynomial numerators satisfy the exact identities

 D(T)+E(T)+2(T^2+T+1)J(T)=0,                            (4)
 D(T)^2+(6T^2+2T-2)D(T)J(T)
                 +(12T^4+8T^3+4)J(T)^2=4P(T).          (5)

All displayed identities follow by literal coefficient multiplication.
They involve only rational integer polynomials, not numerical
approximations or a bounded search.

Equation(4) gives M_1=0 for(1). Equations(2),(5) give
M_2=8P(t)/J(t)^2=0. Finally J(t)d+B(t)=D(t)+B(t)=0,
so (3) gives M_3=0.

## 3. Exact support and the fourth moment

The three quadratic coefficients t^2,2t,1 are nonzero and satisfy the
singular-Hessian equation (2t)^2=4t^2.

Neither linear coefficient vanishes. For a direct proof, suppose d=0.
Equations(2)--(3), together with M_2=M_3=0, give

 K(t)=3t^4+2t^3+1=0,   B(t)=0.

But B(t)/2=2t^2 K(t)-6t^2+3. Hence t^2=1/2, and then
K(t)=t+7/4. This forces t=-7/4 and t^2=49/16, contradicting t^2=1/2.

If e=0, swap X,Y and divide the entire phase by t^2. This preserves
all three zero moment equations and gives the same normalized form
with t'=1/t and zero X-linear coefficient. The preceding argument,
which uses only native identities(2)--(3) and not P or J, again applies.
Thus e!=0 too. The actual support is exactly five.

The phase is nonzero because its Y^2 coefficient is1. Its first three
moments vanish, so the independently proved uniform four-window implies
M_4!=0. This proves Theorem1.

## 4. Scope

The example is algebraic and exact. It proves a three-window fails on
the zero-constant singular-Hessian family at start n=1, and does not
assert anything about the optimal window for all nonsingular quadratics,
about a smaller coefficient field, or about publication or priority.
It requires no computational algebra certificate as an unproved premise:
all necessary polynomial identities and denominator/support arguments
are displayed above.
