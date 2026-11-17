# Complete opposite-weight classification through degree seven (KCUB33)

Let K be a field of characteristic zero. Fix coprime positive integers
m,q with m>=2. Every polynomial Keller map in three variables, equivariant
for source and output weights (m,-q,-m), of ordinary degree at most seven
is an automorphism. After diagonal output normalization it is precisely

    F(a,b,c) = (a, b, c + a^(q-1)*b^m*k(a^q*b^m)),

where every nonzero monomial k_j*v^j satisfies

    m+q-1 + (m+q)*j <= 7.

The inverse subtracts the displayed shear. There is no support or invariant
profile restriction. This is a classification of this weighted family,
not a global Keller theorem or a new counterexample. The condition m>=2
is essential to the argument: m=1 includes the known weight(1,-1,-2)
counterexample after permuting the negative coordinates. No claim is made
here about that exceptional family.

## Complete weighted modules

The three weights are distinct because gcd(m,q)=1 and m>=2. Equivariance
therefore makes the linear part diagonal; the nonzero constant Jacobian
makes its three entries nonzero. Scale the outputs so that the linear
part is the identity and the Jacobian is one.

Put u=ac, v=a^q*b^m, s=m+q. An invariant monomial has b exponent divisible
by m and a exponent equal to the c exponent plus q times that quotient.
Consequently the invariant ring is K[u,v]. More generally the weight-m
module is a*K[u,v]. For weight -q, the b exponent is congruent to 1 modulo
m, hence is at least 1; its module is b*K[u,v]. For weight -m, the b
exponent is divisible by m. If the c exponent is positive the monomial
lies in c*K[u,v]. Otherwise it is a^(q-1)*b^m times a power of v.
The relation u*a^(q-1)*b^m=c*v lets us choose the latter coefficient
independent of u. Thus, uniquely,

    F1=a*f(u,v), F2=b*g(u,v),
    F3=c*h(u,v)+a^(q-1)*b^m*k(v),
    f(0,0)=g(0,0)=h(0,0)=1.

The degree budgets are 2*i+s*j<=6 for monomials u^i*v^j in f,g,h,
and s-1+s*j<=7 for v^j in k. These are complete modules, with no generic
coefficient or nonvanishing assumption. The argument fails at m=1
because weight -q can then have b exponent zero.

## The boundary removes every u-only coefficient

At b=0, since m>=2, the v dependence and the extra generator contribute
no off-block b derivative. Direct differentiation gives

    det JF |_(b=0) = g(u,0) * d/du [u*f(u,0)*h(u,0)] = 1.

Both factors are polynomial units in K[u]. Their values at zero are one,
so g(u,0)=1 and the derivative is one. In characteristic zero, integration
and the factor u imply f(u,0)*h(u,0)=1. Polynomial units and the constant
normalization now give f(u,0)=h(u,0)=1. Substitution u=ac is injective,
so this is an identity, not just a condition at finitely many points.

For s>=9, all v terms and the extra generator exceed degree seven; this
already proves F=identity. For every smaller s, the remaining coefficients
fit the following single universal model, with unavailable coefficients
set to zero:

    f=1+A*v+B*u*v+C*v^2,
    g=1+D*v+E*u*v+G*v^2,
    h=1+H*v+I*u*v+K*v^2,
    k=L+M*v.

Indeed s>=3, so v^3 and u^2*v exceed the budget. All primitive pairs
with 3<=s<=8 give fourteen complete models, with 169 nonlinear source
atoms before the boundary and 43 after it. This finite enumeration checks
the module calculation; the written degree argument covers all larger s.

## Three coefficients remove all mixed terms

Write {P,Q}=P_u*Q_v-P_v*Q_u and R=u*h+v*k. The Jacobian identity is

    J = g*(f*R_u+R*f_u+q*v*{R,f})
        +m*v*(f*{R,g}+R*{f,g}).

For example it follows by introducing U=F1*F3=f*R and
V=F1^q*F2^m=v*f^q*g^m. Chain rule in (a,u,v) coordinates gives
Jac_(u,v)(U,V)=f^q*g^(m-1)*det JF; expanding and cancelling gives the
polynomial identity above. The coordinate calculation is made in the
fraction field and yields a polynomial identity, so no a=0, b=0 or
parameter boundary is discarded. The independent checker also verifies
the identity against the actual source Jacobian in all fourteen models.

In the universal model the three diagonal coefficients of J-1 are

    [u*v]     = (q+2)*B+(m+1)*E+2*I,
    [u^2*v^2] = (m+q+2)*B*E+(q+3)*B*I+(m+2)*E*I,
    [u^3*v^3] = (m+q+3)*B*E*I.

If some of B,E,I are nonzero, take their actual nonempty support and
let n be its size. Set all inactive coefficients to zero in the nth
identity. It becomes the product of all active coefficients times

    1+(q+1)*1_(B active)+m*1_(E active)+1_(I active).

This scalar is a positive integer and hence nonzero in K. A field has no
zero divisors, so the identity cannot vanish. The seven possible nonempty
supports are all covered, including every coefficient-zero boundary.
It follows that B=E=I=0. The field hypothesis matters; this is not an
argument about nilpotent coefficient rings.

## Univariate degree finishes the classification

Now f,g,h depend only on v, and the exact Jacobian factorization is

    J=h(v)*[g(v)*(f(v)+q*v*f'(v))+m*v*f(v)*g'(v)] = 1.

Polynomial units force h=1 and the bracket to be one. If r=deg f and
t=deg g, the bracket has leading coefficient

    (1+q*r+m*t)*lc(f)*lc(g).

The integer factor is nonzero for every r,t>=0. Hence its degree is r+t,
which must be zero. The normalized constants give f=g=h=1, leaving k
arbitrary within the degree budget. This proves the stated triangular
form and its explicit inverse over K and every characteristic-zero
extension. No plane Jacobian theorem or external classification is used.

## Reproduction and trust boundary

Run from this directory:

    python verify_opposite_weight_reflection.py

The [independent checker](verify_opposite_weight_reflection.py) uses only
stdlib rational arithmetic and the existing [sparse polynomial helper](verify.py).
It reconstructs every weighted source monomial for the fourteen finite
models, all boundary factorizations and actual source determinants. It
replays the three universal symbolic-weight coefficients, all seven
support masks, the pure-v factorization and all nine possible degree
pairs, and checks both inverse compositions of fourteen non-affine
triangular controls. Its [receipt](opposite-weight-reflection/validation.json)
pins the models, producers, arithmetic helper and this proof.

The [universal producer](opposite_weight_structure.py) and
[finite-model producer](opposite_weight_reflection.py) use SymPy for
discovery. Five primary Singular unit certificates and the earlier
weight212 screen are retained as unused discovery history; neither the
primary Groebner answers nor their large multipliers are proof premises.
See the [artifact boundary](opposite-weight-reflection/README.md).

The module completeness, polynomial-unit implications, nonempty-support
field split and all-weight degree bound are written mathematics, not a
formal proof or an identified independent external review. Exact replay
checks their finite algebraic premises; it does not by itself formalize
the universal argument.

This closes only primitive weights (m,-q,-m), m>=2, through degree seven.
In particular it closes the proposed (2,-1,-2) gate without any coefficient
search being left unresolved. The m=1 exception, other unequal weight
ratios, unweighted constructions, all twelve remaining strict profiles
and global simplicity minima remain outside this theorem. Whole-profile
counts stay 62 excluded / 12 remaining. No positive witness bound changes.
