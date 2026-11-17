# Complete weight(2,-1,-r) classification through degree seven (KCUB35)

Let K have characteristic zero. For every integer r>=4, every polynomial
Keller map of source and output weights (2,-1,-r), of ordinary degree at
most seven, is an automorphism. There is no support or invariant-flag
restriction. After diagonal output scaling to identity linear part,
the complete normal forms are as follows.

| r | normalized map |
|---|---|
| 4 | (a,b,c+b^4*(L+M*a*b^2)) |
| 6 | (a,b,c+L*b^6) |
| even r>=8 | identity |
| 5 | Z=c+M*b^5; F=(a,b+L*a^2*Z,Z), L,M arbitrary |
| 7 | (a,b+L*a^3*c,c+M*b^7), with L*M=0 |
| 9,11,13 | (a,b+L*a^((r-1)/2)*c,c) |
| odd r>=15 | identity |

The even inverses subtract the third-coordinate shear. In the odd cases,
recover b=Y-L*X^((r-1)/2)*Z, then c=Z-M*b^r, with M=0 where absent.
All displayed coefficient boundaries are included. No smaller counterexample
or unrestricted whole-profile exclusion follows.

The separate r=1,2,3 classes are already excluded by KCUB32, KCUB33 and
KCUB17 respectively. They are not premises of this r>=4 proof, and their
old UNKNOWN/timeouts remain untouched. This result closes the proposed
mixed-weight raw-sparsity route, including the two complete eleven-slot
models r=5,7, without changing a positive witness bound.

## Complete source models and half-cover identity

The weights are distinct, so the linear part is diagonal. Nonzero constant
Jacobian makes it invertible; output diagonal scaling normalizes it to
identity and the determinant to one. On a!=0 take a=s^2,b=u/s,c=v/s^r.
Write F=(s^2*f,s^-1*g,s^-r*h). Each b,c monomial fixes its a exponent
uniquely. Direct differentiation gives

    J=f*{g,h}+g*{f,h}/2+r*h*{g,f}/2=1,                  (1)
    {P,Q}=P_u*Q_v-P_v*Q_u.

The factor one half is essential. The finite checker reconstructs every
original source determinant and verifies descent of each quotient monomial:
a zero-weight monomial u^i*v^j becomes a^((i+r*j)/2)*b^i*c^j.
All exponents in the actual Jacobian are nonnegative integers. Thus this
fraction-field calculation is an identity of source polynomials, including
a=0, not a polynomial invariant-coordinate assumption.

The ten complete source models r=4..13 have respectively
17,11,10,11,9,7,9,7,6,7 nonlinear atoms: 94 in total, with 142 native
Jacobian coefficient equations. These are full bases, with no support cap
or collision equation. The all-r module bounds below cover every larger r;
the theorem is not a finite ratio scan.

## Even r: boundary units and three diagonal coefficients

Write r=2t, t>=2, and use the polynomial invariants U=a*b^2,V=a^t*c.
Exponent congruences give the complete modules

    F1=a*f(U,V), F2=b*g(U,V),
    F3=c*h(U,V)+b^(2t)*k(U).

The last coefficient can be chosen independent of V because
V*b^(2t)=c*U^t. Every monomial in f,g,h satisfies
3*i+(t+1)*j<=6; every U^j in k satisfies 2t+3j<=7.
The three normalized constant terms are one.

At b=0 the Jacobian is

    g0(V)*[h0(V)*(f0(V)+t*V*f0'(V))+V*f0(V)*h0'(V)]=1. (2)

Here f0=f(0,V), and similarly for g0,h0. Off-block b derivatives vanish
because the first/third outputs are even in b and the second is b times
an even polynomial. Polynomial units force g0=1 and the bracket to one.
If degrees of f0,h0 are d,e, the leading scalar in that bracket is
1+t*d+e, a nonzero positive integer. Hence d=e=0 and f0=h0=1.
This removes every V-only coefficient without a parameter restriction.

For t>=3 the remaining f,g,h depend only on U, since U*V exceeds the
budget. For t=2 the full boundary-reduced model is

    f=1+A*U+B*U^2+C*U*V,
    g=1+D*U+E*U^2+G*U*V,
    h=1+H*U+I*U^2+K*U*V, k=L+M*U.

The three coefficients [U^j*V^j] of J-1 are

    4C+3G+2K,
    6CG+5CK+4GK,
    7CGK.

For the actual nonempty support of C,G,K, choose j equal to its size.
The jth coefficient is the product of all active variables times
1+3*1_(C active)+2*1_(G active)+1_(K active), a positive integer.
Thus it cannot vanish in a characteristic-zero field. The seven nonempty
supports give C=G=K=0, covering all coefficient-zero boundaries.

In every even case the determinant now factors as

    J=h(U)*[g(U)*(f(U)+U*f'(U))+2U*f(U)*g'(U)].          (3)

The extra third-coordinate shear does not affect it. Units force h=1.
The bracket has leading scalar 1+deg f+2deg g, so f=g=1. Only k remains.
Its degree budget gives exactly r4: k=L+M*U, r6: k=L, and r>=8: k=0.
This proves all even normal forms without an external theorem.

## Odd r: complete models and two tiny unit identities

Now r>=5 is odd. The half-cover makes f even and g,h odd under
(u,v)->(-u,-v). The complete degree-seven modules are

    f=1+mu*u^2+P*u^4+nu*u*v,
    g=u+B3*u^3+B5*u^5+(A0+A1*u^2)*v,
    h=(1+H2*u^2+H4*u^4)*v+N*u*v^2+M*u^r.              (4)

For r=5,7 all eleven nonlinear slots are present. For r=9,11,13,
nu=A1=N=M=0, leaving seven slots. For r>=15 also A0=0,
leaving only the six U-dependent slots. These cutoffs come from the
ordinary source degree, not quotient degree: nu's source monomial has
degree(r+7)/2, A0's has degree(r+1)/2, and A1,N have degree(r+7)/2.
The b^r term requires r<=7.

Only r5/r7 can have nu!=0. A diagonal conjugacy in c makes nu=1,
without taking roots or losing any nonzero coefficient. It preserves the
complete source class and identity linear part. For each complete nu1
model the retained native coefficient equations I_i satisfy

    sum_i multiplier_i*I_i=1.

Each certificate has just nine multiplier monomials, maximum degree one:
[r5 certificate](mixed-weight2/r5-nu1-unit.json) and
[r7 certificate](mixed-weight2/r7-nu1-unit.json).
The [independent checker](verify_mixed_weight2.py) reconstructs the two
full source maps and all 34 specialized coefficient equations, then
multiplies all 18 rational multiplier terms exactly. No primary Groebner
verdict is trusted. These identities exclude every nu!=0 map; no generic
extra parameter factor or coefficient-zero stratum remains.

## Odd univariate first quotient

We now have f=phi(u)=1+mu*u^2+P*u^4. Suppose it is nonconstant.
Put A=A0+A1*u^2, B=u+B3*u^3+B5*u^5, C=N*u,
D=1+H2*u^2+H4*u^4, E=M*u^r (or zero where unavailable).
Then g=A*v+B and h=C*v^2+D*v+E. The v^2,v,constant coefficients of (1) are

    (2-r)*phi'*A*C/2+phi*(2A'*C-A*C'),
    phi'*[(1-r)*A*D+2B*C]/2+phi*(A'*D+2B'*C-A*D'),
    phi'*(B*D-r*A*E)/2+phi*(B'*D-A*E').                 (5)

A multiple root of phi makes J zero, impossible. Thus phi is squarefree.
Every root rho is nonzero because phi(0)=1. If N!=0, C(rho)!=0.
At rho the first coefficient forces A(rho)=0, the second then forces
B(rho)=0, and the constant coefficient is zero, again impossible. Hence N=0.

If A=0, J=D*(phi*B'+phi'*B/2). The second factor has degree
 deg phi+deg B-1>=2, with nonzero leading scalar deg B+deg phi/2.
It cannot be a polynomial unit. If A!=0, the remaining v coefficient gives

    phi^((r-1)/2)*D=k*A

for a nonzero field constant k, by differentiating the rational quotient.
Since deg phi>=2 and r>=5, its left degree is at least four, while
 deg A<=2. This is impossible. The A=0 branch, repeated roots and every
source-coefficient boundary have been handled separately. The rational
constant-field assertion and all root arguments are in characteristic zero.
Consequently phi=1, so every normalized map has first output a.

## Odd constant-first family and exact normal forms

For f=1, equation (1) is the plane determinant {g,h}=1, but no plane
Jacobian theorem is needed. Its v^2 coefficient is exactly

    2A'*C-A*C'=N*(3A1*u^2-A0).

If N!=0, this forces A0=A1=0, so g=B. The v coefficient is then
2N*u*B', incompatible with the constant determinant one. Thus N=0.
The v coefficient now gives A'*D-A*D'=0. If A=0, the constant coefficient
B'*D=1 gives B=u,D=1. Otherwise D=k*A for a nonzero constant k, and

    A*(k*B'-E')=1

forces A constant. Normalize D(0)=1 to get D=1 in every case. Integration
using B(0)=E(0)=0 gives B=u+A*E. Therefore, writing A=L,E=M*u^r,

    F3=c+M*b^r, F2=b+L*a^((r-1)/2)*F3.

For r5 the product L*M is allowed: its source term has degree seven.
For r7 it exceeds the B degree budget, forcing L*M=0. For r9 and larger
M is absent. The source cross term L*a^((r-1)/2)*c is permitted only
through r13. This yields exactly the odd normal forms in the opening table.
Their explicit inverses prove sufficiency over the original coefficient
field and globally at a=0; no square-root variable enters these inverses.

## Reproduction and trust boundary

Run from this directory:

    python verify_mixed_weight2.py

The [receipt](mixed-weight2/validation.json) pins all ten complete models,
all three producers, arithmetic helper, two exact unit certificates and
this proof. The stdlib Fraction checker verifies all 94 source atoms,
142 native coefficient equations and ten full original determinants;
five complete even boundary factorizations; three even diagonal identities
and all seven active supports; eighteen leading-degree cases; both odd
units; three symbolic-r odd univariate identities and the constant-first
top identity; and all eleven complete generic inverse charts. The r7
charts L=0/M=0 cover its full product-zero locus.

The [primary model producer](mixed_weight2.py),
[structural producer](mixed_weight2_structure.py) and
[bounded Singular certificate producer](certify_mixed_weight2.py) retain
all models, source scripts, logs, hashes and terminal checkpoints under
[mixed-weight2](mixed-weight2/README.md). Both primary units completed
in less than a second; there is no new timeout or unresolved solver handle.
The exact units, rather than primary status labels, are the accepted evidence.

The all-r exponent modules and degree cutoffs, unit/leading-degree
implications, nonzero-support field split, roots, rational derivatives
and diagonal normalizations are written universal mathematics. Exact replay
checks their finite polynomial premises, not a formalization or external
independent full review. No external classification is imported.

All twelve remaining unrestricted strict profiles remain open; the full
partition stays 62 excluded / 12 remaining. Raw3D13, normalized3D[10,48],
identity-linear-part cubic12/36, free-linear-part cubic12/31+17 and
cubic-linear89 remain distinct unchanged positive bounds. Global degree,
support and special-form dimension minima are not settled by this result.
