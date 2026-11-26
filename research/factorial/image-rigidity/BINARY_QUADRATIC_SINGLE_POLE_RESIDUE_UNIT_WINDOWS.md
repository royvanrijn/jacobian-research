# Positive unit windows for a binary matching residue with one pole

This note classifies the window length of the ACTUAL native binary residue series, not arbitrary degree-four rational series. It gives a small-ramification constant-shift corollary and isolates the cubic pole profile. No exact moment counterexample or full binary-quadratic closure is asserted.

Work in a field K of characteristic2. Write

 N(T)=1+(d+e+beta)T+de T^2,
 D(T)=(1+alpha T+d^2T^2)(1+gamma T+e^2T^2),
 R(T)=N(T)/D(T)=sum_(j>=0) R_j T^j.                       (1)

Assume R is NONPOLYNOMIAL and has exactly ONE distinct geometric pole. There is then a nonzero lambda in an algebraic closure such that its REDUCED denominator is (1+lambda T)^r,1<=r<=4. The exceptional profile below has lambda in K. No general pole-parameter membership assertion is needed: passing to an algebraic closure does not change which coefficients are zero.

**Proposition.** Every block of TWO consecutive positive coefficients of R contains a nonzero coefficient, except when

                         R(T)=1/(1+lambda T)^3.           (2)

In case(2), the coefficient pattern is

 R_(4h)=lambda^(4h), R_(4h+1)=lambda^(4h+1),
 R_(4h+2)=R_(4h+3)=0,

so the least positive coefficient window is exactly THREE. The exceptional native residue profile is, up to swapping the variables,

              alpha=e=0, gamma=d=beta=lambda!=0.         (3)

## 1. Native degree bounds after cancellation

Put d0=degD,n0=degN and let c be the degree of gcd(N,D), so the reduced denominator order is r=d0-c and its numerator degree is h=n0-c. Constants of both reduced polynomials can be normalized to1.

If1<=d0<=3, then de=0: if both d and e were nonzero, D would have degree4. Hence n0<=1. If d0=4, then de!=0 and n0=2 EXACTLY. Therefore:

* r=1 implies h<=1;
* r=2 implies h<=1;
* r=3 implies h<=1;
* r=4 implies c=0,d0=4,n0=2.

For r=1 these follow from d0=1,c=0 or d0=2,c=1; larger c would exceed n0. For r=2 the only additional possible case is d0=4,c=2,h=0. For r=3 one may have d0=3,c=0,h<=1 or d0=4,c=1,h=1. These observations explicitly allow raw poles canceled by N. It would be incorrect to assume the unreduced D has only the one surviving pole before this case analysis.

## 2. Each reduced order

For r=1 write R=(1+aT)/(1+lambda T). Polynomiality occurs precisely at a=lambda and is excluded. At every positive index j,

                     R_j=(a+lambda)lambda^(j-1)!=0.

For r=2 write R=(1+aT)/(1+lambda T)^2. Since the denominator is1+lambda^2T^2, every even coefficient R_(2h)=lambda^(2h) is nonzero. A positive two-window therefore suffices, whether or not a is zero.

For r=3 write R=(1+aT)/(1+lambda T)^3 and multiply numerator and denominator by1+lambda T:

 R=[1+(a+lambda)T+a lambda T^2]/[1+lambda^4T^4].             (4)

If a!=0, the coefficients in residues0 and2 modulo4 are nonzero, giving a positive two-window. If a=0, (4) gives exactly(2) and its asserted pattern, hence a three-window and no two-window.

For r=4, there is no cancellation. Every geometric root of each quadratic factor of D is the same pole. Both d,e are nonzero. Each constant-one quadratic must therefore equal (1+lambda T)^2. Consequently alpha=gamma=0 and d=e=lambda (squaring is injective in a field of characteristic2). In particular N_2=de=lambda^2!=0, and

                 R=[1+aT+lambda^2T^2]/[1+lambda^4T^4].    (5)

Again residues0 and2 have nonzero coefficients, so a two-window suffices. This proves the proposition's window assertion.

## 3. The exceptional profile is exactly the one-axis mixed residue

If the reduced R is(2), it cannot arise from d0=4,c=1: in that case n0=2 and the reduced numerator has degree1, while(2) has constant numerator. Thus d0=3,c=0,N=1 and D=(1+lambda T)^3.

Degree3 implies exactly one of d,e is nonzero and the opposite diagonal coefficient is nonzero. Swap variables so d!=0,e=0,gamma!=0. The native denominator is

                   (1+alpha T+d^2T^2)(1+gamma T).

The two factors must respectively be (1+lambda T)^2 and1+lambda T. Thus gamma=lambda,alpha=0,d=lambda. The identity N=1 gives beta=d. This proves(3), and direct substitution proves the converse.

When an odd-order reduced denominator is present, lambda belongs to K because its coefficient of T is r*lambda=lambda. In the even cases r=2 or4 the asserted even coefficients contain lambda^2 or lambda^4 in K; although a purely inseparable extension could be needed to name lambda in an arbitrary imperfect K, no such membership is needed for coefficient nonvanishing. In the native r=4 case the preceding factor argument gives lambda=d in K. In the exceptional r=3 case(3) also gives lambda=d in K. Accordingly ONLY these native-profile uses require the membership statement; the general even reduced-order window proof remains valid over the algebraic closure without asserting that every pole parameter lies in K.

## 4. Consequence on a coefficient-value grid with denominator twice-even

Use the [discrete-grid multiplicity theorem](QUADRATIC_CONSTANT_APPELL_MULTIPLICITIES_FROM_DISCRETE_VALUE_GRIDS.md): if the native integral coefficients R_j lie in a characteristic-zero valued field E with v(E^*) contained in (1/e_E)Z and v(2)=1, then at delta=1-2^(-a), a unit R_(m mod2^a) permits root multiplicity at most gcd(2^a,e_E). An assumed positive unit-window w therefore excludes a zero-moment block of length w+2^min(a,v_2(e_E)); all other constant valuations have bounds no larger than max(w,2).

Suppose v_2(e_E)=2 and the nonpolynomial native residue has exactly one geometric pole. If the residue is not(2), w<=2 and the critical multiplicity bound is at most4. Thus every constant shift has no SIX consecutive zero moments, at ALL positive starts. This gives actual-support SFC for dense six-term binary quadratics in that sector, not all sparse actual-support bounds.

For the exceptional cube(2), w=3 and the general transfer alone yields at most seven, so a six-window is not concluded. If a six-zero block at a critical value with a>=2 starts at n not equal to2 modulo4, one of its last two indices m=n+4,n+5 is0 or1 modulo4. Since2^a is divisible by4, R_(m mod2^a) is a unit. The block supplies P_m,P_(m-1),...,P_(m-4)=0, so its constant is a root of multiplicity at least5, contradicting the bound4. Therefore a surviving six-block in this cubic profile must start at

                              n=2 mod4.                  (6)

At a=1 the multiplicity bound is at most2 and the three-window already excludes a six-block. At a=2, the separate exact quarter-critical residue analysis sharpens(6) to n=2 mod8 and requires the critical constant to cancel the pole completely. No statement here proves a native zero block at any of these residual indices.

