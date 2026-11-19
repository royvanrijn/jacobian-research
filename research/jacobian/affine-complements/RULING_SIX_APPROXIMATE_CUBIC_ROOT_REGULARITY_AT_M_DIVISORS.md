# The cubic square coefficient is regular at every M-divisor

PF2RULED167. 17 November 2026. Actual complex polynomial plane Keller maps
in the complete quadratic ruling-six class of PF165.

Write its regular square completion as

    P=T^2+R4, T=(h*w+K)*z, K=A2/2,
    R4=d*z^2+a1*w*z+a0*w+b1*z+b0,
    d=b2-N*h^2-K^2, w^2=N+M*z.                            (1)

Here b3=M*h^2 and a2=2hK are the polynomial normal form proved
in [PF165](RULING_SIX_SIMPLE_ESCAPE_GRAPHS_AND_QUARTIC_FACE_REDUCTION.md).
All coefficients in(1) are polynomials in u. At EVERY root c of M,

    ord_c d >= ord_c h + ord_c M.                         (2)

Consequently the next cubic square coefficient

    B1=d/(2hM)                                            (3)

has no pole at any M-root, including the central root. It also has no
pole at an off-M h-root where K(c)!=0. At the actual A pole divisors,
the retained weighted face gives ord d=ord h, so B1 is regular there
too. This note leaves the off-M common zeros of h,K which carry no
coframe pole unresolved; it does not assert global polynomiality of B1.
The [global quartic-divisibility theorem](RULING_SIX_GLOBAL_QUARTIC_DIVISIBILITY_AND_CUBIC_REMAINDER.md)
closes those remaining coefficient poles and makes B1 polynomial.
No first-derivative integrality, target-field membership or JC2 closure
follows from the local statement.

## 1. Exact equation in a regular translated variable

Keep the original fields and full finite normalizations of PF165:

    B=C[P,Q], F=C(P,Q), E=C(u,w), L=C(x,t), u=t^2,
    N=u*Sigma^2, M=u*A, w=t*(Sigma+A*x),
    dP wedge dQ=du wedge dw/M.                             (4)

At a noncentral M-root, N(c)!=0 and both reached values
w0=+/-sqrt(N(c)) satisfy b2(c)+2h(c)K(c)w0!=0. At the central
M-root N(0)=0 and b2(0)!=0. Every actual M-root has multiplicity
at least two, though the argument below only needs it positive.

Put U=h*w+K and C(U)=(U-K)^2-N*h^2. Substitution in the WHOLE
polynomial P, with all lower terms retained, gives

    F0(U)=(U^2+d)*C(U)^2
       +h*M*a1*(U-K)*C(U)+h^2*M*b1*C(U)
       +h^3*M^2*a0*(U-K)+h^4*M^2*b0,
    F0(U)=p*h^4*M^2 on P=p.                               (5)

Thus F0 is polynomial in U over C[u], not an initial substitute.
On a generic P fibre, h and K are p-independent, and its coframe is

    dQ=-(partial_p U)*du/(hM).                             (6)

These identities allow poles in w but use a translated U with finite
limit at the points constructed below.

## 2. A constant simple branch gives a forbidden full graph

Fix an M-root c and write v=u-c, j=ord_c h>=0 and m=ord_c M>=1.
Suppose d(c)!=0. Modulo v, every lower term after the first line of(5)
vanishes, so the equation reduces to

    (U^2+d(c))*C_c(U)^2=0.                                (7)

Both nonzero roots U0 of U^2+d(c) avoid C_c=0. If h(c)!=0,
a simultaneous root would give w0=(U0-K(c))/h(c) with w0^2=N(c),
and the identity

    U^2+d=h^2*(w^2-N)+2hK*w+b2

would contradict the reached-line unit condition in Section1.
If h(c)=0, then C_c(U)=(U-K(c))^2. A simultaneous root U0=K(c)
would give b2(c)=0, again impossible since a2(c)=0 and the reached
line has degree two. Thus each U0 is a SIMPLE root in(7).

All-order implicit lifting gives Uhat in C[p][[v]] with constant term
U0. The recursion divides only by the nonzero constant derivative in(7).
Equation(5) gives the exact orders

    ord_v(U_p)=4j+2m, ord_v(dQ)=3j+m>=0.                  (8)

The source escapes. If j>0, U0-K(c)!=0 and w=(U-K)/h has a pole.
If j=0, its finite limit has w0^2!=N(c), so z has a pole. For c!=0,
either original t choice is a unit and w0/t-Sigma(c) is nonzero,
so x is polar. At c=0, w0!=0 and the factor 1/t makes x polar.
If w is already polar the same original-source conclusion is immediate.

The actual Q is polynomial in w over C(u), so its evaluation lies in
C[p]((v)). The holomorphic differential in(8) forbids all negative
Laurent powers; its residue is q0(p) in C[p]. The valuation has residue
C(p), and its centre on the original full finite normalization S_E
has closed image Q=q0(P): its residue ring is integral over C[P]
and lies inside C(P), hence equals C[P]. Passing to the original source
field t^2=u preserves that residue; at c=0 use the required quadratic
uniformizer. All resulting source primes are omitted because x or t
is polar. The graph is consequently an irreducible component of the
actual nonproper-value set, excluded by the exact Chau component
remark retained in PF165 Section1. No finite graph point is deleted.

This contradiction proves d(c)=0.

## 3. The whole reached-line condition makes C_c(0) a unit

After d(c)=0, one has

    b2(c)=N(c)*h(c)^2+K(c)^2.                              (9)

If j>0 then b2(c)!=0 and K(c)!=0, so C_c(0)=K(c)^2!=0.
If j=0 and c=0, N(0)=0 and b2(0)!=0 give the same conclusion.
If j=0 and c!=0, write n0^2=N(c). The two reached-line units are

    b2(c)+/-2h(c)K(c)n0=(h(c)n0+/-K(c))^2.

They imply K(c)^2-N(c)h(c)^2!=0. Thus in EVERY case

    C_c(0)!=0.                                           (10)

This is where retaining both noncentral source lines, and separately
the central source line, is essential.

## 4. Every insufficient vanishing order gives two simple Puiseux graphs

If d is identically zero, (2) is immediate. Otherwise put

    delta=ord_c d>=1.

Suppose delta<j+m. Choose alpha=1 if delta is even and alpha=2
if it is odd, put r=alpha*delta/2, and substitute

    v=s^alpha, U=s^r Y.

Divide(5) by s^(alpha*delta). Its ENTIRE reduction is

    C_c(0)^2*(Y^2+d_delta), d_delta!=0.                    (11)

Every lower term has order at least alpha*(j+m)>alpha*delta;
the parameter term has order alpha*(4j+2m)>alpha*delta.
The two nonzero roots in(11) are simple constants. Hensel recursion
therefore gives Yhat in C[p][[s]], retaining every p-value.
For odd delta this is a normalized quadratic Puiseux branch: u-c
has order two and U has odd order delta, so their valuations have
gcd one. For even delta u-c=s itself is a uniformizer.

Along either branch F0_U has EXACT order r, with leading coefficient
2*C_c(0)^2*Y(0), a nonzero constant. Differentiating(5) and using(6)
therefore gives

    ord_s(U_p)=alpha*(4j+2m)-r,
    ord_s(dQ)=alpha*(3j+m+1)-r-1
             =alpha*(3j+m+1-delta/2)-1>=0.                (12)

For example the final inequality follows from delta<=j+m-1,
which bounds its parenthesis below by (5j+m+3)/2>=2.

The original source again escapes. If j>0, (10) makes K(c) a unit,
so w tends to infinity as -K/h. If j=0, w tends to -K(c)/h(c),
whose square differs from N(c) by(10). The source-lift argument of
Section2 therefore applies, including the central t^2=u cover.
The actual companion evaluates in C[p]((s)); (12) forces its residue
to be polynomial in p. The full-normalization graph argument of
Section2 excludes both branches. This proves delta>=j+m and hence(2).

## 5. Off-M unit K and the known A pole divisors

At an off-M zero of h with K(c)!=0, take j>=1,m=0. If d(c)!=0,
the two roots of U^2+d(c) are opposite nonzero constants; at most one
is K(c), so choose one with C_c(U0)=(U0-K(c))^2!=0. It gives the
same simple graph as Section2, with w polar and coframe order3j.
Thus d(c)=0 and C_c(0)=K(c)^2 is a unit. The proof(11)--(12)
now gives ord_c d>=j. Hence B1 is regular there as well.

At an actual A pole divisor off M, put its slope L0 and nonzero
weighted root c0. The retained initial is

    I=C0*Y^2*(Y-c0)^4.

Its coefficients in degrees6,5,4 are respectively
C0,-4c0C0,6c0^2C0. If h=v^(2L0)*h0+... and M(c)=M0,
comparison gives K=-2c0*h0*v^L0+... and
b2=6c0^2*h0^2*v^(2L0)+.... Therefore

    d=2c0^2*h0^2*v^(2L0)+..., ord_c d=2L0=ord_c h.        (13)

In particular these ACTUAL pole divisors do not cause a pole in B1.
Only off-M common zeros of h,K carrying no coframe pole remain for
a global polynomiality argument. This is a restriction on the actual
remaining class, not a classification or construction of that class.

## 6. Verification and scope

The checker controls the exact translated sextic, the reached-line
unit identities, both constant-root cases, normalized Puiseux parity,
complete symbolic valuation margins and the coframe orders. All-order
lifting, normal full images and the component theorem establish the
unrestricted result; finite checks do not replace those proofs.

Replay with

    .venv/bin/python research/jacobian/affine-complements/verify_ruling_six_M_approximate_root.py

The original target, full finite normalizations, both noncentral reached
lines, central source line and all original primes remain in scope.
The imported Chau/one-line assertions retain their existing assurance
limits. No new external theorem is imported and no independent whole-
proof review, formal verification or external human review is claimed.
[PF188](RULING_SIX_SQUARE_DEFORMATION_AND_MARKED_QUADRATIC_EXCLUSION.md) excludes the retained marked quadratic ruling-six class.
Derivative integrality and JC2 remain open.
