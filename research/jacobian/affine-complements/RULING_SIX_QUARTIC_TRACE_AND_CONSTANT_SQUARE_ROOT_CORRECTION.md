# A quartic trace forces the constant square-root correction to be regular

PF2RULED170. 17 November 2026. Actual complex polynomial plane Keller maps.

At every off-M zero of h whose outer face contains all four escaping
sheets, h divides the cubic remainder coefficient C1 locally.
Equivalently the constant correction C1/(2h) to the polynomial-part
square root is regular there. This includes EVERY surviving A pole
divisor, as well as an outer face with two distinct double roots.

A trace identity in an actual quartic factor proves the result.
A failed divisibility would force a root with holomorphic coframe
and affine first parameter response, giving a forbidden closed graph.

The same local divisibility holds at every M-root and at every
off-M h-root where K is a unit. The remaining unequal-slope
two-cluster loci and the full ruling-six class are not declared
closed here. T and the lower-ruling remainders
are not assumed to lie in the target field or full finite normalization.
The [nested two-cluster theorem](RULING_SIX_REGULAR_CUBIC_APPROXIMATE_ROOT_AND_QUADRATIC_REMAINDER.md)
resolves those remaining coefficient poles and gives a regular quadratic remainder.

## 1. Actual regular cubic after the linear correction

Use the complete quadratic source-image class, original fields and
full normalizations of PF165--169:

    w^2=N+M*z, M=u*A, N=u*Sigma^2, u=t^2,
    b3=M*h^2, a2=2h*K, h,K in C[u],
    r4=b2-N*h^2-K^2, B0=r4/(2hM) in C[u].                (1)

The regularity of B0 is the global conclusion of
[PF169](RULING_SIX_GLOBAL_QUARTIC_DIVISIBILITY_AND_CUBIC_REMAINDER.md).
Keep that result as an input; no new normalization replaces the
actual B=C[P,Q], S_E or S.

Define the regular cubic and cubic remainder by

    W1=h*w*z+K*z+B0*w,
    C1=a1-2B0*K,
    D2=b1-2B0*h*N-B0^2*M,
    D0=b0-B0^2*N,
    P=W1^2+C1*w*z+D2*z+a0*w+D0.                         (2)

All coefficients in (2) are polynomials in u, with arbitrary
degree. The actual companion has arbitrary degree and satisfies
dP wedge dQ=du wedge dw/M. The complete outer-square result
of PF168 is an indispensable input below.

Let c be an off-M zero of h whose outer face consumes four
escaping sheets. Put v=u-c. PF168 gives, with rational ell>0,

    ord_c h=j=2ell,
    in(P*)=(aY^3+bY^2+c1Y)^2, a*c1!=0, Y=v^ell*w,
    P*=v^(-2ell)*in(P*)+terms of greater v order.          (3)

The sign of the square root is chosen to agree with h/M.
The unit c1 need not come from a repeated root: both a quartic
root and two distinct double roots are included.

The source coframe on P=p is

    omega=dQ*=-du/(M*P*_w).                              (4)

## 2. Whole scaled equation near the four-sheet cluster

Assume C1 is not identically zero, and put a1ord=ord_c C1.
Matching the whole initial square in (3), including its cubic
and quadratic coefficients, gives

    a1ord>ell, D2(c)=0.                                  (5)

To see the matching, (1) gives the degree-three, degree-two and
degree-one coefficients of W1* as h/M, K/M and B0-hN/M.
PF168's coefficient identity r4/M^2=2a*c1*v^j+...
makes the last coefficient have initial c1. Therefore the
initial of v^ell*W1* is exactly the square root
aY^3+bY^2+c1Y in (3). In the difference (2), a C1 term at
order ell would contribute a nonzero cubic initial, and a unit
D2 would contribute a nonzero quadratic initial. Both must vanish.

Suppose for contradiction

    ell<a1ord<2ell.

All fractional exponents can be cleared by v=s^alpha; retain
the v-valuation notation until needed. Put

    epsilon=v^(-ell)*C1/M, e=ord_v epsilon=a1ord-ell,
    0<e<ell,
    F(v,Y,p)=v^(2ell)*(P*(c+v,v^(-ell)*Y)-p).             (6)

F has regular coefficients in the cleared parameter, and its
leading coefficient in Y is a unit. Define

    S(v,Y)=v^ell*W1*(c+v,v^(-ell)*Y).

Its constant coefficient is -v^ell*K*N/M and has order at
least 2ell. Therefore

    S=Y*A(v,Y)+O(v^(2ell)),                              (7)

where A is a regular quadratic polynomial in Y with unit leading
coefficient a2(v), and A(0,0)=c1!=0.

Substitution into the WHOLE identity (2) gives

    F=Y^2*(A^2+epsilon*Y+d)
                   +terms divisible by v^ell,
    d=D2/M in C[[v]], d(0)=0.                            (8)

Indeed the omitted terms of degrees one and zero have orders
at least ell and 2ell respectively. The parameter occurs ONLY
as the exact term -p*v^(2ell).

Normalize the unit leading coefficient a2(v)^2. Hensel factorization
separates the two roots at Y=0 from the other four roots:

    F/a2^2=F2*G, F2 monicdegree2, G monicdegree4,
    F2(0,Y)=Y^2,
    G(0,Y)=(A(0,Y)/a2(0))^2.                            (9)

The factors are coprime, because A(0,0)!=0. The recursion is over
C[p][[s]], dividing only by constant-unit resultants. Uniqueness
modulo v^ell gives

    G=(A/a2)^2+(epsilon/a2^2)*Y+d/a2^2
                                      modulo v^ell.      (10)

This remains valid when the four roots consist of two different
double clusters. No degree-two refinement is discarded.

## 3. The actual quartic trace detects a small derivative

For a monic quartic

    G=Y^4+bY^3+cY^2+dY+e0,

write its depressed linear coefficient as

    q(G)=d-b*c/2+b^3/8.

A square of a monic quadratic has q=0 identically. Adding a
constant changes no q, and adding epsilon*Y changes q by epsilon.
Thus (10) gives

    q(G)=epsilon/a2^2 modulo v^ell,
    ord_v q(G)=e<ell.                                   (11)

There is an exact multiplication/field-trace identity

    Tr(G'(Y))=-8*q(G)                                   (12)

in the actual rank-four algebra over the local coefficient field.
For completeness, Newton sums of the roots are
sum y=-b, sum y^2=b^2-2c, sum y^3=-b^3+3bc-3d.
Substitution in sum(4y^3+3by^2+2cy+d) proves (12).
The identity is not a trace in the original finite target
normalization or a nilpotence surrogate.

All four roots of G are regular and nonzero at v=0.
From (11),(12) and the valuation inequality for sums, at least
one root has

    mu=ord_v G'(Y)<=e<ell.                               (13)

At that root, F2 is a unit and G=0. Equations (9),(4) give

    ord_v P*_w=mu-ell<0,
    ord_v Y_p=2ell-mu>ell.                              (14)

The second equality also follows directly by differentiating the
WHOLE polynomial F0(v,Y)-p*v^(2ell). In particular the original
coframe is holomorphic. With v=s^alpha its exact order is

    alpha*(ell-mu)+alpha-1>=0.                           (15)

## 4. The first parameter response must be affine

Let r=ord_v Y_p>ell and let lambda be the first coefficient
of Y outside C, after its finite C-coefficient truncation Ypre.
The constants of the extended partial_p derivation are C, so r
is the exact index of that first coefficient. Normalize the
Puiseux parameter to make all indices integral.

Expand the p-INDEPENDENT whole polynomial F0 at Ypre.
Every coefficient and Ypre are regular, so every Taylor term
of degree at least two has order at least

    2r>2ell.

The derivative at the actual root has order mu, and the difference
from its value at Ypre has order at least r>mu. Hence the linear
term has exact order mu+r=2ell, with nonzero constant leading
slope. The constant term has order at least 2ell, since a smaller
constant coefficient could not cancel. The initial equation is

    A0*lambda+B0const=p, A0 in C*.                       (16)

Thus lambda is affine in p.

The entire polynomial-source companion Q*, evaluated on this
branch, is finite by (15). The full Taylor-leading argument of
PF159 makes its residue polynomial in lambda, hence polynomial
in p. Equivalently, (16) is a constant-unit Hensel equation,
so all later coefficients lie in C[p] and the finite companion
residue q0(p) is in C[p].

Since Y has nonzero limit, w=v^(-ell)*Y is polar. Both original
unit t lifts over c have x=(w/t-Sigma)/A polar. The ORIGINAL
full-normalization closed-image argument of PF165 therefore
gives the whole graph Q=q0(P) as an irreducible nonproper-value
component. The retained published Chau component theorem excludes it.

This contradiction proves

    ord_c C1>=2ell=ord_c h.                              (17)

The zero C1 case satisfies (17) trivially.

## 5. Unit-K roots, including all M-roots

The same conclusion holds whenever h vanishes but K is a unit,
including M-roots. Put j=ord_c h>=1, m=ord_c M>=0, and suppose
a=ord_c C1<j.

There is a p-independent exact zero wo of W1 in C((v)) with

    wo=-K/h+O(v^m), ord_v wo=-j,
    ord_v(W1_w(wo))=-j-m.                               (18)

For an all-order proof put U=hw+K and clear h^2*M in W1=0:

    U*((U-K)^2-N*h^2)+B0*h*M*(U-K)=0.

Its U derivative at v=U=0 is K(0)^2, a unit. The unique U
lift has order at least j+m, proving (18). The derivative
order follows directly from the cubic, and its leading coefficient
is nonzero independent of p.

At wo, the C1*w*z term in (2) has exact order a-3j-m.
It strictly dominates D2*z, a0*w and D0, whose orders are
at least -2j-m, -j and0 respectively. Substitute

    w=wo+v^((a-j+m)/2)*Y.

The correction has strictly larger order than wo. Scaling the
whole P*-p equation gives a nonzero constant plus a nonzero
constant multiple of Y^2. Its roots are simple and nonzero.
After clearing the possible denominator two, Hensel lifting
is over C[p] with constant-unit derivative. Exact orders give

    ord_v w_p=(5j+3m-a)/2,
    ord_v omega=(5j+m-a)/2>0.                            (19)

Thus the actual companion is finite with polynomial p residue,
and the original source x is polar. At c=0 use u=t^2;
at c!=0 use the unit t lifts. The same whole closed graph
contradiction proves ord_c C1>=j.

Every M-root where h vanishes has K a unit. Indeed PF169
gives ord r4>=j+m>0 and (1) gives b2(c)=K(c)^2 there.
The reached-line condition makes b2(c) nonzero, since
a2(c)=2h(c)K(c)=0. This also covers the central root,
where PF162 supplies b2(0)!=0. If h is a unit there,
the desired divisibility has no local condition.

Consequently C1/(2h) is regular at EVERY M-root and every
off-M unit-K h-root, as well as all four-sheet outer faces.

## 6. Scope and the remaining poles of the constant correction

The result (17) applies at every off-M h-root whose outer edge
contains all four escaping sheets. In particular it applies at
EVERY A pole divisor, with ell=L and j=2L. Thus the rational
constant correction

    C0=C1/(2h)

has no pole at any divisor supporting a pole of the actual
primitive Q on the generic P curve.

Any remaining pole of C0 can only lie at an off-M common
zero of h and K with specialization degree two and outer
escaping mass two. The other two escaping sheets have a
smaller pole slope. That unequal-slope case is outside the
four-sheet trace argument and the unit-K argument. No global polynomial
termination or target-field membership is inferred.

The checker verifies the full quartic trace in its actual
companion multiplication algebra, the monic-square linear
invariant, Hensel factor precision with constant cancellation,
and the exact derivative/response inequalities. The checker passes
107 exact controls. Whole-polynomial Taylor expansion and the full closed-graph theorem establish the
unbounded statement, not finite coefficient prefixes.

Replay with

    .venv/bin/python research/jacobian/affine-complements/verify_ruling_six_quartic_trace_correction.py

No new external classification is imported. The full original
normalizations and all their primes remain in scope. [PF188](RULING_SIX_SQUARE_DEFORMATION_AND_MARKED_QUADRATIC_EXCLUSION.md) excludes the retained marked quadratic ruling-six class.
Full derivative integrality, normal-tail termination and HC4/JC2 remain
open. Independent whole-proof
review, formal verification and external human review remain unclaimed.
