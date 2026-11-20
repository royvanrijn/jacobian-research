# A uniform finite purity prefix for every two-height starting profile

20 November 2026. Work over ANY characteristic-zero
field K in the sphere algebra S=K[x,y,t]/(t^2+2xy-1). Its angular
functional is the exact rational polynomial functional

    Lambda(F)=CT_x (1/2) integral_(-1)^1
                F(x,(1-t^2)/(2x),t) dt.

Fix an integer d>=1. Define the explicit integers

    N_d=floor((d+3)^2/4),
    A_d=floor((d+5)^2/4),
    B_d=3*A_d*(N_d-1)-4.                            (U1)

**Theorem AllStartUniformPurityPrefix.** For EVERY integer a>=0 and
EVERY regular two-height phase

    p=t^a(C(x)+t^2D(x))/x, C,D in K[x], C(0)+D(0)=0,

of least sphere degree at most d, the following are equivalent:

1. Lambda(p^m)=0 for every positive integer m;
2. Lambda(p^m)=0 for1<=m<=B_d.

This is a uniform degree-bounded purity test, NOT an assertion about
arbitrary height profiles or arbitrary Image polynomials. It gives
no ordinary moment-ideal stabilization or nilpotence exponent bound.
The identically zero phase is pure and is included trivially.

## 1. Imported odd-start gcd/resultant bound

The companion odd-start finite-gcd proof is the only new input to
the odd case. Its setup requires no rationality theorem for the
selected square root or selected coordinate BEFORE purity.
Normalize a unit center C(0)=-D(0)=1 and suppose a is odd. Put

    f=x^2*(-D)^a/C^(a+2),
    W=(-D/C)^((a-1)/2)/C,
    H=2W/f'.

These are rational, f=x^2+O(x^3),W=1+O(x),H=1/x+O(1).
If n=deg f and h=deg H as rational MAP degrees in reduced form,
that theorem gives the following exact implication:

    if p is NOT pure, its first nonzero EVEN moment has order
    m <= (n-1)*(n+4h)-4.                            (U2)

The proof constructs the unique off-diagonal local branch f(j)=f
through j=-x+O(x^2), translates purity into a rational polynomial
identity on it, removes the FULL common gcd from BOTH polynomials,
and applies a nonzero residual Sylvester resultant. The common gcd
can contain components away from the origin even in a nonpure case;
those components are not omitted from this premise.

The retained formal endpoint formula remains valid before purity.
Only positive even m can be nonzero when a is odd. Its nonzero
rational endpoint factor identifies that first nonzero moment with
the first nonzero even coefficient of the formal endpoint series.

## 2. Explicit odd-height degree bounds

Let c=deg C,e=deg D, and write a=2q+1. Define

    E=2CD+a*x*C*D'-(a+2)*x*C'*D.

It is nonzero because f' is nonzero. Direct differentiation gives

    H=-2*C^(q+3)/(x*(-D)^q*E).                      (U3)

No canceled factor is being assumed absent. This formula gives valid
UPPER bounds even when numerator and denominator share factors:

    n <= max(2+a*e,(a+2)*c),
    h <= max(((a+5)/2)*c, 1+c+((a+1)/2)*e).         (U4)

Here deg E<=c+e was used. In the unit branch D is nonzero, and the
exact least-degree formula gives

    c<=d-a+1=L, e<=d-a-1=L-2,
    a<=d-1, L>=2.

Consequently

    n<=(a+2)*L,
    h<=((a+5)/2)*L.

For the h-bound its second competitor is smaller by d+1; the
n-bound follows similarly, with difference2d. As n>=2,h>=1, the
right side of(U2) increases with either variable. Hence a nonpure
odd-start unit phase has first bad even order at most

    3*(a+4)*L*((a+2)*L-1)-4.                       (U5)

The two elementary concave-product bounds

    (a+2)*L <= floor((d+3)^2/4)=N_d,
    (a+4)*L <= floor((d+5)^2/4)=A_d

give(U5)<=B_d. These are integer upper bounds, not optimization
over root positions or coefficient values. In particular they are
valid over arbitrary K and for an initially nonrational selected X.

## 3. Every even height, including zero

Let a>=0 be even and normalize a unit center as above. The formal
selected X satisfies X^2=-D/C, and

    u=x*(-D/C)^(a/2)/C in K(x), u=x+O(x^2),
    W=X^(a-1)/C,
    W^2=(-D/C)^(a-1)/C^2 in K(x),
    T=W/u'=1+O(x).

The endpoint formula has a nonzero rational factor at EVERY positive
order m in this parity, so purity is equivalent to T=1. Equivalently

    W^2=u'^2.                                      (U6)

Indeed T^2=1 and T(0)=1 force T=1 in the characteristic-zero formal
power-series field. If nonpure, the first nonzero moment has order
m=ord_x(T-1)=ord_x(T^2-1), since T+1 is a unit. Formal substitution
by u=x+O(x^2) preserves this order.

Put n0=deg u and w0=deg W^2 in reduced rational form. Differentiation
gives deg u'<=2n0, so deg u'^2<=4n0. Clearing the unit denominators
at0 in W^2-u'^2 bounds its nonzero numerator degree by w0+4n0.
Consequently

    m<=w0+4n0.                                     (U7)

For even a>=2 and L=d-a+1, the coefficient degree bounds give

    n0<=((a+2)/2)*L,
    w0<=(a+1)*L,
    m<=(3a+5)*L
      <=(d+2)^2/4+(d+3)^2/2.                       (U8)

The last bound separately maximizes (a+1)L and2(a+2)L, so it is
deliberately coarse. No rationality of X is required before testing
the rational identity(U6).

For a=0 the negative exponent in W^2 is interpreted rationally:

    u=x/C, W^2=-1/(CD).

Thus n0<=d+1,w0<=2d and(U7) gives m<=6d+4.

For d>=1 one has N_d>=d+3 and A_d>=2d+7. These follow by expanding
the squares, with remainders (d+3)(d-1)/4 and(d-1)(d+3)/4,
respectively. Hence

    B_d>=6d^2+33d+38,

which dominates both6d+4 and the even bound(U8), as well as4.

## 4. Zero centers and completion of the iff

If C(0)=D(0)=0, the only angular-weight-zero part is
t^a(c1+d1*t^2), c1=C'(0),d1=D'(0). The retained all-height
second/fourth moment theorem proves

    Lambda(p^2)=Lambda(p^4)=0 iff c1=d1=0.

Its exact cleared Sylvester determinant is4096 times a polynomial
in the integer a with strictly positive integer coefficients, so
it remains nonzero over every characteristic-zero K. Thus those
two moments force x^2 to divide C,D. Conversely that divisibility
gives strictly positive angular weights and ALL moments vanish.
Since B_d>=4, the prefix contains both required tests.

Every nonzero regular phase is in exactly one of the unit and
zero-center branches. If it were nonpure, Sections2--4 would give
a nonzero moment within the first B_d orders. This proves the
nontrivial direction of Theorem AllStartUniformPurityPrefix. The
other direction is immediate.

## 5. Reduced moment-ideal stabilization, with exact scope

Fix a and a finite coefficient box whose regular members have least
sphere degree at most d. Let R be its universal Q-coefficient ring
with the linear regularity relation C(0)+D(0)=0 imposed. The exact
moments M_m=Lambda(p^m) are polynomials in R. Put

    I_B=(M_1,...,M_(B_d)),
    I_infinity=(M_m : m>=1).

**Corollary UniformReducedMomentLocus.** One has

    radical(I_B)=radical(I_infinity).               (U9)

Here is a direct algebraic proof using the all-characteristic-zero-
fields scope, with no Nullstellensatz or field-embedding premise.
If some M_j did not lie in radical(I_B), the prime characterization
of a radical would give a prime ideal P containing I_B but not M_j.
The fraction field F=Frac(R/P) has characteristic zero because R is
a Q-algebra. The universal coefficient point over F would have its
first B_d moments zero but M_j nonzero, contrary to the theorem.
Thus I_infinity is contained in radical(I_B), and the reverse radical
inclusion follows from I_B being contained in I_infinity.

Conversely, radical equality implies the universal field-point prefix\ncriterion: for every higher moment M_j, some positive power M_j^e\nlies in I_B. At a field-valued coefficient point annihilating I_B,\nthat power is zero only when M_j is zero. Thus the reduced-ideal\nidentity and the universal field-point criterion imply one another.\n\nThis states a REDUCED, set-theoretic finite moment locus. It does NOT
assert I_B=I_infinity, finite generation by the same prefix as ordinary
ideals, any effective exponent M_j^e in I_B, or a bound on scheme
multiplicities. The coefficient box is fixed; there is no finite
parameterization claim for arbitrary polynomial degrees.

## 6. Degree dependence is necessary

**Proposition DelayedTwoHeightNonpurity.** Fix any integer a>=0, any
positive EVEN integer N, and any nonzero epsilon in a characteristic-
zero field K. The regular unit-center phase

    C=1, D=-1+epsilon*x^N,
    p=t^a(C+t^2D)/x
      =2t^a y+epsilon*t^(a+2)*x^(N-1)

has least sphere degree d=N+a+1. Its FIRST nonzero positive angular
moment occurs at exponent N, with exact value

    Lambda(p^N)
      =epsilon*N*integral_0^1 t^(aN+2)*(1-t^2)^(N-1) dt
      =epsilon*2^(N-1)*N!
           /product_(j=0)^(N-1)(aN+3+2j) !=0.       (DN1)

This is exactly linear in epsilon at that exponent. In particular no
degree-INDEPENDENT finite purity prefix can recognize all these
profiles, even at any one fixed starting height. A prefix recognizing
every profile of degree at most d must be at least N=d-a-1 along
this subsequence. No optimality claim for the quartic upper bound
is made.

**Proof.** Give x weight1,y weight-1,t weight0. The two summands of
p have weights-1 and N-1. The term in p^m using j perturbation
summands has weight

    j*(N-1)-(m-j)=j*N-m.

For1<=m<N this is never0 for an integer0<=j<=m, so every angular
moment vanishes. At m=N it is0 only for j1. Consequently its entire
weight-zero contribution is

    epsilon*N*t^(aN+2)*(1-t^2)^(N-1).

Since N is even, aN+2 is even. The exact angular integral is(DN1).
The rational coefficient is nonzero: every factor in its closed
form is a positive integer, hence nonzero in K. The closed form
also follows by expanding the binomial and integrating termwise,
or by repeated elementary integration by parts. The exact degree
formula max(degC+a-1,degD+a+1) gives N+a+1. QED.

These examples are NOT pure and are not claimed as counterexamples
to eventual safety under a purity hypothesis. Their role is solely
to demonstrate the necessity of degree dependence and the exact
meaning of a finite-prefix certificate.

## 7. Sources and verification boundary

Inputs are the companion odd-start gcd/resultant proof, the retained
formal endpoint identity for all heights, exact least-degree formula
and all-height zero-center second/fourth-moment theorem. The even
rational identity and all explicit degree bounds are proved here.
Classical polynomial-segment moment multiplicity bounds are a related
subject, but none is being applied to this different phase functional.
No historical-priority, external-human-review or formal-verification
claim is made.

The [uniform checker](verify_all_start_two_height_effective_purity_prefix.py),
[source ledger](all_start_two_height_effective_purity_prefix_source_ledger.json)
and [receipt](all_start_two_height_effective_purity_prefix_verification.json)
retain the rational and degree-bound controls. A separate
[delayed-moment checker](verify_two_height_delayed_nonpurity.py) and
[receipt](two_height_delayed_nonpurity_verification.json) test Section6.
