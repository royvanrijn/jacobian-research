# Polynomial-residue ramified boundary forces a line and closes ruling five

PF2RULED161. 17 November 2026. Actual complex polynomial PLANE Keller maps.

This note
now closes the ENTIRE quadratic ruling-five ACTUAL TARGET class, including
its arbitrary bad divisors and PF160's previously exceptional parameters.
The mechanism uses the ACTUAL Jacobian mate and the FULL finite normalization:
a forced ramified prime has residue a polynomial graph, forbidden by the
retained published Chau no-affine-line-COMPONENT theorem.

The universal theorem is stronger than this degree-five corollary. For an
odd ruling degree 2n+1, n>=1, the canonical even z-leading coefficient
b_n MUST vanish at the central multiple fibre. All lower coefficients,
all central pole orders and every companion degree are retained. This is
not an exclusion of all higher odd targets, all cyclic images or trivial
deck. Generic first-derivative integrality and HC4/JC2 remain OPEN.

## 1. Actual full-normalization setup and retained imported input

Let B=C[P,Q] subset R=C[x,t], F=Frac B, L=Frac R, with the actual
constant source Jacobian nonzero. Retain the FULL finite normalization S
of B in L, codifferent Lambda, and EVERY original ramified and unramified
boundary prime. Suppose there is an actual quadratic top factor L/E.
No generic-degree minimality is assumed.

[PF157](CANONICAL_DECK_JET_FIELDS_AND_RATIONAL_POLE_ATLAS_COVERAGE.md) covers EVERY such marked source image by

    V: w^2=N(u)+M(u)*z,
    N=u*Sigma(u)^2, M=u*A(u),
    u=t^2, w=t*(Sigma+A*x), z=2*Sigma*x+A*x^2.                (1)

Here A has order k>=1 at zero, Sigma(0)!=0, and gcd(A,Sigma)=1.
Put m=ord_0 M=k+1>=2. The ACTUAL target functions P,Q lie in O(V),
and their volume can be normalized by a nonzero constant so that

    dP wedge dQ=Omega=du wedge dw/M.                         (2)

The source pullback of Omega is -2 dx wedge dt. Rescaling a target
coordinate is legitimate and does not change B,F,S,Lambda or any prime.
Let S_E be the full finite normalization of B in E, not just O(V).

The imported global input is Chau, *Two remarks on non-zero constant
Jacobian polynomial maps of C2*, Ann. Polon. Math.82.1(2003),39--44,
Remark on printed43 / full PDF page5: a nonsingular polynomial map of
C2 has no irreducible exceptional-value component isomorphic to A1.
For a Keller map the critical-value set is empty, so this concerns its
NONPROPER-value components. Other components MAY be present.

This is the COMPONENT Remark already adopted/audited in PF112 and ME0,
not Corollary3 about the WHOLE exceptional curve and not a prohibition
on arbitrary cuspidal contractible curves. The retained six-page
[published primary PDF](sources/global-boundary-topology/chau-2003-two-remarks.pdf)
is hash pinned; the exact Remark was rechecked in that file and at the
[publisher](https://www.impan.pl/shop/en/publication/transaction/download/product/84523).
It is IMPORTED, not independently reproved. No new classification,
withdrawn stronger theorem or whole-exceptional-set interpretation is used.

## 2. Universal odd-coordinate hypothesis and the two simple branches

Let an actual target coordinate have ruling degree 2n+1, n>=1. Its unique
regular expression is

    P=w*A_P(u,z)+B_P(u,z),
    A_P=sum_(j=0)^n a_j(u)z^j, a_n!=0,
    B_P=sum_(j=0)^n b_j(u)z^j.                              (3)

Assume for contradiction

    b_n(0)!=0.                                             (4)

All coefficient polynomials have arbitrary u-degree. Put a=ord_0 a_n>=0.
On M!=0 substitute z=(w^2-N)/M, giving P* in C(u)[w]. With p an
independent target parameter define the polynomial

    F_n(u,w,p)=M^n*(P*(u,w)-p) in C[u,w,p].                  (5)

Its reduction at u=0 is

    F_n(0,w,p)=w^(2n)*(a_n(0)*w+b_n(0)),                    (6)

because N(0)=M(0)=0. Notice the entire dependence on p is divisible
by M^n. We use ordinary Hensel lifting over C[p], not only over C(p).

If a=0, (6) has the SIMPLE NONZERO root
w0=-b_n(0)/a_n(0). Its w derivative there is the nonzero CONSTANT
a_n(0)*w0^(2n). There is a unique lift

    w_hat in C[p][[u]], w_hat(0)=w0.                        (7)

If a>=1, set y=1/w and use the reciprocal polynomial

    G_n(u,y,p)=M^n*y^(2n+1)*(P*(u,1/y)-p).                  (8)

Its reduction at u=0 is b_n(0)*y. Thus y=0 is a simple root, with
nonzero CONSTANT derivative b_n(0), and its unique lift satisfies

    y_hat in u^a*C[p][[u]],
    leading(y_hat)=-leading(a_n)/b_n(0),
    w_hat=1/y_hat in u^-a*C[p][[u]].                        (9)

The order assertion follows from G_n(u,0,p)=a_n(u) and its unit y
coefficient. Inverting y_hat/u^a divides by a nonzero constant.

The all-order coefficient assertion is essential. Hensel recursion solves
EACH next coefficient by dividing by the same nonzero constant derivative;
therefore every coefficient remains in C[p]. No p-dependent denominator,
exceptional p value, uniform nilpotence or formal termination is assumed.
Earlier calculations of a finite prefix are not substituted for this proof.

## 3. The actual Jacobian makes the polynomial residue regular

Differentiate (5) or (8) at fixed u. Since the relevant derivatives are
units with p-independent constant term, the exact orders are

    ord_u(partial_p w_hat)=n*m+(2n-1)*a,
    leading(partial_p w_hat) in C*.                         (10)

For a=0 this follows from w_p=M^n/F_(n,w). For a>=1,
y_p=M^n*y^(2n+1)/G_(n,y) and w_p=-y_p/y^2. The leading y coefficient
in (9) is p-independent, so the leading constant in (10) is p-independent.

Every regular image function, including the ACTUAL arbitrary-degree
companion Q, belongs to C(u)[w] after the substitution in (1). Thus

    Q_hat=Q*(u,w_hat) in C[p]((u)).                          (11)

This Laurent series has a finite lower exponent and polynomial p
coefficients. Restricting the actual volume (2) to P=p gives

    dQ_hat/du=-partial_p(w_hat)/M.

Its right side has EXACT nonnegative order

    E0=(n-1)*m+(2n-1)*a>=0.                                (12)

A negative Laurent coefficient in Q_hat would differentiate to a pole,
which (12) forbids. Consequently

    Q_hat in C[p][[u]], q0(p)=Q_hat(0) in C[p],
    ord_u(Q_hat-q0(p))=e0=E0+1,
    leading(Q_hat-q0(p)) in C*.                            (13)

The constant q0 is a POLYNOMIAL, not a rational function with finite
poles. This is what makes the target change

    T=Q-q0(P)                                              (14)

a legitimate polynomial TARGET AUTOMORPHISM. It preserves the entire
original B,F,S,Lambda and every boundary prime. It is not an approximate
relation falsely declared a target coordinate.

## 4. The forced prime in the FULL finite normalization

Equations(7)--(9) embed E into C(p)((u)), with valuation ord_u, residue
field C(p), and P=p. Since P,Q are regular in this DVR, the full S_E
is contained in it by integrality. Its centre is a height-one prime D0
over the target graph

    Gamma={Q=q0(P)}={T=0}.                                 (15)

The centre is divisorial because P has transcendental residue. Normality
of S_E makes its localization the actual DVR, normalized by ord_u u=1.
The residue S_E/D0 is integral over C[P] and is a subring of C(P).
Since C[P] is integrally closed, this gives EXACTLY

    S_E/D0=C[P], f_(D0/Gamma)=1,
    e_(D0/Gamma)=(n-1)*m+(2n-1)*a+1.                        (16)

There is no deleted finite p-value and no punctured graph in (16).

Now pass to the ORIGINAL full S in L=E(t), t^2=u. The odd u valuation
forces one prime R0 over D0, with quadratic index TWO and residue degree
ONE. Its normalized valuation has ord_R0(t)=1. From (1),

    x=(w/t-Sigma)/A,
    q=-ord_R0(x)=2*(m+a)-1>0.                              (17)

The term w/t dominates Sigma in both cases; its nonzero leading
coefficient is independent of p. Thus R0 is OUTSIDE the whole source X.
It is an ORIGINAL omitted RAMIFIED prime, not a comparison curve on a
translated open. Its exact original target data are

    e_R=2*((n-1)*m+(2n-1)*a+1),
    f_R=1, S/R0=C[P], image(R0)=Gamma.                      (18)

All other primes over Gamma and every other target curve remain in S.
For n>=2, e_R>=6. The source-volume check is explicit:

    ord_R0(partial_p x)=2*e0-1=e_R-1,

as required by the actual normal ramification of dP wedge dQ. The large
p-independent pole of x has NOT been confused with this volume order.

Alternatively, for EACH p0 in C the unit-root construction converges
by the analytic implicit function theorem, giving w_hat(u,p0) or y_hat.
Put u=t^2 and x=(w_hat(t^2,p0)/t-Sigma)/A. Then x has the fixed pole
q from (17), while P=p0 and Q tends to q0(p0). These actual escaping
arcs exhibit EVERY point of the polynomial graph, independently confirming
that no finite endpoints were silently removed.

## 5. The global contradiction and complete ruling-five closure

The image Gamma of the original omitted prime R0 is an irreducible
component of the nonproper-value set of the actual plane Keller map.
It is the CLOSED polynomial graph in (15), hence algebraically smooth
and isomorphic to A1. The imported Chau COMPONENT Remark excludes it.
Other nonproper components do not affect this contradiction.

Therefore (4) is impossible. We have proved, for ALL actual quadratic
source images and EVERY odd target ruling degree2n+1 with n>=1,

    b_n(0)=0.                                              (19)

This does not exclude all higher odd targets: their surviving central
restriction may use b_(n-1),b_(n-2), etc. The coefficient a_n may have
arbitrary zeros anywhere; no leading coprimality or multiplicity bound
has been imposed. All lower coefficients and companion degrees remain.

For ruling degree FIVE, n=2. On the reached whole line t=0, w=0 and
z=2Sigma(0)*x is affine-linear. PF52's actual whole-line criterion makes
the restriction degree at least TWO for a noninvertible target, so
b_2(0)!=0. This contradicts (19). Hence

    EVERY actual quadratic-image target coordinate has ruling degree>=6.
                                                                    (20)

Degrees below four and degree four were already retired by PF157/PF158.
No generic-degree minimality is required. The ENTIRE PF159 bad-order
(15j,10j,5j) ruling-five residual and every PF160 exceptional parameter
as an ACTUAL quadratic-image target are now retired. PF160's stronger
rational-mate exclusion away its exceptional resultant stays useful;
(20) is NOT a no-rational-mate theorem for arbitrary rational companions.

For nontrivial FULL generic deck order g, PF157 and (20) give the new
reduced target coverage bound

    max(6,2g)<=r<s, neither r nor s divides the other.         (21)

In the quadratic case an odd ruling r=2n+1 has central source-line degree
at most n-1 by (19), and at least2 by PF52. Thus odd rulings now start
at SEVEN; a ruling-seven coordinate has central restriction degree TWO.
Both-even pairs start at SIX. No theorem here closes those survivors.

## 6. Codifferent ledger, scope guards and remaining full objective

For the forced original R0, the full codifferent threshold is 1-e_R.
PF141 gives a generic first source velocity pole q+e_R, with q from(17),
using the actual source volume. For n>=2 one has

    e_R-q=2*(n-2)*m+4*(n-1)*a+3>0,

so this prime contributes2 to the first-velocity target clearing exponent
1+ceil(q/e_R). The FULL exponent still takes the maximum over EVERY
original prime sharing Gamma. Neither other ramified primes nor the
unramified packet is dropped. This local ledger does not itself supply
marked quadratic regularity or derivative integrality; the class is
excluded by its global graph image instead.

Both simple-root hypotheses and n>=1 matter. On S2, P=w+g(u) has the
rational mate Q=1/u with dP wedge dQ=Omega. It has n=0; the coframe is
polar and (12) is false. When b_n(0)=0, (6) need not have any nonzero simple root and may even
vanish identically. No polynomial-residue conclusion is asserted in that
case. The regular ruling7 shape P=u^2*w*z^3+z^2+w has b3(0)=0 and
b2(0)!=0; its cleared polynomial vanishes identically at u=0. It
illustrates the central-degree boundary without asserting a Keller mate.

The actual dimension-two steps are the scalar relative coframe from the
2x2 Jacobian, the common whole-plane quadratic source chart, and Chau's
plane no-A1-component theorem. The known nonproper3D Keller control does
not satisfy these hypotheses. No dimension-free first-derivative lemma
from constant Jacobian, generic finiteness or derivative generation is
inferred.

The uniform proof is the all-order C[p] Hensel recursion, actual mate
identity(12), full finite-normalization residue argument(16)--(18), and
the imported component obstruction. Exact finite-prefix checks support
these arguments and do not replace them. Original S,Lambda and ALL
ramified primes remain; no new external classification or predecessor
checker replay is used. The checker passes 82 exact arithmetic and symbolic
controls. The retained published Chau input was rechecked
at primary sources, not independently reproved.

PF142's mixed minimal case, PF156's fixed-v2 quintic graph branch,
PF158's ruling-four branch and now the ENTIRE quadratic ruling-five
actual-target branch remain closed. Ruling>=6, higher cyclic images,
trivial deck, full marked quadratic regularity, generic derivative
integrality and normal-tail increasing ideal powers remain OPEN.
Full HC4/JC2 is not proved and no counterexample is constructed.
Independent whole-proof replay, formal verification and external review
remain unclaimed.
