# Half-specialization pole bound excludes every quadratic ruling-four target

PF2RULED158. 17 November 2026. Actual complex polynomial PLANE Keller maps.

PF157 supplied complete quadratic source-image coverage. This note now
excludes the ENTIRE ruling-degree-four target-coordinate family in that
class, including ALL bad common leading divisors left by PF151/PF157.
There is no bound on the degree of the companion or on lower coefficients.
The argument is a new local generic-parameter bound using the DEGREE OF
THE SPECIALIZATION, followed by a complete generic-curve pole audit.

The result supplies no automatic integrality on the original full finite
normalization. Higher ruling degrees and the trivial-deck case remain
open; no plane Keller counterexample is constructed.

## 1. A reciprocal first-parameter-coefficient lemma

Put k=an algebraic closure of C(p), with its unique extension of d/dp.
Its constants are C, as proved in [PF149](MONIC_PARAMETER_POLE_BOUNDS_AND_EVEN_LEADING_CYCLIC_TARGETS.md), Section2.
Let A(s,X) in C[[s]][X] have degree D and let its reduction at s=0
have degree d, where 1<=d<D. Put e=D-d. Suppose

    A(s,w)=p, w in k((s)), nu(w)=-h<0.

Then

    nu(partial_p w) >= ceil(D*h/e)-2*h.                         (1)

The coefficient field C[[s]] is p-independent. No simplicity of the
escaping root, restriction on the vanishing of the leading coefficient,
or bound on the p-independent polar part of w is assumed.

Proof. Set Z=1/w and F(s,Z)=Z^D*A(s,Z^-1). The coefficient of Z^e
in F(0,Z) is a nonzero constant; all lower coefficients vanish there.
The equation becomes

    F(s,Z)=p*Z^D, nu(Z)=h>0.                                   (2)

There is a FIRST coefficient of Z outside C; otherwise (2) contradicts
p-transcendence. Write its index r and value lambda, and let Z0 be the
finite sum of earlier coefficients. Thus Z0 in s*C[s], possibly zero,
Z=Z0+s^r*(lambda+O(s)), and lambda not in C.

In the polynomial F(s,Z0+s^r*Y), the coefficient of Y^e has order EXACTLY
e*r: its quotient by s^(e*r) reduces to the nonzero Z^e coefficient
of F(0,Z). This uses Z0(0)=0. If q is the minimum order of all its
coefficients, q<=e*r. Its nonzero initial polynomial is H(Y) in C[Y].

If e*r<D*h, then q<D*h. Taking the initial coefficient of (2) gives
H(lambda)=0. Since C is algebraically closed this contradicts lambda
not in C. Hence r>=ceil(D*h/e). The constants of partial_p are C, so
nu(partial_p Z)=r and nu(partial_p w)=r-2h, proving (1).

In particular, if d>=D/2, the right side of (1) is NONNEGATIVE.
This is the decisive extra information beyond the unrestricted monic
coefficient bound in PF149. It comes from the nonzero LOW reciprocal
coefficient Z^(D-d), whose index counts the roots escaping the affine
chart, rather than from the full monic degree D.

## 2. A complete plane generic-fibre obstruction

Let P in C[t,x] have x-degree D>=2, leading coefficient a(t)!=0,
and b=deg_t a. Assume

    b>=D;
    deg_x P(c,x)>=ceil(D/2) at EVERY root c of a.                (3)

Then P has NO rational mate Q with nonzero constant Jacobian, even in
ANY finite characteristic-zero extension of C(t,x). All lower
coefficients and their t-degrees are arbitrary.

Work on every component of the WHOLE geometric generic curve P=p.
Its regular nonzero relative form for the source volume dx wedge dt is

    omega=(partial_p x)*dt=dt/P_x.                             (4)

The other residue patches give regularity at affine points where P_x=0.
The generic curve is smooth in characteristic zero; the argument does
not assume its geometric irreducibility. Every generic component has
nonconstant t. Indeed a whole vertical component would require all
positive-x coefficients to vanish at a CONSTANT common root, where
P could not equal transcendental p.

Take EVERY point on the projective normalization of each component.
If t is finite and x finite, it is an affine point. If t=c finite and
x has a pole, then a(c)=0, since otherwise a polynomial with unit leading
coefficient cannot have a polar root of P=p. Choose a local parameter
with t=c+s^alpha, alpha>=1. The coefficients of A(s,X)=P(c+s^alpha,X)
are in C[[s]], and its specialization degree is d=deg_x P(c,x).
Write nu(x)=-h<0. By (1),(3),

    ord_s omega >= ceil(D*h/(D-d))-2h+alpha-1
                >= alpha-1 >=0.                              (5)

This retains every branch, any ramification alpha, arbitrary multiple
escaping roots and arbitrarily large leading-coefficient multiplicity.
If d=D there are no escaping roots. No simple-root genericity is used.

At EVERY point with t=infinity write t=s^(-alpha), alpha>=1. The
monic equation P(t,x)/a(t)=p/a(t) has p-independent coefficients in
C((s)) and nu(1/a)=b*alpha. PF149's proved monic coefficient lemma gives

    nu(partial_p x)>=ceil(b*alpha/D),
    ord_s omega>=ceil(b*alpha/D)-alpha-1>=-1.                   (6)

The last inequality is exactly b>=D. In particular arbitrarily large
p-independent polar parts of x have not been discarded.

Thus omega has AT MOST SIMPLE poles on the entire complete generic
curve. It cannot be a nonzero rational exact differential: a pole of a
rational primitive gives a differential pole of order at least two;
a primitive without poles on a complete curve is constant.

A mate with Jac_(x,t)(P,Q)=c0!=0 would satisfy dQ=c0*omega along P=p,
a contradiction. In any finite extension the normalized pullback has
order E*ord(omega)+(E-1)>=-1, so the same argument applies there too.
This is a curve-differential obstruction, not polynomial trace landing
or a proposed stability assertion for an integral lattice.

## 3. Actual low-degree plane Keller consequence

For every NONINVERTIBLE plane Keller map and every polynomial source
frame (x,t), every target polynomial coordinate restricts to degree at
least TWO on every whole source-coordinate line t=c. This is the
retained [PF52 reached-line criterion](SOURCE_SYMMETRY_DESCENT_AND_TRACE_GATE.md), Section3:
degree one makes the map injective on the line; degree zero makes the
complementary target coordinate affine-linear there by the actual
constant Jacobian. The retained one-line theorem then forces inversion.

Consequently, for any such target coordinate P with 2<=deg_x P<=4,

    deg_t(leading_x P) <= deg_x P-1.                            (7)

Indeed degree at least two satisfies the half-degree requirement in
(3) for D=2,3,4. The conclusion uses NONINVERTIBILITY and the WHOLE
source lines. It is not a statement about an arbitrary polynomial P,
a general quasi-finite map, or punctured source lines.

## 4. Close the entire covered quadratic ruling-four family

By [PF157](CANONICAL_DECK_JET_FIELDS_AND_RATIONAL_POLE_ATLAS_COVERAGE.md), EVERY actual quadratic
factor L/E, with no minimality assumption, has marked source image

    V_R0: w^2=u*Sigma(u)^2+u*A(u)*z,
    u=t^2, w=t*(Sigma+A*x), z=2*Sigma*x+A*x^2,                   (8)

where A is monic, A(0)=0, gcd(A,Sigma)=1, and deg A>=1. The form comes
from the reduced fraction Sigma/A=R0+2uR0' for proper rational R0 with
a pole at zero. This includes every order-two subgroup of a larger deck
group. Original target functions P,Q are regular on THIS actual image.

Suppose an actual target coordinate has generic ruling degree FOUR.
Its unique regular expression is

    P=c(u)*z^2+b1(u)*z+b0(u)+w*(a1(u)*z+a0(u)), c!=0.           (9)

All lower coefficient polynomials and c are arbitrary. There is no
coprimality, root-multiplicity, Newton-support or companion-degree bound.
Pull back (9) by the original whole-plane map (8). Its x-degree is four,
and its leading x coefficient is EXACTLY

    a4(t)=c(t^2)*A(t^2)^2,
    deg_t a4=2*deg c+4*deg A>=4.                               (10)

The w terms have x-degree at most three, so none has been omitted from
the leading calculation. This contradicts (7). Therefore

    EVERY actual quadratic-image target coordinate has ruling degree>=5.
                                                                    (11)

Degrees below four were already excluded by PF157's r>=2*g bound
with g=2 (or directly by the central reached-line degree). Equation
(10) now excludes degree four WITHOUT normalizing its bad divisor.
In particular ALL of PF151/PF157's r=4, s=4j+2 common-leading-divisor
residual, including arbitrary multiplicities>=3, is RETIRED.

At any nontrivial full deck order g, combine (11) for g=2 with PF157
for g>=3. After the legitimate elementary target reductions retained
in PF151/PF157, the remaining necessary ruling degrees are

    max(5,2g)<=r<s, neither r nor s divides the other.            (12)

All original target automorphisms preserve B,F,S,Lambda and every prime.
For g=2, the next degree is r=5; both-even pairs now begin at r>=6.
The higher-degree common-leading-divisor problem remains open.

## 5. Sharp local boundaries and a genuine new surface test family

The half-degree threshold cannot be weakened merely because two is a
necessary whole-line degree. Consider local polynomials

    P=w^d*(1-t^a*w)^e, D=d+e.

Put j=gcd(e,a*d), alpha=e/j, h=a*alpha, k0=a*d/j and

    t=s^alpha, Z=1/w=s^h*(1+s^k0*Y),
    Y^e=p*(1+s^k0*Y)^D.                                      (13)

For any lambda^e=p, the last equation has a unique Hensel solution
Y=lambda+O(s^k0), since e*lambda^(e-1)!=0. The parameter is normalized:
gcd(alpha,k0)=1. The first parameter coefficient of Z is at index
r=h+k0=D*h/e, so (1) is SHARP on these branches.

For D=4,d=2, any a>=1, alpha=1 and partial_p w has order zero;
omega is holomorphic, even for arbitrarily high leading multiplicity.
For D=4,d=1,a=3, alpha=1 and omega has order -2. This illustrates
why the actual PF52 degree-two condition matters. For D=5,d=2,a=4,
alpha=3,h=12,r=20 and omega again has order -2. Thus the present
argument DOES NOT close ruling degree five from line degree two.
These are local polynomial controls, not Keller pairs or counterexamples.

A regular family genuinely beyond PF149's leading-multiplicity gate is
on S2, with M=u^2,N=u and the source marking Sigma=-1,A=u:

    P_h=((u-1)^h*z+w)^2+z, h>=3.                               (14)

Its even z-leading coefficient c=(u-1)^(2h) is coprime to M and has
multiplicity 2h>4, outside PF149's n=2 finite-root condition. In the
whole plane marking, degree along t=0 and t=+/-1 is TWO, and along
every other fixed-t line is FOUR. Thus it passes the named line-degree
gate. Formula (3), or its plane pullback (10), excludes any rational
mate. No assertion that these P_h have a companion is made.

## 6. Full normalization, dimension and remaining objective

Keep the actual original B=C[P,Q], R=C[x,t], L/F, FULL finite
normalization S and codifferent Lambda. Every original ramified and
unramified boundary prime remains in the ledger. This proof uses the
actual source volume and all projective generic-curve branches. Its
auxiliary p-derivative is not assumed to preserve S; no replacement of
S by a translated open, comparison normalization or field trace algebra
is made. No unmarked trace is promoted to integral marked traces.

The identifiable plane inputs are the one-dimensional geometric generic
fibre with its scalar coframe, and PF52's injectivity-on-a-whole-line
consequence. A three-dimensional Keller map has a surface generic fibre
for one coordinate, so the complete-curve exact-differential argument
does not apply. The retained known three-dimensional nonproper Keller
control still has constant Jacobian and a polar first velocity. Constant
Jacobian, generic finiteness or derivative generation alone never supplies
our half-specialization/curve hypotheses.

PF142's mixed minimal case and PF156's fixed-v2 quintic graph branch
stay closed. The r=4 quadratic residual is now also closed. The trivial-
deck case, r>=5 quadratic targets, higher cyclic targets, the full marked
quadratic regularity needed for derivative integrality, and normal-tail
increasing-power containment remain OPEN. Full HC4/JC2 and an actual
counterexample remain OPEN. No global nilpotent-tail certificate is claimed.

The proof is the universal first-coefficient argument (1)--(6) and the
actual pullback (10), rather than a bounded coefficient experiment.
The accompanying exact controls replay new identities and retain pinned
predecessors. The checker passes 98 exact arithmetic and symbolic controls;
no predecessor checker or new external classification is invoked. Independent whole-proof replay, formal verification, external
review and novelty priority remain unclaimed.
