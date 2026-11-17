# Odd leading cyclic coordinates and the reciprocal boundary coframe

PF2RULED150. 17 November 2026. ALL degrees over C.

The odd highest-degree alternative left by
[PF149](MONIC_PARAMETER_POLE_BOUNDS_AND_EVEN_LEADING_CYCLIC_TARGETS.md)
can also be excluded on a substantial unbounded family.

On S2, EVERY

    P=w*A(u,z)+B(u,z),
    deg_z A=n>=1, deg_z B<=n,                      (1)

has NO rational Jacobian mate in any finite source-field extension
provided the leading coefficient a(u) of A is coprime to M=u^2
and every zero of a has multiplicity <=2*n+1. ALL lower coefficients
are arbitrary polynomials in u, of unbounded degree. In particular
every nonzero CONSTANT leading a is allowed.

On quadratic MX0 the same conclusion holds if, at each root r of
the model polynomial b0 where

    b(r)^2=N(r)*a(r)^2,                            (2)

one has n>=k-1. Here b(u) is the coefficient of z^n in B, possibly
zero, and M=u^2*b0(u)^k. No restriction is needed at roots where
(2) fails. All leading cancellations, higher contacts and ramified
branches satisfying this condition are retained.

For ACTUAL noninvertible plane Keller targets, n=1 is already
excluded by the retained source-line criterion. The whole n=2
family in (1), with the stated leading a, is also excluded on
EVERY quadratic MX0 model: the remaining leading cancellations
make the restriction to a reached source line affine-linear or
constant. This last extension is a REGULAR source-line argument;
we do NOT assert a universal rational-mate exclusion outside
the displayed local pole bound.

Together with PF149, this excludes all constant-leading coordinates
of either parity on S2 and on quadratic MX0 models with k<=4.
It does not show that every target coordinate has constant or
admissible leading coefficient. Full HC4/JC2 and first-derivative
integrality remain OPEN.

## 1. Setup and the whole generic curve

Keep V={w^2=N(u)+M(u)*z}, Omega=du wedge dw/M, with the complete
smooth source image of S2 or a quadratic MX0 model. At u=0, M has
order TWO and N has a SIMPLE zero. At every nonzero root r of M,
M has order k>=2 and N(r)!=0. These statements follow directly
from the actual MX0 formula

    M=u^2*b0(u)^k,
    N=u*(-b0+2*(1-k)*u*b0')^2,
    k>=2, b0 square-free, b0(0)!=0.

S2 has M=u^2,N=u, with no nonzero roots of M.

Write a=a_n, b=b_n, and D=2*n+1. Assume

    gcd(a,M)=1,
    ord_r a<=D for EVERY zero r of a.              (3)

Work on EVERY component of the geometric generic curve P=p over
an algebraic closure of C(p). As in PF149 it is smooth, the
relative form omega defined by Omega=dP wedge omega is regular
and NONZERO on its WHOLE affine locus, and u is nonconstant.
No component of the generic curve can be contained in a fixed
fibre u=r over C. All affine M=0 patches are retained.

On M*a!=0, put z=(w^2-N)/M. Then

    Q_u(T)=M^n*P*(u,T)/a(u)

is MONIC of degree D=2*n+1, and

    Q_u(w)=p*h(u), h=M^n/a,
    omega=-(partial_p w)*du/M.                    (4)

Here partial_p fixes u; the coframe is the actual one, with all
mixed coefficients retained. Curve degree D is NOT the original
target extension degree or integrality over its target B.

We reuse PF149's proved monic parameter lemma: for any monic Q
over C((s)), Q(W)=p*h with p-independent h gives

    nu_s(partial_p W)>=ceil(nu_s(h)/deg Q).        (5)

It follows from the first Laurent coefficient of W outside C;
no bounded degree, root simplicity or nondegeneracy is assumed.

## 2. Infinity over u=infinity and zeros of a

At u=s^-alpha, alpha>=1, put m=deg M and d=deg a.
Equation (5) and the actual coframe (4) give

    ord_s omega
      >= ceil((((n+1)*m-(2*n+1)+d)/D)*alpha)-1
      >=0.                                       (6)

Indeed m>=2 and the numerator is at least 1+d>0. This covers
EVERY ramification index and arbitrary p-independent polar
parts of w.

At a finite zero r of a, M(r)!=0 by (3). With
u=r+s^alpha and e_r=ord_r a, the same calculation gives

    ord_s omega>=ceil((1-e_r/D)*alpha)-1>=-1.      (7)

Every equality e_r=D is retained. It can give a logarithmic pole.

If M(r)*a(r)!=0, no finite-u infinity point exists: a pole of z
would make w polar of half its order, so w*a*z^n has strictly
higher pole than every other term of (1). It cannot cancel
in P=p. Thus only the M=0 branches remain.

## 3. Exact reciprocal coframe at every M=0 end

At a finite-u infinity point over u=r with M(r)=0, one has
a(r)!=0. A pole of w is again impossible by odd leading
dominance. Therefore w is finite, z is polar, and

    w(0)=-b(r)/a(r).

Put U=u-r and y=1/z. On EACH normalized branch
alpha=nu(U)>=1 and beta=nu(y)>=1. Define the p-independent
formal polynomials

    Atilde=y^n*A(u,1/y), Btilde=y^n*B(u,1/y).

Atilde is a UNIT at (U,y)=(0,0). On the generic curve,

    w=W(U,y,p)=(p*y^n-Btilde)/Atilde,
    H=M-y*(W^2-N)=0,
    Hbar=Atilde^2*H
       =M*Atilde^2-y*(Btilde^2-N*Atilde^2)
             +2*p*Btilde*y^(n+1)-p^2*y^(2*n+1).  (8)

The EXACT relative differential is

    omega=-Atilde*y^(n-1)*dy/Hbar_U
          =-y^(n-1)*dy/(Atilde*H_U).              (9)

For completeness, differentiate H=0 at fixed y and at fixed p:

    H_U*U_p=2*y*W*W_p,
    H_U*U_y=(W^2-N)+2*y*W*W_y,
    W_p=y^n/Atilde.

The determinant of (u,w) with respect to (p,y), divided by
M=y*(W^2-N), is -W_p/(y*H_U), which is (9).
This calculation does not divide by w at the boundary, and
therefore also applies when its limiting value is ZERO.

These are completed GENERIC-CURVE coordinates, not an affine
intermediate plane or a replacement finite normalization.

## 4. Nonzero residual discriminant

Let e=ord_r M and put Delta=b^2-N*a^2. If Delta(r)!=0,
W^2-N has nonzero constant term. Equation (8) forces

    beta=e*alpha.

In H_U=M'-y*(2*W*W_U-N'), the first term has order
(e-1)*alpha and every other term has order at least beta.
Consequently

    nu(H_U)=(e-1)*alpha,
    ord omega=n*beta-1-(e-1)*alpha
             =((n-1)*e+1)*alpha-1>=0.            (10)

No branch is removed by the nonzero discriminant: these
branches generally EXIST and are being proved holomorphic.

## 5. The zero limiting w at the simple N root

Suppose N(r)=0 and Delta(r)=0. On these actual models r=0,
e=2, N has a simple zero, b(r)=0, and W lies in the ideal (U,y).
The following calculation works for any e>=2 with N simple.

Write M=c*U^e+..., N=d*U+..., c*d!=0. Since W in (U,y),
W_U is regular and nu(W)>=min(alpha,beta). Equation H=0 is

    M+y*N-y*W^2=0.                               (11)

If alpha<beta, its lowest terms require
beta=(e-1)*alpha; then M'+y*N' has nonzero leading
coefficient (e-1)*c*U^(e-1), of order beta.
The term y*W*W_U has larger order.

If alpha>beta, y*N' alone leads H_U, again of order beta.

If alpha=beta, (11) is impossible for e>2. For e=2 the
two leading terms force c*U+d*y=0, and
2*c*U+d*y=c*U has nonzero leading coefficient. Its order
is beta. In ALL cases,

    nu(H_U)=beta,
    ord omega=(n-1)*beta-1>=-1.                   (12)

For n>=2 this is nonnegative; n=1 retains possible simple
poles. This proves the full central S2 boundary case, including
all cancellations and all ramified branches.

## 6. Nonzero N and arbitrarily high leading contact

It remains to treat N(r)!=0, Delta(r)=0. Now b(r)!=0.
Write e=ord_r M. In (8), the p-dependent terms have order
at least (n+1)*beta, and

    nu(Hbar_p)=(n+1)*beta.                        (13)

We prove, on EVERY normalized branch,

    nu(U_p)>=ceil((n+1)*beta/e).                  (14)

Choose y=s^beta and use coefficientwise partial_p fixing s.
Let q be the exponent of the first coefficient of U outside C.
If all coefficients were in C, (8),(13) would contradict
transcendence of p. Earlier terms give a finite p-independent
truncation U0(s), and write U=U0+s^q*(lambda+O(s)),
lambda not in C.

The p-independent part H0 of Hbar satisfies

    H0(U,0)=M(U)*a(U)^2,

of U-order e with nonzero leading coefficient. In
H0(U0+s^q*Y,s^beta), the coefficient of Y^e has exact
order e*q and nonzero constant leading coefficient.
Higher U powers and terms divisible by y contribute only
higher order to that coefficient. Hence the minimum
coefficient order is at most e*q.

If e*q<(n+1)*beta, all p-dependent terms of Hbar lie
strictly above this minimum. Its nonzero initial polynomial
has coefficients in C and would annihilate lambda, forcing
lambda in C. Contradiction. Thus e*q>=(n+1)*beta.
The constants of partial_p are C, so nu(U_p)=q. This proves
(14) without a root multiplicity or contact bound.

Differentiate (8) at fixed y. By (13),(14),

    nu(Hbar_U)=(n+1)*beta-nu(U_p),
    ord omega=nu(U_p)-beta-1
      >=ceil((n+1)*beta/e)-beta-1>=-1             (15)

whenever n+1>=e. If n+1>e, the bound is nonnegative.
The equality n+1=e is NOT discarded: it can give actual
nonzero logarithmic residues, as the control below shows.

For MX0 the nonzero roots have e=k, yielding the stated
condition n>=k-1 only at those roots with Delta(r)=0.
Every higher contact and normalized branch under that
condition is covered by (14),(15), not assumed generic.

## 7. Nonexactness, finite extensions and the source-line corollary

Equations (6),(7),(10),(12),(15) cover the WHOLE projective
normalization, componentwise. The nonzero omega is regular
on the affine curve and has at most SIMPLE poles elsewhere.
PF149's elementary complete-curve argument therefore rules
out rational exactness: a pole of a primitive creates a
differential pole of order at least two; a pole-free primitive
is constant.

Every finite ramified pullback preserves this property:
E*ord+(E-1) is nonnegative for ord>=0 and equals -1 for
ord=-1. No index or branch is omitted. Normalized FIELD
trace is an alternative rational-exactness argument only,
not polynomial trace landing.

Now use the retained
[PF52](SOURCE_SYMMETRY_DESCENT_AND_TRACE_GATE.md) criterion
for an ACTUAL noninvertible plane Keller target: every target
coordinate restricts to degree at least TWO on every
source-coordinate line. A constant restriction also fails:
the constant Jacobian forces its companion to be affine
linear on that line, hence injectivity on the line and the
retained inversion theorem.

For n=1, at the reached source line t=0, w=0 and z is a
nonconstant affine-linear function of x. P restricts to
degree at most one. This is already PF52, not a new regular
linear-z frontier.

For n=2 on ANY MX0 model, the only roots not already treated
by the rational theorem are nonzero roots of b0 with
Delta(r)=0 and k>3. There w0=-b(r)/a(r) is one of the two
reached constant w values. The corresponding whole
fixed-t line has z affine-linear. Its z^2 coefficient
a(r)*w0+b(r) vanishes, so the restriction has degree<=1.
PF52 excludes this case too. We do not turn this regular
argument into a no-rational-mate theorem.

Thus all admissible-leading n=2 coordinates are excluded as
actual noninvertible target coordinates, in every k.

Together with PF149 and the n=1/degree-one line cases:
S2 and MX0 with k=2,3,4 exclude EVERY constant-leading
coordinate of both parities in the canonical expression
B+w*A. Higher k only leaves low odd orders n=3,...,k-2
outside this combined gate, and only when their leading
data cancel at a nonzero M root. Leading divisors outside
(3) are a separate open alternative in all models.
No reduction of arbitrary coordinates to these hypotheses
has been supplied.

## 8. Exact controls and a real limitation of the pole method

A genuinely mixed central S2 example is

    P=w*z^3+z^2+u*w*z+u^5+w.

It has degree two on t=0 and degree ten on the other
marked source axis h=0, under the whole S2 plane chart.
At the u=0 infinity point, (8) has branches with leading
u=y^2 and u=-y-y^2. Formula (9) has leading
-y*dy and +y*dy, respectively. Both branches are retained;
no pole at the zero limiting w was silently removed.

For an MX0 high-contact equality control, choose k=4,b0=u-1,
n=3, a=1 and

    b(u)=-u^3/2+17*u^2/4+3*u-3/4,
    P=(w+b(u))*z^3.

Here N=u*(1-7*u)^2, ord_(u=1)(b^2-N)=4, and
b(1)=6. The leading equation at u=1 is

    U^4+12*p*y^4=0.

There are FOUR distinct normalized local branches U=lambda*y+...,
lambda^4=-12*p; each is unramified over y. Equation (9) has
residue -1/(4*lambda^3), NONZERO on every branch. This is
the equality n=k-1 in (15), not an omitted degeneracy.

The same model with

    P=(w+b(u))*z^3+z^2+u^17+w

passes the marked source-line degree gate even on the
leading-cancelled line: its degree there is TWO, at the
other u=1 line THREE, and on ordinary fixed-t lines SEVEN.
Its nonzero constant leading a and all higher contacts are
still excluded by the theorem, without a companion bound.

The low-degree restriction in the pole proof is real.
For k=5,b0=u-1,n=3 and

    b(u)=u^4/4-13*u^3/8+55*u^2/8+25*u/8-5/8,
    P=(w+b(u))*z^3,

one has ord_(u=1)(b^2-N)=5 with N=u*(1-9*u)^2.
The leading equation is U^5+16*p*y^4=0: ONE normalized
branch with nu(U)=4,nu(y)=5 and five conjugate expressions
over y. Its differential has order 3*5-1-4*4=-2.
Therefore the logarithmic conclusion cannot simply be
asserted when n+1<k. This is NOT a counterexample: no
rational primitive or plane target is constructed, and
this particular control already fails PF52 on the
cancelled source line. It tests the exact boundary of
the general rational pole method, not a reopened target.

The odd degree-one case also stays separate: on S2,
P=w+g(u) has rational mate 1/u and is excluded from the
theorem n>=1. Source-line arguments, not a false rational
nonexactness assertion, handle it for actual targets.

## 9. Original finite normalization and remaining goal

Retain the original B,R,F,L, FULL finite normalization S
and trace-dual Lambda from
[PF143](MARKED_QUADRATIC_VELOCITIES_AND_CYCLIC_TRACE_POLES.md),
including EVERY ramified and unramified boundary prime.
The reciprocal curve equation, U_p, and monic polynomial
over C(u) are not S, a finite target order, or a claim of
normalization stability. Field trace was used only for
rational exactness. No Tr(R) subset B, marked quadratic
regularity or first-derivative integrality is inferred.

[PF142](PARTIAL_KUMMER_PLANE_CHARTS_AND_CYCLIC_EQUIVARIANT_FACTORS.md)'s mixed
2+3 minimal alternative remains retired. The full objective
still requires actual integrality, exclusion of every
remaining cyclic target, or an actual counterexample.
It is OPEN.

The universal proof is (4)--(15), with the complete M=0
audit and the written first-coefficient argument. Exact
supporting controls are not a bounded target census.
No new external theorem import or predecessor checker
replay is used. Independent whole-proof replay, formal
verification, external review and novelty priority are
unclaimed.
