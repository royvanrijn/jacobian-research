# Complete twelve-slot weight(1,-2,-4)/(1,-2,-5) closure (KCUB37)

Subsequent [KCUB40](EQUAL_NEGATIVE_FULL_DEGREE7.md) closes full111 through
degree seven by a changed curve/Euler-Wronskian proof. The current
degree<=7 diagonal-weight raw<=12 frontier is empty:only112 can collide,
with exactly13 terms. All larger frontier and full111-unresolved descriptions
below are historical stages;original timeout records are unchanged.

Subsequent [KCUB39](ZERO_WEIGHT_CLASSIFICATION.md) closes every one-zero-weight
three-variable class in all degrees. The current degree<=7 diagonal-weight
raw<=12 frontier contains only full111. All fourteen/twelve/ten-type
frontiers and unresolved zero-weight descriptions below are historical
stages, not current next steps; their solver records remain unchanged.

Subsequent [KCUB38](DISTINCT_NONZERO_WEIGHT_MINIMUM.md) closes312/412, leaving ten current
diagonal-weight improvement types. The twelve-type stage below is retained
with all original timeout evidence.

Over every characteristic-zero field, every three-variable polynomial Keller
map of ordinary degree at most seven, with source/output weights(1,-2,-4)
or(1,-2,-5), is a polynomial automorphism. There is no support, collision
or coefficient-nonvanishing assumption. These are the two complete twelve-slot
models retained by KCUB36. Both close here by a changed boundary proof,
not by increasing the time cap on their old native coefficient queries.

After diagonal output normalization to identity linear part, the exact
families are as follows; all displayed parameters belong to the base field.

For weights(1,-2,-4), either

    S=c+eta*b^2+theta*a^2*b^3,
    F=(a,b+lambda*a^2*S,S),                              (4-I)

or

    Y=b+lambda*a^2*c,
    F=(a,Y,c+kappa*Y^2).                                (4-II)

All parameters are arbitrary. For weights(1,-2,-5), either

    S=c+eta*a*b^3+theta*a^3*b^4,
    F=(a,b+lambda*a^3*S,S), lambda*theta=0,              (5-I)

or, putting Q=kappa*b^2-a*c,

    F=(a,b-lambda*a^2*Q,
         c-2*kappa*lambda*a*b*Q+kappa*lambda^2*a^3*Q^2). (5-II)

All lambda,kappa in (5-II) are allowed. The polynomial Q is invariant under
(5-II), and its inverse replaces lambda by -lambda. The other families
are compositions of two elementary shears; (5-I)'s product-zero locus is
covered by lambda=0 or theta=0. Five generic inverse charts are directly
verified. No new counterexample or new construction mechanism is claimed.

Combining this new closure with the [complete KCUB36 census](WEIGHT_RESONANCE_CENSUS7.md)
reduces the degree-seven diagonal-weight raw<=12 frontier from fourteen
weight types to twelve: (3,-1,-2),(4,-1,-2), full(1,-1,-1), and the nine
(0,p,-q), gcd(p,q)=1,1<=p<=q,p+q<=7. This import is used only for the
frontier corollary; the two complete automorphism proofs below import no
external classification. Counts refer to the displayed weight coordinates;
general source-frame sparsification remains outside this frontier.

## Complete twelve-slot source modules

The three weights are distinct. Nonzero constant Jacobian makes the
necessarily diagonal linear part invertible. Output scaling normalizes it
to identity and the determinant to one, preserving degree and occurrences.
Set u=a^2*b,w=a^r*c, r=4 or5. Write F=(a*f,a^-2*g,a^-r*h).
Complete source enumeration gives

    f=1+mu*u+P*u^2+nu*w,
    g=u+A2*u^2+A3*u^3+L*w+B*u*w,
    h=w+E_(r-2)*u^(r-2)+E_(r-1)*u^(r-1)
         +H1*u*w+H2*u^2*w+T*w^2.                       (1)

Each coefficient multiplies exactly one native source monomial, so the
basis has twelve nonlinear occurrences, not a partial support template.
There are fifteen native determinant coefficient equations in each model.
The independent checker enumerates every compatible source monomial and
compares the full original determinant to

    J=f*{g,h}+2*g*{f,h}+r*h*{g,f}=1,
    {P,Q}=P_u*Q_w-P_w*Q_u.                              (2)

The fraction-field calculation extends to a polynomial source identity;
no a=0 boundary is omitted.

## Nonzero transverse first coefficient: complete graph reduction

If nu!=0, a diagonal c conjugacy makes nu=-1. Write f=N(u)-w,
N=1+mu*u+P*u^2, and v=f. With G(u,v)=g(u,N(u)-v) and similarly H,
(2) becomes

    (2G+v*Gv)*Hu-Gu*(rH+v*Hv)=1.                        (3)

The full native shape gives

    G=G0(u)+v*G1(u),       deg G0<=3, deg G1<=1,
    H=H0(u)+v*H1(u)+T*v^2, deg H0<=4, deg H1<=2.

At v=0 the exact Wronskian is 2G0*H0'-r*G0'*H0=1. Neither polynomial
is zero. If m=deg G0,n=deg H0 and 2n!=rm, the leading degree is m+n-1,
forcing m+n=1. For r5 the bounds permit no positive balanced pair, so
G0,H0 are both affine. For r4 the only additional pairs are(1,2),(2,4).
The latter is impossible by subtracting a leading power: put
H0=k*G0^2+R with deg R<4. The Wronskian is unchanged. If R=0 it is zero;
otherwise its leading coefficient has degree deg R+1>=1 and scalar
2*(deg R-4), so cannot be constant. Thus

    r4: deg G0<=1, deg H0<=2;
    r5: deg G0<=1, deg H0<=1.                            (4)

This removes the balanced quadratic boundary without any root or generic
parameter assumption.

Write G0=A*u-beta and G1=alpha*u+beta, using G(0,1)=0 at the original
source point u=w=0. That point also has H=0 and Hv=-1. The v^3 and v^2
coefficients of (3) are

    -(r+2)*alpha*T=0,
    3*G1*H1'-(r+1)*alpha*H1-(r+2)*A*T=0.                (5)

If alpha=0, G=A*u+beta*(v-1) is affine. For A=0, (3) is
beta*(3v-2)*Hu=1, impossible. For A!=0 take q=G. Its vector field sends
q to2A*q and v to-A*v, so (3) is

    A*(2q*Hq-v*Hv-r*H)=1.

Its polynomial solutions have constant term -1/(rA); every other monomial
q^i*v^j satisfies 2i-j=r, hence i>0. At the source q=0,v=1 this contradicts
H=0. This excludes all affine G in any degree.

Consequently alpha!=0, T=0, and H1 is nonzero because H1(0)=-1.
Leading degrees in (5) now require 3*deg H1=r+1. For r4 this is impossible.
For r5 it gives H1=gamma*G1^2. The value at zero implies beta!=0 and
gamma=-beta^-2. A diagonal b conjugacy normalizes beta=1 without roots.
Now G0=A*u-1,G1=1+alpha*u,H1=-(1+alpha*u)^2. By (4) and H(0,1)=0,
H0=1+e*u.

The constant boundary equation, normalized original g_u=1, and original
h_u=0 respectively give

    -2e-5A=1, A+alpha+mu=1, e-2alpha-mu=0.

This complete linear system has constant determinant -2. Therefore

    e=-(1+5A)/2, alpha=-3*(1+A)/2, mu=(5+A)/2.

Let E0=2G0*H0'-5A*H0-1 and
E1=3G1*H0'+2G0*H1'-5alpha*H0-6A*H1. Their u coefficients are

    Q=3A*(5A+1)/2,
    S=-3*(A+1)*(7A-5)/2.

The exact four-monomial identity

    ((21A+181)/960)*Q+((15A+128)/960)*S=1

excludes the final graph branch. P never needed a nonzero assumption or
an elimination: these necessary graph/source conditions already contradict
J=1 for every P. Every zero coefficient and normalization denominator
boundary has been covered separately. Thus nu=0.

## Nonconstant univariate first quotient

Put f=phi(u)=1+mu*u+P*u^2 and suppose it is nonconstant. Write

    g=A(u)*w+B(u), h=C*w^2+D(u)*w+E(u),
    deg A<=1, deg B<=3, deg D<=2, C constant,
    B(0)=E(0)=0, B'(0)=D(0)=1.

A multiple root of phi makes (2) zero, impossible. Its w^2,w,constant
coefficients are

    (4-r)*phi'*A*C+2phi*A'*C,
    phi'*[(2-r)*A*D+4B*C]+phi*(A'*D+2B'*C-A*D'),
    phi'*(2B*D-r*A*E)+phi*(B'*D-A*E').                  (6)

If A=0, J=(2phi'*B+phi*B')*hw. The first factor has positive degree
 deg phi+deg B-1 and nonzero leading scalar2deg phi+deg B, so cannot
be a polynomial unit. Thus A is nonzero.

For r5, if C!=0, evaluate at any simple root of phi. The first coefficient
forces A=0 there, the second forces B=0 there, and the constant coefficient
becomes zero. This is impossible, so C=0.

For r4 and C!=0, the first coefficient forces A'=0. Thus A is a nonzero
constant. The target shear h->h-(C/A^2)*g^2 preserves (2), exactly because
r=2*2, and removes its w^2 term. Its new w coefficient Dtilde still has
Dtilde(0)=1. The remaining w equation gives (phi^2*Dtilde)'=0, impossible.
The shear may raise degree; no original degree bound on the changed h or
Dtilde is used. The A=0 and C=0 boundaries were handled separately.

It remains to take C=0 in the original native degree bounds. The w equation
implies by a rational derivative that

    phi^(r-2)*D=k*A

for a nonzero field constant k. Since r>=4 and deg phi>=1, its left
degree is at least two; its right degree is at most one. This contradiction
excludes every nonconstant phi. Therefore f=1 in both complete classes.

## Exhaustive linear-first normal forms

Now {g,h}=1. Its w^2,w,constant coefficients are

    2A'*C=0, A'*D+2B'*C-A*D'=0, B'*D-A*E'=1.           (7)

If C=0 and A=0, polynomial units give B=u,D=1. If C=0 and A!=0,
D=k*A and A*(kB'-E')=1 force A constant and D=1. In either case
B=u+A*E. Substitution into the native source budgets gives (4-I), and
(5-I) with lambda*theta=0; its two branches retain every field point.

If C!=0, (7) makes A constant and nonzero. Put t=C/A!=0. Integration
using the normalized values gives

    D=1+2tB, B+tB^2-AE=u.                               (8)

For r4, E starts at u^2 and has degree at most3. Leading degrees in (8)
give deg B<=1, so B=u,E=(t/A)*u^2. Put lambda=A,kappa=C/A^2 to obtain
(4-II). The displayed family extends over lambda=0; its degenerate cases
were already included in (4-I).

For r5, E starts at u^3 and has degree at most4, so deg B<=2. Write
B=u+b2*u^2. The missing u^2 coefficient gives b2=-t, and (8) gives
E=(-2t^2*u^3+t^3*u^4)/A. Put lambda=A,kappa=C/A^2. Expanding back in
the source variables yields exactly (5-II). The displayed source formulas are polynomial in lambda,kappa and need
no roots; A=0 was treated separately. All C!=0 maps have lambda!=0;
the extension to lambda=0 is harmless and gives identity.

Conversely all four displayed families have degree<=7, the stated weights,
determinant one and polynomial inverses. Family I recovers
b=Y-lambda*X^(r-2)*Z and then subtracts its third-coordinate shear using
that b. Family(4-II) recovers c=Z-kappa*Y^2 and then b=Y-lambda*X^2*c.
Family(5-II) preserves Q and has inverse lambda->-lambda. This proves the
full classification over the original field, including a=0.

## Exact evidence, preserved timeouts and remaining objective

Run from this directory:

    python verify_weighted124125.py

The [independent checker](verify_weighted124125.py) uses only stdlib rational
arithmetic and the established sparse polynomial helper. It reconstructs
all24 source atoms and30 native equations, both full graph transformations,
three universal univariate coefficients and the target-square identity,
two rational-power numerators, four integer-boundary leading cases, four
graph coefficients, fifteen affine vector-field eigenvalues, the constant
linear pivot and four-term unit, and all five generic inverse charts.
The quadratic-invariant flow is checked directly. Its
[receipt](weighted-124-125/validation.json) pins every model, producer,
identity, old timeout checkpoint and this proof.

The old KCUB36 native models/scripts/logs and all four five-second timeout
records remain unchanged. Their statuses describe those solver runs, not
the current mathematical state of these two classes. No unchanged native
ideal was rerun, and no time cap was raised. The changed proof uses complete
graph degree reductions and a one-variable rational unit. The source
[producer](weighted124125.py) and [artifact note](weighted-124-125/README.md)
retain all new inputs. No bounded primary Groebner answer is a premise.

The degree, field-unit, rational-derivative, normalization and complete-case
arguments are written, not formalized or externally independently reviewed.
The two weight closures use no external theorem. KCUB36 is an explicit
import only for the twelve-type frontier corollary; its original imports
retain their own boundaries.

The remaining mixed types(3,-1,-2),(4,-1,-2), full111 and nine opposite-sign
zero-weight cases remain unresolved. This is a diagonal-weight-coordinate
frontier, not arbitrary source-frame or unweighted sparsity. All positive
bounds and62 excluded/12 remaining unrestricted profiles are unchanged.
Global degree, support and special-form dimension minima remain open.
