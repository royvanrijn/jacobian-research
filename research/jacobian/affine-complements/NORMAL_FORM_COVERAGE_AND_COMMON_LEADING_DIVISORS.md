# Normal-form coverage, common leading divisors, and dimension-three controls

PF2RULED151. 17 November 2026.

Priority: prove COVERAGE of hypothetical minimal plane counterexamples,
rather than enlarge the list of individually excluded coordinate families.
Global coverage is OPEN. This note proves a conditional normal form inside
the marked quadratic models and identifies its exact remaining leading
divisors. It does not assume that these models cover every counterexample.

## 1. Three global coverage obligations

Retain F subset L for a hypothetical noninvertible complex plane Keller
map of MINIMAL GENERIC DEGREE, its target algebra B, full finite
normalization S and trace-dual Lambda.

The following are distinct obligations.

G1. Exclude the case with NO actual proper source-LND-equivariant
intermediate F subsetneq E subsetneq L, or supply another complete
normal form for it. Derivative generation alone is not integrality.

G2. If such an intermediate exists, PF142 makes it an ACTUAL cyclic
deck fixed field in minimal generic degree. Classify the actual marked
source image Spec(R intersect E), or retain it as an unclassified
cyclic-image residual. MX0 is an explicit family, not a classification
of all such images. Cyclic orders greater than two and unclassified
quadratic images cannot be discarded.

G3. For an actual quadratic image identified with S2 or MX0, prove
that the coupled leading-divisor residual below is impossible or can
be moved by justified source/target coordinate changes into an existing
exclusion. Individual constant-leading exclusions do not prove this.

PF142's mixed 2+3 minimal alternative is RETIRED. It is not a residual
in G2. No automatic source equivariance, proper intermediate, model
coverage or good leading coefficient is asserted.

The conditional normal form below is useful because it exposes an
obstruction to a generic "make the leading coefficient constant"
argument: after valid target reductions, an odd ruling degree forces
BOTH leading coefficients to have a DOUBLE zero at u=0.

## 2. Degrees and hypotheses of the conditional normal form

Assume the actual factor has the marked smooth quadratic model

    V={w^2=N(u)+M(u)*z}

of S2 or MX0. In both cases M has a double zero at u=0, N has a
simple zero there, and the original whole-plane chart has a reached
source-coordinate line t=0 on which

    u=w=0, z=nonconstant affine-linear function of x.

Every other M-root fibre is also reached in its entirety by fixed-t
source-coordinate lines, with constant w and affine-linear z.
For MX0, M=u^2*b0(u)^k, k>=2, b0 square-free, b0(0)!=0.

The actual target functions P,Q belong to O(V) and have

    dP wedge dQ = c*Omega, Omega=du wedge dw/M, c!=0.

Scale one target coordinate to make c=1. On M!=0 substitute
z=(w^2-N)/M. Define r=deg_w P*, s=deg_w Q*, where P*,Q* are
polynomials in w over C(u). These are GENERIC RULING DEGREES;
they are not source total degrees, potential degrees, graph degrees
or the generic degree [L:F].

Write the unique regular expressions as

    P=w*A_P(u,z)+B_P(u,z),
    Q=w*A_Q(u,z)+B_Q(u,z).

For a ruling degree r put epsilon_r=r mod2 and q_r=floor(r/2).
Let c_P(u) be the nonzero leading coefficient of B_P in z if
r is even, and of A_P if r is odd. Then the coefficient of w^r
in P* is

    l_P=c_P/M^q_r.                                  (1)

The analogous statements hold for Q.

PF52's reached-line criterion applies to EVERY polynomial target
coordinate, including every coordinate obtained in the elementary
target changes used below. In particular the restrictions of both
P,Q to t=0 have degree at least TWO for a noninvertible target.

## 3. Plane Jacobian identity forces common leading powers

The actual TWO-dimensional identity is

    P*_u Q*_w - P*_w Q*_u = 1/M.                    (2)

Since r,s are positive, the coefficient of w^(r+s-1) gives

    s*l_P'*l_Q-r*l_P*l_Q'=0.                       (3)

Thus l_P^s/l_Q^r is a nonzero constant. Put

    d=gcd(r,s), r=d*r0, s=d*s0, gcd(r0,s0)=1.

Rational-function divisor orders in C(u) then give

    l_P=alpha*H^r0, l_Q=beta*H^s0,                 (4)

for H in C(u)* and alpha,beta in C*. To verify this without
adjoining a root of u, (3) implies
s0*ord(l_P)=r0*ord(l_Q) at EVERY point of P1. Coprimality
makes these orders r0*h and s0*h for an integral divisor.
The corresponding rational H exists on P1; the remaining
functions have zero divisor and hence are constants. Constants
can be adjusted because C is algebraically closed. Equivalently
a Bezout combination of r0,s0 constructs H after constant scaling.

This is a plane-specific use of (2), not a dimension-free
consequence of constant Jacobian or generic finiteness.

## 4. Valid elementary target reduction terminates

If r divides s, (4) gives l_Q=lambda*l_P^(s/r), lambda in C*.
The polynomial TARGET AUTOMORPHISM

    (P,Q) -> (P,Q-lambda*P^(s/r))                   (5)

strictly lowers the ruling degree of Q. The symmetric operation
applies when s divides r; equality is included. It preserves
the target algebra B, target field F, original generic degree,
FULL S and Lambda, and the constant Jacobian identity.

Repeated application terminates because the sum of the two
nonnegative integer ruling degrees decreases. It cannot end
with a degree at most three: such a function restricts at
u=w=0 to a polynomial in z of degree at most one, contrary
to PF52. The same rules out degree zero.

Consequently every hypothetical compatible noninvertible target
admits an elementary target frame with

    r,s>=4, r<s, r does NOT divide s,
    r0,s0>=2 and gcd(r0,s0)=1.                     (6)

Minimal generic degree is preserved if it was assumed.

The approximate relation Q^r0-lambda*P^s0 also cancels leading
terms, but when r0,s0>=2 it is NOT thereby a target coordinate.
It must not be silently substituted for operation (5).
This is exactly where a further reduction theorem is needed.

## 5. Complete leading-divisor residual inside these models

There are two parity cases.

### 5.1 Both ruling degrees even

Then d is even. Put

    G=H*M^(d/2).

By (1),(4),

    c_P=alpha*G^r0, c_Q=beta*G^s0.                 (7)

Although initially rational, G is a POLYNOMIAL: either of its
positive powers is polynomial, so it has no finite pole.

If G is coprime to M and every root of G has multiplicity at
most d, every root of c_P has multiplicity at most d*r0=r.
PF149 then excludes P as an even-leading target coordinate.
The same applies to Q. Hence EVERY survivor must satisfy

    gcd(G,M)!=1
      OR some root of G has multiplicity >=d+1.   (8)

In particular G is nonconstant. This is a COMMON leading
divisor for the two coordinates, not two independently
chosen bad coefficients.

There is a sharper lowest-degree normal form. If r=4,
a zero of G at ANY M-root kills the z^2 coefficient of P.
On each corresponding reached whole source line, w is
constant and all remaining terms have z-degree at most one.
PF52 excludes that. Therefore r=4 forces

    s=4*j+2>=6, d=2,
    c_P=alpha*G^2, c_Q=beta*G^(s/2),
    gcd(G,M)=1, and some root of G has multiplicity >=3. (9)

These are necessary top forms, not constructed Keller pairs.

### 5.2 At least one ruling degree odd

Let h0=ord_(u=0) H and J=d+h0, an integer. Because ord_0 M=2,

    ord_0 c_P=r0*J-epsilon_r,
    ord_0 c_Q=s0*J-epsilon_s.                     (10)

Polynomiality and the odd component force J>=1.
After (6), an odd r0 or s0 is at least THREE; an even one
is at least TWO. Thus (10) gives

    u^2 divides BOTH c_P and c_Q.                 (11)

In particular neither leading coefficient is admissible
under the coprimality-to-M hypotheses of PF149/PF150.
This is a concrete reason those coordinate exclusions do
not automatically provide coverage.

If one component had ruling degree four, its even
z^2 coefficient would vanish at u=0 by (11), and its
restriction to t=0 would have degree<=1. So this mixed/odd
case has NO component of ruling degree four.

Equations (4),(6),(10),(11), together with all retained
exclusions in every target frame, describe the residual.
No assertion is made that a polynomial pair satisfying
only these leading data is Keller or exists.

## 6. The central first-order data are coupled too

On V put K=N'+M'*z. Since N'(0)!=0,

    Omega=dw wedge dz/K

near the reached central curve u=w=0. Both partial
derivatives of u in the local coordinates (w,z) vanish
on that curve. Restricting (2) intrinsically therefore gives

    A_P(0,z)*B_Q'(0,z)-B_P'(0,z)*A_Q(0,z)
       =1/N'(0).                                 (12)

Also deg B_P(0,z),deg B_Q(0,z)>=2 by PF52.
Thus the residual is subject to a genuine polynomial
Bezout identity and both whole-line degree constraints.
This uses the ACTUAL constant-Jacobian coframe.

A unimodular polynomial row or this first-order identity
does not supply global inverse completeness or full
trace-dual multiplication. The original S is not replaced
by the leading coefficient ring or by a primitive order.

## 7. A mandatory three-dimensional control for integrality

The known Alpoge/Gallagher map, already retained in
HC4MYKO2 and HC4MYKW1, is

    U0=1+x*y,
    A=U0^3*z+y^2*U0*(4+3*x*y),
    B=y+3*x*U0^2*z+3*x*y^2*(4+3*x*y),
    C=2*x-3*x^2*y-x^3*z.

Its Jacobian is -2. Three retained rational points give
a collision, and the explicit source arc

    x=1/t, y=0, z=2*t^2

has target (2*t^2,6*t,0), regular at t=0. It proves
nonproperness directly. Along this arc the differential
matrix is

    [ 0       6*t        1       ]
    [ 6*t^2  13          3/t     ]
    [-4      -3/t^2     -1/t^3  ].                (13)

For a generic constant direction v and generic target
linear function f=lambda*A+mu*B+nu*C, the first velocity
delta_v f has pole order THREE if nu*v3!=0.

This velocity is NOT integral over the target polynomial
algebra: evaluating any hypothetical monic relation on
the arc gives coefficients in C[[t]], whereas its
highest velocity power has a unique strictly lowest
valuation. This is impossible. This reasoning uses
integrality in the FULL normalization, not an arbitrary
primitive order.

Even source derivative generation does not save a
dimension-free lemma. All mixed ordinary partials of
C generate the WHOLE source ring, since

    x=-C_xy/6, z=-C_xxx/6,
    y=-C_xx/6-x*z.                               (14)

Thus nonzero constant Jacobian, generic finiteness,
finite mixed-jet generation and failure of generic
first-velocity integrality occur together. We make
NO claim here about generation by one fixed-direction
jet algebra; that is a different hypothesis.

The plane-specific step (3) visibly fails in this
control if misapplied in dimension three. All three
components have z-degree one, but their leading
coefficients U0^3,3*x*U0^2,-x^3 are not constant
multiples of one another. For example

    partial_y(3*x/U0)=-3*x^2/U0^2 !=0.

A three-by-three determinant with TWO base variables
does not give the scalar one-variable Wronskian (3).
The known weighted two-dimensional quotient also has
Jacobian -2*(2-3*u-v)^2, not a nonzero constant.

Any proposed automatic integrality implication must
state the extra plane hypothesis and identify the step
that uses it. The present leading normal form does
NOT itself prove integrality.

Primary control attribution:
[announced-map audit](https://github.com/shadybrook/jacobian-counterexample-research/blob/main/paper/main.md),
with its links to Alpoge's announcement and Gallagher's
weighted derivation. This is a rechecked CONTROL, not
a new three-dimensional counterexample claim.

## 8. The November preprint and the graph-degree convention

The named source is David Niedbala Giraudin,
[Hexagon 2611.00161v1](https://hexagonmath.org/2611.00161v1),
6 November 2026. It claims HC4 for potentials of TOTAL
degree at most five and explicitly states it is
unrefereed. Five ancillary programs report 247 checks;
the separate ternary-top dependency contributes 29.

The preprint's total-potential-degree hypothesis must be distinguished
from the graph-degree convention used here. No formal verification or
external review of that imported statement is asserted.

The local graph degree is DIFFERENT. In the retained
Meng--Yang family, U=1+x*y and

    A=U^3*p+3*x*U^2*q-x^3*r,
    B=y^2*U*(4+3*x*y)*p
       +(y+3*x*y^2*(4+3*x*y))*q+(2*x-3*x^2*y)*r,
    Phi=L*A^2+C*A+N*B, L*N!=0.

After r=R(x,y,p,q) with EXACT degree five and top R5,

    deg(A|graph)=8,
    top_8(A|graph)=-x^3*R5,
    top_16(Phi|graph)=L*x^6*R5^2 !=0.

Hence the graph potential has EXACT degree SIXTEEN,
not five. The external theorem does not remove those
unmarked graph branches by a degree comparison.
The old standard graph COLLISION route is already
closed in all degrees by HC4MYGALL1; unmarked
constant-determinant graph termination is separate.

Under cotangent doubling Phi=y^T F(x), potential
degree is 1+deg F. So the degree-five potential
theorem supplies plane maps only through degree
four, not a new arbitrary-degree JC2 frontier.

The distinctions are:

| Quantity | Degree in this discussion |
| --- | --- |
| Preprint's potential | total degree <=5 |
| Exact quintic graph R | degree 5 |
| Its Meng--Yang pulled-back potential | total degree 16 |
| Plane F in cotangent doubling | deg Phi-1 |
| r,s in (6) | degrees in w over C(u), not total degrees |
| Minimality in PF142 | generic field degree [L:F] |

## 9. Current coverage boundary

The written conditional normal form (3)--(12) is
proved inside the specified marked quadratic models.
The GLOBAL coverage theorem G1--G3 is NOT proved.

The immediate task is to eliminate or normalize the
common bad divisor in (8) or (11), using the actual
plane Jacobian and full original S. In the all-even
case, the lower-degree element
Q^r0-lambda*P^s0 belongs to B, but it is not known
to be a target coordinate. Treating it as one would
hide the missing normal-form theorem.

The no-proper-factor and unclassified-cyclic-image
alternatives remain explicit global residuals.
No ramified boundary prime is discarded; no S or
Lambda is replaced. PF142's mixed minimal retirement
and all completed proofs remain preserved. Full
HC4/JC2 and an actual counterexample remain OPEN.
