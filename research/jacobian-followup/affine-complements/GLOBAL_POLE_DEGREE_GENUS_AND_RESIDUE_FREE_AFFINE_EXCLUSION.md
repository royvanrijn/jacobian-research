# A global pole-degree genus obstruction closes the residue-free affine control

PF2RULED160. 17 November 2026. Actual complex polynomial PLANE Keller maps.

[PF159](POLYNOMIAL_MATE_RESPONSES_AND_QUADRATIC_RULING_FIVE_ESCAPE_COVERAGE.md) made progress by using the actual polynomial mate and by giving
complete necessary ruling-five escape coverage. Its residue-free affine
control is now EXCLUDED, including rational mates in every finite field
extension. The mechanism is GLOBAL: a primitive's pole budget bounds its
map degree, a zero-trace prime-degree ruling makes the two maps birational,
and the resulting genus upper bound contradicts an exact branch certificate.

The same proof excludes a Zariski-open one-parameter family with an explicitly
defined finite exceptional polynomial. The full ruling-five target class,
generic derivative integrality and HC4/JC2 remain OPEN. No normal-tail
nilpotence or increasing-power containment is inferred.

## 1. The complete-curve pole budget and the genus test

Let C be a smooth complete geometrically integral curve over an algebraically
closed characteristic-zero field k. Let u:C->P1 have PRIME degree D>=2.
Let omega be a nonzero rational differential with

    Tr_(k(C)/k(u))(omega)=0.                                 (1)

If omega=dQ for Q in k(C), a simple pole of omega is impossible. Every
pole of order m>=2 is exactly a pole of Q of order m-1. Therefore

    b=sum_(ord_A omega<0)(-ord_A omega-1)=deg(Q:C->P1)>0.      (2)

The sum uses EVERY complete-curve point. If there are no poles or a simple
pole, exactness is already impossible.

Q cannot belong to k(u): then (1) would give D*dQ=0. Since D is prime,
the intermediate field k(u,Q) is consequently all of k(C). The map
(u,Q) is birational onto its integral image in P1 times P1. Its bidegree
is (b,D), so

    genus(C)<=(D-1)*(b-1).                                  (3)

For completeness the image has arithmetic genus (b-1)(D-1), by
0->O(-b,-D)->O->O_image->0 and chi O(a,c)=(a+1)(c+1).
Normalization subtracts nonnegative delta invariants. This proves (3)
without assuming smoothness of the image or injectivity of the affine pair.

More generally, if primeness is unavailable, gcd(deg u,deg Q)=1 suffices
to make k(u,Q)=k(C). The PRIME condition in the zero-trace version cannot
simply be dropped. On Y^2=t^11-1, genus5, take u=t^2 and Q=t. Then
deg u=4,deg Q=2,Tr(dQ)=0, but k(u,Q)=k(t) is a proper index-two
subfield. The naive genus upper bound3 would be false. Its twelve
hyperelliptic branch points (including infinity) give genus5 directly.
This is only a curve control, not a Keller example.

If omega has no rational primitive, it has none after ANY finite
characteristic-zero function-field extension: normalized FIELD trace of
a primitive would give one downstairs. This is rational exactness, not
integral trace landing on an original polynomial target.

## 2. Actual ruling coframes have zero trace

Let k be an algebraic closure of C(p). On a geometrically integral generic
curve P*(u,w)=p, with P* in C(u)[w] of degree D>=2, let

    Omega=du wedge dw/M(u),
    omega=-du/(M*P*_w)=-partial_p(w)*du/M.                   (4)

The monic equation P*-p has w^(D-1) coefficient independent of p. Thus
Tr(w) is independent of p and

    Tr(partial_p w)=partial_p Tr(w)=0, Tr(omega)=0.            (5)

Equivalently sum_(roots w_i)1/P*_w(u,w_i)=0 follows by expanding the
partial fractions of 1/(P*(u,W)-p) at W=infinity. This uses ACTUAL
multiplication/field trace in the degree-D generic-curve algebra. It is
not a trace computation in the original full finite normalization S.

For ruling degree FIVE, (3) therefore applies to EVERY rational primitive:

    genus(C)<=4*(b-1),                                     (6)

with no coprimality assumption on b. For an actual quadratic top factor
L/E and actual target field F=C(P,Q), b is [E:F]; the original generic
Keller degree is 2b. The ruling degree5 and b are distinct quantities.
No automatic B-integrality follows from this degree identity.

A useful geometric-integrality check is also elementary. Suppose P*-p
has degree D over C(p)(u), is irreducible there, and at a constant u=a
its coefficients are regular and its finite w specialization has degree d
with d simple generic roots.
The number c of geometric generic components divides BOTH D and d:
the constants Galois action is transitive on components, preserves their
u degrees and preserves finiteness of w over u=a. Thus gcd(D,d)=1 makes
c=1. The total degree D and the finite branch count d are being counted
with all their geometric conjugates, not by choosing one sheet.

Every actual ruling-five target from PF159 has a degree-two specialization
off M. Hence this criterion applies. It also applies directly to the
family below, with D=5,d=2 at u=1. Irreducibility over C(p)(u) follows
from the elementary generic polynomial fact f(W)-p is irreducible over
K(p) for f in K[W] nonconstant: a factor independent of p would divide
the coefficient -1 of p. Clearing fixed u denominators introduces no
additional function-field component.

## 3. A one-parameter residue-free affine-response family

On S2, retain PF157's actual whole-plane SOURCE marking

    w^2=u+u^2*z,
    u=t^2,w=t*(-1+u*x),z=-2*x+u*x^2.                         (7)

Put v=u-1, T=1-v^10*w, g0=1-2v and take tau in C*:

    P_tau=w^2*T^3+v^2*g0*w*T+tau*v^41*z^2+v^11*w+v.          (8)

Tau=1 is exactly PF159's residue-free affine control. The parameter
changes none of the leading or local coframe identities there. For
clarity the complete audit is retained as follows.

The generic ruling degree is5, with leading coefficient -v^30, a fifth
power up to constant. At u=1 the two finite roots have w^2=p. There
are exactly THREE further normalized ends, all unramified over v:

    C=0:      w=v^-10+v^-2*Y, initial(P_tau)=-Y;
    C=+/-i:   w=v^-10+C*v^-4-C*v^-3+v^-2*Y,
                         initial(P_tau)=2Y+C.                (9)

All initial slopes are units. Direct differentiation before substitution
gives

    P*_w/v^2=-g0+O(v^2)   for C=0,
    P*_w/v^2=2*g0+O(v^2)  for C=+/-i.                       (10)

The first two coefficients are independent of Y and tau. Therefore
each coframe omega=-du/(u^2*P*_w) has a DOUBLE pole, with ZERO residue,
because u^2*g0=(1+v)^2(1-2v)=1+O(v^2).

There are no other finite leading roots except u=0. At u=0 the canonical
z^2 coefficient is -tau, a unit. The complete M-root audit of PF159,
including leading a2=0 ON M, proves every end there holomorphic. Its
monic-parameter infinity calculation likewise proves every u-infinity
end holomorphic. Every affine generic-curve point is retained and regular.
Thus (9) is the COMPLETE polar divisor of omega, and

    b=3.                                                   (11)

The geometric generic curve is integral by Section2, since D=5 and its
finite specialization at u=1 has degree2. Any rational primitive would
make genus<=8 by (6). The next section contradicts this globally.

The named source-line tests also stay valid for every tau!=0: restriction
degrees are2 on t=0,+/-1, and5 on all other fixed-t lines. The source
t derivative is nonzero on the whole t=0 line, and P_u=1 at fixed w
along u=1. These checks do not assert global submersion or a target pair.
At tau=0 the t=0 restriction is constant, so PF52 already rules it out
as an actual noninvertible quadratic-image target coordinate.

## 4. Exact discriminant and the finite exceptional polynomial

Clear the generic equation by u^4. With u=v+1, define the explicit quintic

    F(v,w,p,tau)=-u^4*v^30*w^5+(3u^4*v^20+tau*v^41)*w^4
      -3u^4*v^10*w^3+[u^4*(1-v^12*(1-2v))-2tau*u*v^41]*w^2
      +u^4*[v^2*(1-2v)+v^11]*w
      +u^4*(v-p)+tau*u^2*v^41.                              (12)

Let Delta=Disc_w(F), using the polynomial discriminant, and define

    R(v,p,tau)=Delta/[v^96*(v+1)^10].                        (13)

This is an INTEGER polynomial, of v-degree215, with identities

    leading_v R=-256*tau^6,
    R(0,p,tau)=-16p,
    R(-1,p,tau)=256*(p+1)^2*tau^6.                          (14)

Here is a universal derivation of the division and degree; it is not
inferred from a modular sample. Use Delta=a5^8 product_(i<j)(w_i-w_j)^2.

At v=0, (9) gives three polar roots with common leading v^-10 and
next distinct coefficients0,+i,-i at order v^-4. There are two finite
roots +/-sqrt(p). The orders in the product are

    8*30 - 6*20 - 3*8 =96;

the leading constant is (4p)*(-4)=-16p. This is valid over the generic
p,tau field and includes all five roots.

At u=0, four roots satisfy

    w=epsilon*sqrt(u)+(epsilon/2)*z0*u^(3/2)+...,
    epsilon=+/-1, z0^2=-(p+1)/tau.

The fifth root is -tau*u^-4+.... The two within-sign squared differences
have product (p+1)^2*u^6/tau^2; the four cross-sign ones give256*u^4;
the four differences with the polar root give tau^8*u^-32; and a5^8
gives u^32. The result is order10 and constant256*(p+1)^2*tau^6.

At v=infinity one root is tau*v^7+.... The other four have

    w=epsilon*sqrt(u)+delta,
    delta^2~epsilon*v^(-11/2)/(4tau), epsilon=+/-1.

These follow by substituting in (8): the competing terms are
-v^30*w^5 and tau*v^37*(w^2-u)^2, and all other terms have smaller
pole order. The two within-sign squared differences have product
-v^-11/tau^2; the cross-sign ones give256*v^4; the four differences
with the large root give tau^8*v^56; and a5^8 gives v^272.
Thus deg_v Delta=321 with leading coefficient -256tau^6. This proves
all of (13),(14). The Puiseux conjugates are counted as roots for the
discriminant product, not mislabeled as unramified normalized points.
Divisibility over the generic parameter field implies exact division
in Z[p,tau,v], since both linear divisors in (13) are monic.

The EXPLICIT exceptional polynomial is the finite resultant

    B(tau)=Res_v(R(v,3,tau), partial_v R(v,3,tau)) in Z[tau].   (15)

Its coefficients need not be expanded to define it exactly. It is a
specified resultant of the quintic discriminant (12),(13), with fixed
integer division. Since Disc of a quintic is degree8 in its coefficients,
deg_tau R<=8 and

    deg_tau B<=8*(215+214)=3432.                             (16)

B is NONZERO. The retained exact certificate over F_101 at p=3,tau=1
computes R of degree215 and gives an explicit polynomial Bezout identity

    U*R+V*partial_v R=1,
    B(1)=12 modulo101.                                     (17)

The recorded leading value is47, R(0)=53 and R(-1)=56 modulo101,
agreeing with (14). These are exact finite-field polynomial operations,
not numerical root approximations or a probabilistic genus estimate.
The stored source recipe and Bezout vector are replay inputs. The checker
passes 70 exact arithmetic and symbolic controls.

For every COMPLEX tau!=0 with B(tau)!=0, R(v,3,tau) has215 distinct
roots, none at0 or-1 by (14). The same holds over C(p): the resultant
in p is nonzero because its specialization at p=3 is nonzero; the degree
and leading coefficient are fixed. At each of these simple discriminant
roots the quintic leading coefficient is nonzero. A simple discriminant
zero gives one simple quadratic ramification: factor off the other
simple roots, and the residual quadratic discriminant has a simple zero.
Thus EVERY such root contributes1 to geometric Riemann--Hurwitz.

For the integral degree-five cover C->P1_u, retaining every other
ramification point as well,

    2g(C)-2=-10+Ramification >= -10+215,
    g(C)>=104.                                             (18)

The other contributions can only increase this lower bound. They are
not omitted or assumed unramified. This contradicts genus<=8 from (11).

Consequently P_tau has NO rational constant-Jacobian mate, also in ANY
finite characteristic-zero field extension, for

    tau!=0, B(tau)!=0.                                     (19)

In particular tau=1 is excluded. Within this one-parameter family the
only remaining possible actual quadratic-image targets have nonzero tau
among the at most3432 roots of the explicitly defined B. This is a finite
necessary exceptional class, not a claim that any exceptional value is
realizable. Tau=0 already fails the retained whole-line criterion.

## 5. Exact source scope and the original full normalization

The certificate expands ONLY the modular specialization in (17). The
universal discriminant and B are defined by exact polynomial operations;
their full coefficient expansions are NOT claimed computed. Their
universal division, degree and boundary values have the written proof
above. A single exact modular nonvanishing proves B is not identically
zero; it is not a bounded experiment promoted to a universal identity.
No predecessor checker is replayed and no external classification is used.

The prime-degree field trace in (5) is taken in the GENERIC CURVE field.
It does not identify that field algebra with the original finite target
normalization, nor assert any integral marked trace landing. Keep the
actual B_original=C[P,Q], whole source R_original, original L/F,
FULL finite normalization S, codifferent Lambda and EVERY original
ramified and unramified boundary prime. Rational field trace is used
only to descend exactness from a finite extension. No stability of S or
Lambda under partial_p is assumed. Generic derivative integrality remains
UNSUPPLIED.

The identifiable two-dimensional inputs are the complete generic CURVE,
its scalar coframe from the actual plane volume, and the two ruling/mate
maps used in (3). A three-dimensional Keller example has a surface generic
fibre, so this genus argument does not apply to it. The retained known
nonproper3D control still guards against dimension-free integrality from
constant Jacobian, generic finiteness or derivative generation.

PF142's mixed minimal case, PF156's fixed-v2 quintic graph branch and
PF158's quadratic ruling-four branch remain closed. PF159's particular
residue-free affine control tau=1 is now closed; its wider necessary
affine-response ruling-five class remains OPEN. The genus criterion
(6) supplies a global test across that class, rather than another local
residue gate. The finite B-root family, other ruling-five and higher
quadratic/cyclic targets, trivial deck, full marked regularity and normal-
tail increasing powers remain OPEN. Full HC4/JC2 is not proved, and no
counterexample is constructed. Independent whole-proof replay, formal
verification and external review remain unclaimed.
