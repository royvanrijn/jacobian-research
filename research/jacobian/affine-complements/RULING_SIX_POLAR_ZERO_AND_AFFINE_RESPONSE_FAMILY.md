# A polar zero excludes a uniform ruling-six affine-response family

PF2RULED163. 17 November 2026. Complex quadratic source images.

On the retained smooth surface S2, an explicit infinite family of regular
ruling-six functions has no rational Jacobian mate, even after any finite
characteristic-zero extension of its function field. All three polar
coframe ends have affine first parameter response and zero residue. The
obstruction comes from a fourth escaping end: its coframe has a zero too
large for a primitive with the available pole degree.

This closes the displayed family, including every integer K>=62 and every
complex tau!=0. It does not close all ruling-six functions or prove HC4/JC2.
The later [simple-escape theorem](RULING_SIX_SIMPLE_ESCAPE_GRAPHS_AND_QUARTIC_FACE_REDUCTION.md)
also excludes its entire B face class; the pole/zero proof here remains
an independent obstruction for the displayed functions.

## 1. The regular family and the full geometric generic curve

Use the complete S2 source marking from
[PF157, Section 7](CANONICAL_DECK_JET_FIELDS_AND_RATIONAL_POLE_ATLAS_COVERAGE.md#7-complete-quadratic-image-presentation-beyond-mx0):

    S2: w^2=u+u^2*z,  Omega=du wedge dw/u^2,
    u=t^2, w=t*(-1+u*x), z=-2*x+u*x^2.                 (1)

Omega is a regular unit volume on the whole smooth affine surface; its
source pullback is -2 dx wedge dt. Put v=u-1, T=1-v^10*w and g=1-2v.
For every integer K>=62 and complex tau!=0 define the REGULAR function

    P_(K,tau)=w^2*T^3+v^2*g*w*T+v^41*z^2
                         +v^11*w+v+tau*v^K*z^3.         (2)

Its ruling degree is six, with canonical even z-leading coefficient
b3=tau*v^K. There is no bound on a hypothetical mate's degree.

Let p be transcendental, k an algebraic closure of C(p), and take the
smooth complete normalization of EVERY geometric generic component P=p.
On the whole affine generic curve define the nonzero coframe by

    Omega=dP wedge omega,  omega=-du/(u^2*P*_w),          (3)

where P* is obtained by substituting z=(w^2-u)/u^2. The affine generic
curve is smooth by generic smoothness; Omega makes omega regular and
nonzero there. No geometric-integrality assumption is needed below:
the final contradiction applies to the component carrying a specified
end, using the total pole count as an upper bound.

Writing P* as a degree-six polynomial in w over C(u), its top coefficients are

    p6=tau*v^K/u^6, p5=-v^30,
    p4=3v^20+v^41/u^4-3tau*v^K/u^5.                    (4)

Every geometric generic component has nonconstant u. Off u=0 the
degree-six leading coefficient can vanish only at u=1, where P=w^2
is nonconstant; on u=0 the restriction is cubic in z. No whole
constant-u fibre can therefore be a generic component.

The whole-plane source restriction has degree three on t=0, degree two
on t=+/-1, and degree six on every other fixed-t line. At t=0 the
original source t derivative of P is -2; at u=1 its u derivative at
fixed w is 1. Thus the family passes those named reached-line and
differential tests. These tests do not assert global submersion or
existence of a Keller target.

## 2. Exactly three double coframe poles over v=0

The first five local branches are controlled by the first five terms
of (2), namely the regular affine-response function of
[PF159, Section 6](POLYNOMIAL_MATE_RESPONSES_AND_QUADRATIC_RULING_FIVE_ESCAPE_COVERAGE.md#6-a-residue-free-affine-response-boundary-control).
We give the complete substitutions needed here.

There are two finite roots at v=0, with w^2=p. They are distinct and
simple over k. The three escaping branches have the following exact
initial substitutions, with Y a response variable:

    C=0:    w=v^-10+v^-2*Y,             initial P=-Y;
    C=+i,-i:
            w=v^-10+C*v^-4-C*v^-3+v^-2*Y,
                                         initial P=2Y+C. (5)

In each case the initial slope in Y is nonzero, so ordinary formal
implicit lifting gives an actual branch over k[[v]] with the displayed
finite p-independent polar truncation and first response. The parameter
v itself is normalized: these are three distinct unramified points,
not three conjugates of one ramified point.

On every substitution in (5), the added term tau*v^K*z^3 has exact
v order K-60>=2, with leading coefficient tau. Thus it changes neither
initial affine response. Its w derivative has order K-50>=12.
For the first five terms of (2), direct differentiation and substitution
give, as identities uniform in the response variable Y,

    C=0:    P*_w/v^2= -g+O(v^2);
    C=+i,-i: P*_w/v^2=2g+O(v^2).                         (6)

The coefficient assertions up to v^1 are independent of Y. Consequently
the later implicit correction to Y cannot change them. The added term
also cannot change them, by its derivative order. Equations (3),(6)
show that each of these three ends has an exact DOUBLE coframe pole.
Its residue is zero because

    u^2*g=(1+v)^2*(1-2v)=1+O(v^2).                       (7)

The escaping triple face is
I(Y)=Y^2*(1-Y)^3, of slope L=10, with the sixth coefficient strictly
above that face for every K>=62. It is the constant-K(Y) case of
PF162's off-M class B. The family therefore passes the affine-response
and local residue conditions at every pole. These conditions do not imply rational
exactness.

## 3. The sixth end has a coframe zero of order 4K-150

Put h=K-30>=32 and substitute w=v^-h*Y into the ENTIRE equation
P*=p. Multiplication by v^(5K-180) gives a polynomial in Y with
coefficients in C[p,tau,tau^-1][[v]] and reduction

    Y^5*(tau*Y-1)=0.                                    (8)

All lower terms have strictly positive v order; for example the
degree-four term of least possible order has scaled order K-40>0.
The nonzero root Y=tau^-1 is simple. Implicit lifting therefore gives
one more normalized branch, unramified over v, with

    w=tau^-1*v^-(K-30)*(1+O(v)).                         (9)

It differs from the three order-ten polar branches and from the two
finite branches. These six branches account for ALL roots of the
generic degree-six equation above v=0.

Using (4), the degree-six and degree-five terms in P*_w have the
same leading order and give respectively constants 6tau^-4 and
-5tau^-4. Every lower term has larger order. Hence on (9),

    P*_w=tau^-4*v^(150-4K)*(1+O(v)),
    ord_(9)(omega)=4K-150>=98.                           (10)

In particular this end is holomorphic and has a very large ZERO.
The cancellation coefficient in (10) is 6-5=1; it is not omitted.

## 4. Every other end is holomorphic

If u has a finite limit different from 0 and 1, the top coefficient
in (4) is a unit and w has no pole. Then z is finite as well, so the
point belongs to the whole affine curve and is already covered by (3).

At u=0, the canonical coefficient b3=tau*(-1)^K is nonzero. An end
with z polar would have w pole order at most half that of z, by
w^2=u+u^2*z. The term b3*z^3 would then strictly dominate every
lower z term and every w*z^j term, j<=2, a contradiction. Finite z
also makes w finite. Thus there is no omitted affine-curve point
over u=0. This retains all six local sheets and all central
ramification; no rational-patch point is substituted for the whole
surface.

At u=infinity, write u=s^-alpha on any normalized branch, alpha>=1.
Apply the monic-parameter lemma of
[PF149, Sections 2--4](MONIC_PARAMETER_POLE_BOUNDS_AND_EVEN_LEADING_CYCLIC_TARGETS.md#2-a-universal-monic-parameter-lemma):
for a monic degree-D polynomial F over C((s)) with F(w)=p*h(s),
the first nonconstant coefficient of w has order at least
ceil(ord_s(h)/D). Its derivative with respect to p has exactly that
first coefficient order. The proof compares the initial polynomial
after a finite p-independent truncation; its nonzero initial cannot
vanish at a coefficient transcendental over C. It retains every
ramification index and all multiple leading roots.

Here D=6 and h=u^6/(tau*v^K), so ord_s h=(K-6)*alpha.
Since ord_s(du/u^2)=alpha-1, the resulting uniform bound is

    ord_s omega >= ceil(K*alpha/6)-1 >=0.                 (11)

Equations (5)--(11) and affine regularity exhaust the complete curve,
componentwise. The ONLY poles of omega are the three double poles
in (5). In particular the total prospective primitive pole budget
over all components is three.

## 5. A local degree of at least 99 contradicts the entire pole budget

On any smooth complete curve in characteristic zero, if a nonzero
rational differential is dQ, a pole of Q of order m has differential
pole order m+1. Thus (5)--(11) would make every primitive Q have at
most three simple poles on the component containing (9). It follows
that the nonconstant map Q to P1 has degree at most three.

At the point (9), equation (10) is holomorphic. A rational primitive
there is finite, and characteristic-zero termwise differentiation
gives the EXACT local order

    ord_(9)(Q-Q((9)))=4K-149>=99.                         (12)

This local multiplicity cannot exceed the degree of the entire map Q.
The inequalities deg Q<=3 and deg Q>=99 contradict one another.
This proof retains every end; a high zero alone, without the complete
pole budget, would not give the contradiction.

Therefore (3) is not rationally exact on the component containing
(9). There is no rational Jacobian mate in C(u,w). A mate in any
finite characteristic-zero extension would give one downstairs by
normalized FIELD trace, which commutes with the relative differential.
Hence no such mate exists after any finite extension either.

An actual whole-plane Keller companion is polynomial in the source
and belongs to the finite quadratic extension C(x,t)/C(u,w). It
would provide the forbidden primitive. Thus every function (2) is
excluded as an actual target coordinate on the marked source.

## 6. Scope and verification

The proof is uniform in all integers K>=62 and all complex tau!=0.
No genericity restriction or finite exceptional tau polynomial is
needed. The coframe zeros and complete pole budget, rather than a
bounded coefficient search, establish the exclusion. Tau=0 drops
the ruling degree to five and is outside this statement.

The checker verifies the full substituted identities for K=62,63,70,
the unchanged affine responses and first two derivative coefficients,
the new simple root and derivative constants, the whole-source
restrictions, and the symbolic exponent bounds used for all K>=62.
The checker passes 96 exact controls. Those finite controls support
the written uniform proof. Replay from the repository root with

    .venv/bin/python research/jacobian/affine-complements/verify_ruling_six_polar_zero.py

The default writes a new receipt to a temporary directory; the retained
receipt records the executed run and source hashes.

This supplies an explicit excluded family inside the ruling-six
affine-response boundary, not complete ruling-six closure. The original
full normalization, all its primes and codifferent remain unchanged
in an actual-map application; the curve calculation does not assert
full-normalization derivative integrality or normal-tail termination.
No new external classification is imported. Independent whole-proof
review, formal verification and external human review remain unclaimed.
