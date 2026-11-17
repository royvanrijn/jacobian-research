# Every global degree-four coordinate conductor is excluded

Let K=(P,Q) be an actual complex polynomial plane Keller map, in a
source frame where P is monic of positive degree in x. Let

    (C[P,Q,y]:C[x,y])=kappa*C[x,y]

be its exact coordinate-order conductor. Then

    deg_x kappa != 4.                                    (1)

This includes special source lines on which the degree drops from
four to two. Together with the
[unrestricted nodal-action exclusion](UNRESTRICTED_NODAL_ACTION_AND_DEGREE_TWO_CONDUCTOR_EXCLUSION.md)
and the positive even-degree theorem, every noninvertible Keller
map in such a frame must satisfy

    deg_x kappa>=6.                                      (2)

This is a necessary constraint on actual Keller sources.
Higher conductor degrees, target derivative integrality, HC4 and
JC2 remain open.

The new whole-plane input is elementary but essential: since P is
monic in x, a conductor point with finite x-coordinate cannot have
an escaping mate over a finite y-value while retaining the same P.
This makes the generic two-node midpoint polynomial through every
degree drop. Two moving action periods then give incompatible
integer orders at a zero of the quartic leading coefficient.

## 1. Fibre degrees, singularities and the generic two-node model

Use the exact specialization of the
[coordinate-order theorem](GENERIC_COORDINATE_ORDER_CONDUCTOR_AND_CLOSED_COLLISIONS.md).
On every source line y=c, the normalization of
A_c=C[P_c,Q_c] is C[x], its exact conductor is kappa_c*C[x], and

    deg_x kappa_c=2*dim_C(C[x]/A_c).

For a noninvertible Keller map this degree is positive and even
on EVERY source line. Assume for contradiction deg_x kappa=4.
Every fibre therefore has degree four or TWO.

The [constant-degree-four proof](CONSTANT_DEGREE_FOUR_CONDUCTORS_AND_PERIOD_ODE.md)
excludes degree four on all fibres. Consequently the leading
coefficient a_4(y) is NONCONSTANT and has a complex zero. At every
zero of a_4 the fibre degree is two, and its order is the exact
one-node order. In particular its two conductor roots are distinct.

Every fibre normalization is an immersion, since P_x,Q_x cannot
both vanish. The local conductor-order argument in the linked proof
therefore classifies a degree-four fibre as two ordinary nodes or
one ordinary tacnode. No cusp or ramified normalization branch occurs.

The generic fibre must be two-nodal. If it were persistently
tacnodal, every generic conductor root would be repeated. But a
simple finite conductor root on any degree-two fibre lifts by
Hensel to a root with unit x derivative, and hence stays generically
simple. This is a contradiction.

The audited genus-two classification then makes the two generic
node pairs centrally symmetric. Their common midpoint is the
rational root centroid

    h(y)=-a_3(y)/(4*a_4(y)),                              (3)

where a_3 is the coefficient of x^3 in kappa. The reflection
x -> 2*h(y)-x sends each generic conductor point to its actual
same-(P,Q,y) mate.

The retained primary classification is Neumann, Invent. Math.98(1989),
Section8, page482, together with corrected Neumann-Rudolph Lemma7.1.
Its exact hypotheses, source PDFs and page locators are pinned in
[the source audit](sources/neumann-two-node-genus-two-source-audit.json).
The genus-two applicability and the node/tacnode models are proved
in the constant-degree-four note. No uniform family of target
automorphisms is assumed.

## 2. Monicity makes the centroid polynomial through every drop

Fix a zero c of a_4 and put v=y-c. Choose one simple root x_0 of
kappa(x,c). The implicit recursion gives

    x(v) in C[[v]], x(0)=x_0, kappa(x(v),c+v)=0.

Its generic mate is

    z(v)=2*h(c+v)-x(v) in C((v)).

The mate identities hold in the function field. For example,
kappa divides P(2h-x,y)-P(x,y) in C(y)[x], because the latter
vanishes at all four distinct generic conductor roots.
Substitution of the Hensel root therefore gives

    P(z(v),c+v)=P(x(v),c+v).                              (4)

The right side is regular. Equation(4) is a MONIC polynomial
equation for z over C[[v]], because P is monic in x. Thus
z belongs to C[[v]]. Directly, if z had negative order, the
leading monomial of P(z,c+v) would strictly dominate all lower
monomials with regular coefficients and could not cancel.

Consequently h=(x+z)/2 is regular at c. This works at EVERY zero
of a_4; all other finite points were already regular in(3).
Hence

    h in C[y].                                           (5)

Make the single polynomial source shear u=x-h(y).
It preserves the nonzero constant source Jacobian. Generic central
symmetry now makes the entire conductor even:

    kappa=a(y)*u^4+b(y)*u^2+c_0(y).                      (6)

The constant coefficient c_0(y) is nonzero on every fibre. For a
degree-four fibre, the hyperelliptic node or tacnode normalization
has nonzero radii: a zero radius would ramify the normalization,
contradicting immersion. For a degree-two fibre, the even quadratic
has two distinct conductor roots; if c_0 vanished they would be
the single point u=0, again impossible.

Thus c_0 is a nonzero CONSTANT, which we denote by c. Also a and b
cannot vanish together, since that would give a constant conductor
fibre. We have proved the exact survivor shape

    kappa=a(y)*u^4+b(y)*u^2+c,
    c in C*, gcd(a,b)=1, a nonconstant.                  (7)

The use of the ACTUAL monic target coordinate in(4) is essential.
This is not a consequence of arbitrary plane-open completion
geometry or of an arbitrary regular source coordinate.

## 3. Two periods with variable product

On the generic source fibre normalize(7) by its nonzero leading
coefficient. Put

    S=-b/a, K_0=c/a,
    A+B=S, A*B=K_0, A!=B,
    r^2=A, s^2=B.

The four conductor parameters are {-r,r,-s,s}. The canonical
genus-two plane generators are

    p_0=u^2,
    q_0=u*(u^2-A)*(u^2-B).

The original centered affine-line normalization parameter differs
from the canonical one only by a nonzero scale; a diagonal target
scaling absorbs it. Thus a polynomial target automorphism on each
generic fibre carries the actual generator pair to this canonical
pair in the SAME source parameter u.

Matching-endpoint action integrals are all multiplied by the SAME
constant nonzero Jacobian under a polynomial target automorphism.
The canonical periods give a common algebraic scale Lambda:

    I_r=integral_(-r)^r P*Q_u*du=Lambda*r*C_A,
    I_s=integral_(-s)^s P*Q_u*du=Lambda*s*C_B,

    C_A=3*A^3-7*K_0*A,
    C_B=3*B^3-7*K_0*B.                                  (8)

The constant 8/105 in the canonical integrals is absorbed into
Lambda. Neither r nor s is zero and C_A,C_B cannot both vanish.
Interchange the labels if necessary so that C_A is not identically
zero; its nonvanishing open set is used for the first quotient below.

These are identities of algebraic functions. Define Lambda by
the first quotient on a nonempty open set; pointwise target
equivalence gives the second identity there. This does not posit
a target automorphism depending algebraically on y.

Write J=P_u*Q_y-P_y*Q_u in C*. The exact moving-endpoint identity is

    I_r'=-2*J*r, I_s'=-2*J*s.                             (9)

The boundary term is P(r)*(d/dy)Q(r)-P(-r)*(d/dy)Q(-r),
which vanishes because both coordinates match. Derivatives are
taken in the separable algebraic function field.

## 4. The exact two-scale equations

Differentiate A^2-S*A+K_0=0. Combining this with(8),(9) gives

    Lambda'*C_A
       +(7/2)*Lambda*(3*A^2*S'-5*A*K_0')=-2*J,           (10)

and the same identity with B. The useful derivative identity is

    C_A'+(A'/(2*A))*C_A
       =(7/2)*(3*A^2*S'-5*A*K_0').

Let N=3*S^2-10*K_0. Since

    C_A-C_B=(A-B)*N,
    C_A+C_B=S*(3*S^2-16*K_0),

subtracting(10) and then eliminating Lambda' from its sum yields

    2*N*Lambda'+7*(3*S*S'-5*K_0')*Lambda=0,
    105*K_0*Lambda*(2*K_0*S'-S*K_0')=-4*J*N.              (11)

All factors and constants are retained. Substitute S=-b/a,K_0=c/a
and define the polynomials

    D=3*b^2-10*c*a,
    H=b*a'-2*a*b',
    E=3*a*b*b'-(3*b^2-5*c*a)*a'.

Then(11) is exactly

    2*a*D*Lambda'+7*E*Lambda=0,
    105*c^2*Lambda*H=-4*J*a^2*D.                         (12)

These identities are initially generic and continue as rational
function identities wherever the expressions are defined.

## 5. A zero of the leading coefficient would have order 2/5

Choose ANY zero of a. With v=y-c_1, write

    a=a_m*v^m+..., b=b_0+...,
    m>=1, a_m*b_0!=0.

The nonzero b_0 follows from gcd(a,b)=1. Directly,

    D=3*b_0^2+O(v),
    H=m*b_0*a_m*v^(m-1)+...,
    E=-3*m*b_0^2*a_m*v^(m-1)+....

In particular H and D are nonzero rational functions. The second
identity in(12) makes the scale EXPLICITLY RATIONAL:

    Lambda=-4*J*a^2*D/(105*c^2*H).

Therefore Lambda has exact order m+1, with nonzero leading
coefficient lambda_0. The first identity in(12), at order 2m,
now has the nonzero common factor
b_0^2*a_m*lambda_0 and forces

    6*(m+1)-21*m=0,
    m=2/5.                                               (13)

This contradicts the positive INTEGER order of a polynomial zero.
The argument includes m=1 without a derivative-order exception.

Thus a has no zero and is constant, contradicting(7).
This proves the entire degree-four exclusion(1).

For an additional algebraic check, eliminating Lambda from(12)
also gives

    2*a*D*H'
       =H*[33*a*b*b'-(9*b^2+25*c*a)*a'].

The shorter first-order valuation argument(13) is sufficient for
the theorem and avoids unnecessary exceptional-case divisions.

## 6. Exact scope and replay

The checker verifies the variable-product period derivatives,
both equations(11), their polynomial normalizations(12), the
optional second-order identity, and the all-integer local
coefficient obstruction. It checks exact Hensel/mate controls
showing where monicity is used.

A non-Keller control

    P=u^2, Q=u*((y-1)*u^4+u^2+1)

has conductor degree four generically and degree two at y=1.
Its normalization is immersed on every fibre, its centroid and
constant-term constraints pass, and its Jacobian is 2*u^6.
Thus the geometric normal form(7) alone is not claimed impossible.
The actual constant nonzero Jacobian supplies the contradiction
through the two period equations.

No bounded target or source degree, uniform family target change,
trace-integrality assumption or discarded original source prime
is used. The remaining noninvertible class has generic conductor
degree at least SIX. This does not prove the absence of that class
or automatic trace perfection, HC4 or JC2.

Exact controls support the written proof and the audited primary
classification. Independent replay, formal verification and
external human review are not claimed.
