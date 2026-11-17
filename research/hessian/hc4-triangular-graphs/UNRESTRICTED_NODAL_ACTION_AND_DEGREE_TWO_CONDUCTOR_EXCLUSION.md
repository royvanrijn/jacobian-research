# Unrestricted nodal action and degree-two conductor exclusion

Let

    A_0={f in C[t]: f(1)=f(-1)}
       =C[t^2-1,t*(t^2-1)].

For EVERY polynomial generator pair A_0=C[p,q],

    I(p,q)=integral_(-1)^1 p(t)*q'(t)*dt != 0.              (1)

There is no bound on either degree. The essential imported result is
Neumann's unconditional classification of polynomial curves with one
node, checked in the original source. Combining(1) with the endpoint
variation identity excludes the ENTIRE degree-two conductor case for
a noninvertible plane Keller map. Automatic trace integrality and
the higher-conductor cases remain unresolved.

This resolves the unrestricted action obligation left in the
[low-degree proof](LOW_DEGREE_NODAL_ACTION_AND_CONDUCTOR_EXCLUSION.md)
and the [quartic proof](QUARTIC_NODAL_GENERATORS_AND_PAIR_NORMS.md).
Those elementary arguments remain independent low-degree proofs.

## 1. Exact imported statement and its applicability

Neumann, *Complex algebraic plane curves via their links at infinity*,
Invent. Math.98(1989),445-489, Section8, pages482-483, proves the
following particular classification:

> A reduced irreducible polynomial curve whose normalization is the
> affine line and whose only singularity is one ordinary node is
> carried to Y^2=X^3+X^2 by a polynomial automorphism of C^2.

This is a paraphrase of the conclusion, not a quotation. Page482
specifies polynomial automorphisms in the low-genus classification;
page483 explicitly applies that classification to the singular nodal
curve. The original source, including its relevant proof, is retained
[here](sources/neumann-1989-plane-curves-links-infinity.pdf).
Its exact page locators and dependency audit are retained in the
[source audit](sources/neumann-nodal-classification-source-audit.json).

Every plane embedding supplied by A_0=C[p,q] has the stated properties.
Indeed C[t] is the finite normalization, its only identification is
the transverse gluing of t=1 and t=-1, and Spec(A_0) has no other
singularity. A nonconstant polynomial p also shows that the
normalization map is finite. At infinity its projective normalization
is P^1 with just the single omitted point t=infinity, so the curve has
one place at infinity.

The source argument does not assume the singular curve itself is a
regular fibre. The one place at infinity makes its link a knot.
The CORRECTED Neumann-Rudolph Lemma7.1 says that a reduced fibre
with a knot at infinity has a defining polynomial all of whose
fibres are regular at infinity. Its statement and proof were checked
on pages349 and351 of the
[1988 corrigendum](sources/neumann-rudolph-1988-unfoldings-corrigendum.pdf).
Smoothing the one node then gives a smooth once-punctured genus-one
fibre, to which Neumann's low-genus table applies.

The table gives a cubic Weierstrass equation after polynomial target
changes. For the given nodal fibre, its cubic has a double root r and
a distinct simple root s. Thus its equation has the form

    Y^2=(X-r)^2*(X-s), r!=s.

Translation and nonzero complex scaling put it into
Y^2=X^2*(X+1). This elementary factorization justifies the final
normal form without using a mistaken uniqueness sentence on page483:
X^3-Y^2+X^2 has TWO singular fibre values, 0 and4/27. The prescribed
nodal fibre nevertheless has the displayed normal form.

No conditional classification from the later parametric-lines paper
is used.

## 2. Polynomial target changes scale the action

Let Phi=(F,G) be a polynomial automorphism of C^2 and let its
constant nonzero Jacobian be lambda. The polynomial one-form

    F*dG-lambda*X*dY

is closed. It is polynomial-exact: a closed polynomial form
a*dX+b*dY has the polynomial primitive

    H(X,Y)=integral_0^1
           [X*a(tau*X,tau*Y)+Y*b(tau*X,tau*Y)]*dtau.

The ordinary formal differentiation identity verifies dH=a*dX+b*dY.
Consequently, if p and q have matching values at the two endpoints,

    I(F(p,q),G(p,q))=lambda*I(p,q).                        (2)

The exact correction H(p,q) has equal endpoint values.

Choose Phi from the imported classification, and write

    (s(t),z(t))=Phi(p(t),q(t)).

The target curve is z^2=s^2*(s+1). Its normalization parameter tau
satisfies s=tau^2-1,z=tau*(tau^2-1). The induced isomorphism of
normalizations is an automorphism of C[t], so tau=a*t+b with a!=0.
The two preimages of the node must match the unordered sets
{1,-1} in the two normalizations. Hence b=0 and a=epsilon in{1,-1}:

    s=t^2-1, z=epsilon*t*(t^2-1).

Direct integration gives

    I(s,z)=epsilon*integral_(-1)^1(3*t^4-4*t^2+1)*dt
          =8*epsilon/15.

Equation(2) proves the unrestricted formula

    I(p,q)=8*epsilon/(15*lambda) != 0.                     (3)

The automorphism may depend on the individual generator pair.
No uniform target automorphism for a family is required.

## 3. The entire degree-two conductor case is excluded

Use the actual plane-Keller setup and the exact conductor theorem in
[the generic-coordinate proof](GENERIC_COORDINATE_ORDER_CONDUCTOR_AND_CLOSED_COLLISIONS.md).
If its global source conductor polynomial has directional degree TWO,
a polynomial source shear and a constant rescaling put the order into

    C[P,Q,y]=C[y,u^2-1,u*(u^2-1)].                         (4)

The Jacobian

    J=P_u*Q_y-P_y*Q_u

is a nonzero constant in these source coordinates. Each specialization
y=c of(4) is EXACTLY A_0=C[P(u,c),Q(u,c)], so(1) applies on EVERY fibre.
In particular,

    I(y)=integral_(-1)^1 P(u,y)*Q_u(u,y)*du

is a polynomial with no complex zero.

On the other hand the two target functions and their y-derivatives
match at u=1,-1. Differentiating and integrating by parts gives

    I'(y)=[P*Q_y]_(u=-1)^(u=1)
           -integral_(-1)^1(P_u*Q_y-P_y*Q_u)*du
          =-2*J.                                         (5)

Since J!=0, I is affine with exactly one complex zero. This contradicts
unrestricted nodal-action nonvanishing. Therefore the degree-two
conductor case cannot occur.

Together with the earlier evenness and positive-degree result, a
noninvertible plane Keller map must have global directional conductor
degree at least FOUR. This conclusion uses the generic source-finite
coordinate and conductor hypotheses of the linked proof. It does not
assert that arbitrary degree-two unramified normalization orders are
Keller maps, or that higher conductor degrees are impossible.

## 4. Exact evidence and remaining scope

The checker pins the two original source PDFs and the source audit.
It verifies the standard normalization and action, target-action
scaling for explicit polynomial automorphisms, degree-six examples
with unbounded added shear degree, the endpoint variation identity,
and the elementary double-root cubic normalization. Controls retain
a zero-action polynomial immersion with TWO nodes, which is outside
the exact A_0 hypothesis, and check the source's two critical values.

The imported topological classification is established by the cited
written proof; the checker does not prove it. These are internal
source inspection and exact symbolic controls. Independent replay,
formal verification and external mathematical review are not claimed.
