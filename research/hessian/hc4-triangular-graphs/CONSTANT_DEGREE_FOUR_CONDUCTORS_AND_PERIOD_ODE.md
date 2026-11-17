# Constant degree-four conductors and the two-period equation

Let K=(P,Q) be an actual complex polynomial plane Keller map, in a
source frame where P is monic of positive x-degree. Let kappa generate
the coordinate-order conductor

    (C[P,Q,y]:C[x,y])=kappa*C[x,y].

It is IMPOSSIBLE that

    deg_x kappa(x,c)=4 for EVERY c in C.                  (1)

Repeated conductor roots are allowed in this statement. Thus every
global degree-four conductor survivor must have a degree drop on
some source line. The subsequent
[global degree-four theorem](GLOBAL_DEGREE_FOUR_CONDUCTOR_EXCLUSION.md)
excludes those degree-drop cases too. Larger conductor degrees remain open.

The proof strengthens the
[four-simple-root result](SIMPLE_FOUR_ROOT_CONDUCTORS_AND_TWO_NODE_ACTION.md).
It uses the same audited genus-two classification, but allows moving
node pairs and their tacnode collisions. A polynomial differential
equation rules out the moving case.

## 1. Exact fibre conductors and smooth normalization branches

The [coordinate-order theorem](GENERIC_COORDINATE_ORDER_CONDUCTOR_AND_CLOSED_COLLISIONS.md)
specializes its monogenic bases on every source line. Thus

    A_c=C[P(u,c),Q(u,c)] subset C[u]

has finite birational normalization C[u], its EXACT conductor is
kappa(u,c)*C[u], and

    dim_C(C[u]/A_c)=deg_u(kappa(u,c))/2=2.                 (2)

The same theorem's polynomial-immersion argument says that P_u,Q_u
cannot both vanish. Every branch of each fibre curve is therefore
smooth; its normalization map has nonzero tangent.

For completeness, the conductor orders determine the possibilities
even when its roots are repeated. At a singular image point, choose
a constant linear target coordinate xi whose derivative is nonzero
on every branch above that point. The branches are graphs

    eta=phi_i(xi).

Their reduced local equation, up to a unit, is the product of
eta-phi_i(xi). The monogenic trace-dual conductor formula on branch i
therefore gives conductor order

    c_i=sum_(j!=i) ord(phi_i-phi_j)
       =sum_(j!=i) intersection_multiplicity(i,j).        (3)

The derivative of xi on the normalization is a unit. Hence a point
with at least three branches consumes conductor length at least SIX.
Every singular point with two branches consumes twice their positive
intersection multiplicity.

The total conductor length is FOUR. The fibre therefore has either
TWO ordinary nodes or ONE ordinary tacnode: two smooth branches with
intersection multiplicity two. A cusp, ramphoid cusp or further
unibranch singularity would contradict the normalization immersion.
There are no other singularities.

## 2. The genus-two classification gives an even quartic

The [retained primary-source audit](sources/neumann-two-node-genus-two-source-audit.json)
pins Neumann's genus-two table and the corrected Neumann-Rudolph
one-place-at-infinity lemma. The latter makes the defining polynomial
good at infinity because its reduced fibre has affine-line
normalization with just one place at infinity.

In the two-node case, the pinched affine line has Euler characteristic
-1, and smoothing two nodes changes this to -3. In the tacnode case,
the pinched affine line has Euler characteristic zero. The local
tacnode smoothing is a double cover of a disk with four branch points;
its Euler characteristic is -2, compared with one for the unsmoothed
local curve, so the total again becomes -3. The unchanged end is one
puncture. Thus a nearby smooth fibre has genus TWO in both cases.

Neumann's sole realizable genus-two row is hyperelliptic quintic,
up to polynomial ambient automorphisms. A two-node fibre consequently
has a normalization of the form

    X=tau^2+gamma,
    Y=tau*(tau^2-R_1)*(tau^2-R_2),
    R_1*R_2*(R_1-R_2)!=0.                               (4)

For a tacnode the same formula holds with R_1=R_2!=0. The nonzero
conditions also follow from immersion: if R_i=0, both derivatives
of X,Y vanish at tau=0.

The two affine-line normalization parameters differ by
tau=a*u+b. Hence the conductor roots in either case are centrally
symmetric, with a common midpoint. In the tacnode case each of the
two roots has multiplicity two.

Under(1), the leading coefficient of kappa has no complex zero and
is a nonzero constant. Make kappa monic. Its root centroid is the
polynomial h(y)=-[u^3]kappa/4. The one source shear u=x-h(y)
makes EVERY fibre conductor even, including the tacnode fibres.
Consequently the whole polynomial has the form

    kappa=u^4-S(y)*u^2+C(y), S,C in C[y].                 (5)

Every root pair in(4), repeated or not, has nonzero radius.
Thus C(c)!=0 on EVERY fibre, so C is a nonzero CONSTANT. Write

    C=K_0 in C*.                                         (6)

This unit conclusion uses immersion, not just conductor degree.

## 3. Common pointwise action scaling on the generic node cover

First suppose S is nonconstant. Over an algebraic function-field
extension choose a,b,r,s satisfying

    a+b=S, a*b=K_0, a!=b,
    r^2=a, s^2=b.

The generic fibre has nodes at {-r,r} and {-s,s}. Its canonical
plane generators are

    p_0=u^2,
    q_0=u*(u^2-a)*(u^2-b)
       =u*(u^4-S*u^2+K_0).                              (7)

The polynomial target automorphism giving(4) can be followed by a
translation and nonzero diagonal scaling to give EXACTLY(7) in the
original centered normalization parameter u. This is the same
argument as in the four-simple-root proof, applied fibrewise.
A different target automorphism is allowed on each fibre.

Polynomial target automorphisms multiply ALL matching-endpoint
action integrals by their same constant nonzero Jacobian. Therefore
the actual periods have a COMMON pointwise scale:

    I_r=integral_(-r)^r P*Q_u*du=Lambda*r*C_a,
    I_s=integral_(-s)^s P*Q_u*du=Lambda*s*C_b,             (8)

where the harmless constant 8/105 has been absorbed into Lambda and

    C_a=3*a^3-7*K_0*a, C_b=3*b^3-7*K_0*b.                (9)

Indeed the canonical periods are (8/105)*r*C_a and
(8/105)*s*C_b. These cannot both vanish because K_0!=0.

The periods in(8) are algebraic functions: evaluate a polynomial
antiderivative at the algebraic endpoints. Define Lambda=I_r/(r*C_a)
on the nonempty open set where its denominator is nonzero.
The pointwise classification then makes the second identity in(8)
an identity of algebraic functions. No target automorphism varying
algebraically with y is asserted or needed.

Let J=P_u*Q_y-P_y*Q_u be the nonzero constant source Jacobian,
unchanged by the centroid shear. The endpoints now MOVE, but
integration by parts still gives

    I_r'=-2*J*r, I_s'=-2*J*s.                             (10)

For example, the endpoint contribution for I_r is
P(r)*(d/dy)Q(r)-P(-r)*(d/dy)Q(-r), which vanishes because BOTH
target coordinates match at these endpoints. All derivatives in
(8)--(10) are taken in the separable algebraic function field.

## 4. Exact elimination gives a polynomial differential equation

Put

    d=a-b, N=3*S^2-10*K_0.

Since a*b=K_0,

    a'/a=S'/d, b'/b=-S'/d.

Divide the derivatives of(8) by r and s, respectively. Equations(10)
become

    Lambda'*C_a+Lambda*(C_a'+S'*C_a/(2*d))=-2*J,
    Lambda'*C_b+Lambda*(C_b'-S'*C_b/(2*d))=-2*J.           (11)

The following exact identities use only a+b=S,a*b=K_0:

    C_a-C_b=d*N,
    C_a+C_b=S*(3*S^2-16*K_0).

Subtracting(11) and dividing by d gives

    2*N*Lambda'+21*S*S'*Lambda=0.                         (12)

Their sum is

    S*(3*S^2-16*K_0)*Lambda'
       +(21/2)*(S^2-2*K_0)*S'*Lambda=-4*J.

Eliminate Lambda' using(12). The coefficient simplification is

    N*(S^2-2*K_0)-S^2*(3*S^2-16*K_0)=20*K_0^2.

Hence

    210*K_0^2*Lambda*S'=-4*J*N.                          (13)

For nonconstant S, N is not the zero polynomial. Differentiate(13)
and use(12), on the nonempty open set where the divisions are
defined. This gives the exact identity

    2*(3*S^2-10*K_0)*S''=33*S*(S')^2.                    (14)

Both sides are polynomials in y, so the identity holds globally.
A nonconstant polynomial S of degree n>=1 cannot satisfy it.
Comparison of the leading coefficients gives

    6*n*(n-1)=33*n^2,
    equivalently 3*n*(9*n+2)=0,

which is impossible. This includes n=1, for which the left side
of(14) vanishes and the right side does not.

## 5. Constant configurations are also impossible

It remains that S is constant. The two radii in(4) are then fixed,
with either two distinct node pairs or one repeated tacnode pair.
For distinct pairs the canonical periods

    (8/105)*r^5*(3*r^2-7*s^2),
    (8/105)*s^5*(3*s^2-7*r^2)

cannot both vanish. In the repeated case r^2=s^2, the tacnode
period is

    integral_(-r)^r u^2*d[u*(u^2-r^2)^2]
       =-32*r^7/105 != 0.                               (15)

Select a fixed nonzero canonical period. Pointwise target-action
scaling makes the actual corresponding period nonzero on EVERY
complex y. But the fixed-endpoint version of(10) makes it an affine
polynomial with nonzero slope, hence with a complex zero.
This is a contradiction.

Together with Section4, this proves that(1) never occurs.

## 6. Controls and exact remaining scope

The checker verifies the conductor and index of the whole canonical
family, both period formulas, the moving-endpoint derivative identity,
the polynomial elimination(11)--(14), and the degree contradiction.
It also retains a non-Keller immersion family

    P=u^2, Q=u*(u^4-y*u^2+1),

whose conductor has degree four on every fibre and passes through
tacnodes at y=2,-2. Its source Jacobian is -2*u^4. Thus the geometry
of such families alone supplies no contradiction; the constant
nonzero Keller Jacobian in(10) is essential.

The theorem now excludes degree-four conductors with repeated
roots as well, PROVIDED no fibre degree drops. A global degree-four
survivor must therefore have a special degree drop; by the retained
evenness and positivity theorem this drop is to degree TWO.
The subsequent global degree-four theorem excludes this remaining
degeneration, using monicity through the degree drop and the variable-product
period equations. Higher conductor degrees, automatic trace integrality,
HC4 and JC2 remain open.

The genus-two classification is an audited primary-source input.
Exact symbolic controls do not prove it or the unrestricted analytic
arguments. Internal algebraic audits and exact controls are recorded;
independent replay, formal verification and external human review
are not claimed.
