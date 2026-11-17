# Monic parameter pole bounds exclude cyclic target coordinates of even leading degree

PF2RULED149. 17 November 2026. ALL degrees over C.

On every quadratic source-image model in the stated class and on S2,
consider

    P=B(u,z)+w*A(u,z),
    deg_z B=n>=1, deg_z A<=n-1,
    B=c(u)*z^n+sum_(i<n) b_i(u)*z^i.               (1)

ALL lower coefficients of A,B are arbitrary polynomials. There is NO
bound on their u-degree or on the degree of a hypothetical companion.
If c is nonzero, coprime to M, and every zero of c has multiplicity
at most 2*n, then P has NO rational Jacobian mate, even in ANY finite
extension of the source-image function field.

In particular EVERY nonzero CONSTANT c is allowed. This includes ALL
above-face mixed u,z terms, arbitrary lower w coefficients and arbitrary
base terms within (1). It is not a reduction of all coordinates to (1).
Odd highest w degree, other leading coefficient divisors and full
first-derivative integrality remain open. HC4/JC2 is not closed.

## 1. General theorem with explicit leading coefficient conditions

Let V={w^2=N(u)+M(u)*z} be a smooth affine surface over C, with
M a nonzero polynomial and N ANY polynomial. Its residue form

    Omega=du wedge dw/M

is a regular unit on all V, using the usual smooth residue patches
also at M=0. Put m=deg M, b=deg c. For r in C let

    m_r=ord_r M, b_r=ord_r c.

For P as in (1), the following two conditions suffice:

    b+n*m >= 2*n;                                 (2)
    b_r+n*m_r <= 2*n for EVERY root r of c.        (3)

There is no restriction on N's degree. When m>=2, (2) is automatic.
The conclusion is absence of a rational mate dP wedge dQ=Omega,
also in every finite characteristic-zero function-field extension.

For the actual quadratic MX0 family, M=u^2*b0(u)^k with k>=2,
b0 square-free, b0(0)!=0, deg b0>=1. Every zero of M has multiplicity
at least two. Condition (3) is therefore EXACTLY the stated sufficient
condition: gcd(c,M)=1, and b_r<=2*n at every zero of c. The same
description applies to S2, where M=u^2,N=u.

We prove the general assertion, including equality in BOTH (2),(3).

## 2. A universal monic parameter lemma

Let k be an algebraic closure of C(p), and let partial_p be its unique
extension of d/dp. Its constants are C. Indeed if a nonconstant
lambda in k had derivative zero, p would be algebraic and separable
over C(lambda); differentiating its minimal polynomial would give
a zero derivative of that separable polynomial at p, a contradiction.

Let Q(T) be ANY monic polynomial of degree D>=1 over C((s)).
Let h in C((s))* be independent of p. On ANY branch suppose

    Q(W)=p*h, W in k((s)).                         (4)

Write nu for s-order. Then

    nu(partial_p W) >= ceil(nu(h)/D).              (5)

Here partial_p fixes s and acts coefficientwise; no hypothesis on
simple roots or root multiplicities is needed.

Proof. There is a FIRST coefficient of the Laurent expansion of W
that is not in C. Otherwise both Q(W) and h would be in C((s)),
contradicting (4). Let its exponent be r and its coefficient lambda.
All earlier terms form a finite Laurent polynomial W0 in C[s,s^-1].
Thus

    W=W0+s^r*(lambda+O(s)), lambda not in C.

The monic polynomial Q(W0+s^r*Y) has coefficients in C((s)) and
leading Y coefficient s^(D*r). Let q be the minimum order of all
its coefficients. Its nonzero initial polynomial H(Y) is in C[Y],
and q<=D*r.

If D*r<nu(h), then q<nu(h). Taking the leading coefficient in (4)
gives H(lambda)=0. Algebraic closedness of C forces lambda in C,
a contradiction. Hence D*r>=nu(h). Since the earlier coefficients
are constant and partial_p lambda!=0, one has
nu(partial_p W)=r. This proves (5).

This is a coefficient argument for EVERY normalized curve branch.
It permits arbitrarily large p-independent polar terms in W and
unbounded degree in the coefficients of Q. Such polar terms do not
contribute to partial_p W.

## 3. Actual generic coframe and the monic curve equation

Work on the WHOLE geometric generic fibre P=p over k, componentwise
if disconnected. It is smooth: on a smooth characteristic-zero
surface every critical curve maps to one value, and the other
critical points are finite. The unit Omega defines a regular,
NONZERO relative differential omega there by

    Omega=dP wedge omega.

This also covers every finite affine point with M=0, not just the
rational coordinate patch below.

On M*c!=0, substitute z=(w^2-N)/M and write P*(u,w) for the
result. The polynomial in w

    Q_u(T) = [sum_(i=0)^n b_i(u)*M(u)^(n-i)
                            *(T^2-N(u))^i
              +T*sum_(i=0)^(n-1) a_i(u)*M(u)^(n-i)
                            *(T^2-N(u))^i] / c(u)  (6)

is MONIC of degree D=2*n over C(u). Its coefficients depend only
on u, not on p. The generic equation is exactly

    Q_u(w)=p*h(u), h(u)=M(u)^n/c(u).               (7)

The degree 2*n is a degree over the CURVE parameter u. It is NOT
the original target extension degree, and (6) is NOT a monic
equation proving source integrality over the target B.

On this patch the actual coframe is

    omega=-du/(M*P*_w)=-(partial_p w)*du/M,         (8)

where partial_p fixes u. Differentiate P*(u,w)=p to see
P*_w*partial_p w=1. All mixed partials of A,B are therefore retained;
no additive coframe, omitted F_u term or target symmetry is assumed.

Every generic component has nonconstant u. A component with constant
u would lie in a fibre u=r where the highest coefficients vanish,
forcing r in C; a whole component of that fibre can only map to a
constant in C, not the transcendental p. Equivalently, on M*c!=0
the leading coefficient in (6) precludes a vertical component;
no component of the generic fibre can be contained in the fixed
divisors M*c=0. Thus the local curve parameters used below exist.

## 4. Every infinity point over u=infinity

At any normalized projective generic-curve point over u=infinity,
choose parameter s, absorbing a formal unit root, so

    u=s^(-alpha), alpha>=1.

The constant field is k, but u now has coefficients in C and is
fixed by partial_p. Equation (7) is of the form (4), with

    nu(h)=(b-n*m)*alpha.

Using (5) in (8) gives the exact expression and lower bound

    ord_s omega=(m-1)*alpha-1+nu(partial_p w)
      >= ceil(((m-2)/2+b/(2*n))*alpha)-1 >= -1.    (9)

The last inequality is (2). If (2) is strict, the bound is
nonnegative at every such point. Arbitrary leading cancellations,
multiple roots, large polar parts of w and EVERY ramification
index alpha are included.

This argument needs no degree gap between N and M, no Newton-face
condition on lower coefficients and no bound on the order at which
p first enters a Puiseux expansion.

## 5. Every finite-u point at infinity

If u has finite value r and c(r)!=0, there is NO missing infinity
point. To verify this, a pole of z of order beta>0 would give w
pole order at most beta/2 from w^2=N+Mz, because M,N stay finite.
The term c*z^n would have pole order n*beta. Every lower B term
has smaller order, and every w*A term has order at most
(n-1)*beta+beta/2<n*beta. They cannot cancel it in P=p.
If z is finite, w is also finite, so the point belongs to the
whole affine curve.

Thus every finite-u infinity point lies over a root r of c,
which is a CONSTANT in C. On each normalized branch choose

    u=r+s^alpha, alpha>=1.

Equation (7) again has p-independent Laurent coefficients and

    nu(h)=(n*m_r-b_r)*alpha.

Since ord du=alpha-1 and ord M=m_r*alpha, (5),(8) yield

    ord_s omega=(1-m_r)*alpha-1+nu(partial_p w)
      >= ceil((1-m_r/2-b_r/(2*n))*alpha)-1 >= -1.  (10)

This is (3), including EQUALITY. If the inequality in (3) is
strict at r, the lower bound is nonnegative. This treats all
branches over that value, with finite OR polar w, any ramification
and arbitrary degeneracy of the specialized curve polynomial.
In particular no original ramified boundary prime is removed
because a curve equation has a vanishing leading coefficient.

## 6. At most simple poles forbid every rational primitive

The entire differential omega is now regular on the affine generic
curve and has AT MOST SIMPLE POLES on its projective normalization,
by (9),(10).

A nonzero rational exact differential on a smooth complete curve
must have a pole of order at least TWO: any pole of its primitive
of order r>=1 gives differential order -r-1. If the primitive has
no poles, it is a regular function on a complete curve and hence
constant. Thus omega is not rationally exact. This argument
includes both holomorphic omega and omega with logarithmic poles;
it does not assert holomorphicness when (9) or (10) is sharp.

The conclusion persists in EVERY finite characteristic-zero source
extension. A ramified pullback has order

    E*ord(omega)+(E-1).

For ord>=0 this is nonnegative, and for ord=-1 it is exactly -1.
The pulled-back nonzero differential still has no pole of order
at least two and hence no rational primitive. Alternatively,
normalized FIELD trace of a hypothetical primitive gives one
downstairs. Neither argument asserts polynomial trace landing.

All geometric generic components, all normalized branches and
all ramification indices are retained. Three Puiseux conjugates
of one branch are not counted as three unramified points.

## 7. Applications and exact boundary controls

Every quadratic MX0 model and S2 has the complete smooth affine
source image and whole-plane cover retained in
the [complete quadratic source-image presentation](CANONICAL_DECK_JET_FIELDS_AND_RATIONAL_POLE_ATLAS_COVERAGE.md#7-complete-quadratic-image-presentation-beyond-mx0).
The pulled-back volume is a nonzero constant multiple of the source
plane volume, as retained in
the explicit source pullback in that presentation.
Consequently (1)--(3) exclude ACTUAL target coordinates, even
when their hypothetical companion does not descend from the
whole plane source field to the source-image field.

For constant c, the conclusion now allows EVERY polynomial B
monic up to a nonzero constant in z, and EVERY polynomial
A with deg_z A<n. No Newton support restriction remains.
This strictly enlarges the displayed forms of PF147/PF148;
no claim that all surface-automorphism transports have been
classified is made. The regular linear-z cases were already
subject to PF52's source-line exclusion; the no-rational-mate
conclusion here is stronger and must not be advertised as
a new regular linear-z frontier.

A genuinely above-face S2 example is

    P=w+z^2+u*z+u.

It has degree two on BOTH named source axes, under
u=t^2,w=t*(1+t^2*h),z=2*h+t^2*h^2.
The term u*z is outside PF148's n=2,m=2 region.
The present all-branch pole bound, not an axis census, excludes
a rational mate.

On the actual MX k=2,b0=u-1 model, use

    M=u^2*(u-1)^2, N=u*(1-3*u)^2,
    P=z^2+(u^9+1)*z+u^15
        +w*(u^7*z+u^3+1).                        (11)

It has arbitrary high mixed terms, a nonconstant w coefficient
depending on z, and constant leading c=1. At every marked
fixed-t source line its degree is two at a(t)=0 and four
otherwise. This passes the old necessary line-degree gate;
it asserts no test of all polynomially recharted source curves.

For genuinely nonconstant leading c, any square-free polynomial
coprime to M is permitted in ANY degree. Repeated zeros up to
order 2*n are also allowed. There is no bound on their number.

Logarithmic equality is REAL, not a defect to erase. On S2,

    P=(u-1)^2*z

has n=1,c=(u-1)^2, with b_1=2=2*n. At its finite-u infinity
points the required differential has residues
-1/(2*sqrt(p)) and +1/(2*sqrt(p)), both nonzero.
So (10) can give a simple pole, while nonexactness still follows.

For a nonlinear, geometrically disconnected control, put

    H=z+w+u^2/4,
    P=H^2
     =z^2+(3*u^2/2)*z+u+u^4/16
         +w*(2*z+u^2/2).                         (12)

This has constant leading c and fits n=2. Its two geometric
generic components H=+sqrt(p), H=-sqrt(p) inherit the nonzero
logarithmic form divided by 2*H. Both are retained. A general
holomorphicness claim on S2 would be false.

The finite leading-coefficient condition cannot simply be dropped:
on S2, P=w+u*z has the rational mate Q=1/w. Here n=1,c=u,
m_0=2,b_0=1, so (3) fails, while (2) holds.
This polar pair is NO whole-plane counterexample.

The global condition is also material: for M=1,N=u the plane
pair P=w+z^2,Q=z has unit Jacobian. Its c is constant and
(2) fails. These controls prevent promoting the theorem to
arbitrary functions or arbitrary smooth affine surfaces.

## 8. Original normalization and remaining full goal

Retain the original target B, whole plane source R, FULL finite
normalization S and trace-dual Lambda, with every ramified and
unramified boundary prime, from
[PF143](MARKED_QUADRATIC_VELOCITIES_AND_CYCLIC_TRACE_POLES.md).
The monic polynomial (6) is a generic-curve device over C(u).
It is NOT integrality over B or a new finite target normalization.
The coefficientwise parameter derivative is NOT a source action
asserted to preserve S or Lambda. Field trace is used only for
rational exactness, never Tr(R) subset B or full marked regularity.

[PF142](PARTIAL_KUMMER_PLANE_CHARTS_AND_CYCLIC_EQUIVARIANT_FACTORS.md)'s mixed
2+3 minimal alternative remains retired. The present result
excludes a substantial surviving cyclic target family, including
all displayed even monic forms and the stated nonconstant
leading coefficients. No reduction of arbitrary target coordinates
to this family is supplied. Odd highest degree after clearing M,
leading divisors failing (3), and the actual first-derivative
integrality implication remain open. Full HC4/JC2 and an actual
counterexample remain OPEN.

The universal proof is the monic coefficient lemma (4),(5) and
the COMPLETE valuation audit (6)--(10). Supporting symbolic
controls are not substituted for this proof. No new external
theorem import or predecessor checker replay is used.
Independent whole-proof replay, formal verification, external
review and novelty priority remain unclaimed.
