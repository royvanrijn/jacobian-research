# Formal parameter coefficients force primitive equations

PF2RULED173. 17 November 2026. Actual complex polynomial plane Keller maps,
within the complete marked quadratic source model.

The leading coefficient of a ruling-six coordinate imposes an exact
polynomial differential equation, independently of the degree of its
actual companion. Every survivor must admit rho in C[u] satisfying

    rho_u - (G_u/(3G))*rho = 1, G=M^2*h.                  (1)

Consequently D=deg b3 is divisible by three, so the retained lower bound
D>=22 strengthens to D>=24. A multiplicity-three root of G is impossible.
A polynomial G with at least two roots cannot be a perfect cube or have
an affine rotational symmetry of order other than three. In particular
the entire numerical degree D=27 is excluded, and the Kummer extension
has degree exactly three in every actual survivor.
The equation excludes coefficient choices which pass the former numerical
ledger. With PF171, the next rational-character coefficient also forces
A*du/(hM) to have a rational primitive, where A is the coefficient of z
in the quadratic remainder. Neither condition closes the entire class.

## 1. Actual inputs and the formal inverse

Use the actual fields, original source and volume form of PF165:

    P,Q in C(u)[w], deg_w P=6,
    dP wedge dQ=du wedge dw/M,
    p6=(h/M)^2, M,h in C[u] nonzero,
    G=M^2*h, m=deg M, H=deg h, D=m+2H.                 (2)

All companion degrees and all coefficient degrees are unrestricted.
The [square leading form](RULING_SIX_SIMPLE_ESCAPE_GRAPHS_AND_QUARTIC_FACE_REDUCTION.md)
supplies (2) and at least one actual A pole root away from M.
There is also the central root of M. Hence G has at least TWO distinct
roots in every actual survivor considered here. The
[retained-source bound](RULING_SIX_RETAINED_SOURCE_PACKETS_AND_STRICT_INFINITY_BOUNDS.md)
supplies D>=22. The regular square root of PF171 is compatible with the
argument but is not needed for (1).

Let Z be a new formal variable, independent of u, and choose a with

    a^3=M/h, K=C(u,a), y=M/a, y^3=G.

The extension K/C(u) has degree one or three. The scaled equation

    Z^-6 * P(u,Z*V)=1

has reduction p6*V^6=1 at Z^-1=0. Its derivative at V=a is
6*p6*a^5=6/a, a unit in K. Formal implicit recursion gives

    w=a*Z+c0+c1*Z^-1+... in K((Z^-1)), P(u,w)=Z^6.        (3)

Every coefficient belongs to K. Substitution into the ACTUAL polynomial
Q gives a Laurent series with finite positive part,

    Q(u,w(u,Z))=sum_j q_j(u)*Z^j, q_j in K.

Extend the u derivation uniquely to K, and differentiate with Z fixed.
Using the actual Jacobian identity in (2),

    partial_u Q(u,w(u,Z)) = -1/(M*P_w(u,w(u,Z))).          (4)

The right side starts with

    -(a/(6M))*Z^-5 = -(1/(6y))*Z^-5.

Thus the coefficient q_-5 satisfies

    partial_u(-6*q_-5)=1/y.                              (5)

This is an algebraic primitive on the cubic Kummer function field.
It is forced by a coefficient of the full actual companion, rather than
by a truncation or a proposed target change.

## 2. Descent to a polynomial differential equation

Put rho0=-6*y*q_-5 in K. Since y_u/y=G_u/(3G), equation (5) gives

    partial_u rho0 - (G_u/(3G))*rho0=1.

Normalized field trace K/C(u) commutes with the extended derivation.
Applying it produces rho in C(u) satisfying (1). This also covers the
case that G is a cube and K=C(u); no irreducible Kummer hypothesis is
required. Conversely, a rational solution of (1) makes rho/y an
algebraic primitive of du/y.

In fact rho has no finite pole. At a root c of G of multiplicity n>0,
a hypothetical rho of negative order k has leading term in the left
side of (1) with coefficient k-n/3, which is nonzero, and order k-1<0.
At a point off G the same argument has coefficient k. Both contradict
the regular right side. Therefore rho is a POLYNOMIAL.

At a root c of G, rho(c)=0: otherwise the logarithmic coefficient in
(1) would give a pole. If its zero order exceeded one, both terms on
the left would vanish at c. Thus each G-root is a simple rho-root,
and comparison of constant coefficients gives

    rho_u(c)=3/(3-n).                                    (6)

In particular n=3 is impossible. Roots of rho off G are simple as
well, with rho_u(c)=1. These assertions follow locally without an
upper bound on any multiplicity.

## 3. Degree congruence and the new lower bound

Let N=deg G=2m+H and e=deg rho. Since G has at least two distinct
roots, rho has degree at least two. At infinity the leading term of
the left side of (1) is

    (e-N/3)*leading(rho)*u^(e-1).

It must vanish because the right side is constant. Hence

    e=N/3 is an integer, H congruent to m modulo 3,
    D=m+2H congruent to 0 modulo 3.                       (7)

Combined with D>=22 this proves the unrestricted necessary bound

    D>=24.                                               (8)

The two-root hypothesis matters. For G=u^n with n!=3,
rho=3u/(3-n) solves (1), even when n is not divisible by three.
These single-root controls are outside the actual A-plus-central
configuration.

## 4. Rotational symmetry obstruction

Suppose, after translating u=c+v, that G has a nontrivial rotational
symmetry

    G(c+zeta*v)=G(c+v), zeta of exact order ell>1.         (9)

The logarithmic coefficient transforms with weight minus one. Thus
for any solution rho, each polynomial

    zeta^-j * rho(c+zeta^j*v)

still solves (1) in v. Their average rho1 is a solution containing
only exponents congruent to one modulo ell. The two-root argument
forces its degree to be e=N/3, so e is congruent to one modulo ell.
But (9) also gives ell dividing N=3e. Consequently ell divides three:

    ell=3.                                               (10)

In particular NO reflection-invariant G with at least two roots can
occur. This includes reflection after any affine translation; no
reality assumption on the coefficients is used.

More generally, if G(c+zeta*v)=lambda*G(c+v) for a nonzero constant
lambda, the same projection proves N/3 congruent to one modulo ell.
The stronger conclusion ell=3 uses the invariant form (9).

## 5. Proper exclusions and surviving controls

The new equation is stronger than the former numerical ledger.
Take

    M=u^2, h=u^6*(u-1)^10*(u-2)^10,
    G=[u*(u-1)*(u-2)]^10.

Then m=2, H=26, D=54. Two candidate A slopes five give S=10,
n_A=2 and the old numerical bounds meet at b=10:

    ceil((D+3m)/6)=10=2S-5*n_A.

Nevertheless G is reflection-invariant around u=1, so (1) has no
rational solution. This excludes these leading coefficients for an
actual Keller map. This coefficient control is not asserted to be an
actual realization of the other local equations.

The equation has nontrivial positive controls. For G=(u^3-1)^10,

    rho=u-(5/2)*u^4+(15/7)*u^7-(9/14)*u^10

solves (1). Its order-three rotational symmetry is allowed by (10).
After translation, one can assign a central M-root among its three
roots, so this is not a single-root exception.

The new minimum D=24 is NOT excluded by (1). Indeed

    M=u^2, h=u*(u-1)^10, G=u^5*(u-1)^10,
    rho=(-42u+420u^2-945u^3+810u^4-243u^5)/28

satisfies it exactly. Its numerical data m=2, D=24, one A slope five
and b=5 meet the former ledger bounds. No actual companion or Keller
realization is claimed for either positive control.

## 6. Perfect cubes, the degree-27 ledger and the two-root case

### A multi-root cube is impossible

Suppose G=J^3 up to a nonzero constant, absorbed in J over C.
Then (1) makes F=rho/J a rational function with

    F_u=1/J.                                            (11)

If J has a simple root, the right side has a nonzero residue and
cannot be a rational derivative. Otherwise write the distinct roots
c_1,...,c_s of J with multiplicities n_i>=2; here s>=2 because G
has at least two roots. At c_i a rational primitive F has a pole
of exact order n_i-1, and it has no other finite poles. Since
1/J decays to order at least two at infinity, subtracting a constant
makes F vanish there. In lowest terms F=A0/B0 is then proper, with

    deg B0=deg J-s.

The exact leading order of (11) at infinity gives

    deg A0-deg B0-1=-deg J, so deg A0=1-s<0.

This is impossible for a nonzero numerator. Thus G is not a cube.
In particular K/C(u) in Section 1 has degree EXACTLY THREE in every
actual survivor. The general trace proof above still covers degree
one algebraically, without assuming this conclusion in advance.

### The entire degree D=27 is excluded

Suppose D=27. The marked source model has m>=2, and D=m+2H
makes m odd, so m>=3. Let n_A>=1 count the distinct A-supporting
base values and let S be the sum of their slopes. The complete ledger gives

    (9+m)/2 <= b <= 2S-5*n_A <= H-5*n_A,
    H=(27-m)/2.

Consequently m+5*n_A<=9, forcing m=3, n_A=1, H=12.
Now S<=H/2=6 and b>=6 force S=6. There is one A root c off M,
with slope six, which uses the full degree of h:

    h=constant*(u-c)^12.

Every M-root has multiplicity at least two in the marked model,
and the central root has multiplicity at least two. A degree-three
M therefore has only its central root, so M=constant*u^3.
These multiplicity facts are part of the complete model retained in
[PF162](RULING_SIX_FINITE_END_COVERAGE_AND_QUADRATIC_BOUNDARY_ALGEBRA.md).
It follows that G=constant*u^6*(u-c)^12 is a cube with two roots,
contradicting the preceding theorem. Thus

    D!=27.                                               (12)

This excludes a whole degree, not only one coefficient control.
The D=24 control in Section 5 remains allowed by the leading ODE.

### Exact leading-ODE classification when G has two roots

After affine translation and scaling, let

    G=u^a*(u-1)^b, a,b>=1, a+b=3e>=6.

For these two-root polynomials, the congruence and noncube conditions
are sufficient as well as necessary for (1). To see sufficiency,
assume 3 does not divide a; it then does not divide b. The Kummer
field is parametrized by

    u=t^3/(t^3-1), y=t^a/(t^3-1)^e,
    du/y=-3*t^(2-a)*(t^3-1)^(e-2)*dt.

This field is C(t): t^3 and t^a belong to C(u,y), and gcd(3,a)=1.
The last differential is a Laurent polynomial with no t^-1 term,
because every exponent is 2-a+3k and 3 does not divide a. It has
a Laurent-polynomial primitive. Multiplication by y gives the
following polynomial solution of (1) directly:

    rho=-3*sum_(k=0)^(e-2)
       binom(e-2,k)*(-1)^(e-2-k)
       *u^(k+1)*(u-1)^(e-k-1)/(3k+3-a).                  (13)

All denominators are nonzero. Necessity is the degree theorem and
the multi-root cube obstruction: when a+b is divisible by three,
G is a cube exactly when a and b are both divisible by three.
This is a classification of the LEADING ODE only. It constructs
neither a target companion nor a solution of the other residue
conditions.

## 7. The next rational character constrains the quadratic remainder

Use the additional exact regular normal form of
[PF171](RULING_SIX_REGULAR_CUBIC_APPROXIMATE_ROOT_AND_QUADRATIC_REMAINDER.md):

    P=H3^2+A(u)*z+B(u)*w+C(u), z=(w^2-N)/M,
    deg_w H3=3, leading_w(H3)=h/M.                     (14)

Here A is the quadratic REMAINDER coefficient; it is unrelated to the
marked coefficient used in M=u*A in earlier notes. All coefficients
in (14) belong to the actual source-image algebra, and A,B,C are
polynomials in u. No target-field membership is assumed.

Write the inverse (3) as w=sum_j a_j*Z^j, with a_1=a. We claim the
exact negative coefficient is

    a_-3=-A/(6h) in C(u).                                (15)

First invert the PURE cubic H3(w0)=Z^3 with the same leading term aZ.
Its three formal roots are w0(Z), w0(zeta3*Z), w0(zeta3^2*Z).
Their sum is the constant-in-Z cubic trace. Consequently the coefficient
of Z^-3 in w0 is zero. This trace argument is valid over C(u), without
any coefficientwise regularity assumption in C[u,w].

The remainder in (14) has leading term (A/M)*w^2. Relative to w0,
it changes w first in degree at most Z^-3: the derivative of H3^2
starts with 6*(h/M)*a^2*Z^5, whereas the quadratic remainder starts
with (A/M)*a^2*Z^2. Its correction in degree Z^-3 is therefore
-A/(6h). Linear and constant remainder terms begin later. Higher
Taylor terms cannot contribute in degree Z^2, since the correction
has order at least Z^-3. This proves (15), including A=0.

Differentiating the entire inverse with respect to P=Z^6 gives the
coefficient identity, for every integer j,

    partial_u q_(j-6)=-j*a_j/(6M).                       (16)

For j=-3, equations (15),(16) give

    partial_u q_-9=-A/(12hM).                            (17)

The right side is rational in u. Normalized trace from K to C(u)
therefore shows that the rational differential

    A*du/(hM)                                           (18)

has a rational primitive. Equivalently all its residues vanish.
Indeed partial fractions of a rational function have a rational
primitive exactly when every simple-pole coefficient is zero; the
polynomial part and every higher-pole term integrate rationally.
The residue at infinity then vanishes as well.

In particular, at every root c of hM with n=ord_c(hM)>0,

    A!=0 implies ord_c A != n-1.                        (19)

That order would make (18) have a simple pole with nonzero residue.
At an actual off-M A pole of slope L, ord h=2L, so
ord_c A=2L-1 is excluded. This supplements the half-order bound
ord_c A>=L from PF174; orders below 2L-1 can have higher poles and
must satisfy the full residue equation, rather than a leading-order
rule alone. If A=0, condition (18) is vacuous and does not close the
linear-remainder branch.

The rational residue condition is proper even after the leading ODE
and the half-order bound. Keep the D=24 positive leading data from
Section 5 and choose

    A_bad=u^3*(u-1)^9=hM/(u-1).

Then A_bad has order nine at the slope-five A root, so it passes the
half-order bound, but A_bad/(hM)=1/(u-1) has residue one. These
quadratic remainder coefficients cannot occur with an actual companion.
In contrast

    A_good=-u^3*(u-1)^8=-hM/(u-1)^2

passes both conditions: A_good/(hM) is the derivative of 1/(u-1).
Neither coefficient control is asserted to realize the other equations
or an actual Keller map.

## 8. Scope and verification

The exact checker verifies the leading formal coefficient, the primitive
identity, the polynomial ODE, local derivative constraints, degree and
symmetry obstructions, an exact inconsistency certificate for the
reflection example, the cube obstruction, the complete degree-27 integer
ledger, two-root polynomial primitives, the positive controls, the
full generic inverse coefficient (15), and both rational-residue controls.
The formal implicit construction, pole argument and degree/projection
proofs establish the unrestricted necessary condition. Finite controls
do not classify all solutions of (1) or of the residue equations.

Replay:

    .venv/bin/python research/jacobian/affine-complements/verify_ruling_six_leading_primitive.py

No target membership or target integrality is inferred from (1).
[PF188](RULING_SIX_SQUARE_DEFORMATION_AND_MARKED_QUADRATIC_EXCLUSION.md) excludes the retained marked quadratic ruling-six class.
JC2, all-prime derivative integrality and normal-tail termination remain
open. Independent whole-proof replay,
formal verification and external human review are not claimed.
