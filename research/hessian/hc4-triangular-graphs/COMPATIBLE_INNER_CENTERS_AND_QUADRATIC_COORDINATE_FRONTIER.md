# Inner critical points reduce to the four fixed source-orbit tests

15 November 2026. Claims **HC4MYJN0/1**. Write, over C,

    h=q_A+H(g),     g=q_B+2v^T A x+c,
    A=A^T invertible, B=B^T, T=A^(-1)B, deg H=d>=2.

For every affine inner quadratic g with a critical point,
the universal four-variable source-orbit inverse question
reduces to the same four fixed tests already retained:
the coordinate and quadratic orbits, h0 and h1.
The same holds whenever an admissible isotropic zero
direction exists, even if g has no critical point.

Beyond those tests, an additional affine-image mechanism
must have g a degree-two polynomial coordinate, singular
rank-two or rank-three T, and no admissible null direction.
Two explicit global coordinate normal forms describe
that remaining class. They are not asserted to be closed.
Neither the four tests nor unrestricted HC4/JC2 are solved.

## 1. Statements and exact retained inputs

**JN0.** Let 1<=k<=3, let Q be a nondegenerate quadratic
in 4-k additional variables, and let

    f=q_C+H(q_D+2a^T C z+c),  C=C^T invertible,
    D=D^T,  C^(-1)D invertible.

Every polynomial source reparameterization of Q+f with
unit Hessian has polynomially invertible gradient.
All offsets, nonlinear outer degrees and source degrees
are admitted, including arbitrary mixing of the quadratic
and composite variables.

This is a theorem for these specific stabilized families.
It is proved below from actual critical mass or whole-fibre
profiles, not from bare inverse closure or nonexistence
for an unstabilized function.

**JN1.** In four variables suppose at least one condition holds:

1. v belongs to im T, equivalently Crit(g) is nonempty.
2. There is a nonzero u with Tu=0, A(u,u)=0, A(v,u)=0.

The universal inverse question for this family is exactly
the conjunction of the same four fixed source-orbit tests.
The first-condition family alone already contains all four.
No new fixed test is created by compatible inner centers
or admissible isotropic zero directions.

After the prior closures, any additional affine mechanism
beyond the four tests must satisfy

    det T=0, rank T in {2,3}, v not in im T,
    ker T intersect v^perp contains no nonzero A-isotropic vector.
                                                               (1)

Its inner g is a polynomial coordinate of degree two.
An explicit invertible polynomial coordinate map gives
one of the following source normal forms:

    h=H(t)+(x/a)t+q(z)-(x/a)p(x,z),                     (2)

with two z variables, p of degree<=2, a!=0; or

    h=H(t)+q(z)+(t-p(z))^2/(4a^2),                    (3)

with three z variables and p of degree<=2.
Thus the part outside H has total degree at most three
or four respectively. No degree is bounded for H or
an eventual full source reparameterization.

Retained foundations:

- [HG0/1/2](GLOBAL_VALUE_ORBITS.md): native critical schemes,
  their constant critical-Hessian factor, whole generic
  profile L^3-NL under hypothetical gradient nonautomorphy,
  split one-variable exclusion, and the two old endpoints.
- [HC0](HODGE_TRACE_CLOSURE.md): additive/product Hodge
  invariant applied to explicit strata.
- [JH0](SIMPLE_STATIONARY_POLES_AND_NONCYCLIC_KRYLOV_CLOSURE.md):
  complete reciprocal-Hessian mass with compatible simple
  stationary groups and the literal infinity condition.
- [JI0](AFFINE_POLAR_CURVES_AND_COMPLETE_INVERTIBLE_SEMISIMPLE_CLOSURE.md):
  exact all-offset distinct and scalar plane profiles.
- [JL0/1](JORDAN_POLAR_CURVES_AND_ALL_INVERTIBLE_PENCIL_CLOSURE.md):
  exact plane Jordan genus/punctures, whole eigenvector
  fibre identity, and complete invertible four-variable closure.
- [JM0/1](ZERO_EIGENVECTOR_SHEARS_AND_RANK_ONE_PENCIL_CLASSIFICATION.md):
  the global null shear, stabilized actual critical mass,
  and complete rank-at-most-one four-test reduction.
- [JD0](NILPOTENT_QUADRATIC_PENCILS_AND_ENDPOINT_ORBITS.md)
  supplies the proved constant self-adjoint primary frames.
  [HP0](QUADRATIC_ISOMETRY_SATURATION.md) and
  [HM0](AMBIENT_QUADRIC_GRAPHS.md) supply literal ambient transfers.

Use L=uv, U=u+v, epsilon=d mod 2,

    a_d=(d+epsilon-2)/2, rho=2-epsilon,
    M=L-a_d U+1-rho.

M is the exact generic double-curve profile from JL0.
The hyperelliptic curve w^2=kappa-P(z), deg P=2d,
has d-1 compact genus and two infinity punctures:
the generic polynomial has 2d simple roots, so the
retained characteristic-zero Riemann--Hurwitz formula
gives

    B_d=L-(d-1)U-1.                                   (4)

In particular d-1>a_d for every d>=2.
We do not substitute a zero eigenvalue into JI0's
distinct-nonzero-eigenvalue theorem.

## 2. Plane composites with a quadratic hyperbolic stabilization

Here k=2; normalize Q=2xy. For every plane f,
the x!=0 stratum solves y, while x=0 retains the
whole plane curve and free y. Thus

    E(2xy+f=kappa)=L^3-L^2+L E(f=kappa).                (5)

If the invertible plane pencil has two distinct eigenvalues,
JI0 gives strictly positive total curve genus in every
affine offset. Its U term in (5) forbids HG0's profile.

For a scalar nonzero plane pencil, the exact JI0 profiles
are d(L-1) for zero offset, B_d for nonisotropic offset,
and L-d for nonzero isotropic offset. The first gives
L^3+(d-1)L^2-dL, with a nonzero L^2 coefficient.
The second has positive genus. The last gives L^3-dL.
For that last case the complete native plane critical set
has at most one point: its stationary inner function is
constant, the regular equation is linear, and the only
determinant parameter is inconsistent because of the
nonzero isotropic scalar offset. The quadratic
stabilization retains that full critical count.
Since d>=2, HG0 cannot have coefficient d.

For a plane Jordan pencil, translate its invertible
inner quadratic to center. JL0's marked plane has
positive total genus, again sufficient in (5).
The unmarked plane has exact profile M-d.
For d>=3 the genus a_d is positive. For d=2,
(5) is L^3-3L; its full Morse critical count is at
most two, as follows explicitly.

In the JL0 frame

    f=2zw+alpha z+H(2lambda zw+z^2),

the regular stationary point lies on the isotropic
eigenline, so its inner value is zero and the equation
is s-H'(0)=0. There is at most one regular point.
At s=-1/lambda stationarity fixes z=lambda alpha/2
and leaves w free. Its inner value is affine-linear
in w if alpha!=0, and constant if alpha=0.
For d=2 there is at most one extra point or a
forbidden whole critical line. Thus N<=2, not three.
No singular stationary group is dropped.

This proves inverse closure for k=2. The argument uses
actual profiles and counts; it does not claim every
unstabilized nonexistence statement is stable.

## 3. Ternary composites with a quadratic-line stabilization

Let k=3 and normalize Q=w^2. The inner pencil on
the other three variables is invertible. Translate it
to center, retaining

    h=w^2+q_3+ell_3+H(q_B3).                           (6)

The original affine-inner form is also available by the
inverse translation. Use that form for the semisimple
Witt cases and (6) for the unmarked-chain cuts below.

### Complete critical-mass cases

For a semisimple ternary pencil, assume each nonzero
projected affine offset is nonisotropic on its eigenspace.
Its marked lines form a nondegenerate orthogonal Krylov
space. At an active eigenvalue R has a double pole
and the exact finite budget is 2(d-1)>=k_lambda.
An inactive repeated eigenspace has a positive-dimensional
stationary critical quadric, already forbidding unit Hessian.
An inactive simple root supplies its whole JH0 group.
Thus the finite mass conditions hold except possibly
for d=2 and a scalar primary of multiplicity three.
Infinity has R bounded, P of degree three and F~s,
so W=O(s^(-4)). The complete mass then stabilizes
by JM0 and excludes unit Hessian in four variables.

For a ternary length-three Jordan block with marked
head, the exact active R pole order is four and its
determinant multiplicity three. For a marked length-two
block plus one line, its active pole order is three
and its multiplicity at most three. Their top cyclic
pairings are nonzero in the JD/JH frames. A distinct
remaining line is active or a compatible simple root.
The same finite and infinity conditions hold, so the
actual mass excludes unit Hessian after stabilization.

### Semisimple isotropic and large-primary cases

Suppose a nonzero projected offset is isotropic. It
has eigenspace dimension two or three. A full-eigenspace
Witt pair gives the literal JI whole-fibre strata.
With an additional quadratic zero-inner variable w,
the residual plane has the explicit form

    f0=z^2+w^2+linear(z)+H(mu z^2+linear(z)+constant),
    mu!=0.

Its whole generic curve is w^2=kappa-P(z), deg P=2d,
and has profile B_d from (4).

For an isotropic eigenspace of dimension two, the
subtracted base has rank two, profile L^2+(d-1)L;
the whole profile is

    E=L^3-L^2-(d-1)L+L B_d.

It has a nonzero genus coefficient d-1.
For a scalar dimension-three isotropic primary, the
subtracted base has rank one, profile L M; the whole
profile is L^3-L M+L B_d. Its genus coefficient is
d-1-a_d>0. All closed Witt strata are retained.

The remaining scalar dimension-three nonisotropic offset
also has a global fibre model. Choose

    q_A=2uv+x^2+w^2,
    g=2lambda uv+lambda x^2+2a x+c, a!=0.

On h=kappa,t=g solve polynomially

    x=[t+lambda H(t)-lambda kappa+lambda w^2-c]/(2a),
    2uv=kappa-H(t)-x(t,w)^2-w^2.

This is a whole-fibre hyperbolic isomorphism, with its
displayed polynomial inverse. Its zero curve is
x^2+w^2+H(lambda x^2+2a x+c)=kappa,
again of profile B_d. The full profile is
L^3-L^2+L B_d, with positive genus d-1.
This includes the quadratic outer degree where the
finite active multiplicity budget can fail.

### Nonsemisimple unmarked heads

The only ternary nonsemisimple partitions are 3 and 2+1.
Their marked heads were covered by complete critical mass.
For the unmarked length-three head, use

    q_A=2xz+y^2+w^2,
    q_B=lambda(2xz+y^2)+2xy,
    ell=A x+B y.

The JL whole-fibre cut along partial_z has ternary
quadric 2xy-lambda w^2-lambda(Ax+By), of rank three.
The cut base is forced if B!=0 and rank one if B=0.
Thus E(Y)=L^3-L^2+L M and E(Y0)=L^2 or L M.
The remaining plane is

    y^2+w^2+B y+H(lambda y^2),

of profile B_d. Its genus d-1 exceeds any subtracted
a_d, while the Y genus only adds a nonnegative term.
The full HG0 profile is impossible.

For an unmarked length-two block plus one line use

    q_A=2xy+z^2+w^2,
    q_B=2lambda xy+x^2+mu z^2,
    ell=A x+B z, lambda*mu!=0.

The cut along partial_y has ternary quadric
x^2+(mu-lambda)z^2-lambda w^2-lambda(Ax+Bz).
If mu!=lambda its rank is three and the cut rank two.
If mu=lambda and B!=0, both are forced.
If mu=lambda and B=0, their profiles are
L^3+(d-1)L^2 and L M.
In all cases the plane
z^2+w^2+Bz+H(mu z^2) has profile B_d.
Again d-1>a_d prevents every possible genus cancellation.

These alternatives exhaust every invertible ternary
pencil and every affine offset. All cut-stratum
identities are the literal JL/JI fibre maps, with the
new zero-variable curve computed by (4), not by
extending a nonzero-eigenvalue curve assertion.

For k=1, the composite is a one-variable polynomial of
degree 2d>=4, and Q has rank three. HG1's critical
second-derivative argument excludes unit Hessian.
This completes JN0 for k=1,2,3.

## 4. Every admissible null kernel reduces to the same fixed tests

Suppose u satisfies condition 2 of JN1. The actual
JM0 shear gives

    h ~ 2xy+q(z1,z2)+H(g0(z1,z2)).

The quadratic pairing q is nondegenerate. Its reduced
inner pencil S is the self-adjoint operator induced by T
on u^perp/ C u: u^perp is T-stable by self-adjointness,
and T u=0 makes the quotient well-defined. Setting x=0
and suppressing y in the explicit Witt frame gives
exactly those quadratic and affine coefficients.

If rank S=2, the reduced plane pencil is invertible
and JN0 applies to the full four-variable suspension.
If rank S<=1, the new full pencil is 0 on the hyperbolic
pair and S on the plane, so has rank at most one.
JM1 then gives precisely the old/fixed four-test
universal reduction. No bare plane conclusion is
assumed to stabilize. This covers every admissible
null kernel without a Jordan or degree enumeration.

## 5. Every inner critical point is covered

The inner stationary equation is T x=-v, so Crit(g)
is nonempty exactly when v belongs to im T.
Self-adjointness gives

    im T=(ker T)^perp.                                 (7)

If T is invertible, JL1 already applies. Otherwise,
if ker T contains a nonzero isotropic u, (7) makes
A(v,u)=0 automatically. Section 4 applies.

If ker T contains no such vector, over C it must be
a one-dimensional nondegenerate space: every complex
quadratic form in at least two variables has a nonzero
isotropic vector, including degenerate forms. A
one-dimensional degenerate kernel is itself isotropic.
Thus

    C^4=ker T orthogonal-direct-sum (ker T)^perp,

with both summands nondegenerate. The restriction of T
to its three-dimensional complement is invertible,
and v lies in that complement by (7). The potential
is a nondegenerate quadratic line plus an invertible
ternary affine composite. JN0 applies directly.

Consequently every compatible inner critical center
reduces to the same four universal tests.

### All four tests already occur in this compatible-center family

For the coordinate test take

    q_A=2xz+y^2+w^2,  g=2xy+2y,
    H(t)=-t^2/4+t/2.

Here T has a zero Jordan block of length three and a
zero line, v=(0,1,0,0) belongs to im T, and x=-1,y=0
is an inner critical point. The null shear along
partial_z gives exactly 2xZ+w^2+y, a coordinate.
The old quadratic test q_A=2xy+z^2+w^2,g=x^2,H=t^2
has an inner critical point and the explicit quadratic
source orbit already retained in JM.

For h0,h1 use

    q_A=2xz+y^2+w^2,
    g_a=2xy+2y+w^2+2a w+a^2,
    H(t)=-t^2/4+t/2.

Its pencil has a zero length-three block and a
nonzero line. The affine inner has critical points
x=-1,y=0,w=-a, with z free. The same global null
shear gives

    2xZ+y^2+w^2+H(2y+w^2+2a w+a^2).

This is JM's resonant rank-one example, with the
plane roles interchanged. The explicit JM inverse
normalization gives h0 at a=0 and h1 at a=-1/2,
up to constants and variable names.

Every example is the actual nondegenerate ambient
graph Q=q_A+a0 X^2+b0 X+c0, X=g, with a0!=0.
Projection to the retained source is its polynomial
parameter inverse. Thus necessity as well as
sufficiency of the four tests holds even in the
first-condition family alone. This does not identify
each individual compatible potential with a fixed
source orbit; the other alternatives are inverse-closed.

## 6. The remaining affine inner must be a quadratic coordinate

Outside the four fixed tests and the inverse-closed
branches above, v is not in im T and no admissible
u exists. By (7), choose u in ker T with
a=A(v,u)!=0. Its directional derivative is

    D_u g=2A(Tx+v,u)=2a.

Hence in coordinates with u=partial_y,

    g=2a y+p(other variables),
    t=g,   y=(t-p)/(2a)                                (8)

is an actual global polynomial coordinate map and inverse.
If u is isotropic, a constant Witt frame gives
q_A=2xy+q(z1,z2), and (8) gives exactly (2).
If u is nonisotropic, normalize A(u,u)=1, split
q_A=y^2+q(z1,z2,z3), and obtain exactly (3).

If ker T contains an isotropic vector, the absence
of an admissible null vector makes its affine pairing
nonzero, so the cubic form (2) is available.
Otherwise the kernel is the nondegenerate line from
Section 5 and the quartic form (3) is available.
The affine changes and (8) have global polynomial
inverses; the full source-orbit problem is retained.

The prior JL/JM closures restrict any additional
affine mechanism to rank two or three and det T=0.
Here g has a nonzero quadratic part, so its degree
is exactly two. Further, ker T intersect v^perp
has no nonzero isotropic vector. Its complex
dimension is at most one; since the affine pairing
on ker T is nonzero, this is the corresponding
zero-primary constraint in (1).

This is a necessary frontier and a normal-form
reduction, not unit-Hessian nonexistence or gradient
inverse closure for the remaining coordinate-inner
cubic/quartic forms. Native Hessians are not invariant
under (8); only the entire source orbit is identified.

## 7. Scope and source audit

All source, graph and outer degrees are unrestricted.
The two fixed h0,h1 full source orbits remain open;
JK excludes their base-preserving maps and retains
their entire passing special-value Hodge profiles.
Neither old actual endpoint, a universal two-point
Keller equivalence, the remaining quadratic-coordinate
frontier nor unrestricted HC4/JC2 is solved.

The primary
[Riemann--Hurwitz foundation](https://stacks.math.columbia.edu/tag/0C1B)
and retained [Soprunov record](https://arxiv.org/abs/math/0511684v2)
were rechecked. No new unproved external theorem,
stable cancellation, arbitrary finite-cover E
multiplicativity or generic source inverse is imported.
