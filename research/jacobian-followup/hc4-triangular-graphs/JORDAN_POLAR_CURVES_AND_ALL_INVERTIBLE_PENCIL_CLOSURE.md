# Every invertible affine quadratic pencil is inverse-closed

15 November 2026. Claims **HC4MYJL0/1**. In four variables let

    h=q_A+H(q_B+2v^T A x+c),   T=A^(-1)B,
    A=A^T invertible, B=B^T,  deg H=d>=2.

**JL1.** If T is invertible, every polynomial source
reparameterization of h having nonzero constant Hessian determinant
has a polynomial gradient inverse. Diagonalizability, cyclicity
of the offset and nondegeneracy of its Krylov space are not
required. The affine offset, outer degree and source degree
are unrestricted.

Together with JI, this closes every nonsemisimple nonzero-eigenvalue
branch. It is inverse closure, not nonexistence of every
unit-Hessian source orbit. Singular T, the two fixed quadratic
multiplier models h0,h1, both old plane endpoints and unrestricted
HC4/JC2 remain open.

The actual graphs X=P0(q_B+2v^T A x+c) in
Q=q_A+a0 X^2+b0 X+c0, a0!=0, P0 nonconstant, have
outer H=a0 P0^2+b0 P0+c0 of degree>=2. The result therefore
applies in every admitted graph/source degree, with the retained
certified quadratic-similitude and literal Meng--Yang transfers.

## 1. Retained foundations and normalization

Use the exact whole-fibre obstruction [HG0](GLOBAL_VALUE_ORBITS.md):
under unit Hessian and gradient nonautomorphy, generic E=L^3-NL,
where N is the full critical count. A positive-dimensional or
degenerate native critical scheme already forbids unit Hessian.
The additive/product Hodge invariant and affine-curve formula
are explicitly retained in [HC0](HODGE_TRACE_CLOSURE.md).
Characteristic-zero Riemann--Hurwitz is retained in
[IY0](CONTRACTIBLE_RADICAL_CURVES_AND_PENCIL_DEFECT.md).

The constant self-adjoint primary frames are proved in
[JD0](NILPOTENT_QUADRATIC_PENCILS_AND_ENDPOINT_ORBITS.md);
adding lambda q_A gives the corresponding nonzero-eigenvalue
blocks. The cyclic-offset exclusion is
[JF1](AFFINE_QUADRATIC_PENCIL_RESIDUES_AND_ENDPOINTS.md),
and the nondegenerate proper-Krylov budget is
[JH1](SIMPLE_STATIONARY_POLES_AND_NONCYCLIC_KRYLOV_CLOSURE.md).
All semisimple cases and the two-variable distinct/scalar
curve profiles are [JI0/1](AFFINE_POLAR_CURVES_AND_COMPLETE_INVERTIBLE_SEMISIMPLE_CLOSURE.md).
The remaining singular two-point tests are
[JK0/1](SUSPENDED_MULTIPLIER_ORBITS_AND_BASE_PRESERVING_CLOSURE.md).
Ambient transfers are [HP0](QUADRATIC_ISOMETRY_SATURATION.md)
and [HM0](AMBIENT_QUADRIC_GRAPHS.md).

Because T is invertible, translate the inner quadratic to its
center. Explicitly u=T^(-1)v and x=z-u give
q_B+2v^T A x+c=q_B(z)+c-v^T A T^(-1)v.
Absorb this constant in H. Up to a constant in the potential
we may work with

    h=q_A+ell+H(q_B),       ell any linear form.           (1)

This is a literal source translation with polynomial inverse,
not nonlinear Hessian invariance.

For the rest of the proof write

    epsilon=d mod 2,  a=(d+epsilon-2)/2,  rho=2-epsilon,
    M=L-a U+1-rho,   U=u+v.

Here M is the E-polynomial of a generic smooth double curve
Y^2=t+lambda H(t)-lambda kappa+constant, lambda!=0.
It has genus a and rho actual infinity points. Every occurrence
has this same profile even if its constant differs. Generic
kappa avoids only finitely many values depending on the fixed
polynomials and coefficients.

## 2. Plane Jordan curves (JL0)

Let lambda!=0 and

    f(x,y)=2xy+alpha x+beta y+H(2lambda xy+x^2),
    t=2lambda xy+x^2,  P(t)=t+lambda H(t)-lambda kappa.

### A cyclic marked direction: beta!=0

On f=kappa one solves polynomially

    y=[x^2-lambda alpha x-P(t)]/(lambda beta).

The resulting exact affine curve has equation

    2x^3+(beta-2lambda alpha)x^2
       -(2x+beta)t-2lambda x(H(t)-kappa)=0.               (2)

It is monic of degree three in x. Thus t is finite free
of total degree three, with basis 1,x,x^2.
Its polar is the smooth parabola

    J(t,f)=4x^2+2(beta-lambda alpha)x+2lambda beta y=0.

On it

    t=-4x^3/beta+(2lambda alpha/beta-1)x^2.

Therefore f restricted to the polar has exact degree 3d:
H(t) dominates its degree-at-most-three remaining terms.
A generic kappa meets this parabola transversely in 3d
distinct points. Those are all the simple affine
ramification points of t on the smooth generic curve.

At t=infinity, (2) has one small branch
x~[-beta/(2lambda h_d)]t^(1-d), unramified over t.
The two large leading roots are
x~+-sqrt(lambda h_d)*t^(d/2). For even d these give
two unramified points; for odd d they pair into one
point of index two. The scaled initial roots are simple:
the small-root initial equation is linear, and at either
large nonzero root the derivative is 4lambda h_d!=0.
The finite degree three accounts for all branches.

Consequently the actual puncture count is 3-epsilon
and infinity ramification is epsilon. If the curve has
c_C components and total compact genus G_C, componentwise
Riemann--Hurwitz gives

    2G_C-2c_C=-6+3d+epsilon,
    G_C=(3d-6+epsilon)/2+c_C
        >=(3d-4+epsilon)/2 > a >=0.                     (3)

No connectedness theorem is assumed. In particular the
total genus is positive in every d>=2, including d=2.

### An unmarked head: beta=0

The exact double-curve presentation is

    x^2-lambda alpha x=P(t),    y=(t-x^2)/(2lambda x).

Its compact genus is a and its infinity punctures rho.
For generic kappa the d points with x=0 have distinct
t values, none equal to zero. They are not actual plane
points, since the original inner t would be zero at x=0.
All other points are recovered by the displayed inverse.
Thus the whole generic plane profile is

    E(f=kappa)=M-d.                                    (4)

If alpha=0 these d omitted points are branch points,
one point per root, not two. This does not change (4).
These two statements comprise JL0 for every affine offset.

## 3. One whole-fibre identity and its rank profiles

Suppose a nonzero isotropic eigenvector with eigenvalue lambda
annihilates ell. A constant Witt frame writes

    q_A=2xy+q(z1,z2),   q_B=2lambda xy+g0(x,z1,z2),
    ell=A x+ell0(z1,z2).

On h=kappa put t=q_B. Define Y and Y0 by

    Y:  g0-lambda q-lambda ell=P(t),
    Y0=Y intersect {x=0},
    f0=q+ell0+H(g0(0,z)).

On x!=0, Y recovers y=(t-g0)/(2lambda x);
the equation of Y then makes h=kappa. This is a
regular isomorphism with regular inverse on that open set.
At x=0 retain the whole level f0=kappa and the free y line.
Hence

    E(h=kappa)=E(Y)-E(Y0)+L E(f0=kappa).                 (5)

Y is an affine quadric in three variables equal to P(t).
After constant completion of squares, its exact profiles are:

| Quadratic rank and remaining affine force | E(Y) |
|---|---|
| Nonzero linear force on the quadratic kernel | L^3 |
| Rank three | L^3-L^2+L M |
| Rank two, no force on its kernel | L^3+(d-1)L^2 |
| Rank one, no force on its kernel | L^2 M |

For Y0, a quadric in two variables, the corresponding profiles
are L^2 for a kernel force, L^2+(d-1)L for rank two,
L M for rank one without force, and d L^2 for rank zero
without force. These formulas are literal products and
hyperbolic strata: uv=P(t) contributes
(L-1)L+dL=L^2+(d-1)L, and
2uv+w^2=P(t) contributes (L-1)L^2+L M.
No finite-cover E multiplicativity is used.

For the plane f0 use JL0, JI0, or the scalar profiles:
zero offset d(L-1), nonisotropic offset
L-(d-1)U-1, nonzero isotropic offset L-d.
For a distinct nonzero semisimple plane its total genus
is at least 2d-3+epsilon>a.
For a cyclic plane Jordan curve (3) its genus also exceeds a.
Thus a subtracted L M cannot cancel either of these larger
genera. All curve genus coefficients are nonnegative.

The remaining sections exhaust the fixed partitions of four.
They are spectral rank cases, independent of outer/source
degree, not a growing polynomial-degree calculation.

## 4. One block of length four

Use

    q_A=2xw+2yz,
    q_B=lambda q_A+2xz+y^2.

If ell has nonzero w coefficient its marked vector is cyclic,
so JF1 already excludes unit Hessian. Otherwise write
ell=A x+B y+C z and use (5) with the eigenvector partial_w.
The ternary quadric is

    2xz+y^2-lambda(Ax+By+Cz).

It has rank three, so E(Y)=L^3-L^2+L M.
If C!=0, Y0 has a kernel force and E(Y0)=L^2.
The plane f0=2yz+By+Cz+H(2lambda yz+y^2)
is cyclic Jordan, of positive genus by (3). Both M and
f0 enter with positive signs in (5); their negative genus
terms cannot cancel.

If C=0, E(Y0)=L M and f0 is the unmarked Jordan curve
of (4). Formula (5) is

    E=L^3-L^2+L(M-d).

For d>=3 its genus a is positive. For d=2 it is L^3-3L.
Under a hypothetical Morse critical set that exceptional
case has at most one point, as proved in Section 8.
Neither profile can be HG0's L^3-NL.

## 5. A block of length three and one line

Use

    q_A=2xz+y^2+w^2,
    q_B=lambda(2xz+y^2)+2xy+mu w^2,
    lambda*mu!=0.

A nonzero z coefficient of ell gives a cyclic marked
length-three primary. Its cyclic restricted pairing is
nondegenerate by the retained top-moment argument.
JH1 covers it: its active multiplicity is at most four
and m(3+1)>=4 for m=d-1>=1; any different unmarked
one-dimensional primary is simple. Thus assume
ell=A x+B y+C w and cut along partial_z.

The ternary quadric is

    2xy+(mu-lambda)w^2-lambda(Ax+By+Cw).

If mu!=lambda, Y has rank three, while Y0 is forced
or rank one. The plane f0 is a distinct nonzero
semisimple pencil, with genus greater than a. Its genus
cannot be canceled by Y0; the contribution of Y only
adds a nonnegative genus. This excludes HG0.

Let mu=lambda. If C!=0, both Y and Y0 are forced,
so E=L^3-L^2+L E(f0). The plane is scalar. A nonisotropic
offset has genus d-1>0. For the nonzero isotropic offset
the profile is L^3-dL and the Morse critical count is
at most one (Section 8), also incompatible with HG0.
If C=0,B!=0, Y has rank two, Y0 is forced and the scalar
plane offset is nonisotropic, again giving positive genus.
If B=C=0, direct substitution in (5) gives

    E=L^3+(d-1)L^2-dL,

whose L^2 coefficient is nonzero for every d>=2.

## 6. One length-two block and two semisimple lines

Use

    q_A=2xy+z^2+w^2,
    q_B=lambda(2xy)+x^2+mu z^2+nu w^2.

First suppose the y coefficient of ell is zero:
ell=A x+B z+C w. Cut along partial_y. The ternary
quadric is

    x^2+(mu-lambda)z^2+(nu-lambda)w^2
                                    -lambda(Ax+Bz+Cw).

Let r count the remaining eigenvalues unequal to lambda.

- r=2: Y has rank three and Y0 rank two. A distinct
  semisimple plane or a nonisotropic scalar plane has
  positive genus, adding to that of M. For a zero scalar
  offset the L^2 coefficient in (5) is d-1, hence nonzero.
  For a nonzero isotropic scalar offset, a>0 excludes
  d>=3; at d=2 the exact profile is L^3-4L with at most
  two Morse critical points (Section 8).
- r=1: if ell has nonzero coefficient on the remaining
  lambda line, Y,Y0 are forced and f0 has positive genus.
  Otherwise Y has rank two, Y0 profile L M, and the
  distinct plane genus exceeds a, so cancellation is
  impossible.
- r=0: a nonzero remaining offset forces both Y,Y0.
  Its nonisotropic scalar plane has positive genus.
  Its isotropic scalar plane gives L^3-dL with at most
  one Morse critical point. If the remaining offset is
  zero, the exact profile is L^2 M-dL: either a>0 or,
  at d=2, its L^2 coefficient is -1.

Now let the y coefficient B be nonzero. The marked
length-two primary has nondegenerate cyclic pairing.
JH1 covers all cases except these:

1. The other primary is a distinct scalar eigenvalue of
   multiplicity two with a nonzero isotropic offset.
   A Witt pair in that scalar primary gives the same
   whole-fibre strata as JI: E=L^3-L^2-(d-1)L+L E(f0).
   Here f0 is the cyclic Jordan plane, with positive genus.
2. All eigenvalues equal lambda and d=2, where the active
   budget 3(d-1)>=4 can fail. The following whole-fibre
   model works in every d and closes it.

For case 2 write ell=A x+B y+C z+D w, B!=0 and put
K=(C^2+D^2)/B^2. On the entire fibre solve

    Y0(t,x)=[x^2-lambda A x-P(t)]/(lambda B),
    y=Y0-Cz/B-Dw/B,
    Z=z-Cx/B,  W=w-Dx/B,
    Z^2+W^2=S(t,x):=(t-x^2)/lambda-2xY0+Kx^2.           (6)

These formulas and their reverse give a global polynomial
fibre isomorphism. Thus E=L^3-L^2+L E(S=0).
The curve S=0 is the actual plane section Z=W=0.
Its quadratic forms are a constant shear of
2xy and 2lambda xy+x^2; its y-linear coefficient is
the same B!=0. It is a cyclic Jordan plane, so (3)
gives positive genus and the HG0 profile is impossible.

## 7. Two blocks of length two

Use q_A=2xy+2zw and

    q_B=2lambda xy+x^2+2mu zw+z^2.

If lambda!=mu and both head coefficients ell_y,ell_w
are nonzero, the marked vector is cyclic and JF1 applies.
Otherwise choose a block whose head coefficient vanishes
and cut along its eigenvector. Y has rank three and Y0
rank two. If the other plane is cyclic Jordan, its genus
is positive. If both heads vanish, use (4) to obtain

    E=L^3-2L^2+2L M-(2d-1)L.

For d>=3 its genus coefficient 2a is positive. For d=2
this is L^3-5L, with at most three Morse critical points
by Section 8.

Let lambda=mu. The head coefficient vector is
(ell_y,ell_w); its squared norm is ell_y^2+ell_w^2.

- If the head vector itself is zero, both heads are unmarked. Y has rank
  two, Y0 profile L M and f0 profile M-d. Therefore
  E=L^3+(d-1)L^2-dL, with nonzero L^2 coefficient.
- If the vector is nonisotropic, a constant orthogonal
  change in the two block multiplicity coordinates makes
  one head zero and the other nonzero. Cutting the
  unmarked block forces Y,Y0 via the other dual coefficient.
  The remaining cyclic Jordan plane has positive genus.
- If the vector is nonzero isotropic, interchange the
  blocks if needed and write it (beta,i beta), beta!=0.
  Put r=x+iz, s=x-iz,
  p=(y-iw)/2, q=(y+iw)/2. Then

      q_A=2rp+2sq,   q_B=lambda q_A+rs,
      ell=a0 r+b0 s+2beta q.

  Cutting along partial_p forces Y,Y0 because of q.
  The plane f0 is scalar with offset b0 s+2beta q.
  If b0!=0 it is nonisotropic, of genus d-1.
  If b0=0 it is isotropic and gives L^3-dL with at
  most one Morse critical point. Both are excluded
  from the nonautomorphic HG0 profile.

The first bullet means that the vector itself is zero,
not merely that its squared norm is zero; the third
bullet retains exactly that nonzero isotropic case.

## 8. The low-genus exceptions keep their full critical counts

The few profiles without a genus or L^2 obstruction above
must not be disposed of by generic point counting.
Here are their complete bounds, under the Morse assumption
necessary for unit Hessian.

For (1), put u=(1/2)A^(-1)ell. At a regular stationary
parameter s=H'(q_B), the unique point is

    x(s)=-(I+sT)^(-1)u.

In the low-genus cases every moment
u^T A T^j u is zero. Expanding the resolvent and using
Cayley--Hamilton gives q_B(x(s))=0 identically.
Thus the regular equation is s-H'(0)=0 and there is
at most one regular critical point. The moment vanishing
is immediate in the displayed constant frames:

- Length four with ell=A x+B y has u in its bottom
  two chain positions, all pairings with its T powers zero.
- The scalar-isotropic cases in Sections 5,6,7 have zero
  squared offset norm and zero higher nilpotent moments.
  Their only determinant root is blocked by a nonzero
  linear coefficient on the stationary kernel. They
  have at most one critical point in total.
- In Section 6,r=2's scalar-isotropic case, u is the sum
  of isotropic eigenvectors in two distinct primary
  spaces, so all moments vanish.
- With two distinct unmarked length-two blocks, u is
  a sum of their isotropic eigenvectors, again with
  every moment zero.

At a singular parameter the following elementary
stationary calculations supply all omitted source points.

For length four with ell=A x+B y, s=-1/lambda forces
x=0,y=lambda B/2,z=lambda A/2, with w free. The
inner q_B is constant along that line. Hence either
no point is critical, or the whole line is. A Morse
critical set contains no such group: N<=1 in the
L^3-3L exception.

For an unmarked length-two block q_A=2xy,
q_B=2lambda xy+x^2, ell=A x, stationarity at
s=-1/lambda fixes x=lambda A/2 and leaves y free.
When other primary eigenvalues differ, their coordinates
are uniquely fixed. In the low-genus cases their inner
quadratic values are zero. The remaining inner value is
affine-linear in y, or constant if A=0. For d=2 the
equation H'(q_B)=-1/lambda therefore gives at most
one point, or a whole forbidden critical line.

Accordingly Section 6,r=2 has at most one such group
at lambda; its scalar-isotropic primary is blocked.
Together with the regular point N<=2, incompatible
with coefficient four. Two distinct unmarked
length-two blocks give at most one point at each
singular parameter plus the regular point, N<=3,
incompatible with coefficient five.

Every possible full critical point is regular or lies
at one of those determinant roots. None is discarded.
If any line or degeneracy occurs, HG0's native critical
scheme obstruction already applies. These bounds are
not a properness assertion or a multivariable residue theorem.

## 9. Exhaustion and exact remaining boundary

An invertible self-adjoint four-variable pencil is either
semisimple, already covered by JI1, or has one of the
fixed Jordan configurations in Sections 4--7:
4, 3+1, 2+1+1, or 2+2, with all eigenvalues nonzero.
Their marked-head alternatives are exhausted by JF/JH
or the whole-fibre arguments above. Thus every
hypothetical unit-Hessian member with nonautomorphic
gradient is excluded. This proves JL1.

No degree of H, source automorphism, graph polynomial
or affine coefficient has been bounded. The proof uses
one whole-fibre surgery, two plane Jordan profiles,
constant quadratic ranks and complete low-genus counts.

Invertibility remains essential to these arguments:
the inner translation, nonzero-eigenvalue plane curves
and generic double-cover degrees all use it. A remaining
affine quadratic-composite image must have singular T.
The singular h0,h1 full source orbits are not closed,
and the old coordinate and quadratic endpoints are
not solved. Arbitrary image geometry, full w-dependent
global recovery and unrestricted HC4/JC2 remain open.

We rechecked the retained primary
[Riemann--Hurwitz foundation](https://stacks.math.columbia.edu/tag/0C1B)
and the related [Arzhantsev--Zaitseva suspension record](https://arxiv.org/abs/2309.06170v3).
No homogeneous-space, cancellation, connectedness or
suspension-classification theorem is imported. Every
fibre isomorphism used here is explicit.
