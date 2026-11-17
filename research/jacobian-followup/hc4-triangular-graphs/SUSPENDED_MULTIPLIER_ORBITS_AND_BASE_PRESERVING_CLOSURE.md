# Suspended multipliers reduce to four source-orbit tests

15 November 2026. Claims **HC4MYJK0/1**. This treats

    h=2xy+w C(z)+D(z),     C,D in C[z],                   (1)

in arbitrary multiplier and remainder degrees. Under arbitrary
polynomial source automorphisms the universal gradient-inverse
question reduces to four fixed source-orbit tests: the old
coordinate and nondegenerate-quadratic orbits, and

    h_0=2xy+w(z^2-1),       h_1=2xy+w(z^2-1)+z.          (2)

The last two full source orbits remain open. They have no
unit-Hessian member under **any** source automorphism preserving
the base algebra C[x,z] setwise, including nonlinear base and
nonlinear dual automorphisms. An arbitrary source automorphism
may change that base algebra; no such alignment is assumed.

## 1. Statements and exact retained inputs

**JK0.** Every potential (1) outside the four source-orbit tests
is inverse-closed: any polynomial source reparameterization
having nonzero constant Hessian determinant has polynomially
invertible gradient. More precisely:

- A squarefree multiplier C of degree at least three has no
  unit-Hessian member in its full polynomial source orbit.
- A multiplier with repeated roots is either excluded by a
  degenerate critical locus, or is inverse-closed by its whole
  generic Hodge profile and actual critical count.
- Nonzero constant C gives the coordinate orbit; linear C
  gives the nondegenerate-quadratic orbit.
- Squarefree quadratic C reduces, by an explicit polynomial
  source automorphism and constant translation of h, to
  h_alpha=2xy+w(z^2-1)+alpha z. For alpha=0 this is h_0;
  for alpha!=0 its entire source orbit differs from h_1's
  only by a nonzero scalar factor. Both remain two-critical-point
  tests, not proved Keller endpoint equivalences.
- If C=0, a nonconstant linear D gives a coordinate; all other
  D give a positive-dimensional critical locus.

These four tests also suffice for, and each occurs in, the genuine
closed graph family

    Q=2xy+z^2+w^2-X^2,       X=w+A(z),                  (3)

with arbitrary A in C[z]. The ambient quadratic form is
nondegenerate. A graph's actual inverse parameter map is its
projection to x,y,z,w. The resulting potential has
C=-2A and D=z^2-A^2, and every parameterization of this fixed
graph factors by a source automorphism. Certified quadratic
similitudes and the literal original Meng--Yang transfer apply.

**JK1.** No source automorphism preserving C[x,z] setwise
makes (1) have unit Hessian if deg C>=2. There is no bound
on the base automorphism, dual polynomial degrees or D.
The more general orientation test proved in Section 4 says:
for f=D(b)+p(b) eta with a finite reduced base zero set Z
of at least two points and Jp nonzero on Z, such a
base-preserving unit-Hessian source member can exist only
if Jp takes one common nonzero value on all of Z.
Equality of its square alone is insufficient.

The retained inputs are:

- [HG0](GLOBAL_VALUE_ORBITS.md): actual critical schemes and
  full fibres survive source automorphisms. Unit Hessian in
  a full source orbit requires all native critical Hessian
  determinants to be one common nonzero number. Under
  gradient nonautomorphy the whole profile is
  E=L^3-NL+N_kappa L^2, L=uv.
- [HG2](GLOBAL_VALUE_ORBITS.md#3-the-two-known-endpoints-are-full-source-orbits-hc4myhg2):
  the coordinate and nondegenerate-quadratic source orbits
  are exactly the two old plane endpoint tests.
- [HR0](RELATIVE_CRITICAL_GRAPH_CLASSIFICATION.md#1-two-relative-critical-scheme-obstructions-hc4myhr0):
  the full-Hessian Jacobian-rank argument, the three
  projections of Tang's constant jointly affine plane, and
  the exact relative graph and rank-one formulas.
- [HC0](HODGE_TRACE_CLOSURE.md): additive/product Hodge
  invariant, applied only to explicit product strata.
- [HP0](QUADRATIC_ISOMETRY_SATURATION.md) and
  [HM0](AMBIENT_QUADRIC_GRAPHS.md): actual scalar/constant
  source-potential transfers.

The derivative rigidity and all polynomial inverse formulas below
are proved directly. No connectedness, cancellation, properness,
row normalization of every two-point Keller fibre, or invariance
of Hessians under nonlinear source maps is imported.

## 2. Critical schemes and every whole value fibre

Suppose first C is nonconstant and nonzero. Let r be its number
of distinct complex roots and s the number of its simple roots.
At a simple root rho, the full critical point is exactly

    x=y=0,   z=rho,   w=-D'(rho)/C'(rho).

Its critical value is D(rho), and its native critical Hessian is

    det Hess h=4 C'(rho)^2.                              (4)

At a repeated root rho the z derivative is D'(rho). If it
is nonzero there is no full critical point above that root.
If it is zero, the entire w line is critical and is forbidden
in every unit-Hessian source orbit. If all repeated roots
are of the former kind, the full critical set consists of
exactly s reduced Morse points.

The complete plane curve w C(z)+D(z)=kappa projects to
the z line. Off the r roots of C, w is uniquely solved.
At a root with D(rho)=kappa, the entire w line is free.
Writing a_kappa for the number of such distinct roots,
including repeated roots only once, gives

    E(w C+D=kappa)=L-r+a_kappa L.

For (1), the x!=0 stratum solves y and is C* times A2.
At x=0 retain that whole plane curve and the free y line.
Thus for every kappa,

    E(h=kappa)=L^3-rL+a_kappa L^2.                       (5)

This includes all special values, not only a generic fibre.
No multiplicities are substituted for distinct roots in (5).

If C has repeated roots with D' nonzero on all of them,
then r>s. Under hypothetical nonautomorphy HG0 would require
the generic profile L^3-sL, contradicting (5). Therefore
every unit-Hessian source member in this case has
automorphic gradient. If s!=1 even that alternative is
impossible, since a gradient automorphism has exactly one
zero. We do not infer that all s=1 source orbits are empty.
There is no need to classify the repeated-root degrees.

For squarefree C, N=r=s and a_kappa is exactly the count
of critical points at that value. Hence the **entire**
HG0 whole-value profile, including every special-value
jump, is satisfied. Repeating that test cannot close the
quadratic cases in (2).

### A short all-degree derivative rigidity

Let C be squarefree of degree n>=1. If (4) is one common
nonzero constant at every root, then

    C divides C'^2-c,       c!=0.

For n>=2 write C'^2-c=C Q, where deg Q<=n-2.
Differentiate and reduce modulo C:

    C'(2C''-Q)=C Q'.

Since gcd(C,C')=1, C divides 2C''-Q. Its degree is
strictly below n, so Q=2C'' and

    C'^2-2 C C''=c.                                     (6)

If C has leading coefficient a, the coefficient of
z^(2n-2) on the left is a^2*n*(2-n). For n>=3 it is
nonzero, contradicting (6). Thus n<=2.
The linear case is immediate. A squarefree quadratic
satisfies (6) with its nonzero discriminant.

This proves full source-orbit nonexistence for every
squarefree multiplier of degree at least three at once.
It is not a test of a growing sequence of degrees.

## 3. Exact normal forms and the four universal tests

If C is a nonzero constant, solve w to make h one
component of a polynomial automorphism: it is a coordinate.
If C=0 and D is nonconstant linear, solve z instead.
For C=0 with D constant or of degree at least two,
D' has a zero (or vanishes identically), giving a full
critical line or plane.

For C=a(z-rho), let

    Z=z-rho,  W=a w+[D(z)-D(rho)]/(z-rho).

The quotient is polynomial, the map has a polynomial inverse,
and h=2xy+ZW+D(rho) is nondegenerate quadratic.

For squarefree quadratic C, choose m,q,a with q*a!=0 and

    C(z)=a((z-m)^2-q^2),   D(z)=C(z)E(z)+ell z+d0.

Set

    Z=(z-m)/q,     W=a q^2(w+E(z)).

The inverse is z=m+qZ, w=W/(a q^2)-E(m+qZ), and

    h=2xy+W(Z^2-1)+alpha Z+beta,
    alpha=ell q,     beta=ell m+d0.                      (7)

For alpha!=0,

    h_alpha(x,y,z,w)
       =alpha*h_1(x/alpha,y,z,w/alpha).                  (8)

Both the map in (8) and its inverse are constant linear.
Nonzero scalar and constant changes preserve unit-Hessian
and gradient-inverse conclusions for the corresponding
functions. Equations (7),(8) identify full source orbits;
they do not assert native gradient conjugacy under a
nonlinear shear.

Each h_alpha has exactly two Morse critical points, both
of determinant 16. Thus a unit-Hessian member of either
full orbit in (2) would necessarily have nonautomorphic
gradient. Universal gradient invertibility on these two
orbits is therefore equivalent to absence of unit-Hessian
members there. Neither absence is proved in this note.

All other cases in Section 1 are either excluded or
inverse-closed. Conversely the four tests occur in (1).
This proves its exact reduction to four source-orbit tests.
It establishes no equivalence to a universal two-point
plane Keller problem and does not prove either old endpoint.

For the actual graphs (3), A=1 gives a coordinate and A=z
gives 2xy-2zw. A=z^2-1 gives h_0 up to a polynomial
source map and a constant; A=z(z-1) gives h_alpha with
alpha!=0 and hence the scalar source orbit of h_1.
Indeed D modulo C is the remainder of z^2 modulo A:
at two distinct roots rho1,rho2 its values are rho1^2,rho2^2.
They agree exactly when rho1+rho2=0.
Every other A is covered by the previous all-degree alternatives.
The four-test reduction is therefore exact for this closed
nondegenerate ambient graph family too.

## 4. All base-preserving source maps: orientation rather than its square

Use a two-variable base b and two free variables eta. Let

    f=D(b)+p(b) eta,

where p=(p1,p2), its base zero scheme Z is reduced and finite,
#Z=N>=2, and Jp is nonzero at every point of Z.
Its full critical set has one point over each b in Z:
the base critical equations solve eta using Dp(b)^T.

Let a polynomial source automorphism preserve C[b] setwise.
It has the form

    Phi(b,eta)=(sigma(b),tau(b,eta)),

where sigma is a plane polynomial automorphism and tau
is a relative polynomial automorphism, with no affine-linearity
assumption. Write delta=det Dsigma and
epsilon=det D_eta tau. Both are nonzero constants.
Put g=f composed with Phi and E=span(partial_eta).
Its full relative critical scheme and vertical Hessian satisfy

    Z_E(g)=sigma^(-1)(Z) times A2,
    (H_E g)|Z_E(g)=0.                                   (9)

The first identity follows from
grad_eta g=(D_eta tau)^T p(sigma(b)), using the global
polynomial inverse of that matrix. The second follows
because each second derivative of tau is multiplied by
one component of p. These are scheme identities.

Suppose g has unit Hessian. Its N full critical points
exclude an automorphic gradient, so apply the retained Tang
constant jointly affine plane W. Its projection to the
retained base has rank zero, one or two. Here are all cases.

- Rank two: constant coordinates preserving the retained base
  give g=A(z)+b^T K(z), det DK constant nonzero. The relative
  critical equations solve
  b=-DK(z)^(-T) grad A(z). Thus Z_E(g) is a graph
  isomorphic to A2, contrary to (9)'s N disconnected components.
- Rank one: the retained HR0 form is
  g=A(t,U)+s B(t,U)+V C0(t,U), with (B,C0) a plane Keller pair.
  The relative critical equations are
  C0=0 and A_U+s B_U+V (C0)_U=0.
  If (C0)_U were nonzero at any point of C0=0, solve V
  to obtain a relative critical point with nonzero vertical
  Hessian entry, contrary to (9).
  Hence (C0)_U vanishes on C0=0. The Keller component C0
  is squarefree, so C0 divides (C0)_U. Degree in U forces
  (C0)_U=0. The Keller identity then makes C0 affine-linear
  in t and B affine-linear in U with constant nonzero
  coefficient. The pair is triangular and invertible,
  which makes the full cotangent gradient invertible,
  the same contradiction.
- Rank zero: W=E. Thus g=A(b)+q(b) eta and Jq is a
  nonzero constant. This is the only remaining possibility.

In the last case differentiate

    q(b)=(D_eta tau)^T p(sigma(b)).

At a point b0 in sigma^(-1)(Z), every term differentiating
D_eta tau is multiplied by p=0. Taking determinants gives

    Jq(b0)=delta*epsilon*Jp(sigma(b0)).                   (10)

Because delta,epsilon,Jq are constants, Jp must have
one common nonzero value at every point of Z. This proves
the stated general orientation test. It controls the value
of Jp, not merely Jp^2.

For (1), p=(2x,C(z)), so Jp=2C'(z).
If C is squarefree of degree n>=2, its n roots cannot
all have C'(rho)=c0: the polynomial C'-c0 has degree at
most n-1 and n roots, forcing C' constant and n=1.
This contradicts (10), including the squarefree quadratic
where the two derivative values are opposite.

If C has a repeated root, the relative critical scheme
of (1) is V(2x,C(z)) times A2 and is nonreduced.
Every base-preserving Phi retains that scheme property.
But two selected rows of an invertible full Hessian are
everywhere independent. Their common zero scheme is smooth
of dimension two whenever nonempty, by the same Jacobian-rank
argument already retained in HR0. This contradicts
the nonreduced scheme. It excludes every repeated C as well.
Therefore every deg C>=2 is excluded under the entire
base-preserving automorphism class, proving JK1.

This includes arbitrary nonlinear base maps sigma and
arbitrary nonlinear dual maps tau over that base. It is
stronger than a fibre-affine matrix calculation. It still
does not constrain an arbitrary automorphism changing
the base algebra C[x,z].

## 5. What remains

There is no longer an open search in the degrees of C,D,
A or a base-preserving automorphism. For the suspended
univariate-multiplier branch the only additional full-source
tests beyond the old endpoints are the two fixed functions
in (2). Their entire special-value Hodge profiles pass,
and every base-preserving source route is excluded.
Any possible unit-Hessian source member must change the
base algebra; that is a structural requirement, not an
example or proof of existence.

General singular/nonsemisimple affine pencils, other image
geometry, full w-dependent source recovery, both actual
old endpoints and unrestricted HC4/JC2 remain open.
For standard Jacobian-rank context we rechecked the primary
[Stacks smooth presentation](https://stacks.math.columbia.edu/tag/00T6)
and [standard smooth criterion](https://stacks.math.columbia.edu/tag/00T7).
That criterion is already explicit in HR0; no new unproved
external theorem or formal verification is claimed.
