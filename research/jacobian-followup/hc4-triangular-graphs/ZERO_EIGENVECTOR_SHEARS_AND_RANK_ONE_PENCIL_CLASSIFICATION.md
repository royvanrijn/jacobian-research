# Zero-eigenvector shears classify every rank-one affine pencil

15 November 2026. Claims **HC4MYJM0/1**. The inner quadratic
rank, affine offset and nonlinear outer polynomial are treated
without degree bounds. Every rank-at-most-one affine pencil
reduces its universal source-orbit inverse question to the same
four fixed tests already retained: the coordinate and quadratic
orbits, h0=2xy+w(z^2-1), and h1=h0+z.

For a rank-one pencil with nonzero eigenvalue, the only surviving
source orbits have a precisely stated quadratic outer resonance.
All other cases have no unit-Hessian source member. The two
resonant full source orbits and both old endpoint problems
remain open. Rank-two/rank-three singular pencils and
unrestricted HC4/JC2 are not classified here.

## 1. Two general reductions (JM0)

### Polynomial elimination of a null paired coordinate

For any polynomial K(x,z), independent of y, the function

    h=2xy+K(x,z)

has the exact polynomial source normal form

    h=2xY+K(0,z),
    Y=y+[K(x,z)-K(0,z)]/(2x).                            (1)

The quotient is a polynomial by evaluation at x=0. It is
independent of y. The inverse subtracts that same quotient;
the Jacobian determinant is one. No degree restriction is
placed on K. This is equality of full source-potential
orbits, not Hessian invariance or gradient conjugacy.

In particular let A=A^T be nondegenerate, T=A^(-1)B
self-adjoint, and g=q_B+2v^T A x+c. Suppose there is
a nonzero u with

    Tu=0,      A(u,u)=0,      A(v,u)=0.                  (2)

Choose u' with A(u,u')=1 and replace it by
u'-A(u',u')u/2. The resulting constant hyperbolic pair
has Gram matrix [[0,1],[1,0]], with nondegenerate
orthogonal complement. In its coordinates u=partial_y,

    q_A=2xy+q(z),       g=g(x,z) independent of y.

Indeed B u=A T u=0 and A(v,u)=0 kill its quadratic
and linear y derivatives. Equation (1), with
K=q+H(g), gives the actual global source equivalence

    q_A+H(g) ~ 2xy+q(z)+H(g(0,z)).                       (3)

Every inverse is polynomial; a generic or local coordinate
claim has not been substituted for a global source map.

### The complete critical-mass obstruction survives stabilization

Suppose a polynomial f has a finite nonempty Morse critical
set and

    sum_Crit(f) 1/det Hess f=0.                          (4)

For any nondegenerate quadratic Q in additional variables,
Crit(Q+f) is exactly {0} times Crit(f), up to a constant
translation if Q is affine quadratic. Each determinant is
det(Hess Q) times the corresponding determinant of f.
Thus (4) still holds, with a common nonzero factor.

If any arbitrary polynomial source reparameterization of
Q+f had unit Hessian, the native critical determinants
would all be one common C!=0. Their sum would be N/C!=0,
contradicting (4). Native degenerate critical points are
also preserved as obstructions. Consequently every full
JH0 critical-mass certificate remains an exclusion after
any nondegenerate quadratic stabilization, including source
maps mixing all old and new variables. This does not follow
merely from nonexistence of an unstabilized source member;
it follows from the actual complete critical-mass identity.

The retained [JH0](SIMPLE_STATIONARY_POLES_AND_NONCYCLIC_KRYLOV_CLOSURE.md)
proves (4) with all compatible simple stationary groups,
finite regular points and infinity accounted for.

## 2. Rank-one nonzero eigenvalue: the invariant classification (JM1)

In four variables let

    h=q_A+H(q_B+2v^T A x+c),
    A=A^T invertible, B=B^T, T=A^(-1)B,
    rank T=1, lambda=tr T!=0, d=deg H>=2.

A rank-one operator satisfies T^2=(tr T)T. Hence
E=T/lambda is the projection onto its one-dimensional
nondegenerate eigenspace, and I-E projects onto the
nondegenerate three-dimensional zero eigenspace. Put

    v1=Ev,        v0=(I-E)v,
    eta=A(v0,v0), nu=A(v1,v1),
    beta=c-nu/lambda.                                  (5)

Here eta is the norm on the zero primary space, not the
norm A(v,v) on the whole source.

There is no unit-Hessian member in the full polynomial
source orbit unless all four conditions hold:

    eta!=0,  d=2,
    H''=-1/(2eta),    gamma:=H'(beta)!=0.                (6)

If (6) holds, the entire source orbit is, up to a
nonzero scalar and a constant, exactly h0's if v1=0,
and exactly h1's if v1!=0. No existence of a
unit-Hessian member is asserted in either orbit.

In particular every d>=3 is excluded for every affine
offset. In the actual nonconstant graph
X=P0(g), Q=q_A+a0 X^2+b0 X+c0, a0!=0, one has
d=2deg P0. Thus every graph degree deg P0>=2 is
excluded in this nonzero-eigenvalue rank-one family,
throughout every polynomial source parameterization.

## 3. The zero-space norm and the genuine plane suspension

Normalize a coordinate z on the nonzero eigenspace so
q_A=z^2+q_0 and q_B=lambda z^2. Write v1=a partial_z,
so a^2=nu.

If eta!=0, v0 has a nondegenerate orthogonal complement
of dimension two inside the zero eigenspace. Over C
this complement contains a hyperbolic pair, constructed
as in Section 1. Choose w along v0 and normalize it so
v0=b partial_w with b^2=eta. The pair is orthogonal
to v0, so without any nonlinear change,

    h=2xy+f(z,w),
    f=z^2+w^2+H(lambda z^2+2az+2bw+c).                  (7)

If eta=0 and v0!=0, choose the isotropic zero eigenvector
u=v0 in (2). After (3), its component has no surviving
zero-line contribution; the residual inner polynomial is
lambda z^2+2az+c. If v0=0 choose any hyperbolic pair
in the zero eigenspace, with the same conclusion. The
reduced potential is

    2xy+w^2+[z^2+H(lambda z^2+2az+c)].

The bracket has degree 2d>=4. It cannot have a
unit-Hessian source member: if its derivative has a
multiple root there is a degenerate native critical
point; otherwise its second derivative would need one
common nonzero value at all 2d-1 critical roots.
Subtracting that constant leaves a polynomial of degree
at most 2d-2 with 2d-1 roots, forcing the bracket to
be quadratic. This is impossible. This is the retained
[HG1](GLOBAL_VALUE_ORBITS.md#2-two-all-degree-nonexistence-consequences-hc4myhg1)
argument, including arbitrary source mixing.

Thus eta=0 is completely excluded. For eta!=0 use
the actual suspension (7), with b!=0.

## 4. All nonresonant degrees are excluded by complete critical mass

For the plane f in (7), the exact stationarity data are

    P(s)=1+lambda s,
    z(s)=-a s/(1+lambda s),   w(s)=-b s,
    R(s)=beta-2b^2 s+a^2/[lambda(1+lambda s)^2],
    F(s)=s-H'(R(s)),        W(s)=1/[P(s)F(s)].            (8)

At alpha=-1/lambda, if a!=0 then stationarity is
inconsistent, R has a double pole, and W is regular
with zero order 2(d-1)-1>=1. If a=0 the determinant
root is simple and compatible. The full JH0 group of
2(d-1) critical points is retained; that pole is not
discarded or counted as one point.

At infinity, R=-2b^2 s+beta+O(s^(-2)). If d>=3,
F has polynomial growth of degree d-1 and
W=O(s^(-d)); in particular it has the required decay
of order at least two. F is not identically zero.

For d=2 write H'(r)=u r+v, u!=0. Then

    F=(1+2u b^2)s-(u beta+v)
                         -u a^2/[lambda(1+lambda s)^2].

If 1+2u b^2!=0, W=O(s^(-2)) and F is nonzero.
The finite alternatives and infinity condition are
exactly the complete JH0 boundary certificate.

A hypothetical unit-Hessian source member of the full
four-variable suspension forces its plane critical set
to be Morse. JH0 then gives a finite nonempty set with
zero complete reciprocal-Hessian mass. Section 1's
stabilization proves the contradiction in four variables.
No unit-Hessian assumption on the plane f itself is used.

Thus all d>=3 and every nonresonant quadratic are
excluded. The only remaining leading coefficient is
u=H''=-1/(2b^2)=-1/(2eta).

## 5. Every quadratic resonance has a literal multiplier normal form

Assume that leading resonance. Put u0=a/lambda,
Z=z+u0, beta=c-a^2/lambda. The inner becomes
lambda Z^2+2bw+beta. Write

    H(r)=-(r-beta)^2/(4b^2)+gamma(r-beta)+delta,
    gamma=H'(beta), delta=H(beta).

Expansion, with its inverse translation z=Z-u0, gives

    h=2xy+w C(Z)+D(Z)+u0^2+delta,
    C(Z)=-(lambda/b) Z^2+2gamma b,
    D(Z)=-lambda^2 Z^4/(4b^2)
                     +(1+gamma lambda)Z^2-2u0 Z.        (9)

These identities hold in the original source potential;
no gradient conjugacy is inferred.

If gamma=0, C has one double root at zero. If u0=0
there is a full critical w line. If u0!=0 there is
no full critical point, whereas the whole generic
profile of (9) is L^3-L. Unit Hessian would then have
nonautomorphic gradient and require L^3 by HG0,
a contradiction. Both cases are excluded. This is
also the explicit repeated-root case of
[JK0](SUSPENDED_MULTIPLIER_ORBITS_AND_BASE_PRESERVING_CLOSURE.md).

If gamma!=0, choose r!=0 with r^2=2gamma b^2/lambda.
The quadratic C has simple roots +-r. Divide
D=C E0+ell Z+d0. Its remainder has ell=-2u0,
since its even part is constant modulo C.
The global source map and inverse are

    zeta=Z/r,   Omega=-(lambda/b)r^2(w+E0(Z)),
    Z=r zeta,  w=-b Omega/(lambda r^2)-E0(r zeta).

They yield

    h=2xy+Omega(zeta^2-1)+alpha zeta+constant,
    alpha=-2u0 r.                                      (10)

Since the nonzero eigenspace is a nondegenerate line,
v1=0 iff a=0 iff u0=0. Thus alpha=0 exactly for
v1=0; otherwise the exact scalar source equivalence
h_alpha=alpha h1(x/alpha,y,zeta,Omega/alpha)
reduces to h1. This proves all assertions in (5),(6).

The two full source orbits remain open; JK1 excludes
all their base-preserving source maps and JK0 records
every passing special-value Hodge profile. This note
does not align a general full source map with that base.

## 6. All rank-at-most-one pencils give exactly the four old/fixed tests

If rank T=1 but tr T=0, the rank-one identity gives
T^2=0. If T=0 it is nilpotent too. The complete
[JF2 nilpotent affine-family reduction](AFFINE_QUADRATIC_PENCIL_RESIDUES_AND_ENDPOINTS.md#4-jf2-every-nilpotent-affine-pencil-has-only-the-two-old-endpoints)
already makes their universal inverse question exactly
the conjunction of the two old plane endpoints. This
is not an assertion that every individual nilpotent
potential is source-equivalent to a coordinate or quadric.

Combined with the classification above, the universal
gradient-inverse question for all rank-at-most-one
affine quadratic-composite graph images, over all
polynomial source parameterizations, is exactly the
conjunction of the old coordinate/quadratic tests and
absence of unit-Hessian members in h0,h1.

Each test occurs within this actual graph family:

- Q=x^2+y^2+z^2+w^2-X^2, X=x+1 gives a coordinate
  and has rank T=0.
- Q=2xy+z^2+w^2+X^2, X=x^2 gives the quadratic orbit
  by y->y+x^3/2, with rank-one nilpotent T.
- Put H(t)=-t^2/4+t/2 and
  h_a=2xy+z^2+w^2+H(z^2+2az+2w+a^2).
  This has rank-one semisimple T. The actual source
  automorphism Z=z+a, W=-w-Z^2/4+5/4, with inverse
  z=Z-a, w=-W-Z^2/4+5/4, gives

      h_a=2xy+W(Z^2-1)-2a Z+a^2+5/4.

  The choices a=0 and a=-1/2 realize h0 and h1.
  This is the nondegenerate ambient graph
  Q=q_A-X^2/4+X/2, X=z^2+2az+2w+a^2.

The constant ambient changes and source/scalar identities
are literal, as required by [HP0](QUADRATIC_ISOMETRY_SATURATION.md)
and [HM0](AMBIENT_QUADRIC_GRAPHS.md).
The inverse map of every stated graph is its retained-source
projection, so closed embedding and actual recoverability
have not been inferred just from a quadratic identity.

## 7. Exact scope after this reduction

[JL1](JORDAN_POLAR_CURVES_AND_ALL_INVERTIBLE_PENCIL_CLOSURE.md)
closes every invertible affine pencil. JM adds the entire
rank-at-most-one family to the same four fixed tests.
Any further affine-image mechanism beyond those tests
must have singular T of rank two or three, and still
escape the earlier centered/cyclic/nilpotent and general
critical-mass/geometric certificates.

This does not solve any of the four tests. The full
h0,h1 source orbits, both actual old endpoints, higher
rank singular pencils, full w-dependent embedding
recovery and unrestricted HC4/JC2 remain open.

For residue context we rechecked the retained primary
[Soprunov record](https://arxiv.org/abs/math/0511684v2).
No multivariable Euler--Jacobi, properness, stable
cancellation or new unproved external theorem is
imported. The stabilization uses the complete native
critical mass, and the zero-eigenvector shear has its
explicit polynomial inverse.
