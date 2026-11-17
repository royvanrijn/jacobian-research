# Additive base perturbations preserve the holomorphic cyclic obstruction

PF2RULED147. 17 November 2026. ALL degrees over C.
For EVERY quadratic MX0 surface and classical S2, every function

    P=a*w+f(z)+g(u), a in C*, deg f>=2,

has NO rational Jacobian mate, even in any finite extension of its
function field. The polynomial g is ARBITRARY in degree and coefficients.
This includes genuine u,z dependence and targets breaking fibre scaling.
No degree bound on the companion or target-field invariance is assumed.
General mixed coordinates, full finite-lattice integrality and HC4/JC2
remain OPEN. PF142's mixed 2+3 minimal alternative is not reopened.

## 1. General hypersurface statement and the actual coframe

Let V={w^2=N(u)+M(u)*z} be smooth over C, with

    m=deg M > ell=deg N >=1,
    Omega=du wedge dw/M.

The residue expression makes Omega a regular unit form on all V.
Choose a!=0, f,g polynomial, n=deg f>=2. Set d=deg g if g is
nonconstant, and d=0 if g is constant, including zero.

As in PF146, the geometric generic fibre P=p is the entire smooth
affine curve

    G=N+M*z-(p-f(z)-g(u))^2/a^2=0.                  (1)

Generic smoothness follows directly: every critical curve of a regular
function on a smooth characteristic-zero surface maps to one value;
there are finitely many components and isolated critical points.
Work componentwise if the geometric generic curve is disconnected.

Put K=N'+M'*z. Differentiating the ORIGINAL surface equation and P gives

    G_u=K+(2*w/a)*g'(u),
    Omega=dP wedge omega,
    omega=dz/(a*G_u).                               (2)

On the fibre use omega=-du/(a*G_z) where G_u=0. Smoothness of (1) makes
omega a regular NONZERO differential on the WHOLE affine generic curve.
A rational mate would force dQ=omega there.

We prove omega holomorphic at EVERY point of its smooth projective
normalization. Only after that complete audit do we use nonexactness.

## 2. Every finite-z infinity branch

At infinity, z cannot have a pole while u stays finite: the pole order
of f(z)^2 is greater than that of N+M*z. Thus u has a pole; put
alpha=-nu(u)>=1.

Suppose z is finite at this point. If 2*d<m its residue is zeta=0.
If 2*d=m, its residue MUST be the nonzero constant

    zeta=g_d^2/(a^2*M_m).

If 2*d>m, no finite-z branch exists. These follow by comparing the
highest poles of w^2, N and M*z; m>ell. Constant g belongs to the
first case.

Set Z=z-zeta and form the polynomial over C(p)

    R(u)=N(u)+M(u)*zeta-(p-f(zeta)-g(u))^2/a^2.

Its degree L satisfies

    1<=L<m, and L>=d when d>=1.                     (3)

The leading m-degree terms cancel in the equality case 2*d=m.
The generic coefficient p*g prevents cancellation below degree d.
For d=0, L=ell>=1.

Expanding (1), the other terms besides R+M*Z are divisible by Z and
have u-degree at most d. Since L>=d, they have strictly smaller pole
than R at a branch with nu(Z)>0. Consequently

    nu(Z)=(m-L)*alpha,
    nu(G_u)=-(L-1)*alpha.

Indeed the leading derivative coefficient is (L-m)*R_lead, nonzero.
Since dz=dZ, the EXACT differential order is

    ord(omega)=(m-1)*alpha-1 >=0.                   (4)

This includes branches where z approaches a NONZERO constant; retaining
only the vanishing-z branch from PF146 would miss them.

## 3. Polar-z branches away from a leading balance

Now let beta=-nu(z)>0. The term M*z strictly dominates N, so w has
pole order

    W=(m*alpha+beta)/2 >0.

Write J=d*alpha (zero for constant g) and H=n*beta.
The equation a*w+f+g=p gives a COMPLETE set of alternatives:
if J!=H, W=max(J,H); if J=H, W<=J.

In (2), K has pole order 2*W-alpha, and the g' term, when nonzero,
has pole order W+J-alpha.

* If J<W, K dominates. The exact order is again (4).
* If J>W, then J=H>W: f and g have cancelled their leading poles.
  The g' term dominates and

      ord(omega)=W-alpha+(n-1)*beta-1 >=0.          (5)

  Here W-alpha is a POSITIVE integer because
  2*(W-alpha)=(m-2)*alpha+beta>0.
* If J=W>H, the two derivative terms have equal order, but their
  combined leading coefficient is (m-2*d)*w^2/u. It is nonzero,
  since beta=(2*d-m)*alpha>0. Formula (4) still holds.

No bound on d was used. The only unresolved branch in this list is
J=H=W, treated next rather than excluded by a genericity assumption.

## 4. The resonant leading branch has at most a double root

In the balanced case J=H=W,

    m=(2*n-1)*t, d=n*t, t a POSITIVE INTEGER.        (6)

The integrality of t follows from gcd(n,2*n-1)=1.
Let v=1/u and use the local rescalings

    Z=v^t*z, T=v^d*w, M0(v)=v^m*M(1/v).

At the branch, Z,T are nonzero units. The surface equation is

    T^2=M0(v)*Z+v^(2*d)*N(1/v).

The last term vanishes at v=0. After choosing its nonzero leading T,
T(v,Z) is a formal regular unit, with T_Z=M0/(2*T) a unit. Define

    A(v,Z)=a*T(v,Z)+v^d*f(v^-t*Z)+v^d*g(v^-1).

It is a formal power series INDEPENDENT OF p, and the generic-fibre
equation is A(v,Z)=p*v^d. At v=0 its leading function is

    A0(Z)=a*T0(Z)+f_n*Z^n+g_d, T0^2=M_m*Z.

Every zero is simple or DOUBLE. At a multiple zero,

    a*T0=-2*n*f_n*Z^n,
    g_d=(2*n-1)*f_n*Z^n,
    A0''=n*(2*n-1)*f_n*Z^(n-2)/2 !=0.              (7)

Thus a higher leading multiplicity is impossible, not merely absent
in sampled coefficient choices.

The volume and required differential in these local coordinates are

    Omega=-v^(m-d-2)*(T_Z/M0)*dv wedge dZ,
    omega=v^(m-2)*(T_Z/M0)*dv/A_Z.                  (8)

At a simple leading zero A_Z is a unit and the order is nonnegative.
At a double zero solve A_Z(v,Zc(v))=0 coefficient by coefficient.
Taylor expansion gives

    A(v,Z)=B(v)+Y^2,

where B is in C[[v]], Y=(Z-Zc(v))*sqrt(U(v,Z)), and U is a formal UNIT.
This is a local identity obtained by a formal unit square root,
not a claim that a global comparison surface is affine space.

The curve now has the COMPLETE local equation

    Y^2=p*v^d-B(v).

If q is the v-order of its right side, then

    1<=q<=d.                                        (9)

All coefficients of B are in C, so the transcendental generic p cannot
cancel the coefficient at v^d. Normalize this binomial after a unit-root
extension. Every branch has
e=2/gcd(2,q), nu(v)=e, nu(Y)=e*q/2; for even q both branches are retained.
By (8), omega is a UNIT times v^(m-2)*dv/Y, of exact order

    e*(m-1)-1-e*q/2
      >= e*(m-1-d/2)-1
      >= e-1 >=0.                                  (10)

The final inequality uses n>=2 and (6):
m-1-d/2=((3*n-2)*t)/2-1>=1.

All formal rescalings were used only to compute normalized curve
valuations. No formal action or local coframe has been promoted to a
global target. Equations (4), (5), (10) cover EVERY infinity branch,
including leading cancellations, double roots and ramified branches.

## 5. Nonexactness, finite source extensions and cyclic applications

The required omega is therefore a NONZERO holomorphic differential
on every complete projective generic component. It is not rationally
exact: a pole of a rational primitive would create a differential pole;
with no poles the primitive is constant.

Normalized field trace rules out a rational mate in ANY finite
characteristic-zero extension of C(V). It traces a hypothetical exact
differential to a rational exact differential downstairs. It does NOT
assert polynomial trace landing. Equivalently every ramified pullback
has order e'*ord(omega)+(e'-1)>=0, with every index retained.

Apply the theorem to ALL quadratic MX0 models:

    M=u^2*b(u)^k,
    N=u*(-b+2*(1-k)*u*b')^2,
    k>=2, b square-free, b(0)!=0, r=deg b>=1.

Here m=k*r+2>ell=2*r+1, and PF146 retains the complete smooth source
image and q^*Omega=-2*dx wedge dt. Every
a*w+f(z)+g(u), deg f>=2, has no rational mate even in the WHOLE plane
source field. The companion need not descend to the source-image field.

For classical S2 use M=u^2,N=u, m=2,ell=1. The resonant branch (6)
cannot occur when n>=2, but the finite nonzero-z and leading-cancellation
branches still matter. This now excludes genuinely mixed S2 coordinates
that were not handled by PF146's h=0 line restriction.

For example P=w+z^2+u pulls back under PF144's S2 plane chart to

    t+t^3*h+(2*h+t^2*h^2)^2+t^2.

On t=0 and h=0 its degrees are BOTH two, and on an ordinary vertical
line its degree is four. No test of every polynomial source line is
claimed. The actual generic-fibre differential provides the certificate.
This P breaks the original fibre-scaling pencil.

## 6. Exact resonant controls and scope boundaries

The proof does not discard the balanced coefficient locus. A concrete
MX0 resonant example is k=4,b=u-1 and

    M=u^2*(u-1)^4, N=u*(1-7*u)^2,
    P=4*w+z^2+3*u^4-8*u^3+(20/3)*u^2-(40/27)*u.

Here n=2,m=6,d=4,t=2, leading Z=1,T=-1 is a double zero.
The formal critical value has
B(v)=(10/81)*v^4+O(v^5), so q=d=4 and the two local branches are
unramified. Formula (10) gives order TWO, despite the maximal generic
delay in (9).

The lower universal control M=u^3,N=u,
P=4*w+z^2+3*u^2 has m=3,d=2,t=1, q=2 and order ZERO.
It is included to test the sharp holomorphic boundary; it is not
asserted to be another original Keller normalization.

The nonlinear hypothesis is retained. On S2, P=w+z+u^2/4 is a polynomial
surface translate of z, and its required relative differential has a
nonzero logarithmic residue. Thus holomorphicness must not be asserted
for arbitrary linear f after adding g. Constant f gives actual rational
mates: P=w+g(u) has Q=1/u on S2.

Additive separation is ALSO a real boundary: on S2,
P=w+u*z has the rational mate Q=1/w. This multiplicative mixed control
lies outside the theorem; it is not a regular target pair because Q
is polar. General products in u,z, nonconstant w coefficients and
higher w dependence remain separate. No reduction of all target
coordinates to the tested family is supplied.

The exclusion survives all surface automorphisms of these cyclic models
and all polynomial target rechartings, since unit volumes differ by
constants. It rules out any ACTUAL target coordinate in those transported
families, with no bound on its companion.

## 7. Original finite normalization and remaining goal

Retain the original B,R,F,L, FULL finite normalization S and trace-dual
Lambda, including every ramified and unramified boundary prime.
The local Y,v model in Section 4 is a completed GENERIC CURVE branch,
not S, an affine intermediate plane, or a replacement order. Field
trace was used for rational exactness only, never Tr(R) subset B or
Tr(Lambda*Lambda) subset B.

PF142's mixed 2+3 minimal case remains retired. PF143's actual generic
first-velocity marked quadratic integrality is still UNSUPPLIED.
The present theorem supplies an actual coframe obstruction for an
unbounded genuinely mixed cyclic target-coordinate family. It does
not close all cyclic targets, HC4/JC2, or construct a counterexample.

The universal proof is (1)--(10), including the formal double-root
normalization reproduced above. Supporting exact controls check the
coframe identity, full valuation cases, the double-root calculation,
the resonant examples and boundary countercontrols. No bounded target
census, new external theorem import or predecessor checker replay is
used. Independent whole-proof replay, formal verification, external
review and novelty priority remain unclaimed.
