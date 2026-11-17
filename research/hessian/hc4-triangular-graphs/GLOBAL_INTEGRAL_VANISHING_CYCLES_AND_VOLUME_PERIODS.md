# Global integral vanishing cycles and the volume-period barrier

Theorems HC4MYKA0/1/2, over C. The complete middle-homology volume
periods of every remaining normalized model pass the necessary cotangent
test. This includes the integral lattice, its monodromy, and the
point-labelled local vanishing cycles, not only separate local expansions.

It does not identify Seifert pairings, thimble linking, actual singular
fibre gluing, torsor projections or source orbits. Those data are not
included in the period object considered here.

## 1. Integral suspension cycles and their periods (HC4MYKA0)

Let p(z,w) have smooth zero curve C, and let

    X={a x y+p(z,w)=0}, a in C*, Omega=dx dy dz dw.

The topological realization of the retained smooth-suspension/purity
construction in [JZ1](SMOOTH_VALUE_FAMILIES_AND_THE_MOTIVIC_HOMOTOPY_BARRIER.md)
gives an integral suspension isomorphism

    tau:H1(C,Z) -> H3(X,Z).                              (1)

Here is the cycle realization and normalization of (1). Bound a compact
one-cycle gamma on C by a real two-chain D in the contractible plane.
Over its interior take the circle in axy=-p, and collapse that circle
where p=0. The resulting tube over D represents tau(gamma). Tubular
excision identifies its normal circle with the positively oriented
complex normal circle in the purity construction. Thus it is the
integral suspension map, with multiplicity one, not merely a rational
homology map. Different bounding chains give homologous tubes. This
also follows by applying the same circle construction to a bounding
three-chain for their difference in the plane.

On x!=0 the Gelfand--Leray form satisfying dh wedge omega=Omega for
h=a x y+B(z,w) is

    omega= -dx/(a x) wedge dz wedge dw.                  (2)

For a primitive lambda with d lambda=dz wedge dw, circle integration
and Stokes therefore give

    integral_tau(gamma) omega
       = -(2 pi i)/a * integral_gamma lambda.           (3)

One may equivalently take lambda=-w dz. Formula (3) is independent
of the primitive and is an exact contour formula, not an expansion.
The collapsed fibres have measure zero in this computation; remove
small neighbourhoods and take their limit. The form is regular on X,
so the limit agrees with the cycle period.

The construction is natural in smooth families. In a local plane node
AB=epsilon the positive A-circle is the local vanishing circle; its
tube is the generator of the suspended four-variable Morse vanishing
sphere. In local quadratic coordinates this is S1*S1=S3, with degree
one. Hence (1) carries local plane vanishing circles to the corresponding
local three-cycles, up to their independent orientation choices.
This compatibility is used below; no classification of singularities
at infinity is substituted for it.

For a plane Keller pair (F,G) with N zeros b_i, consider

    H=d(b)+P F(b)+Q G(b), j=Jac(F,G) in C*.

On every smooth value s avoiding c_i=d(b_i), projection H=s to
U=A2 minus {b_i} is the actual affine-line torsor from JU/JZ.
Its contractible fibres identify H3 with H3(U,Z). Small three-sphere
links of the N points are an integral basis: the punctured plane
retracts to a wedge of N three-spheres. The torsor projection is over
the same constant U on the entire noncritical value open, so these
classes are flat and monodromy is the identity.

Near b_i use u=F,v=G and absorb d-c_i into P,Q as in JU2. The local
potential is c_i+uP'+vQ', and the volume is j^-1 du dv dP' dQ'.
Its Morse vanishing sphere projects with degree one to the link of
b_i. Consequently the N local vanishing cycles themselves are that
integral basis, with exact periods

    (2 pi i)^2 * (s-c_i)/j,                             (4)

after orientations are chosen. Quadratic integration fixes the common
factor in (4). Thus the entire period functional on H3 is specified
by N exactly linear functions with one common nonzero slope. Coincident
critical values keep distinct basis vectors, even when their period
functions coincide.

Multiplying the ambient constant volume form by a scalar multiplies
all periods by that one scalar. Actual polynomial source automorphisms
do exactly this. This is therefore an appropriate necessary global
period object for full source equivalence, not an invariant of one
chosen quotient alone.

## 2. Complete native cycle bases and exact residues (HC4MYKA1)

### Squarefree suspended multipliers in every degree

Let

    h=2xy+C(z)w+D(z),

with C nonconstant squarefree and D arbitrary. Let r_1,...,r_N be
all the roots of C, and restrict s away from c_i=D(r_i). The whole
plane curve Cw+D=s is the punctured line

    w=(s-D(z))/C(z), z in A1 minus {r_i}.

Let e_i be the positive puncture loops. They are an integral basis
of H1. At s=c_i the plane node has coordinates A=z-r_i and

    E=(C(z)w+D(z)-c_i)/(z-r_i).

The derivative of E in w at the critical point is C'(r_i)!=0,
so these are actual local analytic coordinates and h_plane-c_i=A E.
An A-circle with radius tending to zero and |s-c_i|/radius tending
to zero is the local vanishing circle. It is exactly e_i. This remains
true when several c_i coincide: the nodes lie at different r_i.

By (1), all N local three-cycles tau(e_i) are an integral basis of
the whole H3. The punctures are fixed and their labels give identity
monodromy. From (3) with primitive -w dz, their exact periods are

    (2 pi i)^2 * (s-D(r_i))/(2 C'(r_i)).                (5)

This is an arbitrary-degree formula for the full lattice. In the
remaining C=z^2-1 models, choose delta_r=r*tau(e_r), r=+-1.
Then both slopes in (5) become the same 1/4.

### The fixed K2 model

Use the exact whole-family JS presentation

    2xy+F_s(u,v)=0, F_s=v+(s+uv)^2-1,
    z=s+uv, u=w, v=1-2xy-z^2.

The whole exceptional curve is the same punctured line z!=+-1, with

    v=1-z^2, u=(z-s)/(1-z^2).
    s=z+u(z^2-1).                                     (6)

Thus its total curve family is actually the plane (u,z) with value
function z+u(z^2-1). For r=+-1 the exact local factorization is

    s-r=(z-r)[1+u(z+r)].

At its node z=r,u=-1/(2r), the two factors are local coordinates.
The vanishing circle is e_r, the puncture loop around z=r.
Therefore their two tubes form the full integral H3 basis, again
with identity monodromy. This curve-family identity does not assert
a four-variable K2/h1 source equivalence.

The native four-volume must be kept, since the presentation uses s.
On the smooth K2 fibre, direct substitution into Omega/dh gives

    omega_K2=dx/(2x) wedge du wedge dv.                 (7)

Indeed dz wedge dw=-u du wedge dv on a fixed s fibre, and the native
dy coefficient of dh is 2ux. The apparent u denominator cancels.
For lambda=u dv the contour formula therefore has a plus sign.
On the curve (6),

    lambda=2z(z-s)/(z^2-1) dz,
    Res_(z=r) lambda=r-s.

Choose delta_r=-tau(e_r). Formula (7) gives the exact global periods

    (2 pi i)^2 * (s-r)/2.                              (8)

No local square-root completion is needed to obtain this global result.

### The complete three-point family, including collisions

In the actual suspension frame put

    S_rho=2PV+B_rho(u,v),
    T=uv+rho, B_rho=u T^2-u+uv+rho.

Its noncritical value open is s notin {-1,1,rho}. The whole plane
curve is parametrized by

    t in A1 minus {1,-1,s},
    u=(s-t)/(t^2-1), v=(t-rho)(t^2-1)/(s-t).

Let e_+,e_-,e_s denote the three positive puncture classes.
They form a global flat integral basis: moving the labelled puncture
s around a fixed puncture is a pure braid and fixes these homology
classes. This statement concerns homology, not the fundamental group.

The point-labelled vanishing circles are

    gamma_+=e_+ + e_s, gamma_-=e_- + e_s, gamma_rho=e_s. (9)

Here is a local check that also handles rho=+-1. The exact identities

    B_rho-1=(T-1)[u(T+1)+1],
    B_rho+1=(T+1)[u(T-1)+1],
    B_rho-rho=u(T^2+v-1)

give local product coordinates at the three distinct Morse points.
For the first two, the first factor t-r is a local annulus coordinate.
Its small circle encloses both punctures r and s when s approaches r,
so its class is e_r+e_s.

For the third point choose u as the local annulus coordinate. The
equation u(t^2-1)+t=s solves t analytically near u=0, since its t
derivative there is one. On a small disk in u, t-s=-u(t^2-1) has a
simple zero at u=0 for s!=+-1, while t never equals either fixed
puncture. Its annulus circle therefore winds once around s and
around no fixed puncture. This proves gamma_rho=e_s also when
rho=+-1. At a collided critical value the corresponding nodes and
local annuli are still distinct.

The row matrix in (9) is

    [[1,0,1],[0,1,1],[0,0,1]], determinant one.

Thus the three local cycles form a unimodular basis of the entire
H3 after suspension. There is no hidden integral index or additional
middle-homology period to test.

For the native source order (P,V,u,v), equivalently (u,v,P,V),
formula (3) uses the primitive -v du and coefficient 2. Directly,

    v du = [(rho-1)/(t-1)+(rho+1)/(t+1)
                         +(s-rho)/(t-s)-1] dt.         (10)

Its residues in the puncture basis are rho-1,rho+1,s-rho.
Applying (9) therefore gives the three exact suspended periods

    (2 pi i)^2/2 * (s-1, s+1, s-rho).                  (11)

The identical formula holds for every fixed rho, including rho=+-1.
For example at rho=1 the first and third basis vectors have the same
period, and their nonzero difference has identically zero period.
This is allowed and occurs equally for cotangent data with repeated
critical values. One must not mistake that zero period for a missing
cycle or a failure of primitivity.

## 3. The complete global period comparison passes (HC4MYKA2)

The coordinate endpoint has H3=0. The nondegenerate quadratic endpoint
has one primitive Morse cycle and its familiar exactly linear quadratic
period. Equations (5), (8) and (11) treat every other normalized test.

After choosing orientations, their finite-point-labelled vanishing
cycles are integral bases with identity monodromy. Their period vectors
are one common nonzero constant times

    (s-c_1,...,s-c_N).

This is precisely (4), allowing the single volume scalar permitted by
an actual polynomial source automorphism. The local first coefficient
is fixed by the native critical Hessian determinant; its common value
was already the HG necessary condition. No separate scaling at different
critical points is used here.

Every global integral cycle is an integer combination of this basis.
Its period is therefore the corresponding integer combination of the
displayed linear functions. Analytic continuation creates no additional
periods or monodromy mixing. Consequently the integral H3 local system
together with its point-labelled vanishing classes and its native-volume
period functional supplies no new full-source exclusion of the
remaining normalized models. It is stronger than the separate local
JU2 statement and adds information not asserted by an unmarked
homotopy-type comparison alone.

This conclusion is limited to this stated period object. It does not
compute a Seifert form, linking of chosen thimbles, a compactification
at infinity, relative homology marking an external boundary, or actual
singular-fibre gluing. Nor does equality of periods supply an actual
quotient action, globally Keller coefficients or polynomial coordinates.
The two old endpoints, E2, the full three-point family, individual
full-source orbit questions, actual graph embedding recovery and
unrestricted HC4/JC2 remain open.

## Reproducibility and assurance

[The checker](verify_global_vanishing_periods.py) verifies exact residue
identities, native-volume chart conversion, all three product-coordinate
Jacobians, the unimodular cycle matrix and collision controls.
[The source context](sources/global-vanishing-periods-context.json)
pins the retained homotopy/suspension inputs, local-period proof, whole
curve families and current frontier. The integral tube construction and
cycle identification are written proofs, not conclusions from symbolic
samples. [Receipt](receipt-global-vanishing-periods.json).
No independent replay, formal verification, external review or
completeness of more refined marked invariants is claimed.

