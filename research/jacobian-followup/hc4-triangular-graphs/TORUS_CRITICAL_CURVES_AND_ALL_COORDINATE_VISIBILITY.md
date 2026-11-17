# Torus critical curves and complete one-coordinate visibility barriers

Theorems HC4MYJY0/1/2, over C. The previous turn contracted the universal
affine frontier to the two old endpoint inversion questions, absence of
an actual two-point Keller fibre, and the whole S_rho family. This note
does not solve those questions. It closes the remaining individual
coordinate-preserving classes by an all-degree global curve argument and
the retained relative Hessian invariant.

## 1. A two-point torus critical curve cannot survive (HC4MYJY0)

First consider a plane Keller map K=(F,G). Suppose the **whole** G=0
scheme is isomorphic to Gm. Suppose also F has exactly two zeros on it.
This is impossible.

Choose its actual Laurent coordinate t, so O(G=0)=C[t,t^-1]. The Keller
differential makes the restriction of F etale: on the smooth G=0 curve
its differential never vanishes. The units of the Laurent ring are c t^k.
Thus dF/dt is a monomial. Integrating a Laurent polynomial in characteristic
zero gives

    F|_(G=0) = b + a t^m,  a!=0, m in Z minus {0}.

The logarithmic case m=0 cannot be a Laurent polynomial. Exactly two
zeros force b!=0 and |m|=2. Every fibre K^-1(alpha,0), alpha!=b, therefore
has exactly two points. No extra point or component is allowed, because
G=0 was the whole fibre, not one component.

We import the precise statement in the Remark after Theorem 1 of
[Nguyen Van Chau, math/0408048v1](https://arxiv.org/html/math/0408048v1):
the exceptional value curve of a complex nonsingular polynomial plane map
has no irreducible component isomorphic to A1. For a Keller map its critical
value set is empty, so this is its nonproper-value set. Consequently the
target line G=0 is not contained in that set. A generic point of this line
lies outside it, where the map has the generic field degree number of
preimages. That degree must be two. This passage from a whole line of
two-point fibres to generic degree uses the no-line theorem; a single
two-point fibre does not determine generic degree.

But a characteristic-zero extension of degree two is abelian Galois,
contradicting retained [HW0](ABELIAN_CLASSGROUP_OUTPUT_OBSTRUCTION.md).
This proves the plane assertion. The primary paper also states that the
nonproper set is an algebraic curve and the map off its exceptional set
is a covering; no unknown global Keller-surjectivity statement is imported.

Now let f(s,z1,z2,z3) have exactly two full reduced Morse critical points.
Assume its whole three-relative critical scheme

    C_f = V(f_z1,f_z2,f_z3)

is Gm and its projection to s is nonconstant. Then no polynomial source
automorphism preserving C[s] setwise can produce unit Hessian.

Normalize a base affine automorphism so s is fixed. Relative critical
schemes and their s projections transport by the actual inverse, as in
[GK0](RELATIVE_CRITICAL_LOCUS.md). Suppose a unit member exists. Its two
full critical points prevent its gradient from being an automorphism.
The constant Tang plane projects to the retained s direction with rank
zero or one, precisely as in [GR1](CRITICAL_COMPONENTS.md).

Rank zero gives d(s,U)+P F(s,U)+Q G(s,U), with a Keller coefficient pair.
Its relative critical curve projects to the finitely many s values at
the common zeros of F,G. This contradicts the nonconstant projection.

Rank one gives, after constant affine changes retaining s,

    g = d(U,Z) + s F(U,Z) + V G(U,Z), J(F,G) in C*.

The whole three-relative critical scheme projects isomorphically onto
G=0: its other two equations solve uniquely and polynomially for s,V.
Thus G=0 is exactly Gm. The full critical scheme is obtained by also
setting g_s=F=0, so F has exactly two zeros on this torus. The plane
assertion rules this out. This proves the four-variable assertion in
every source degree. The retained coordinate need not be the torus
parameter, only a nonconstant function on the whole torus.

The word whole is essential. An additional component of G=0 could supply
other points over the target line and invalidate the generic-degree count.
The theorem neither excludes an arbitrary two-point Keller fibre nor
changes the exact open E2 gate in JX1.

## 2. All four native K2 coordinates must become nonlinear (HC4MYJY1)

For K2=w(2xy+z^2-1)+z, the relative scheme retaining w is exactly

    C_w = V(x,y,2wz+1) = Gm_w.

Indeed the relative derivatives are 2wy,2wx,2wz+1, and the last equation
makes w a unit. Projection to w is nonconstant. Its two full critical
points are (x,y,z,w)=(0,0,r,-1/(2r)), r=+-1, both Morse. Section 1
excludes every w-preserving source change. The native determinant on
this entire curve is the passing constant 4; this is a genuine global
obstruction beyond the relative determinant test.

Retaining z instead, the relative equations are

    2wy=0, 2wx=0, 2xy+z^2-1=0.

They contain the two-dimensional closed locus w=0,2xy+z^2=1. Under unit
Hessian the three relative derivative rows are independent rows of the
invertible Hessian, so the entire nonempty relative scheme must be smooth
of pure dimension one. The dimension-two locus is impossible, and it
transports under any z-preserving source automorphism.

Together with the x/y exclusion in
[JW2](GENERIC_COORDINATE_SPECIALIZATION_AND_K2_ONE_COORDINATE_CLOSURE.md),
this excludes unit Hessian whenever **any one** of Phi^*x,Phi^*y,Phi^*z,
Phi^*w is affine-linear. To see this extension, normalize that nonconstant
affine pullback to the corresponding native coordinate by an actual
constant affine source precomposition. Only this affine normalization
uses full Hessian congruence away from a critical locus. Every remaining
K2 unit candidate must move all four native pullbacks nonlinearly.

For h1=2xy+(z^2-1)w+z, retaining w gives the same whole relative torus
V(x,y,2zw+1), so section 1 excludes every w-preserving map as well.
For h0=2xy+(z^2-1)w, its w-relative scheme is V(x,y,zw), two lines
meeting at the origin. This singular scheme cannot be the smooth
one-dimensional relative scheme of a unit Hessian. Thus h0 also has
no w-preserving unit source. Together with JR's x/y exclusions, a unit
member of h0 or h1 must make x,y,w pullbacks nonlinear. In contrast their
z-preserving unit-existence question remains exactly E2 by
[JX1](GLOBAL_SIGN_TWISTS_AND_TWO_POINT_KELLER_EQUIVALENCE.md).

K2 remains redundant in the universal affine test list by JX2. Closing
its last coordinate-preserving classes does not resolve its individual
full-source orbit or identify it with h1 or C2.

## 3. All four coordinates in the S_rho frame must become nonlinear (HC4MYJY2)

Use the actual polynomial suspension frame of
[JQ0](TWO_NORMAL_MODELS_AND_RELATIVE_RANK_ONE_ALIGNMENT.md):

    S_rho = 2PV + B_rho(u,v),
    T=uv+rho, B_rho=u T^2-u+uv+rho.

It has three full Morse critical points, each full Hessian determinant 4,
for every complex rho, including rho=+-1.

Retaining u, the whole relative critical scheme is

    P=V=0, u(2uT+1)=0.

Its two factors generate comaximal ideals. One component is the actual
closed affine line u=0,v free. The other is exactly Gm_u, with
v=-1/(2u^2)-rho/u, and projects nonconstantly to u. GR1 applies and
excludes unit Hessian: a unit member would have an invertible gradient,
contradicting its three full critical points. This uses the whole
line-plus-torus geometry, not a rational-curve substitute for an affine
line component.

Retaining v, it suffices to use GK0's determinant class on the relative
scheme. Its equations are P=V=0,B_u=0, where

    B_u = 3u^2v^2+4rho uv+v+rho^2-1.

For arbitrary rho it contains the following entire open curve, with
regular inverse t=uv+rho:

    Y(t)=1-3t^2+2rho t,
    v=Y(t), u=(t-rho)/Y(t), Y(t)!=0.                    (1)

On it the full Hessian determinant is

    4 H(t)/Y(t)^2,
    H(t)=-3t^4+4rho t^3+6t^2-12rho t+4rho^2+1.        (2)

If a v-preserving unit member existed, GK0 would force the native
determinant class on the whole relative scheme to be one constant.
The full critical points belong to that scheme and have determinant 4,
so the constant must be 4. Equation (2) on the open curve would force
H=Y^2 as polynomial identities. Their leading coefficients are -3 and
9, respectively, an immediate contradiction for every rho. No generic
parameter exception, ramification classification or finite rho sample
is needed; Y is never the zero polynomial and its complement is dense.
The computation chart (1) is used only on this actual relative curve,
not asserted to be a global source automorphism.

The earlier JR1 theorem excludes all source maps retaining P or V in
this same actual suspension frame. Combining these four exclusions and
constant affine normalization, a unit member S_rho composed with Phi
would require all four Phi^*u,Phi^*v,Phi^*P,Phi^*V to be nonlinear.
The original K3-to-S map retains w=u and y=v, so the new u/v results
also exclude all K3 source maps with either native w or y pullback
affine. The P,V frame coordinates are nonlinear in the original K3
coordinates; we do not turn them into an unproved native x/z exclusion.

These are visibility barriers, not a theorem forcing an arbitrary source
map to retain an affine pullback. The global three-point family, E2, both
old endpoints, arbitrary HC4 embeddings and unrestricted HC4/JC2 remain
open. The next step is intrinsic global alignment or a genuinely global
Keller/embedding obstruction, not further coordinate-preserving searches.

## Reproduction and assurance

[The checker](verify_torus_coordinate_visibility.py) binds
[the primary and retained inputs](sources/torus-coordinate-visibility-context.json)
and checks exact native relative ideals, the torus and curve inverse
identities, the parameter-uniform determinant formula and leading-term
contradiction, all critical points, and the passing torus determinant
control. [Receipt](receipt-torus-coordinate-visibility.json).
The generic-degree argument is proved above and uses the explicitly
imported no-line statement with its plane/nonsingular hypotheses checked.
The receipt is a local symbolic replay, not independent replay, formal
verification or external review. No arbitrary two-point fibre exclusion
or full-source rectification is inferred from these checks.
