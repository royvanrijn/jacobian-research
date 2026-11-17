# Two-affine geometry: relative submersions and the exact endpoint branches

14 November 2026. Below the three-affine projection obstruction, two
global tests separate complete endpoint branches from new embedding
questions. Neither test requires a degree cutoff or a classification
of polynomial plane embeddings in three-space.

## 1. Two directional partials with no common zero force the cotangent plane (HC4MYHO0)

Let f be a complex polynomial on A4 with unit Hessian. Fix a constant
two-plane E of source directions and affine source coordinates
(s,t,u,v) with E=span(partial_u,partial_v). Then

\[
 (f_u,f_v)=(1)
 \quad\Longleftrightarrow\quad
 f=D(s,t)+uF(s,t)+vG(s,t),\quad
 J(F,G)\in C^\times,\quad(F,G)=(1).                         \tag{1}
\]

The forward implication is the substantive rigidity statement.
The right side supplies an actual relative coordinate completion:
choose aF+bG=1 and set

\[
 V=f,\quad R=-bu+av,\qquad
 u=a(V-D)-GR,\quad v=b(V-D)+FR.                             \tag{2}
\]

Thus a relative submersion in these two free variables is automatically
a polynomial coordinate over C[s,t], under the unit-Hessian hypothesis.
This is a global assertion: no common zero means the ideal in the
whole four-variable source ring is the unit ideal. Generation on
one chosen fibre or on a generic chart does not suffice.

Since (f_u,f_v)=(1), grad f misses zero, so it is not an automorphism.
Tang gives a constant jointly affine plane W. Its projection to
the chosen base (s,t) has rank 0, 1 or 2. Constant affine changes
preserving that base exhaust the following cases.

In rank two the normal form is

\[
 f=D(U,V)+sA(U,V)+tC(U,V),\quad J(A,C)\in C^\times.
\]

At any fixed (U,V), invertibility of D(A,C) solves for s,t making
both free derivatives f_U,f_V zero. This contradicts the hypothesis.

In rank one the normal form is

\[
 f=D(t,U)+sA(t,U)+V C(t,U),\quad J_{t,U}(A,C)\in C^\times.
\]

The nonconstant polynomial C has a zero (t0,U0). At it, A_U and C_U
cannot both vanish. Hence one can choose s,V solving
D_U+sA_U+V C_U=0. Together with f_V=C=0 this again produces a
common zero of the two free derivatives.

Only rank zero remains. Then W equals E and f is globally affine
in u,v in the actual source coordinates. Its Hessian determinant
is J(F,G)^2, and the free gradient has no zero precisely when
(F,G)=(1) in C[s,t]. This proves (1), with the explicit inverse (2).

There is an intrinsic recognition consequence. Under unit Hessian,
there exists a constant two-plane with unit directional-partial
ideal if and only if the full critical scheme of f is empty.
For the reverse direction Tang supplies W. Any zero of its coefficient
pair would solve the remaining two gradient equations uniquely in
the W variables, because the coefficient Jacobian is invertible,
and give a full critical point. Thus the pair has no common zero.
Such a plane is **unique**, by the rank argument
above applied to any other candidate plane. If the critical scheme
is nonempty, every directional-partial ideal vanishes there and none
is the unit ideal. In particular two distinct planes with unit
directional-partial ideals forbid unit Hessian. No unit-ideal plane
Keller pair is supplied, and its existence remains unresolved.

## 2. Constant-fibre geometry for a cotangent-aligned retained base (HC4MYHO1)

Consider a closed embedding

\[
 j(s,t,u,v)=(s,t,H_1,H_2,H_3),\quad
 C[s,t,u,v]=C[s,t,H_1,H_2,H_3],
\]

in a nondegenerate quadratic ambient A5. Suppose its actual potential
is affine in the free variables:

\[
                      f=D(s,t)+uF(s,t)+vG(s,t),\qquad
                      J(F,G)\in C^\times.                  \tag{3}
\]

Thus the retained base is the quotient by the cotangent plane;
this is an explicit alignment hypothesis, not a conclusion for
arbitrary two-affine embeddings.

The unretained three-dimensional quadratic kernel has rank r=1,2
or 3. After constant affine normalization and subtracting the ambient
critical value, the exhaustive forms and constant-fibre restrictions
are:

| r | Ambient form | Restriction on a zero of (F,G) |
|---|---|---|
| 3 | q2(s,t)+XY+Z^2, q2 nondegenerate | No zero is possible |
| 2 | sZ+t^2+XY | On s=0 a zero requires D(0,t)=t^2; its image fibre is one plane component |
| 1 | X^2+sY+tZ | The origin (s,t)=(0,0) is necessarily a zero; other zeros are not excluded |

The kernel rank is at least one because a nondegenerate form on
five-space has no totally isotropic three-plane. The displayed
normal forms follow by splitting its nondegenerate part, pairing
its radical with retained base directions, and absorbing cross
terms by affine coordinate changes.

At a zero of (F,G), the whole source fibre A2_(u,v) has constant
f value and embeds as a closed surface in the corresponding ambient
quadratic level in A3. If the level is irreducible, equality of
dimension makes the image the whole level.

For rank three, normalize the fibre equation to XY+Z^2=c. At c=0
it is an irreducible singular cone, so cannot be A2. At c!=0 its
compactly supported Hodge--Deligne polynomial, with T=uv, is

\[
                 (T-1)T+2T=T^2+T\ne E(A^2)=T^2.            \tag{4}
\]

Indeed X!=0 gives C* times A1 and X=0 gives two affine lines.
Only the retained additive/multiplicative Hodge invariant is used;
no quadric classification theorem is imported. Thus (F,G)=(1).

For rank two on s=0, a nonzero level XY=c is C* times A1, whose
ring has nonconstant units, unlike A2. At level zero the surface
is the union of two planes; the closed irreducible source image
must equal one of them. Away from s=0 the Z coefficient is nonzero
and the level is A2, so this argument supplies no exclusion there.

For rank one, at s=t=0 the source potential is affine in u,v but
equals X0(u,v)^2. Polynomial degree makes X0 constant, forcing
F(0,0)=G(0,0)=0. This proves attainment, not uniqueness of that
zero fibre. The singleton conclusion of the special p-graph family
must not be transferred to arbitrary two-affine embeddings.

The rank-three aligned branch is **exactly** the existing unit-ideal
endpoint, including existence. For any (F,G)=(1) and arbitrary D,
use V=f and R from (2), and take the nondegenerate ambient quadratic
Q=st+XY+Z^2. The explicit embedding

\[
                      j=(s,t,1,V-st-R^2,R)                 \tag{5}
\]

has potential V and polynomial source recovery by (2). Its unretained
kernel is nondegenerate of rank three, and it keeps s,t affine.
If (F,G) is Keller, the source Hessian is a unit and its gradient
misses zero. Conversely every embedding in the aligned rank-three
branch gives such a unit-ideal Keller pair. Thus unit-Hessian
existence in this branch, not just failure of gradient invertibility,
is equivalent to that unresolved endpoint. Alignment and closed
embedding are essential; arbitrary substitutions do not suffice.

## 3. Every affine-hyperplane image has only the two old endpoint possibilities (HC4MYHO2)

This classification allows arbitrary source parameterizations and
requires no retained affine source coordinates. Write the ambient
quadratic in centered normal coordinates as

\[
 Q(z)=\lambda_*+\tfrac12 z^T H z,\quad H\text{ invertible}.
\]

Suppose a closed embedding j:A4->A5 has image in the affine hyperplane
ell(z)=c, ell a nonzero constant row. Its closed image has dimension
four, so is the entire hyperplane; j is a polynomial isomorphism
onto it. Put

\[
                     n=H^{-1}\ell^T,\qquad \nu=\ell n.
\]

The exhaustive classification is:

| Hyperplane type | Restricted potential up to affine coordinates | Unit-Hessian source orbit |
|---|---|---|
| nu!=0 | Nondegenerate Q4 plus lambda_*+c^2/(2nu) | Exactly the singleton-zero Keller endpoint for nonautomorphic members |
| nu=0,c=0 | Rank-three quadratic with a critical line | No unit-Hessian member under any source automorphism |
| nu=0,c!=0 | Rank-three quadratic plus a nonzero linear radical coordinate | Exactly the unit-ideal plane Keller endpoint |

For nu!=0, the restriction to ker ell is nondegenerate and its
critical point is (c/nu)n, giving the displayed critical value.
The full source-orbit equivalence with the singleton-zero endpoint
is retained HC4MYHG2.

For nu=0, n is tangent to the hyperplane and is its quadratic
radical. The derivative of Q in direction n is ell(z)=c. If c=0,
the restricted critical scheme is the entire line Cn. Its pullback
cannot have unit Hessian. If c!=0, the radical coordinate occurs
linearly with a nonzero constant coefficient, so the restriction
is a polynomial coordinate with an explicit triangular completion.
Its full source orbit is all polynomial coordinates, giving exactly
the unit-ideal endpoint by retained HC4MYHG2.

Under the literal ambient Phi inverse, the same classification applies
to all Meng--Yang embeddings whose normal-coordinate image has an
affine linear relation. That entire branch supplies no new endpoint.
Only the isotropic hyperplane through the critical point is excluded
unconditionally; neither of the other two endpoint problems is solved.

## Remaining scope and evidence

The [checker](verify_two_affine_geometry.py) and
[receipt](receipt-two-affine-geometry.json) keep exact rank, critical-point,
quadric-stratum and source-inverse controls separate from these written
all-degree proofs. The only geometric inputs beyond elementary algebra
are the reproduced Tang alternative and the retained basic Hodge
invariant. The source-orbit endpoint equivalences are reused explicitly.

The two retained source coordinates alone do not force alignment in
Section 2. That section does not settle other positions of the intrinsic
plane or the remaining degenerate-kernel families. The no-common-zero
partial test in Section 1 does force alignment, but its existence is
already exactly the unit-ideal endpoint.

A new general embedding direction beyond the known endpoints must
avoid an affine relation among the five normal coordinates, have at
most two independent affine retained normal source coordinates, and
retain genuine global ring generation. Nonlinear image hypersurfaces,
the unaligned two-affine cases, the unrestricted SL2 row equation,
both plane endpoints and unrestricted HC4/JC2 remain open. No new
external classification, formal verification, independent replay,
external review or novelty priority is claimed.
