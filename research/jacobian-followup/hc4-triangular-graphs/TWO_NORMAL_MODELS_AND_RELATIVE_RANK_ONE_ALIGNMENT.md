# Two normal models and a relative rank-one alignment obstruction

15 November 2026. Claims **HC4MYJQ0/1**. The three resonant cubic
families in JP are not three separate universal tests. An explicit
quadratic-preserving polynomial shear identifies the length-two
family with the scalar family. All their parameters then normalize
to one fixed function. The length-three family normalizes to one
parameter and has an explicit global hyperbolic suspension form.

Thus the entire affine quadratic-composite universal inverse
question reduces to the four earlier tests together with

    K2=w(2xy+z^2-1)+z,
    K3_rho=w(2xy+z^2-1)+2x(rho-z)+z, rho in C.          (1)

These full source-orbit questions remain open. Neither model is
identified with h0,h1. For the three-critical-point family,
no equivalence with a one- or two-critical-point model is possible.

There is also a new all-degree alignment lemma: a relative critical
surface of vertical Hessian rank one that is not an affine plane
forces gradient inverse under unit Hessian. It excludes every
source reparameterization of K2 preserving its natural base algebra.
The hyperbolic suspension gives the corresponding base-preserving
exclusion for every K3_rho. A possible unit-Hessian source member
must genuinely change these base algebras.

## 1. Statements and retained foundations

Work over C. **JQ0** asserts the exact universal reduction to (1)
and the four earlier fixed tests, with unrestricted original outer
degree and unrestricted subsequent source degree. It also asserts
the global polynomial equivalence

    K3_rho ~ 2PV+B_rho(u,v),
    B_rho=u(uv+rho)^2-u+uv+rho.                         (2)

The generic plane curve B_rho=kappa is exactly A1 minus the three
distinct points -1,1,kappa, for kappa outside {-1,1,rho}.
Its three Morse points have critical values -1,1,rho, including
the coincident-value parameters rho=+-1. There is no assertion
that rho is a complete invariant of full source equivalence.

**JQ1.** Let g in C[b1,b2,eta1,eta2] have nonzero constant full
Hessian determinant, and let

    Z_E=V(g_eta1,g_eta2), E=span(partial_eta1,partial_eta2).

Suppose Z_E is nonempty, its vertical Hessian has rank at most
one everywhere and rank one somewhere, and Z_E is not isomorphic
to A2. Then grad g is a polynomial automorphism. It suffices
to establish E(Z_E)!=L^2 using the retained additive invariant.
No closedness assumption on its image in the retained base is
needed. This is an assertion about the constant plane E in the
actual unit-Hessian source coordinates, not a nonlinear distribution.

For K2 the retained base is C[x,w] and the free variables are y,z.
Every base-preserving source automorphism is excluded under unit
Hessian. For (2), the retained base is C[u,v], with P,V free;
every source automorphism preserving that base is likewise excluded.
Equivalently K3_rho in (1) has natural base C[w,y], which (9)
preserves exactly.
Both assertions allow arbitrary nonlinear base and relative maps.

Retain [JP2](ALL_AFFINE_PENCILS_TO_THREE_RESONANT_CUBIC_FAMILIES.md),
which gives the three starting cubic families and every reverse
inclusion in the actual affine quadratic-composite family.
[JM0](ZERO_EIGENVECTOR_SHEARS_AND_RANK_ONE_PENCIL_CLASSIFICATION.md)
supplies the literal global hyperbolic null shear; all maps below
also have their inverses written explicitly.

The critical chain rule and whole-value profile are
[HG0](GLOBAL_VALUE_ORBITS.md); the additive invariant is
[HC0](HODGE_TRACE_CLOSURE.md). Use the retained constant Tang
two-plane theorem [HCTANG2](../tang-hc4-audit/STRUCTURAL_REVIEW.md),
the three projection-rank normal forms and graph/envelope theorem
[HR0](RELATIVE_CRITICAL_GRAPH_CLASSIFICATION.md), and the relative
chain-rule argument already used in
[JK1](SUSPENDED_MULTIPLIER_ORBITS_AND_BASE_PRESERVING_CLOSURE.md).
HR0 uses the retained HC2 theorem, not a four-variable assertion
from de Bondt's low-dimensional theorem.

## 2. The length-two family equals the scalar family

Write q=2xy+z^2 and

    F_N=q+alpha x+e z+t(lambda q+N-c),
    lambda*c*e!=0, N=0,x^2,2xz.                         (3)

For any polynomial p=p(x,t), the actual source automorphism

    S_p: x->x, t->t,
         z->z+xp,
         y->y-zp-xp^2/2                               (4)

preserves q identically. Its inverse is S_-p, and its Jacobian
is one. The dependence on x,t is unrestricted, because those
coordinates are fixed. Let F_scalar=(1+lambda t)q+ez-ct.
Then, literally,

    F_scalar composed with S_(alpha/e)=F_0,
    F_scalar composed with S_((tx+alpha)/e)=F_(x^2).     (5)

These are global source equivalences, not a statement about
the native Hessians off the critical schemes. In particular
every length-two resonant source orbit is already a scalar
resonant source orbit. No cancellation theorem is used.

Choose r!=0 with r^2=c/lambda. The invertible constant map

    x=rX, y=rY, z=rZ,
    t=eW/(lambda r)-1/lambda                           (6)

gives

    F_scalar=e r K2(X,Y,Z,W)+r^2.

Its inverse is X=x/r,Y=y/r,Z=z/r,
W=lambda r(t+1/lambda)/e. Combining (4)--(6) handles every
alpha and every admitted nonzero parameter in both families.

## 3. The length-three family has just one remaining parameter

For N=2xz set r^2=c/lambda and

    x=e lambda X, z=rZ,
    y=r^2 Y/(e lambda)-rZ/lambda,
    t=eW/(lambda r)-1/lambda,
    rho=alpha lambda/(2r).                              (7)

This is an invertible constant affine source map. Its inverse is

    X=x/(e lambda), Z=z/r,
    Y=e lambda(y+z/lambda)/r^2,
    W=lambda r(t+1/lambda)/e.

Direct substitution into (3) gives

    F_(2xz)=e r K3_rho(X,Y,Z,W)+r^2.                    (8)

Thus all lambda,c,e are removed, not bounded or searched.
Every rho occurs by lambda=c=e=1,alpha=2rho. The old four
tests and both new normalized models are already contained
in the admitted affine family by JP's reverse inclusions.
Hence the universal reduction is exact. Literal
[HP0](QUADRATIC_ISOMETRY_SATURATION.md) and
[HM0](AMBIENT_QUADRIC_GRAPHS.md) transfers are retained.

Choosing the other square root r only changes a normal-form
representative. We make no claim of a full classification
of rho under source and value equivalences.

## 4. The coupled cubic is a genuine hyperbolic suspension

For K3_rho in (1) define the new source coordinates

    u=w, v=y, V=wy+rho-z,
    P=x+(wV-2w^2y-2w rho-1)/2.                          (9)

Their polynomial inverse is

    w=u, y=v, z=uv+rho-V,
    x=P-(uV-2u^2v-2u rho-1)/2.                         (10)

Writing T=uv+rho, substitution gives exactly (2).
Indeed the remaining part before the last shear is
B_rho+V(uV-2u^2v-2u rho-1), so (9) absorbs its entire
V-divisible part. No division by u, V or any other variable
occurs in this global source map.

The plane gradient is

    (B_rho)_v=u(2uT+1),
    (B_rho)_u=T^2-1+v(2uT+1).

Its complete critical set is therefore

    (u,v)=(0,1-rho^2),                 value rho,
          (-1/2,2(rho-1)),             value 1,
          (1/2,-2(rho+1)),             value -1.         (11)

The plane Hessian determinant is -1 at every point. The
suspended full Hessian determinant is four at all three.
The points stay distinct for every rho, including rho=+-1.

On a generic plane level set put T=uv+rho. Its exact inverse
parameterization is

    u=(kappa-T)/(T^2-1),
    v=(T-rho)(T^2-1)/(kappa-T),
    T in A1 minus {-1,1,kappa}.                         (12)

If kappa is outside {-1,1,rho}, every actual level point has
u!=0 and T avoids exactly those three values. Thus (12) is
an isomorphism with regular inverse T=uv+rho, not a merely
generic rational description. Its curve profile is L-3.

All exceptional levels can also be retained without covers:

    B_rho-1=(T-1)(u(T+1)+1),
    B_rho+1=(T+1)(u(T-1)+1),
    B_rho-rho=u(T^2+v-1).                              (13)

For rho!=+-1 each critical level is the union of two smooth
curves crossing at its one Morse point, with profile 2L-3.
At rho=1 its coincident level factors as
u v(u^2 v+2u+1); at rho=-1 it factors as
u v(u^2 v-2u+1). The two line components and the punctured-line
component form a chain with two distinct transverse intersections,
giving 3L-3. The other critical level still has profile 2L-3.
Thus every plane level has L-3+n_kappa L, with n_kappa the
number of its critical points. The hyperbolic suspension has
the retained whole profile L^3-3L+n_kappa L^2. These are
passing controls, not inverse closure of (2).

## 5. A relative rank-one alignment lemma without a closed base image

We prove JQ1. Under unit Hessian the two selected Hessian rows
are independent, so Z_E is smooth of dimension two. Suppose
the gradient is nonautomorphic. Tang supplies a constant jointly
affine two-plane W in these actual source coordinates.

If W projects to the retained base with rank zero, W=E and
the vertical Hessian is identically zero. This contradicts
rank one somewhere on Z_E.

If that projection has rank two, HR0's constant normal form is

    g=D(z)+b^T K(z), det DK in C*,
    Z_E={b=-DK(z)^(-T) grad D(z)}.

This is an actual polynomial graph isomorphic to A2. It
contradicts the assumed abstract nonisomorphism, for example
the sufficient certificate E(Z_E)!=L^2. No statement about
the image of that graph in the retained base is required.

For rank one the retained normal form is

    g=A(t,U)+s B(t,U)+V C(t,U), J(B,C) in C*,
    Z_E={C=0, A_U+s B_U+V C_U=0},
    det H_E=-C_U^2 on Z_E.

If C_U were nonzero at any point of C=0, the second equation
could be solved for V, giving a relative critical point with
rank-two vertical Hessian. This contradicts rank at most one.
Hence C_U vanishes everywhere on C=0. A Keller component is
squarefree: otherwise all its partial derivatives vanish on
a repeated factor, contradicting J(B,C) in C*. Thus C divides
C_U. Degree in U forces C_U=0. Now C=C(t), and the Keller
identity makes C affine-linear in t and B affine-linear in U
with constant nonzero coefficient. The plane Keller pair is
an explicit triangular automorphism. Its cotangent gradient
is then invertible, contradicting the assumed nonautomorphy.
This exhausts all three constant ranks and proves JQ1.

This lemma uses only the actual relative critical scheme and
vertical Hessian, not an unproved rectification of a varying
distribution. If the source automorphism changes the base
algebra, those relative data need not refer to a constant
plane in the actual unit-Hessian coordinates.

## 6. All natural-base source maps are excluded

For K2 use base (x,w) and free variables (y,z). Its relative
gradient ideal is

    (2wx,2wz+1)=(x,2wz+1),

since x=x(2wz+1)-z(2wx). Hence its full relative critical
surface is exactly Gm_w times A1_y, with z=-1/(2w),x=0.
It is smooth and has profile L^2-L, so is not A2.
Its vertical Hessian is diag(0,2w), of rank one everywhere.

Let Phi be any polynomial source automorphism preserving
C[x,w] setwise. It has an arbitrary plane base automorphism
and an arbitrary relative polynomial automorphism in y,z.
On the entire relative critical scheme the vertical chain
rule gives Hessian congruence; all extra terms have a vertical
gradient coefficient and vanish. Thus Phi preserves that
surface up to isomorphism and preserves its rank-one condition.
JQ1 forces gradient inverse for any unit-Hessian member.
But K2 has exactly two full critical points, so no source
member can have an automorphic gradient. This proves actual
unit-Hessian nonexistence throughout this base-preserving class.

For (2) use base (u,v) and free variables (P,V). Its relative
critical scheme is exactly P=V=0, a global section over A2.
Any base-preserving source automorphism retains that scheme
isomorphism over the base. HR0's global-section theorem forces
gradient inverse under unit Hessian. The three points (11)
exclude that possibility for every rho. Hence every source
map preserving this suspension base is excluded too.
In the cubic coordinates this is exactly preservation of C[w,y]:
the relative equations solve z=wy+rho and x=w^2y+w rho+1/2,
so the same global-section conclusion also holds directly.

The two and three full critical points rule out equivalence
between K2 and any K3_rho. They also rule out equivalence of
K3_rho with either old endpoint orbit or h0,h1. Equal counts
alone do not identify K2 with an old model. Its two distinct
critical values also exclude equivalence with h0, whose two
points share one value. Equivalence with h1 remains open;
matching its count and complete profiles supplies no source map.
The four earlier tests, the fixed K2
full source orbit, the K3_rho full source family, actual global
embedding recovery and unrestricted HC4/JC2 remain open.

We rechecked the retained primary
[standard smooth presentation](https://stacks.math.columbia.edu/tag/00T6),
[smooth criterion](https://stacks.math.columbia.edu/tag/00T7)
and [de Bondt v3](https://arxiv.org/html/1203.6605v3).
The exact HR0 envelope argument retains only its HC2 input.
No general cancellation, stable equivalence, local Morse
rectification or four-variable theorem is imported from those
sources. Native Hessians are not invariant under the nonlinear
maps above. The checker controls the written global identities
and scheme boundaries, not an independent/formal proof audit.
