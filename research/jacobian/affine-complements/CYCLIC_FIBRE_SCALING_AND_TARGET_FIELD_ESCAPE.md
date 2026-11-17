# Cyclic fibre scaling cannot preserve an actual plane target field

PF2RULED144. 16 November 2026. ALL degrees over C.
This excludes every fibre-scaling-equivariant plane target in the ENTIRE
MX0 cyclic source-image family, including nonlinear target actions and
arbitrary target degrees. Section 6 also treats every classical S_g pseudo-plane.
Target-field invariance is NOT automatic.
The unrestricted cyclic target and full HC4/JC2 goal remain OPEN.

## 1. The regular fibre action is different from the generic deck action

Retain MX0's g,k>=2 and a positive-degree square-free b(u), b(0)!=0.
Use different names for its TWO commuting finite actions. Put

    u=t^g, a(t)=t^2*b(u)^k,
    sigma(t)=-b(u)+g*(1-k)*u*b'(u),
    W={a(t)*z=v^g-sigma(t)^g}.

The DECK group D=mu_g acts freely on W by

    d_epsilon(t,v,z)=(epsilon*t,epsilon*v,epsilon^-2*z).

Its quotient V=W/D is smooth affine. MX0's selected whole plane chart

    v=sigma+a*x, z=(v^g-sigma^g)/a

gives q:X=A2_(x,t)->V SURJECTIVE ETALE of generic degree g.
Let R=C[x,t], L=C(x,t), C=O(V)=R intersect E, E=C(V). The entire line
L0={t=0} maps ISOMORPHICALLY onto the reduced multiple-fibre line C0
of V->A1_u. Exact unit volume Omega on V pulls back to dx wedge dt.

There is ANOTHER regular action on W:

    rho_gamma(t,v,z)=(t,gamma*v,z), gamma in mu_g.     (1)

It preserves the equation because gamma^g=1, commutes with D, and
therefore induces a regular action on V. In the invariant rational slice

    s=v/[g*t^(g+1)*b(u)^k], E=C(u,s),

it is exactly (u,s)->(u,gamma*s). It scales Omega by gamma. It is NOT the
deck group D: D fixes E pointwise, whereas rho acts nontrivially on E.

On the selected source function field it has the expression

    t -> t, x -> gamma*x+(gamma-1)*sigma/a.           (2)

This has genuine poles at t=0 and at the b-roots. We never treat (2) as
a polynomial automorphism of the chosen source, nor apply affine-space
fixed-point theory to that rational source action.

## 2. The COMPLETE fixed schemes on V

For any nonidentity gamma in mu_g, put

    Gamma = Spec C[u,1/(u*b(u))].

It is a CLOSED smooth curve in V, obtained from v=0 in W. The equation
forces a!=0 there, because sigma is nonzero at every root of a. Taking
the free D quotient gives exactly the displayed punctured line.

There is also one distinguished point o on C0, represented by
t=0,z=0,v^g=sigma(0)^g. Then

    V^rho_gamma =
       Gamma disjoint-union C0,    if gamma=-1;
       Gamma disjoint-union {o},   if gamma^2!=1.    (3)

This is an exhaustion, not a partial list of fixed points. A point of V
is fixed exactly when rho_gamma(w)=d_epsilon(w) for some epsilon in D.
For t!=0, this forces epsilon=1, hence v=0 and gives Gamma.
For t=0, v!=0, so epsilon=gamma and (1-gamma^-2)*z=0.
This gives precisely the remaining sets in (3). Finite-order fixed
schemes of smooth complex varieties are smooth, so these pointwise
descriptions identify the reduced fixed SCHEMES as well.

Here are explicit tangent checks. Near C0, v is a unit, and the
D-invariant functions

    n=t/v, w=v^2*z

are local parameters at o. Their rho characters are gamma^-1,gamma^2.
Their differentials form a basis at o: on W, (t,z) are local parameters,
dn=dt/v and dw=v^2*dz there. Thus o is an ISOLATED fixed point for every
odd-prime-order gamma. For gamma=-1 the fixed tangent is the C0 direction.

Near Gamma the invariant local coordinates are u and r=v/t. The action
fixes u and sends r to gamma*r, so its fixed tangent has dimension ONE.

The mixed fixed dimensions in (3) are compatible with this source
geometry. V is NOT affine space; no Smith connectedness is asserted for V.

## 3. Target-field stability would give a regular target action

Suppose f:V->Y=A2 is an ACTUAL ETALE morphism, with target generators
P,Q in C. Its composition K=f*q is an actual polynomial plane Keller
map of degree

    d=g*b0, b0=[E:F]>=1, F=C(P,Q).                    (4)

It is noninvertible because g>=2. Rescale a target coordinate to have
Jac_(x,t)(P,Q)=1.

Assume a nontrivial subgroup H of the regular fibre group in (1)
preserves F AS A FIELD. The actual plane intersection theorem PF52 gives

    F intersect R=B=C[P,Q],
    F intersect C=B.                                 (5)

For every rho in H, rho(P),rho(Q) belong to C intersect F=B.
Applying rho^-1 gives rho(B)=B. Hence H acts by POLYNOMIAL automorphisms
of the ENTIRE target Y, and f is equivariant. This includes nonlinear
target actions; it does not assume homogeneous chosen generators.

At a fixed source point, equivariant etaleness identifies the invariant
tangent spaces. More precisely the restriction f^H:V^H->Y^H is etale,
by uniqueness of lifts along an etale map. This is the fixed-restriction
argument already retained in UASMITHFIX1.

The induced action of any prime-order subgroup is nontrivial. Otherwise
f*rho=f, and the etale local inverse at a fixed point makes rho the
identity on a nonempty open, hence on connected V. Alternatively, the
volume gives the exact target determinant

    Jac_Y(rho|B)=gamma !=1                           (6)

for rho=rho_gamma: f^*(dP wedge dQ)=Omega and rho^*Omega=gamma*Omega.

## 4. Odd-prime and involution exclusions

Choose a prime-order subgroup of H.

If its order p is ODD, (3) gives both a one-dimensional fixed component
Gamma and an isolated fixed point o. Etaleness forces Y^H to have local
dimension one at f(Gamma) and zero at f(o). But the ENTIRE fixed scheme
of a p-group on A2_C is smooth and connected, hence has constant dimension.
This is the already adopted Smith input, rechecked in
[Haution, Proposition3.2.1](https://arxiv.org/html/1507.04582#S3.SS2).
Contradiction. We use only target fixed-locus connectedness and smoothness;
there is no claim that V is acyclic or that a composite-order fixed locus
is automatically connected.

If p=2, (3) says that C0=A1 is fixed pointwise. Its image lies in the
smooth one-dimensional component T of Y^H. The target Smith theorem
also makes T the entire connected fixed scheme. It is A1: its smooth
compactification has genus zero and one point at infinity, by the
acyclic curve argument retained in UASMITHFIX1.

Consequently

    L0=A1 --q isomorphism--> C0=A1 --f etale--> T=A1

is a polynomial etale endomorphism of the affine line, hence has degree
ONE in characteristic zero. Thus K is injective on the WHOLE source
line L0. The retained Gwozdziewicz one-line theorem in PF52 forces K
to be an automorphism, contradicting (4).

The involution argument can even avoid target Smith: a smooth affine
curve dominated by A1 has genus zero and scalar units, hence is A1.
The nonconstant etale restriction C0->T then has degree one. The odd
prime argument still needs the retained connectedness input.

Every nontrivial finite subgroup has a prime-order subgroup. Therefore

    NO nonidentity rho_gamma preserves F.             (7)

In particular, there is NO etale plane target on this V equivariant
for ANY nontrivial subgroup of its regular fibre-scaling mu_g.
This is an actual all-degree exclusion of that target architecture,
not merely a necessary local numerical ledger.

## 5. Field escape and finite normalization

If a hypothetical target exists, its stabilizer under the regular
fibre group is TRIVIAL. There are exactly g distinct translated target
fields rho_gamma(F) inside E. For every gamma!=1,

    rho_gamma(P) notin F OR rho_gamma(Q) notin F.      (8)

Indeed membership of both would give rho_gamma(B) subset B and hence
rho_gamma(F) subset F. Its finite-order iterations force equality,
contradicting (7). This is a field obstruction, stronger than the failure
of the chosen outputs to be eigenfunctions.

It excludes every target pair that becomes fibre-scaling-equivariant
after ANY polynomial target rechart. It includes all pairs of fibre
eigenfunctions and every compatible polynomial target action, not just
a fixed diagonal ansatz. It does NOT imply that every target pair has
that symmetry. Averaging its components need not preserve its constant
Jacobian, so (7) cannot be promoted to an exclusion of all targets.

Retain the original full normalization S of B in L, its trace-dual Lambda
and EVERY ramified boundary prime. If the excluded field stability held,
(2) would extend to an automorphism of S by integral closure, and the
action on E would likewise extend to the FULL S_E. No boundary prime is
discarded to prove (7). Its contradiction already occurs on the actual
smooth V and target Y; S_E is not identified with V or C. The fibre action,
generic deck group and arbitrary monodromy remain separate.

PF142 still retires the mixed 2+3 and every non-power minimal equivariant
factor. PF143's full-S-marked quadratic integrality gate remains open.
This result closes the symmetry-preserving part of the surviving cyclic
target family. The remaining actual target must break EVERY nontrivial
regular fibre scaling; nonequivariant cyclic targets, original full
finite-lattice landing and full HC4/JC2 remain OPEN.

## 6. The same mechanism excludes equivariant targets on the classical pseudo planes

There is a useful general form of the involution argument. Let a smooth affine
V admit an etale plane map q of degree m>1, taking a whole source line
isomorphically to a curve C=A1. If a nontrivial prime-order automorphism of V
fixes C pointwise and preserves the field of an etale plane target, that
target is impossible. Equations (5) and the fixed restriction put C etale
over a smooth affine target fixed curve T. A1 domination forces T=A1,
and the one-line theorem contradicts the product degree. This version
needs NO Smith theorem.

Apply it to the ENTIRE classical family already retained in the
[acyclic surface proof](../universal-keller-acyclic-boundary/PROOF.md),

    S_g={u+u^2*v=w^g}, g>=2.

For plane coordinates (t,h), put

    u=t^g, w=t*(1+t^g*h),
    v=sum_(j=1)^g binom(g,j)*t^(g*(j-1))*h^j.        (9)

These polynomial formulas land on S_g. The smooth surface has the
global unit residue volume Omega=du wedge dw/u^2; near u=0 its regular
expression is -dv wedge dw/(1+2uv). Its pullback by (9) is
g*dt wedge dh. Thus (9) is ETALE on the ENTIRE plane, has generic
degree g, and maps t=0 isomorphically onto C={u=w=0}=A1_v.
The generic field degree follows by t^g=u and
h=(w/t-1)/u; t^g-u is irreducible over C(u,w).

The regular fibre action w->gamma*w fixes C pointwise. Therefore NO
plane target on S_g can have function field preserved by ANY nontrivial
subgroup of this action, in any target degree. This includes the
canonical-trivial exceptional pseudo-plane family; it is not an exclusion
of its unrestricted targets. Neither a torus action on the target nor
linearization of its finite action is assumed.

## 7. Sources and verification boundary

The source surface and full plane charts are MX0's retained constructions.
PF52 supplies the field/algebra intersection and its previously adopted
one-line theorem. UASMITHFIX1 supplies the fixed-restriction proof and
the primary Smith input. Haution's complete Proposition3.2.1 statement
and proof were reread; no new classification, finite-action linearization
or external theorem is imported.

The unbounded proof is (1)--(9). New exact controls check commuting
actions, local invariant coordinates, complete fixed-point equations,
character dimensions, the explicit quotient and selected source line.
They do not independently prove Smith theory or the one-line theorem.
No predecessor mathematical checker is rerun. Independent whole-proof
replay, formal verification, external review and full-goal completion
remain unclaimed.
