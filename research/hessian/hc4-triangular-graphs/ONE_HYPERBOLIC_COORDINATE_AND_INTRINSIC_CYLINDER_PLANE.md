# One hyperbolic coordinate and the intrinsic cylinder plane

15 November 2026. Claims **HC4MYJR0/1**. The source obstruction for
h0,h1 and the normalized three-point suspension extends beyond
automorphisms preserving a two-variable base. Keeping just one
hyperbolic coordinate algebra is already impossible under unit
Hessian, even when the other three source variables mix arbitrarily.

The new input is an intrinsic plane algebra inside the exceptional
three-variable trace. Its recovery from locally nilpotent derivations
forces a single constant scaling on the free cylinder coordinate.
The unsquared Keller orientation then contradicts a polynomial
having the same nonzero derivative at all its distinct roots.

Full source maps that change both hyperbolic coordinate algebras
remain open. The fixed K2/h1 equivalence is not proved or disproved.

## 1. Precise statements

For a polynomial B in C[z,w], define its relative invariant algebra

    ML_B(C[z,w,y])= intersection ker delta,

where delta ranges over all locally nilpotent C-derivations of
C[z,w,y] with delta(B)=0. This definition is used only through
the elementary proof below; no classification of additive actions
or stable equivalences is imported.

**JR0.** Suppose B has geometrically integral generic fibre whose
geometric coordinate ring is that of P1 with at least two points
removed. Then

    ML_B(C[z,w,y])=C[z,w].                              (1)

If two such plane polynomials B and B' become equal by an actual
three-variable source automorphism after adjoining the free
coordinate y or S, that automorphism identifies the two plane
subalgebras. The free coordinate has the form

    y=a S+b(t,R), a in C*,                              (2)

in the other cylinder coordinates. Thus its linear coefficient
is one constant on all critical components, not a separately
chosen unit on each component. This is the specific cylinder
rigidity needed here, not a new general stable-equivalence theorem.

There is also a first-normal-slope invariant. If the exceptional
trace of a potential is B(z,w) with y free and its first normal
derivative is a(z,w)y+L(z,w), then a source equivalence preserving
C[x] multiplies all slopes a(p), p in Crit(B), by one common
nonzero constant and permutes the points. This follows from (2)
and the critical chain rule. It does not allow independent scalings
on separate critical lines. The trace x=0 must correspond under
the equivalence; this is automatic in the K2/h1 application below.

**JR1.** Let B in C[z,w] have a finite Morse critical scheme with
N>=2 points and satisfy (1). Suppose a(z,w) takes one common
nonzero value A at every critical point of B. For arbitrary
polynomials L(z,w) and R(x,y,z,w), put

    h=B(z,w)+x(a(z,w)y+L(z,w))+x^2 R(x,y,z,w).           (3a)

No polynomial source automorphism Phi preserving C[x] setwise
can make h composed with Phi have nonzero constant full Hessian
determinant. All higher normal jets R and all source degrees are
unrestricted. In particular this applies to

    h(x,y,z,w)=2xy+B(z,w).

The same holds for preservation of C[y] in this hyperbolic case,
by symmetry.
Only one coordinate algebra is retained; all other source mixing
and all polynomial degrees are unrestricted.

In particular this excludes every C[x]- or C[y]-preserving
source map for

    h0=2xy+w(z^2-1),
    h1=2xy+w(z^2-1)+z.                                 (3)

More generally it covers 2xy+w C(z)+D(z) for every squarefree
nonconstant C with deg C=N>=2 and arbitrary D. The previously
proved stronger full-source obstructions at higher multiplier
degree remain retained; no old case is reopened.

For the exact [JQ](TWO_NORMAL_MODELS_AND_RELATIVE_RANK_ONE_ALIGNMENT.md)
suspension

    2PV+B_rho(u,v), B_rho=u(uv+rho)^2-u+uv+rho,

it excludes every source map preserving C[P] or C[V], separately,
for every rho in C, including rho=+-1. In the normalized cubic
coordinates these two actual polynomial coordinates are

    V=wy+rho-z,
    P=x+(wV-2w^2y-2w rho-1)/2.                          (4)

They are genuine coordinates by JQ's explicit global inverse.
These one-coordinate exclusions supplement JQ's exclusion of
all maps preserving the different natural base C[w,y].
The new assertion is made in the suspension representative:
the actual function is (2PV+B_rho) composed with a source map
retaining one of its actual affine coordinate algebras. For
K3_rho composed with Phi, write that same function using the
left-composed map J_rho composed with Phi, where J_rho is (4)
with u=w,v=y. Apply JR1 only when this composed map retains an
actual affine hyperbolic coordinate. Mere preservation by Phi
of a nonlinear pulled-back algebra is not the hypothesis.

Retain the constant jointly affine plane theorem
[HCTANG2](../tang-hc4-audit/STRUCTURAL_REVIEW.md), the critical-chain
rule in [HG0](GLOBAL_VALUE_ORBITS.md), and the actual affine-line
component theorem
[GR0](CRITICAL_COMPONENTS.md#1-a-plane-keller-fibre-with-an-affine-line-component).
The two projection ranks for a retained single coordinate are
the same constant forms already used in GR1. The new argument
closes their rank-zero case using (1), not a nonlinear distribution.

## 2. An elementary proof of the intrinsic plane algebra (JR0)

If delta is locally nilpotent on a characteristic-zero domain,
put ord_delta(a)=max{n:delta^n(a)!=0} for nonzero a. The map
exp(T delta) is a ring homomorphism to the polynomial ring in T.
Degrees add under products, so

    ord_delta(ab)=ord_delta(a)+ord_delta(b).

In particular every unit is in ker delta, since a unit and its
inverse multiply to one. Localizing by elements in ker delta
retains local nilpotence.

Let K=C(B) and extend a B-annihilating delta to the generic
fibre and then to an algebraic closure Kbar, acting trivially
on Kbar. By hypothesis its coefficient ring is

    Kbar[T,1/((T-a1)...(T-ar))], r>=1,

after placing one omitted completion point at infinity. In the
cylinder over this curve, T-a1 is a unit. Hence delta(T)=0,
and every geometric curve function is fixed. Faithful scalar
extension shows delta fixes z and w already before extension.
This proves C[z,w] is contained in every kernel.

Conversely partial_y is a B-annihilating locally nilpotent
derivation and its kernel is exactly C[z,w]. This proves (1).
There is no assumption that an arbitrary derivation is nilpotent.

An actual cylinder source equivalence conjugates precisely the
sets of B- and B'-annihilating locally nilpotent derivations,
so identifies their intersections. After the induced plane
isomorphism, it is an automorphism of a one-variable polynomial
algebra over a plane polynomial ring. Degree in the free variable,
applied also to the polynomial inverse, forces degree one.
Its leading coefficient is a unit of the plane ring, hence a
single nonzero constant. This proves (2).

For B=w C(z)+D(z), with C nonconstant, the generic fibre over
K=C(s) has the exact coordinate ring

    K[z,1/C(z)], w=(s-D(z))/C(z).

Indeed C and s-D are coprime over K[z], so C becomes a unit
on the actual generic level. Thus (1) holds in every multiplier
degree, with no requirement on the degree of D.

For B_rho, JQ supplies the actual generic parameter T=uv+rho
and its inverse, giving A1 minus {-1,1,s}. Over C(s) these
three points are distinct for every rho. Therefore (1) holds
for B_rho as well, including the critical-value collisions.

## 3. Retain one coordinate and apply the two constant Tang ranks

Let Phi preserve C[x]. It sends x to a x+b, a!=0. A constant
affine change of the actual source normalizes this to Phi^*x=x.
It preserves the property of unit Hessian and gradient inverse.
Write g=h composed with the normalized Phi.

At x=0 the native full critical equations give Crit(B) and
y=-L(p)/A. Their full Hessian determinants are
-A^2 det Hess B(p), hence nonzero. There are at least N>=2 full
critical points, regardless of additional points away from x=0.
Thus a hypothetical unit-Hessian g cannot have an automorphic
gradient. Tang gives
a constant jointly affine two-plane W. Its projection to x
has rank one or zero. Only constant changes preserving x are
used in passing to the following normal forms.

The critical scheme of the native exceptional trace h|x=0 is

    Crit(B) times A1_y.

It is exactly N disjoint affine lines, as a reduced scheme.
Each is also a closed irreducible component of the entire
three-relative critical scheme C_h=V(h_y,h_z,h_w). Indeed
h_y=x(a+xR_y), and the factor in parentheses is the nonzero A
along each line. Locally the relative ideal is therefore
(x,B_z,B_w), a smooth reduced line. Additional relative components
away from x=0 are allowed. Because Phi preserves x, the vertical
first chain rule identifies the entire corresponding C_g with
C_h and identifies their exceptional traces. In particular C_g
has actual closed A1 components, not merely normalizations.

### Projection rank one

The constant normal form is

    g=D(z1,z2)+x F(z1,z2)+Q G(z1,z2), J(F,G)=c in C*.

The three-relative equations are G=0 and
grad D+x grad F+Q grad G=0. Their coefficient determinant
is c, so x and Q are polynomial functions of z1,z2. Therefore
the whole relative critical curve is isomorphic to G=0.
It has an actual closed affine-line component. GR0 makes
(F,G) invertible, and cotangent reconstruction makes grad g
invertible. This contradicts N>=2. This is the retained GR1
rank-one argument, with no new affine-line rectification input.

### Projection rank zero

Use base (x,t) and free variables (P,Q):

    g=A(x,t)+P F(x,t)+Q G(x,t), J_(x,t)(F,G)=c in C*.

At x=0 write F0(t)=C(t)a(t),G0(t)=C(t)b(t), with a,b coprime
and C their greatest common divisor. C is nonconstant and has
exactly N distinct simple roots. To justify every part: the
critical scheme of the trace g0=g|x=0 is precisely C_g intersect
{x=0}, so is
N disjoint affine lines. Its last two equations are F0=G0=0;
at a common root J(F,G)!=0 makes (F_t,G_t) nonzero. Hence
the common multiplicity is one, and the remaining t derivative
cuts exactly one affine line in (P,Q). An identically zero
pair or no common root is incompatible with this scheme.

Choose polynomials a1,b1 with a b1-b a1=1, and put

    R=aP+bQ, S=a1P+b1Q,
    P=b1R-bS, Q=-a1R+aS.

This is an actual source automorphism with determinant one,
independent of x. The whole exceptional trace becomes

    g0=A0(t)+C(t)R,

with S free. Its generic plane curve is A1 minus the N roots
of C, by the exact solved R formula. Thus JR0 recovers the
plane subalgebra C[t,R] intrinsically on this trace.

## 4. The free cylinder scale forces unsquared orientation

The restriction of Phi to x=0 is an actual three-variable
source equivalence between B(z,w) with y free and
A0(t)+C(t)R with S free. Their intrinsic plane algebras agree
by JR0. Consequently the old free coordinate has the form

    Phi^*y|x=0=d S+E(t,R), d in C*.                    (5)

On Crit(g0), the vertical derivatives of h vanish after Phi.
Differentiating g=h composed with Phi in x therefore gives

    g_x|Crit(g0)=Phi^*(a(z,w)y+L(z,w))|Crit(g0).

Terms differentiating the other three coordinates have coefficients
h_y,h_z,h_w and vanish on that entire scheme. Since Phi^*x=x,
there is no nonconstant normal factor. Since a=A on every
native critical point, equation (5) says the coefficient of S
in g_x is the same nonzero A d on every line. All terms x^2 R
have disappeared from this identity; no degree bound on R is used.

On the line over a root t_i of C, the constant Tang form gives
the coefficient directly:

    coeff_S(g_x)=-b(t_i)F_x(0,t_i)+a(t_i)G_x(0,t_i)
                 =-c/C'(t_i).                         (6)

Here F_t=C'a,G_t=C'b at the root, and the displayed sign
uses J_(x,t)=F_xG_t-F_tG_x. The derivative of the (R,S) change
adds no x term, since that change was independent of x.

Thus C'(t_i)=-c/(A d), the same nonzero constant at all N roots.
The polynomial C'+c/(A d) has degree at most N-1 and N distinct
zeros. It vanishes identically, forcing C' constant and N=1.
This contradicts N>=2 and closes rank zero. JR1 follows.

This uses the orientation itself, rather than its square in
the critical Hessian determinant. Crucially the constant free
scale comes from the intrinsic plane algebra, even though Phi
was allowed to mix all three nonretained variables. A separate
choice of scale at each critical line would destroy the argument.

## 5. Applications and exact remaining boundary

For the family (3) and every squarefree C of degree N>=2,
the plane critical points are C(z)=0,w=-D'(z)/C'(z), with
det Hess B=-C'(z)^2!=0. JR0 and JR1 apply in all degrees.
In particular every source map preserving just C[x] or just
C[y] is excluded for h0,h1. This is stronger than merely
preserving their two-variable base algebra, but does not
classify arbitrary source maps.

These exclusions also retain arbitrary x^2 R perturbations of
the potential, and arbitrary first-normal terms xL(z,w), under
the stated constant-slope condition. Higher normal jets do not
reopen the obstruction. In the general perturbed family only
C[x] preservation is asserted; the C[y] statement uses the
symmetric hyperbolic representative.

For B_rho, JQ gives exactly three distinct Morse points for
every rho. JR1 excludes preservation of C[P] or C[V] separately
in the exact suspension coordinates, including rho=+-1. In
the original normalized cubic model use (4) only to express
the very same source function via its left-composed source map,
not to conjugate away a nonlinear Hessian change. JQ's
exclusion of maps preserving C[w,y] remains independently valid.
None of these statements requires bounding rho or a source degree.

N>=2 is essential to the nonexistence conclusion: for a nondegenerate
quadratic B there is one critical point, and 2xy+B itself has unit
Hessian and automorphic gradient. Generic-fibre rigidity is also
a stated hypothesis, not a conclusion for every plane polynomial.
For B=z the generic plane fibre is A1 and the cylinder derivation
partial_w prevents recovery of the full plane algebra in (1).

There is also a direct new limit on the K2/h1 equivalence question.
No source equivalence between them, even up to nonzero scalar and
constant, can preserve C[x] or C[y]. For C[x], both full three-relative
gradient ideals equal (x,z^2-1,2wz+1). Their critical schemes are
exactly the two lines x=0,z=+-1,w=-1/(2z), with y free.
Thus any retained-x equivalence must carry the trace
x=0 to itself, so its affine x translation is zero. JR0 forces
one common free-y scale on its two critical lines. The first
normal derivative for h1 has slope 2 on both lines, whereas
for K2 it has slopes 2w=-1,1. These two unordered slope pairs
cannot differ by one nonzero common multiplier. Source and value
scalars do not change that distinction. The C[y] case is symmetric.
This does not exclude a global equivalence changing both algebras.
Nor does JR1 exclude a one-coordinate-preserving unit-Hessian
route for K2 itself: its unequal normal slopes fail JR1's stated
constant-slope hypothesis. JQ's natural-base exclusion for K2
remains valid, and its remaining one-coordinate route is open.

The proof uses a constant Tang plane only in the actual unit-Hessian
source coordinates. It does not apply that rank argument to a
nonlinear pulled-back distribution. The source algebra C[x] is
assumed preserved; it is not recovered intrinsically from the full
four-variable potential. In particular neither full-source h0,h1
nor full-source K3_rho is closed here. The K2/h1 equivalence,
both old endpoints, actual global embedding recovery and unrestricted
HC4/JC2 remain open.

For primary context we inspected
[Dubouloz--Heden--Kishimoto v1](https://arxiv.org/abs/1707.08768v1)
on additive actions and punctured-surface torsors. We do not import
their classification, assert properness of a new action or identify
an arbitrary quotient. We rechecked the retained affine-line
rectification context in
[Jelonek v3, Theorem 2.4](https://arxiv.org/html/2011.03472v3#S2.Thmdefi4).
The needed affine-line theorem is already retained in GR0.
All new unit, kernel and cylinder arguments are proved above.
The finite checker controls these written symbolic identities and
scope boundaries, not an independent or formal proof audit.
