# Exact collision action and first-derivative reconstruction

Date: 2026-11-16. Claims HC4MYLC0 and HC4MYLC1.

FULL THEOREM TARGET: prove HC4/JC2 or construct an actual counterexample.
CURRENT LAST UNSOLVED IMPLICATION: the ACTUAL polynomial plane source forces
Tr(Lambda_L*Lambda_L) subset B, equivalently regular inverse Gram entries;
the alternative full marked-energy integrality implication is also OPEN.
Neither implication is proved here. External review/outreach remains DEFERRED.

We strengthen the ACTUAL compatibility in the
[LB collision construction](GLOBAL_COLLISION_CHORDS_AND_PENCIL_DEFECT_DROP.md).
A canonical exact action generates its ordered chord cover, and its square
generates the unordered cover. First chord derivatives recover the actual
source points, including as regular functions across chord ramification.
Primitive action orders can nevertheless be singular; no order is silently
substituted for its normalization.

## 1. A global polynomial primitive, without an H1 assumption

Let K=(P,Q):A^2->A^2 be complex Keller, j=Jac(P,Q)!=0 constant.
Use source coordinates (x,y). The polynomial one-form

    P*dQ-j*x*dy

is closed because dP wedge dQ=j*dx wedge dy. It has a POLYNOMIAL
primitive S, unique up to a constant:

    dS=P*dQ-j*x*dy.                                      (1)

This needs no cohomology hypothesis on the collision surface. Directly,
write the form A dx+B dy; integrate A in x and then integrate the residual
B-S_y, which is independent of x, in y. Characteristic zero permits the
polynomial antiderivatives.

Retain the entire ordered collision surface C and its swap quotient D
from LB0. For p=(x1,y1),q=(x2,y2) in C put

    u=(p+q)/2, v=p-q,
    A(p,q)=(S(q)-S(p))/j,
    F(p,q)=A(p,q)-u_x*v_y.                               (2)

These are regular functions on the ACTUAL affine C; F is polynomial in
the four source-point coordinates restricted to C. The additive constant
of S cancels. Since K(p)=K(q), pulling back (1) gives

    dA=x1*dy1-x2*dy2,
    dF=v_x*du_y-v_y*du_x.                                (3)

The swap sends F to -F. Thus G=F^2 is a global regular function on D.
There is also a direct chord interpretation:

    F=-(1/j)*integral_(straight segment q->p) P*dQ.         (4)

The source integral of x dy on the straight segment is u_x*v_y, proving
(4) from (1). This is an integral of a polynomial in the segment parameter,
so (4) is algebraic and independent of a choice of analytic continuation.
Changing the target primitive by an exact form dh(P,Q) does not change
the endpoint action, since the two target endpoints coincide. No real
positivity, nonzero area at EVERY finite pair, or nonlinear source-coordinate
invariance is inferred.

## 2. LC0: exact chord derivatives and universal primitive elements

On the nonvertical chord chart write

    y=m*x+b, r=v_x!=0, u=u_x.

Its pullback is exactly C_r={v_x!=0}, an affine open of C. The chord
coordinates are regular there and u_y=b+m*u, v_y=m*r. Formula(3) becomes

    dF=r*db+r*u*dm.                                     (5)

The chord map is generically etale by LB0. The rational derivations
partial_b,partial_m on the finite separable component fields over
k0=C(m,b) therefore give

    F_b=r, F_m=r*u.                                     (6)

In particular BOTH first derivatives are regular on the ENTIRE C_r,
including chord ramification points, and F_b is a UNIT there. This does
not say that these rational derivations preserve the entire coordinate
ring or that second derivatives are regular. From (6) one reconstructs

    r=F_b, u=F_m/F_b,
    x1=u+r/2, x2=u-r/2,
    y1=b+m*x1, y2=b+m*x2.                               (7)

There is no missing vertical direction. On the other chart x=n*y+c put
s=v_y!=0. Formula(3) gives

    F_c=-s, F_n=-s*u_y,
    u_y=F_n/F_c, s=-F_c,
    y1=u_y+s/2, y2=u_y-s/2,
    x1=c+n*y1, x2=c+n*y2.                               (8)

These charts cover the entire ordered collision surface and use the SAME
global F. No branch, ramification point, or action-zero locus was deleted
from this ORDERED first-derivative reconstruction.

Let E_C and E_D be the finite etale GENERIC PRODUCT algebras of the
ordered and unordered chord maps over k0. Their dimensions are 2N and N,
where N is LB0's total generic chord degree. No irreducibility or generic
deck-group hypothesis is needed. Then

    E_C=k0[F],  E_D=k0[G].                               (9)

Here is a proof that also handles product algebras and separates their
components. In a finite separable product extension, each base derivation
extends uniquely. For ANY element z, the square-free minimal polynomial
M_z gives

    partial_i z=-(partial_i M_z)(z)/M_z'(z) in k0[z],      (10)

because M_z'(z) is a unit in the finite reduced algebra k0[z]. Thus
k0[F] contains F_b,F_m. Equation(6) makes F_b nonzero in EVERY component;
it is therefore a unit of k0[F]. Formulas(7) recover all source-point
coordinates, which generate E_C over k0, proving the first equality.
It also shows F is nonconstant and not identically zero on ANY component.

For G, its derivatives satisfy

    G_b=2*F*r, G_m=G_b*u,
    u=G_m/G_b, r^2=G_b^2/(4*G).                          (11)

The elements G and G_b are nonzero in every generic component, hence
units in k0[G]. The invariant source-point coordinates are the midpoint
u_x,u_y and the three quadratic products v_x^2,v_x*v_y,v_y^2. They are
recovered from (11), u_y=b+m*u_x and v_y=m*v_x. They generate E_D,
proving the second equality. This does not assume that a primitive element
of each individual field automatically separates a product: reconstruction
proves that separation for the WHOLE product algebra.

Consequently the square-free monic minimal polynomial H(V;m,b) of G has
degree EXACTLY N, and that of F is

    Phi(T;m,b)=H(T^2;m,b), degree EXACTLY 2N.               (12)

Indeed swapping the two endpoints pairs F and -F; the reconstruction
separates all ordered generic points. None has F=0 generically. Thus the
even polynomial in (12) is square-free and already has the full generic
degree. No global monicity or integrality over C[m,b] is asserted: these
are polynomials over its FIELD k0. Clearing denominators can produce a
nonmonic primitive action hypersurface.

The derivative argument (10) is standard, not a new primitive-element
theorem. A particularly relevant primary precedent is Adamus--Hajto,
arXiv:2610.36736v1, Theorem4.15, printed/PDFp6: the Legendre potential of
a generically finite gradient map is primitive because its Nambu derivatives
recover the source coordinates. Its proof, Theorem4.17 pp6-7 and the
cotangent discussion pp7-8 were read. That paper explicitly leaves HC4/JC2
open. We use it as RELATED PRIMARY CONTEXT, not an imported theorem or a
new solution. The collision action, whole product-algebra separation and
the following exact ring recoveries are proved directly here. No priority
claim is made for the derivative mechanism.

## 3. LC1: complete chart rings from first derivatives

Write A_C=O(C_r). Formula(7) is a global coordinate-ring statement:

    A_C=C[m,b,F,F_b,F_m,F_b^-1].                          (13)

All right-side elements are REGULAR on C_r by (6), and F_b^-1=r^-1.
Conversely C_r is the ACTUAL closed ordered collision locus in the source
four-space localized at r. Its ring is generated by x1,y1,x2,y2,r^-1,
and (7) expresses all these in the displayed ring. This proves equality,
including finite special points and ramification. Equation(8) gives the
same statement in the vertical chart, with F_c invertible.

On the unordered chart D_r=C_r/S2 one has

    A_D=O(D_r)=C[m,b,u,r^2,r^-2],                         (14)

as a quotient by its actual relations: the free swap negates r and fixes
m,b,u. The invariant square G and its first derivatives are regular on
ALL of D_r, but the ratios in (11) need not be regular expressions at
G=0. On the precise open D_r minus {G=0}, they give

    A_D[G^-1]=C[m,b,G,G_b,G_m,G^-1,G_b^-1].               (15)

Indeed G_b^2=4*G*r^2 makes G_b a unit after inverting G, and (11) recovers
the generators of (14) and their necessary inverses. Action-zero points
are NOT covered by (15); the global ordered formula (13) still covers
their preimages. No unproved removal of that locus is used.

These are identities in the coordinate rings of the normal collision
surfaces, not proofs that the SMALL orders C[m,b,F] or C[m,b,G] are normal
or smooth. For a polynomial primitive relation Psi(m,b,T)=0 obtained by
clearing denominators in (12), differentiation gives on C_r

    Psi_T*F_b+Psi_b=0, Psi_T*F_m+Psi_m=0.                 (16)

These relations must not be read as globally inverting Psi_T. They can
give regular quotients even at Psi_T=0.

In fact EVERY ramification point of the ordered chord map maps to a
SINGULAR point of its primitive action hypersurface:

    chord ramification ==> Psi=Psi_T=Psi_m=Psi_b=0.       (17)

If Psi_T were nonzero at such an image, the hypersurface would be smooth
and etale over the chord base there. By (9) the map from C_r to its
hypersurface is birational componentwise; it is quasi-finite since the
chord map is. Zariski Main over that normal smooth neighbourhood makes
it an open immersion, contradicting chord ramification. Hence Psi_T=0,
and (16), whose quotients are genuinely regular on C_r, forces both
remaining partials zero. This also works for the reduced whole product
relation: Psi_T!=0 isolates one smooth hypersurface component. The converse
to (17) is NOT asserted; unramified finite branches can collide in an order.

For orientation, suppose N=2 and D is CONNECTED. Its chord map must
ramify. Otherwise, in its normal finite degree-two completion, a branch
value would have one local fibre factor of length2. Surjectivity forces
that point to be retained, contradicting etaleness of the open source.
Thus the completion would be finite etale over L_line. But L_line is
the total space of O_(P1)(1), hence simply connected, so that completion
would split into two copies. Its connected open D could not have total
generic degree2. Therefore (17) forces finite singularities of the
ordered primitive action hypersurface in this minimal connected case.
It does NOT exclude that case: the full collision surface is smooth and
the first-derivative ring (13) can resolve precisely those order singularities.

The original target algebra can also be recovered as an EXACT INTERSECTION
inside the coordinate ring O(C) of the WHOLE ordered collision surface:

    p^*C[x,y] intersect q^*C[x,y] = C[P,Q].               (20)

The common ring on the right is embedded by the common K-value. A
noninvertible plane Keller map has d>=3: degree2 would be Galois and is
already excluded by KMFI1. Suppose f(p)=g(q) on C. After passing to a
generic splitting field, the equality says f_i=g_j for EVERY ordered pair
of DISTINCT generic sheets i,j. Given j!=k choose a third sheet i; then
g_j=f_i=g_k. Thus all conjugates of g are equal, g belongs to the target
field, and the same is true of f with f=g. The retained PF52 identity
F_target intersect C[x,y]=C[P,Q] proves (20) as a ring equality, not just
a field equality. Conversely every target polynomial is plainly common.

ALL ordered collision components are required for this argument unless
pair-transitivity already makes C irreducible. The claim is not made for
one selected orbital component. The d>=3 hypothesis is essential: in a
quadratic extension the sole off-diagonal component can identify the two
entire source fields, and the two polynomial subrings can have a much
larger intersection than the target. Equation(20) supplies another exact
actual-target compatibility condition for the first-derivative reconstruction;
it does not prove that an abstract action function produces a polynomial
intersection ring or generators with unit Jacobian.

## 4. Exact nontrivial controls and what they do not prove

The sharp pseudo-plane family of LA1 supplies a genuine, nonvacuous
exact-source control. In its k=2 case the ACTUAL etale surjective map
from A2_(x,t) to the pseudo-plane T2 has generic deck formula on x!=0

    sigma(x,t)=(-x,-t-2/x^2).

This is NOT an etale endomorphism of an affine-plane TARGET. It comes from
the free quotient of x^2*y=z^2-1 and the chart z=1+x^2*t,
y=2*t+x^2*t^2, already credited to Dubouloz--Palka in LA1.
The descended exact primitive from LA1 pulls back with polynomial S

    S=-t*x*(t^2*x^4+3*t*x^2+4).

For p=(x,t),q=sigma(p), formulas(2) give

    F=-4/x, b=-1/x^2, m=(t+1/x^2)/x,
    F_b=2*x, F_m=0, G=-16*b.                             (18)

Equations(7) recover BOTH actual points. The ordered action relation is
F^2+16*b=0; its squared action is primitive of unordered degree1. This
does not contradict LB1: the chord map misses, for example, the source
line x=0, which is injective for this pseudo-plane-target map. Thus
surjectivity onto ALL affine source lines, specific to the plane Keller
one-line theorem, is essential. The control demonstrates that canonical
primitive action and first-derivative reconstruction ALONE do not close
the required target-plane compatibility.

Exactness in (3) cannot be inferred merely from the equality of symplectic
two-forms. On Gm times A1, the rational symplectic shear
(x,y)->(x,y+1/x) changes x dy by -dx/x. This closed form has no regular
algebraic primitive there: differentiating a Laurent polynomial has no
x^-1 coefficient. Our proof supplies the actual polynomial S on the
WHOLE original A2 and uses equality of K-values; it does not assume H1(C)=0.

Regular UNIT first derivatives do not imply regular second derivatives.
On b=z^2 take F=z^2+z^3, independent of m. Then

    F_b=1+(3/2)z, F_m=0, F_bb=3/(4z),
    (F-b)^2-b^3=0.                                      (19)

Near z=0, the first derivative is a unit and reconstructs z regularly,
while the second derivative has a pole. The primitive plane equation
is singular there: all three derivatives in m,b,F vanish. The regular
ratio -Psi_b/Psi_T equals 1+(3/2)z although both numerator and denominator
vanish. Reconstruction(7) gives two distinct source-point maps with the
SAME nonzero Jacobian -3/8 at z=0, in coordinates(m,z). Thus this local
model retains unit first derivative, chord ramification and local etale
source-point maps; it still supplies NO whole polynomial Keller K with
their common value. It verifies the necessity of distinguishing the
normal first-derivative ring from its singular primitive order.
Likewise primitive orders
may have finite collisions even at unramified points, as already proved
by KW1; no old order checker or whole-proof replay is needed here.

## 5. Remaining gap

LC0/1 replaces an abstract two-cover compatibility proposal by a specific
global exact action, primitive generic equations and COMPLETE ordered
first-derivative chart rings. It does not make Phi/H monic over the chord
coordinate ring, make all chord derivatives regular, make primitive orders
normal, force chord degree N=1, or exclude the connected N=2 case. It
does not turn this action into the marked HC4 energy or make its traces
regular over the original Keller target. Those are different covers and
different functions. The universal trace-dual/marked-energy gap remains OPEN.

The checker supports the exact differential identities, reconstruction,
nontrivial pseudo-plane control, source locks and scope. It does NOT prove
the global algebraic-geometric arguments or decide HC4. Independent replay,
formal verification, external review and goal completion remain false.
