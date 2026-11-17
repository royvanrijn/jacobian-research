# Exact quotient holes, torsor residues and the local period barrier

15 November 2026. Claims **HC4MYJU0/1/2**. The necessary action now
has an exact plinth ideal and an exact quotient-image defect. Its free
part carries a torsor whose annihilator and point residues are determined
in every degree. These are additional necessary global tests, not an
action classification or a recovered source coordinate system.

A separate possible obstruction, the complete local volume period at a
Morse point, is investigated here and exhausted for the remaining models.
Every remaining normalized model passes it to all orders. Computing more
coefficients of these local periods cannot close the full-source gap.

## 1. Statements and retained inputs

Work over C. Retain [JT2](PRIMITIVE_QUOTIENTS_PLINTH_AND_GLOBAL_ACTION_CERTIFICATES.md):
if h has N>=2 finite Morse critical points and its source orbit contains a
unit-Hessian member, there are actual polynomial source coordinates
(b1,b2,P,Q), their polynomial inverse, and

    h=d(b)+c H, H=P F(b)+Q G(b), c in C*, J_b(F,G)=j in C*.
    D=G partial_P-F partial_Q, A=ker D=B[H]=B[h], B=C[b1,b2].

Put I=(F,G) in B and Z=V(I)={b_1,...,b_N}. Each zero is reduced,
since j is a unit. The unique critical point above b_i has value
c_i=d(b_i). The old notation b_i for a point is distinct from b1,b2
for coordinate functions. No nonautomorphic Keller pair is constructed.

**JU0.** For every such Keller translation the exact plinth ideal is

    pl(D)=I A.                                             (1)

It is a radical complete intersection, and A/pl(D)=product_i C[h].
In quotient coordinates (b,h), its zero set consists of N disjoint lines
L_i={b_i} times A1_h. The quotient morphism has the exact set-theoretic
image

    (A3 minus union_i L_i) union {(b_i,c_i):1<=i<=N}.        (2)

Its fixed ideal is I R. Its contraction to A is

    I R intersect A=(F,G,h-d(b)),                          (3)

so each fixed plane is contracted to the single indicated point on its
own plinth line. All other points of that line are absent from the image.
This proves the exact formula deliberately left open by JT0; the old
proof, checker and receipt remain unchanged.

**JU1.** Set U=A2_b minus Z. On every value h=kappa the free part above
U is a Ga torsor with transition function

    (kappa-d)/(c F G)                                    (4)

on the cover D(F),D(G), with the sign fixed by the slice ordering below.
In H^1(U,O_U)=B_FG/(B_F+B_G), its class has exact annihilator

    ann_B[(kappa-d)/(cFG)]=I:(kappa-d)
                         =intersection_{c_i!=kappa} m_(b_i). (5)

The empty intersection means B. In particular it is I at every
noncritical value. Its local Grothendieck residues, measured against
db1 wedge db2, are exactly

    (kappa-c_i)/(c j).                                   (6)

The class belongs to the I-annihilated socle: no higher pole orders are
present in these actual coefficient coordinates. The same global j
occurs at every point. Abstract isomorphisms of quotient algebras do not
retain these embeddings, residues or the source coordinate inverse.

**JU2.** A necessary local volume-period test for an N>=2 unit-Hessian
source orbit is that at each Morse point the vanishing three-cycle
period of the Gelfand--Leray form Omega/dh is exactly a nonzero constant
times h-c_i, with no higher powers. Polynomial source changes multiply
all these periods by one common nonzero constant, with independent
orientation signs for the chosen cycles.

The coordinate endpoint has no critical point. The quadratic endpoint,
h0=2xy+w(z^2-1), h1=h0+z, fixed K2=w(2xy+z^2-1)+z, and every member
K3_rho=K2+2x(rho-z), including rho=+-1, all pass this local period test
exactly, not just to a finite order. This does not assert completeness
of period invariants for analytic equivalence, simultaneous analytic
equivalence of different points, polynomial gluing, a new action, or
unit Hessian in any actual global source coordinates.

Retained sources are [JS](GLOBAL_ADDITIVE_ACTIONS_AND_THE_K2_TRANSLATION_BARRIER.md),
[JQ](TWO_NORMAL_MODELS_AND_RELATIVE_RANK_ONE_ALIGNMENT.md),
[HG](GLOBAL_VALUE_ORBITS.md),
[JK](SUSPENDED_MULTIPLIER_ORBITS_AND_BASE_PRESERVING_CLOSURE.md), and
[Tang](../tang-hc4-audit/STRUCTURAL_REVIEW.md). All new action, residue
and period arguments are given below. The primary
[Francoise article record](https://www.numdam.org/item/AST_1978__59-60__119_0/)
and [Stacks etale-presentation lemma](https://stacks.math.columbia.edu/tag/00U9)
were inspected as context. No uninspected isochore classification or
flatness theorem is required in the proofs here.

## 2. An elementary highest-weight calculation proves the exact plinth

First use independent variables U,V,P,Q and derivations

    D0=V partial_P-U partial_Q,
    E0=P partial_V-Q partial_U,
    W=U partial_U+V partial_V-P partial_P-Q partial_Q.

Then [D0,E0]=W, [W,E0]=-2E0, and both D0,E0 kill
H0=UP+VQ. A monomial has weight deg_(U,V)-deg_(P,Q).
D0 raises weight by two. E0 is locally nilpotent, since it differentiates
only U,V and fixes P,Q.

Suppose D0 r=p(H0). Its weight-zero part is D0 r_(-2)=p(H0).
Write r_(-2)=r'. Since E0 p(H0)=0 and W r'=-2r', induction gives

    D0 E0^m r'=-m(m+1) E0^(m-1) r', m>=1.               (7)

If r' were nonzero, apply (7) with m one greater than the last nonzero
E0 iterate of r'. The left side is zero and the right side is nonzero.
Thus p=0. In particular no nonzero polynomial in H0 alone is a D0 image.

Return to B and a=sum_k a_k(b)H^k in pl(D). At any b_i, the completed
local ring of B is C[[F,G]]: the invertible linear Jacobian constructs
the formal inverse recursively. This is a local formal chart only.
Write a=D r in that ring, keeping r polynomial in P,Q. Extract weight
zero on the right and weight minus two from r. These are finite
polynomials in F,G,P,Q: the P,Q degree is bounded, and a fixed weight
then bounds the F,G degree. The preceding calculation applies unchanged
and gives

    sum_k a_k(b_i)H^k=0.

Consequently every a_k vanishes at every point of Z. The ideal I is
radical: its zero-dimensional local quotients have independent linear
parameters F,G and hence length one. Therefore a_k is in I. This proves
pl(D) is contained in IA. The reverse containment follows immediately
from D(P H^k)=G H^k and D(-Q H^k)=F H^k. This proves (1), without
transferring a three-variable principal-plinth theorem to dimension four.

The kernel equality A=B[H] is retained from JS/JT. Formula (1) and the
Chinese remainder theorem give A/IA=product_i C[H]. On the source fixed
ideal IR, H vanishes; conversely the map A -> R/IR=product_i C[P,Q]
has kernel (I,H). Replacing H by (h-d)/c gives (3). If b is outside Z,
the nonzero linear form PF+QG assumes every value. If b=b_i, it is zero
and h=c_i for every P,Q. This proves (2) and its exact fibre description:
an affine line off the plinth locus, an affine plane at each image
point on it, and an empty fibre at every other point of it.

These statements concern the quotient by the necessary action, not the
image of the original plane Keller map. Its nonproperness is not assumed
away, and quotient surjectivity is not asserted.

## 3. The torsor class has an exact annihilator and affine residues

On D(G) use the slice t_G=P/G. On D(F) use t_F=-Q/F.
Both satisfy D(t)=1, and t_G-t_F=H/(FG). For h=kappa this is (4).
These charts and their inverse formulas are polynomial after the stated
localizations, so this is an actual Zariski Ga torsor on U.

Here is an elementary useful injection, including its kernel proof:

    B/I -> B_FG/(B_F+B_G), [a] |-> [a/(FG)].              (8)

It is well defined. If a/(FG)=A/F^m+B'/G^n, take m,n>=1 and clear
denominators:

    a F^(m-1)G^(n-1) in (F^m,G^n).

For the coprime pair F,G in the UFD B,

    (F^m,G^n):(F^(m-1)G^(n-1))=(F,G).                  (9)

For completeness reduce a putative relation modulo F^m. Coprimality
makes G a nonzero divisor modulo every F^m, by unique factorization,
so cancel G^(n-1), obtaining a F^(m-1)=G C modulo F^m.
Reduction modulo F^(m-1) makes C divisible by F^(m-1). Cancel that
power in the resulting equality to obtain a in (F,G). This proves
(9) and injectivity of (8). The same argument handles a unit in the
pair, when both sides of (8) are zero.

The class (4) is the image of (kappa-d)/c under (8). Multiplication
by a kills it exactly when a(kappa-d) is in I. Since B/I is the product
of its N residue fields, this gives precisely (5), including repeated
critical values. The double Cauchy residue in local parameters F,G is
the numerator at their common zero; db1 wedge db2=j^-1 dF wedge dG.
This proves (6). The factor j is globally constant, not a separately
chosen rescaling at each missing point.

This gives a further necessary global obstruction: a proposed quotient
whose plinth is nonreduced, has the wrong missing-line image, or whose
torsor class fails these annihilator and residue requirements cannot be
the JT2 action. Passing all tests still does not supply b1,b2 as part of
a complete source coordinate system or establish the coefficient Keller
condition. No classification of all actions fixing h follows.

## 4. Why local volume periods are necessary

Let Omega be the native constant four-form. On a smooth local level of
h, Omega/dh denotes a three-form eta satisfying dh wedge eta=Omega;
its restriction to the level and its period are independent of choices.
A polynomial source automorphism has constant nonzero Jacobian, so
pullback changes Omega by a single constant. The period condition in
JU2 is therefore an actual full-source invariant up to that constant.

For the required Tang form, use u=F,v=G near b_i. This analytic base
change has constant Jacobian j. Write
d(b(u,v))-c_i=u alpha(u,v)+v beta(u,v), and absorb alpha,beta in P,Q.
After also accounting for c, the function becomes exactly

    h-c_i=u P'+v Q'

and the volume form becomes a constant multiple of du dv dP' dQ'.
These are local analytic coordinates, not a polynomial global inverse.
For a nondegenerate quadratic in four variables, scaling coordinates
by sqrt(epsilon) makes Omega/dh scale by epsilon. Its vanishing-cycle
period is thus a nonzero constant times epsilon, with no higher terms.
This proves necessity. The period argument does not claim that every
unit-Hessian polynomial in every dimension has this property: the
retained four-variable Tang reduction and N>=2 hypothesis are used.

## 5. The complete remaining local periods pass

For h=2xy+w C(z)+D(z) at a simple root a of C, put c_i=D(a).
For its plane part, the period of dw dz/dh around the local vanishing
circle is, up to orientation, the contour integral

    integral dz/C(z)=2 pi i / C'(a).                    (10)

It is independent of the nearby value. Indeed solving for w on the
level gives the Gelfand--Leray form dz/C(z); a small contour surrounding
a is the vanishing annulus generator. The plane node is Morse since
C'(a)!=0. Adding the hyperbolic pair 2xy suspends this circle: its
three-cycle period is a fixed nonzero factor times the integral of
the plane period from c_i to the level. This assertion follows directly
by slicing the suspended cycle along 2xy and applying Fubini; the
hyperbolic annulus has constant residue period. Thus (10) gives an
exactly linear four-dimensional period. This covers h0,h1, in fact
every simple-root suspended multiplier in arbitrary degree.

For K2 the critical points are x=y=0,z=+-1,w=-1/(2z), so w0^2=1/4.
Locally w is a unit. The translation Z=z+1/(2w) has Jacobian one and

    K2=w(2xy+Z^2)+f(w), f(w)=-w-1/(4w),
    f(w)-f(w0)=-(w-w0)^2/w.

Choose a local square root r=sqrt(w) and put

    X=r x, Y=r y, T=r Z, t=(w-w0)/r.

Then K2-f(w0)=2XY+T^2-t^2 and the volume density in these coordinates
is rho(t)=w^-3/2 dw/dt. Put S=r+w0/r, so S^2=t^2+4w0 and S is an
even analytic function of t on the branch through t=0. A direct identity is

    rho(t)=2/(w+w0)=(1/w0)(1-t/S).                      (11)

Its even part is exactly 1/w0. Choose a quadratic vanishing sphere
in coordinates making 2XY+T^2-t^2 a sum of four squares. Reflection
t -> -t kills the odd part of rho in its period (equivalently in the
associated thimble integral). Only the constant part remains. Thus
the K2 period is exactly linear, at both points, to all orders.

For K3_rho the two critical points with w!=0 have the local
Jacobian-one change y -> y+(rho-z)/w, reducing exactly to K2.
For its third point use the actual global JQ suspension

    K3_rho=2PV+B_rho(u,v),
    B_rho=u(uv+rho)^2-u+uv+rho.

The JQ source change has constant nonzero Jacobian. Its third plane
critical point is u=0,v=1-rho^2 with value rho and plane Hessian
determinant -1, even when rho=+-1. Put L=uv+rho. On the nearby level
B_rho=s, solve the branch of u L^2+L-u=s through L=rho:

    L=(-1+sqrt(1+4u(u+s)))/(2u).

The plane Gelfand--Leray form is, up to sign,

    du/[u sqrt(1+4u(u+s))].                            (12)

Its residue at u=0 is identically one for every s near rho. A small
u-circle represents the local vanishing annulus, just as in (10).
Its plane period is therefore constant, and hyperbolic suspension
makes the three-cycle period exactly linear. No division by rho^2-1
was used, so the colliding critical values rho=+-1 are included.

The quadratic endpoint passes by scaling and the coordinate endpoint
has no critical point to test. This proves JU2. In particular the
new local test does not decide K2 versus h1, any S_rho source orbit,
or either old endpoint. It does not exclude different global periods,
behaviour at infinity, or quotient/torsor obstructions.

## 6. Consequence for the research frontier and assurance

Retain the JT exclusion of every standard embedded quotient and its
value-preserving conjugates. For a genuinely different prospective
action, require the reduced plinth, exact missing quotient lines,
contracted fixed-plane images, and the socle torsor class with (5),(6).
Then recover the complete polynomial coordinate inverse and the SAME
base pair in JT2, or prove that this is impossible. Neither local
quadratic normal forms nor more terms of the local volume periods
address that missing global recovery.

The checker provides exact finite identity, residue, degeneracy and
scope controls for the written all-degree arguments. It is not
independent replay, formal verification or external review. No
nonstandard necessary action, global coordinate certificate or new
source member is constructed. The normalized full-source questions,
both old endpoints, general graph embedding recovery and unrestricted
HC4/JC2 remain open. The original goal is not complete.
