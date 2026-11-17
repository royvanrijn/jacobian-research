# Canonical deck jet fields and rational pole-atlas coverage

PF2RULED157. 17 November 2026. Actual complex polynomial PLANE Keller maps.

This is a coverage refinement for the ENTIRE nontrivial generic-deck
case, with every cyclic order retained. It does not assume MX0 or a
quadratic hypersurface model. The no-proper-equivariant-factor residual
of a minimal-degree counterexample is exactly the TRIVIAL-DECK case.
Quadratic images receive a complete explicit presentation in Section7.
The nontrivial-deck source images are described by one trace-free rational
primitive, its actual translated-plane atlas, and computable principal
part partitions. All source and target maps remain the original ones.

The same normal form gives a lower bound 2g on both reduced ruling
degrees and a polar relative quadratic first-velocity trace for every
cyclic order. Neither field reconstruction nor this polar trace proves
integrality or excludes the cyclic target. Full finite normalization,
codifferent and ALL original ramified primes remain in scope.

## 1. Retained geometry and notation

Let B=C[P,Q] subset R=C[x,t], Jac_(x,t)(P,Q)=1, F=Frac(B), L=Frac(R).
Assume the actual map is noninvertible. Let S be the FULL finite
normalization of B in L and Lambda its codifferent. The ENTIRE etale
completion U is the open retained in MP0/LA0. It is not all of Spec(S).

The retained [MP0 common-pencil theorem](../../hessian/hc4-triangular-graphs/COMMON_UNRAMIFIED_PENCIL_AND_AUTOMATIC_ADDITIVE_COMPLETION.md)
makes the generic rational deck group D cyclic. If |D|=g>1, there are
polynomial source coordinates (x,t), centered at the faithful base fixed
point, such that a generator acts by

    tau(t)=zeta*t, tau(x)=zeta^(-1)*x+b(t), b in C(t),             (1)

where zeta has exact order g. The multiplier follows from preservation
of the ACTUAL two-form dx wedge dt=dP wedge dQ.

By [MQ0](../../hessian/hc4-triangular-graphs/CYCLIC_DECK_SOURCE_IMAGE_AND_THE_UNIQUE_MULTIPLE_FIBRE.md),

    E=L^D, W=union_j tau^j(X), V=W/D,
    X=A2 -> V SURJECTIVE ETALE, W->V finite etale degree g,
    W,V smooth AFFINE, O(V)=R intersect E.                       (2)

The t=0 fibre of W has exactly g disjoint reduced A1 components,
including the original source line. Its quotient is one A1 of
multiplicity g in u=t^g. All other fibres are reduced. W may be smaller
than the entire U, and V may be smaller than the entire etale locus of
the finite intermediate normalization. Neither replaces S or Lambda.

## 2. Exact reconstruction by finite source jets

For delta=p(t)*partial_x, let

    J_delta=C[delta^j P,delta^j Q : 0<=j<=M],
    M=max(deg_x P,deg_x Q), E_delta=Frac(J_delta).

This is a finite list; higher derivatives vanish. The retained
[MR0 replica-index theorem](../../hessian/hc4-triangular-graphs/REPLICA_DEGREE_AND_SOURCE_DERIVATIVE_FIELD_RECONSTRUCTION.md)
gives [L:E_delta] dividing deg(p)+1. For monomial p=t^k it also makes
every proper E_delta an actual cyclic deck fixed field of an order
m dividing k+1. The E_delta=F case is excluded by the retained
source-LND target obstruction. No normalization of J_delta is assumed.

For every k>=0, conjugation of the field derivation gives

    tau*(t^k partial_x)*tau^(-1)=zeta^(k+1)*t^k partial_x.

Let h=gcd(g,k+1), and let H<=D be its unique order-h subgroup.
H commutes with delta and fixes P,Q, hence fixes every finite jet:
E_delta subset L^H. Thus [L:E_delta]>=h. If E_delta is proper, MR0
makes it a deck fixed field of order m dividing both g and k+1, so
m<=h. Equality follows, including h=1. Consequently

    E_(t^k partial_x)=L^H,
    [L:E_(t^k partial_x)]=gcd(g,k+1)  for ALL k>=0.                (3)

In particular

    E_(t^(g-1) partial_x)=L^D=E.                                 (4)

This does not require degree minimality. The actual generic quotient
has [E:F]>1 by MQ0's retained Galois-Keller obstruction, so (4) is a
PROPER intermediate F subsetneq E subsetneq L.

There is a stronger formula at MINIMAL GENERIC DEGREE. For any nonzero
p(t)=sum p_n t^n put

    h=gcd(g,{n+1 : p_n!=0}).

The order-h subgroup H is exactly the subgroup commuting with delta.
Again E_delta subset L^H. If E_delta is proper, PF142 makes it L^K for
an actual subgroup K<=D. Since delta preserves E_delta, every element
of K commutes with delta: a derivation on E_delta has a UNIQUE extension
to the finite separable L. Thus K subset H. The reverse inclusion
follows because H fixes E_delta. Hence, at minimal generic degree,

    E_delta=L^H, [L:E_delta]=h.                                  (5)

Formula (5) is stated in the CANONICAL common-pencil frame. It does not
claim every source LND in every frame preserves that pencil.

PF142 already makes every proper source-LND-equivariant factor cyclic
at minimal generic degree, in every source frame. Combining that theorem
with (4) gives the exact global dichotomy

    no proper source-LND-equivariant factor  <=>  D is trivial.   (6)

When D is trivial, every finite source-LND derivative field is L, in
all source frames. This is a FIELD condition; that trivial-deck
residual is not excluded. Equations (3)--(6) do not make the raw jet
algebra normal, finite over B, or equal to the full S.

## 3. A single rational primitive gives the actual deck action

Let Pi=(1/g)Tr_(L/E) and start with s0=x/(g*t^(g-1)). From (1),
Pi(s0)-s0 is a rational function of t. Set

    s=Pi(s0)=x/(g*t^(g-1))+beta(t), Pi(beta)=0.

Then

    u=t^g, E=C(u,s), ds wedge du=dx wedge dt.                    (7)

Indeed L=C(t,s) and its degree over C(t^g,s) is exactly g. Also E=F(u)
by MQ0's coordinate-primitive argument; it need not be a polynomial plane.

The actual target primitive P*dQ and the full source polynomial
Poincare argument of
[MV0](../../hessian/hc4-triangular-graphs/TRACE_FREE_RATIONAL_SLICES_AND_MIXED_SOURCE_IMAGE_GRAPH_MODELS.md)
give a UNIQUE trace-free G in C(t), and a polynomial F0 in R, with

    dG=beta*du, Pi(G)=0, F0-G in E.                              (8)

For clarity, pull the regular target primitive to X and write it as
x*dt+dF0. Taking normalized field trace and subtracting gives

d(F0-Pi(F0))=beta*du. Its primitive is rational in t because the
differential has zero dx coefficient. Set G=F0-Pi(F0).
The condition F0-G in E is retained: rational exactness alone does not
supply a global polynomial representative for arbitrary proposed data.

Let H(t) be the polynomial part of G. Averaging respects polynomial
parts, so Pi(H)=0. The legitimate polynomial source shear

    x_new=x+H'(t), G_new=G-H, F0_new=F0-H

preserves (7)--(8), B,F,S,Lambda and the common source pencil.
We may therefore require G to be STRICTLY PROPER rational. Its finite
principal parts determine it completely. It is not set equal to zero.

Writing X0=x+G'(t), we obtain the exact normal form

    s=X0/(g*t^(g-1)), tau(X0)=zeta^(-1)*X0,
    tau(x)=zeta^(-1)x+zeta^(-1)G'(t)-G'(zeta*t).                  (9)

The rational pole data, including all character components, are retained.
For any subgroup H of order h>=2, the same action is described by
G_H=G-Pi_H(G), where Pi_H averages over H. The same polynomial F0 works:
F0-G_H=(F0-G)+Pi_H(G) belongs to L^H. Thus a single full-deck primitive
supplies all the subgroup pole normal forms selected by (3)--(5).

## 4. Exact pole-atlas and image algebra coverage

Index the translated charts as U_j=tau^(-j)(X), and use the global t
and normalized fibre coordinate

    w_j=zeta^j*tau^j(x)=x+xi_j(t),
    xi_j=(G(t)-G(zeta^j*t))', 0<=j<g.                            (10)

These are WHOLE planes. Their ACTUAL intersections are

    U_i intersection U_j = D(d_ij(t)) times A1,
    w_j=w_i+xi_j-xi_i,                                         (11)

where d_ij has exactly the finite poles of xi_j-xi_i as its roots.
If the difference is regular at a, its transition extends over the
whole fibre; separatedness identifies the two chart maps there. A pole
prevents overlap at every point of that fibre. This is the complete
atlas argument of PF138, now with all transitions determined by G.

In the common function field, the ACTUAL image algebras are

    T_G=O(W)=intersection_(j=0)^(g-1) C[t,w_j],
    C_G=O(V)=T_G^D=R intersect C(t^g,s).                         (12)

Affineness and finite generation in (12) come from the actual W,V in
(2). No claim is made that arbitrary rational G automatically gives an
affine surface or a plane target. Formula (12) is an exact normal form
for EVERY actual nontrivial-deck source image, not a classification
of admissible arbitrary pole data up to isomorphism.

Let PP_a denote the finite Laurent principal part at t=a. Because a
rational derivative has a pole exactly where its primitive has a pole,
the whole fibre components in (11) are indexed by the equivalence classes

    i~_a j iff PP_a G(zeta^i*t)=PP_a G(zeta^j*t).                 (13)

Thus all fibre component counts and source sheet counts are computable
from the finite principal parts. At a!=0 the number of classes is the
number of reduced V-components over u=a^g, and their cardinalities are
the source fibre cardinalities on those components. This includes
intermediate clusters, not just fully split MX0-type fibres. MQ0's
Euler identity becomes chi(V)=1+sum_(nonzero pole orbits)(#classes-1).

At zero write the finite principal part of G as

    PP_0 G=sum_(r in A0) a_r t^(-r), a_r!=0.

The g central translated components are distinct by freeness. Hence

    A0 is nonempty, gcd(g,{r:r in A0})=1.                       (14)

Indeed the subgroup fixing this principal part has order equal to that
gcd; a nontrivial stabilizer of a central A1 component would have a fixed
point. For j=1,...,g-1 define

    nu_j=max{r in A0 : g does NOT divide j*r}.

All these numbers exist by (14). In chart U_j the original x=w_j-xi_j
has EXACT pole order nu_j+1 along its t=0 component. Normalized trace
zero removes all r divisible by g, and max_j nu_j=max A0=nu_*>=1.

This is the explicit rational-pole residual replacing an unspecified
cyclic image model. It includes all cyclic orders and mixed characters.
It does not assert that MX0 or S2 covers it, nor that an arbitrary member
has the required affine geometry, polynomial representative or target.

## 5. A coverage degree bound on every target coordinate

For any polynomial target coordinate H (including after a target
automorphism) put

    r_H=deg_s H=deg_x H, b_H=deg_x H(x,0).

The first equality holds because H in E is a polynomial in s over C(u):
substitution (7) is affine linear in x, and invariance fixes its coefficients.
PF52's actual reached-line theorem gives b_H>=2 for a noninvertible map.

More generally, for ANY invariant source polynomial H with b_H>0,

    r_H >= g*b_H.                                              (15)

Use the finite field extension L/C(t,H), of degree r_H. Its FULL DVR
normalization above the prime t=0 has the g central W-components among
its places. Each has ramification index one and residue degree b_H.
The fundamental degree inequality, retaining all other places, gives
r_H>=sum_(central places) e*f=g*b_H. This comparison normalization is
used only for degree bookkeeping; it does not replace the ORIGINAL S.

Therefore every polynomial target coordinate has r_H>=2g. The actual
plane identity in coordinates (u,s) is

    P_s Q_u-P_u Q_s=1.

Its highest s-coefficient gives the scalar one-variable Wronskian
sdeg*l_P'*l_Q-rdeg*l_P*l_Q'=0. The same legitimate elementary target
reductions as PF151 terminate with

    2g<=r<s, neither r nor s divides the other,                  (16)

and common leading powers l_P=alpha*H0^(r/d), l_Q=beta*H0^(s/d),
d=gcd(r,s). These are RULING degrees, not total or field degrees.
The bad common leading-divisor problem is not thereby eliminated.

## 6. Generic first velocities have a polar relative quadratic trace

For any invariant target function f use coordinates (t,X0) in (9) and set

    X_f=f_(X0)=f_x, A_f=f_t|_(X0).

Then tau(X_f)=zeta*X_f, tau(A_f)=zeta^(-1)*A_f, and

    f_t|x=A_f+G''*X_f,
    a_lambda=(lambda*partial_x+partial_t)f
                  =lambda*X_f+A_f+G''*X_f.                     (17)

The actual pair has X_P*A_Q-A_P*X_Q=1. Put
M_ij(lambda)=Tr_(L/E)(a_lambda,i*a_lambda,j). Since Pi(G)=0 and
Theta=t*partial_t commutes with averaging,

    Pi(t^2 G'')=Pi(Theta(Theta-1)G)=0,
    Tr(G''*X_i*X_j)=0.                                         (18)

For the second identity use that X_i X_j/t^2 belongs to E. Therefore

    [lambda]M_ij=g(X_i A_j+A_i X_j),
    det([lambda]M)=-g^2.                                      (19)

This is the determinant of the LINEAR COEFFICIENT MATRIX, not an
assertion that det M(lambda) is constant. All its constant terms remain.
The quadratic lambda coefficient is zero for g>=3. For g=2 it is
2X_i X_j, regular on V because it is invariant and polynomial on X.

For generic nonzero constant target linear f, X_f(x,0) is not identically
zero. Since G'' has pole order nu_*+2, (17) makes A_f have precisely
that pole on the original source line; f_t itself is polynomial.
Thus [lambda]Tr(a_lambda^2)=2g X_f A_f has exact pole nu_*+2 on
the RETAINED central curve of V. For all but at most one lambda the
relative quadratic trace is nonregular there. This is an all-order
pole identity for arbitrary admissible G, including mixed characters,
not a finite unmarked-moment prefix and not just the MX0 subclass.

By PF143's full codifferent argument there is an integral marker in
the actual intermediate normalization S_E with

    Tr_(L/F)(s_marker*a_lambda^2) not in B, s_marker in S_E subset S. (20)

A polar relative trace on the retained curve is outside the codifferent
of S_E/B; transitivity proves (20) without guessing cancellation of sheets.
No marker in O(W) is silently declared integral over B.

The explicit pole packet also strengthens its quantitative ledger.
On the central omitted component in U_j, q_j=nu_j+1 and the ORIGINAL
ramification index is e_j=1, because it lies in the actual etale locus.
All these components map to the same original target curve Gamma.
For a generic first velocity the exact FULL-S marked quadratic exponent
m_Gamma of PF143 consequently satisfies

    m_Gamma >= 2*nu_*+4 >=6.                                   (21)

Other omitted primes above Gamma, including ramified ones outside W,
remain in the full maximum and can increase this bound. Equality is
not asserted from the central packet alone.

## 7. Complete quadratic-image presentation, beyond MX0

Every actual quadratic top factor L/E is Galois, so it is the fixed field
of the unique order-two deck subgroup. By (3), it is already the finite
jet field for t*partial_x in the canonical frame. No minimality is needed.

For g=2, trace-free G is odd. In the proper normalization it has exactly

    G(t)=t*R0(u), u=t^2,

where R0 in C(u) is STRICTLY PROPER and has a pole at u=0. Write

    H0(u)=R0(u)+2u R0'(u)=Sigma(u)/A(u),                         (22)

in reduced form with A monic. Then

    gcd(A,Sigma)=1, A(0)=0, Sigma(0)!=0,
    deg Sigma<deg A.                                           (23)

If R0 has pole order k0>=1 at zero and order m_a>=1 at each nonzero
finite pole a, the denominator is exactly

    A=u^k0 product_(a!=0)(u-a)^(m_a+1).                         (24)

At zero the factor 1-2k0 is nonzero; at a!=0 differentiation raises the
pole order by one. At infinity R0~c*u^(-l), l>=1, and the factor 1-2l
is nonzero. These observations prove (23)--(24), with no constant
pole-order or square-free-denominator restriction.

The ACTUAL quadratic source image has the explicit affine presentation

    V_R0={w^2=u*Sigma(u)^2+u*A(u)*z},                            (25)
    u=t^2, w=t*(Sigma(u)+A(u)*x),
    z=2*Sigma(u)*x+A(u)*x^2.                                   (26)

To prove coverage, put v=Sigma+A*x. The two charts of its finite cyclic
cover are v=+Sigma+A*x and v=-Sigma+A*y on

    W={A(t^2)*z=v^2-Sigma(t^2)^2}.

They are whole planes, with actual overlap D(A), because A and Sigma
are coprime. Their transition is y=x+2Sigma/A=x+2G'(t), exactly
(10). Thus they recover the ACTUAL W from its full atlas, not just its
function field. The involution (t,v,z)->(-t,-v,z) is free, including
at t=0 where v!=0. Its invariant ring is generated by u=t^2,w=t*v,z,
with the sole relation (25). This proves the exact marked-image
presentation for every actual quadratic subgroup, including order-two
subgroups of larger D. In the latter case first use G_H=G-Pi_H G.

Conversely EVERY proper rational R0 with a pole at zero gives this
smooth affine W and free involution. Off A=0 smoothness uses the z
derivative; at A=0 it uses v!=0. The quotient V_R0 is smooth affine.
The two whole charts cover W, so (26) is SURJECTIVE ETALE of generic
degree two, and O(V_R0)=C[x,t] intersect C(t^2,s). The source volume is

    (du wedge dw)/(uA) = -2 dx wedge dt.

The normalized ds wedge du is therefore its regular unit volume.
No etale plane TARGET is constructed by this converse.

There is an explicit polynomial primitive representative as well.
Since A*R0 is polynomial, choose C0 in C[u] with

    C0*Sigma=-A*R0 modulo A,
    D0=R0+C0*Sigma/A in C[u].

Such C0 exists by (23). Then

    F0=t*C0(u)*x+t*D0(u), F0-G=2u*C0(u)*s in E.                 (27)

This verifies ALL of (8), not only rational exactness. It also proves
that every model in this explicit quadratic class has exact unit volume.
Distinct R0 are not asserted to classify isomorphism classes injectively.

The retained S2 model is R0=1/u. MX0's quadratic family is the special
choice R0=1/(u*B(u)^(k-1)), with monic square-free B, B(0)!=0:

    A=u*B^k, Sigma=-B+2(1-k)u B'.

Arbitrary leading constants of B are handled by constant normalization.

Thus (25) is a COMPLETE quadratic-image coverage class, rather than a
claim that MX0 exhausts it. Central pole order and all nonzero pole
orders may vary independently through (24).

For every such model M=uA has degree at least two. Hence PF149's
retained even-leading criterion excludes ANY actual target coordinate
w*A_P+B_P with deg_z B_P=n>=1, deg_z A_P<n and constant nonzero leading
coefficient of B_P. Its infinity inequality is n*deg M>=2n and there
are no leading-coefficient roots to check. This applies to ALL quadratic
images covered here. No theorem makes that leading coefficient constant.

The central multiplicity is ord_0 M=k0+1, not always two. In a coupled
reduced target frame let d=gcd(r,s), and l_P=alpha*H^(r/d),
l_Q=beta*H^(s/d). The correct general central valuations are

    ord_0 c_P=(r/d)ord_0 H+(k0+1)floor(r/2),
    ord_0 c_Q=(s/d)ord_0 H+(k0+1)floor(s/2).                      (28)

Here c_P,c_Q are the canonical z-leading coefficients of B or A according
to parity, as in PF151. Its double-zero odd conclusion was for k0=1;
it is not silently imposed on all of (25). The central packet of this
quadratic image alone gives m_Gamma>=4k0+2 by (21). The common bad
leading-divisor normalization problem remains G3.

There is a sharper consequence at the lowest ruling degree. If r=4,
PF52 makes c_P nonzero at EVERY M-root: otherwise P has degree<=1
on the corresponding whole source line, where w is constant and z is
affine linear. If s were odd, gcd(4,s)=1. Writing m0=ord_0 M>0,
(28) would give 4ord_0 H+2m0=0 and then
ord_0 c_Q=-m0/2<0, impossible. Therefore EVERY actual quadratic image
with reduced r=4 has

    s=4j+2, d=2, c_P=alpha*G0^2, c_Q=beta*G0^(s/2),
    gcd(G0,M)=1, some root of G0 has multiplicity>=3.              (29)

Here G0=H*M is polynomial by its two polynomial positive powers.
If all its roots had multiplicity<=2, PF149 would exclude c_P; its
infinity condition again follows from deg M>=2. This extends PF151's
lowest-degree coverage to ALL quadratic images in (25), including k0>1.
The bound (16) also forces the FULL deck order to be two when r=4.
It does not eliminate the bad common divisor in (29).

## 8. Controls, dimension two and the exact remaining gaps

For the proper primitive G=t^(1-g), a whole source invariant polynomial is

    Z=((t^g*x-(g-1))^g-(-(g-1))^g)/t^g.

It has deg_x Z=g and central degree one, so Z^d attains equality in
(15) for every d>=1. Also F0=t*x/(g-1) is a polynomial representative:
F0-G=(g/(g-1))*u*s. These are controls for the source pole geometry and
degree bound; no polynomial Keller target pair is constructed.

The finite arithmetic checks use ACTUAL multiplication in C(u)[t]/(t^g-u)
for field traces, and include mixed central pole characters and a nonzero
pole orbit with cluster cardinalities (2,4). They do not substitute this
field algebra for S or Lambda. Outside the explicit quadratic construction,
arbitrary pole controls are not asserted to be admissible affine source
images. None is asserted to have a plane Keller target.

The dimension-two steps are identifiable: the retained common SOURCE
A1-pencil; the scalar volume multiplier in (1),(7); and the FINITE
extension L/C(t,H) in (15). The latter would have positive transcendence
degree for a three-variable source and two comparison functions. The
known nonproper three-dimensional Keller control therefore supplies no
counterexample to these hypotheses, and still rules out dimension-free
integrality claims from constant Jacobian or derivative generation.
No automatic regularity or integrality is inferred here.

The G1 residual is now exactly trivial deck, with full finite-jet FIELD
reconstruction for every source LND at minimal degree. The G2 actual
nontrivial-deck images are covered by (7)--(14); its quadratic images have
the complete presentation (25)--(27). Higher-order pole-data admissibility
and target realization remain unresolved. G3's common bad leading
divisors remain, now alongside the all-order bound (16). The full marked
quadratic regularity needed for generic integrality is still UNSUPPLIED.
Original non-equivariant factors, the entire finite normalization and
every ramified boundary prime remain relevant. PF142's mixed minimal
case stays retired; PF156's fixed-v2 quintic graph branch stays closed.
Full HC4/JC2 and an actual counterexample remain OPEN. No new external
classification or predecessor checker replay is used. Independent whole-
proof replay, formal verification and external review are unclaimed.
The narrow checker passes 118 exact arithmetic and symbolic controls.
