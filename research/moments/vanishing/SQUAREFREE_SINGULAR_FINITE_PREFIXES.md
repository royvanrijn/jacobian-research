# Explicit pure prefixes for concurrent, cuspidal and tangent leading cubics

Current quantitative refinement, 15 November 2026: the [exact Euler representative proof](EULER_REPRESENTATIVES_GLOBAL_297.md) gives global pure297, odd296, conic/full-W-linear258 and all-homogeneous258. Ordinary-node complements are297/271/244, with homogeneous nonisotropic-node101; cusp255/tangent284 are unchanged. Current intervals are global[20,297], odd[20,296] and all-homogeneous[14,258]. Earlier353/336 bounds below are historical valid coarser evidence. Leastness,20/14 sufficiency, expanded I297 cofactors and independent review remain open.

Current quantitative refinement, 15 November 2026: the [joint-degree actual-seed proof](JOINT_DEGREE_GLOBAL_353.md) gives global pure353, odd352, full W-linear336, rank-one140, mixed-parent255 and cusp255/tangent284. Current least-prefix intervals are global[20,353], odd[20,352] and all-homogeneous[14,336]. Earlier590/888 bounds and rank-one425/mixed-parent320 estimates below are historical valid coarser evidence. Leastness,20/14 sufficiency, expanded I353 cofactors and independent review remain open.

Current quantitative refinement, 15 November 2026: the [balanced logarithmic seed proof](BALANCED_LOGARITHMIC_GLOBAL_590.md) gives global pure590, all-homogeneous336, conic336 and full W-linear425. Current intervals are global/odd[20,590] and all-homogeneous[14,336]; smooth[20,94] and homogeneous-smooth[14,56] are unchanged. Earlier888/1855 bounds and frame-plus888 statements below retain their historical scopes. Leastness,20/14 sufficiency, expanded I590 cofactors and independent review remain open.

Current global refinement,15 November2026: the [rescaled singular seed proof](RESCALED_SINGULAR_GLOBAL_888.md)
now proves888 sufficient as a global pure-only prefix, without an extra
frame condition. Current intervals are global[20,888], odd[20,888] and
all-homogeneous[14,888]. Earlier1855 bounds and reservations about888
alone below are historical. The new proof remains internal and
independently unreviewed; leastness and expanded global cofactors are open.

Current quantitative status,15 November2026: the [global pure-prefix proof](GLOBAL_TERNARY_CUBIC_PURE_PREFIX.md)
now gives a sufficient1855 pure-only gate and least-prefix interval[10,1855].
Earlier absence-of-global-upper or unquantified-node statements below are
historical. Their earlier class/certificate scopes remain unchanged; the
new proof remains internal and independently unreviewed.

15 November 2026. Internal written quantitative deductions from retained
analytic exclusions and exact filtered reductions. The checks do not
certify the analytic input. No independent specialist replay, external
review or formal verification. No moments are expanded to the stated gates.

## Results

For complex ternary polynomials of degree3 with arbitrary lower terms:

| Claim | Leading projective curve | Sufficient initial pure-zero gate |
|---|---|---:|
| **G3CONCURRENTFINITE20261115** | Three distinct concurrent lines | 104 |
| **G3CUSPFINITE20261115** | Irreducible cuspidal cubic | 720 |
| **G3TANGENTFINITE20261115** | A line tangent to a smooth conic | 935 |

At the stated gate, all positive pure moments vanish and mixed moments
vanish for m>3*deg Q. These sufficient pure gates are not claimed optimal.
On each class, the complementary non-frame locus has a nonzero pure
moment by the listed bound; the frame locus uses the smaller retained
92 or425 gate as described below. No global arbitrary-cubic pure-only
upper bound is asserted.

Together with the retained smooth-leading188 and repeated-line99 bounds,
these results leave the three ordinary-node leading types unquantified:
irreducible nodal cubics, transverse line-smooth-conic curves, and three
independent lines. The numerical least global prefix remains
**OP-GMC3-PURE-PREFIX**, open with lower bound10 and a terminating exact
search in principle. The ternary cubic GMC decision remains internally
closed; analytic review, scheme/term questions and lower-power GVC remain separate.

## 1. A common exact quadratic elimination, including an isotropic direction

Put the unique leading singular point at the T-axis. In the concurrent
case, a binary linear change and constant rescaling put its squarefree
binary cubic into u^3+v^3: over C a projective linear transformation sends
any three distinct points of P^1 to any other three, and an overall
scalar can be absorbed in the coordinate change. For the other cases use
the retained canonical leading forms. Thus write

    P=alpha*T^2+T*(H2(u,v)+ell(u,v))+G3(u,v)+g2(u,v),
    Q=m*T^2/2+T*J(u,v)+q2(u,v), det Hess Q!=0,       (1)

where ell,J have degree<=1, g2,q2 degree<=2, and Q is the full affine
Gaussian exponent. All means and mixed terms remain arbitrary.

| Type | H2 | G3 |
|---|---|---|
| Concurrent | 0 | u^3+v^3 |
| Cusp | u^2 | v^3 |
| Tangent | u^2 | u*v^2 |

Set D=m-2*alpha*t and L=t*(H2+ell)-J. If D!=0 in C[t], completing T
is exact over C(t), with normalized seed

    I_Gamma=constant*D^(-1/2)*integral exp(U/D) du dv,
    U=D*(t*(G3+g2)-q2)+L^2/2.                    (2)

In every type deg_t D<=1 and deg_t U<=2. This identity is a generic-t
twisted-module calculation, retaining the original formal Wick functional
by Gaussian integration by parts. When D(0)=0 it is not a conditional
t-adic integral in T, and no convergence is assumed.

If D=0 identically, alpha=m=0. The vertex is isotropic and F2 vanishes
there, so it supplies a W-linear Gaussian frame. Section5 treats this
branch by the retained finite W-linear results. Every D!=0 member is
pure-excluded by the retained class proofs: alpha!=0 is gradient-coercive;
alpha=0,m!=0 is the complete nonisotropic formal-frame exclusion. In
particular the normalized seed has no nonzero constant actual period.
We use that stronger period conclusion, not just a claim about pure nullity.

## 2. The fixed leading algebras and weighted bases

Write U=U_top+R. With the positive integer weights below, R is strictly
lower than U_top. Let a=t^2/2, b=tD when those names occur.

| Type | Weights(u,v) | U_top | Bound on weight R | Basis | rank |
|---|---|---|---:|---|---:|
| Concurrent | (1,1) | tD*(u^3+v^3) | 2 | 1,u,v,uv | 4 |
| Cusp | (3,4) | a*u^4+b*v^3 | 10 | 1,u,u^2,v,uv,u^2v | 6 |
| Tangent | (2,3) | a*u^4+b*u*v^2 | 7 | 1,u,u^2,v,v^2 | 5 |

All listed top coefficients are nonzero over C(t). Their gradients form
a weighted homogeneous regular sequence. The first two rows have pure
power gradients. For the last row,

    g_u=4a*u^3+b*v^2, g_v=2b*u*v,

are coprime and have weighted degrees6,5. Their quotient has dimension
6*5/(2*3)=5. The following exact principal divisions prove the listed
basis explicitly; monomials already in the basis are unchanged:

    u^i*v^j, i,j>=1: (u^(i-1)*v^(j-1)/(2b))*g_v;
    u^3: -(b/(4a))*v^2+(1/(4a))*g_u;
    u^i, i>=4: (u^(i-3)/(4a))*g_u-(u^(i-4)*v/(8a))*g_v;
    v^j, j>=3: (v^(j-2)/b)*g_u-(2a*u^2*v^(j-3)/b^2)*g_v. (3)

There are no omitted parameter opens: a,b are already nonzero by D!=0.
The associated weighted Koszul complex of the ordinary twisted derivatives
partial_i+U_i/D has this regular symbol sequence. Its leading-term
exactness makes these representatives a full quotient basis. This is
the same filtered regular-sequence argument used in the retained smooth
and fixed D5 proofs. Derivative corrections must be retained; ordinary
critical-ideal reduction alone is not the quotient computation.

## 3. Uniform connection degrees from the actual divisions

For a principal decomposition f=NF+sum h_i*(U_top)_i, integration by
parts gives, in the twisted quotient,

    [f]=[NF-sum h_i*R_i-D*sum partial_i h_i].       (4)

Each correction strictly decreases integer weight. Choose a common
principal divisor A as follows:

| Type | A, up to nonzero constants | deg_t A<=q | Minimum weight of a reducible monomial |
|---|---|---:|---:|
| Concurrent | tD | 2 | 2 |
| Cusp | a*b | 4 | 8 |
| Tangent | a*b^2 | 6 | 5 |

In every row A times each principal cofactor in(3), or in pure-power
division, is polynomial in t. Its degree is at most q-2. The normalized
principal remainder has degree<=q. Since deg_t R_i<=2 and deg_t D<=1,
each entire right side of(4), after multiplication by A, has degree<=q
per scalar coefficient. This includes the same-weight basis remainder
of u^3 in the tangent row. All other remainders have smaller weight.

The normalized connection source is

    C=(U_t*D-U*D')/D^2-D'/(2D).                  (5)

Its numerator has t-degree<=2 and weight at most that of U_top. A
connection column starts with source numerator times a basis element.
The maximum weights are5,22,14 respectively. Weight drops by at least
1,2,1 per non-basis reduction, because top/lower gaps are1,2,1;
the derivative corrections drop by the entire top weight. Hence the
maximum number r of successive corrections along a branch is

    concurrent: (5-2)+1=4;
    cusp: floor((22-8)/2)+1=8;
    tangent: (14-5)+1=10.                       (6)

Branching does not increase depth. Basis terms are padded by the
remaining power of A. Thus the full connection matrix has a common
denominator H=D^2*A^r and polynomial numerators N with

    deg_t H<=B=2+q*r, deg_t N_ij<=B.             (7)

The resulting B values are10,34,62. Coefficient specializations can
lower these degrees, but cannot increase them. Normal-form coefficients
and constant coordinate changes affect only constants in t. We do not
assume an ordinary Morse critical set, distinct critical values or a
nonvanishing Hessian of P.

## 4. Nonzero zeroth coefficient and finite pure exclusion

Apply the retained [seed cofactor/initial-order argument](REPEATED_LINE_FINITE_PREFIX.md).
At minimal seed rank s<=mu, the iterated connection columns have
numerators of degree<=jB over H^j. Polynomial cofactors produce

    sum_(j=0)^s q_j(t)*I^(j)=0,
    deg q_j<=B*s*(s+1)/2.

On D!=0, the complete retained no-constant-period theorem forces q0!=0:
otherwise the minimal seed equation would admit the constant solution,
and the seed-submodule/actual-period bridge would produce a forbidden
nonzero constant normalized period. Seed cyclicity in the full ambient
module is not assumed. This is the analytic dependency of the bounds.

If I=1+O(t^(K+1)), the derivative terms have order greater than deg q0
when K=B*s*(s+1)/2+s. They cannot cancel the nonzero polynomial q0.
The worst bounds at ranks4,6,5 are therefore

    10*4*5/2+4=104;
    34*6*7/2+6=720;
    62*5*6/2+5=935.                             (8)

This works at D(0)=0 as well: the final relation is polynomial and acts
on the original formal Wick series, which has constant term1. We do not
evaluate a singular eliminated integral at t=0.

## 5. The isotropic zero-F2 frame branches

When alpha=m=0, use a circular Gaussian frame with its W direction at
the leading singular point. For concurrent distinct lines, F3 is
independent of W, so its W coefficient is affine linear. This lies in
the d=0 face of the retained fourteen-coefficient theorem, with pure
gate92, below104.

At a cusp or an ordinary tangency, the projective tangent cone at the
singular point is a nonzero rank-one quadratic q2(Z,T). In the W-linear
expression P=W*(a+bZ+cT+q2)+G, write q2=dZ^2+eZT+fT^2. Rank one says
e^2-4df=0. If f!=0 the retained Gaussian null rotation puts it into
f*T^2, with sufficient rank-one pure gate425. If f=0 then e=0, and
the fourteen-coefficient family has sufficient gate320. Both are below
720 and935. These gates imply all-order pure nullity, without asserting
that every frame is pure-null or that the whole frame locus is empty.

The retained W-linear mixed theorem gives m>3*deg Q on all pure-null
members. No sharpness of coefficient3 within each of these three
squarefree-leading subclasses is claimed.

## Checks, source scope and remaining work

The [checker](verify_squarefree_singular_finite_prefixes.py),
[receipt](squarefree_singular_finite_prefixes_verification.json) and
[source audit](squarefree_singular_finite_source_audit.json) distinguish
universal completion/weight/division identities, exact rational connection
models and arithmetic from written basis/completeness and analytic input.

The inspected [Euler-Mellin source](https://arxiv.org/html/2505.12458v1)
uses generic exponents and a regular holonomic proper-direct-image setting;
it is not a degree bound for the present exponential Gaussian seed.
[Solvable-algebra Groebner degree bounds](https://arxiv.org/abs/0710.4945)
likewise do not by themselves supply the missing specialized integration
presentation or a nonzero zeroth seed coefficient. Neither theorem is
imported here. The source audit records actual reading and version limits.

The global least pure prefix and every explicit arbitrary-cubic numerical
upper bound remain open. Ordinary-node leading types require a quantitative
reduction retaining all Newton faces; their pure exclusions are already
proved internally and are not reopened as qualitative GMC questions.

The next [exact ordinary-node checkpoint](ordinary_node_phase_shift_checkpoint.json)
and its [script](ordinary_node_phase_shift_checkpoint.py) retain a useful
normalization without claiming a prefix. After completing T in
P=alpha*T^2+T*(uv+h)+a*u^3+d*v^3+lower terms, translate
old_u=u+jv/t, old_v=v+ju/t. The mixed cubic phase coefficients cancel.
The new phase is V/(t^2*D), with deg_t V<=4,deg_t(t^2*D)<=3 and
normalized connection-source numerator degree<=6. Required nonzero
axis coefficients distinguish nodal, transverse and triangle types.
Full twisted divisions and denominator bounds remain to be constructed;
a commutative critical-ideal Groebner basis is insufficient.
