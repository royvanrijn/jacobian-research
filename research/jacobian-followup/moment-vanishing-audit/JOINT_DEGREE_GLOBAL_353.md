# Joint-degree actual seeds: global353, W-linear336 and rank-one140

15 November 2026. Internal written refinement of the
[balanced logarithmic590 proof](BALANCED_LOGARITHMIC_GLOBAL_590.md),
[rescaled singular seeds](RESCALED_SINGULAR_GLOBAL_888.md), and retained
[rank-one](RANK_ONE_WT2_CLASSIFICATION.md) and
[mixed-parent](MIXED_PERIOD_AT_INFINITY.md) analytic exclusions.

| Coefficient scope | New sufficient pure gate or exclusion |
|---|---:|
| Every complex ternary polynomial of total degree at most three | **353** |
| Odd ternary cubics | **352** |
| All sixteen W-linear cubic coefficients | **336** |
| Rank-one quadratic W coefficient | **140** |
| Fourteen-coefficient mixed parent | **255** |
| Irreducible cusp, arbitrary lower terms | **255**, complement225 |
| Tangent line-conic, arbitrary lower terms | **284** |

The global test uses only the initial pure moments, without an extra frame
hypothesis. The retained cubic theorem then gives the Gaussian W-linear
frame and sharp mixed cutoff m>3*deg Q. The least global prefix is in
**[20,353]**, odd-only **[20,352]**. All-homogeneous **[14,336]**,
smooth-leading **[20,94]** and homogeneous-smooth **[14,56]** are unchanged.
Both isolated lower witnesses remain complete. Neither leastness nor20/14
sufficiency is proved.

The [checker](verify_joint_degree_seed_bounds.py),
[certificate](joint_degree_seed_bounds_certificate.json),
[receipt](joint_degree_seed_bounds_verification.json) and
[source audit](joint_degree_seed_bounds_source_audit.json) bind exact auxiliary
algebra. Analytic exclusions, all-coefficient reduction inductions, complete
class coverage and nonzero zeroth coefficients remain written mathematics.
No independent replay, external review, formal verification, expanded global
radical cofactors or full coefficient-ring search is claimed.

## 1. Cusp and tangent: charge joint degree, not parameter degree alone

Retain the exact binary elimination and rescaling, with

    delta=t^2*D, deg_t D<=1, D!=0 in C[t],
    I=constant*t^-2*D^-1/2*integral exp(W/delta) dx dy,
    deg_(t,x,y) W<=4.

The source is still

    S=W_t*delta-W*delta_prime-2*t^3*D^2-t^4*D*D_prime/2,
    deg_(t,x,y) S<=6.

Both the volume factor and its -2/t term are retained. Actual derivatives
have R_j/delta^(2j), where

    R_0=1,
    R_(j+1)=delta^2*R_j_prime-2j*delta*delta_prime*R_j+S*R_j,
    deg_(t,x,y) R_j<=6j.

Use a=1/2,b=D. The weighted top phases and bases remain

    cusp: a*x^4+b*y^3, weights(3,4),
          basis 1,x,x^2,y,xy,x^2*y;
    tangent: a*x^4+b*x*y^2, weights(2,3),
             basis 1,x,x^2,y,y^2.

Their common principal divisors are A=a*b and A=a*b^2, of maximum
t-degree k=1 and2. The retained explicit principal divisions show that,
for a monomial of ordinary amplitude degree M, the cleared cofactors have
joint degree at most M+k-3. For example the tangent y^j division has

    A*h_x=a*b*y^(j-2),
    A*h_y=-2*a^2*x^2*y^(j-3),

both of joint degree at most j-1=M+k-3. The mixed-monomial and pure-x
families satisfy the same estimate. For x^3 the exceptional principal
remainder is -b^3*y^2/4; its joint degree5 equals M+k. Thus the complete
cleared remainder, including derivative corrections, has joint degree
at most M+k: W gradients have joint degree at most3 and delta has degree
at most3, with differentiation lowering degree. The old parameter-only
charges4/5 discarded the compensating decrease in amplitude degree.

Weight still strictly decreases by2 for cusp or1 for tangent, except for
already-normal basis terms. Consequently the old reduction depths remain

    r_j=6j-3 (cusp), r_j=8j-4 (tangent).

Padding branches by unused powers of A preserves their joint-degree
bounds. For N_j/H_j with H_j=delta^(2j)*A^r_j, the row of ordinary
basis degree d_i now satisfies

    deg_t N_(i,j)<=b_j-d_i,
    cusp: b_j=6j+r_j=12j-3;
    tangent: b_j=6j+2r_j=22j-8.

The denominator degrees are no greater than b_j and their valuations are
at least4j. Applying the retained actual-seed/cofactor argument gives

    cusp: 12*21-3*6-(0+1+2+1+2+3)=225;
    tangent: 22*15-8*5-(0+1+2+1+2)=284.

Every smaller seed rank and row subset containing the unit row is covered;
the full rank gives the largest bound. These are D!=0 complement
exclusions. The D=0 frame branches are included in section4, rather than
silently treating225/284 as full-class numbers.

## 2. Rank-one quadratic W coefficient: rescale the ordinary volume

For

    P=W*(a+bZ+cT+fT^2)+G(Z,T), f!=0, deg G<=3,
    q=1-bt, A(T)=a+cT+fT^2,

retain the exact phase and normalization

    Phi=t*G(t*A(T)/q,T)-T^2/2,
    I=constant*q^-1*integral exp(Phi) dT.

Set y=t*T and Ahat=a*t^2+c*t*y+f*y^2. The transformed normalization and
phase numerator are

    I=constant*(t*q)^-1*integral exp(W/delta) dy,
    delta=t^2*q^3,
    W=sum_(i+j<=3) g_ij*t^(3-i-j)*q^(3-i)*Ahat^i*y^j
           -q^3*y^2/2.

Here deg_(t,y) W<=6, so deg_t W_k<=6-k, and deg_t delta<=5.
The new -1/t is required because dT=dy/t is an ordinary volume form;
this differs from the unchanged logarithmic form dx/x in the conic proof.
The normalized source numerator is polynomial:

    S=W_t*delta-W*delta_prime-delta^2/t+b*delta^2/q,
    deg_(t,y) S<=10.

Relative to the old source pulled back to y, the source difference is
-Dop(y)/(t*delta), where Dop(h)=delta*h_y+W_y*h. Thus the same original
formal Gaussian seed is retained; no singular t=0 integral is introduced.

Let d be the actual y-degree of W, between2 and6. Degree cannot drop below2
because the original phase has Gaussian quadratic coefficient -1/2 at t=0.
Put A=d*W_d, nonzero, with deg_t A<=6-d. The full polynomial twisted basis
is1,y,...,y^(d-2). A leading cancellation of a monomial of degree M has
cleared remainder of joint degree at most M+(6-d), while lowering amplitude
degree strictly. This follows directly from deg_t W_k<=6-k and from the
strictly smaller derivative contribution delta*h_y.

Actual seeds R_j/delta^(2j) have joint degree at most10j and amplitude
degree at most d*j. Thus r_j=d*(j-1)+2 cancellations suffice and

    H_j=delta^(2j)*A^r_j,
    deg_t N_(i,j)<=b_j-i,
    b_j=10j+(6-d)*(d*(j-1)+2).

Higher scalar coefficients have valuation at least4j at t=0. The full-rank
zeroth cofactor bound, subtracting the basis row degrees, is

| Actual d | 2 | 3 | 4 | 5 | 6 |
|---|---:|---:|---:|---:|---:|
| Sufficient exclusion outside the pure-null locus | 18 | 50 | 93 | 132 | **140** |

All smaller ranks and row subsets containing0 are included. Outside
b=0,G=0, the retained complete shrinking/escape/flattened-Gaussian proof
excludes a constant normalized period and hence makes the zeroth
coefficient nonzero. Therefore the first140 zero pure moments force
b=0,G=0. These are exactly P=W*A(T), whose sharp mixed cutoff remains
deg Q. The new bound does not replace that analytic classification with
a generic saddle calculation.

## 3. The mixed-parent interior: use its actual seeds

For the normalized c*d!=0 interior of

    P=W*(a+bZ+cT+dZ^2)+G(Z,T),

retain the old phase in coordinate z,

    F=z/t-a-z^2, Phi=t*G(z,F)-F^2/2,
    I=constant*t^-1*integral exp(Phi) dz.

Write V=t^2*Phi. Every coefficient of V has t-degree at most3; the same
holds for S=t*V_t-2V-t^2. Actual phase degree d is4,5 or6, with polynomial
basis1,z,...,z^(d-2). Unreduced seeds now have R_j/t^(3j), with

    R_(j+1)=t^3*R_j_prime-3j*t^2*R_j+S*R_j,
    coefficient t-degree R_j<=3j, amplitude degree<=d*j.

Each leading cancellation has denominator A=d*V_d and coefficient
numerator degree at most3. Using r_j=d*(j-1)+2 gives numerator and
common-denominator degree b_j=3j+3*r_j. There is no row-grade subtraction
here: only parameter degree is being tracked. The resulting gates are
72,144,255 for d=4,5,6. Higher scalar valuations at least3j remove the
old additional derivative loss. The retained analytic proof makes the
zeroth coefficient nonzero throughout the interior. The c=0 and d=0
faces retain35/92, so255 suffices for the full fourteen-coefficient family.

## 4. Add back the frame branches and the global cover

Cusp and tangent D=0 branches have rank-one quadratic W coefficient.
If its T^2 coefficient is nonzero, section2 gives140 after the retained
Gaussian null rotation. Otherwise the mixed-parent255 gate applies.
Thus full cusp has max(225,140,255)=255 and full tangent has
max(284,140,255)=284. These are sufficient gates, without optimality.

The nondegenerate quadratic W branch remains336, including its trapped
component. The f=0 face is max(336,255)=336. Together with rank-one140,
all sixteen W-linear coefficients therefore have sufficient gate336.
The sharp general mixed cutoff3*deg Q is unchanged.

| Leading type, arbitrary lower terms | Full-class sufficient gate/exclusion |
|---|---:|
| Degree at most two | 7 |
| Smooth cubic | 94 |
| Irreducible nodal cubic | 353 |
| Irreducible cuspidal cubic | 255 |
| Transverse line-smooth-conic | 336 |
| Tangent line-smooth-conic | 284 |
| Three independent lines | 336 |
| Three distinct concurrent lines | 104 |
| Repeated leading line(s) | 99 |

The ordinary-node values retain the prior353/334/304 complement bounds
and conic336 frame exclusion. The maximum353 proves the global pure-only
test. For an odd polynomial all odd pure moments vanish automatically,
so the first352 include every needed condition through353. This parity
step is valid because353 is odd; it was unavailable for the older even590.
The separate all-homogeneous336 proof remains intact.

In the rational twenty-coefficient ring, radical(I353)=radical(I_infinity).
This is not ordinary ideal equality or an expanded cofactor computation.
At most ceil(log2(353-19))=9 exact monotone search rounds would determine
the least global K. No such full coefficient-ring computation is performed.
The original twenty/fourteen-moment witness searches are complete and are
not repeated.

## Evidence and remaining work

The checker verifies63 symbolic joint-degree weighted divergence identities,
replays all22 connection and8 direct seed identities in the four retained
full-Gaussian weighted models, checks the universal rank-one rescaling and
source identity, and certifies polynomial connections and the first three
actual seeds in eleven models. They include b=0, b!=0, actual degrees2--6,
the flattened Gaussian chart, the pure-null boundary, and mixed-parent
degrees4--6. The first run passed algebra and failed only at JSON emission
of a symbolic integer; its exact input and terminal failure are preserved.

Current arXiv version metadata for Crew--Feldbrugge2610.22767v2 and
Helmer--Tellander2505.12458v1 was rechecked. No new body review or theorem
import is claimed; their previously audited generic formulas do not replace
the analytic exclusions used here.

The ordinary-node353 complement is now the largest class bound. Improving
its actual source/reduction degrees or proving smaller sufficiency on
complete coefficient strata is justified next work. Exact leastness,
20/14 sufficiency, independent review, expanded global radical cofactors,
exceptional coefficient schemes, term-count minima and lower-power GVC
remain open.
