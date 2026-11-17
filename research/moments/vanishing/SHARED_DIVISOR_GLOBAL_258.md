# Shared divisors and weighted Euler seeds: a global258 pure prefix

15 November 2026. Internal written refinement of the
[exact Euler297 proof](EULER_REPRESENTATIVES_GLOBAL_297.md), with its full
original Wick functional, case cover and analytic exclusions retained.

For every complex polynomial P of total degree at most three in three
independent standard real Gaussian variables,

    E(P^j)=0 for j=1,...,258
      iff every positive pure moment vanishes.        (1)

No separate Gaussian-frame condition is required. Pure nullity then gives
the retained W-linear Gaussian frame and sharp mixed cutoff m>3*deg Q.
The least global and odd-only pure prefixes lie in **[20,258]**. The
all-homogeneous interval stays **[14,258]**, smooth-leading **[20,94]** and
homogeneous-smooth **[14,56]**. Both exact isolated lower witnesses remain
complete. Neither leastness nor20/14 sufficiency is proved.

The [checker](verify_shared_divisor_weighted_seeds.py),
[certificate](shared_divisor_weighted_seeds_certificate.json),
[receipt](shared_divisor_weighted_seeds_verification.json) and
[source audit](shared_divisor_weighted_seeds_source_audit.json) bind exact
auxiliary algebra. Complete rank/class coverage, all-coefficient induction
and analytic nonzero zeroth coefficients remain written mathematics.
Independent replay, external review, formal verification and expanded global
radical cofactors remain outstanding.

## 1. Keep the shared D factor in ordinary-node principal coefficients

Use the full ordinary-node normalization, including every affine Gaussian
coefficient and lower term. Its rescaled binary phase is

    W=C*x^2*y^2+A*x^3+B*y^3+E*x^2+F*y^2
         +J*xy+L*x+N*y+K0,
    C=1/2, delta=t^2*D, deg_t D<=1, D!=0.

The [retained exact normalization](ordinary_node_phase_shift_checkpoint.json)
proves the stronger coefficient identities

    A=a*D, B=d*D, E=D*Le, F=D*Lf,
    Le=t*b+3a*jv-q_uu/2,
    Lf=t*c+3d*ju-q_vv/2.                            (2)

Thus Le,Lf have t-degree at most one. The complete node/frame split
makes each needed principal coefficient nonzero. It excludes every member
having any isotropic ordinary node with F2=0 before applying this complement
calculation. Checking only one chosen node is still insufficient.

The old common principal divisors C*A*B, C*A*F and C*E*F repeat D. Inspecting
the complete principal divisions shows that the following suffice instead:

| Leading type | Required nonzero opens | New common divisor, up to constants | Maximum degree k |
|---|---|---|---:|
| Irreducible node | a*d*D!=0 | D | 1 |
| Transverse line-conic | a*D*Lf!=0, d=0 | D*Lf | 2 |
| Three independent lines | D*Le*Lf!=0, a=d=0 | D*Le*Lf | 3 |

For example, nodal x^3 reduction has basis coefficient B/A=d/a, with no
D denominator, while its divergence coefficient1/(3A) needs only D.
For x^i, i>=4, the additional coefficient B/(2AC) likewise cancels D.
Pure-y reduction needs1/B and hence a single D. The transverse x- and
y-axis divisions need1/(aD) and1/(D*Lf), whose common denominator is
D*Lf, not D^2*Lf. For triangles, F/E=Lf/Le and the other divisions have
denominators D*Le, D*Lf or constant multiples thereof. Their common
denominator is D*Le*Lf, not D^2*Le*Lf.

No coprimality, distinct-root or nonvanishing-at-zero condition is imposed
on D,Le,Lf. Coincident factors only make the common divisor conservative.
Nonzero constants a,d are elements of the coefficient field and do not
contribute t-degree.

Retain all four actual divergence operators from the earlier proof,
including the derivative of the fourth operator's coefficients. The new
cleared cofactors and remainders are the old ones divided by one shared
D and nonzero constants. They are polynomial in D,Le,Lf. Give those symbols
weight one, J weight2, L,N weight3, K0 weight4 and delta weight3.
The old complete projection estimates therefore become

    cleared cofactor joint degree<=M+k-3,
    cleared whole remainder joint degree<=M+k

for an input monomial of ordinary amplitude degree M. Substitution of
polynomials of degrees at most those weights preserves the bounds, also
when D is constant. Every non-basis correction still strictly lowers
amplitude degree. Equal-degree remainders are already basis terms.

## 2. Smaller ordinary-node complement gates

Use the unchanged exact radial Euler representative and its complete
actual-seed recurrence from the297 proof. It preserves t^-2*D^-1/2,
keeps the Euler amplitude derivative and volume divergence, and gives

    joint degree R_j<=6j, amplitude degree R_j<=3j,
    derivative denominator delta^(2j).

The same depths3j-2,3j-2,3j-1 now use principal divisor costs1,2,3.
For the unchanged full bases and their ordinary grades d_i,

    H_j=delta^(2j)*A0^r_j,
    deg_t N_(i,j)<=b_j-d_i,
    b_j=9j-2,12j-4,15j-3, respectively.

Higher scalar coefficients have valuation at least4j at t=0. The retained
minimal-seed/cofactor and analytic nonzero-zeroth argument gives

| Type | Full rank | Zeroth degree bound | Complement exclusion |
|---|---:|---|---:|
| Irreducible node | 7 | 9*28-2*7-11 | **227** |
| Transverse line-conic | 6 | 12*21-4*6-8 | **220** |
| Three independent lines | 5 | 15*15-3*5-6 | **204** |

Every smaller seed rank and selected row subset containing the unit row
is included. The conic frame branch retains258, so the full-class bound
for each ordinary-node type is258. The complement numbers are not silently
promoted to full-class gates.

## 3. Weighted Euler representatives cancel both weighted top terms

For cusp and tangent use the same delta=t^2*D, normalization t^-2*D^-1/2
and original source

    S=delta*W_t-delta_prime*W-2*t^3*D^2-t^4*D*D_prime/2.

Their top phases are x^4/2+D*y^3 with weights(3,4), and
x^4/2+D*x*y^2 with weights(2,3). Their lower phase weights are at most10
and7, respectively. Put

    beta_x=delta_prime/4,
    beta_y=2*t*D/3 (cusp),
    beta_y=t*D-delta_prime/8 (tangent).

For every amplitude f, add the exact image

    delta^-2*(beta_x*G_x(x*f)+beta_y*G_y(y*f)).

The same quotient connection then has representative

    nabla_beta(f)=f_t+(beta_x/ delta)*x*f_x
                         +(beta_y/delta)*y*f_y+(S_beta/delta^2)*f,
    S_beta=S+delta*(beta_x+beta_y)
                   +beta_x*x*W_x+beta_y*y*W_y.      (3)

The volume term delta*(beta_x+beta_y) is required. This is not a new
normalization or a logarithmic dy/y connection. The choice4*beta_x=delta'
cancels the x^4 source. For cusp,3*beta_y=2*delta/t; for tangent,
beta_x+2*beta_y=2*delta/t. Together with delta=t^2D these cancel the y^3
or xy^2 source exactly. Thus

    joint degree S_beta<=6,
    weighted degree S_beta<=10 (cusp), <=7 (tangent).

The full actual-seed recurrence is

    R_(j+1)=delta^2*R_j_prime-2j*delta*delta_prime*R_j+S_beta*R_j
               +delta*(beta_x*x*(R_j)_x+beta_y*y*(R_j)_y). (4)

Its difference from the old recurrence on the same R_j has exact witnesses
beta_x*x*R_j and beta_y*y*R_j. Keeping both derivative terms proves that
the original moment series is still being differentiated.

## 4. Weighted complement and full-class bounds

From(4), joint degree R_j<=6j, while weighted degrees are at most10j and7j.
The unchanged reductions lower weight by2 or1, with minimum reducible
weights8 and5. Consequently the sufficient depths decrease to

    cusp: r_j=5j-3;
    tangent: r_j=7j-4.

The common divisors remain D and D^2, up to nonzero constants, with
joint costs1/2. Padding unused reduction depths keeps each row estimate:

    H_j=delta^(2j)*A0^r_j,
    deg_t N_(i,j)<=b_j-d_i,
    cusp: b_j=11j-3;
    tangent: b_j=20j-8.

Higher scalar valuations remain at least4j. With basis grade sums9 and6,
the retained cofactor argument gives complement exclusions

    cusp: 11*21-3*6-9=204;
    tangent: 20*15-8*5-6=254.

All smaller actual seed ranks and unit-row subsets are included. The D=0
frame branches use retained rank-one140 or mixed-parent255, so the full
cusp and tangent gates are both255. No optimality or new mixed-cutoff
sharpness within either leading subclass is claimed.

## 5. The global258 cover

| Leading type, arbitrary lower terms | Full-class sufficient gate/exclusion |
|---|---:|
| Degree at most two | 7 |
| Smooth cubic | 94 |
| Irreducible nodal cubic | max(227,258)=258 |
| Irreducible cuspidal cubic | max(204,140,255)=255 |
| Transverse line-smooth-conic | max(220,258)=258 |
| Tangent line-smooth-conic | max(254,140,255)=255 |
| Three independent lines | max(204,258)=258 |
| Three distinct concurrent lines | 104 |
| Repeated leading line(s) | 99 |

Every coefficient point is assigned by the retained complete leading-type
and frame cover. Its maximum258 proves(1). Since258 is even, parity does
not lower the odd-only bound to257. The completed odd witness20 still
gives the lower bound20; homogeneous258 and its witness14 are unchanged.

In the rational twenty-coefficient ring,
radical(I258)=radical(I_infinity), without ordinary ideal equality or
expanded radical-membership cofactors. At most ceil(log2(258-19))=8 exact
monotone search rounds against I258 would determine the least global K.
No such full coefficient-ring search or high-order moment expansion is run.

## Evidence, sources and remaining work

The checker verifies135 symbolic shared-divisor divergence identities and
four universal axis factorizations. It replays all36 connection and18 seed
certificates from six retained full-Gaussian ordinary-node models after
exact polynomial division by the redundant factor. It also verifies30
symbolic weighted Euler connection identities and22 connection/12 seed
certificates in four full-Gaussian weighted models. Both regular and
D(0)=0 cases are included. The initial run passed; no failure is discarded.

A fresh literature search and primary arXiv version checks of
[Long's Gaussian counterexamples](https://arxiv.org/abs/2607.18186),
[Wilson's two-variable proof](https://arxiv.org/abs/2607.23887), and
[Long's compact-group result](https://arxiv.org/abs/2609.16178) found their
retained v1 versions. This is a metadata/abstract check, not a fresh body
review or a claim of exhaustive literature coverage. No new external
theorem or producer implementation is imported. Prior analytic assurance
limits remain in force.

The [historical proof index](PREFIX_PROOF_HISTORY.md) identifies the retained
coarser numerical frontiers. Hash-bound353 and297 proof bodies stay unchanged.
Current mathematical authority is the ledger and this note. The current
bottlenecks are conic258 and mixed-parent255. Exact leastness,20/14
sufficiency, expanded global radical certificates, independent review,
exceptional coefficient schemes, term-count minima and lower-power GVC remain
open; neither lower-witness search nor the qualitative cubic case cover is
reopened.
