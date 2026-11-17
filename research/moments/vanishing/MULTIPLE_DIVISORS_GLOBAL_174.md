# Several divisor filtrations: global174 and homogeneous146

15 November 2026. Internal written refinement of the [pole-filtered199 proof](POLE_FILTERED_GLOBAL_199.md). The original Wick functional, full quotient bases, complete coefficient/frame cover, period comparison and analytic nonzero zeroth exclusions are retained.

Every complex ternary polynomial P of total degree at most three satisfies

    E(P^j)=0 for j=1,...,174 iff all positive pure moments vanish.

No extra Gaussian-frame hypothesis is required. Pure nullity gives the retained W-linear frame and sharp mixed cutoff m>3*deg Q. Least global and odd-only prefixes lie in **[20,174]**; the all-homogeneous prefix lies in **[14,146]**. Smooth-leading **[20,94]** and homogeneous-smooth **[14,56]** are unchanged. The first-nonzero20/14 witnesses are complete. Both new upper bounds are even, so parity does not lower them by one. No leastness or20/14 sufficiency is claimed.

The [checker](verify_multiple_divisor_seed_bounds.py), [certificate](multiple_divisor_seed_bounds_certificate.json), [receipt](multiple_divisor_seed_bounds_verification.json) and [source audit](multiple_divisor_seed_bounds_source_audit.json) support exact auxiliary algebra. Universal all-degree inductions, complete class/rank coverage and analytic exclusions remain written inputs. No independent replay, external review, formal verification or expanded global radical cofactors are claimed.

## 1. Work with independent divisors before specializing them

The ordinary-node phase has the retained form

    W=u^2*v^2/2+D*(a*u^3+d*v^3+Le*u^2+Lf*v^2)
                                      +J*u*v+L*u+N*v+K,
    delta=t^2*D, h=t^3*D^2.

For nodal type a,d are nonzero constants; for transverse type d=0,a!=0; for triangle type a=d=0. The retained principal opens are A0=D, D*Lf and D*Le*Lf, respectively. Le,Lf are nonzero on the types requiring them. The coefficients have joint grades D,Le,Lf=1, J=2, L,N=3, K=4, giving joint degree W<=4.

Initially treat t,D,Le,Lf as independent variables of grade1, and each other coefficient/derivative jet as an independent variable with its retained degree bound. The formal time derivative maps each linear divisor to its constant derivative and lowers its degree by one. This represents the actual polynomial phase after specialization, including constant divisors and coincident roots. Ordinary principal division is polynomial in this ring after clearing the same A0, with constants a,d invertible where already required.

Several separate valuation estimates can now be combined: distinct formal variables are relatively prime in this polynomial ring, so their forced powers divide simultaneously. Only after this product divisibility and its graded degree bound are established do we substitute actual linear polynomials in t. The specialized divisors may be proportional, share roots with each other or with t, or become nonzero constants. No coprimality assertion about those specialized divisors is used. Dividing by a specialized constant alone would not prove the degree saving.

## 2. Fractional axis grades complement the D grade

For any one divisor F and nonnegative rational amplitude weights (rho_u,rho_v), define pole grade

    pole_F(t^a*F^b*u^i*v^k)=rho_u*i+rho_v*k-b.

All other independent coefficients have pole grade0. Time differentiation raises this grade by at most1. Use the following grades on the indicated complete types:

| Type | Divisor | rho_u,rho_v | Actual seed pole bound |
|---|---|---|---|
| All three ordinary types | D | 1,1 | 2j |
| Transverse | Lf | 1/3,1 | 2j |
| Triangle | Le | 1,1/2 | 2j |
| Triangle | Lf | 1/2,1 | 2j |

The complete common-time radial Euler source is

    Sbar=t*D*W_t-(2D+t*D_prime)*(W-E(W)/4)-t^2*D^2,
    E=u*partial_u+v*partial_v.

It has grade<=2 for every listed filtration. On the transverse/triangle types the absent cubic axes are essential to the axis estimates. The Euler derivative coefficient delta*delta_prime/(4t) has grade<=2; h times a derivative and h_prime have no larger cost. Therefore the full actual recurrence, with its denominator derivative and vertical Euler field retained, gives grade<=2j by induction.

For an input exponent e=(i,k), put rho(e)=rho_u*i+rho_v*k. The unchanged A0-cleared principal formulas give, for each listed divisor separately,

    pole_F(remainder), pole_F(Cu), pole_F(Cv)<=rho(e)-1.  (1)

The principal formulas have finitely many exponent branches; adding extra monomial powers adds their amplitude grades on both sides, and differentiated amplitudes only improve the bound. Thus checking each explicit branch proves(1) for arbitrary exponents. The checker additionally expands all identities through ordinary degree12. Nonbasis remainders strictly drop ordinary amplitude degree; equal-degree normal terms lie in the retained basis. Padding a completed basis row by A0 also lowers each pole grade by one, since each listed divisor occurs once in A0.

At the retained depth r_j, (1) gives normal row divisibility

    F^(r_j-2j+ceil(rho(e))).

The witness bound in(1), propagated with its A0 padding, gives common witness divisibility F^(r_j-2j). Negative displayed exponents at the first seed mean multiplication, rather than removal of a polynomial factor. All normalized numerators and witnesses are polynomial after this operation, as follows from the valuation inequalities; no singular specialization is used to infer it.

A tempting integral transverse grade, v exponent minus Lf exponent, fails on input u^3: the cleared normal contains Lf^2*v^2, of pole grade0, while the required bound is-1. The fractional u weight1/3 corrects exactly this obstruction. The accepted certificate retains that counterexample; the failed grade is not used.

## 3. Combine the formal valuations and take cofactors

For ordinary types put b_j=r_j-2j. Cancel the product of divisor powers F^b_j from the completed numerators, denominators and exact divergence witnesses. Rescale each completed row by the product F^ceil(rho_F(e)); the actual seeds have already been computed. We do not differentiate a new time-dependent basis and omit its derivative terms.

| Type | Depth r_j | Common factor | New H_j | Row degree bound |
|---|---|---|---|---|
| Node | 3j-2 | D^(j-2) | t^(3j)*D^(4j) | 7j-2(i+k) |
| Transverse | 3j-2 | (D*Lf)^(j-2) | t^(3j)*D^(4j)*Lf^(2j) | 9j-2(i+k)-ceil(i/3+k) |
| Triangle | 3j-1 | (D*Le*Lf)^(j-1) | t^(3j)*D^(4j)*Le^(2j)*Lf^(2j) | 11j-2(i+k)-ceil(i+k/2)-ceil(i/2+k) |

These follow from the old joint bounds5j+c*r_j-|e|, c=1,2,3, subtracting the formal degree of every cancelled common and row factor before specialization. The full row-grade sums are22,23,23 at ranks7,6,5. The zeroth cofactor bounds are

    node:       7*28-22=174,
    transverse: 9*21-23=166,
    triangle:  11*15-23=142.

Every smaller actual seed rank and row subset containing the unit has a smaller bound. Higher scalar coefficients contain the new H_j with valuation>=3j, so the original initial-order contradiction loses no extra derivative order. Analytic nonzero zeroth coefficients come from the same rational quotient seed and retained exclusion, not from the polynomial factor checks alone.

## 4. The cusp has a D-pole filtration as well

Use the common-time weighted Euler connection, phase W=u^4/2+D*v^3+R with old weights(3,4), lower weight<=10 and joint degree<=4. Work first with independent formal D and the compatible derivative as above. The grade

    pole_D=u exponent+2*v exponent-D exponent

gives source grade<=4 and actual seed grade<=4j. The weighted Euler derivative coefficients and denominator derivative respect that induction. The unchanged D-cleared principal divisions, using W_u leading2u^3 and W_v leading3D*v^2, satisfy (1) with rho=i+2k. The old weighted termination depth is r_j=5j-3.

Thus row coefficients contain D^(j-3+rho), and witnesses contain D^(j-3). Cancel D^(j-3), scale completed rows by D^rho, and get

    H_j=t^(3j)*D^(6j),
    deg N_(e,j)<=9j-|e|-rho(e).

For basis1,u,u^2,v,u*v,u^2*v the combined row grades are0,2,4,3,5,7, summing21. At full rank6 the sufficient complement gate is9*21-21=168; every smaller seed rank is covered. Parent92 and rank-one140 frame branches are below168, giving **full cusp168**. Constant D and D(0)=0 are included by the formal argument; no coefficient open is added.

## 5. Conic displacement estimates give full146

Retain the structured Laurent phase, negative module and Euler recurrence from the preceding proofs. High and low paths do not feed outside each other's retained boundary: high cancellation at n>=p produces exponents at least-3; low cancellation at n<-3 produces exponents at most p-1. The unchanged basis is y^-3,...,y^(p-1). Denominator and volume derivative terms remain in the actual seeds.

For p=q=3 and deg K1, use weights(t,y)=(2,1). The unreduced seed has weighted degree<=8j. A high step displaces its exponent downward; its weighted cost is at most the magnitude of that displacement. A low step, after its polynomial K^3 division, has cost at most its upward displacement. The endpoint coefficients are nonzero constants after removing K^3, so no additional denominator is needed.

Starting a high path at n<=3j and ending in row k gives weighted degree<=11j-k, hence coefficient degree<=floor(11j/2)-k. Starting a low path at n>=-2j gives weighted degree<=10j+k, hence coefficient degree<=5j. Both are at most floor(11j/2)+3 for every row -3,...,2. The denominator h^j has degree4j and is below this bound. Summing these safe row-independent bounds at full rank6 gives **132**, and every smaller rank is below132. This deliberately conservative estimate is sufficient here.

For p2,q3 and constant nonzero K, the positive coefficient A2 has degree<=1. Retain only A2^(2j-1) as common clearing denominator; no low endpoint padding is needed. High paths have the preceding joint-degree bound5j-1-k. On low paths the coefficient degree increase is at most exponent displacement; starting with deg R_(n,j)<=3j-n and n>=-2j gives degree<=7j+k before positive-endpoint padding. Including A2^(2j-1) gives9j-1+k. At j1 no low input exists, so its bound is4-k. For j>=2 the low bound dominates the high one on all rows-3,...,1. The full determinant bound is

    4+sum_(j=2)^5(9j-1)+sum_rows k-2*min_rows k=127.

All smaller ranks/row choices containing0 are checked. Shared A2/K/t factors remain allowed by the retained separate-path proof. The denominator has degree<=5j-1 and valuation>=3j.

The complete conic endpoint cover is now

    p3q3, deg K0/1/2: 146/132/108,
    p2q3, deg K0/1/2: 127/136/105,
    p<=1,q3: <=130; q<3: <=145.

Thus **conic146** is uniform for cubic G. Parent92 and unrestricted rank-one140 give **full W-linear146**. This is not an arbitrary-degree G claim. The homogeneous cover smooth56,cusp10,reducible10,nonisotropic node72,isotropic-node conic146 gives **homogeneous146**. Since146 is even, no145 parity reduction follows.

## 6. Global174 and evidence limits

| Leading type with arbitrary lower terms | Full sufficient gate/exclusion |
|---|---:|
| Degree at most two | 7 |
| Smooth cubic | 94 |
| Irreducible node | max(174,146)=174 |
| Irreducible cusp | 168 |
| Transverse line-smooth-conic | max(166,146)=166 |
| Tangent line-smooth-conic | 132 |
| Three independent lines | max(142,146)=146 |
| Three distinct concurrent lines | 104 |
| Repeated leading line(s) | 99 |

The maximum174 proves the global pure-only test. The odd gate also stays174. In the rational twenty-coefficient ring radical(I174)=radical(I_infinity), without ordinary ideal equality or expanded global cofactors. At most ceil(log2(174-19))=8 monotone exact radical-search rounds against I174 would determine the least global K; no full coefficient-ring round has run.

The checker verifies formal source and principal/divergence filtrations, including fractional grades and the rejected integral-grade counterexample. It exactly replays24 binary seed identities in eight full-Gaussian models and60 Laurent seed identities in twenty retained conic models. It also checks individual linear-K displacement costs and all determinant/rank/row arithmetic. The accepted run completes in about25 seconds. This is a new auxiliary checker over retained certificate inputs, not independent reproduction from a separate implementation or formal verification of the analytic proof.

All earlier proof bodies and receipts remain unchanged. No fresh literature-coverage claim, external theorem or producer package is imported in this algebraic refinement; the retained primary-source audit remains in force. The controlling bounds are now ordinary node174 and cusp168. Exact leastness,20/14 sufficiency, expanded global radical cofactors, independent analytic review, exceptional schemes, global term minima and lower-power GVC remain open. Qualitative cubic coverage and both lower witnesses stay closed internally.
