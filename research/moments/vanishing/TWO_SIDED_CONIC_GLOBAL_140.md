# Two-sided conic degree bounds: global140, homogeneous130

15 November 2026. Internal written refinement of the [normal-module146 proof](NORMAL_MODULE_GLOBAL_146.md). The original Wick functional, full quotient bases, logarithmic volume normalization, complete coefficient/frame cover and analytic exclusions are retained.

For every complex ternary polynomial P of total degree at most three,

    E(P^j)=0 for j=1,...,140 iff all positive pure moments vanish.

No extra Gaussian-frame hypothesis is required. Pure nullity gives the retained W-linear frame and sharp mixed cutoff m>3*deg Q. Least global and odd-only prefixes lie in **[20,140]**; the all-homogeneous interval is **[14,130]**. Smooth-leading **[20,94]** and homogeneous-smooth **[14,56]** stay separate. The first-nonzero20/14 witnesses are complete. Both upper bounds are even; parity gives no139/129 reduction. Exact leastness and20/14 sufficiency remain open.

The [checker](verify_two_sided_conic_seed_bounds.py), [certificate](two_sided_conic_seed_bounds_certificate.json), [receipt](two_sided_conic_seed_bounds_verification.json) and [source audit](two_sided_conic_seed_bounds_source_audit.json) support auxiliary exact algebra. All-degree degree/module/path inductions, complete class/rank coverage and analytic nonzero zeroth coefficients remain written inputs. There is no independent replay, external review, formal verification or expanded global radical certificate.

## 1. Use both endpoint coefficient bounds

The retained structured cubic conic phase is

    W=A0+sum_(k=1)^3(A_k*y^k+K^k*B_k*y^(-k)),
    deg A0<=3, deg A_k,deg B_k<=3-k, deg K<=2.

If deg K0, both endpoints obey deg W_k<=3-|k|. If deg K1, they obey deg W_k<=3-max(k,0). These strengthen the old one-sided bound3-k. The functions |k| and max(k,0) are subadditive, so their coefficient bounds propagate through products without selecting generic coefficients.

For actual negative order q3, retain the previous negative Euler representative

    B=t*K_prime-2*K/3, h=t^3*K,
    S_K=K*(t*W_t-2W-t^2)+B*y*W_y,
    R_(j+1)=h*R_j_prime-j*h_prime*R_j
                           +t^2*B*y*(R_j)_y+S_K*R_j.

The source has support[-2,p]. All denominator and vertical derivative terms are retained. For deg K0 the full recurrence gives deg R_(n,j)<=3j-|n|. For deg K1 it gives deg R_(n,j)<=4j-max(n,0). Negative coefficients remain divisible by K^(-n). No new singular integral or volume factor is introduced.

## 2. Sharpen the q3 high/low paths

Keep the full basis y^-3,...,y^(p-1). High and low paths do not feed outside one another's retained boundary. With p3, positive leading A3 is a nonzero constant; low division by K^3 is polynomial in the retained negative module. With p2, clear only A2^(2j-1) at the end. Shared factors among A2,K,t are permitted by the same separate-path proof.

For p3 and constant K, a high path ending in row k has degree<=3j-k; a low path has degree<=3j+k. The second bound follows by starting at n<0 with degree<=3j+n and adding at most k-n along a low path. Therefore every row has degree<=3j+|k|. The full-rank gate is3*21+9=72.

For p3 and linear K, high paths have degree<=4j-k. Each low cancellation divides K^3, of degree3, and every nonleading phase coefficient has degree<=3; hence no low step increases coefficient degree. Low rows have degree<=4j. Thus the uniform row bound4j+max(-k,0) gives4*21+6=90.

Along rational endpoint paths, coefficient degree means numerator degree minus denominator degree at infinity. Once the bounded common endpoint power is cleared, the normal coefficients are polynomial and this becomes their ordinary degree.

For p2, split by the actual degree of A2 rather than assuming it is linear. If A2 is linear, raw high coefficient growth is at most exponent displacement. After common A2^(2j-1) padding, the constant-K row bound is5j-1+|k|, and the linear-K bound is6j-1+max(-k,0). If A2 is constant, raw high growth is at most twice the displacement, and the corresponding safe bounds are5j-2k and6j-2k. At j1 there is no low input; this matters for the constant-K constant-A2 bound at the positive row. Every higher seed and every smaller actual rank is included.

| p2,q3 | deg A2=0 | deg A2=1 |
|---|---:|---:|
| Constant K | 85 | 77 |
| Linear K | 100 | 91 |

The unchanged quadratic-K bounds are108 for p3 and105 for p2. The retained p<=1,q3 bounds remain at most130. Higher scalar coefficients retain h^j times the stated endpoint powers, hence valuation>=3j.

## 3. The remaining p3,q2 case: positive Euler and a smaller module denominator

The only old q<3 gate above108 was p3,q2 at145. Here A3 is a nonzero constant and the negative endpoint is K^2*B2, with nonzero B2 of degree<=1. Use the exact positive Euler representative

    nabla_+(f)=f_t+(2/(3t))*y*f_y
                  +(t*W_t-2W-t^2+(2/3)*y*W_y)*f/t^3.

Its difference from the original normalized connection is(2/(3t))*Dlog(f), where Dlog=t^(-2)*(t^2*y*partial_y+y*W_y). The logarithmic measure dy/y has no additional ordinary-volume divergence constant. This cancels the positive y^3 source coefficient.

Clear the whole operator with h=t^3*K, and put

    S_+=K*(t*W_t-2W-t^2+(2/3)*y*W_y),
    R_(j+1)=h*R_j_prime-j*h_prime*R_j
                           +(2/3)*t^2*K*y*(R_j)_y+S_+*R_j.

The source support is[-2,2], and seeds have support[-2j,2j]. Multiplication by K in h restores any K factor lost under time differentiation. The negative coefficient module is therefore retained. High division by3*A3 introduces no denominator. Low division by -2*K^2*B2 removes the K^2 already in the module, leaving only B2 denominators. At most2j-2 high and2j-2 low steps are needed. Thus the normal clearing denominator is just B2^(2j-2).

For this module argument, allow B2 inverses formally while retaining K as an independent divisor. Clear the bounded B2 power before substituting actual coefficients. It remains valid if actual B2 shares a root with K or t; no coprimality open is imposed. Exact divergence witnesses stay in C(t), and every expanded identity uses the original logarithmic twisted operator.

## 4. All K-degree strata of p3,q2

When K is constant, the new source has the two-sided bound3-|k| and gives seed coefficients3j-|k|. If B2 is linear, raw low degree growth is at most upward exponent displacement; high paths obey the matching downward bound. Adding B2^(2j-2) padding gives normal row degree<=5j-2+|k|, hence71 on rows-2,...,2. If B2 is constant, low growth is bounded by twice the displacement; the safe uniform bound5j+2|k| gives87. All smaller rank/row choices containing the unit are included.

When K is linear, use weights(t,y)=(2,1). The phase has weight<=6 and source/seed cost8. Each high step has whole cost<=3, each low step after K^2 division and B2 clearing has cost<=4. Padding a completed row by B2 costs at most2, below4. Thus

    2*deg N_(k,j)+k<=8j+3*(2j-2)+4*(2j-2)=22j-14.

The full determinant degree is129 after rounding each polynomial row bound down individually. The coarser sum-before-rounding estimate130 is valid but not the value asserted by the exact row arithmetic. If K is quadratic, weights(1,1) give phase weight3, seed cost5, high cost0 and low whole cost1, so deg N_(k,j)<=7j-2-k and the gate is95. Both B2-degree strata are covered by these latter estimates. Every new H_j=h^j*B2^(2j-2) has valuation>=3j and obeys the same safe degree bounds.

The complete conic cover is therefore

    p3,q3, deg K0/1/2: 72/90/108;
    p2,q3, deg K0/1/2: <=85/100/105;
    p<=1,q3: <=130;
    p3,q2, deg K0/1/2: <=87/129/95;
    every other q<3 case: <=108.

This proves **uniform conic130** for cubic G. It does not cover arbitrary-degree G. The full W-linear family has gate max(conic130,parent92,rank-one140)=**140**. Homogeneous smooth56,cusp10,reducible10,nonisotropic-node72 and isotropic-node conic130 give **homogeneous130**. Even130 gives no129 parity deduction.

## 5. Complete global140 cover and assurance

| Leading type with arbitrary lower terms | Full sufficient gate/exclusion |
|---|---:|
| Degree at most two | 7 |
| Smooth cubic | 94 |
| Irreducible node | max(129,130)=130 |
| Irreducible cusp | max(131,92,140)=140 |
| Transverse line-smooth-conic | max(132,130)=132 |
| Tangent line-smooth-conic | 132 |
| Three independent lines | max(118,130)=130 |
| Three distinct concurrent lines | 104 |
| Repeated leading line(s) | 99 |

The maximum140 proves the global pure-only test; odd140 is unchanged by parity. In the rational twenty-coefficient ring radical(I140)=radical(I_infinity), without ordinary ideal equality or expanded global cofactors. At most ceil(log2(140-19))=7 exact monotone radical-search rounds against I140 would find the least global K; none has run.

The checker verifies seven universal source strata,27 positive logarithmic Euler identities, endpoint degree/all-rank arithmetic and69 direct seed identities in23 models. Three models are realized Gaussian conic coordinate phases; synthetic enlarged structured cases include constant/linear endpoints, simple/double K roots and common endpoint/K/t factors. Full all-degree module/path and degree proofs, complete coefficient/rank coverage and the analytic nonzero zeroth arguments remain written inputs. The scalar connection and original Wick functional are unchanged, so those analytic exclusions remain applicable.

The first run passed all source and model identities, then stopped because the expected linear-K p3q2 gate was130 while individual integer row rounding gives129. The failed input and terminal receipt are retained. The accepted run records129; the global/conic maxima are unchanged by this arithmetic correction. No previous proof body or receipt is rewritten, and no new external theorem or producer package is imported. The latest selected-source refresh remains the preceding primary-source audit.

The controlling bound is now unrestricted rank-one140; the conic family has130. Exact leastness,20/14 sufficiency, expanded I140 radical cofactors, independent analytic review, exceptional coefficient schemes, global term minima and lower-power GVC remain open. Qualitative cubic coverage and both lower witnesses stay closed internally.
