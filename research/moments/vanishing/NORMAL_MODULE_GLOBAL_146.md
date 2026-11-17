# Polynomial normal modules: global146

15 November 2026. Internal written refinement of the [multi-divisor174 proof](MULTIPLE_DIVISORS_GLOBAL_174.md), retaining the original Gaussian Wick functional, quotient bases, complete coefficient/frame cover, period comparison and analytic exclusions.

For every complex ternary polynomial P of total degree at most three,

    E(P^j)=0 for j=1,...,146 iff all positive pure moments vanish.

No extra Gaussian-frame hypothesis is imposed. Pure nullity gives the retained W-linear frame and sharp mixed cutoff m>3*deg Q. Least global and odd-only prefixes lie in **[20,146]**, and the all-homogeneous interval stays **[14,146]**. Smooth-leading **[20,94]** and homogeneous-smooth **[14,56]** remain separate. The first-nonzero20/14 witnesses are complete. Even146 gives no145 parity reduction. Leastness and20/14 sufficiency remain open.

The [checker](verify_ordinary_normal_seed_module.py), [certificate](ordinary_normal_seed_module_certificate.json), [receipt](ordinary_normal_seed_module_verification.json) and [source audit](ordinary_normal_seed_module_source_audit.json) support exact auxiliary algebra. Universal all-degree module/filtration inductions, complete class/rank coverage and analytic nonzero zeroth coefficients remain written inputs. No independent replay, external review, formal verification or expanded global radical cofactors are claimed.

## 1. A structural factor left unused in the previous bound

The rescaled, translated ordinary phase has the retained form

    W=u^2*v^2/2+D*(a*u^3+d*v^3+Le*u^2+Lf*v^2)
                                      +J*u*v+L*u+N*v+K,
    delta=t^2*D, h=t^3*D^2.

In the actual full-Gaussian normal form, both L and N are divisible by D. This is universal, not a generic model observation. The retained symbolic phase-shift checkpoint gives

    L/D=-A*jv-C*ju+3*a*jv^2+2*b*jv*t+e*t^2-lu*t,
    N/D=-B*ju-C*jv+2*c*ju*t+3*d*ju^2+f*t^2-lv*t.

Both quotients have degree<=2 in t. The letters A,B,C here are the original quadratic-form coefficients in that checkpoint, not the reduction divisor. The identity includes arbitrary affine Gaussian means and mixed terms, and all nondegenerate metrics already covered there. It remains true for constant D, D(0)=0 and the lower cubic-axis specializations.

Write L=D*L0,N=D*N0. Work first with independent graded D,Le,Lf of degree1, J of degree2, L0,N0 of degree2 and K of degree4. This preserves the joint phase bound4. The earlier enlargement that allowed unrelated L,N was valid for coarser bounds, but unnecessary for these actual phases.

## 2. Polynomial normal coefficients do not require polynomial witnesses

Use the unchanged complete ordinary principal formulas and bases. Previously their common cleared divisors were D,D*Lf,D*Le*Lf for node, transverse and triangle types. After substituting L=D*L0,N=D*N0, every normal remainder has one more factor D. Dividing that whole one-step identity by D gives an exact reduction with normal clearing divisors

    node:1, transverse:Lf, triangle:Le*Lf.

The normal remainders are polynomial in the independent parameter ring. This follows directly from the explicit principal branches. In pure-axis reductions the only potentially unpaired D denominators in the normal coefficients multiply L or N; the displayed universal factors remove them. Mixed reductions already have polynomial normal coefficients. The u^3 to v^3 equal-degree normal relation divides only the already nonzero constant a. All other nonbasis terms strictly drop ordinary amplitude degree. Thus the reduction terminates on the retained full bases without D padding in their normal coefficients.

The exact divergence witnesses can still contain 1/D. They are retained as rational time-coefficient amplitudes in the same quotient over C(t). The resulting identity is

    Anew^r*R-normal
       =delta*(partial_u Cu+partial_v Cv)+W_u*Cu+W_v*Cv,
    Cu,Cv in C(t)[u,v].

There is no requirement that these witnesses be polynomial in t for an equality in the rational quotient. We do not cancel D from the whole connection, alter the time/volume normalization, or assert that the quotient at D=0 is a polynomial fibre. The scalar normal coefficients, whose polynomial degrees control the cofactor bound, are distinguished from the witnesses. Their rational equality is checked using the full original twisted divergence operator.

At the same depths3j-2,3j-2,3j-1, each completed normal numerator is therefore divisible by D^r. The witnesses divided by D^r remain rational and must not be silently discarded. The exact checker retains and expands them.

## 3. Keep the earlier axis filtrations and reduce the complement gates

Actual seed numerators R_j still have joint degree<=5j and amplitude degree<=3j, with the complete common-time Euler recurrence and h=t^3*D^2 unchanged. On the transverse/triangle types, retain the previous independent Le/Lf pole filtrations with fractional amplitude weights and seed bound2j. Their normal row factors combine with D^r formally before specializing actual linear divisors, so coincident roots and constant divisors require no new coprimality open.

Cancel D^r and the common axis powers, then rescale the already computed seed rows as follows:

| Type | Additional common factor after D^r | Completed row scale | New row degree bound |
|---|---|---|---|
| Node | 1 | 1 | 5j-i-k |
| Transverse | Lf^(j-2) | Lf^ceil(i/3+k) | 7j-i-k-ceil(i/3+k) |
| Triangle | (Le*Lf)^(j-1) | Le^ceil(i+k/2)*Lf^ceil(i/2+k) | 9j-i-k-ceil(i+k/2)-ceil(i/2+k) |

Negative early exponents mean multiplication. All normal coefficients are polynomial after these operations. Degree savings are established in the independent formal grading before substitution, which includes constant divisors. We do not differentiate a new rescaled basis while omitting its derivative terms: the actual seeds were computed first.

The new denominators are respectively

    t^(3j)*D^(2j),
    t^(3j)*D^(2j)*Lf^(2j),
    t^(3j)*D^(2j)*Le^(2j)*Lf^(2j).

They have degree bounds5j,7j,9j and valuation>=3j. The full basis row-grade sums are11,15,17 at ranks7,6,5. Their zeroth cofactor bounds are

    node:       5*28-11=129,
    transverse: 7*21-15=132,
    triangle:  9*15-17=118.

Every smaller actual seed rank and row subset containing the unit is included. Higher scalar coefficients contain the displayed H_j, so the retained initial-order contradiction loses no additional derivative order. The same analytic exclusion supplies a nonzero zeroth coefficient, since each actual rational seed class is unchanged.

## 4. Halving the cusp pole weights gives complement131

Use the common-time cusp phase W=u^4/2+D*v^3+R with old weights(3,4), lower weight<=10 and joint degree<=4. In the formal independent-D ring use pole grade

    u exponent/2+v exponent-D exponent.

The complete weighted Euler source has grade<=2. Both denominator-derivative terms and vertical Euler fields respect the actual-seed induction, giving grade<=2j. The unchanged D-cleared principal remainder and its witnesses have grade<=rho(e)-1, where rho(e)=i/2+k. For arbitrary exponents this follows from the two principal branches and the lower phase bounds; the checker additionally expands all identities through old weight30.

At depth r_j=5j-3, the normal row coefficient is divisible by

    D^(3j-3+ceil(i/2+k)),

and the exact witnesses contain the common factor D^(3j-3). Cancel the common factor and scale completed rows by D^ceil(i/2+k). The new denominator is t^(3j)*D^(4j), with degree<=7j and valuation>=3j. Its row bounds are

    deg N_(i,k;j)<=7j-i-k-ceil(i/2+k).

For basis1,u,u^2,v,u*v,u^2*v the row grades are0,2,3,2,4,5, summing16. Thus the complement gate is7*21-16=131, including every smaller seed rank. The retained parent92 and unrestricted rank-one140 frame branches give **full cusp140**. The generic rank-one140 bound is unchanged.

## 5. Global146 and the remaining work

The conic/full-W-linear146 result from the preceding proof now controls the maximum:

| Leading type, arbitrary lower terms | Full sufficient gate/exclusion |
|---|---:|
| Degree at most two | 7 |
| Smooth cubic | 94 |
| Irreducible node | max(129,146)=146 |
| Irreducible cusp | max(131,92,140)=140 |
| Transverse line-smooth-conic | max(132,146)=146 |
| Tangent line-smooth-conic | 132 |
| Three independent lines | max(118,146)=146 |
| Three distinct concurrent lines | 104 |
| Repeated leading line(s) | 99 |

The maximum146 proves the global pure-only test, without a separate frame assumption. Odd and homogeneous gates also equal146; evenness gives no145 deduction. In the rational twenty-coefficient ring radical(I146)=radical(I_infinity), without ordinary ideal equality or expanded global cofactors. At most ceil(log2(146-19))=7 exact monotone radical-search rounds against I146 would determine the least global K; no full coefficient-ring round has run.

The checker verifies two universal axis factors,273 exact ordinary principal normal/divergence identities, the half-weight cusp source and principal bounds, and24 actual seed identities in eight retained full-Gaussian models. Ordinary witnesses are checked over QQ(t), explicitly permitting their D poles; every normal coefficient is checked polynomial with its stated degree bound. The accepted run takes about15 seconds. These are auxiliary checks of a new written internal refinement, not independent reproduction or formal verification of the complete analytic proof. All earlier proof bodies and receipts stay unchanged.

The fresh targeted literature search is not exhaustive. Primary arXiv records for Long2607.18186 and Dvorsky2608.07338 still list v1; their results are already retained. Selected reading of Wang--Zhang2610.16157v1, Introduction equation(1.1) and Theorem statements1.1/1.3, strengthens the existing scope check: it studies high-frequency spatial eigenfunction moments, rather than exact Gaussian polynomial pure/mixed moment identities. No theorem is imported into this polynomial argument. See the linked source audit for exact URLs and read coverage.

The largest remaining sufficient bounds are conic146 and unrestricted rank-one140. Exact leastness,20/14 sufficiency, expanded I146 radical cofactors, independent analytic review, exceptional schemes, global term minima and lower-power GVC remain open. Qualitative cubic coverage and both completed lower witnesses stay closed internally.
