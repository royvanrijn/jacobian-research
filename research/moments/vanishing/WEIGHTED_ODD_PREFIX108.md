# Weighted centered odd sources: sufficient prefix108

15 November 2026. Internal written refinement of the [endpoint/odd126 proof](CONIC108_ODD126_GLOBAL130.md). The original Wick functional, complete quotient bases, centered coefficient/frame cover and analytic exclusions are retained.

Every complex odd ternary polynomial P of total degree at most three satisfies

    E(P^j)=0 for j=1,...,108 iff all positive pure moments vanish.

The least odd prefix is therefore in **[20,108]**. The unrestricted global interval stays **[20,130]**, and the homogeneous interval stays **[14,108]**. Smooth-leading **[20,94]** and homogeneous-smooth **[14,56]** remain distinct. The first-nonzero20/14 witnesses are complete. Even108 gives no107 parity deduction. Exact leastness and20/14 sufficiency remain open.

The [checker](verify_weighted_odd_prefix108.py), [certificate](weighted_odd_prefix108_certificate.json), [receipt](weighted_odd_prefix108_verification.json) and [source audit](weighted_odd_prefix108_source_audit.json) support exact auxiliary algebra. All-degree source/reduction inductions, complete centered coefficient/frame coverage and analytic nonzero zeroth coefficients remain written inputs. There is no independent replay, external review, formal verification or expanded global radical certificate.

## 1. Centered odd phases and their full Euler sources

Retain the preceding exact odd normalizations. P2=P0=0, the original Gaussian remains centered under the linear leading-form normalization, and arbitrary Gaussian cross coefficients ju,jv are kept. D=m is a nonzero constant on each nonisotropic complement. The D=0 and required-axis-zero cases remain on their already covered frame faces.

For transverse and triangle types the ordinary phase is

    W=D*(a*(u+jv)^3+e*t^2*(u+jv)+f*t^2*(v+ju)
                                      -q2(u+jv,v+ju))
                                  +(u*v-ju*jv+ell*t^2)^2/2,

with a!=0 for transverse and a=0 for triangle. On the retained complement the leading quadratic axis coefficients needed by the projection are nonzero constants. For the tangent type,

    W=D*(u*v^2+e*t^2*u+f*t^2*v-q2(u,v))
                                  +(u^2-ju*u-jv*v+ell*t^2)^2/2.

Here q2 is a general constant-coefficient quadratic form. Each identity follows from the same retained Gaussian completion, shift and rescaling; the tangent completion identity is expanded again by the checker. No Gaussian cross coefficient is set to zero.

Keep delta=t^2*D and h=t^3*D^2. Ordinary types use the retained radial Euler connection; tangent uses the retained unsheared weighted Euler connection. Their complete operators are

    nabla(f)=f_t+(V(f)+Sbar*f)/h.

For ordinary types the Euler coefficients are u/2 and v/2 times D after common-time normalization. For tangent they are u/2 and3v/4 times D. The complete source includes the ordinary binary volume divergence. The30 exact gauge identities check it explicitly against the original normalized source, including its scalar volume term. No logarithmic-measure divergence convention is substituted here.

## 2. Weights adapted to each restricted phase

Give t weight1 and use the following amplitude weights and source bounds:

| Odd complement | weight(u),weight(v) | Source/seed cost c | Normal projection cost |
|---|---|---|---:|
| Transverse line-smooth-conic | 4/3,2 | 16/3 | 0 |
| Three independent lines | 2,2 | 6 | 0 |
| Tangent line-smooth-conic | 1,3/2 | 4 | 0 |

The constant quartic source term is removed by the full Euler representative. For transverse the remaining coefficients have the structure

    a*u^3+B*u^2+F*v^2+J(t)*u*v+L(t)*u+N(t)*v+K(t),

where a,B,F are constant in t, deg J,L,N<=2 and deg K<=4. The largest source weight is2+4/3+2=16/3, from the possible t^2*u*v term. For triangle the cubic axis term is absent; the same calculation gives6. In the odd tangent phase, top u^4/2+D*u*v^2 has weight4, while its lower derivatives obey weight(R_u)<=3 and weight(R_v)<=2. Its complete weighted Euler source has weight<=4.

In all three cases h*partial_t, h_prime and the complete Euler vertical field add at most2 to total weight. Consequently the full actual seed recurrence

    R_(j+1)=h*(R_j)_t-j*h_prime*R_j+V(R_j)+Sbar*R_j

has weighted degree<=c*j. The denominator derivative and every amplitude derivative are retained. These weights apply to the centered odd coefficient subfamilies, not arbitrary lower polynomial terms.

## 3. Every principal normal reduction preserves the new weight

All required leading divisors are constant in time. Use the unchanged full bases and principal branches, dividing by those constants in the coefficient field. The transverse branch estimates are seen at their smallest inputs: u^3 has weight4 and its normal v^2 and possible t^2*v terms have weight<=4; v^3 has weight6 and its possible t^2*v^2 term has weight6. Mixed branches are lower. Multiplication by additional amplitude powers adds the same weight to both sides, and derivative-of-amplitude terms have no larger weight because delta has weight2.

For triangle, u^2 has weight4; the possible t^2*u and t^2*v terms also have weight4, as does its retained v^2 row. The v^3 branch has weight6, as do the possible t^2*u*v and t^2*v^2 terms. The mixed and special u^2*v branches obey the same bound. These are the explicit retained principal formulas, including their divergence corrections, not a commutative critical projection.

For tangent, the u^3 to v^2 normal relation preserves weight3. A mixed input has its v-derivative witness of weight input-5/2, followed by R_v of weight<=2. The pure-u branch has u-witness weight input-3 followed by R_u<=3. In the pure-v branch the u-witness has weight input-3 followed by R_u<=3, while the v-witness has weight input-5/2 followed by R_v<=2. Their delta-divergence terms are smaller. Thus every branch preserves total weight for arbitrary exponents.

The original ordinary-degree or(2,3)-weight termination arguments remain in force. Constant padding has zero weight, so complete normal rows satisfy

    deg_t N_(e,j)<=floor(c*j-weight(e)).

The retained actual seed depths3j-2,3j-1,7j-4 for transverse,triangle,tangent are sufficient. There is no new time-dependent basis and no missing basis derivative. Exact witnesses are retained in the original twisted quotient, and H_j=h^j has valuation3j.

## 4. Zeroth cofactor bounds and every rank

At seed rank s select any rows including the unit, and use

    deg a0<=floor(c*s*(s+1)/2-sum_selected weight(e)).

The checker includes every smaller actual seed rank and row subset. At full rank the bounds are

    transverse: floor((16/3)*21-40/3)=98,
    triangle:   6*15-12=78,
    tangent:    floor(4*15-15/2)=52.

The same analytic exclusions supply nonzero zeroth coefficients on the relevant complements, since the operator and original Wick seed are unchanged. Higher scalar coefficients contain H_j with valuation>=3j, so the initial-order contradiction loses no additional derivative order. The three even bounds are sufficient; no numerical optimality is claimed.

## 5. Add back every odd frame face

Odd transverse and triangle classes combine their new complements98/78 with conic108, giving full108. Odd tangent combines complement52,parent92 and rank-one95; oddness makes moment95 automatic, so94 suffices on that rank-one face and the full tangent class. Odd node108 and odd cusp94 remain from the preceding proof.

| Leading type | Full odd sufficient gate/exclusion |
|---|---:|
| Degree at most two | 7 |
| Smooth cubic | 94 |
| Irreducible node | 108 |
| Irreducible cusp | 94 |
| Transverse line-smooth-conic | 108 |
| Tangent line-smooth-conic | 94 |
| Three independent lines | 108 |
| Three distinct concurrent lines | 104 |
| Repeated leading line(s) | 98 |

The maximum108 proves the odd-only test. It is a complete restricted coefficient/class cover, not a parity deduction from the even unrestricted130 bound. Homogeneous108 and full W-linear108 remain unchanged. No extra Gaussian-frame condition is imposed on the retained unrestricted global130 test.

The checker verifies three universal weighted source strata,30 exact source-gauge identities,262 principal projection/divergence identities,18 actual seed identities in six full nondegenerate Gaussian models, and all cofactor/rank/row arithmetic. Each model retains Gaussian cross terms and a nonzero metric determinant. The accepted run completes in about6 seconds. Universal all-degree and full coefficient/frame coverage and the analytic exclusions remain written, independently unreviewed inputs.

All previous proof bodies and receipts stay unchanged. No fresh literature-coverage claim, external theorem or producer package is imported; the preceding selected primary-source audit remains in force. Exact least global[20,130],odd[20,108],homogeneous[14,108],20/14 sufficiency, expanded radical cofactors, independent analytic review, exceptional schemes, global term minima and lower-power GVC remain open. Qualitative cubic coverage and both lower witnesses stay closed internally.
