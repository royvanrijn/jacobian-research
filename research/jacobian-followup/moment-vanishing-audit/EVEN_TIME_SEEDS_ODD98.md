# Even-time seed coefficients: complete odd prefix98

15 November 2026. Internal written refinement of the [odd concurrent36/odd100 cover](ODD_CONCURRENT36_PREFIX100.md), retaining its full centered coefficient/frame cover and original Wick functional.

Every complex odd ternary polynomial of total degree at most three has all positive pure Gaussian moments zero if and only if its first **98** pure moments vanish. The least odd prefix is therefore in **[20,98]**. The general interval remains **[20,130]**, homogeneous **[14,96]**, smooth-leading **[20,94]** and homogeneous-smooth **[14,56]**. Conic/full-W-linear96 and both exact first-nonzero20/14 witnesses remain unchanged. Even98 gives no97 parity deduction.

The improvement uses parity of the *completed polynomial normal coefficients*, not only automatic vanishing of odd moments. It lowers odd node100 to98, cusp complement74 to72 and tangent complement52 to50. Transverse98 and triangle78 complements remain unchanged. No corresponding arbitrary-lower-term improvement is asserted.

The [checker](verify_even_time_odd_seeds.py), [certificate](even_time_odd_seeds_certificate.json), [receipt](even_time_odd_seeds_verification.json) and [source audit](even_time_odd_seeds_source_audit.json) bind auxiliary exact algebra. Universal all-degree parity/reduction induction, complete coefficient/frame choices, original period comparison and analytic nonzero zeroth coefficients remain written, independently unreviewed inputs. There is no independent replay, external review, formal verification or expanded global radical certificate.

## 1. The five centered odd phases lie in C[t^2,u,v]

Retain the exact centered phases and original ordinary or weighted Euler representatives from the [endpoint proof](CONIC108_ODD126_GLOBAL130.md) and [weighted source proof](WEIGHTED_ODD_PREFIX108.md). In each nonisotropic complement D=m is a nonzero constant. Gaussian cross coefficients remain arbitrary; leading-form normalizations are linear, so the original zero mean is retained. The nodal phase is

    W=D*(a*(u+jv)^3+d*(v+ju)^3+e*t^2*(u+jv)
              +f*t^2*(v+ju)-q2(u+jv,v+ju))
                      +(u*v-ju*jv+ell*t^2)^2/2.

Set d=0 for transverse and a=d=0 for triangle. Cusp and tangent phases are respectively

    W=D*(v^3+e*t^2*u+f*t^2*v-q2(u,v))
                      +(u^2-ju*u-jv*v+ell*t^2)^2/2,
    W=D*(u*v^2+e*t^2*u+f*t^2*v-q2(u,v))
                      +(u^2-ju*u-jv*v+ell*t^2)^2/2.

Every coefficient is a polynomial in t^2. Their complete normalized sources Sbar are likewise even in t, including the full ordinary-volume divergence. The common scalar denominator is h=t^3*D^2, which is odd in t; h_prime is even. The full Euler vertical fields have the form

    V=D^2*t^2*(bu*u*partial_u+bv*v*partial_v),
    (bu,bv)=(1/2,1/2) on node/transverse/triangle,
             (1/2,2/3) on cusp, (1/2,3/4) on tangent.

No shear representative is substituted on the tangent stratum. These are the same normalized operators used by the retained actual seed certificates.

## 2. An exact recurrence over tau=t^2

Write R_j(t,u,v)=F_j(tau,u,v), tau=t^2, and retain the full actual recurrence

    R_(j+1)=h*(R_j)_t-j*h_prime*R_j+V(R_j)+Sbar*R_j.

The chain rule gives exactly

    F_(j+1)=2D^2*tau^2*(F_j)_tau-3j*D^2*tau*F_j
         +D^2*tau*(bu*u*(F_j)_u+bv*v*(F_j)_v)
         +Sbar(tau,u,v)*F_j.

Starting with F0=1, this proves every unreduced R_j is even, for all j. It keeps the derivative of h and every vertical derivative. The computation is a polynomial recurrence identity; it does not replace the original scalar differentiation by an unaccounted derivative in tau. The actual normalized seeds still have denominators H_j=t^(3j)*D^(2j), with valuation3j.

Each retained principal branch has coefficients in C(t^2): its phase coefficients and delta=t^2*D are even, all its required leading divisors are constant in time, and ordinary amplitude differentiation preserves that coefficient field. Its exact divergence witnesses and remainder therefore lie in C[t^2,u,v] after division by those constant coefficient-field units. The old termination and degree proofs still apply. The fixed quotient bases are unchanged. Induction through each completed reduction gives

    N_(e,j) in C[t^2].

There is no new time-dependent basis, scalar connection or missing basis derivative. The D=0 and required-axis-zero cases remain on their retained frame faces, rather than being treated as nonzero constant divisors.

## 3. Round each polynomial coefficient before taking the determinant

Retain the previous row bound deg_t N_(e,j)<=c*j-w(e). Evenness strengthens it to

    deg_t N_(e,j)<=2*floor((c*j-w(e))/2).

For a selected set of s basis rows containing the unit, the zeroth-cofactor determinant uses columns j=1,...,s. Its degree is at most

    max_permutations sum_(j=1)^s
                      2*floor((c*j-w(e_permutation(j)))/2).

This formula bounds every determinant term; it does not assume the leading determinant coefficient is nonzero. Include every smaller actual seed rank and every row set containing the unit. For rational source weights, the rounding need not separate into independent column and row terms, so the checker computes the maximum assignment exactly.

| Complement | Source cost c | Basis weights w(e) | Even-coefficient cofactor bound |
|---|---:|---|---:|
| Node | 4 | 0,1,1,2,2,2,3 | **98** |
| Cusp | 4 | 0,1,2,1,2,3 | **72** |
| Transverse | 16/3 | 0,4/3,2,8/3,10/3,4 | 98 |
| Triangle | 6 | 0,2,2,4,4 | 78 |
| Tangent | 4 | 0,1,2,3/2,3 | **50** |

For the node, each of the three odd-weight rows loses one additional degree:4*28-(11+3)=98. For cusp the analogous value is4*21-(9+3)=72. Tangent rounds its row penalties to0,2,2,2,4, giving4*15-10=50. The transverse fractional assignment still permits98; no further decrease is inferred from rounding the total alone. Triangle rows already have the required parity, so78 remains.

The same retained nonzero-zeroth-coefficient analytic exclusions apply to the original quotient connection and seed. All scalar relation coefficients remain polynomial, and higher ones retain H_j with valuation at least3j. Consequently each displayed even cofactor bound is a sufficient pure-prefix exclusion on its complement. This is stronger information than the earlier step that used odd-moment vanishing to turn odd node101 into100 and cusp75 into74. No assertion of numerical optimality is made.

## 4. Reassemble every odd class and frame face

Odd node combines its complement98 with conic96, giving full98. Odd cusp combines72,parent92 and rank-one95-to94, retaining94. Odd tangent combines50 with the same92/94 frame bounds, retaining94. Transverse98 and triangle78 still combine with conic96 to give full98 and96. Concurrent36 and repeated98 remain from the previous complete cover.

| Leading type | Full odd sufficient prefix |
|---|---:|
| Degree at most two | 7 |
| Smooth cubic | 94 |
| Irreducible node | **98** |
| Irreducible cusp | 94 |
| Transverse line-smooth-conic | 98 |
| Tangent line-smooth-conic | 94 |
| Three independent lines | 96 |
| Three distinct concurrent lines | 36 |
| Repeated leading line(s) | 98 |

The maximum98 proves the full odd-only test. Node, transverse and repeated-line cases now share that sufficient upper bound. Homogeneous96, conic/full-W-linear96 and unrestricted130 stay unchanged. The unrestricted cover retains generic node129,cusp130,transverse127,tangent120,triangle118,concurrent104,repeated99,quadratic7 and smooth94. There is no added frame hypothesis on the original polynomial.

## Exact evidence and remaining assurance

The checker verifies five generic even phase/source strata;120 symbolic-order tau recurrence identities;225 complete principal/divergence parity identities;30 original twisted seed identities directly over tau from ten retained nondegenerate centered Gaussian models;174 even normal rows; all scalar denominator valuations; and every rank/row determinant assignment. All original Gaussian cross coefficients are retained in the generic identities. The first run passes in about6 seconds. No prior checker main is run and no prior receipt or proof body is rewritten.

The finite checks support the all-degree branch and parity argument; they are not formal certification of the full analytic theorem. Universal reduction termination, field/axis/frame coverage, period comparison and analytic exclusions remain independently unreviewed written inputs. No fresh literature-coverage claim, new external theorem or producer package is imported; the prior selected primary-source refresh stays retained.

Least global[20,130],odd[20,98],homogeneous[14,96],20/14 sufficiency, expanded I130 radical cofactors, independent analytic review, exceptional coefficient schemes, global term minima and lower-power GVC remain open. Generic cusp130/node129, odd node/transverse/repeated98 and homogeneous conic96 are the next controlling estimates. Both lower witnesses and qualitative cubic coverage stay closed internally.
