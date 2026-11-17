# Endpoint degree strata and odd prefixes: conic108, odd126, homogeneous108

15 November 2026. Internal written refinement of the [rank-factor/global130 proof](RANK_FACTORS_AND_REDUCTION_GLOBAL_130.md). The original Wick functional, full quotient bases, complete coefficient/frame covers and analytic exclusions are retained.

| Scope | Sufficient pure prefix | Least-prefix interval |
|---|---:|---:|
| Every complex ternary polynomial of total degree<=3 | 130, unchanged | [20,130] |
| Odd ternary cubics | **126** | **[20,126]** |
| Homogeneous ternary cubics | **108** | **[14,108]** |
| Complete sixteen-coefficient W-linear family | **108** | Not computed |
| Cubic conic family | **108**, an exclusion | Not computed |

No extra Gaussian-frame condition is imposed on the global130 test. Pure nullity gives the retained W-linear frame and sharp mixed cutoff m>3*deg Q. The smooth-leading[20,94] and homogeneous-smooth[14,56] intervals and completed first-nonzero20/14 witnesses remain unchanged. All displayed upper bounds are even; no one-step parity reduction follows. Leastness and20/14 sufficiency remain open.

The [checker](verify_endpoint_degrees_and_odd_prefixes.py), [certificate](endpoint_degrees_and_odd_prefixes_certificate.json), [receipt](endpoint_degrees_and_odd_prefixes_verification.json) and [source audit](endpoint_degrees_and_odd_prefixes_source_audit.json) support auxiliary exact algebra. Universal degree/module/path inductions, full coefficient/rank coverage and analytic nonzero zeroth coefficients remain written inputs. Independent replay, external review, formal verification and expanded global radical cofactors remain outstanding.

## 1. Positive order1, negative order3

Retain the structured phase

    W=A0+sum_(k=1)^3(A_k*y^k+K^k*B_k*y^(-k)),
    deg A_k,deg B_k<=3-k, deg K<=2.

On p1,q3, A2=A3=0 and A1 is a nonzero polynomial of actual degree ell0,1 or2. The same negative Euler representative and logarithmic measure dy/y give h=t^3*K, vertical field t^2*(t*K_prime-2K/3)*y*partial_y and source support[-2,1]. Retain every denominator derivative. Unreduced seed coefficient bounds are

    deg K0: deg R_(n,j)<=3j-|n|,
    deg K1: deg R_(n,j)<=4j-max(n,0),
    deg K2: deg R_(n,j)<=5j-n.

The first two were already proved in the preceding two-sided audit; the third is the retained quadratic-K joint grade. Negative coefficients contain K^(-n). The full basis is y^-3,y^-2,y^-1,1. At most j high cancellations are needed, and low division removes K^3 already in the negative module. Thus only A1^j is required as a common normal denominator. High paths never create exponents below-3, and low paths never create exponents above0. This remains true when A1,K,t share roots: high outputs are kept rational until A1^j is cleared, and no unsupported polynomial-module assertion is applied to them.

Until endpoint clearing, degree denotes rational degree at infinity. Let n>=1 start a high path and k in[-3,0] be its final row. The largest coefficient-growth per unit downward displacement is as follows; phase and Euler terms are both included.

| deg K | ell0 | ell1 | ell2 |
|---|---:|---:|---:|
| 0 | 2 | 1 | 0 |
| 1 | 2 | 1 | 1/2 |
| 2 | 2 | 3/2 | 1 |

For example the Euler step has degree increase2-ell and displacement1; the negative phase steps have displacements2,3,4 and coefficient degree increases deg(W_-r)-ell. Taking the maximum of these ratios proves the table for every path. Combine it with n<=j and then add ell*j for A1^j padding. Low paths have coefficient degree increase at most upward displacement for constant K, and never increase coefficient degree for linear or quadratic K after their K^3 division. The resulting safe completed row bounds are

| deg K | ell0 | ell1 | ell2 |
|---|---|---|---|
| 0 | 4j-2k | 4j-k | 5j |
| 1 | 5j-2k | 5j-k | 6j+ceil(-k/2) |
| 2 | 7j-2k | 8j+ceil(-3k/2) | 9j-k |

They also bound initially retained basis terms and low outputs, not just high paths. The new H_j=h^j*A1^j has degree at most(3+deg K+ell)j and valuation>=3j. All actual seed ranks and row subsets containing1 give gates

    deg K0: 52/46/50,
    deg K1: 62/56/64,
    deg K2: 82/90/96,

in ell0/1/2 order. Hence the entire p1,q3 locus has sufficient exclusion96, replacing130. The p0,q3 regular-end case retains69 and its unit-containing full basis. Identically zero K belongs to a lower actual negative order, already covered there.

## 2. The linear-K p3,q2 locus admits the same one-sided seed grade

Retain the positive logarithmic Euler representative from the preceding proof, with h=t^3*K, source K*(t*W_t-2W-t^2+(2/3)*y*W_y), vertical field(2/3)*t^2*K*y*partial_y and support[-2,2]. For linear K its source satisfies deg S_k<=4-max(k,0), so the whole recurrence gives deg R_(n,j)<=4j-max(n,0). The negative endpoint is K^2*B2 and the only normal clearing is B2^(2j-2), with no coprimality open.

If B2 is linear, a low step divides a degree3 endpoint while every other coefficient has degree<=3, so low coefficient degree never increases. High paths have degree<=4j-k before padding. Clearing B2^(2j-2) gives the uniform row bound6j-2+max(-k,0). Its full rank5 cofactor degree is83.

If B2 is constant, low coefficient growth is at most upward exponent displacement. Its outputs have degree<=6j+k, while high outputs remain<=4j-k. The safe bound6j+|k| on rows-2,...,2 gives96. Every smaller seed rank and row choice is included. The former129 estimate is superseded by these two whole-degree strata, not by a generic sample.

The complete conic cover is now

    p3,q3, deg K0/1/2: 72/90/108;
    p2,q3, deg K0/1/2: <=85/100/105;
    p1,q3: <=96; p0,q3: <=69;
    p3,q2, deg K0/1/2: <=87/96/95;
    every other q<3 case: <=108.

Thus **conic108** is uniform for cubic G, not arbitrary-degree G. Parent92 and rank-one95 give **full W-linear108**. Homogeneous smooth56,cusp10,reducible10,nonisotropic-node72 and isotropic-node conic108 give **homogeneous108**. Even108 gives no107 parity deduction.

## 3. Odd node and cusp complements have constant D

For an odd polynomial P, P2=P0=0. The leading-form normalizations are linear changes, so the transformed Gaussian remains centered, with arbitrary nondegenerate covariance. In the nonisotropic-node/cusp complement alpha=0 and D=m!=0 is constant. Gaussian cross coefficients ju,jv need not vanish; they are not affine mean terms. Gaussian linear terms and j0 vanish, while the linear polynomial coefficient ell multiplying T remains arbitrary.

With q2 a general constant-coefficient quadratic form, the exact rescaled phases are

    nodal:
      W=D*(a*(u+jv)^3+d*(v+ju)^3
             +e*t^2*(u+jv)+f*t^2*(v+ju)-q2(u+jv,v+ju))
                                  +(u*v-ju*jv+ell*t^2)^2/2,
    cusp:
      W=D*(v^3+e*t^2*u+f*t^2*v-q2(u,v))
                                  +(u^2-ju*u-jv*v+ell*t^2)^2/2.

The node identity is checked as an exact specialization of the universal retained phase-shift checkpoint; the cusp identity is checked from the retained completion and rescaling. Both have joint degree<=4. Their full common-time Euler sources have joint degree<=4 as well. The same actual numerator recurrence therefore gives deg_joint R_j<=4j, including all denominator and amplitude derivative terms.

For node, the polynomial-normal module already proved has joint cost0; constant D makes its normal denominators constant in time. For cusp the lower phase R satisfies joint_degree R_u<=3 and joint_degree R_v<=2. Its u^3 and v^2 principal divisions therefore also preserve joint degree: their witnesses have amplitude degrees input-3 and input-2, respectively, while delta has degree2. This proves cost0 for every principal branch. No actual seed or quotient is replaced by a commutative critical algebra.

At the same retained full bases and depths, the normal row bounds are4j-|e|, with H_j=t^(3j)*D^(2j). The cofactor bounds are

    odd nonisotropic node:4*28-11=101,
    odd nonisotropic cusp:4*21-9=75.

Each bound is odd, so the original centered Gaussian moment at that order vanishes automatically for an odd P. The sufficient complement prefixes are therefore100 and74. All smaller actual seed ranks and row subsets containing the unit are included. D=0 is assigned to the retained frame branches; it is not used as a constant nonzero complement.

## 4. Complete odd126 cover

The odd nodal full class combines complement100 and conic frame108, giving108. For odd cusp combine complement74, parent92, and rank-one95; oddness makes moment95 automatic, so94 suffices on the rank-one face and on the full odd cusp. Generic transverse complement127 similarly gives126 on odd P, while its conic frame has108. Repeated-line99 becomes98 on odd P. The complete odd cover is

| Leading type | Odd sufficient gate/exclusion |
|---|---:|
| Degree at most two | 7 |
| Smooth cubic | 94 |
| Irreducible node | 108 |
| Irreducible cusp | 94 |
| Transverse line-smooth-conic | 126 |
| Tangent line-smooth-conic | 120 |
| Three independent lines | 118 |
| Three distinct concurrent lines | 104 |
| Repeated leading line(s) | 98 |

Its maximum126 proves the odd-only sufficient prefix. This is a complete restricted coefficient/class cover, not a parity inference from the even global130 bound.

The unrestricted cover is now quadratic7,smooth94,nodal129,cusp130,transverse127,tangent120,triangle118,concurrent104,repeated99. Its maximum remains130, controlled solely by the cusp complement. In the rational twenty-coefficient ring radical(I130)=radical(I_infinity) remains valid, without ordinary ideal equality or expanded global cofactors. At most7 exact monotone radical-search rounds against I130 would determine the least global K in principle; none has run.

## Evidence and boundary

The checker verifies four generic conic source strata, the exact two centered odd phase/source identities,78 actual seed reductions in22 conic and four odd binary models, and all endpoint-degree/rank/row arithmetic. Four conic models are realized Gaussian coordinate phases; enlarged synthetic cases include constant/linear/quadratic endpoints and shared endpoint/K/t factors. All odd binary metrics have nonzero determinant and retain Gaussian cross terms. The accepted run takes about14 seconds. The initial input retained an obsolete assertion restricting endpoint degree to0/1; it stopped when the new quadratic endpoint was encountered, after universal source checks. That failed input and terminal receipt remain preserved; the accepted assertion explicitly allows degree2 on p1.

Universal all-degree degree/module/path and centered Gaussian coefficient/frame coverage, period comparison and analytic nonzero zeroth arguments remain written inputs, independently unreviewed. No fresh literature-coverage claim or new external theorem/producer package is imported; the previous selected primary-source audit stays in force. All earlier proof bodies and receipts remain unchanged. Exact leastness,20/14 sufficiency, expanded radical cofactors, independent analytic review, exceptional schemes, global term minima and lower-power GVC remain open. Qualitative cubic coverage and both lower witnesses stay closed internally.
