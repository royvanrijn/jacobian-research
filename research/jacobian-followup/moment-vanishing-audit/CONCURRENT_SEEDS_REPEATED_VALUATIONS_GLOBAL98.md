# Actual concurrent seeds and repeated-line valuations: general prefix98

15 November 2026. Internal written refinement of the [general104 cover](GAUSSIAN_PROJECTION_INFINITY_GLOBAL104.md), retaining the original Gaussian Wick functional, complete leading-type/frame cover and analytic nonzero-zeroth exclusions.

Every complex ternary polynomial of total degree at most three has all positive pure Gaussian moments zero if and only if its first **98** pure moments vanish. The least general prefix is in **[20,98]**. The odd least interval is also **[20,98]**; homogeneous **[14,96]**, smooth-leading **[20,94]**, homogeneous-smooth **[14,56]**, conic/full-W-linear96 and both complete first-nonzero20/14 witnesses remain unchanged. Even98 gives no97 parity deduction.

The new full concurrent-line gate is **92**, from complements62/52 and the retained isotropic frame92. The repeated-line gate falls99 to **96**, because the actual cleared connection denominator always vanishes at zero and already compensates for differentiation in the initial-order argument. Triangle98 now controls the general cover, followed by transverse97 and node/repeated/conic96. This is sufficient, not optimal.

The [checker](verify_concurrent_seeds_repeated_valuations.py), [certificate](concurrent_seeds_repeated_valuations_certificate.json), [receipt](concurrent_seeds_repeated_valuations_verification.json) and [source audit](concurrent_seeds_repeated_valuations_source_audit.json) bind exact auxiliary algebra. Universal all-degree seed/reduction/minor induction, complete coefficient/frame choices and analytic nonzero-zeroth period arguments remain independently unreviewed written inputs. No independent replay, external review, formal verification or expanded radical certificate is asserted.

## 1. Keep the full concurrent completion and all lower terms

Normalize three distinct concurrent leading lines to u^3+v^3, keeping their vertex on the T-axis. With arbitrary affine Gaussian means and cross coefficients write

    P=alpha*T^2+T*(ell_u*u+ell_v*v+ell_0)+u^3+v^3+g2(u,v),
    Q=m*T^2/2+T*(J_u*u+J_v*v+J_0)+q2(u,v),
    D=m-2alpha*t.

The full metric is nondegenerate; no condition m!=0 is imposed. On D not identically zero, complete T exactly and rescale the old u,v to u/t,v/t. Retain the normalized period prefactor t^(-2)D^(-1/2), ordinary binary volume and phase W/(t^2D), where

    W=D*G+F^2/2,
    F=(t*ell_u-J_u)*u+(t*ell_v-J_v)*v+t^2*ell_0-t*J_0,
    G=u^3+v^3+t*g2_2+t^2*g2_1+t^3*g2_0
                                      -q2_2-t*q2_1-t^2*q2_0.

Subscripts denote homogeneous amplitude degree. Thus F is affine in u,v, with linear coefficients of time degree<=1 and constant coefficient<=2. G has quadratic coefficients of time degree<=1, linear coefficients<=2 and constant<=3. Giving t,u,v weight1, W has joint degree<=4, including linear D. Gaussian means and every lower polynomial term are retained.

## 2. A radial Euler representative cancels the binary cubic source

Put delta=t^2D,h=t^3D^2,eps=D_prime. After the common time cancellation the original source is

    S_old=(delta*W_prime-delta_prime*W)/t
                                  -2t^2D^2-t^3D*eps/2.

Add the exact ordinary twisted derivatives with witnesses(2D/3)u*f and(2D/3)v*f. The full operator becomes

    L f=h*f_prime+(2/3)t^2D^2*(u*f_u+v*f_v)+S*f,
    S=tD*W_prime-(2D+t*eps)*W+(2D/3)*(u*W_u+v*W_v)
                                  -(2/3)t^2D^2-t^3D*eps/2.

The cubic D(u^3+v^3) cancels. Hence S has amplitude degree<=2. Its joint degree is<=5 for linear D and<=4 for nonzero constant D. The scalar volume divergence in this gauge is(4/3)t^2D^2. In particular it is an ordinary-volume calculation, not a logarithmic-measure one. All determinant,denominator,time and Euler derivatives are retained.

For the actual numerators use

    R0=1, R_(j+1)=L R_j-j*h_prime*R_j.

Their amplitude degree is<=2j. Their joint degree is<=5j for linear D and<=4j for constant D. Indeed derivative and Euler terms add at most4 or2 respectively; source multiplication gives the asserted largest cost. This is an all-j degree induction on the actual numerators.

## 3. Stop the reductions by amplitude degree

The complete rank-four twisted basis is1,u,v,uv. Its leading gradients are3D*u^2,3D*v^2. Their regular sequence and strict lower-degree corrections give the retained full filtered basis, not merely a critical-point sample.

For a nonbasis monomial f of amplitude degree n, choose C_u=u^(i-2)v^k/3 when i>=2, or C_v=u^i v^(k-2)/3 otherwise. With W_low=W-D(u^3+v^3), the exact cleared remainder is

    D*f=-C_u*(W_low)_u-C_v*(W_low)_v
                   -delta*(partial_u C_u+partial_v C_v)
                   +delta*div C+W_u*C_u+W_v*C_v.

The remainder has amplitude degree<=n-1 and joint degree<=n+1. The same statement holds termwise after multiplying by time coefficients. Thus each D-cleared principal step costs at most1 joint degree and strictly decreases amplitude. A seed of amplitude<=2j needs depth at most2j-1. Pad shorter branches by the remaining powers of D.

Consequently completed polynomial rows N_(e,j)=D^(2j-1)*NF_e(R_j) satisfy

    linear D: deg_t N_(e,j)<=7j-1-|e|,
    constant D: deg_t N_(e,j)<=6j-1-|e|.

For constant D, multiplying by D contributes no degree; the direct source joint bound4 and the same at-most-one joint cost per reduction give the second estimate. This is not subtraction of a constant divisor's degree. The scalar denominator is

    H_j=t^(3j)*D^(4j-1), valuation at least3j.

At actual seed rank s choose s rows containing the unit. The zeroth cofactor degree is at most c*s*(s+1)/2-s-sum_selected|e|, with c=7 or6. Every smaller rank and row set is included. Full rank4 gives gates70-4-4=62 and60-4-4=52. Higher scalar coefficients retain H_j and thus valuation>=3j. No additional derivative-order allowance is needed.

The [complete concurrent exclusion](CONCURRENT_LINE_EXCLUSION.md) supplies the nonzero zeroth coefficient on all D!=0 members, including alpha!=0 and the nonisotropic formal-frame branches. At D(0)=0 the final relation is polynomial and acts on the original Wick series; it is not a conditional Gaussian integral at zero.

When D vanishes identically,alpha=m=0. The vertex is isotropic with zero quadratic value there, so the retained W-linear frame has an affine W coefficient, the d=0 fourteen-coefficient face, with gate92. All lower terms and means are included. Therefore the full concurrent gate is max(62,52,92)=92. The odd concurrent36 result remains separate and unchanged.

## 4. The actual repeated-line denominator vanishes at zero

Retain the [repeated-line proof](REPEATED_LINE_FINITE_PREFIX.md). In each pure-exclusion branch, elimination gives a phase U(t,u)/D(t), determinant prefactor D^(-1/2), and actual phase degree d=3 or4, with deg_t U<=3,deg_t D<=2. Its cleared connection denominator is

    H=D^2*a^d, a=d*[u^d]U,
    B=4+3d, deg H<=B, deg N_ij<=B.

The old bound allowed arbitrary degree-complete rational models with a(0)!=0 and D(0)!=0. Actual eliminated phases have an additional property: H(0)=0.

On the nonisotropic repeated-form chart, the transverse Gaussian matrix is I-2tA after a constant Gaussian normalization. With a possibly affine Gaussian cross term J(u),a quadratic polynomial vector H0(u),and deg g<=3,deg q<=2,the full numerator has the form

    U=D*(t*g-q)+(t*H0-J)^T*adj(I-2tA)*(t*H0-J)/2.

At t=0 it is -q+J^T J/2, of amplitude degree<=2. Thus the actual cubic or quartic leading coefficient a vanishes at zero. D(0)!=0 and a!=0 identically still give H(0)=0. The statement permits all cancellations that lower the actual degree from4 to3.

On the isotropic nonzero-W^2 chart, the retained elimination gives D=t*(A*t-1/2),up to a nonzero constant. It vanishes at zero, regardless of the leading coefficient. Hence H(0)=0 on this chart as well. Identically zero leading coefficients that remove both cubic and quartic phases remain on the separate retained frame branches; they are not asserted to lie in this exclusion calculation.

The conclusion H(0)=0 therefore covers every actual repeated-line pure-exclusion chart, rather than the larger unrestricted model space of the old degree checker.

## 5. Use the retained scalar valuations before adding an order allowance

At minimal seed rank s<=d-1 the retained cleared columns h_j have degree<=jB. Polynomial cofactors C_j give the same relation

    sum_(j=0)^s A_j I^(j)=0, A_j=H^j*C_j,
    deg A_j<=C=B*s*(s+1)/2.

The retained analytic exclusion ensures A0!=0. Since H has positive valuation, every higher A_j has valuation>=j. If I=1+O(t^(C+1)), then for j>=1,

    val(A_j I^(j)) >= j+(C+1-j)=C+1.

But A0 I=A0+O(t^(C+1)),and the nonzero polynomial A0 has degree<=C. This is impossible. Therefore the sufficient exclusion gate is C, without the old extra s. No factor has been divided from the scalar relation and no valuation has been lost.

At d=3 this gives13*3=39; at d=4 it gives16*6=96. All lower actual seed ranks give smaller gates. The adapted repeated-line frames retain c=0 gate35,d=0 gate92 and the restricted c*d interior quartic gate93. All are below96,so the full repeated-line class has sufficient gate96. Its qualitative classification and sharp cubic mixed cutoff stay unchanged.

## 6. Complete general98 cover

| Leading type with arbitrary lower terms | Full sufficient prefix |
|---|---:|
| Degree at most two | 7 |
| Smooth cubic | 94 |
| Irreducible node | 96 |
| Irreducible cusp | 95 |
| Transverse line-smooth-conic | 97 |
| Tangent line-smooth-conic | 95 |
| Three independent lines | 98 |
| Three distinct concurrent lines | 92 |
| Repeated leading line(s) | 96 |

The maximum98 gives the general test with no extra frame hypothesis. General and odd sufficient prefixes now coincide at98; this does not prove equality of their least prefixes. The full homogeneous gate96 is unchanged. Pure nullity retains the internally closed ternary cubic GMC conclusion and sharp mixed cutoff m>3deg Q.

In the rational twenty-coefficient ring radical(I98)=radical(I_infinity),without ordinary ideal equality or expanded cofactors. At most7 exact monotone radical-search rounds against I98 determine the least global prefix in principle; none has been run in the full coefficient ring. Equality of consecutive moment ideals is insufficient.

The checker verifies the full concurrent completion/rescaling,28 source-gauge and62 principal identities,20 original twisted seed identities and five complete scalar relations in nondegenerate Gaussian models. These include linear D,D(0)=0,nonzero constant D and zero ell. It also verifies actual repeated elimination valuations,10 exact connection columns and four full scalar relations for cubic/quartic models with regular and zero D(0). All actual-rank/row degree and higher scalar valuation checks pass. Old checker mains are not rerun. No failure occurred in the accepted new input.

Universal seed/reduction/minor induction,complete leading-type/metric/frame coverage,the seed-to-period bridge and nonzero zeroth analytic exclusions remain independently unreviewed written inputs. The finite algebra does not certify the complete analytic theorem. All earlier proof bodies and receipts stay byte-identical. No fresh literature-coverage claim,external theorem or producer package is imported.

Least global[20,98],odd[20,98],homogeneous[14,96],20/14 sufficiency,expanded I98 cofactors,independent analytic review,exceptional schemes,global term minima and lower-power GVC remain open. General triangle98,transverse97 and node/repeated/conic96 are the next numerical frontier; odd node/transverse/repeated98 remain separate. Qualitative cubic coverage and both lower witnesses stay closed internally.
