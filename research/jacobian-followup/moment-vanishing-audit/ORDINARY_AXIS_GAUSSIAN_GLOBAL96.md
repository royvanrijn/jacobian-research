# Exact quadratic axis projection: general and odd sufficient prefix96

15 November 2026. Internal written refinement of the [general98 cover](CONCURRENT_SEEDS_REPEATED_VALUATIONS_GLOBAL98.md), retaining the original Gaussian Wick functional, complete coefficient/frame choices and analytic nonzero-zeroth exclusions.

Every complex ternary polynomial of total degree at most three has all positive pure Gaussian moments zero if and only if its first **96** pure moments vanish. The least general and odd prefixes each lie in **[20,96]**. The homogeneous interval remains **[14,96]**; smooth-leading **[20,94]**, homogeneous-smooth **[14,56]**, conic/full-W-linear96 and both complete first-nonzero20/14 witnesses remain unchanged. Even96 gives no95 parity deduction. Equality of these upper bounds does not establish equality of least thresholds.

The new uniform complementary bounds are transverse77 and triangle82. Their retained conic/frame gate96 gives full96 for both types. The full cover is quadratic7,smooth94,node96,cusp95,transverse96,tangent95,triangle96,concurrent92,repeated96. Its maximum96 now also implies odd96. Numerical leastness remains open.

The [checker](verify_ordinary_axis_gaussian_projection.py), [certificate](ordinary_axis_gaussian_projection_certificate.json), [receipt](ordinary_axis_gaussian_projection_verification.json) and [source audit](ordinary_axis_gaussian_projection_source_audit.json) bind exact auxiliary algebra. The universal formal Laurent quotient/regular axis connection,restored-direction/complementary-minor argument,joint-degree specialization,polynomial extension and complete analytic/frame cover remain independently unreviewed written inputs. No independent replay,external review,formal verification or expanded radical certificate is asserted.

## 1. Complete the quadratic variable exactly

On the retained transverse and triangle complements the actual full phase is

    W=D*(a*u^3+Le*u^2+Lf*v^2+J*u*v+L0*u+N0*v+K)
                                                    +(u*v+r)^2/2,
    delta=t^2D, h=t^3D^2.

Here a!=0 for transverse and a=0 for triangle. The original Gaussian means and all lower/cross terms remain included. Transverse requires Lf not identically zero; triangle also requires Le not identically zero. The excluded identically zero axes and D remain on their separately retained frame charts.

As a polynomial in v write

    W=A*v^2/2+B*v+C,
    A=u^2+2D Lf,
    B=(r+D J)*u+D N0,
    C=D*(a*u^3+Le*u^2+L0*u+K)+r^2/2.

For the following formulas one may absorb D J into a new symbol R=r+D J in B and use C as displayed. At Lf=0 the residual phase contains an additional regular polynomial/constant contribution from J; it does not change its endpoints or either comparison determinant below. For explicit endpoint formulas below, use R in B and C, so C's constant r^2/2 is retained rather than silently replaced by R^2/2.

Exact quadratic completion gives Phi=C-B^2/(2A). Conditional formal Gaussian moments in v have mean -B/A and variance -delta/A. For every polynomial amplitude this is a finite algebraic operation; no asymptotic integration or convergence is assumed. The full period normalization changes from t^(-2)D^(-1/2) to t^(-1)A^(-1/2). In particular the determinant derivative is retained.

For arbitrary coefficients A,B,C and every amplitude f the normalized Gaussian operation Jv obeys

    Jv((A*v+B+delta*partial_v)f)=0,
    Jv((W_u+delta*partial_u)f)
        =delta*partial_u Jv(f)+(Phi_u-delta*A_u/(2A))*Jv(f).

The corresponding time identity also differentiates delta:

    Jv((W_prime-delta_prime*W/delta)*f+delta*f_prime)
      =delta*(Jv(f))_prime
       +(Phi_prime-delta_prime*Phi/delta+delta_prime/2
                                           -delta*A_prime/(2A))*Jv(f).

Combining delta_prime/(2delta)=1/t+D_prime/(2D) with the original prefactor gives the stated new normalization and its complete time derivative. These identities keep the original twisted operator and Gaussian functional.

## 2. A finite Laurent quotient regular along the axis

Work first over the generic coefficient field with independent formal z=Lf and with t,D,N0 and required Le,a invertible. Since u is a Laurent variable,

    1/A=u^(-2)*(1+2D z/u^2)^(-1)

is coefficientwise analytic in z. Factor A^(-1/2)=u^(-1)*(1+2D z/u^2)^(-1/2). Absorb the latter analytic unit into the amplitude and retain volume du/u.

At z=0 put R=r+D J. Then

    Phi0=D*(a*u^3+Le*u^2+L0*u+K)
        +(r^2-R^2)/2-D*R*N0/u-D^2*N0^2/(2u^2).

Its Laurent endpoints are (3,-2) for transverse and (2,-2) for triangle. The full one-variable twisted operator with volume du/u is delta*u*partial_u+u*Phi0_u. Its derivative term has exponent shift0. Both endpoint multiplication coefficients are units on this generic chart. Thus the exact twisted quotient has bases u^(-2),u^(-1),1,u,u^2 for transverse and u^(-2),u^(-1),1,u for triangle,of ranks5 and4.

For completeness, reduce an exponent k>=p using the upper endpoint u^p: the lowest new exponent k-p-2 is still>=-2. Reduce k<-2 using the lower endpoint: the highest new exponent k+2+p is<=p-1. Both procedures terminate. If a nonzero Laurent amplitude had twisted image entirely in[-2,p-1],its highest exponent would force its input maximum<=-1 and its lowest exponent would force its input minimum>=0,an impossibility. Hence the displayed representatives are independent as well as spanning. This is the full twisted quotient,including divergence,not a commutative critical-ideal sample.

Correct this division coefficient by coefficient in z. Each coefficient is a finite Laurent polynomial and endpoint units remain invertible,so the resulting complete formal quotient is free of the same rank and introduces no negative z powers. Differentiate all Gaussian,phase,volume,basis and coefficient terms to obtain its induced time connection. Its coefficients remain analytic in z,including z_prime and D_prime. Differentiating formal analytic series does not create negative z powers. No extra assertion of constant period or generic seed cyclicity is made.

Denote this projection by pi. The complete old numerator operator intertwines with h times that regular finite connection. In particular every pi(R_j) is analytic in z. The old Euler fields are retained exact integration-by-parts terms. Only this analyticity,not a stronger valuation claim about original R_j,is used below.

## 3. Restore the missing direction and bound its determinant

Retain the actual old polynomial columns and denominators. On transverse,

    E=(1,u,v,u^2,uv,v^2),
    N_j=Lf^(2j-2)*n_j,
    H_j=t^(3j)D^(2j)Lf^(2j-2).

On triangle let ai=ceil(i+k/2),bi=ceil(i/2+k) on the old basis1,u,v,uv,v^2. Then

    E_(i,k)=Le^ai Lf^bi u^i v^k,
    N_j=Le^(2j)Lf^(2j)*n_j,
    H_j=t^(3j)D^(2j)Le^(2j)Lf^(2j).

These are exactly the [retained polynomial normal modules](NORMAL_MODULE_GLOBAL_146.md) and [stopped transverse clearing](RANK_FACTORS_AND_REDUCTION_GLOBAL_130.md). Every rational divergence witness remains in the same original quotient. Consequently the finite projected columns have axis orders at least2j-2 and2j respectively.

Add one unrestricted analytic row selecting the v^2 coordinate. This restores a square comparison matrix Q. At z=0 the Gaussian mean is

    mu=-R/u-D*N0/u^2.

For transverse,removing the v^2 column leaves finite columns1,u,mu,u^2,u*mu. They are independent in the five Laurent basis positions; their determinant is a nonzero sign times D^2 N0^2. Thus Q has determinant order0.

For triangle the four finite columns are1,Le*z*u,Le*z*mu,Le^2*z^2*u*mu. Their determinant is a nonzero sign times Le^4*z^4*D^2 N0^2. The order is4. Contributions from higher z coefficients cannot lower these orders or cancel the indicated nonzero leading terms.

Complementary-minor identities express each selected minor of Q^(-1) as an analytic complementary minor of Q divided by det Q. Cauchy--Binet therefore loses at most0 or4 axis orders once per selected minor. The added coordinate remains unrestricted; the missing direction has not been discarded.

At actual seed rank s put w0=0,and wj=2j-2 for transverse or2j for triangle. A raw cofactor C_l omits column l from N0,...,Ns. Its axis order is at least the sum of the smallest s-1 weights among the remaining indices,minus the comparison loss,clipped at0. Write this bound b_l. At full ranks6/5 the zeroth bounds are

    transverse b0=0+2+4+6+8=20,
    triangle b0=2+4+6+8-4=16.

For every smaller rank and l>=1,the axis exponent of H_l plus b_l is>=b0. The complete certificate includes these inequalities and all row selections containing the unit.

## 4. Combine independent divisor powers before specializing degrees

The [previous full Gaussian projection proof](FORMAL_GAUSSIAN_PROJECTION_GLOBAL109.md) already gives raw cofactor D orders2T_s for C0 and max(0,2T_s-2l) for C_l,where T_s=s(s-1)/2. The bounded cleared columns are polynomial in independent D,Lf and the original parameters. D and Lf are distinct prime variables in that ring; their divisibilities therefore combine multiplicatively.

Divide the whole scalar relation by

    D^(2T_s)*Lf^b0.

Higher scalar coefficients H_l*C_l remain divisible by t^(3l),since their separate D and Lf exponents meet both common powers. All other axis factors have nonnegative exponents. This remains true when physical D,Lf,t share a root. No r,N0,J,Le or trace factor is cancelled. The retained analytic nonzero zeroth coefficient stays nonzero.

To use this cancellation uniformly for constant divisors,establish the degree bound in the independent weighted polynomial ring before specialization. Give t,u,v,D,Le,Lf,J weight1,L0,N0,r weight2 and K weight3,with their time derivatives assigned one lower weight and D_prime weight0. W has joint degree4,the complete ordinary source has joint degree5 and amplitude<=3,and actual R_j has joint degree<=5j. The existing polynomial-normal reductions preserve this grading with their bounded axis clearings. They give exactly

    transverse: row joint degree<=7j-2-|e|,
    triangle: row joint degree<=9j-|e|-ai-bi.

Thus raw full-rank zeroth joint degrees are127 and118. A homogeneous monomial factor D^a Lf^b of weight a+b lowers polynomial joint degree by a+b upon division. This proves

    transverse:127-30-20=77,
    triangle:118-20-16=82.

All smaller actual ranks and selected row sets have no larger bounds. Specializing D,Lf,Le to affine or nonzero constant time polynomials and all other coefficients to their retained degree-bounded time polynomials cannot increase weighted joint degree. Hence77/82 hold uniformly,including constant D or constant axes. This is not subtraction of the degree of a constant factor from a bound known only after univariate specialization: polynomial division and weighted degree are established first in the independent ring.

Generic inversions of N0,D,t,Le and leading endpoint coefficients were used only to prove polynomial divisibilities. Their low divisor coefficients vanish over the generic field and therefore identically. Polynomial extension includes N0=0,zero axis slopes,shared roots and degenerate leading coefficients. Identically zero required axes and D stay on their separately retained frame charts. No generic chart is imposed on the final full class statement.

## 5. Complete general/odd96 cover and assurance

The uniform complements77/82 lie below the retained conic/frame96. Full transverse and triangle gates are therefore96. Together with the preceding cover this gives

    quadratic7,smooth94,node96,cusp95,transverse96,tangent95,
    triangle96,concurrent92,repeated96.

Its maximum96 is the general pure-only gate. It also gives odd96 by inclusion of the odd coefficient subfamily. The homogeneous gate was already96. No extra Gaussian-frame hypothesis is added. All-order nullity retains the internally closed ternary cubic GMC classification and sharp mixed cutoff m>3deg Q; odd mixed cutoff stays m>deg Q.

In the rational twenty-coefficient ring radical(I96)=radical(I_infinity),without ordinary ideal equality or expanded cofactors. At most7 exact monotone radical-search rounds against I96 determine the least general prefix in principle; none has been run in the full coefficient ring. Equality of consecutive moment ideals is insufficient.

The checker verifies21 exact quadratic Gaussian connection identities,two comparison determinants,31 finite polynomial-normal joint-degree identities and10 complete scalar cofactor families with65 coefficient checks. Six families are retained full Gaussian models;four fresh full Gaussian constant-axis models additionally verify22 original twisted seed identities in the [preserved producer input](constant_axis_gaussian_probe_input.py),[model data](constant_axis_gaussian_models.json) and [terminal receipt](constant_axis_gaussian_probe_receipt.json). The receipt states its original tmp replay location. Models include D(0)=0,shared D/Lf roots,constant D and constant axes. No previous checker main is rerun;no failed new input occurred.

Universal coefficientwise Laurent division and regular axis connection,the restored-direction/complementary-minor argument,independent-variable divisibility,joint-degree specialization,polynomial extension and complete analytic/frame cover remain independently unreviewed written inputs. Finite identities support auxiliary algebra and do not certify the complete analytic theorem. All earlier proof bodies and receipts remain byte-identical. No fresh literature-coverage claim,external theorem or producer package is imported.

Least general[20,96],odd[20,96],homogeneous[14,96],20/14 sufficiency,expanded I96 cofactors,independent analytic review,exceptional schemes,global term minima and lower-power GVC remain open. Conic/full-W-linear96,node/repeated96 and rank-one95/smooth94 now control the common numerical frontier. Qualitative cubic coverage and both lower witnesses remain closed internally.
