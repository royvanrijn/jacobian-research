# Actual Euler seeds lower the repeated-line pure gate to92

15 November 2026. Internal written refinement of the [repeated96 proof](CONCURRENT_SEEDS_REPEATED_VALUATIONS_GLOBAL98.md),retaining the original Gaussian Wick functional,actual elimination charts,adapted frame coverage and analytic nonzero-zeroth exclusions.

Every complex ternary cubic whose leading form is L^2 M or L^3 has all positive pure Gaussian moments zero if its first **92** vanish. The complementary cubic/quartic exclusions improve39/96 to **33/78**. The restricted mixed quartic frame also improves93 to **78**,leaving the affine-W frame92 as the full repeated-line maximum. This is sufficient,not optimal.

General,odd and homogeneous sufficient prefixes stay **96**. Least general/odd intervals remain **[20,96]**,homogeneous **[14,96]**. Conic/full-W-linear96,rank-one95,smooth-leading[20,94],homogeneous-smooth[14,56] and both complete first-nonzero20/14 witnesses stay unchanged. The [checker](verify_repeated_euler_seeds.py),[certificate](repeated_euler_seeds_certificate.json),[receipt](repeated_euler_seeds_verification.json) and [source audit](repeated_euler_seeds_source_audit.json) support auxiliary exact algebra. Universal induction,complete elimination/frame coverage,the seed-to-period bridge and analytic nonzero-zeroth arguments remain independently unreviewed written inputs.

## 1. Cancel the highest amplitude before reducing actual seeds

The retained pure-exclusion charts have period prefactor D^(-1/2) and phase U(t,u)/D(t),with actual amplitude degree d=3 or4,deg_t U<=3,deg_t D<=2. Put

    a=d*[u^d]U, L=a_prime*D-a*D_prime,
    C=D*U_prime-D_prime*U-D*D_prime/2.

Here a and D are nonzero polynomials;their roots at zero are allowed. The original connection numerator is D^2*f_prime+C*f. Multiply the entire operator by a and add the exact ordinary twisted derivative with witness -(L/d)*u*f. This gives

    h=a*D^2,
    S=a*C-(L/d)*(u*U_u+D),
    V f=-(D*L/d)*u*f_u,
    Lnew f=h*f_prime+V f+S*f.

The exact difference is D*partial_u(-(L/d)u*f)+U_u*(-L/d)u*f. Thus the ordinary-volume divergence -(L/d)D and Gaussian determinant term in C are retained. Every derivative of a,D and h remains present. The leading u^d coefficient of a*C is a*L/d;the same coefficient is subtracted by (L/d)uU_u. Consequently S has amplitude degree<=d-1.

The time degrees satisfy

    deg a<=3, deg L<=4, deg h<=7,
    deg S<=7, deg(D*L)<=6.

For actual numerators use R0=1 and

    R_(j+1)=h*(R_j)_prime-j*h_prime*R_j+V R_j+S R_j.

Time derivative and vertical terms add at most6 to coefficient degree;source multiplication adds at most7. Their amplitude degrees do not exceed the source allowance. Hence,for every j,

    deg_t R_j<=7j, deg_u R_j<=(d-1)j.

These are bounds on actual unreduced numerators. The new Euler operator has been used before seed formation;old seed columns are not reinterpreted after deleting a field or basis derivative.

## 2. Stop the original twisted reductions at the basis

Keep the complete one-variable twisted basis1,u,...,u^(d-2). One principal step on u^n,n>=d-1 is

    a*u^n=-u^(n-d+1)*(U_u-a*u^(d-1))
                  -(n-d+1)D*u^(n-d)
                  +D*partial_u(u^(n-d+1))+U_u*u^(n-d+1).

Its remainder has amplitude degree<=n-1 and coefficient time degree<=3. The all-degree estimate follows directly from these displayed terms;the finite principal checks are auxiliary evidence. Starting with amplitude at most(d-1)j,depth

    r_j=(d-1)j-(d-2)

suffices. Pad shorter branches by a. Completed polynomial normals N_(k,j) represent a^r_j R_j and have coefficient degree bounded by

    7j+3r_j=(3d+4)j-3(d-2).

Their unchanged period derivative column is obtained with scalar denominator

    H_j=h^j*a^r_j=D^(2j)*a^(d*j-d+2).

Rational quotient equality uses the full original operator D partial_u+U_u;all cleared divergence witnesses are retained and checked.

At actual seed rank s<=d-1,choose any s rows containing the unit. Its zeroth cofactor degree is at most

    (3d+4)*s*(s+1)/2-3(d-2)*s.

Every smaller rank/row choice is included. Full rank gives d3:13*3-3*2=33 and d4:16*6-6*3=78.

## 3. Retain the actual zero-valuation condition

The [preceding actual-chart argument](CONCURRENT_SEEDS_REPEATED_VALUATIONS_GLOBAL98.md) proves that either a(0)=0 or D(0)=0 on every repeated-line exclusion chart. On a nonisotropic repeated-form chart,U(0) has amplitude degree<=2,so a(0)=0. On the isotropic nonzero-W^2 chart,D=t(A*t-1/2),up to a nonzero constant. This is an actual elimination property,not a condition imposed on all arbitrary rational phase models.

The exponent d*j-d+2 is at least j for j>=1,d>=3. Therefore every H_j has valuation>=j on either chart. Signed polynomial cofactors give the same scalar relation,whose higher coefficients retain H_j. The analytic nonzero-zeroth period exclusion is unchanged by multiplication of the operator and the exact Euler gauge.

If I=1+O(t^(K+1)) with K equal to the zeroth polynomial degree bound,each higher term has order at least j+(K+1-j)=K+1. It cannot cancel the nonzero zeroth polynomial. Thus33/78 are sufficient complement gates without an extra derivative-order allowance. No scalar factor is cancelled.

## 4. The adapted frame quartic is covered by the same calculation

The surviving repeated-line frames are the retained affine-W face92,the c=0 face35,and the c*d!=0 restricted quartic with

    z(u,t)=(u-t*(a0+b*u+d0*u^2))/(c*t),
    Phi=t*(alpha*u^3+G2(u,z))-z^2/2,
    normalized prefactor=constant/t.

Put U=t^2 Phi,D=t^2. Then deg_t U<=3 and the actual amplitude degree is4. Its quartic coefficient is

    [u^4]U=t^2*d0^2*(t*[T^2]G2-1/2)/c^2,

which is never identically zero. The prefactor is exactly D^(-1/2) up to a constant,and D(0)=0. The full preceding Euler calculation therefore applies to this actual frame family with the same gate78. The retained mixed-phase analytic exclusion supplies its nonzero zeroth scalar coefficient. It is not inferred from degree alone.

All frame faces are retained:affine-W92,c=0 gate35,and restricted interior78. Therefore the full repeated-line gate is

    max(33,78,92,35,78)=92.

The sharp repeated-class mixed cutoff m>3deg Q remains unchanged. No claim is made that every frame is null or that92 is the least pure prefix.

## 5. Current cover,literature refresh and assurance

The full cover is now quadratic7,smooth94,node96,cusp95,transverse96,tangent95,triangle96,concurrent92,repeated92. Its maximum96 and all least-prefix intervals remain unchanged. Odd repeated92 follows by inclusion;general/odd/homogeneous upper96 does not determine their least thresholds. The current radical-search target remains I96,with at most7 exact monotone rounds in principle and no full coefficient-ring search executed.

The checker verifies14 full time-dependent Euler gauge identities,21 principal identities,16 original twisted seed identities,six full scalar cofactor relations with22 coefficient checks,and every actual-rank/row degree and scalar valuation bound. Four phase families retain the earlier cubic/quartic elimination models;two additional actual adapted quartics include both nonzero and zero quadratic T^2 coefficient. All new inputs passed. No old checker main is rerun.

The [selected primary-record refresh](repeated92_primary_record_refresh_20261115.json) rechecked abstracts and submission histories for [Long's Gaussian counterexamples](https://arxiv.org/abs/2607.18186),[Long's xz/SU(2) counterexamples](https://arxiv.org/abs/2607.19012),[Wilson's two-variable theorem](https://arxiv.org/abs/2607.23887),[Dvorsky's GVC counterexample](https://arxiv.org/abs/2608.07338) and [Wang--Zhang's spatial eigenfunction moments](https://arxiv.org/abs/2610.16157). All five still list v1. The targeted search is not exhaustive and this refresh is not a fresh body/proof read. Their earlier body audits remain separate;no new external theorem or package is imported into this refinement.

Universal all-j degree/reduction/minor induction,complete actual elimination and adapted frame coverage,the seed-to-period bridge and analytic nonzero-zeroth arguments remain independently unreviewed written mathematics. Finite exact checks do not certify the complete analytic theorem. Every earlier proof body and receipt remains byte-identical.

Least general[20,96],odd[20,96],homogeneous[14,96],20/14 sufficiency,expanded I96 cofactors,independent analytic review,exceptional schemes,global term minima and lower-power GVC remain open. Conic/full-W-linear96 and its node/transverse/triangle frame branches now control the numerical maximum,followed by rank-one95 and smooth94. Qualitative cubic coverage and both lower witnesses remain closed internally.
