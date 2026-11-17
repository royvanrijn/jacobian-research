# Two regular rank-one branches: rank-one57 and general94

15 November 2026. Internal written refinement of the [general95/conic94 proof](CONIC_CONSTANT_TERM_GLOBAL95.md). Retain the original Gaussian Wick functional,complete lower-degree/frame coverage and analytic nonzero-zeroth exclusions.

Every complex ternary polynomial of total degree at most three has all positive pure Gaussian moments zero if and only if its first **94** pure moments vanish. Odd and homogeneous sufficient prefixes remain **94**. Least intervals are general **[20,94]**,odd **[20,94]**,homogeneous **[14,94]**;smooth-leading **[20,94]**,homogeneous-smooth **[14,56]** and both complete first-nonzero20/14 witnesses stay unchanged. These are sufficient bounds,not optimality.

The rank-one W-linear family improves95 to **57**. With conic94 and parent92,full W-linear is **94**. Smooth94 and conic94 now control the general maximum. The [checker](verify_rank_one_two_branch_projection.py),[certificate](rank_one_two_branch_projection_certificate.json),[receipt](rank_one_two_branch_projection_verification.json) and [source audit](rank_one_two_branch_projection_source_audit.json) bind auxiliary exact algebra. Universal formal branch quotients and regular connections,complementary-minor/joint-degree proof,polynomial extension and complete analytic/frame coverage remain independently unreviewed written inputs.

## 1. Retain the original fixed-point phase and its moving quadratic

For P=Wc*(a+bZ+cT+fT^2)+G(Z,T),f!=0,deg G<=3,the retained Gaussian--Lagrange identity has q=1-bt and normalized one-variable phase

    I=(1/q)*integral exp(t*G(t*A(T)/q,T)-T^2/2) dT,
    A(T)=a+cT+fT^2.

This is the same formal Wick period used in the [rank-factor proof](RANK_FACTORS_AND_REDUCTION_GLOBAL_130.md),not a new Gaussian average. After y=tT it has prefactor1/(tq),phase W0/(t^2q^3),with

    ahat=a*t^2+c*t*y+f*y^2,
    W0=sum_(i+j<=3) g_ij*t^(3-i-j)*q^(3-i)*ahat^i*y^j-q^3*y^2/2.

Only its actual degree6 case needs a new bound. Then g30!=0 and the y^6 coefficient f^3*g30 is a nonzero constant in time and q. Its full ordinary twisted rank is5. Keep the already normalized complete operator

    h=t^3*q^4,
    S=t*q*(W0)_prime-(2q-3bt)*W0-t^2*q^4+b*t^3*q^3,
    R0=1,R_(j+1)=h*(R_j)_prime-j*h_prime*R_j+S*R_j.

All q,denominator and normalization derivatives remain included. Work first with independent q and compatible time derivation partial_t-b partial_q.

## 2. Two formal branches remove the apparent q pole

On the generic coefficient field let r1,r2 be distinct roots of A. The equations

    A(T_i(q,Z))=q*Z/t, T_i(0,Z)=r_i

have unique coefficientwise formal solutions. Their coefficients are polynomials in Z with coefficient-field denominators involving t,f and r1-r2. They introduce no negative q powers. In the original phase the variable t*A(T)/q is exactly Z,so the branch phase becomes

    t*G(Z,T_i(q,Z))-T_i(q,Z)^2/2.

It is analytic in q. The normalized volume also loses its apparent pole:

    (1/q)*(partial T_i/partial Z)=1/(t*A_prime(T_i)).

Its leading coefficient is the nonzero unit1/(t*A_prime(r_i)). Keep this complete Jacobian and the Gaussian T_i^2 term. Every time derivative,including q_prime=-b and implicit branch derivatives,has coefficients analytic in q. A contour or convergence assertion is not required.

At q=0 each phase is t*G(Z,r_i)-r_i^2/2,a cubic in Z with nonzero leading coefficient t*g30. The ordinary one-variable twisted quotient therefore has basis1,Z and rank2. Division by its quadratic leading derivative includes the lower derivative correction and has no q denominator. Correct the quotient coefficient by coefficient in q:each correction is a finite polynomial in Z,so this gives a complete free formal rank-two quotient on each branch.

Equivalently,factor the phase exponential into its q=0 cubic exponential and a formal analytic perturbation. At each q order the perturbation and full Jacobian have polynomial Z amplitudes. Their original ordinary divergence identities map to the corresponding full cubic twisted relations. Differentiating these identities gives a time connection preserving analytic formal vectors. Thus the two branches provide four finite projected rows with a regular time connection. They do not replace the full rank-five quotient by a four-dimensional sample:the missing direction is restored below.

For each branch projection pi_i,the full original seed operator intertwines with h times this regular connection. Exact cancellation of the h_prime terms yields

    pi_i(R_j)=h^j*nabla_i^j(pi_i(1)).

Consequently every one of these four projected seed coordinates has q order at least4j. The original normal columns N_j in1,y,...,y^4 are polynomial in q;their optional endpoint-clearing factors are constants,so the same projected orders apply to them.

## 3. The comparison determinant loses precisely two q orders

Normalize each branch projection by its leading Jacobian/phase unit. At q=0 its constant row on T^k is r_i^k,while its Z row starts at positive q order. To first order,

    T_i=r_i+q*Z/(t*A_prime(r_i))+O(q^2),
    pi_i(T^k)=r_i^k*pi_i(1)
                 +q*k*r_i^(k-1)*Z/(t*A_prime(r_i))+O(q^2).

An analytic row operation subtracts the first-order Z coordinate of pi_i(1) times its constant row. It leaves the leading Z row displayed by the derivative term. All such normalizations are analytic units;none divides by q.

Add a fifth unrestricted analytic row selecting the y^4 coordinate. Remove that column to compute the leading comparison determinant on k=0,...,3. With y^k=t^k T^k,the remaining four rows are the two value rows at r1,r2 and their two derivative rows multiplied by q/(t A_prime(r_i)). The confluent Vandermonde determinant is

    plus/minus q^2*t^4*(r1-r2)^2/f^2.

It is nonzero generically,so the complete comparison matrix Q has determinant q order2. Higher q terms cannot cancel this leading coefficient. Complementary-minor identities charge this loss once per selected minor. The unrestricted fifth row has not been discarded or assumed to be a formal eigenvector.

## 4. Cancel the common q factor at every actual rank

At actual seed rank s<=5,selected minors of columns N1,...,Ns have four slow rows of orders4j and at most one unrestricted row. Put the latter on the largest column index. The determinant order is at least

    b0=max(0,4*T_s-2), T_s=s(s-1)/2.

For the cofactor omitting column l from N0,...,Ns,set w0=0,wj=4j. A safe bound b_l is the sum of the smallest s-1 weights among the remaining indices,minus2,clipped at0. This also handles every smaller actual rank and every unit-containing row selection through Cauchy--Binet and complementary minors.

The scalar coefficients are C0 and H_l*C_l,where H_l=h^l times optional constant endpoint factors. For every s,l>=1,

    b_l+4l>=b0.

Thus the whole scalar relation contains q^b0. Divide only that factor. Higher coefficients retain t^(3l) and nonnegative q powers. The original analytic nonzero zeroth coefficient remains nonzero. At full rank5,b0=4*10-2=38. No t,root-discriminant,f or g30 factor is cancelled from the scalar relation.

The generic branch proof inverted the quadratic discriminant. After bounded normal clearing,the raw cofactors are polynomial in independent q and all original coefficient parameters. Their low q coefficients vanish on the dense distinct-root coefficient field,and hence identically. Divisibility therefore extends to merging roots,c^2-4af=0,as well as every special lower coefficient. Degree drops g30=0 belong to the separately retained lower actual-degree strata. The f=0 face is the retained parent family. Physical q=1-bt is never identically zero,and q(0)=1. The b=0 case follows by the same polynomial extension and joint-degree specialization below;it is not assigned a nonexistent q root.

## 5. Divide in the joint-graded ring before setting q constant

Give t,y,q weight1,with original coefficients weight0 and q_prime=-b. The retained exact W0 has joint degree6,and the full source S has joint degree7. The actual R_j has joint degree<=7j. The leading y derivative coefficient is a nonzero constant,so complete ordinary twisted reduction preserves joint degree and gives normal row bounds7j-k.

At rank s the raw zeroth cofactor joint degree is at most7*s*(s+1)/2-sum_selected k. Its common q^b0 is a homogeneous monomial in this ring. Polynomial division therefore lowers joint degree by b0 before q is specialized. At full rank5,

    7*15-10-38=57.

Every smaller actual rank and selected row set has no larger value. Substitution q=1-bt,including constant q=1 at b=0,cannot increase joint degree. This is not subtraction of the degree of a constant factor from a univariate bound known only after specialization. The constant-q exact model is included in the certificate.

The actual degree2/3/4/5 strata keep their earlier gates9/21/51/54. They already cover all degree/leading-coefficient boundaries and the b=0,G=0 pure-null face,without applying a nonzero-zeroth exclusion to that null family. Consequently the full rank-one family has sufficient gate max(9,21,51,54,57)=57.

## 6. Complete general94 cover and assurance

Conic94,parent92 and rank-one57 give full W-linear94. Combine this with the retained complements and frame covers:

    quadratic7,smooth94,node94,cusp92,transverse94,tangent92,
    triangle94,concurrent92,repeated92.

The maximum94 proves the general pure-only gate with no additional frame hypothesis. Odd94 and homogeneous94 were already retained and remain unchanged. Even94 gives no93 parity deduction. Sharp mixed cutoffs and the qualitative cubic GMC classification remain unchanged. In the rational twenty-coefficient ring radical(I94)=radical(I_infinity),without ordinary ideal equality or expanded cofactors. At most7 exact monotone radical-search rounds against I94 determine the least general prefix in principle;none has run on the full coefficient ring.

The checker verifies two formal branch equations through q^4 and their normalized Jacobians through q^3,the exact original phase rescaling and joint source bounds,the leading Hermite determinant,15 original twisted seed identities and three full scalar relations with18 coefficient checks. Models include nonconstant q,constant q and a merged-root quadratic. Their zeroth cofactors are nonzero and every post-cancellation scalar valuation and actual-rank/row degree inequality passes.

The [initial input](rank_one_two_branch_projection_initial_input.py) and [failure receipt](rank_one_two_branch_projection_initial_failure.json) are preserved. It compared expanded rational expressions syntactically;the apparent residual had every q^0,...,q^4 coefficient identically zero after exact common-denominator cancellation. The [diagnostic](rank_one_two_branch_projection_initial_diagnostic.json) records this. The accepted assertions apply exact rational cancellation before comparison. No mathematical hypothesis was weakened and no counterexample was discarded.

Universal coefficientwise branch quotients and regular time connections,the comparison/minor argument,joint-degree specialization,polynomial extension,complete lower-degree/frame coverage and analytic nonzero-zeroth exclusions remain independently unreviewed written inputs. Finite checks support auxiliary algebra and do not certify the complete analytic theorem. All earlier proof bodies and receipts remain byte-identical. No fresh literature-coverage claim,external theorem or package is imported;the preceding selected primary-record refresh remains retained.

Least general[20,94],odd[20,94],homogeneous[14,94],20/14 sufficiency,expanded I94 cofactors,independent analytic review,exceptional schemes,global term minima and lower-power GVC remain open. Smooth94 and conic94 now control the common numerical maximum. Qualitative cubic coverage and both lower witnesses remain closed internally.
