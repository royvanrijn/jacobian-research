# Formal Gaussian projection at infinity: general sufficient prefix104

15 November 2026. Internal written refinement of the [general109 proof](FORMAL_GAUSSIAN_PROJECTION_GLOBAL109.md), retaining the original Wick functional, complete coefficient/frame cover and analytic nonzero-zeroth exclusions.

Every complex ternary polynomial of degree at most three has all positive pure Gaussian moments zero if and only if its first **104** pure moments vanish. The least general prefix is in **[20,104]**. Odd **[20,98]**, homogeneous **[14,96]**, smooth-leading **[20,94]**, homogeneous-smooth **[14,56]**, conic/full-W-linear96 and both complete first-nonzero20/14 witnesses remain unchanged. Even104 gives no103 parity deduction for arbitrary polynomials.

The new nonzero constant-D complements are node80,transverse91,triangle93,cusp71,tangent61. Retaining the previous linear-D bounds and every frame face gives full node96,transverse97,triangle98,cusp95,tangent95. Concurrent104 now controls the general cover, followed by repeated99. No optimality is asserted.

The [checker](verify_gaussian_projection_at_infinity.py), [certificate](gaussian_projection_at_infinity_certificate.json), [receipt](gaussian_projection_at_infinity_verification.json) and [source audit](gaussian_projection_at_infinity_source_audit.json) bind auxiliary exact evidence. All-order formal Morse integration, regular time connection, comparison/minor proof, polynomial extension and complete analytic/frame cover remain written, independently unreviewed inputs. No independent replay, external review, formal verification or expanded radical certificate is asserted.

## 1. Rescale the full phase, including all lower terms

Throughout this new argument D is a nonzero constant and D_prime=0. Put x=1/t, u=U/x, v=V/x. For the actual phase W=F^2/2+D G define

    Fhat=x^2 F(1/x,U/x,V/x), Ghat=x^3 G(1/x,U/x,V/x).

All coefficients are analytic polynomials in x. This follows from the retained degree bounds: p,q have degree<=1,r<=2, quadratic coefficients of G have degree<=1, linear coefficients<=2 and its constant<=3. The full means, cross coefficients and lower terms are retained. The exponent becomes

    W/(t^2D)=Fhat^2/(2D x^2)+Ghat/x.

Generic leading coefficients are inverted temporarily. Ordinary Fhat=UV+rhat, with rhat=x^2r(1/x); cusp/tangent Fhat=U^2+phat U+qhat V+rhat, where phat=xp(1/x),qhat=xq(1/x). Set z=Fhat and H=D Ghat in these coordinates. The exponent is

    (z^2/2+x H(U,z,x))/(D x^2).

The coordinate volume is dU dz/U for ordinary cases and dU dz/qhat for cusp/tangent. The original prefactor t^(-2)D^(-1/2) and the scaling Jacobian x^(-2) leave D^(-1/2). These identities include the original ordinary-volume and Gaussian-determinant terms.

## 2. Center before applying the formal Gaussian heat series

A naive heat expansion of exp((H-H0)/(D x)) at z=0 has infinitely many terms of the same x order. We instead center the transverse critical point formally. The equation

    phi+x H_z(U,phi,x)=0

has a unique phi in x K[U][[x]] (Laurent U in ordinary cases), by coefficient recursion with unit linear coefficient. Its coefficients are polynomial/Laurent in U. Define

    Phi(U,x)=H(U,phi,x)+phi^2/(2x), Phi(U,0)=H(U,0,0).

It is analytic in x. Since H has z degree at most3, after z=phi+w the phase numerator has the exact form

    z^2/2+xH=x Phi+(1+xH_zz(U,phi,x))*w^2/2
                                      +xH_zzz(U,phi,x)*w^3/6.

The formal Morse coordinate

    y=w*sqrt(1+xH_zz(U,phi,x)+xH_zzz(U,phi,x)*w/3)

has Jacobian1 at x=0 and has a unique inverse w=y+O(x). For each x coefficient these changes and their Jacobians are polynomials in y and polynomial/Laurent in U. The transformed exponent is exactly

    Phi/(D x)+y^2/(2D x^2).

Apply the normalized formal Gaussian heat series of variance -D x^2 to the transformed amplitude and its full coordinate Jacobian. It is now well-defined coefficientwise: at each x order only finitely many amplitude coefficients and Gaussian moments contribute. The Jacobian1/U or1/qhat is retained. The transverse normalized factor contributes x*D^(1/2), leaving prefactor x times that Jacobian. Its logarithmic time derivative is regular; in particular d(log x)/dt=-x.

This is coefficientwise formal integration, not a convergence or contour assertion. Gaussian integration by parts follows term by term. Coordinate changes retain divergence with the stated volume. In the original z coordinate it gives exactly

    pi((z+xH_z+D x^2 partial_z)f)=0.

In particular pi(z)=-x pi(H_z). This identity will determine the tangent comparison determinant; the centered shift must not be dropped.

## 3. The finite quotient connection is regular for d/dt

The residual one-variable twisted operators are

    D x partial_U+Phi_U,
    D x U partial_U+U Phi_U  (with volume dU/U).

Their x=0 multiplication terms have respectively degree5 for cusp and4 for tangent, and Laurent endpoint spans6/5/4 for node/transverse/triangle. The leading terms of Phi0 are D times those in the [previous projection proof](FORMAL_GAUSSIAN_PROJECTION_GLOBAL109.md): cusp -U^6/qhat0^3,tangent U^5/qhat0^2,ordinary endpoints(3,-3),(3,-2),(2,-2), with nonzero generic endpoint coefficients.

Consequently these formal quotients have bases1,...,U^4 for cusp,1,...,U^3 for tangent and Laurent intervals[-3,2],[-2,2],[-2,1] for the ordinary cases. To construct the projection, divide at order x^0 by the monic leading multiplication operator, retain the displayed basis, and correct the higher x coefficients recursively. At each order all remainders are finite polynomials/Laurent polynomials. Multiplication endpoints are units and the derivative term has positive x order. This gives a free complete formal module of ranks5,4,6,5,4 without introducing negative x powers. No finite commutative sample replaces the original twisted quotient: the multiplication quotient is used only to initialize this full formal deformation.

The physical time derivation is -x^2 partial_x. Its action on the residual phase is

    d(Phi/(D x))/dt=(Phi-x Phi_x)/D,

which is regular in x. All normalized-prefactor,coordinate-Jacobian,basis and coefficient derivatives are regular as well. Differentiation of the integration-by-parts identities, or their direct commutator, shows that this descends to a finite connection nabla_t=-x^2 partial_x+B(x), with B analytic. Thus nabla_t preserves analytic formal vectors. This uses the physical time derivation; a partial_x connection alone would have an irregular pole.

The full old numerator operator L intertwines with h*nabla_t under the projection. All Euler fields and the tangent shear are retained exact integration-by-parts terms. Original constant-D h has degree3 in t in all cases: t^3D^2 for ordinary/cusp and the separately normalized unsheared t^3D for the retained tangent constant-D seeds. The exact seed recurrence R_(j+1)=L R_j-j h_prime R_j therefore gives

    pi(R_j)=h^j nabla_t^j pi(1).

The h_prime contributions cancel exactly. Hence the projected original numerator has x pole order at most3j. Neither the Euler terms nor the ordinary volume/determinant derivatives can be omitted from this identity. Simultaneously rescaling u,v changes the period Jacobian as above; its derivative is included in nabla_t.

## 4. Compare with the unchanged polynomial seed rows

Retain every old constant-D polynomial column N_j and scalar denominator H_j. Write its row-degree estimate as

    deg_t N_(i,j)<=c*j-o-g_i,
    X_j=x^(c*j-o) diag(x^(-g_i)) N_j(1/x).

Then X_j is analytic at x=0. The old row parameters are:

| Case | c | o | g_i |
|---|---:|---:|---|
| Node | 4 | 0 | 0,1,1,2,2,2,3 |
| Transverse | 6 | 2 | 0,1,1,2,2,2 |
| Triangle | 8 | 0 | 0,3,3,6,5 |
| Cusp | 6 | 1 | 0,1,2,1,3,4 |
| Tangent | 8 | 0 | 0,2,4,3,6 |

For ordinary cases the completed row basis has t degree g_i. The old column clearing has degree0 for node,2j-2 for transverse and4j for triangle. Combining these factors with the projected numerator pole3j gives slow projected order at least j in every ordinary case.

For cusp the row basis is (1,u,u^2,v,Duv,Du^2v). Its actual t scaling is the ordinary monomial degree; the difference g_i-|e_i|=(0,0,0,0,1,1) supplies comparison columns1,U,U^2,V,xUV,xU^2V, up to constant-D row units. Thus the projected X_j has order at least (6j-1)-3j=3j-1.

For tangent the analogous difference is rho=(0,1,2,2,4). Comparison columns are1,xU,x^2U^2,x^2V,x^4V^2, up to constant-D units. Its projected order is at least8j-3j=5j. Only constant row units are suppressed here; the actual columns and scalar relation remain unchanged.

Add one unrestricted analytic row to restore the missing quotient direction: select the unit coordinate for ordinary,v for cusp,v^2 for tangent. The comparison determinants have x orders0,2,6. Ordinary follows by restriction UV=-rhat0 and the independent Laurent bases. Cusp follows from V0=-(U^2+phat0 U+rhat0)/qhat0: after removing the v column the determinant of1,U,U^2,xUV0,xU^2V0 has leading term x^2/qhat0^2.

For tangent replace the V column by its combination with the first three giving x^2 z/qhat. The exact Gaussian identity pi(z)=-x pi(H_z) makes it -x^3 pi(H_z)/qhat. Here H_z(U,0,0) has leading term -2D U^3/qhat0^2, while all full quadratic lower terms have degree at most2. The determinant has leading term2D x^6/qhat0^3. The centered transverse shift is precisely what produces this nonzero fourth finite column.

Complementary-minor identities charge the comparison determinant order once per selected minor. This is the same complete-matrix argument as in the previous proof; it retains the unrestricted missing direction rather than deleting it.

## 5. Degree savings, valuations and specialization

For actual seed rank s, at most one row of a transformed determinant can use that unrestricted direction. Put it on column s. The s-1 slow columns supply an infinity-order saving of at least

    node/transverse/triangle: T_s=s(s-1)/2,
    cusp: max(0,3T_s-(s-1)-2),
    tangent: max(0,5T_s-6).

The same lower bound applies to every selected original row minor by Cauchy--Binet and complementary minors. Subtract it from the old raw degree c*s*(s+1)/2-o*s-sum_selected g_i. All smaller ranks and row selections containing the unit are included. At full rank the bounds are:

| Complement | Raw degree | Infinity saving | New constant-D gate |
|---|---:|---:|---:|
| Node,rank7 | 101 | 21 | 80 |
| Transverse,rank6 | 106 | 15 | 91 |
| Triangle,rank5 | 103 | 10 | 93 |
| Cusp,rank6 | 109 | 38 | 71 |
| Tangent,rank5 | 105 | 44 | 61 |

These are degree bounds on the original zeroth cofactor. No t factor is cancelled from any scalar relation. Higher scalar coefficients and their valuations>=3l at t=0 are unchanged, as is the retained nonzero zeroth coefficient. Thus an argument at infinity has not weakened the initial-order contradiction at zero.

The formal argument first inverts generic leading coefficients, including qhat0,rhat0 and the ordinary axis slopes. Every coefficient of the old cleared polynomial cofactor above the new bound vanishes over that generic field, hence identically in the coefficient ring. This polynomial extension retains vanishing leading coefficients,constant or zero q/r,zero axis slopes and shared roots. Axes that vanish identically and identically zero D still use their separate frame cases. No new leading-coefficient denominator is introduced into the actual scalar relation. Constant D is never treated as a degree saving.

## 6. Full cover and evidence boundary

Combine new constant-D gates80/91/93/71/61 with retained linear-D87/97/98/77/66. The full frame gates are conic96 in ordinary cases and parent92/rank-one95 in cusp/tangent cases. The resulting full cover is

    quadratic7,smooth94,node96,cusp95,transverse97,tangent95,
    triangle98,concurrent104,repeated99.

Its maximum104 proves the general sufficient prefix with no added frame condition. Odd98,homogeneous96,conic/full-W-linear96 and sharp mixed cutoffs stay unchanged. In the rational twenty-coefficient ring radical(I104)=radical(I_infinity), without ordinary ideal equality or expanded cofactors. At most7 exact monotone radical-search rounds against I104 determine least K in principle; none has been executed in the full coefficient ring.

The checker verifies the generic transverse critical equation through x^5,the formal Morse identity through x^4,the regular residual phase derivative,five comparison determinants,eight retained original constant-D Gaussian families with53 full scalar-cofactor checks,and all rank/row degree estimates. Original twisted seed certificates and previous sources/receipts are hash-bound; their checker mains are not rerun. No accepted-input failure occurred. Finite checks support auxiliary algebra and do not certify the universal formal or analytic theorem.

All-order coefficientwise formal Morse/Gaussian integration,regular physical-time connection,missing-direction comparison/minor argument,polynomial extension and complete analytic/frame coverage remain independently unreviewed written inputs. Every earlier canonical proof and receipt stays byte-identical. No fresh literature-coverage claim,external theorem or producer package is imported.

Least global[20,104],odd[20,98],homogeneous[14,96],20/14 sufficiency,expanded I104 cofactors,independent analytic/formal-connection review,exceptional schemes,global term minima and lower-power GVC remain open. Concurrent104 and repeated99 now control the general numerical frontier. Qualitative cubic coverage and both lower witnesses remain closed internally.
