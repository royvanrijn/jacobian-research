# Formal Gaussian projection and a general sufficient prefix109

15 November 2026. Internal written refinement of the [general115 proof](TANGENT_POLE_RESIDUE_GLOBAL115.md), retaining the original Gaussian functional, complete coefficient/frame coverage and analytic nonzero-zeroth arguments.

Every complex ternary polynomial of total degree at most three has all positive pure Gaussian moments zero if and only if its first **109** pure moments vanish. The least general prefix is in **[20,109]**. Odd **[20,98]**, homogeneous **[14,96]**, smooth-leading **[20,94]**, homogeneous-smooth **[14,56]**, conic/full-W-linear96 and both complete first-nonzero20/14 witnesses remain unchanged. General109 does not give108 by parity for arbitrary non-odd polynomials.

The improved linear-D complements are node87,transverse97,triangle98,cusp77,tangent66. Their separately retained nonzero constant-D bounds101,106,103,109,105 now dominate. In particular the full cusp gate109 controls the general cover. This is a sufficient bound, not optimality.

The [checker](verify_gaussian_regular_projection.py), [certificate](gaussian_regular_projection_certificate.json), [receipt](gaussian_regular_projection_verification.json) and [source audit](gaussian_regular_projection_source_audit.json) bind auxiliary exact evidence. The all-order formal projection and regular finite connection, complementary-minor proof, polynomial extension, full frame cover and analytic exclusions remain written, independently unreviewed inputs. No independent replay, external review, formal verification or expanded radical certificate is asserted.

## 1. Project the completed square formally, while retaining the missing direction

All five complements have

    W=F^2/2+D*G, delta=t^2*D, D_prime=eps.

Work over the generic coefficient field K, with t and the indicated nonzero chart coefficients inverted, and with independent formal D. Coefficient derivatives remain included. The ordinary cases have F=u*v+r. The cusp and tangent cases have F=u^2+p*u+q*v+r. Their full G, Gaussian means and cross coefficients are exactly those in the retained proofs.

Set z=F. For ordinary cases v=(z-r)/u and the volume is du*dz/u. For cusp/tangent v=(z-u^2-p*u-r)/q and the volume is du*dz/q. Write H(u,z)=G(u,v(u,z)), H0=H(u,0). Here H is polynomial of z degree at most3, with Laurent-polynomial u coefficients in the ordinary cases and polynomial u coefficients in the other two.

For any amplitude f define the formal Gaussian projection

    J_D(f)=sum_(k>=0) (-t^2 D/2)^k/k!
               * [partial_z^(2k)(exp((H-H0)/t^2)*f)]_(z=0).

Each coefficient of D is a finite polynomial or Laurent polynomial in u. Thus this is well-defined coefficient by coefficient; it requires neither convergence nor an analytic contour. The sign is fixed by the original positive exponent exp(z^2/(2t^2D)): its formal variance is -t^2D.

Formal Gaussian integration by parts gives, for every amplitude,

    J_D((z+D*H_z+t^2D*partial_z)f)=0,
    J_D((H_u+t^2*partial_u)f)
                           =(H0_u+t^2*partial_u)J_D(f).

These follow term by term from the defining heat series. For example, differentiation of the exponential in z supplies H_z/t^2, and the Gaussian relation moves multiplication by z to -t^2D partial_z. Differentiation in u commutes with each finite coefficient and includes every derivative of H. For ordinary Laurent coordinates the volume du/u is retained: the remaining exact Euler operator is t^2*u*partial_u+u*H0_u. Replacing that volume by du would give the wrong divergence.

The two-variable twisted relations with delta=t^2D therefore map to the corresponding one-variable twisted relations after coefficientwise projection. The transverse relation maps to zero; the other relation maps to D times the one-variable relation, with the coordinate-Jacobian divergence retained. This remains an identity for Laurent-series witnesses obtained by the old bounded denominator clearing.

The normalized transverse Gaussian factor cancels the D^(-1/2) from the original completion. More explicitly, the original period prefactor is t^(-2)D^(-1/2), and transverse formal integration contributes (t^2D)^(1/2). The remaining factor is t^(-1), together with the Jacobian1/u or1/q. This cancellation includes the determinant and ordinary-volume derivatives already present in the full source. It does not change the original Gaussian functional.

One can also check the cancellation directly, for H independent of D:

    partial_D J_D(f)=-J_D(f)/(2D)-J_D(z^2*f)/(2*t^2*D^2)

for D-independent f. The apparent negative powers cancel coefficient by coefficient. D-dependent amplitudes add J_D(partial_D f). Physical time derivatives additionally differentiate t,H,p,q,r and the Jacobian; all their coefficients are regular in D on the generic chart. No denominator derivative or Euler/shear term is omitted.

## 2. The finite one-variable quotient has a regular connection

Reduce J_D coefficient by coefficient in the one-variable twisted quotient. At D=0 the relevant leading terms of H0 are:

| Complement | H0 endpoints or leading term | Finite quotient rank |
|---|---|---:|
| Node | a*u^3+d*(-r)^3*u^(-3) | 6 |
| Transverse | a*u^3+Lf*r^2*u^(-2) | 5 |
| Triangle | Le*u^2+Lf*r^2*u^(-2) | 4 |
| Cusp | -u^6/q^3 | 5 |
| Tangent | u^5/q^2 | 4 |

On this generic field the endpoint coefficients are units. In the polynomial cases, division by t^2 partial_u+H0_u reduces to bases1,...,u^4 or1,...,u^3. The derivative term has lower degree than the leading multiplication term. In the Laurent cases, the operator t^2 u partial_u+u H0_u reduces to any consecutive interval of p+q exponents for nonzero endpoints(-q,p); the upper and lower reductions strictly shorten the exponent range. Bases used here are [-3,2],[-2,2],[-2,1], respectively. These explicit degree divisions establish a finite free quotient of the displayed rank. They are applied coefficientwise in D, so no D denominator is introduced.

Denote the resulting projection by pi. The residual phase H0/t^2 and normalized prefactor have coefficients regular in D. Differentiating them and then performing these same one-variable divisions defines a finite connection nabla that preserves K[[D]]^r. This assertion means that differentiation of an analytic formal vector again has no negative D powers; it does not mean differentiation preserves its positive order. The commutator of the twisted derivative with the phase time derivative proves the connection descends to the quotient, or equivalently one differentiates the preceding integration-by-parts identities coefficientwise.

The original Euler fields and the tangent shear are exact two-variable integration-by-parts terms. Thus, for the complete induced seed operator L and its actual h,

    pi(L f)=h*nabla(pi(f)).

Ordinary/cusp h=t^3D^2. For the sheared tangent operator h=t^3D^3: the additional whole factorD is retained. This is essential to the tangent exponent below.

The actual numerator recurrence is R0=1 and

    R_(j+1)=L R_j-j*h_prime*R_j.

By induction the h_prime terms cancel exactly:

    pi(R_j)=h^j*nabla^j(pi(1)).

All iterated connection vectors are analytic formal series. Hence pi(R_j) contains D^(2j) for ordinary/cusp and D^(3j) for tangent. The argument uses the complete differential connection, rather than eigenvalues of a commutative Jacobi algebra.

This projection has rank one less than the original twisted quotient. We do not discard its missing direction. In the determinant argument a further analytic row restores that direction, with a precisely bounded determinant loss.

## 3. The comparison determinant costs0,2,6 powers, once per minor

The actual polynomial columns N_j and old scalar denominators H_j are unchanged:

| Complement | Completed columns | Scalar H_j | Slow projected order |
|---|---|---|---:|
| Ordinary | old axis-cleared n_j | h^j times old axis clearing | 2j |
| Cusp | D^(2j-1)*n_j | t^(3j)D^(4j-1) | 4j-1 |
| Tangent | D^(3j)*n_j | t^(3j)D^(6j) | 6j |

The old ordinary axis factors are polynomial and are units at generic D=0; they cannot lower these orders. No change is made to the earlier normal reductions or actual twisted identities.

Make a square comparison matrix Q by adding one unrestricted analytic row to the finite projection matrix pi(E), where E is the old completed basis. The extra row selects the unit coordinate in the ordinary cases, v in the cusp case and v^2 in the tangent case.

For ordinary cases at D=0, restriction to uv=-r identifies only uv with -r among the retained basis entries. In the displayed Laurent bases the other entries are independent. The extra unit row distinguishes that relation, so det Q is a D-adic unit. All old axis row scales are generic units as well. The determinant loss is0.

For cusp E=(1,u,u^2,v,Duv,Du^2v). Remove the v column using the extra row. At D=0 put v0=-(u^2+pu+r)/q. The five finite columns have leading terms

    1,u,u^2,D*u*v0,D*u^2*v0.

Their determinant in1,u,u^2,u^3,u^4 has leading term D^2/q^2. Higher D terms cannot change this nonzero leading coefficient. The loss is2.

For tangent E=(1,Du,D^2u^2,D^2v,D^4v^2). Remove the v^2 column. Replace the v column by its combination with the first three that gives D^2*z/q. Gaussian integration by parts gives exactly pi(z)=-D*pi(H_z). At z=0 the full phase uv^2+G2 has

    H_z=-2*u^3/q^2 + terms of u degree at most2.

Thus the four finite columns have independent leading degrees0,1,2,3 with D orders0,1,2,3, and det Q has leading term2D^6/q^3. The full G2 changes only lower u degrees. The loss is6. Here the extra row remains unrestricted and may have any analytic seed values; it is not assumed to be a formal eigenvector.

For every selected s-row minor, undoing Q loses at most kappa=ord_D det Q, once. Indeed complementary-minor identities express every s-minor of Q^(-1) as an analytic complementary minor of Q divided by det Q. Cauchy--Binet then transfers the transformed lower bound. A bound that charged kappa for every row would be unnecessarily weak.

## 4. Stronger common scalar factors retain every higher valuation

Put w0=0. For j>=1 set wj=2j,4j-1,6j in the ordinary,cusp,tangent cases, respectively. For a raw cofactor C_l omitting column l from N0,...,Ns, let B_l be the sum of the smallest s-1 weights among the remaining column indices, minus kappa, clipped at0. In every transformed minor at most one row can use the unrestricted added coordinate. All other rows have the stated slow orders. Consequently ord_D C_l>=B_l for every selected row set.

Combine this with the retained triangular bounds T_s=s(s-1)/2 for C0 and max(0,T_s-l) for higher C_l. Write b_l for the maximum of the two bounds. These include all smaller ranks, even where the projection loss makes its estimate weaker.

The common scalar factor is D^(b0). The scalar coefficient of order l>=1 is H_l*C_l. For every rank and l,

    b_l+ord_D H_l >= b0.

For ordinary cases b0=2T_s. For full cusp rank6 b0=4T6-5-2=53. For full tangent rank5 b0=6T5-6=54. The complete smaller-rank/row inequalities are included in the certificate. Dividing this factor from the entire relation leaves every higher coefficient divisible by t^(3l), because only D factors were cancelled. This remains true when D(0)=0 or an axis has a common t/D root. The analytic nonzero zeroth coefficient remains nonzero.

All formal work above takes place on a generic coefficient chart. The old columns, cleared cofactors and their scalar coefficients are polynomial in D and the original coefficient parameters. Their low D coefficients vanish over the generic field, hence identically. The divisibilities therefore extend to q=0,r=0,degenerate comparison determinants,zero traces,D(0)=0 and shared roots. No q,r,t or trace factor is cancelled from the actual scalar relation. Identically zero axes and identically zero D stay on their separately retained frame charts. Generic inversion alone would not justify those cases without this polynomial extension.

## 5. Complete cover and assurance

Subtract the new common powers from the retained raw full-rank degrees:

| Complement | Raw degree | New D power | Linear-D gate | Retained constant-D gate | Full class |
|---|---:|---:|---:|---:|---:|
| Node | 129 | 42 | 87 | 101 | 101 |
| Transverse | 127 | 30 | 97 | 106 | 106 |
| Triangle | 118 | 20 | 98 | 103 | 103 |
| Cusp | 130 | 53 | 77 | 109 | 109 |
| Tangent | 120 | 54 | 66 | 105 | 105 |

Every smaller actual rank and selected row set is checked; its maximum is the displayed linear-D gate. No degree saving is claimed for constant D. All parent/conic/rank-one and zero-axis frame gates remain at most96 and below these full-class maxima.

The full cover is quadratic7,smooth94,node101,cusp109,transverse106,tangent105,triangle103,concurrent104,repeated99. Its maximum109 gives the general test with no extra Gaussian-frame condition. Pure nullity retains the internally closed classification and sharp mixed cutoff m>3deg Q. Odd98 and homogeneous96 remain separate.

In the rational twenty-coefficient ring radical(I109)=radical(I_infinity). This asserts neither ordinary ideal equality nor expanded radical cofactors. At most7 exact monotone radical-search rounds against I109 determine the least global prefix in principle. None has been executed on the full coefficient ring; equality of consecutive ideals would not suffice.

The checker verifies119 formal Gaussian coefficient identities through D order4, five generic comparison determinants, and15 retained original Gaussian cofactor families with100 full scalar-coefficient checks. These include D(0)=0,zero q and shared q/D or axis/D roots. Retained original twisted identities, source audits and constant-D proofs are hash-bound; old checker mains are not rerun. Initial cusp order and cross-family probes are preserved as motivation. No checker failure occurred in the accepted new input.

The universal coefficientwise formal integration, regular finite connection, comparison determinant/complementary-minor argument, polynomial extension and complete analytic/frame cover are written inputs. Finite checks support their algebra but do not certify the complete analytic theorem. All earlier proof bodies and receipts remain byte-identical. No fresh literature-coverage claim or external theorem is imported.

Least global[20,109],odd[20,98],homogeneous[14,96],20/14 sufficiency,expanded I109 cofactors,independent analytic review,exceptional schemes,global term minima and lower-power GVC remain open. Constant-D cusp109,transverse106 and tangent105 now control the general numerical frontier. Qualitative cubic coverage and both lower witnesses remain closed internally.
