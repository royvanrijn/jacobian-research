# Tangent pole residue and constant-D normalization: general prefix115

15 November 2026. Internal written refinement of the [triangular divisor/general120 proof](TRIANGULAR_DIVISOR_FACTORS_GLOBAL120.md), retaining the original Wick functional, full coefficient/frame cover and analytic exclusions.

Every complex ternary polynomial of total degree at most three has all positive pure Gaussian moments zero if and only if its first **115** pure moments vanish. The least general prefix is therefore in **[20,115]**. Odd **[20,98]**, homogeneous **[14,96]**, smooth-leading **[20,94]** and homogeneous-smooth **[14,56]** remain unchanged. Conic/full-W-linear96 and both exact first-nonzero20/14 witnesses remain retained. Odd115 does not imply114 for arbitrary non-odd P.

The new full tangent-line/smooth-conic gate is **110**:110 on linear D,105 on nonzero constant D, and retained92/95 frame faces. Cusp115 now controls the general cover.

The [checker](verify_tangent_rank_one_lattice.py), [certificate](tangent_rank_one_lattice_certificate.json), [receipt](tangent_rank_one_lattice_verification.json) and [source audit](tangent_rank_one_lattice_source_audit.json) support auxiliary exact algebra. Universal pole-lattice/seed/triangular-minor induction, polynomial extension across degenerate traces, the constant-D all-degree path estimate, full coefficient/frame coverage and analytic nonzero zeroth coefficients remain written, independently unreviewed inputs. No independent replay, external review, formal verification or expanded global radical certificate is asserted.

## 1. Keep the full tangent Gaussian square and weighted Euler shear

The actual rescaled tangent phase is

    W=D*(u*v^2+G2)+F^2/2,
    F=u^2+p*u+q*v+r, delta=t^2*D,
    G2=a*u^2+b*u*v+c*v^2+e*u+f*v+g.

Here D=m-2alpha*t, p,q have degree at most1,r at most2,a,b,c at most1,e,f at most2,g at most3. Full Gaussian means and cross coefficients are retained. Leading normalizations and all zero-divisor/frame charts are those of the previous tangent proof.

Retain bx=delta_prime/4,by=tD-delta_prime/8 and the complete unsheared weighted Euler source S_beta. Put C=[u^2*v]S_beta. The normalized sheared source and denominator are

    Sbar=(D*S_beta-C*u*W_v/2)/t, h=t^3*D^3.

Its full vertical field includes both Euler fields and the shear -delta*C*u*partial_v/(2t). No ordinary-volume divergence or denominator derivative is dropped.

## 2. A pole scaling exposes the rank-one residue

Work first with independent D,D_prime=eps and full derivatives of p,q,r and G2. Set u=U/D,v=V/D^2. The exact leading limits are

    D^4*W -> Psi=U*V^2+(U^2+q*V)^2/2,
    D^3*Sbar -> Sigma=(5*t*eps*q^2/16)
                                     *(U^3+q*U*V-2*V^2).

Both rescaled twisted gradients have derivative coefficient t^2*D^5, which vanishes in the residue. Their multiplication limits are Psi_U,Psi_V. Thus the residue projection is the full commutative Jacobi quotient of Psi in the basis1,U,U^2,V,V^2. The retained weighted principal sequence and strict lower-weight divisions give rank5 throughout the parameter family; no smaller sampled quotient is substituted.

The factor5/16 in Sigma includes by at D=0: by=-t^2*eps/8. That contribution does not vanish even though D does. Omitting it would give the incorrect1/4 coefficient.

Multiplication by Sigma has residue matrix A0=t*eps*A_base. All rows of A_base are zero except its fifth row,

    (-25*q^2/32,125*q^4/256,-625*q^6/2048,
                              -625*q^5/512,-15625*q^8/8192).

Consequently

    trace A0=-15625*t*eps*q^8/8192,
    A0^2=(trace A0)*A0,
    rank A0<=1.

All150 two-by-two minors of the augmented first-seed/residue matrix vanish. The trace is not identically zero, so the generic nonzero-trace hypothesis of the retained formal differential splitting lemma is satisfied. No q factor is cancelled from the actual scalar relation.

## 3. Include every derivative before taking the residue

Keep rho(i,k)=i+2k and the completed row lattice

    E=(1,D*u,D^2*u^2,D^2*v,D^4*v^2).

The previous source pole-grade bound3 and normal filtration show that N_j=D^(3j)*n_j is polynomial in D, where n_j is the actual numerator coordinate vector in E. Its scalar denominator is

    H_j=t^(3j)*D^(6j).

Applying the same filtration to each E column, whose pole grade is0, proves D^3 times the full induced source matrix analytic in D. This supplies an analytic matrix, not just a formal residue candidate.

Basis derivatives contribute h*rho*D_prime/D before scaling; the Euler and shear fields are also included. After the D^3 matrix scaling every such contribution has D order at least5. The derivative coefficient in the N_j recurrence is D^3*h=t^3D^6, and its additional scalar j terms

    -j*D^3*h_prime-3j*D^2*h*D_prime

have D order at least5. The only matrix residue is therefore the A0 above. N1 mod D=A0*unit, so it lies in the residue image.

The derivative/scalar orders6/5 satisfy the weaker2/1 hypotheses of the [triangular divisor lemma](TRIANGULAR_DIVISOR_FACTORS_GLOBAL120.md). Its full formal block gauge includes all derivative corrections. At every actual seed rank s, the zeroth raw minor has D^T_s and the higher ones D^max(0,T_s-l), where T_s=s(s-1)/2.

Generic nonzero trace is used only to prove these polynomial divisor identities before specialization. The polynomial completed columns and cofactors extend them to q=0,zero trace at D(0)=0 and shared q/D/t roots. Such a degenerate member is not directly diagonalized by dividing by its zero trace. Identically zero D remains a retained frame branch.

Cancel the same D^T_s from every scalar coefficient. Since H_l contains t^(3l)D^(6l), each higher coefficient still contains at least t^(3l)D^(5l). Its valuation at t=0 is at least3l, including D(0)=0. The retained nonzero zeroth coefficient stays nonzero, so there is no extra derivative-order loss.

The old linear-D row bounds are9j-g_i with g=(0,2,4,3,6). Every actual rank/row choice is included. At full rank5 the new gate is

    9*15-15-T5=120-10=110.

## 4. Constant D needs a different whole-operator normalization

A constant D gives no polynomial-degree saving from a D power, so the linear-D estimate alone would not improve the previous uniform120 gate. Treat this entire coefficient stratum separately.

Use the *unsheared* weighted Euler representative and D_prime=0. Its common-time denominator t^3D^2, source S_beta/t and both Euler fields have a common factor D. Cancelling that whole factor gives

    h_new=t^3*D,
    S_new=t*W_t-2W+(1/2)u*W_u+(3/4)v*W_v-(3/4)t^2D,
    V_new=t^2*D*((1/2)u*partial_u+(3/4)v*partial_v).

The divergence scalar is retained. This is a different representative on a justified constant-D stratum, not deletion of the shear term from a previously computed seed. Its actual seeds are recomputed and checked in the original twisted quotient. The whole old connection is D times this normalized one; because D_prime=0, the old numerators are D^j times the new ones for every j.

For the uniform degree proof retain independent D of formal degree1 but with derivative0, as in the prior specialization argument. The full phase has joint degree<=4 and pole grade rho-D exponent<=4. The new source has joint degree<=4,pole grade<=4 and old(2,3)-amplitude weight<=7. The h, h_prime and vertical derivatives respect the induction, so R_j has joint degree<=4j,pole grade<=4j and old weight<=7j.

The unchanged complete principal reductions, cleared by D^2, have joint cost<=2 and preserve the pole grade after clearing. At depth r_j=7j-4, cancel the common power

    D^(2r_j-4j)=D^(10j-8)

and retain the same completed row factors D^rho. Stronger witness grades ensure their common factor is also present; the exact identities are checked in the original rational twisted quotient. The new normal rows and scalar denominator satisfy

    deg_t N_(i,j)<=8j-|e_i|-rho(e_i)=8j-g_i,
    H_j=t^(3j)*D^(5j), valuation>=3j.

At full rank5 the cofactor degree is8*15-15=105; every smaller actual rank and row set is included. Specializing D to any nonzero constant cannot increase the formal joint bound. No degree saving is inferred just by dividing by a constant.

## 5. Reassemble every tangent frame and the complete cover

Identically zero D and required alternate charts retain parent92 and rank-one95. Hence full tangent max(110,105,92,95)=110. Its odd94 class remains separately retained from the stronger centered odd argument.

| Leading type with arbitrary lower terms | Full sufficient prefix |
|---|---:|
| Degree at most two | 7 |
| Smooth cubic | 94 |
| Irreducible node | 108 |
| Irreducible cusp | 115 |
| Transverse line-smooth-conic | 112 |
| Tangent line-smooth-conic | **110** |
| Three independent lines | 108 |
| Three distinct concurrent lines | 104 |
| Repeated leading line(s) | 99 |

The maximum115 gives the new general pure-only test, with no added Gaussian-frame hypothesis. Pure-null cubics retain GMC and the sharp general mixed cutoff m>3deg Q. Odd98 and homogeneous96 remain unchanged.

In the rational twenty-coefficient ring radical(I115)=radical(I_infinity), without ordinary ideal equality or expanded global radical cofactors. At most7 exact monotone radical-search rounds against I115 determine the least prefix in principle; no full coefficient-ring round has been run.

## Evidence, failed normalization and remaining scope

The checker expands the generic pole limits of phase, full sheared source and both gradients; computes the full rank-five Jacobi residue matrix, trace/projector identities and150 augmented minors; verifies the constant-D whole-factor/source/pole and principal-divergence identities; and checks complete original twisted seed and full scalar cofactor relations in nondegenerate Gaussian models through rank5. Models include linear/constant D,D(0)=0,q=0 and a shared q/D root. The source audit records exact counts, hashes and scalar valuations.

The [first input](tangent_rank_one_lattice_initial_input.py) and [terminal failure](tangent_rank_one_lattice_initial_failure.json) are preserved. Its pole-source assertion used the wrong1/4 scalar because the by term at D=0 had been omitted in the proposed limit. The full symbolic source exposed the missing contribution. The corrected5/16 coefficient and its complete residue matrix pass; no failed formula is promoted. This normalization correction changes the residue scalar but not the degree gates.

Universal pole/seed/minor induction, the polynomial extension across zero traces, all-degree constant-D path estimates, every coefficient/frame choice and retained analytic nonzero zeroth arguments remain written, independently unreviewed inputs. The finite exact checks do not formalize the complete analytic theorem. No fresh literature-coverage claim, external theorem or producer package is imported; preceding selected primary-record/body boundaries remain retained. Every earlier canonical proof and receipt remains byte-identical.

Least general[20,115],odd[20,98],homogeneous[14,96],20/14 sufficiency, expanded I115 radical cofactors, independent analytic review, exceptional coefficient schemes, global term minima and lower-power GVC remain open. Cusp115 and transverse112 now control general, while odd node/transverse/repeated98 and homogeneous conic96 retain separate frontiers. Qualitative cubic coverage and both lower witnesses stay closed internally.
