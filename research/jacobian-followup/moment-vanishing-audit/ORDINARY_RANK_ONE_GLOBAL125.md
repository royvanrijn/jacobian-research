# Ordinary Gaussian-square residues: general sufficient prefix125

15 November 2026. Internal written refinement of the [cusp125/general129 proof](CUSP_RANK_ONE_LATTICE_GLOBAL129.md), retaining its complete coefficient/frame cover, original Wick functional and analytic exclusions.

Every complex ternary polynomial P of total degree at most three has all positive pure Gaussian moments zero if and only if its first **125** pure moments vanish. Its least sufficient prefix is therefore in **[20,125]**. Odd **[20,98]**, homogeneous **[14,96]**, smooth-leading **[20,94]** and homogeneous-smooth **[14,56]** remain unchanged. Conic/full-W-linear96 and both exact first-nonzero20/14 witnesses are retained. Odd125 does not imply124 for arbitrary non-odd P.

The three ordinary complement bounds are now node123,transverse122,triangle114 for linear D, and101,106,103 for nonzero constant D. The previously established cusp125 now controls the general cover.

The [checker](verify_ordinary_rank_one_lattice.py), [certificate](ordinary_rank_one_lattice_certificate.json), [receipt](ordinary_rank_one_lattice_verification.json) and [source audit](ordinary_rank_one_lattice_source_audit.json) support auxiliary exact algebra. Universal all-order polynomial-module/seed/minor induction, full coefficient/frame coverage, original period comparison and nonzero zeroth coefficients remain written, independently unreviewed inputs. There is no independent replay, external review, formal verification or expanded global radical certificate.

## 1. Retain the actual translated Gaussian square

The [universal translated phase](ordinary_node_phase_shift_checkpoint.json) has

    W=D*G+F^2/2, F=u*v+r(t), D=m-2alpha*t,
    G=a*u^3+d*v^3+Le*u^2+Lf*v^2+J0*u*v+L0*u+N0*v+K0,
    delta=t^2*D, h=t^3*D^2.

The coefficient degrees of G are at most3 minus amplitude degree, and deg r<=2. Nodal a,d are both nonzero; transverse d=0,a!=0; triangle a=d=0. The required Le,Lf axes on the latter types remain nonzero polynomials on their complement charts. Full affine Gaussian means and cross coefficients remain present in the checkpoint. No centered odd restriction is imposed.

Work first with independent D, D_prime=eps, independent polynomial coefficients and their derivatives. Retain the full ordinary-volume radial Euler source

    Sbar=t*D*W_prime-(2D+t*eps)*(W-E(W)/4)-t^2D^2,
    E=u*partial_u+v*partial_v.

Since E(F)=2uv, its exact value at D=0 is

    S0=-(t*eps*r/2)*F.

The source error Sbar-S0 is divisible by D. Keep the full denominator derivative and vertical field

    V=delta*delta_prime*(u*partial_u+v*partial_v)/(4t)

in every actual seed computation. V has a factor D; its coefficients vanish in the residue after forming the connection, rather than being omitted from the operator.

## 2. An exact identity gives the rank-one residue

The [polynomial normal module](NORMAL_MODULE_GLOBAL_146.md) has normal coefficients polynomial in D: for node no axis denominator is needed; transverse and triangle allow only the retained Lf or Le/Lf inverses before bounded clearing. Exact witnesses may have rational D coefficients, and are retained in the original twisted quotient.

Let N be this normal projection and f be any nonconstant amplitude monomial. If u divides f, put C=f/u. The full twisted derivative gives the exact identity

    F*f=(delta*partial_v+W_v)(C)
                                      -D*(G_v*C+t^2*partial_v C).

If only v divides f, use C=f/v and the u derivative instead. Thus

    N(F*f)=-D*N(G_v*C+t^2*partial_v C)

or its u version. Polynomiality of the retained normal module proves divisibility by D. This is an all-exponent consequence of one explicit operator identity, not a numerical rank observation. It includes the full delta-divergence term.

Every complete retained basis contains the unit and uv. Therefore N(F)=r*unit+uv exactly, while multiplication by F kills every non-unit basis element modulo D. The residue of multiplication by Sbar is

    A0=b*e_unit^T,
    b=-(t*eps*r/2)*(r*unit+uv).

All two-by-two minors of the augmented matrix[b,A0] vanish identically. Neither r nor t*eps is inverted, so zero residues remain included. No extra coefficient open is introduced.

## 3. Extend the residue to every actual seed and all scalar minors

Use the same completed row lattices as the previous normal-module and stopping-budget proofs. Write n_j for normal coordinates of the actual numerator R_j in that lattice, before its scalar column clearing. Any row rescaling by powers of Le,Lf must be differentiated in the induced connection. Those basis-derivative contributions contain h=t^3D^2 and vanish modulo D over the formal axis-localized field. The h derivative has a factor D and the Euler vertical field has a factor D. Consequently the actual recurrence gives

    n1 mod D=b, n_(j+1) mod D=A0*(n_j mod D).

For the node lattice this is literal in its fixed monomial basis. On the axis-scaled transverse/triangle lattices, b is expressed in those same coordinates; the rank-one conclusion is unchanged. Column j is subsequently multiplied by its already bounded axis clearing polynomial. Such column scalars preserve residue rank at most1.

Thus all completed polynomial non-unit seed columns N_j lie in one line modulo D. The formal Le,Lf inverses are cleared before specialization; the resulting polynomial D divisibility extends to shared axis/D/t roots. It is not inferred by substituting a shared root into an uncleared rational basis. Identically zero required axes remain on their retained frame cases.

At actual seed rank s select any rows containing the unit. The zeroth raw cofactor C0 has D^(s-1); a higher raw cofactor C_l, which includes the unit column, has D^max(0,s-2). The scalar relation coefficients are C0 and H_l*C_l. All of them therefore have a common factor D^(s-1). Divide that same factor from the entire scalar relation.

For s>=2 each remaining higher coefficient still contains at least

    t^(3l)*D^(2l-1)

times its nonnegative axis-clearing factor. It has valuation at least3l, including D(0)=0. Rank1 requires no common-factor cancellation and keeps its original valuation. The retained analytic nonzero zeroth coefficient remains nonzero after division. Therefore the same initial-order contradiction applies without additional derivative-order loss.

## 4. Linear-D degree bounds

Keep the previous complete polynomial row estimates and subtract s-1 from the zeroth determinant degree:

| Ordinary complement | Retained row bound | New linear-D sufficient gate |
|---|---|---:|
| Node | 5j minus ordinary row degree | **123** |
| Transverse | 7j-2 minus ordinary row degree | **122** |
| Triangle | 9j minus combined row grade | **114** |

The nodal rank7 ordinary row degrees sum11, so5*28-11-6=123. Transverse rank6 row degrees sum8, giving7*21-12-8-5=122. Triangle rank5 combined row grades sum17, giving9*15-17-4=114. Every smaller actual seed rank and every selected row set containing the unit is included in the arithmetic.

The triangle combined grade is i+k+ceil(i+k/2)+ceil(i/2+k); its old Le/Lf completed row scales and column clearing are unchanged. Transverse keeps Lf^(2j-2), and node has no column-axis denominator. No D factor is removed from a single coefficient in isolation.

## 5. Nonzero constant D is a separate degree stratum

A constant divisor gives no polynomial-degree saving upon cofactor cancellation. Instead the complete source has joint degree at most4, rather than the safe5 bound used for linear D. Indeed W has joint degree at most4 and the constant-D radial source preserves that degree. The actual numerator recurrence gives joint degree R_j<=4j, with all derivatives retained.

The same normal module and axis filtrations apply. Thus the row slopes decrease by1, from5/7/9 to4/6/8; offsets and row grades stay unchanged. Full-rank sufficient complement gates are

    node:4*28-11=101,
    transverse:6*21-12-8=106,
    triangle:8*15-17=103.

All smaller ranks/rows are included. These arbitrary-lower-term constant-D estimates are separate from the stronger centered odd bounds; no oddness is assumed here. All scalar denominators and valuations remain those of the retained actual seed computations.

## 6. Add every frame face and the complete general cover

Conic96 is below every displayed maximal ordinary complement gate. Hence the full ordinary classes have bounds123,122,114. Retain full cusp125, tangent120, concurrent104, repeated99, smooth94 and quadratic7.

| Leading type with arbitrary lower terms | Full sufficient prefix |
|---|---:|
| Degree at most two | 7 |
| Smooth cubic | 94 |
| Irreducible node | **123** |
| Irreducible cusp | 125 |
| Transverse line-smooth-conic | **122** |
| Tangent line-smooth-conic | 120 |
| Three independent lines | **114** |
| Three distinct concurrent lines | 104 |
| Repeated leading line(s) | 99 |

The maximum125 proves the new general pure-only test, with no added Gaussian-frame hypothesis. Pure nullity retains the complete cubic GMC conclusion and sharp uniform mixed cutoff m>3deg Q. Odd98 and homogeneous96 stay unchanged.

In the rational twenty-coefficient ring radical(I125)=radical(I_infinity), without ordinary ideal equality or expanded radical cofactors. At most7 exact monotone radical-search rounds against I125 determine least global K in principle; no full coefficient-ring round has been run.

## Evidence and assurance

The checker verifies132 universal F-monomial operator identities, three exact source residues and1053 rank-one augmented minors. Full Gaussian models check complete original twisted seed identities through every class rank and full scalar cofactor relations, common D factors, constant-D source/normal degree estimates and scalar valuations after cancellation. Models retain nonzero metric determinant, linear D, nonzero constant D, D(0)=0 and shared Le/Lf/D/t roots. The source audit records exact receipt counts and hashes.

Universal all-order polynomial-module/seed/minor induction, every normal/frame coefficient choice, original period comparison and analytic nonzero zeroth arguments remain written, independently unreviewed inputs. Exact auxiliary checks do not formalize the full analytic theorem. A selected primary-record refresh rechecks [Long2607.18186](https://arxiv.org/abs/2607.18186), [Dvorsky2608.07338](https://arxiv.org/abs/2608.07338) and [Wang--Zhang2610.16157](https://arxiv.org/abs/2610.16157): all still list v1. This record/abstract check does not extend the previous selected body-read coverage or establish exhaustive literature coverage. No new external theorem or producer package is imported; the earlier primary-body audit remains retained. Every earlier canonical proof body and receipt stays byte-identical.

Least global[20,125],odd[20,98],homogeneous[14,96],20/14 sufficiency, expanded I125 radical cofactors, independent analytic review, exceptional coefficient schemes, global term minima and lower-power GVC remain open. Cusp125 and node123 now control the general cover; odd node/transverse/repeated98 and homogeneous conic96 keep separate frontiers. Both lower witnesses and qualitative cubic coverage remain closed internally.
