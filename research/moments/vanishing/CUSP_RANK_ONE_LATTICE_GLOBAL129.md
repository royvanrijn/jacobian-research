# A rank-one cusp residue lowers the general pure prefix to129

15 November 2026. Internal written refinement of the [phase-factor and stopped-budget proof](RANK_FACTORS_AND_REDUCTION_GLOBAL_130.md), retaining the [odd98 cover](EVEN_TIME_SEEDS_ODD98.md), conic/full-W-linear96, complete coefficient/frame choices and original Gaussian Wick functional.

Every complex ternary polynomial of total degree at most three has all positive pure Gaussian moments zero if and only if its first **129** pure moments vanish. Its least sufficient prefix is therefore in **[20,129]**. Odd **[20,98]**, homogeneous **[14,96]**, smooth-leading **[20,94]** and homogeneous-smooth **[14,56]** remain unchanged. Both exact first-nonzero20/14 witnesses are complete. Oddness of129 does not give a128 bound for arbitrary, non-odd polynomials.

The new cusp complement bounds are **125** when D is linear and **109** when D is a nonzero constant. Together with its parent92 and rank-one95 frame faces, the full cusp class has sufficient gate125. The unrestricted nodal129 estimate now controls the general cover.

The [checker](verify_cusp_rank_one_lattice.py), [certificate](cusp_rank_one_lattice_certificate.json), [receipt](cusp_rank_one_lattice_verification.json) and [source audit](cusp_rank_one_lattice_source_audit.json) support auxiliary exact algebra. Universal all-degree lattice, seed and minor induction, full coefficient/frame coverage, original period comparison and analytic nonzero zeroth coefficients remain written, independently unreviewed inputs. There is no independent replay, external review, formal verification or expanded global radical certificate.

## 1. Keep the full Gaussian square structure

The actual rescaled cusp phase has more structure than an arbitrary lower weighted polynomial:

    W=D*(v^3+G2)+F^2/2,
    F=u^2+p*u+q*v+r,
    G2=a*u^2+b*u*v+c*v^2+e*u+f*v+g.

Here D=m-2alpha*t is nonzero on the complement, p,q have degree at most1, r degree at most2, a,b,c degree at most1, e,f degree at most2 and g degree at most3. These statements follow directly from the retained full Gaussian completion and rescaling. Means and cross coefficients are retained. No new change of the original Gaussian functional is made.

Work first with independent D and the derivation D_prime=eps, also retaining all time derivatives of p,q,r and G2. Set delta=t^2*D, bx=delta_prime/4, by=2tD/3. The complete common-time weighted Euler source is

    Sbar=(delta*W_prime-delta_prime*W-2t^3D^2-t^4D*eps/2
                  +delta*(bx+by)+bx*u*W_u+by*v*W_v)/t.

At D=0 its exact value is

    S0=-(t*eps)*F*(p*u/4+q*v/2+r/2).

Both Sbar and S0 have amplitude half-weight rho<=2, where rho(u^i*v^k)=i/2+k; Sbar-S0 is divisible by D. This retains the ordinary-volume divergence and Gaussian determinant terms. The highest cubic v^3 source cancels, as in the previous weighted Euler representative.

## 2. The top stopping-budget reduction depends only on F^2/2

Use the old full basis1,u,u^2,v,u*v,u^2*v and stopping potential

    B(e)=max(0,ceil(rho(e)-1)),
    B_i=0,0,0,0,1,1.

Every principal step, cleared by D, has budget grade at most B(input)-1, with grade defined as B(output)-D exponent. Compare it with the base phase

    Wbase=D*v^3+F^2/2,

and omit delta-divergence only for this *associated graded calculation*. The exact full phase and delta-divergence differ from the base cleared remainder by

    -D*Cu*partial_u G2-D*Cv*partial_v G2
                                      -t^2D*(partial_u Cu+partial_v Cv).

This error has the strictly smaller grade B(input)-2. For the u^3 branch, its D-cleared witness contains D and the derivatives of G2 have rho at most1; for the v^2 branch the cleared witness has no D, its lower derivative has rho at most1, and the displayed extra D supplies the strict order. The divergence corrections are smaller as well. These estimates hold for every exponent by the two explicit branch formulas and the ceiling inequality. They are not inferred solely from the finite checked range.

Thus full reduction and base commutative reduction have the same top stopping-budget coefficient. The source error Sbar-S0, which is divisible by D, likewise loses one D order. The Euler vertical field and every denominator/basis derivative term contribute only positive D order to the residue below. They are retained in all actual seed computations; none is discarded from the connection itself.

## 3. An explicit rank-one residue, including every specialization

Let n_j be the completed normal coordinates of the actual numerator R_j in the rational basis

    E=(1,u,u^2,v,D*u*v,D*u^2*v),

and put N_j=D^(2j-1)*n_j for j>=1. Equivalently, component i is the old completed coefficient multiplied by D^(2j-1-B_i). The stopped-budget proof guarantees polynomial N_j and scalar denominator

    H_j=t^(3j)*D^(4j-1), valuation at zero at least3j.

This is the same completed row scaling already used in the130 proof. Basis derivatives must be included in an induced connection calculation. After multiplying that connection by D^2, they vanish modulo D: h=t^3D^2, so their worst h*D_prime/D contribution becomes a multiple of D^3. The actual-seed recurrence similarly shows that all h derivative and coordinate derivative terms vanish in the residue. The only surviving matrix comes from multiplication by S0 in the base quotient.

Define

    b=(q^3*r/6, p*q^3/6, q^3/6, q^4/6, -p*q/2, -q/2)^T,
    A_base=-(q^2/3)*b*e_v^T,

where e_v selects the fourth basis row. Direct full base reduction gives

    N1 mod D=(t*eps)*b,
    A0=(t*eps)*A_base,
    N_(j+1) mod D=A0*(N_j mod D), j>=1.

Every column N_j therefore lies in the same line modulo D. All315 two-by-two minors of the augmented matrix[b,A_base] vanish as polynomial identities. The formulas divide by neither q, p, r, t nor eps, so q=0, a common q/D root, D(0)=0 and rank-zero residues remain included. Constant D is handled separately for parameter degree below.

The all-j recurrence follows from the same full connection and the strict associated-graded comparison above. The finite six-column computation establishes the explicit residue operator; the written lattice/seed induction extends it to every actual seed. No sample rank observation is promoted to a universal statement.

## 4. Divide the common factor from all scalar cofactors

At actual seed rank s choose any s rows containing the unit. Polynomial cofactors of columns N0=unit,N1,...,Ns give raw minors C_l. Since all non-unit columns have rank at most1 modulo D,

    D^(s-1) divides C0.

For l>=1, the minor includes the unit and s-1 non-unit columns, so its residue rank is at most2 and

    D^max(0,s-2) divides C_l.

These are the usual determinant divisibility bounds over the D-adic field; the polynomial identities extend across every zero-parameter specialization. The scalar relation coefficients are A0=C0 and A_l=H_l*C_l for l>=1. Consequently they all contain D^(s-1). Divide the whole relation by that factor.

This division preserves the nonzero zeroth coefficient from the retained analytic exclusion. For s>=2 and l>=1, the remaining higher coefficient still contains

    t^(3l)*D^(4l-2),

so its valuation is at least3l even when D vanishes at t=0. Rank1 needs no common-factor division and retains its original scalar valuation. Thus cancellation creates no derivative-order loss.

For linear D, the old completed row degrees are7j-1-g_i with g=(0,1,2,1,3,4). The new zeroth degree is bounded by

    7*s*(s+1)/2-s-sum_selected g_i-(s-1).

Every smaller actual rank and row selection is included. At full rank6 it is147-6-11-5=125. This gives the linear-D cusp complement bound125, without assuming generic seed rank or generic nonzero q.

## 5. Nonzero constant D has a smaller direct bound

If D is constant, its degree contributes no saving when a common D factor is removed. Use the actual constant-D source instead: joint degree W<=4 and the complete normalized source has joint degree<=4. The actual unreduced numerator therefore has joint degree<=4j and amplitude rho<=2j.

The u^3 principal branch has joint cost at most0. A v^2 branch has joint cost at most1 and decreases B by at least1. At starting seed j the budget is at most2j-1, so a path ending in row i has at most2j-1-B_i such v branches. This gives

    deg_t N_(i,j)<=6j-1-|e_i|-B_i=6j-1-g_i.

Multiplication by constant D row factors changes no degree. All ranks/rows give the full-rank bound6*21-6-11=109. H_j and its valuation are unchanged. Hence109 covers this entire coefficient stratum, with arbitrary lower terms and full Gaussian means/cross coefficients.

Identically zero D is not a nonzero constant complement. Its retained parent92 and rank-one95 frame branches remain covered. The full cusp gate is therefore max(125,109,92,95)=125.

## 6. Complete general129 cover

| Leading type with arbitrary lower terms | Full sufficient prefix |
|---|---:|
| Degree at most two | 7 |
| Smooth cubic | 94 |
| Irreducible node | 129 |
| Irreducible cusp | **125** |
| Transverse line-smooth-conic | 127 |
| Tangent line-smooth-conic | 120 |
| Three independent lines | 118 |
| Three distinct concurrent lines | 104 |
| Repeated leading line(s) | 99 |

The maximum129 supplies the new general test. No extra Gaussian-frame condition is imposed. Its nullity still yields the retained full W-linear classification and sharp mixed cutoff m>3deg Q. The separate odd98 and homogeneous96 conclusions stay unchanged.

In the rational twenty-coefficient ring, radical(I129)=radical(I_infinity); this is not ordinary ideal equality or an expanded radical certificate. At most7 exact monotone radical-search rounds against I129 determine the least global prefix in principle; none has been run. The current least interval is[20,129].

## Evidence and assurance

The exact checker supports the structured source identity and strict branch filtration comparison, the explicit rank-one residue and its315 augmented minors, complete original twisted seed identities in full nondegenerate Gaussian models, full scalar cofactor relations through rank6 and the common D factors, every smaller rank/row estimate and the scalar valuations after cancellation. Models include linear D, D(0)=0, constant D, q sharing a D root and q identically zero. Retained sources and historical proof bodies are hash-bound separately.

Universal all-degree associated-graded reduction, lattice/seed/minor induction, constant-D path estimates, complete coefficient/frame choices and analytic exclusions remain written, independently unreviewed inputs. The checks do not formalize the complete analytic theorem. No fresh literature-coverage claim, external theorem or producer package is imported; the selected primary-source refresh is retained.

Least global[20,129],odd[20,98],homogeneous[14,96],20/14 sufficiency, expanded I129 radical cofactors, independent analytic review, exceptional coefficient schemes, global term minima and lower-power GVC remain open. Generic node129 and transverse127 now control the general cover; odd node/transverse/repeated98 and homogeneous conic96 retain their separate frontiers. Qualitative cubic coverage and both lower witnesses stay closed internally.
