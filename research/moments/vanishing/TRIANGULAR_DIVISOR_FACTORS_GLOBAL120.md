# Triangular divisor factors in scalar cofactors: general prefix120

15 November 2026. Internal written refinement of the [ordinary/general125 proof](ORDINARY_RANK_ONE_GLOBAL125.md) and [cusp residue proof](CUSP_RANK_ONE_LATTICE_GLOBAL129.md), retaining the original Wick functional, complete coefficient/frame cover and analytic exclusions.

Every complex ternary polynomial of total degree at most three has all positive pure Gaussian moments zero if and only if its first **120** pure moments vanish. The least general prefix is therefore in **[20,120]**. Odd **[20,98]**, homogeneous **[14,96]**, smooth-leading **[20,94]** and homogeneous-smooth **[14,56]** remain unchanged, as do conic/full-W-linear96 and both exact first-nonzero20/14 witnesses. Even120 gives no119 parity deduction.

The new linear-D complement bounds are node108,transverse112,triangle108,cusp115. Their nonzero constant-D gates remain101,106,103,109. All frame gates remain smaller, so these are also the full class bounds. Tangent120 now controls the general cover.

The [checker](verify_triangular_divisor_factors.py), [certificate](triangular_divisor_factors_certificate.json), [receipt](triangular_divisor_factors_verification.json) and [source audit](triangular_divisor_factors_source_audit.json) bind exact auxiliary evidence. The general formal block/seed/minor lemma below, polynomial extension across degenerate residue specializations, all coefficient/frame choices and retained analytic exclusions remain written, independently unreviewed inputs. No independent replay, external review, formal verification or expanded global radical certificate is asserted.

## 1. A formal differential lemma gives more than first-order rank counting

Let K be a characteristic-zero differential coefficient field, with a formal divisor D and derivation delta satisfying ord_D(delta f)>=ord_D f-1. Suppose analytic column vectors X_j in K[[D]]^r obey, for j>=1,

    X_(j+1)=a*delta X_j+(A+c_j*Id)*X_j,
    ord_D a>=2, A in Mat_r(K[[D]]), ord_D c_j>=1.

Assume A0=A mod D has rank1 and nonzero trace lambda0, and X1 mod D lies in its image. X0 may be any analytic vector. Put T_s=s*(s-1)/2.

Then every selected s-by-s minor of X1,...,Xs has D order at least T_s. For the minors omitting column l from X0,...,Xs, l>=1, a safe order bound is

    max(0,T_s-l).

These estimates require no assertion that the selected determinant itself is nonzero. They apply at every actual seed rank and every row set, including those containing the unit used by the scalar relation.

To prove them, rank1 and trace lambda0!=0 give A0^2=lambda0*A0, hence a direct sum of its one-dimensional image and zero kernel. In that decomposition a formal invertible gauge P block-diagonalizes the differential operator a*delta+A into one fast scalar block and an (r-1)-dimensional block B in D*Mat(K[[D]]). Construct P order by order: at order D^n the off-diagonal commutator with A0 is invertible with eigenvalues plus/minus lambda0. The derivative correction a*delta(U*D^n) has order at least n+1, so it does not obstruct solving the current order. The coefficient-field derivatives are retained in this construction. The scalar c_j*Id is unchanged by this single gauge.

Write Y_j=P^(-1)X_j. Its fast coordinate stays analytic. Its slow coordinates start in D and have order at least j by induction: B and c_j increase order by at least1, and a*delta also increases order by at least1. At most one row in a transformed determinant can use the fast coordinate. Put it on the largest column index; the other indices1,...,s-1 contribute1+...+(s-1)=T_s. Cauchy--Binet and the analytic invertible gauge preserve this lower bound for every original selected minor.

For a minor containing X0 and omitting X_l, put the fast coordinate on the largest remaining index. If l<s the slow sum is T_s-l. If l=s it is T_(s-1), which is even stronger than T_s-s. Negative bounds are replaced by0. This proves the lemma with the original derivation fully accounted for.

## 2. Verify its assumptions on the existing four residue lattices

On the node/transverse/triangle normal modules, before their bounded axis column clearing, retain h=t^3D^2 and the complete induced recurrence

    n_(j+1)=h*n_j_prime+(M-j*h_prime*Id)*n_j.

The previous proof includes every basis derivative and gives M analytic in D over the formally axis-localized coefficient field. Its residue is

    M0=b*e_unit^T,
    b=-(t*D_prime*r/2)*(r*unit+uv),
    trace M0=-(t*D_prime*r^2)/2.

The trace is not identically zero in the full coefficient family. h has D order2, h_prime has order at least1 and n1 mod D=b. The lemma applies on that generic coefficient field. Multiplying column j by its old axis-clearing polynomial does not reduce its formal D order.

For cusp use the already polynomial rescaled columns N_j=D^(2j-1)*n_j in the completed basis(1,u,u^2,v,Duv,Du^2v). Its derivative coefficient is D^2*h=t^3D^4. The additional j-dependent scalar terms have D order at least3. The same source/potential estimate used to form the existing residue shows the full induced matrix analytic in D, not merely a finite residue formula. The retained residue is

    A0=-(t*D_prime*q^2/3)*b_base*e_v^T,
    trace A0=-(t*D_prime*q^6)/18,
    N1 mod D=(t*D_prime)*b_base.

This trace too is not identically zero. Thus its derivative order4 and scalar order3 meet the lemma's weaker2/1 hypotheses, with N1 in the residue image.

Formal diagonalization may invert the generic trace during its proof. It is not applied directly to a zero-trace member. The resulting minors, after the already bounded axis clearing, are polynomial in D and the original coefficients. Their low D coefficients vanish on the generic coefficient field and therefore vanish identically. Polynomial divisibility consequently extends to r=0,q=0,zero traces at D(0)=0 and shared axis/D/t roots. No trace, q, r or t factor is cancelled from the actual scalar relation. Identically zero required axes and D remain on their separately retained frame charts.

This polynomial extension is essential; a generic formal splitting alone would not justify dropping a degenerate coefficient stratum.

## 3. Cancel a common triangular power without losing scalar valuations

At actual seed rank s choose rows containing the unit and form raw cofactors C_l of columns N0,...,Ns. The new lemma gives

    D^T_s divides C0,
    D^max(0,T_s-l) divides C_l, l>=1.

The scalar coefficients are A0=C0 and A_l=H_l*C_l. Ordinary H_l contains t^(3l)D^(2l) and its nonnegative axis-clearing factors. Cusp H_l=t^(3l)D^(4l-1). Both provide at least D^l, so every scalar coefficient is divisible by the same D^T_s. Divide that common factor from the entire relation.

After cancellation, each ordinary higher coefficient still contains at least t^(3l)D^l times its nonnegative axis factor; cusp retains at least t^(3l)D^(3l-1). These inequalities also hold when T_s<l, since the raw cofactor is polynomial. Hence every higher coefficient retains valuation at least3l at t=0, including D(0)=0. The original nonzero zeroth coefficient stays nonzero. Rank1 has T1=0 and requires no cancellation.

No ordinary or Gaussian functional is changed. The initial-order contradiction remains valid with the same denominator-derivative and period assumptions.

## 4. Every actual rank and row set, not only the maximal generic rank

Subtract T_s from the old raw cofactor degree at each actual rank s:

| Complement | Old raw full-rank degree | T_s at full rank | New linear-D gate |
|---|---:|---:|---:|
| Node,rank7 | 129 | 21 | **108** |
| Transverse,rank6 | 127 | 15 | **112** |
| Triangle,rank5 | 118 | 10 | **108** |
| Cusp,rank6 | 130 | 15 | **115** |

The checker includes every smaller rank and selected row set. Their maxima equal the full-rank displayed bounds. This is a sufficient degree calculation, not a computation of the least prefix or of expanded global radical cofactors.

Nonzero constant D has no degree saving from a D power. Retain node101,transverse106,triangle103,cusp109 from the separately proved constant-D source/path estimates. Conic96,parent92,rank-one95 and all zero-axis/frame cases remain below the respective full class maxima. None is discarded by the generic trace argument.

## 5. Complete general120 cover

| Leading type with arbitrary lower terms | Full sufficient prefix |
|---|---:|
| Degree at most two | 7 |
| Smooth cubic | 94 |
| Irreducible node | **108** |
| Irreducible cusp | **115** |
| Transverse line-smooth-conic | **112** |
| Tangent line-smooth-conic | 120 |
| Three independent lines | **108** |
| Three distinct concurrent lines | 104 |
| Repeated leading line(s) | 99 |

The maximum120 proves the general pure-only test, with no added Gaussian-frame hypothesis. Pure-null cubics retain the complete GMC conclusion and sharp mixed cutoff m>3deg Q. The odd98 and homogeneous96 covers and their sharp mixed cutoffs remain separate and unchanged.

In the rational twenty-coefficient ring radical(I120)=radical(I_infinity), without ordinary ideal equality or an expanded radical certificate. At most7 exact monotone radical-search rounds against I120 would determine least general K in principle; no full coefficient-ring round has been run.

## Evidence and failure boundary

The checker verifies the four exact residue traces and projector identities, six abstract differential recurrences of ranks3/5/7 with derivative orders2/4, and recomputes all full cofactor relations of the ten retained linear-D Gaussian models with the larger D^T_s cancellation. Every scalar valuation after cancellation and every rank/row degree bound is checked. Previous original twisted seed certificates, constant-D checks and all prior receipts remain retained; their mains are not rerun.

The [first abstract input](triangular_divisor_factors_initial_input.py) and [terminal receipt](triangular_divisor_factors_initial_failure.json) are preserved. That synthetic family produced a zero full-rank zeroth cofactor, so it did not supply the intended full cyclic-rank example; the divisibility checks reached before that assertion had passed. The accepted input uses distinct diagonal first-order coefficients and explicitly checks every full determinant is nonzero. No failed run is promoted, and no counterexample to the general formal lemma is discarded.

The general formal splitting/seed/minor argument and its polynomial extension are written inputs, not certified by the finite examples. Complete normal/frame coverage, period comparison and analytic exclusions remain independently unreviewed written inputs. No fresh literature-coverage claim, new external theorem or producer package is imported; the preceding selected primary-record/body audit boundaries remain retained. All earlier canonical proofs and receipts stay byte-identical.

Least general[20,120],odd[20,98],homogeneous[14,96],20/14 sufficiency, expanded I120 radical cofactors, independent analytic review, exceptional coefficient schemes, global term minima and lower-power GVC remain open. Tangent120 and cusp115 now control general, while odd node/transverse/repeated98 and homogeneous conic96 retain separate frontiers. Qualitative cubic coverage and both lower witnesses stay closed internally.
