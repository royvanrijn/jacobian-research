# Phase factors and stopping pole budgets: global130

15 November 2026. Internal written refinement of the [two-sided conic140 proof](TWO_SIDED_CONIC_GLOBAL_140.md), retaining the original Wick functional, full quotient bases, complete coefficient/frame cover and analytic exclusions.

Every complex ternary polynomial P of total degree at most three satisfies

    E(P^j)=0 for j=1,...,130 iff all positive pure moments vanish.

There is no extra Gaussian-frame hypothesis. Pure nullity gives the W-linear frame and sharp mixed cutoff m>3*deg Q. Least global and odd-only prefixes lie in **[20,130]**, while the homogeneous interval stays **[14,130]**. Smooth-leading **[20,94]** and homogeneous-smooth **[14,56]** are unchanged. The first-nonzero20/14 witnesses remain complete. Even130 gives no129 parity reduction. Exact leastness and20/14 sufficiency are open.

The [checker](verify_rank_factors_and_reduction_budgets.py), [certificate](rank_factors_and_reduction_budgets_certificate.json), [receipt](rank_factors_and_reduction_budgets_verification.json) and [source audit](rank_factors_and_reduction_budgets_source_audit.json) support exact auxiliary algebra. All-degree inductions, full coefficient/rank coverage, period comparison and analytic nonzero zeroth coefficients remain written inputs. This is not independent, external or formal verification, nor an expanded global radical certificate.

## 1. Rank-one phase degree forces common q powers

For the complete rank-one family P=W*(a+bZ+cT+fT^2)+G(Z,T), f!=0, retain q=1-bt, y=tT and ahat=a*t^2+c*t*y+f*y^2. Its normalized period has prefactor1/(tq) and phase W0/delta0, where

    W0=sum_(i+j<=3) g_ij*t^(3-i-j)*q^(3-i)*ahat^i*y^j-q^3*y^2/2,
    delta0=t^2*q^3.

Here W0 is a phase numerator, not the circular Gaussian coordinate. Its actual y degree d lies between2 and6. The highest coefficients, successively setting higher ones to zero, are

    y^6: f^3*g30,
    y^5: f^2*q*g21,
    y^4: f*q*(f*t*g20+q*g12),
    y^3: q^2*(f*t*g11+q*g03).

Since q(0)=1 and f!=0, vanishing of y^4 forces g12=0 and then g20=0; vanishing of y^3 similarly forces g03=g11=0. These implications include b=0. At degree2 the remaining leading coefficient divided by q^2 is f*g10*t^2+g02*t*q-q/2, nonzero at t0. Therefore W0 contains q^m with

| d | 2 | 3 | 4 | 5 | 6 |
|---|---:|---:|---:|---:|---:|
| m | 2 | 2 | 1 | 1 | 0 |

Cancel q^m simultaneously from phase numerator and delta: Wm=W0/q^m, deltam=t^2*q^(3-m). This changes neither the phase nor the prefactor. Establish graded bounds first with independent q of degree1 and compatible derivative q_prime=-b. Then deg_joint Wm<=6-m and its leading derivative coefficient has degree<=6-m-d. Actual q=1-bt, including constant q1, is a specialization that cannot increase these bounds.

## 2. Cancel the common factor of the whole rank-one connection

For n=3-m, the original normalized source written with Wm,deltam is

    Sm=deltam*(Wm)_t-deltam_prime*Wm-deltam^2/t+b*deltam^2/q.

It has the polynomial factor t*q^(n-1). Cancelling it from the entire connection gives

    h=t^3*q^(n+1),
    Sbar=t*q*(Wm)_t-(2q-n*b*t)*Wm
                                      -t^2*q^(n+1)+b*t^3*q^n.

Both h and Sbar have joint degree<=7-m. In terms of the original source S0 and h0=delta0^2, their common factor is F=t*q^(m+2), so h0=F*h,S0=F*Sbar. The actual numerator recurrence is

    R_(j+1)=h*(R_j)_t-j*h_prime*R_j+Sbar*R_j.

The old numerator is F^j*R_j for every j. The derivative of F cancels exactly against the derivative of F*h in this identity; h_prime cannot be omitted. Thus the actual moment seed is unchanged, with joint degree R_j<=(7-m)j and amplitude degree<=dj.

The reduced leading coefficient is Anew=Aold/q^m. At the retained depth r=d*(j-1)+2, divide completed normal numerators and H by F^j*q^(m*r). Exact witnesses are divided by the same factor while retaining the original twisted operator. They can have rational time coefficients, which are explicitly checked in the original quotient over C(t). Normal coefficients are polynomial and satisfy

    deg N_(i,j)<=(7-m)j+(6-m-d)*r-i,
    H_j=h^j*Anew^r, valuation(H_j)>=3j.

All actual seed ranks and row subsets containing1 give the following sufficient rank-one gates:

| d | 2 | 3 | 4 | 5 | 6 |
|---|---:|---:|---:|---:|---:|
| Gate | 9 | 21 | 51 | 54 | 95 |

Consequently **rank-one95** replaces140. The analytic classification outside b=0,G=0 is unchanged; its nonzero zeroth coefficient is not proved by these polynomial checks. The pure-null b=0,G=0 boundary remains included, without applying that exclusion to it.

## 3. A pole budget can stop at the basis

The earlier pole grade rho(e)-divisor exponent is preserved by actual rational reductions. A sharper nonnegative potential is available when a simple divisor pole appears only on sufficiently large monomials. Define

    B(e)=max(0,ceil(rho(e)-tau)).

Suppose each individual correction has pole order at most1, a pole can occur only when rho(input)>=tau+1, and the old pole grade does not increase. A term with no pole cannot increase B after accounting for any positive divisor powers. For a term with one pole, rho(output)<=rho(input)-1, and hence B(output)+1<=B(input). This proves the stopping estimate for every explicit reduction branch, not just sampled exponents. Coefficient powers propagate it by linearity. B need not be a multiplicative filtration: actual seed bounds are imported from the already proved linear pole grade before this reduction potential is applied.

For the transverse phase use rho=i/3+k, divisor Lf and tau2. In the polynomial-normal module from the preceding proof, only the pure-v branch v^n,n>=3 introduces 1/Lf. The other branches have no Lf pole. Thus the hypothesis holds and all six retained basis monomials have B0. The earlier seed pole bound2j, together with polynomial unreduced seed coefficients, gives starting budget<=2j-2. Complete normals therefore need at most Lf^(2j-2), rather than Lf^(2j).

At the old depth r=3j-2, cancel D^r*Lf^j from the old normal numerators and H, retaining exact rational witnesses. No row rescaling is needed. The new bounds are

    H_j=t^(3j)*D^(2j)*Lf^(2j-2),
    deg N_(i,k;j)<=7j-2-i-k.

The full rank6 zeroth degree is7*21-2*6-8=127; every smaller seed rank/row choice is included. This is **transverse complement127**. Formal independent D,Lf factors are combined before shared-root or constant specialization, so no coprimality open is added.

## 4. The same stopping estimate gives cusp130

For cusp use rho=i/2+k, divisor D and tau1. Only the principal v^2 branch introduces 1/D; its input has rho>=2. The principal u^3 branch has no D pole. The half-weight pole estimate from the preceding proof thus supplies the stopping inequality. The retained seed pole bound2j gives starting budget<=2j-1.

The basis1,u,u^2,v,u*v,u^2*v has B values0,0,0,0,1,1. At depth r=5j-3, the old normal coefficient contains D^(3j-2+B_i). Cancel common D^(3j-2), and scale completed rows by D^B_i. Witnesses remain explicitly checked over the original rational quotient. Then

    H_j=t^(3j)*D^(4j-1),
    deg N_(i,k;j)<=7j-1-i-k-B_i.

The combined row grades are0,1,2,1,3,4, summing11. At full rank6 the bound is7*21-6-11=130, with all smaller ranks/row subsets included. This is **cusp complement130**. It is a sufficient gate, not a least-threshold computation. As before, degree savings are formal before D becomes a constant or has a zero at t0. Row scaling occurs after the actual seeds have been computed, with no omitted basis derivative terms.

## 5. Global130 and remaining assurance

The retained conic130 and mixed-parent92, together with rank-one95, give full W-linear130. The full cusp gate is max(130,92,95)=130. The tangent complement120 now dominates its frame faces92/95, so full tangent120 follows. The complete cover is

| Leading type with arbitrary lower terms | Full sufficient gate/exclusion |
|---|---:|
| Degree at most two | 7 |
| Smooth cubic | 94 |
| Irreducible node | max(129,130)=130 |
| Irreducible cusp | 130 |
| Transverse line-smooth-conic | max(127,130)=130 |
| Tangent line-smooth-conic | 120 |
| Three independent lines | max(118,130)=130 |
| Three distinct concurrent lines | 104 |
| Repeated leading line(s) | 99 |

The maximum130 proves the global pure-only test. The all-homogeneous gate remains130; even130 gives no129 parity deduction in either odd or homogeneous scope. In the rational twenty-coefficient ring radical(I130)=radical(I_infinity), without ordinary ideal equality or expanded global cofactors. At most ceil(log2(130-19))=7 exact monotone radical-search rounds against I130 would find the least global K; none has run.

The checker proves the five rank-one hierarchy/source/factor identities, symbolic-order recurrence normalization, exact33 seed identities in seven rank-one and four transverse/cusp models, finite principal/divergence budget identities supporting the explicit all-degree branch lemma, and all rank/row cofactor arithmetic. Every model retains the original rational twisted operator. Models include b0/1, actual degrees2 through6, the flattened and pure-null boundaries, D(0)=0 and vanishing axis constants. All earlier proof bodies and receipts remain unchanged.

No external theorem, producer package or fresh literature-coverage claim is imported by this algebraic refinement; the previous selected primary-source audit remains in force. Universal all-degree and full coefficient/rank coverage and analytic period exclusions remain written, independently unreviewed inputs. Exact leastness,20/14 sufficiency, expanded I130 radical cofactors, independent analytic review, exceptional schemes, global term minima and lower-power GVC remain open. Qualitative cubic coverage and both lower witnesses stay closed internally.
