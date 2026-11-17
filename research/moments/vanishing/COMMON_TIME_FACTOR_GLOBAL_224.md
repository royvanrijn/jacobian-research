# Cancel the common time factor: global224 and ordinary complements199/199/189

15 November 2026. Internal written refinement of the
[structured239 proof](STRUCTURED_SEEDS_GLOBAL_239.md), with its original
Gaussian Wick functional, full quotient bases and analytic exclusions retained.

Every complex ternary polynomial P of total degree at most three satisfies

    E(P^j)=0 for j=1,...,224
      iff every positive pure moment vanishes.

No extra Gaussian-frame hypothesis is imposed. Pure nullity then gives the
retained W-linear frame and sharp mixed cutoff m>3*deg Q. Least global and
odd-only prefixes lie in **[20,224]**. The all-homogeneous interval remains
**[14,216]**, smooth-leading **[20,94]**, homogeneous-smooth **[14,56]**.
Both isolated algebraic lower witnesses with first nonzero moments20 and14
remain complete. No leastness or20/14 sufficiency is proved.

The [checker](verify_common_time_seed_factor.py),
[certificate](common_time_seed_factor_certificate.json),
[receipt](common_time_seed_factor_verification.json) and
[source audit](common_time_seed_factor_source_audit.json) bind exact auxiliary
algebra. Full class/rank coverage, all-degree degree inductions and analytic
nonzero zeroth coefficients remain written mathematics. Independent replay,
external review, formal verification and expanded radical cofactors remain open.

## 1. Remove a factor from the whole cleared operator

The retained binary connections have the form

    nabla(f)=f_t+(V(f)+S*f)/h,

where V differentiates only the amplitude variables. For ordinary-node and
cusp phases h=delta^2, delta=t^2*D; for the tangent shear h=delta^2*D.
Every coefficient of V and S is divisible by t. Explicitly, the radial
ordinary-node source divided by t is

    S_E/t=t*D*W_t-(2D+t*D_prime)*(W-E(W)/4)-t^2*D^2,
    E=x*partial_x+y*partial_y.

It is polynomial, has joint t/x/y degree at most5 and amplitude degree at
most3. Its scalar Euler coefficient delta*delta_prime/(4t) is polynomial
of degree at most4. The ordinary cleared denominator becomes t^3*D^2.

For cusp, the weighted Euler source has factor t because delta has factor
t^2 and every beta coefficient has factor t. Dividing source and vector
coefficients by t gives joint source degree5 and unchanged weight bound10.
The cleared denominator again becomes t^3*D^2.

For the tangent shear, its source is D*S_beta-C*x*W_y/2, where
C=[x^2*y]S_beta. Since t divides S_beta it also divides C. Dividing the
whole operator by t gives joint source degree6, unchanged weight bound6
and cleared denominator t^3*D^3. Every Euler and shear derivative coefficient
is divided as well. There is no change to either the original volume
normalization or the quotient connection.

Put hbar=h/t, Vbar=V/t, Sbar=S/t. The actual seeds can be written as
Rbar_j/hbar^j, with the full recurrence

    Rbar_0=1,
    Rbar_(j+1)=hbar*(Rbar_j)_t-j*hbar_prime*Rbar_j
                                     +Vbar(Rbar_j)+Sbar*Rbar_j. (1)

For the old recurrence, R_j=t^j*Rbar_j at every j. Indeed differentiating
t^j*Rbar_j produces the extra j*h/t^2 term in the normalized recurrence,
which is exactly part of -j*(h/t)'. Dropping that derivative contribution
would be wrong. The checker verifies this identity with symbolic j and
independently checks each replayed seed against the old recurrence.

The complete polynomial reductions are C[t]-linear: vertical derivatives
do not differentiate t. Their already-cleared numerators and divergence
witnesses therefore also contain t^j. Cancelling t^j from the numerator,
denominator and witness is exact. This does not require any principal
divisor to be coprime to t and includes D(0)=0. It makes no new parameter
open or generic-rank assumption.

## 2. Retain the reductions and lower every column degree by j

The principal divisors, amplitude supports and reduction depths stay as in
the preceding proofs. The normalized seed joint degrees decrease from6j
to5j for ordinary-node/cusp cases, and from7j to6j for tangent shear.
Thus the new row-degree bounds and complement exclusions are

| Type | Row coefficient t-degree bound | Full-rank zeroth degree/gate |
|---|---|---:|
| Irreducible node | 8j-2-d_i | 8*28-2*7-11 = **199** |
| Transverse line-conic | 11j-4-d_i | 11*21-4*6-8 = **199** |
| Three independent lines | 14j-3-d_i | 14*15-3*5-6 = **189** |
| Irreducible cusp | 10j-3-d_i | 10*21-3*6-9 = **183** |
| Tangent line-conic | 18j-8-d_i | 18*15-8*5-6 = **224** |

The d_i are ordinary basis grades. All smaller actual seed ranks and row
subsets containing the unit row are included; their bounds are smaller.
Every new scalar higher coefficient contains Hbar_j=H_j/t^j, whose valuation
is still at least3j. The original initial-order contradiction therefore
retains no extra derivative-order loss. Nonzero zeroth coefficients follow
from the same retained analytic exclusions, since every derivative seed is
unchanged as a rational quotient class.

In the homogeneous nonisotropic-node chart, W is t-independent and D is
constant, with zero Gaussian means. Its divided source has joint degree3,
and constant principal divisors preserve joint degree. The resulting
cofactor degree is3*28-11=73. A homogeneous cubic is odd, so its moment73
vanishes automatically; the sufficient homogeneous complement prefix is
therefore **72**. The checker records the cofactor-derived sufficient73;
the written oddness deduction supplies72. The full homogeneous gate216
continues to be controlled by the conic branch. Even216 is not reduced to215.

## 3. Add back every frame branch

The retained conic216, mixed-parent210 and rank-one140 gates are unchanged.
For ordinary nodes the conic frame branch dominates each new complement,
giving full-class216. The cusp complement183 is below its mixed-parent210
frame branch. The tangent complement224 is above both210 and140.
The complete cover is

| Leading type, arbitrary lower terms | Full-class sufficient gate/exclusion |
|---|---:|
| Degree at most two | 7 |
| Smooth cubic | 94 |
| Irreducible nodal cubic | max(199,216)=216 |
| Irreducible cuspidal cubic | max(183,210,140)=210 |
| Transverse line-smooth-conic | max(199,216)=216 |
| Tangent line-smooth-conic | max(224,210,140)=224 |
| Three independent lines | max(189,216)=216 |
| Three distinct concurrent lines | 104 |
| Repeated leading line(s) | 99 |

Its maximum224 proves the global pure-only test. Parity does not reduce an
even224 gate to223, so the current odd interval is[20,224]. In the rational
twenty-coefficient ring, radical(I224)=radical(I_infinity), without ordinary
ideal equality or expanded global radical cofactors. At most
ceil(log2(224-19))=8 exact monotone search rounds against I224 would find
the least global K. No such full coefficient-ring search has run.

## Evidence and remaining work

The checker verifies three universal polynomial source/operator factorizations,
the recurrence cancellation at symbolic j and the homogeneous source degree3.
It replays68 complete connection-column and36 direct seed identities in twelve
retained full-Gaussian models, dividing every stored numerator, denominator
and divergence witness exactly. These include D(0)=0 and zero shear coefficient.
It retains the original quotient operator in each expanded identity.

The initial run stopped at a structural equality assertion between expanded
and factored expressions, before model replay. That input and terminal failure
are preserved. The corrected assertion expands their difference and proves
it zero; the accepted run completes all checks. No previous receipt is
rewritten and no new external theorem or producer implementation is imported.

Literature reconciliation: selected reading of [Dvorsky2608.07338v1](https://arxiv.org/html/2608.07338v1), Theorem1.1, Corollary2.1 and the remarks, matches the retained DVG1 five-variable third-order example. The paper uses Q=c, while DVG1 uses its negative; nonvanishing is unchanged. Its older open-dimension remarks do not replace the current ledger. A secondary listing associates a GVC title with2610.02911, but the [primary record](https://arxiv.org/abs/2610.02911) is Taiheng Pan's machine-learning paper. That attribution is rejected, and no theorem is imported from it. Neither source closes the lower-power ternary or ordinary-Laplacian compression questions.

Historical proof bodies remain hash-identical. The present numerical
bottlenecks are tangent224 and conic216. Exact leastness,20/14 sufficiency,
expanded radical certificates, independent review, exceptional coefficient
schemes, term-count minima and lower-power GVC remain open. The qualitative
cubic case cover and both lower-witness constructions stay closed internally.
