# Pole filtrations: global199, odd198 and homogeneous198

15 November 2026. Internal written refinement of the [common-time224 proof](COMMON_TIME_FACTOR_GLOBAL_224.md), retaining the original Gaussian Wick functional, complete quotient bases, coefficient/class cover and analytic exclusions.

For every complex polynomial P of total degree at most three in three standard real Gaussian variables,

    E(P^j)=0 for j=1,...,199 iff all positive pure moments vanish.

There is no extra Gaussian-frame hypothesis. Pure nullity gives the retained W-linear frame and sharp mixed cutoff m>3*deg Q. The least global prefix is in **[20,199]**, the odd-only prefix in **[20,198]**, and the all-homogeneous prefix in **[14,198]**. Smooth-leading **[20,94]** and homogeneous-smooth **[14,56]** are unchanged. The first-nonzero20/14 witnesses are complete. These are sufficient bounds; leastness and20/14 sufficiency remain open.

The [checker](verify_pole_filtered_seed_bounds.py), [certificate](pole_filtered_seed_bounds_certificate.json), [receipt](pole_filtered_seed_bounds_verification.json) and [source audit](pole_filtered_seed_bounds_source_audit.json) support auxiliary exact algebra. All-degree inductions, complete chart/rank coverage, period comparison and analytic nonzero zeroth coefficients remain written inputs. This is not independent, external or formal verification.

## 1. Row scaling after computing actual seeds

In a fixed full quotient basis b_i including1, write an actual derivative seed as N_(i,j)/H_j. If every row numerator contains a known power D^(e_j+rho_i), while H_j and every exact-divergence witness contain D^e_j, cancel the common factor and express the same class in the row basis D^rho_i*b_i. Over C(t), D is already invertible on the retained complement. This is an exact row operation on completed seed columns, not computation of new seeds by differentiating a rescaled basis. The latter would require extra derivative terms.

The same principle applies to powers of t. The scalar cofactor relation concerns the unchanged actual seeds. Its analytic nonzero zeroth coefficient therefore comes from the same retained exclusion. The higher coefficients contain the new H_j. We check their valuations dominate derivative-order loss, so the existing initial-order contradiction applies without an extra rank penalty.

## 2. Tangent complement: a second filtration gives120

Use the common-time normalized tangent connection, delta=t^2*D, h=t^3*D^3, and phase

    W=u^4/2+D*u*v^2+R,
    deg_joint(t,u,v) R<=4, weight_(2,3) R<=7.

To prove degree savings uniformly when D is constant, first work in the formal graded polynomial ring with independent t,D, each of degree1, and derivation partial_t+eps*partial_D, eps of degree0. Allow R to have arbitrary coefficients in this larger ring subject to the same joint bound4. Actual D=m-2*alpha*t and constant nonzero D are specializations; such specialization cannot increase joint degree. Division by a specialized constant alone would not prove the saving.

For a monomial t^a*D^b*u^i*v^k define its pole grade F=i+2k-b. Then F(W)<=4. With beta_u=delta_prime/4, beta_v=tD-delta_prime/8 and C=[u^2*v]S_beta, the source after shear and common-time cancellation is

    Sbar=(D*S_beta-C*u*W_v/2)/t.

It has joint degree<=6, old weight<=6 and F<=3. Indeed the only old-weight<=6 monomial of pole grade4 is v^2; its source coefficient is divisible by D, and the shear u*W_v cannot contribute to v^2. Every other monomial has pole grade<=3. The Euler-vector coefficients have D valuation>=2; the shear coefficient has valuation>=1 and u*partial_v lowers pole grade by1. Also F(h)<=-3 and F(h_prime)<=-2. The complete actual recurrence therefore gives F(R_j)<=3j, joint degree<=6j, and old weight<=6j for every j.

Use the unchanged principal tangent division, cleared by A0=D^2. For an input exponent e=(i,k), rho(e)=i+2k, its remainder and witnesses obey

    F(remainder)<=rho(e)-2,
    F(Cu)<=rho(e)-5, F(Cv)<=rho(e)-4,
    joint_degree(remainder)<=i+k+2.

These follow from the retained explicit principal formulas: each replacement decreases the original (2,3) weight, multiplication by a lower derivative R_u or R_v has pole grade at most3 or2, and the delta-divergence terms are smaller. Thus the same inequalities propagate through every reduction, including padding basis rows by D^2. The old weight termination bound supplies depth r_j=6j-4. The cleared normal coefficient in row e is divisible by

    D^(2*r_j-3j+rho(e))=D^(9j-8+rho(e)).

The divergence witnesses contain at least D^(9j-8); their stronger pole bounds give this directly. Cancelling e_j=9j-8 and scaling rows by D^rho(e) gives

    Htilde_j=t^(3j)*D^(6j),
    deg Ntilde_(e,j)<=9j-|e|-rho(e).

The full basis is 1,u,u^2,v,v^2. Its combined grades |e|+rho(e) are0,2,4,3,6, summing15. At full rank5 the zeroth cofactor degree is9*15-15=**120**; every smaller actual seed rank and row subset containing the unit is bounded by120. Each higher scalar coefficient retains valuation>=3j, including D(0)=0. No coprimality or new coefficient open is imposed. The formal graded argument covers constant D as well.

## 3. Mixed parent: t-adic filtration gives full92

Keep the entire fourteen-coefficient mixed family, with interior normalized phase

    F=z/t-a-z^2, V=t^2*(t*G(z,F)-F^2/2), delta=t^2.

Its actual degree d is4,5 or6. For d5/6 retain the exact ordinary-volume Euler representative from the structured proof:

    S_E=t*V_t-2V-t^2-t^2/d-z*V_z/d,
    R_(j+1)=t^3*R_j_prime-3j*t^2*R_j
                            -t^2*z*(R_j)_z/d+S_E*R_j.

For d4 use the original source S=t*V_t-2V-t^2 and omit both Euler terms. The phase, every source coefficient and its leading derivative coefficient have t-degree<=3. Define pole grade as z exponent minus t exponent. Since H=z-a*t-t*z^2 has pole grade1, each original cubic contribution t^(3-k)*z^i*H^k has grade i+2k-3. The Gaussian H^2 has grade2. It follows that phase/source grades are at most a_d=2,2,3 for d4/5/6, and the full seed recurrence gives grade<=a_d*j.

Let A=leading_z(V_z). Its t valuation ell_d is2,3,3. For d4, A=t^2*L(t), L(0)!=0 and deg L<=1; for d5/6, A is a nonzero constant times t^3. Leading reduction preserves pole grade: the derivative of V has grade<=a_d-1, and its leading monomial attains d-1-ell_d=a_d-1. Clearing by A at each step changes the grade bound by -ell_d. Hence at depth r_j, the coefficient in row z^i is divisible by t^(ell_d*r_j-a_d*j+i). The exact witness has at least the common factor without the extra row exponent, by the same division estimate.

| d | depth r_j | common exponent e_j | new numerator bound in row i | zeroth cofactor gate |
|---|---|---|---|---:|
| 4 | 4j-2 | 6j-4 | 9j-2-i | 45 |
| 5 | 4j-3 | 10j-9 | 5j-i | 44 |
| 6 | 5j-4 | 12j-12 | 6j-i | 80 |

Divide all witnesses and H_j=t^(3j)*A^r_j by t^e_j and scale the completed rows by t^i. The degree bounds follow from the old3j+3r_j numerator bound after these exact divisions. The denominators are t^(5j)*L^(4j-2), constant*t^(5j), and constant*t^(6j), respectively. Their valuations are at least5j and their degrees satisfy the corresponding bounds. The cofactor bounds use every row subset containing0; at full ranks3/4/5 they are9*6-2*3-3=45,5*10-6=44,6*15-10=80.

The c=0 face retains35 and the d=0 linear-WT face retains92. Thus the **full mixed parent has sufficient gate92**, including all coefficient/rank boundaries. No new integral, fibrewise basis or normalization is used.

## 4. Conic: count each reduction step by its displacement

Retain the structured phase

    W=A0+sum_(k=1)^3(A_k*y^k+K^k*B_k*y^(-k)),
    deg A_k,deg B_k<=3-k, deg K<=2,
    B=t*K_prime-2*K/3, h=t^3*K,
    S_K=K*(t*W_t-2W-t^2)+B*y*W_y.

The actual recurrence includes t^2*B*y*partial_y and gives support[-2j,pj], with negative coefficients divisible by K^(-k). The scalar and derivative Euler terms must both be retained. K identically zero belongs to a lower actual endpoint. Cases q<3 retain maximum145; p<=1,q3 retain130. The p3q3 cases deg K1/2 retain199/108. Two improvements complete the cover.

When p=q=3 and K is nonzero constant, normal coefficients are polynomial and deg_t R_(k,j)<=3j-k. High reductions preserve joint t/y grade and never produce an exponent below-3. A low cancellation of exponent n<-3 replaces it by exponents k>n. Its coefficient t-degree increase is at most k-n: the negative phase terms have degrees1/2 for displacements1/2; the Euler term has degree2 for displacement3; positive phase terms have degrees2/1/0 for displacements4/5/6. Following a low path from n to a final row k therefore gives degree<=3j-n+(k-n)<=7j+k. Low reductions never produce exponents above2. There is no low input at j1. Thus the first column has row degree3-k and all later columns have row degree7j+k. For selected rows containing0 at seed rank s>=2, the determinant degree is at most

    3+7*sum_(j=2)^s j+sum_selected k-2*min_selected k.

For all six rows -3,...,2 this is146; every smaller row choice is below146. The denominator h^j still has valuation3j.

For p2,q3, split high and low paths. A high step divides by2*A2, where deg A2<=1, but its output is never below-3. Original low inputs have the structured K divisibility, so low division by -3*K^3*B3 is polynomial and never creates exponent above1. Consequently only A2^(2j-1) is required as a common normal denominator; there is no negative-endpoint padding. This remains valid when A2 shares a factor with K or t: the high part is treated rationally and cleared at the end, and is never subjected to low polynomial division. No coprimality is assumed.

Use the old integer weights and the safe whole-step costs below. For a chosen row k, weight_t*deg N_(k,j)+weight_y*k<=b_j. Here high=2j-1, low=max(0,2j-3), and b_j=seed_cost*j+high_cost*high+low_cost*low. High cost includes common A2 padding; the deg K1 value4 is a conservative bound.

| deg K | weight(t),weight(y) | seed cost | high cost | low cost | full gate |
|---|---|---:|---:|---:|---:|
| 0 | 1,1 | 3 | 1 | 6 | 171 |
| 1 | 2,1 | 8 | 4 | 3 | 136 |
| 2 | 1,1 | 5 | 1 | 0 | 105 |

All actual seed ranks and row subsets of -3,...,1 containing0 are included in floor((sum b_j-weight_y*sum k)/weight_t). H_j=h^j*A2^(2j-1) satisfies the same bound and has valuation>=3j. Thus the **uniform cubic conic exclusion is199**. Together with parent92 and rank-one140, the **full sixteen-coefficient W-linear gate is199**. This is not an arbitrary-degree G statement.

## 5. Tangent frame geometry and the complete cover

On a rank-one tangent frame chart write the leading cubic

    F3=f*W*T^2+g30*Z^3+g21*Z^2*T+g12*Z*T^2+g03*T^3, f!=0.

The tangent cone at the W-axis is f*T^2. For a line and a smooth conic tangent there, the line is T=0, so g30=0. Conversely the remaining conic has determinant -f^2*g21/4, and smoothness forces g21!=0. Thus in the actual conditional Gaussian phase Phi=t*G(t*A(T)/q,T)-T^2/2 the absent Z^3 term implies T-degree<=5. The retained rank-one degree table18/50/93/132/140 can therefore stop at5, giving132 on this tangent chart. This restriction is for genuine line-smooth-conic tangency; the unrestricted rank-one gate remains140. The other tangent frame face uses parent92. The **full tangent gate is max(120,132,92)=132**.

| Leading type with arbitrary lower terms | Full sufficient gate/exclusion |
|---|---:|
| Degree at most two | 7 |
| Smooth cubic | 94 |
| Irreducible node | max(199,199)=199 |
| Irreducible cusp | max(183,92,140)=183 |
| Transverse line-smooth-conic | max(199,199)=199 |
| Tangent line-smooth-conic | 132 |
| Three independent lines | max(189,199)=199 |
| Three distinct concurrent lines | 104 |
| Repeated leading line(s) | 99 |

The maximum199 proves the global pure-only test. For odd cubics moment199 vanishes automatically, so198 suffices. The full homogeneous cover is max(smooth56,cusp10,reducible10,nonisotropic-node72,isotropic-node conic199)=199; homogeneity makes the cubic odd, giving198 there too. The generic conic and generic W-linear gates remain199, since their polynomials need not be odd. In the rational twenty-coefficient ring radical(I199)=radical(I_infinity), without ordinary ideal equality or expanded global cofactors. At most ceil(log2(199-19))=8 exact monotone search rounds against I199 would find the least K; none has run.

## Evidence and boundary

The checker verifies formal-D source/pole estimates, finite principal identities supporting the explicit all-degree division formulas, generic parent source and leading-valuation algebra, exact seed reductions in five full-Gaussian tangent models (including constant D and D(0)=0), three parent models and eleven structured conic models, and all determinant/rank/row arithmetic. Shared A2/K/t factors are included in conic models. Every new seed certificate uses the full original twisted divergence operator.

The first input incorrectly substituted an Euler representative that did not cancel the negative conic endpoint. It failed the support assertion after the tangent/parent checks. The failed input and terminal receipt are retained. The accepted checker uses the exact retained B=t*K_prime-2*K/3 representative; both its source and vertical derivative are included. A second run rejected a shared-factor fixture whose A2 had degree2, outside the stipulated degree<=1. Valid shared-factor fixtures replace it. A third run passed all algebra and stopped at JSON serialization of a SymPy integer; the accepted serializer converts those exact integers. Both failed inputs and receipts are retained. No prior proof body or receipt is rewritten. No fresh literature claim or external theorem is imported by this algebraic refinement; the selected primary-source reconciliation in the preceding audit remains the literature boundary.

Leastness,20/14 sufficiency, expanded global radical cofactors, independent analytic review, exceptional coefficient schemes, term-count minima and lower-power GVC remain open. The qualitative cubic case cover and the completed lower witnesses stay closed internally.
