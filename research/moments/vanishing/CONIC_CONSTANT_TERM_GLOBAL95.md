# A formal constant-term factor: conic94,general95,odd/homogeneous94

15 November 2026. Internal written refinement of the [general96 proof](ORDINARY_AXIS_GAUSSIAN_GLOBAL96.md),with [repeated-line92](REPEATED_EULER_SEEDS_92.md) retained. The original Wick functional,full coefficient/frame cover and analytic nonzero-zeroth exclusions are unchanged.

Every complex ternary polynomial of total degree at most three has all positive pure Gaussian moments zero if and only if its first **95** pure moments vanish. For odd or homogeneous cubics the first **94** suffice. Least intervals are general **[20,95]**,odd **[20,94]**,homogeneous **[14,94]**. Smooth-leading **[20,94]**,homogeneous-smooth **[14,56]** and both complete first-nonzero20/14 witnesses remain unchanged. These are upper bounds,not optimality.

The cubic conic exclusion improves96 to **94**. Retained parent92/rank-one95 then give full W-linear **95**. General95 follows from the complete cover below. Odd/homogeneous94 uses parity of the original centered Gaussian functional;95 does not imply94 for arbitrary non-odd polynomials.

The [checker](verify_conic_constant_term_factor.py),[certificate](conic_constant_term_factor_certificate.json),[receipt](conic_constant_term_factor_verification.json) and [source audit](conic_constant_term_factor_source_audit.json) bind auxiliary algebra. Universal formal constant-term/seed induction,full structured-phase and endpoint/frame coverage,and analytic nonzero-zeroth arguments remain independently unreviewed written inputs. No independent replay,external review,formal verification or expanded radical certificate is asserted.

## 1. Retain the full coefficient module,not just its endpoint

Only the p=q=3 conic case needs a new estimate. Its full phase is

    W=A0+sum_(k=1)^3 A_k*y^k+sum_(k=1)^3 K^k*B_k*y^(-k),
    deg A0<=3, deg A_k,deg B_k<=3-k,
    A3*B3!=0, K not identically zero, deg K=d in{0,1,2}.

Both endpoints are nonzero constants in time. The complete negative coefficients contain K,K^2,K^3 respectively. The [retained completed-row proof](CONIC_ROW_FACTORS_96_ODD104.md) preserves this module throughout actual seed formation and every twisted reduction. An endpoint-only condition would not justify what follows;its earlier counterexample and failed input stay preserved.

Set E_n=y^n for n>=0 and E_n=K^(-n)y^n for n<0. These are the polynomial completed representatives. Their multiplication is

    E_n E_m=K^min(abs(n),abs(m))*E_(n+m) if n*m<0,
    E_n E_m=E_(n+m) otherwise.

The full rank-six twisted basis is E_-3,E_-2,E_-1,1,y,y^2. The logarithmic operator is t^2*y*partial_y+y*W_y,with measure dy/y. Its high and low reductions divide only the constant endpoints3A3 and-3B3 in this completed module. Hence the completed normal coordinates are polynomial in K and the coefficient parameters. Rational witnesses in the original y coordinates remain permitted and retain the full divergence identity.

## 2. A coefficientwise formal constant-term functional selects the unit

Treat K as independent and work over the coefficient field with t,A3,B3 inverted. Define

    pi(f)=CT_y[ exp((W-A0)/t^2)*f ].

This is coefficientwise well-defined as a formal K series on the completed module. The positive exponential contains nonnegative y powers. At each fixed K order,the negative exponential contains only finitely many negative y powers,so only finitely many positive coefficients can contribute to the constant term. No contour or convergence assertion is used.

For an explicit formula write exp(sum A_k*y^k/t^2)=sum P_a*y^a and exp(sum B_k*z^k/t^2)=sum Q_b*z^b. Then

    pi(E_n)=sum_(l>=max(0,n)) P_(l-n)*Q_l*K^(l+max(-n,0)).

Thus pi(E_n) has K order at least abs(n),pi(1)=1+O(K),and at K=0 this functional selects exactly the unit row of the completed basis. Formal constant terms kill the full logarithmic twisted derivative:

    pi(t^2*y*partial_y f+y*W_y*f)=0.

This follows by taking the constant term of y*partial_y(exp((W-A0)/t^2)*f). It remains valid for the Laurent-series witnesses produced by the original bounded divisions. It is not a commutative critical-ideal reduction.

## 3. Keep every time and normalization derivative in the seed recurrence

Retain the actual negative Euler representative

    h=t^3*K, B=t*K_prime-2K/3,
    L f=h*f_prime+t^2*B*y*f_y+S*f,
    S=K*(t*W_prime-2W-t^2)+B*y*W_y,
    R0=1, R_(j+1)=L R_j-j*h_prime*R_j.

The original normalized scalar functional is t^(-1)exp(A0/t^2)*pi. The Euler difference is the exact logarithmic derivative B*(t^2*y*partial_y+y*W_y). Taking pi therefore intertwines L with h times the scalar time connection obtained by differentiating that normalized functional.

This scalar connection preserves formal K series:coefficient and K derivatives are all retained,and differentiation of K^r,r>=0,introduces no negative power. The normalization contributes -1/t+A0_prime/t^2-2A0/t^3,also regular in K. In particular no determinant,denominator,vertical-field or h_prime term has been omitted.

Exact cancellation of the h_prime terms gives

    pi(R_j)=h^j*nabla^j(pi(1)).

Hence pi(R_j) is divisible by K^j. For j>=1,its constant coefficient in K is zero. Write the actual completed normal form as sum n_(k,j) E_k. All n_(k,j) are polynomial in K. Since pi at K=0 selects the unit row,

    K divides n_(0,j), for every j>=1.

The statement follows from the universal recurrence and formal functional,not from the finite sample of row divisions. Time coefficients may subsequently specialize to arbitrary allowed polynomials,and K may be constant,linear,quadratic or have repeated/shared roots. Polynomial identities established before specialization retain the same divisibility. Identically zero K or vanishing endpoints stay in the already covered lower actual-order strata.

## 4. Divide one common factor from the scalar relation

At actual seed rank s<=6 choose any s rows containing the unit and form polynomial cofactors C_l of columns N0=unit,N1,...,Ns. Every entry of the unit row in N1,...,Ns has a K factor,so

    K divides C0.

For l>=1 no additional raw-minor factor is required. The scalar coefficient is A_l=H_l*C_l,where H_l=h^l=t^(3l)K^l in this endpoint-constant case. Thus every scalar coefficient has the common K factor. Divide it from the entire relation.

The retained analytic nonzero zeroth coefficient stays nonzero. Higher coefficients retain t^(3l)K^(l-1),so their valuations are at least3l even if K(0)=0 or K shares a root with t. No other factor is cancelled and no derivative-order allowance is lost. No time-dependent basis is silently differentiated:the completed rows are those of the retained actual seeds,while the formal functional was differentiated with the complete connection above.

The retained completed-row degree bounds are3j+abs(k),4j,5j-abs(k) for K degrees0,1,2. Dividing the unit row factor lowers each selected zeroth cofactor degree by d=deg K. Every actual rank and unit-containing row subset is included,so the p3q3 gates improve

    72,84,96  ->  72,83,94.

For nonzero constant K there is no degree saving;the72 gate is retained. The other conic strata retain their existing bounds:p2q3<=94,p1q3<=84,p0q3<=69,p3q2<=93,p2q2<=93,and all other cases<=63. Therefore the full cubic conic exclusion is94. This is a cubic-G statement,not an arbitrary-degree claim.

## 5. Reassemble the complete frame and parity cover

Full W-linear now has gate max(conic94,parent92,rank-one95)=95. Combine it conservatively with the retained complement bounds:node87,transverse77,triangle82,cusp77,tangent66,with their constant-D strata and every zero-axis/frame face retained. A valid complete cover is

    quadratic7,smooth94,node95,cusp95,transverse95,tangent95,
    triangle95,concurrent92,repeated92.

Its maximum95 proves the general test with no extra frame hypothesis. For odd P,the original centered Gaussian has E(P^95)=0 automatically. Thus the first94 zero moments supply all95 required by the general theorem. A homogeneous cubic is odd,so94 applies there as well;homogeneous degrees at most two already have gate7. No oddness is imposed after an affine change of coordinates:parity is used on the original polynomial and Gaussian functional.

Pure nullity retains the internally closed ternary cubic GMC classification and sharp mixed cutoff m>3deg Q;odd cutoff stays m>deg Q. Conic94/full-W-linear95 do not establish optimality. In the rational twenty-coefficient ring radical(I95)=radical(I_infinity),without ordinary ideal equality or expanded radical cofactors. At most7 exact monotone radical-search rounds against I95 determine the least general prefix in principle;none has run on the full coefficient ring. Equality of consecutive moment ideals is insufficient.

## Evidence and remaining assurance

The checker verifies36 formal logarithmic Ward coefficients through K order3,33 unit-row divisions in11 retained structured families,and30 complete original logarithmic seed identities through rank6 in five new families. The new families include constant/linear/quadratic K,simple zeros and a double zero. All five full scalar relations with35 coefficient checks have nonzero zeroth cofactors and retain the stated valuations after common cancellation. All actual-rank/row degree arithmetic is checked. No old checker main is rerun and no new input failed.

Universal coefficientwise formal constant-term and all-j seed induction,full structured coefficient module,complete endpoint/frame/metric cover,the seed-to-period bridge and analytic nonzero-zeroth exclusions remain independently unreviewed written inputs. Finite identities support auxiliary algebra and do not certify the complete analytic theorem. Every earlier proof body and receipt remains byte-identical. No fresh literature-coverage claim,external theorem or package is imported;the preceding selected primary-record refresh remains retained.

Least general[20,95],odd[20,94],homogeneous[14,94],20/14 sufficiency,expanded I95 cofactors,independent analytic review,exceptional schemes,global term minima and lower-power GVC remain open. Rank-one/full-W-linear95 and smooth/conic94 now control the numerical frontier. Qualitative cubic coverage and both lower witnesses remain closed internally.
