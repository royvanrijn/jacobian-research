# Quantitative claims and bounds

<!-- Generated from MATH_STATUS.json; run make render-navigation. -->

Bounds retain their exact domains and assumptions. Numeric comparisons use identical quantities, domains and parameter restrictions; conditional bounds stay separate. Symbolic criteria are displayed without guessing implication or asymptotic order.

[Dependency graph](DEPENDENCIES.md) · [Source coverage](COVERAGE.md) · [Exact scopes](SCOPES.md)

| Claim | State | Quantity | Domain | Bound | Assumptions |
|---|---|---|---|---|---|
| DIX-QCR1 | proved | total-degree | classical fourth image completing the fixed retained Long triple | lower: 23 |  |
| DIX-QCR1 | proved | momentum-order | classical fourth image completing the fixed retained Long triple | exact: 5 |  |
| DIX-QPH1 | proved | h-coefficient-total-spatial-degree | specific positive scalar branch of quintic phase Long lift, coefficient h^n for n>=2 | exact: 10*n-1 |  |
| DIX-QPH1 | proved | h-coefficient-momentum-order | specific positive scalar branch of quintic phase Long lift, coefficient h^n for n>=2 | exact: n |  |
| JAC-PF153 | proved | ideal-nilpotence-exponent | 30-generator cone ideal in actual fixed-parent degree5 first-jet trace quotient | upper: 4096 |  |
| JAC-PF156 | proved | ideal-nilpotence-exponent | residual obstruction(c) in actual fixed-parent degree5 first-second-jet quotient | exact: 3 |  |
| JAC-PF157 | open | reduced-ruling-degree | actual complex plane Keller target in order-g nontrivial generic cyclic-deck rational pole atlas | lower: 2*g |  |
| JAC-CUBIC12 | proved | nonlinear-monomial-occurrences | integral cubic Keller collision with fixed dimension12; {"dimension": 12} | upper: 36 |  |
| JAC-CUBICLINEAR89 | proved | dimension | integral cubic-linear Keller map with rational collision | upper: 89 |  |
| JAC-CANONICAL-FILTER | proved | source-bundle-rank | canonical invariant-filter receiver for rankN finite etale algebra | exact: 3*N+3 |  |
| JAC-CANONICAL-FILTER | proved | generic-map-degree | canonical invariant-filter receiver for rankN finite etale algebra | exact: 3^(N+1) |  |
| JAC-CANONICAL-FILTER | proved | coordinate-degree | canonical invariant-filter receiver in split coordinates | exact: 4*N |  |
| FOUND-F1 | proved | dimension | ordinary characteristic-zero polynomial Keller map with rational collision | upper: 3 |  |
| DIX-RANK3-FINITE-COUNTEREXAMPLE | proved | weyl-generator-total-degree | finite nonsurjective A3 cotangent endomorphism from normalized F1; {"rank": 3} | upper: 12 |  |
| DIX-LEADING-POISSON-FINITENESS | open | corrected-fourth-image-degree | rank2 actual finite Weyl completion with leading-symbol independence and module-finite Poisson centralizer; {"rank": 2} | upper: max(e2,e1+e2-2+d1) |  |
| ARITH-AMBIENT-HASSE-LOWER-BOUND | proved | ambient-hasse-target-count | all rational targets of the displayed fixed Phi in actual projective target height | lower: (Cphi-epsilon)*B^2/log(B) eventually for every epsilon>0 |  |
| ARITH-AMBIENT-HASSE-LOWER-BOUND | proved | scaled-count-liminf | all rational targets of the displayed fixed Phi, A_Phi(B)*log(B)/B^2 | lower: Cphi |  |
| JAC-CANONICAL-FACTOR-LENGTH | proved | maximum-geometric-nonautomorphism-factor-length | canonical rankN invariant-filter receiver under identity stabilization | exact: N+1 |  |
| JAC-CANONICAL-FACTOR-LENGTH | proved | canonical-volume-bundle-rank | canonical rankN invariant-filter receiver plus determinant-dual line | exact: 3*N+4 |  |
| GAUSSIAN-POLYNOMIAL-SEED-CONSTRUCTOR | proved | phase total degree | Exact seven-variable P_(f,M),deg f=L>=1,M>=0,commutative Q-algebra | upper: L+2M+5 |  |
| GAUSSIAN-PURE-MARKER-DELAY | proved | least bad original-marker degree | Exact five-variable pure Q_(M,k) construction,integers k>M>=0,characteristic-zero fields | exact: M+1 |  |
| FACT-ARITHMETIC-FRONTIER-TERMINAL-ORDERS-THM-MINUSONE | proved | certified terminal orders | Ordinary lower-one gap d=p^s,p prime,r,s>=1,1<=m<p,N>=2 | range: N=m p^r-1 |  |
| FACT-ARITHMETIC-FRONTIER-TERMINAL-ORDERS-THM-NATIVEPRIME | proved | ordinary binomial coprimality criterion | Integers a>=0,d>=1,N>=2;ordinary factorial moments | criterion: Exists prime p>dN and j in1..d with v_p(aN+j)=1;in particular a>=d and aN+1 prime |  |
| FACT-ARITHMETIC-FRONTIER-TERMINAL-ORDERS-THM-COMMONQUOTIENT | proved | ordinary binomial coprimality criterion | Gap d=p^s,p prime,s>=1,a>=0,N>=2 | criterion: floor(aN/d)=floor(a(N-1)/d)=q and q+N<p |  |
| FACT-ARITHMETIC-FRONTIER-TERMINAL-ORDERS-THM-COMMONQUOTIENT | proved | resultant p-adic valuation | Same common-quotient criterion; b=(d-1)/(p-1) | exact: (N-1)v_p((aN)!)+N v_p((a(N-1))!)+bN(N-1) |  |
| FACT-ARITHMETIC-FRONTIER-TERMINAL-ORDERS-THM-PRECISIONOBSTRUCTION | proved | maximum coefficient-comparison precision over unit scalings | N=9h+4,h>=0,3 does not divide h+1,gamma=3h+2;the exact gap3 comparison family | exact: 3 |  |
| FACT-ARITHMETIC-FRONTIER-TERMINAL-ORDERS-THM-PRECISIONOBSTRUCTION | proved | comparison discriminant valuation | Same precision-obstruction family N=9h+4,h>=0,3 does not divide h+1 | lower: N+1>3 |  |
| FACT-GAP-THREE-HURWITZ-THM-BLOCK | proved | whole infinite H3 total-nonnegativity terminal range | Actual W_(n,n+beta-1),every real beta>=0 | range: 2<=n<=16 |  |
| FACT-GAP-THREE-HURWITZ-THM-INITIAL | proved | uniform positive initial Euclidean pivots | Actual gap-three family,every integer n>=4 and real beta>=0 | lower: 9 |  |
| FACT-LOWER-EXPONENT-INITIAL-AND-BOUNDARY-THM-GAPBOUND | proved | sufficient binomial gap | Ordinary factorial M_(j,a,d),integer a>=1,1<=j<=N,N>=2 | upper: D_N=2(N-1)ceil(log2(4N(N-1))) |  |
| FACT-LOWER-EXPONENT-INITIAL-AND-BOUNDARY-THM-LOWERBOUND | proved | sufficient lower exponent | Ordinary factorial M_(j,a,d),fixed N>=2,d>=2,1<=j<=N | upper: A(N,d) as defined in thm:lowerbound |  |
| FACT-LOWER-EXPONENT-INITIAL-AND-BOUNDARY-COR-RECTANGLE | proved | remaining failure domain | Ordinary factorial binomial pairwise coprimality through fixed N>=2,integer a>=1,d>=2 | range: 1<=a<A_N;2<=d<D_N;A_N=(N D_N^2+D_N)(4N+1)^N |  |
| FACT-LOWER-EXPONENT-INITIAL-AND-BOUNDARY-THM-SPARSE | proved | binomial squarefree and pairwise-coprime criterion | Ordinary factorial binomials,integers a>=1,d>=2,N>=2,orders1..N | criterion: d^2>=4(a+d)log(8N) |  |
| MOM-MANUSCRIPT-THM-CUBIC | proved | uniform Gaussian mixed cutoff slope | Pure complex ternary phases of total degree<=3,arbitrary fixed polynomial marker Q | exact: 3 |  |
| IMAGE2-CUBIC-ODD-CLASSIFICATION | proved | uniform Image mixed cutoff slope | Pure odd sphere cubic stratum,lifted to balanced two-pair Image phases,arbitrary multiplier Q | exact: 1 |  |
| IMAGE2-CUBIC-TEN-PARAMETERS | proved | sufficient pure-moment prefix | Exact ten-parameter sphere cubic slice,all complex coefficients including a=0 | upper: 8 |  |
| IMAGE2-CUBIC-ELEVEN-PARAMETERS | proved | sufficient pure-moment prefix | Exact eleven-parameter sphere cubic slice,all complex coefficients including a=0,h=0 | upper: 10 |  |
| IMAGE2-CUBIC-ALL-MULTIPLIERS | proved | uniform Image mixed cutoff slope | All pure bihomogeneous balanced two-pair phases of bidegree(3,3),arbitrary polynomial multiplier Q,characteristic-zero field | exact: 3 |  |
| IMAGE2-MINIMUM-BALANCED-DEGREE | proved | least failing balanced degree | Two-pair Special Image problem,characteristic-zero fields,bihomogeneous bidegree(d,d) | exact: 4 |  |
| IMAGE2-QUADRATIC-SIX-MOMENTS | proved | sufficient pure-moment prefix | Complete nine-parameter sphere quadratic family over C,including all defective branches | exact: 6 |  |
| IMAGE2-SEVEN-PARAMETER-SIX-MOMENTS | proved | sufficient pure-moment prefix | Exact seven-parameter sphere phase sliceS1 | upper: 6 |  |
| IMAGE2-OPPOSITE-QUADRATIC-PROFILES | proved | sufficient pure-moment prefix | Exact opposite-quadratic-height profile sliceP1 | upper: 10 |  |
| IMAGE2-REPEATED-EXCEPTIONAL-FIELDS | proved | sufficient pure-moment prefix | Repeated-profileE1 at r=0 or roots of P2*P3*P4,normalized cubic coefficientzero | upper: 12 |  |
| MRC-PERMUTATION-INDEX-BOUND | proved | meridional rank | Knot with a transitive degree-r permutation quotient and meridian index delta=r-number_of_cycles | lower: ceil((r-1)/delta) |  |
| MRC-ALL-CHARACTERISTIC-BETTI-GAIN | proved | group generator rank | Finitely generated group G,index-r subgroup H,r>=2,any coefficient field | lower: 1+ceil((b1(H)-b1(G))/(r-1)) |  |
| MRC-WHITEHEAD-TORUS-APPLICATION | open | meridional rank and bridge number | Cited fixed untwisted Whitehead pattern iterated N times on coprime torus knot T(p,q),2<=p<q | exact: p*2^N |  |

## Comparable improvements

No weaker registered numeric bound has an available stronger bound in an identical domain.
