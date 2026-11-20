# Ordinary lower-one Strong Factorial at every gap and terminal

20 November 2026. Internal written proof with complete exact finite-field replay and a separate full alternative-formula polynomial Bezout replay. All assurance remains internal.

**Theorem OrdinaryLowerOneStrongFactorial.** For every characteristic-zero field K, every integers d>=1,N>=2, and every A,B in K not both zero, the nonzero polynomial f=X(A+B X^d) has a nonzero ordinary factorial moment in the pair N-1,N. Here L(X^j)=j!. Thus no such nonzero phase has two consecutive positive moments both zero.

The proof combines a uniform complete-contour tail, the older sparse coefficient theorem, exact retained small-gap results, and one exhaustive finite-field computation. Both complete finite computations have passed; their exact mathematical scope and evidence are stated in Section8. This theorem concerns lower exponent ONE and ordinary shape beta=1; it does not prove Strong Factorial for arbitrary lower exponent, all positive Gamma shapes, arbitrary polynomial support, or the unrestricted Factorial conjecture.

**Theorem PracticalUniformGapTail.** Let d>=10 and N>=max(5000,128(d+1)) be integers. Then

    A_N(z;N)=sum_(k=0)^N binom(N,k)(N)_(dk) z^(N-k)

has only simple complex roots. Equivalently, every nonzero ordinary factorial binomial X(A+B X^d) has a nonzero moment in the consecutive pair N-1,N.

**Corollary AllLowerOneGapsAtLeastFortyFour.** The ordinary consecutive assertion holds at EVERY terminal N>=2 for EVERY gap d>=44.

**Corollary PracticalFiniteLowerOneEquivalence.** The all-terminal ordinary lower-one Strong Factorial assertion, over every positive integer gap, is equivalent to the actual finite universal-coefficient tests

    10<=d<=43, 2<=N<max(5000,128(d+1)).                (P0)

There are171812 pairs(d,N) in this deliberately redundant box, and every terminal is less than5632. The analytical argument in Sections1-7 does not execute this box; Section8 gives its separate exhaustive exact computation. No arbitrary lower exponent, all-positive-Gamma-shape, or unrestricted Strong Factorial conclusion is asserted. Rational resultant nonvanishing transfers each stated conclusion to every characteristic-zero coefficient field.

## 1. What is reused, and what is recomputed

Use the exact integral, stationary branches, deformations, curvature/height bounds and Gaussian error formula in the retained FACTORIAL_UNIFORM_LARGE_GAP_SADDLE_SIMPLICITY.md, as extended to d>=10 by FACTORIAL_GAP_TEN_SADDLE_EXTENSION.md. Put

    D=d+1>=11, H=D^D exp(-D), W=eH, c=3D/8.

The actual integral is I_N=J_N+K_N, split at c. The upper packet has domain |w|<=3W, saddle |s_J-D|<1/10, phase Phi=D log s+log(1+w/s^d)-s. The lower packet K_N/w^N has domain |w|>=W/3, saddle |s_K-1|<1/(100D), phase Psi=log s+log(1+s^d/w)-s. Set S_J=RePhi(s_J),S_K=RePsi(s_K). The inherited estimates independent of large N include

    |S_J-logH|<1/100, |S_K+1|<1/100,
    Rea_J>4/(5D), |a_J|<6/(5D),
    Rea_K>9/10, |a_K|<11/10.

The exact leading terms L_J,L_K have modulus prefactors in [2,3]/sqrt(ND) and [2,3]/sqrtN respectively. The second includes w^N. All contour estimates used below hold for N>=2 before any old enormous threshold is imposed. We do NOT apply the old theorem below its stated onset: we reuse its raw inequalities and redo the entire relative-error and derivative budget.

## 2. Smaller local upper error

The gap-ten proof actually establishes that throughout the upper local radius1/2 disk, the whole third derivative has modulus <5/D^2: its base term is <3/D^2 and its logarithmic correction is <14D(3/8)^d<2/D^2. Therefore the cubic Taylor remainder coefficient may be taken as

    B_J=1/D^2, rho_J=1/2, c_J=1/(3D),                (P1)

because Rea_J/2-B_J*rho_J >=2/(5D)-1/(2D^2)>1/(3D) for D>=11. The retained exact local error formula is e^(NS)/N times

    B/((v-rho)c^2)+1/(v(v-rho)c)+1/(v c rho),

where v=|s_*|. Using v>=D-1/10, the upper bracket is at most

    9/(D-3/5)+3D/[(D-1/10)(D-3/5)]+6D/(D-1/10)
       =6+63/(5D-3) <=375/52 <29/4.                 (P2)

This is an exact rational identity and a decreasing function of D. The lower local data remain B_K=1,rho_K=1/10,c_K=1/3; its bracket is <46. Dividing by the exact leading-term lower bounds gives local relative errors at most

    (29/8)*sqrt(D/N),       23/sqrtN.                (P3)

## 3. Integrate the entire central Gaussian tails

It is essential to retain the amplitude1/s instead of bounding the central intervals only by their length. On the entire upper central horizontal interval 3D/4<=x<=5D/4, the retained real curvature inequality gives

    RePhi(x+iIms_J)<=S_J-(x-Res_J)^2/(6D).

Also 1/|s|<=4/(3D). For b>0,rho>0,

    2 int_rho^infinity exp(-b u^2)du
      <=exp(-b rho^2)/(b rho),

by replacing1 by u/rho in each positive tail. With b=N/(6D),rho=1/2, all upper central pieces outside the local interval therefore cost at most

    16/N *exp(NS_J-N/(24D)).                         (P4)

The TWO central pieces are included, even if the finite central endpoints are not symmetric about Res_J. The whole-real tail is an upper bound on their positive modulus integrals. Thus the central relative error is 8sqrt(D/N)exp(-N/(24D)). Along N=128D this is constant-small, not falsely claimed to tend to zero with D.

Similarly the lower central interval1/2<=x<=3/2 has RePsi<=S_K-(x-Res_K)^2/6 and amplitude<=2. The retained local radius is1/10. Its two off-local pieces cost at most

    120/N *exp(NS_K-N/600).                           (P5)

Their relative error is60/sqrtN*exp(-N/600).

## 4. Every remaining contour piece

The gap-ten left upper bound remains exp(NS_J-N/10), including the vertical endpoint at c and the full horizontal[3D/8,3D/4]; their integrated amplitude is <1. The far upper ray remains at most

    10/(DN)*exp(NS_J-ND/200).

No segment is omitted. Dividing by |L_J| gives, together with(P3)-(P4),

    |epsilon_J| <=sqrt(D/N)[29/8+8 exp(-N/(24D))]
       +(1/2)sqrt(ND)exp(-N/10)
       +5/sqrt(ND)exp(-ND/200).                     (P6)

For the lower far interval[3/2,c] and its right vertical, sharpen only the final height margin. Their real modulus envelope is at most9/25; the imaginary correction is at most exp(1/(45000D))<exp(1/1000). Since S_K>=-1-1/100, a height loss1/100 follows from

    exp(1+1/100+1/1000+1/100)
       <(87/32)*(1000/979)=87000/31328<25/9.          (P7)

The first inequality uses the retained e<87/32 and exp(t)<=1/(1-t),t=21/1000. The last rational comparison has numerator margin200 after clearing denominators. Thus the old1/2000 margin can be replaced by1/100 without altering the envelope itself.

The lower origin interval and left vertical already had loss1/100. The far horizontal integrated amplitude is at most log(c/(3/2))=log(D/4)<D, and the endpoint vertical has amplitude <1. Including the origin and left vertical, all lower noncentral pieces cost at most

    4D exp(NS_K-N/100).

This deliberately loose amplitude is adequate. Consequently

    |epsilon_K| <=23/sqrtN+60/sqrtN exp(-N/600)
                     +2D sqrtN exp(-N/100).          (P8)

The integrable endpoint0 was bounded by the original entire integrand, exactly as in the retained proof; the apparent1/s phase amplitude is not integrated alone there.

## 5. Exact uniform error budget at the new onset

Assume N>=5000 and N>=128D. Since e>8/3,

    exp(-N/(24D)) <=exp(-16/3)<(3/8)^5.

The exact square comparison

    [29/8+8(3/8)^5]^2 <128*(33/100)^2                 (P9)

puts the first term of(P6) below33/100. For its left upper term, sqrt(ND)<N/10, and exp(x)>=x^4/24, so

    (1/2)sqrt(ND)exp(-N/10)<12000/N^3
       <=12000/5000^3<1/1000.

For its far upper term, exp(x)>=x^2/2 and sqrt(ND)>=1 give

    5/sqrt(ND)exp(-ND/200)
       <=400000/(N^2D^2)<=400000/55000^2<1/1000.

Thus |epsilon_J|<83/250<1/3.

For(P8), 23/sqrtN<33/100 by529/5000<(33/100)^2. Also60/sqrtN<1 and N/600>=25/3>8, so its central term is <(3/8)^8<1/2000. Finally D<=N/128 and N^(3/2)exp(-N/100) decreases for N>=5000 (its logarithmic derivative is3/(2N)-1/100<0). Since sqrt5000<71,

    2D sqrtN exp(-N/100)
       <=5000*71/64*(3/8)^50<1/1000.

Therefore |epsilon_K|<663/2000<1/3. Every bound is uniform on the ENTIRE respective analytic packet domain, and hence on neighborhoods of the overlap annulus.

## 6. Zero-free regions and derivative separation

Because |epsilon_*|<1/3, the old coarser magnitude prefactors1..4 still hold: the leading lower prefactor2 is reduced only to4/3>1, and the leading upper prefactor3 is increased only to4. Thus the original positivity comparisons prove I_N nonzero on |w|<=W/2 and on |w|>=2W. Their sufficient bound4sqrtD*2^(-N)exp(N/50)<1 follows here from D<=N/128 and N>=5000: 4sqrtD<N, exp(1/50)/2<2/3, and N(2/3)^N<1 for N>=2. The outer comparison has the smaller factor4/sqrtD and is also covered.

On W/2<=|w|<=2W a Cauchy disk of radius W/12 lies inside BOTH packet domains. The new uniform error bound gives

    |epsilon_*'/(1+epsilon_*)|<6/W.

The retained analytic-saddle/prefactor estimates give an additional <1/W for each packet. These estimates do not use the old enormous N threshold. Therefore

    |J_N'/J_N-N/(w+s_J^d)|<7/W,
    |K_N'/K_N-N/(w+s_K^d)|<7/W.                      (P10)

As in the retained proof the reciprocal main terms differ by more than1/(4W), since the lower one has modulus>1/(3W) and the upper one<1/(100W). Their N multiples are separated by>N/(4W)>14/W, because N>=5000>56. At a zero of I_N, both packets are nonzero and K_N=-J_N; hence I_N'=J_N(J_N'/J_N-K_N'/K_N)!=0. This proves PracticalUniformGapTail. The exact all-gap physical-shape resultant proves the ordinary consecutive assertion on the SAME d,N domain.

## 7. Overlap with the stronger retained sparse theorem

Use the actual older Theorem thm:sparse of lower_exponent_initial_and_boundary.tex: for integers a>=1,d>=2,N>=2, d^2>=4(a+d)log(8N) implies individual squarefreeness and PAIRWISE coprimality of ALL actual ordinary moment polynomials through N. At a=1 this is stronger than the separate1/32 Rouche wedge, and it is the input used here.

At d=44,D=45, N0=max(5000,128D)=5760, the exact inequality

    8*N0=46080<exp(1936/180)                          (P11)

has a rational Taylor lower certificate. For every d>=44 the onset is128D, and

    d/dd [d^2/(4(d+1))-log(1024(d+1))]
       =1/4-1/(4D^2)-1/D>0.

Thus the sparse theorem covers every terminal N<=N0 for every d>=44, and the new saddle theorem covers every N>=N0. This proves AllLowerOneGapsAtLeastFortyFour.

The exact retained gap1 discriminant and growing-shape prime-power theorem cover gaps1,2,3,4,5,7,8,9; the retained strict-resultant theorem covers gap6. Together these prove the reverse implication in(P0): small gaps1..9 are closed, d>=44 is closed, and gaps10..43 split into their actual finite prefix and the proved tail. The forward implication is restriction. Neither implication assumes that the finite tests have been executed. The field-transfer statement follows from rational Bezout/resultant identities, not an inequality notion over an arbitrary field.

## 8. The exhaustive remaining finite-field computation

The complete finite computation verifies all171812 cases with no nonconstant gcd, at p=1000003, in all34 gaps with maximum terminal5631. A separate complete replay constructs the coefficients by an alternative recurrence and verifies every full polynomial Bezout identity over F_1000033. Both full runs pass, with no omitted case or alternate-prime fallback. No incomplete or bounded-prefix run is used to assert the theorem.

Let p=1000003. This is prime: the accompanying finite checker tests division by every integer2 through1000; since1000^2<p<1001^2, that exhaustive trial division suffices. For every d=10,...,43 and N=2,...,max(5000,128(d+1))-1, form

    c_k=binom(N,k)*(N)_(dk), 0<=k<=N,
    T_(N,d)(t)=sum c_k t^k,
    A_N(z;N)=z^N T_(N,d)(1/z).                       (P12)

All arithmetic for the finite bridge is exact in F_p. The factorial formula for each coefficient is

    c_k = N!/[k!(N-k)!] * (N+dk-1)!/(N-1)!.

Every factorial argument is at most44*5631-1=247763<247764<p, so these denominators are units modulo p. The coefficient c_0=1, and c_N is also nonzero modulo p. Consequently A_N is monic of full degreeN, with nonzero constant coefficient; its derivative has full degreeN-1 since2<=N<=5631<p. The same full-degree properties hold for T and its derivative. Reversal interchanges the nonzero roots of A_N and T with multiplicity unchanged. Thus a checker using either coefficient orientation performs the same squarefreeness test. The principal checker constructs c_k, reverses that list to the monic A_N in ascending coefficient order, and checks A_N and its actual derivative.

For every pair in(P0), the finite checker computes

    gcd_Fp[A_N(z;N), dA_N(z;N)/dz]=1.                (P13)

Equivalently it may compute the same test on T. In either orientation, exact Euclidean arithmetic on the full polynomial degrees is used; there is no truncation, random evaluation, floating root separation or assumption about a typical parameter.

The complete range consists of29 gaps d10..38 with4998 terminal orders each, and the five gaps d39..43 with5118,5246,5374,5502,5630 orders. Hence it contains exactly

    29*4998+5118+5246+5374+5502+5630=171812

cases. The source-authenticated receipt records every enumerated pair through its complete per-gap case stream, its actual full degree, all171812 successful gcd tests, no missing case, no nonconstant gcd, and the exact proof/checker/input bytes. Small literal factorial and derivative controls verify normalization; malformed-domain and repeated-root controls are rejected. The separate full replay uses sliding rising-factorial and binomial recurrences without the principal factorial tables, works over the different prime1000033, and checks every returned full polynomial Bezout identity. It covers the same171812 pairs and34 gaps with maximum terminal5631 and no fallback. These are exact internal replays, not outside review or formal verification.

For any completed test(P13), the integer resultant Res(A_N,A_N') reduces to a nonzero element of F_p. Full degrees have been preserved, so specialization is the actual resultant and proves that this integer is nonzero. Therefore A_N is squarefree over Q and C, and the exact retained physical-shape resultant gives ordinary lower-one consecutive nonvanishing. This argument proves universal nonvanishing in the phase coefficients A,B, not only a finite collection of numerical phase examples.

The completed full computation proves every test(P0). PracticalFiniteLowerOneEquivalence therefore proves OrdinaryLowerOneStrongFactorial over C. At each fixed d,N the actual moment polynomials have rational coefficients and are coprime, so a rational Bezout identity proves the same conclusion over every characteristic-zero field. Monomial phases A=0 or B=0 have all positive factorial moments nonzero immediately.

## 9. Evidence boundary and remaining scope

The rational-constant checker certifies the changed budgets, Taylor lower bound(P11), monotonicity margins, symbolic Gaussian-tail/local-bracket algebra and deliberately failing corruptions. The separate finite-field checker authenticates the actual exhaustive(P0) computation. Their roles are distinct: neither the rational inequalities alone nor a finite root scan replaces the written uniform contour proof, and the analytical theorem alone does not certify any omitted finite case.

The closest method antecedents are the retained full saddle proofs and the older sparse coefficient theorem, whose actual hypotheses were read. Primary Rouche/Cauchy/Gaussian formulas and the already inspected trigamma series retain exactly their recorded scope. Assurance is internal mathematical audit and exact replay, not formal verification or outside review. No historical priority claim is made.

The proof closes precisely every ordinary binomial whose least exponent is1. The earlier constant-plus-power theorem treats least exponent0. Arbitrary least exponent at least2, all-positive-Gamma-shape terminal determinants, general polynomial supports, and unrestricted Factorial/Strong Factorial remain separate unresolved goals.

Proof support: [analytical checker](verify_factorial_practical_uniform_saddle_tail.py), [rational constants](factorial_practical_uniform_saddle_constants.json), [analytical source ledger](factorial_practical_uniform_saddle_source_ledger.json), and [analytical receipt](factorial_practical_uniform_saddle_verification.json); [complete finite checker](verify_ordinary_lower_one_finite_bridge.py), [finite inputs](ordinary_lower_one_finite_bridge_inputs.json), [finite source ledger](ordinary_lower_one_finite_bridge_source_ledger.json), and [finite receipt](ordinary_lower_one_finite_bridge_verification.json). The finite checker is the principal replay; the analytical checker certifies the exact numerical ledger supporting the full written contours.