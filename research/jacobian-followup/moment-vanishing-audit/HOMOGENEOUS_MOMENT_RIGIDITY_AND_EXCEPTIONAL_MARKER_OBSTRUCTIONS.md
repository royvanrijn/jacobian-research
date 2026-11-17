# Homogeneous Gaussian moment rigidity and exceptional marker obstructions

16 November 2026. Internal written proofs. External review/outreach remains deferred under the closure-first priority; independent replay and formal verification are not claimed. This note closes two structural questions in full stated scope: a global homogeneous Gaussian growth/degree theorem,and the direct joint-pure torus-compatible marker transfer into either exceptional quotient. It does NOT decide either full exceptional Mathieu--Zhao problem.

## 1. Neutral rigidity forbids a mixed-charge jointly pure algebra

Let A be an integral domain over a characteristic-zero field,graded by a torus character lattice Z^r. Let h:A->K be invariant and normalized. Assume the eventual pure radical on the weight-zero algebra is zero:

    b in A_0 and h(b^m)=0 eventually imply b=0.             (1)

Let f_1,...,f_s be nonzero eigenvectors. Assume every positive joint monomial has zero h-value. Then

    0 is outside conv{weight(f_1),...,weight(f_s)}.          (2)

Otherwise rational convexity gives nonnegative integers n_i,not all0,with sum n_i*weight(f_i)=0. The nonzero product b=prod f_i^n_i lies in A_0,and every h(b^m) is0 by joint purity. This contradicts(1). The integral-domain and NONZERO-generator assumptions are essential.

Conversely strict rational separation makes all positive joint monomials have nonzero weight,so they are h-null. For any fixed marker with finite weight support,the weights of powers of every zero-constant polynomial in the f_i escape that support. The whole augmented algebra is marked-safe. With a compatible degree filtration and a separating integer functional ell,if ell(weight f_i)>=eta>0 and every degree-e marker has ell-weights in[-w e,w e],then m*eta>w*e is sufficient.

For a one-dimensional torus,every such jointly pure family has weights of ONE strict sign. More intrinsically,if a finitely generated torus-stable subalgebra B admits h|B as an augmentation character,choose finite weight components of algebra generators. They belong to B,generate it,and their nonconstant components lie in the augmentation ideal. A weight-zero component in that ideal is zero by(1). The remaining eigenvector generators satisfy(2). This classifies torus-stable jointly pure subalgebras under hypothesis(1);it does not show the torus saturation of an individual pure phase is jointly pure.

## 2. The exact obstruction on BOTH exceptional quotients

For Gamma=2O or2I,put A_Gamma=C[SL2]^Gamma with right Gamma action and Haar h. The retained theorem [GTORUSPURECLOSED20261115](VALUATIVE_BOUNDARY_VOLUMES_AND_TORUS_INVARIANT_CLOSURE.md) proves(1) for EVERY maximal left torus and for every algebraic-automorphism-conjugated such action. It is an internal theorem with its original analytic assurance boundary;that boundary is retained here.

Therefore there are NO two nonzero original invariant functions S,T of opposite torus weights with all positive joint Haar moments zero. In particular ST would be a nonzero torus-fixed pure phase,already forbidden by(1). More generally no mixed-weight torus-stable jointly pure finitely generated algebra can occur. This holds in all degrees,with arbitrary right-invariant multiplicities,and after Haar-preserving algebraic automorphisms.

It rules out the DIRECT transfer of the completed six-variable cotangent marker architecture,for ANY radial moment sequence,when its proposed generators are original exceptional functions and its grading is a left torus or its automorphic conjugate. The failure is prior to choosing a marker:the required joint purity itself is impossible. A new polynomial marker or a different prime-separating sequence cannot repair that failed hypothesis.

This is a precise obstruction,not a proof of either full MZ case. A mixed-isotypic pure P may have both torus signs while its individual weight pieces are NOT jointly pure. Its torus saturation need not preserve the cancellations in h(P^m). The remaining strict-leading,image-interior sector therefore survives this obstruction. Noninvariant generators whose selected polynomial becomes invariant,a different non-torus grading,or a non-jointly-pure presentation need separate proofs;they are not excluded by the statement above.

There is a second obstruction that does NOT require invariant generators. Let Ptilde be ANY pure actual homogeneous four-Gaussian invariant lift. Every Gaussian marker of degree<8 on2O,or<12 on2I,has ALL positive marked moments zero. Reynolds averaging preserves degree. Below that threshold the original quotient has only constants,so each even homogeneous averaged marker is a scalar determinant power;odd components vanish by the binary centre. Radial integration then makes every marked moment a nonzero radial factor times h(P^m)=0. Thus the proposed quadratic bad marker cannot detect any such invariant pure phase,even if one tried to present that phase using noninvariant generators. This is an actual all-marker degree obstruction,not an assumption about an invariant marker of the same displayed formula.

## 3. Exponential denominator bounds for every rational sphere moment

Let N>=1 and let sigma be normalized uniform measure on S^(N-1),including the two-point sphere when N=1. Odd coordinate exponents have moment0. For alpha_i>=0,write m=sum alpha_i. The even sphere moment is

    M_alpha=prod_i(1/2)_(alpha_i)/(N/2)_m.                  (3)

This follows by Gaussian polar integration and the product Wick rule. Put L_k=lcm(1,...,k),with L_0=1. The following explicit,conservative integer clears EVERY denominator in degree2m:

    D_N(m)=4^m*L_(2m+N)^(N+2),
    D_N(m)*M_alpha belongs to Z.                           (4)

Here is a direct p-adic proof,including both parities of N. If N=2s,

    M_alpha=(s-1)!*prod(2alpha_i)! /
               (4^m*(m+s-1)!*prod alpha_i!).

At q=p^a,one has floor(2alpha_i/q)-floor(alpha_i/q)>=floor(alpha_i/q). Their sum is at least floor(m/q)-(N-1). Also floor((m+s-1)/q) is at most floor(m/q)+floor((s-1)/q)+1. Thus each Legendre summand for4^m*M_alpha is at least-N. Only q<=2m+N can contribute. This proves the claimed bound for even N.

For N=2s+1,

    M_alpha=(2s)!*(m+s)!*prod(2alpha_i)! /
                  (s!*(2m+2s)!*prod alpha_i!).

The same product estimate leaves,in each Legendre summand,

    floor(2s/q)-floor(s/q)+floor(m/q)+floor((m+s)/q)
                                     -floor((2m+2s)/q)-(N-1).

Use floor(2s/q)-floor(s/q)>=floor(s/q),and split the remaining expression into the two carry errors for m+s and2(m+s). Each error is at least-1. The summand is therefore at least-(N+1). Again only q<=2m+N contribute. Formula(4) follows,including N=1.

The denominator is at most exponential in m for fixed N. Indeed

    L_(2k) divides L_k*binomial(2k,k).

For each prime,the exponent increase in the LCM is at most1;the newly crossed prime power contributes a nonzero valuation to that binomial coefficient. Dyadic induction gives L_(2^j)<=4^(2^j-1),and hence the sufficient bound L_k<=16^k for all k>=1. This proves exponential denominator growth without a prime number theorem,period regularity theorem or a guessed recurrence.

## 4. A global homogeneous Gaussian growth-gap theorem

Let P!=0 be homogeneous of degree D>=1 in N>=1 standard Gaussian variables,with algebraic coefficients. Let Q be ANY polynomial with algebraic coefficients. Write

    A_n=E[Q P^n], k_n=ceil(Dn/2),
    R_n=2^(k_n)*(N/2)_(k_n), b_n=A_n/R_n.                  (5)

There are constants C>=1 and n0>=1,depending on P,Q,N,with

    b_n=0 OR C^(-n)<=|b_n|<=C^n  for every n>=n0.          (6)

Thus on every infinite nonzero subsequence,

    log|A_n|/(n log n) -> D/2.                              (7)

This is an all-dimension,all-degree statement,not a numerical asymptotic observation. Eventual-zero sequences are allowed and excluded from the limiting assertion. Algebraicity of coefficients is essential to the pointwise LOWER bound;no such bound is asserted for arbitrary Liouville coefficients.

Proof. Let epsilon_n be the parity of Dn. Only marker components Q_(epsilon_n+2j) contribute. Set J large enough for all these finitely many j. Polar integration gives

    b_n=sum_j 2^j*(N/2+k_n)_j *
               integral_sphere Q_(epsilon_n+2j)(u) P(u)^n d sigma. (8)

Every prefactor in(8) is the INTEGER product prod_(a=0..j-1)(N+2k_n+2a). All conjugates of b_n are bounded by a fixed exponential times a polynomial in n,using compactness of the sphere. They are consequently at most C_1^n after enlarging C_1.

Place all coefficients in one number field K. Choose an integer M>=1 making every coefficient of M*P and M*Q an algebraic integer. Formula(4),with m<=k_n+J,shows that

    T_n*b_n is an algebraic integer,
    T_n=M^(n+1)*4^(k_n+J)*L_(2(k_n+J)+N)^(N+2).            (9)

There is a fixed C_2 with T_n<=C_2^n for n>=1. If b_n!=0,its field norm is a nonzero rational integer after multiplying by T_n. Hence

    1<=T_n^[K:Q]*prod_(sigma:K->C)|sigma(b_n)|.

Use the exponential upper bounds for all other embeddings to obtain the lower bound in(6) at the chosen embedding. The upper bound is already proved. Finally log R_n=(D/2)n log n+O(n),by elementary factorial/Pochhammer estimates. This proves(7).

The normalized ordinary series sum b_n z^n has the sharp dichotomy

    it is a polynomial,OR it has finite positive radius.   (10)

In particular it is entire iff the original fixed marked tail is eventually zero. This criterion is not an implemented analytic-continuation algorithm,and it does not decide the all-marker condition from one arbitrary marker. The proof requires no new holonomic regularity assertion.

## 5. Exact factorial laws determine homogeneous degree,even over C

Let P!=0 be homogeneous of degree D in N Gaussian variables,now allowing ARBITRARY complex coefficients,and let Q be a fixed finite polynomial marker. Suppose

    E[Q P^n]=c^n*(ell*n+b)! for EVERY n>=1,
    c!=0, ell>=0,b>=0 integers.                            (11)

Then necessarily

    D=2ell.                                                (12)

In particular ell=0 is impossible for positive-degree P. More generally,for a fixed rational non-eventually-zero target sequence whose logarithmic growth order on its nonzero indices exists and is finite theta,the equality to c^n times that sequence forces D/2=theta. It is the exact algebraic specification of all moments that permits descent below;no pointwise lower estimate for arbitrary complex P is claimed.

Proof. With algebraic P,Q,c,this is immediate from(6)--(7),since the target is nonzero at all indices. For arbitrary complex coefficients,fix the finite degrees/support spaces for P and Q. In the polynomial ring over Q in their coefficients,c,and an inverse variable e,impose ec-1 and ALL polynomial moment identities in(11). Gaussian monomial moments are rational. A given complex solution makes this ideal proper. Its quotient is still a finite-type Q-algebra,even though the defining list was infinite. By the weak [Hilbert Nullstellensatz,Stacks00FV](https://stacks.math.columbia.edu/tag/00FV),a maximal ideal has finite residue field over Q. Thus an algebraic solution exists satisfying EVERY identity,with c still nonzero. P cannot become0 because its marked moments are nonzero. Its homogeneity forces degree exactly D. The algebraic theorem contradicts D!=2ell.

The rational-target version uses the same ideal,retaining at least one nonzero target identity and infinitely many nonzero target indices. No finite moment truncation is used to justify specialization;all equations are kept in the proper ideal. This is an existence argument,not an implemented finite-prefix search or an effective bound.

The degree condition is sharp as a possibility:in one circular pair E(ZW)=1,the homogeneous phase P=(ZW)^ell has degree2ell and E(P^n)=(ell*n)!. Adjoining independent Gaussian coordinates preserves the identity. This does not claim every homogeneous phase of that degree has such a law.

## 6. The actual exceptional Gaussian lift cannot carry the old quadratic radial law

The [actual radial bridge](GAUSSIAN_RADIAL_TRANSFER_AND_EXCEPTIONAL_UNIFORMITY.md) uses FOUR standard Gaussian coordinates and homogeneous right-Gamma-invariant lifts. Let S,T be nonzero such homogeneous invariant generators,with positive degrees d_S,d_T,and assume they are jointly pure as the cotangent marker architecture would require. A homogeneous invariant function restricting to a nonzero constant on SL2 is a scalar determinant power,and cannot be pure;restriction to0 would make it the zero polynomial. Thus each restriction is nonconstant. The retained exact invariant-degree theorem gives

    d_S,d_T>=8 for Gamma=2O,
    d_S,d_T>=12 for Gamma=2I.                              (13)

Even if the torus compatibility in Section2 is dropped,the literal quadratic marked radial law

    E[Q S^n T^n]=(2n)!/2^n, n>=1                          (14)

is impossible for EVERY fixed polynomial marker Q,including inhomogeneous ones. Apply(12) to P=ST:law(14) forces d_S+d_T=4,contradicting(13). Extra homogeneous marker components cannot change this degree constraint. More generally an exact factorial radial certificate of slope ell would require d_S+d_T=2ell;the old slope2 law cannot simply be relabelled as a quotient law.

Section2 is stronger for the direct torus-compatible,jointly pure invariant presentation:it forbids that presentation for ANY radial law,not only(14). Section5 supplies a separate global Gaussian obstruction that does not depend on the exceptional analytic foundation. Genuinely different presentations with noninvariant generators,non-torus gradings,nonhomogeneous generators,or non-jointly-pure cancellations remain separate,but the low-marker obstruction still forbids transferring the old quadratic witness into ANY pure actual invariant lift. Each surviving route needs its own actual invariant identity. Neither transfers an ambient Gaussian failure into the quotient.

The last missing implication for a full exceptional MZ theorem remains controlling the mixed-isotypic,image-interior phase with strict leading nullcone and opposite lower weights. Its pieces cannot be assumed jointly pure or replaced by an independent Gaussian law. This note closes the proposed DIRECT marker transfer and the global homogeneous growth/degree theorem,not that last full MZ implication.

## 7. Verification and next action

The [checker](verify_homogeneous_moment_marker_obstructions.py) verifies explicit sphere denominator and LCM identities at small parameters,Gaussian radial normalization in odd/even phase parity,and sharp exact homogeneous factorial examples. It also checks the algebraic weight-balance certificate for the joint-pure obstruction. The [receipt](homogeneous_moment_marker_verification.json) and [source audit](homogeneous_moment_marker_source_audit.json) bind these finite controls. The universal arithmetic,norm,Nullstellensatz and neutral-rigidity arguments are written proofs,not inferred from those controls.

No large prime campaign,degree search,generic invariant solver,external review or new exceptional failure is performed. Next pursue the genuinely non-jointly-pure mixed sector through its actual periods or a different complete global criterion. Retain the strict leading,finite cyclic stabilizer and image-interior conditions. All earlier field94,least intervals,closed cotangent classification,ordinary scaling ideals and radial bridge statements remain unchanged.
