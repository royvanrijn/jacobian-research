# Trifonov's ternary conic estimate and the square-full kernel targets

## 1. The September 2026 source and its exact scope

Ognian Trifonov's [arXiv:2609.17715v1](https://arxiv.org/abs/2609.17715v1), dated15September2026, is *Lattice points close to a smooth curve and squarefull numbers in short intervals. II*. Its Theorem1 gives a short-interval asymptotic for individual square-full counts at theta=8/65+epsilon. Its Theorem2 concerns integer points close to a C^3 graph, with derivative and thin-strip hypotheses. These are context, not a fixed-sum representation bound. [Actual introduction and Theorems1--2](https://arxiv.org/html/2609.17715v1#S1).

The mathematical input below is Lemma10's coefficient-uniform ternary conic box estimate. It is credited to Heath-Brown's2024 communication and proved publicly in Section4.1 via published Browning--Heath-Brown lemmas. The latter paper's Lemma2.6 already contains the squarefull-minor refinement. [Trifonov, Section4.1](https://arxiv.org/html/2609.17715v1#S4.SS1), [Browning--Heath-Brown, exact2018 primary](https://arxiv.org/pdf/1801.00979v2).

## 2. Integral-matrix convention and the uniform conic theorem

For nonzero m let m# be the product of the WHOLE prime powers p^e dividing m exactly with e>=2. This is the positive squarefull part, not generally the largest square dividing m.

Let M be an INTEGER SYMMETRIC nonsingular matrix, q(X)=X^t M X, Delta=abs(det M), and Delta0 the positive gcd of its2by2minors. This literal convention avoids half-integral off-diagonal ambiguity. An application using another polynomial convention must clear and track denominators first. The diagonal specialization below uses precisely its displayed integer matrix, including at2.

**Theorem 1 (classical coefficient-uniform conic box bound).** For R1,R2,R3>0, the number of primitive integer triples with q(X)=0 and abs(Xi)<=Ri is

    O_epsilon(d3(Delta)
      *[1+(R1 R2 R3 (Delta0#)^(3/2)/Delta)^(1/3)]).           (2.1)

This is Trifonov Lemma10. Here is the complete coordinate interpretation of its proof. Browning--Heath-Brown Lemma2.6 covers all primitive zeros by O(C(q)) full-rank integral lattices Lambda satisfying

    det Lambda >> Delta/(Delta0#)^(3/2),   C(q)<=d3(Delta).

Put D_R=diag(R1,R2,R3) and Lambda*=D_R^-1 Lambda. Lemma2.3 supplies a basis g_j with product_j norm(g_j)>=det Lambda* and abs(c_j)<<norm(x*)/norm(g_j) when x*=sum c_j g_j. The unit cube has norm bounded by sqrt3. Restore the ORIGINAL integer basis h_j=D_R g_j in Lambda and set

    q1(c)=q(sum_j c_j h_j).

This is integral and nonsingular. A primitive original point has primitive integer coefficient vector c, because a common divisor of c divides the original point. Lemma2.5 bounds these c in their coefficient box by O(1+(product_j norm(g_j))^(-1/3)). The determinant estimate and number of lattices prove(2.1). The finer C(q) can vanish through local nonsolubility; replacing it by d3(Delta) is only an upper bound.

Three source-transcription details need this precise reading. Trifonov reuses q after scaling; the integral form is q1 on the original integer basis as above. Its printed Lemma13 omits the box inequalities explicit in the original Lemma2.5; the application does use that box. The final basis-product citation is to the basis lemma (Lemma12 there), though labeled Lemma11. None licenses an unbounded count of primitive conic points.

The actual proof of Browning--Heath-Brown Lemma2.6 gives the beta=0 and beta=1 local improvements. Higher beta and dyadic lattice estimates are imported from the cited earlier theorem. Lemma2.3 uses the stated Cassels minima results; Lemma2.5 imports the earlier conic-box result. Those earlier cited proofs are not claimed newly audited here.

## 3. An exact specialization to compatible kernel pairs

Let N>=2 and let b1!=b2 be positive squarefree integers with gcd(b1b2,N)=1. For0<u,v<=H consider

    b1^3 u^2-b2^3 v^2=N h.                               (3.1)

Compatibility additionally asks for a reduced unit fraction u/v; those restrictions may be dropped for an upper bound. The case h=0 is impossible, because prime valuations modulo2 would force b1=b2.

Put g=gcd(b1,b2), X=b1u, Y=b2v and W=1. Then

    q_h(X,Y,W)=b1X^2-b2Y^2-NhW^2=0,                      (3.2)
    M_h=diag(b1,-b2,-Nh),
    Delta=b1b2N abs(h),
    Delta0=gcd(b1b2,b1Nh,b2Nh).                           (3.3)

The point(X,Y,1) is AUTOMATICALLY primitive, even if(u,v) is not. Its box has R1=b1H,R2=b2H,R3=1. No content division or dyadic correction is made.

A prime dividing exactly one kernel occurs in Delta0 to exponent at most1, so contributes nothing to Delta0#. For p dividing both kernels, squarefreeness and gcd(b1b2,N)=1 give

    v_p(Delta0)=min(2,1+v_p(h)).

A prime dividing neither contributes nothing. Therefore, INCLUDING p=2,

    Delta0#=gcd(g,h)^2,
    (Delta0#)^(3/2)=gcd(g,h)^3.                           (3.4)

**Corollary 2 (fixed-shift kernel-pair bound).** For fixed nonzero h, the number R_h of positive integer pairs(u,v) in(3.1) with u,v<=H is

    O_epsilon(d3(b1b2N abs(h))
      *[1+gcd(g,h)*(H^2/(N abs(h)))^(1/3)]).               (3.5)

The injective map(u,v)->(b1u,b2v,1), the stated box and(3.3)--(3.4) prove this by(2.1). The factors b1b2 cancel from its nonconstant term. Extra divisibilities, reducedness and unit conditions only restrict the count.

In a coprime representation box kappa1 a^2b^3+kappa2 c^2d^3=N, every factor of either summand, including the fixed common kappa1,kappa2, is a unit modulo N. Common-neighbor elimination gives

    kappa1(b1^3u^2-b2^3v^2)=N(c2^2-c1^2).

Cancellation modulo N is valid because gcd(kappa1,N)=1. Thus(3.1) holds with INTEGER h and the additional identity kappa1 h=c2^2-c1^2 only restricts witnesses. No hidden multiplier remains in matrix(3.2). Unequal coefficients attached separately to b1,b2, or a nonunit common coefficient, require their own normalization; formula(3.4) is not asserted in those settings. Unrestricted targets retain the existing separate coprime/divisor normalization.

## 4. Shift range, literal summation and its limitation

On the balanced scale b1,b2 in[T,2T), N comparable to T^5 and H=O(T^2),

    0<abs(h)<=8T^3 H^2/N=O(T^2).                         (4.1)

The determinant in(3.3) has polynomial size in T, so its divisor factor is O_epsilon(T^epsilon). For coprime kernels g=1 the nonconstant term in(3.5) is O(T^(-1/3)); the ADDITIVE1 remains. This counts completions for fixed kernels and a shift, not how often a kernel pair is nonempty.

For a positive integer J on the relevant polynomial-height range, literal summation gives

    sum_(0<abs(h)<=J) R_h <<_epsilon T^epsilon
      *[J+(H^2/N)^(1/3)J^(2/3)tau(g)].                    (4.2)

Indeed gcd(g,h)=sum_(d|g,d|h)phi(d), so summing h^(-1/3) gives O(J^(2/3)sum_(d|g)phi(d)/d)<=O(J^(2/3)tau(g)). When H is comparable to T^2, N to T^5 and J=O(T^2), the two displayed costs are at most T^2 and T^(1+epsilon). Large discriminant does not remove the first term.

The current rational-reconstruction lemma already bounds ALL shifts for a FIXED pair more strongly: O_epsilon(N^epsilon(1+H^2/N)) reduced unit fractions realize its modular square roots, which is O_epsilon(T^epsilon) here. Thus(4.2) does not improve that completion count. The remaining target is a saving in the number of DISTINCT compatible pairs. This displayed conic specialization and summation do not by themselves provide it; this is not a claim that no other use of the method can work.

## 5. Current targets and the near-curve comparison

The [dual-graph proof](ERDOS_943_DUAL_KERNEL_GRAPHS.md) proves exponent3/10 and has a unique critical logarithmic point A=B=C=D=N^(1/5). The [pair target](POWERFUL_SUM_GRAPHS.md#square-full-kernel-pairs-at-the-balanced-scale) is K_N(T)<<T^(2-eta+epsilon),0<eta<=1; its current mean is O(T^(1+epsilon)), and infinitely many targets attain orderT. The [six-cycle target](ERDOS_943_SIX_CYCLES.md) is O(T^(2+epsilon)) pointwise, compared with current pointwise O(T^(3+epsilon)) and mean O(T^(1/2+epsilon)). These are distinct quantities.

In Trifonov Theorem2 use M for its curve-domain parameter and T_curve for its height parameter. On a compact balanced fixed(a,d) mixed-cubic branch c=f(b), both parameters are comparable to Q=N^(1/5), delta=0, with abs(f'') comparable to Q^-1 and abs(f''') comparable to Q^-2. Its leading terms then have exponents16/25 and9/14; the latter dominates. This O(Q^(9/14+epsilon)) estimate is weaker there than the existing exact affine-cubic Bombieri--Pila O(Q^(1/3+epsilon)). The improved single-squarefull short-interval exponent8/65 is not a fixed-sum representation exponent.

The conic estimate also does not count the simultaneous Cramer perfect-cube tests for varying six-cycle root labels. Another reduction or averaging estimate would have to control those labels and shifts with the cube, squarefree and coherent modular-color conditions. No such transfer is established here, and neither current gap is closed by this source addition.

## 6. Primary scope and attribution

[T] Ognian Trifonov, *Lattice points close to a smooth curve and squarefull numbers in short intervals. II*, arXiv:2609.17715v1,15September2026. [Retained primary](sources/trifonov-lattice-curve-squarefull-II-2609.17715v1.pdf), SHA2569348f6453a8c0b44a5cea7b429c29c89c6b106bcc7c57f53d5b11b355caf9c38. Read: introduction/Theorems1--2 and WHOLE Section4.1, Lemmas10--13 and the complete Lemma10 proof, printed7--9, including visual inspection. Theorem1/2 are comparison statements; their complete proof sections are not claimed freshly audited.

[BH] T.D.Browning and D.R.Heath-Brown, *Counting Rational Points on Quadric Surfaces*, Discrete Analysis2018:15,29pp., published7September2018; exact [arXiv:1801.00979v2](https://arxiv.org/pdf/1801.00979v2). [Retained primary](sources/browning-heath-brown-quadric-surfaces-1801.00979v2.pdf), SHA256edad12f63a7c07309c113d76a54e9953be753fb08fbe82cd287bab48638fb511. Read: Lemma2.3 statement/full proof, Lemma2.5 statement and its proof reference/comment, and Lemma2.6 statement/ENTIRE displayed proof, printed6--9; printed7--9 visually inspected. Cited earlier minima/conic/p-adic proofs remain imports. No whole29-page proof audit is claimed.

Sections3--4 record an elementary specialization of these classical uniform estimates, including the exact prime-adic computation(3.4). No new pointwise exponent or historical-priority claim for square-full sums follows. Existing arithmetic principals, gap states and assurance flags are unchanged.
