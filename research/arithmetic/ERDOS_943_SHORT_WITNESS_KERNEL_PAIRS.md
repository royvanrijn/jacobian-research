# Short witnesses for compatible square-full kernels

The short-witness part of the balanced kernel-pair problem admits a pointwise power saving. The proof uses signed kernel graphs and a square-full gcd normalization that preserves the original root-height bound. A precise complement reduction isolates the remaining long-witness problem.

## 1. The truncated count

Fix positive constants c,C,C0. Let T>=2, let N be an integer with c T^5<=N<=C T^5, and let 1<=H<=C0 T^2. Let R_N(T;H) count ordered tuples (b1,b2,u,v) such that

* b1,b2 are distinct squarefree integers in [T,2T), with gcd(b1 b2,N)=1;
* 1<=u,v<=H, gcd(u,v)=1, and gcd(uv,N)=1;
* b1^3 u^2 == b2^3 v^2 modulo N.

Let K_N(T;H) count the distinct ordered kernel pairs admitting such a tuple.

**Theorem 1 (short witnesses).** For every epsilon>0,

    K_N(T;H) <= R_N(T;H)
      <<_(epsilon,c,C,C0) T^epsilon
          * min(T^2, H^2/T^(1/2) + H^3/T^2).                 (1.1)

The trivial term in the minimum for R includes a subpower factor, as each fixed kernel pair has at most a subpower number of reduced witnesses by rational reconstruction. In particular, for every fixed 0<=alpha<1/4,

    K_N(T;T^(1+alpha)) << T^(3/2+2alpha+epsilon).            (1.2)

For H below a sufficiently small constant times T there are no witnesses at all: the nonzero difference b1^3u^2-b2^3v^2 is divisible by N and has absolute value at most 8T^3H^2. The power-saving cutoff and its precise long complement are stated in Section 5.

The full count permits H=4T^2. Formula (1.1) then has no saving over T^2. The [full pair target](POWERFUL_SUM_GRAPHS.md#square-full-kernel-pairs-at-the-balanced-scale) and the [pointwise six-cycle target](ERDOS_943_SIX_CYCLES.md) remain open. This theorem does not improve the current all-target exponent 3/10.

## 2. Signed kernel graphs with arbitrary vertex sets

We first prove the signed version of the existing kernel-graph argument. Let kappa1,kappa2,n be fixed positive integers, n>=2. Count coprime positive terms in

    kappa1 a^2 b^3 - kappa2 c^2 d^3 = n,                    (2.1)

where b,d are squarefree and range over specified finite sets B,D of cardinality O(T), and 1<=a,c<=H. Assume all the data have polynomial size in T, n>=c T^5 and H<=C0 T^2. The sets B,D need not be intervals; their numerical members can be larger than T.

**Lemma 2 (signed finite-kernel graph).** Under the hypotheses above, the simple edge count and the representation-tuple count in (2.1) are both bounded by

    O(T^(3/2+epsilon) * (1+H^2/sqrt(n))).

**Proof.** The two terms being coprime implies that every factor on either side is a unit modulo n, and gcd(a,c)=1. For each fixed b,d, the reduced fraction a/c is a unit solution of

    (a/c)^2 == kappa2 d^3/(kappa1 b^3) modulo n.

There are at most 2^(omega(n)+1) such residue classes. Two reduced positive fractions of height <=H in one class have cross-difference divisible by n. Since H^2<n for all sufficiently large T under the stated assumptions, there is at most one such fraction per class. Each fraction is already the coprime pair (a,c). Thus each simple kernel edge has O(T^epsilon) representations. For bounded T the constants absorb the remaining finite range. This argument does not use a ternary conic theorem or imaginary-quadratic finiteness.

For two distinct b1,b2 sharing d, with chosen representation labels (a_i,c_i), eliminate the negative term to obtain

    kappa1[b1^3(a1 c2)^2-b2^3(a2 c1)^2]
       = n(c2^2-c1^2).                                    (2.2)

Put a_i=g alpha_i, c_i=h gamma_i, with gcd(alpha1,alpha2)=gcd(gamma1,gamma2)=1. Coprimality within each representation gives gcd(alpha_i,gamma_i)=1. Hence the fraction

    u0/v0 = a1 c2/(a2 c1) = alpha1 gamma2/(alpha2 gamma1)

is reduced as displayed. It has height <=H^2 and belongs to a modular square-root class of b2^3/b1^3. The angular packing lemma gives

    O(T^epsilon * (1+H^4/n))                              (2.3)

possible fractions: distinct positive fractions of height Z in the same residue have angular separation at least n/(2Z^2), so there are O(1+Z^2/n), and CRT contributes only a subpower number of square-root classes.

For each fraction, its four reduced factors have divisor-many choices. They give

    kappa1 g^2(b1^3 u0^2-b2^3 v0^2)
       = n(gamma2^2-gamma1^2).                             (2.4)

The coefficient of g^2 is nonzero: equality of b1^3 u0^2 and b2^3 v0^2 forces b1=b2 by parity of prime valuations. Thus g has at most one positive value. Substitution gives

    h^2 d^3 = [kappa1 b1^3 g^2 alpha1^2-n]
                 /(kappa2 gamma1^2),                      (2.5)

which has at most one positive square-times-squarefree-cube decomposition. This reconstructs both representations. The codegree bound (2.3) follows. Interchanging the terms, allowing the sign of the fixed right-hand side to change, proves the same bound on the other side. Positivity is used to retain positive roots, not to require a plus sign in (2.1).

For a bipartite graph with both sides O(T) and maximum codegree Q, summing binomial(deg(v),2) and using Cauchy--Schwarz gives e<<T+sqrt(Q)T^(3/2). Therefore both the number of simple edges and the number of representation tuples in (2.1) are

    << T^(3/2+epsilon) * (1+H^2/sqrt(n)).                  (2.6)

For tuples we use the fixed-edge argument above. For distinct kernel pairs the simple-edge bound alone suffices. All divisor factors in this proof have polynomial-size arguments, so a smaller intermediate epsilon absorbs their product.

## 3. Square-full gcd normalization preserves root height

**Lemma 3 (height and kernel recovery).** Every witness can be normalized as below with n>=N, square roots no larger than the original roots, and transformed squarefree kernel sets of cardinality O(T). The normalized tuple and the finite normalization data recover the original tuple.

**Proof.** A witness in Section 1 has a nonzero integer shift

    b1^3 u^2-b2^3 v^2 = N h,
    0<|h|<=J:=floor(8T^3 H^2/N).                          (3.1)

The case h=0 would imply b1=b2. Put m1=b1^3 u^2, m2=b2^3 v^2 and d0=gcd(m1,m2). Since both terms are units modulo N, gcd(d0,N)=1; consequently d0 divides |h|. Moreover d0 is square-full: every positive prime valuation in either m_i is at least 2, hence so is their minimum.

Extract kappa_i as the product of primes occurring to exponent exactly 1 in m_i/d0. Each kappa_i is a squarefree divisor of d0. The remaining square-full part has a unique expression

    m1/d0=kappa1 a^2 b^3,
    m2/d0=kappa2 c^2 d^3,                                 (3.2)

with b,d squarefree. These two terms are coprime. After swapping the sides if h<0, their difference is n=N|h|/d0>=N.

The essential point is

    a<=u,   c<=v.                                         (3.3)

Here is the complete prime-by-prime check. At one prime the original valuation is 2e+3epsilon, epsilon in {0,1}. Let t be the valuation removed in d0. Either t=0 or t>=2. Let delta in {0,1} be its valuation in kappa, and let epsilon' in {0,1} be its valuation in the new squarefree cube kernel. If e' is its new root exponent, then

    2(e'-e)=3epsilon-t-delta-3epsilon'.                    (3.4)

For t=0 the canonical decomposition is unchanged. For t>=2 the right side is at most 1 and is even, so it is at most 0. This proves e'<=e. The square-full condition on d0 is necessary for this proof: dividing 2^3 by the non-square-full divisor 2 produces 2^2 and increases the canonical root from 1 to 2. Such a divisor cannot be the gcd of two square-full numbers.

For fixed d0,kappa_i, the transformed kernel set still has O(T) elements. Write sf(m) for the product of primes of odd valuation in m, and s_i=sf(d0 kappa_i). Parity in (3.2) gives

    b=sf(b1 d0 kappa1)
      =b1 s1/gcd(b1,s1)^2,                               (3.5)

and the analogous formula for d. Multiplication of square classes by s_i is an involution on squarefree integers. Thus (3.5) is injective on the original O(T) kernels, even though the new numerical kernels need not remain in [T,2T). This is why Section 2 was stated for arbitrary finite vertex sets.

Given the normalization data and the normalized tuple, the original m1,m2 and their unique canonical decompositions recover the original witness. Thus the map is injective. There are at most tau(d0)^2 choices for kappa1,kappa2. All normalized terms are at most 8T^3 H^2/d0, so the polynomial-size hypothesis in Section 2 holds. No coefficient bound of the form kappa_i b^3<=n is needed.

## 4. Summing the short shift

Apply (2.6) to (3.2) for each nonzero h, square-full d0 dividing |h|, and the possible kappa_i. Up to T^epsilon factors, this gives

    R_N(T;H) << T^(3/2+epsilon)
      * sum_(1<=j<=J) sum_(d0|j, d0 square-full)
           [1+(H^2/sqrt(N))*sqrt(d0/j)].                  (4.1)

The two signs and coefficient divisor counts are absorbed in the constant and subpower factor. If J=0 the original count is zero.

**Lemma 4 (square-full divisor sums).** For J>=1 the two divisor sums in (4.2) are respectively O(J) and O(sqrt(J)log(2J)).

**Proof.** Square-full integers have a unique expression d0=A^2 B^3 with B squarefree. Consequently

    sum_(d0 square-full) 1/d0 < infinity,
    sum_(d0<=J, d0 square-full) 1/sqrt(d0) << 1+log(2J).

These follow by summing A^-2 B^-3 and A^-1 B^(-3/2), respectively. Reordering the divisor sums therefore yields

    sum_(j<=J) sum_(d0|j, d0 square-full) 1 << J,
    sum_(j<=J) sum_(d0|j, d0 square-full) sqrt(d0/j)
      =sum_(d0 square-full<=J) sum_(ell<=J/d0) ell^(-1/2)
      << sqrt(J) log(2J).                                 (4.2)

Thus

    R_N(T;H) << T^(3/2+epsilon)
          * [J+(H^2/sqrt(N))*sqrt(J)].                    (4.3)

Since J<<H^2/T^2, the terms in (4.3) are at most H^2/T^(1/2) and H^3/T^2. This proves the nontrivial part of (1.1). For a fixed original kernel pair, the same angular/CRT argument at modulus N and fraction height H gives O(T^epsilon) witnesses, since H^2<N for H<=C0T^2 and large T. There are O(T^2) pairs. This proves the other term of (1.1).

The saving is therefore a frequency bound across varying kernel pairs. It does not result from counting O(1) completions for one fixed pair and then assuming that most pairs have no completion.

## 5. A power-saving cutoff and its exact complement

**Corollary 5 (short-cutoff saving).** Fix delta in (0,1/4], and put H0=T^(5/4-delta). Uniformly in the range of Section 1,

    K_N(T;H0) <<_(epsilon,delta,c,C) T^(2-2delta+epsilon).    (5.1)

Indeed, the two terms in (1.1) have exponents 2-2delta and 7/4-3delta, and the second is no larger than the first. The same estimate holds for reduced witness tuples. The cutoff exponent is a sufficient range for a saving; no necessity or sharpness of that exponent is asserted.

Write K_full=K_N(T;4T^2) and K_short=K_N(T;H0). Define K_long to count distinct ordered kernel pairs that have at least one reduced unit witness of height at most 4T^2 but have NO such witness of height at most H0. Thus a pair with both a short and a long witness is counted only in K_short. The word long refers to this complement, or equivalently to the minimum admissible witness height, and not to all tuples having a long witness.

**Proposition 6 (exact complement partition).** For each admissible N,T and the fixed delta,

    K_full = K_short + K_long.                            (5.2)

The short-witness pairs form a subset of the full pairs. Every full pair either belongs to this subset or has no short witness, exclusively. This proves (5.2).

**Theorem 7 (both quantitative directions).** Interpret every power-saving bound here uniformly over the integers N and real T in Section 1, for every epsilon>0, with the constants allowed to depend on the stated fixed parameters.

1. If, for a fixed eta>0,

       K_long <<_epsilon T^(2-eta+epsilon),

   then

       K_full <<_epsilon T^(2-min(eta,2delta)+epsilon).      (5.3)

2. If, for a fixed eta>0,

       K_full <<_epsilon T^(2-eta+epsilon),

   then the SAME estimate holds for K_long.                (5.4)

The first direction follows by adding (5.1) and the assumed complement bound in (5.2). The second follows from 0<=K_long<=K_full. These are proved conditional implication theorems; neither displayed antecedent is claimed to be established.

**Corollary 8 (positive-saving equivalence).** For each one fixed delta in (0,1/4], a uniform positive power saving for K_full exists if and only if a uniform positive power saving for its precise complement K_long exists. The exponent can decrease from eta to min(eta,2delta) in the complement-to-full direction. No statement that every long-witness tuple must have a power-saving count is inferred.

## 6. The remaining pointwise obligation

For the concrete fixed cutoff delta=1/8, the new complement question is:

**Open question.** Is there an absolute eta>0 such that the number of distinct ordered squarefree kernel pairs with a reduced unit witness of height at most 4T^2, but with NO such witness of height at most T^(9/8), is O_epsilon(T^(2-eta+epsilon)) whenever T is comparable to N^(1/5)?

By Theorem 7, a proof would give the full pair-frequency saving with exponent min(eta,1/4). Conversely the full saving would bound this complement. The equivalence does not establish either bound, and the existing full-pair and six-cycle obligations remain open.

The signed-layer proof explains its present limit. At H of order T^2 there can be O(T^2) distinct shifts. Even layers carrying only one new kernel edge each could cover O(T^2) pairs while satisfying all the within-layer graph bounds used above. The original modular reconstruction also bounds each pair's TOTAL number of witnesses across all shifts by a subpower number. That bounds multiplicity rather than the number of nonempty pairs. Additional arithmetic control across the long-height layers is still needed.

## 7. Mathematical inputs and scope

The only reused mathematical ingredients are the elementary unit square-root/CRT reconstruction in [the kernel-graph note](ERDOS_943_KERNEL_GRAPH.md) and the angular packing lemma in [the dual-graph note](ERDOS_943_DUAL_KERNEL_GRAPHS.md). Sections 2--4 give the signed equation, reconstruction, height preservation and divisor summation in full. Arbitrary finite transformed kernel sets are allowed; their members are not silently assumed to remain in [T,2T).

No ternary conic, Pell, Thue, determinant-method or graph supersaturation theorem is needed for Theorem 1 or the complement equivalence. In particular the [retained Trifonov completion estimate](TRIFONOV_2026_TERNARY_CONIC_SOURCE_SCOPE.md) is not a mathematical premise here. Counting completions for one fixed kernel pair does not by itself prove a saving in the frequency of compatible pairs. No finite enumeration is used to prove the unrestricted counting result.
