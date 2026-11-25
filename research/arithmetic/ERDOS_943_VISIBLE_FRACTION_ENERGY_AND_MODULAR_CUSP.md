# Primitive modular energy, the rational cusp and the exact long-complement reduction


The full balanced kernel-pair target admits an exact visible-fraction energy formula. The large raw diagonal collapses after primitive normalization, and all characters of order dividing6 are already below the target scale. This identifies a signed high-order covariance as an equivalent target, while keeping the known short-height saving separate.

**Theorem 1 (primitive energy and cusp identities).** Under the explicit size hypotheses in Sections1-6, equations(2)-(8) below hold for all integer moduli, with the prime-modulus exactness stated in Section6. They are identities, multiplicity bounds and norm bounds; none is a new full-height power saving.

**Theorem 2 (same-exponent spectral reduction and exact long complement).** Fix positive constants c,C, and let c*T^5<=N<=C*T^5. Use all squarefree unit kernels in [T,2T), the full height H_full=floor(4T^2), and the signed covariance W defined in Section3. For every fixed0<eta<=1, the uniform bound K_full<=C_epsilon*T^(2-eta+epsilon) is equivalent to the SAME one-sided upper bound for W. For each fixed0<delta<=1/4, define K_long as the pairs with some full-height witness but NO witness of height at most floor(T^(5/4-delta)). Existence of any uniform positive power saving is equivalent for K_full, K_long and W. Section7 proves both directions and the precise exponent loss min(eta,2delta) from the long complement. These are conditional equivalence theorems, not proofs of any of their power-saving predicates.

## 1. Exact setup

Let N>=2, let B be any finite set of positive squarefree integers in [T,2T), all coprime to N, and put M=|B|. Assume 2T<N, and let H be a positive integer with H^2<N. The balanced application has N comparable to T^5 and H<=C*T^2, with fixed positive constants, so all the further size conditions below hold for sufficiently large T.

Let E(h) count ordered tuples (b1,b2,u,v) with b1,b2 in B, 1<=u,v<=h, gcd(uv,N)=1 and

    b1^3 u^2 = b2^3 v^2 mod N.                               (1)

Let R(H) count the same tuples with gcd(u,v)=1, and let R_off(H) additionally require b1!=b2. Let K(H) count distinct ordered pairs b1!=b2 for which a primitive witness exists. These are counts of witnesses and of pairs respectively; they are not interchangeable without the multiplicity argument below.

For Dirichlet characters modulo N, extended by zero off the units, define

    B3(chi)=sum_(b in B) chi(b)^3,
    U2(chi;h)=sum_(1<=u<=h) chi(u)^2,
    C_H(chi)=sum_(u,v<=H, gcd(u,v)=1, gcd(uv,N)=1) chi(u)^2 conjugate(chi(v)^2).

Write Avg for phi(N)^(-1) times the sum over all these characters. All identities below are finite exact identities.

## 2. Removing the diagonal before estimating

Character orthogonality and Mobius inversion give

    E(h)=Avg |B3(chi)|^2 |U2(chi;h)|^2,
    R(H)=sum_(d<=H, gcd(d,N)=1) mu(d) E(floor(H/d)),
    C_H(chi)=sum_(d<=H, gcd(d,N)=1) mu(d) |U2(chi;floor(H/d))|^2.   (2)

The inversion is valid because u,v are units, their common divisor d is a unit, and its square cancels in (1).

For b1=b2, equation (1) is N | u^2-v^2. Since |u^2-v^2|<H^2<N, it forces u=v. In the primitive count, this leaves only u=v=1. Thus the raw diagonal is M*#{u<=h:gcd(u,N)=1}, while the primitive diagonal is exactly M, and

    R_off(H)=sum_d mu(d)[E(floor(H/d))-M*#{u<=H/d:gcd(u,N)=1}].   (3)

Also Avg C_H=1: its orthogonality condition is again N | u^2-v^2, and primitivity leaves (1,1). Therefore

    R_off(H)=Avg (|B3(chi)|^2-M) C_H(chi).                    (4)

The kernel C_H is real, because exchanging u and v conjugates each summand. It need not be nonnegative: for N=13,H=2 and a character of order4 with chi(2)=i, the three primitive pairs give C_H=1-1-1=-1. Replacing its Mobius sum or the covariance in (4) by absolute values discards essential cancellation.

If (2T)^3<N, then Avg |B3|^2=M as well, since the congruence b1^3=b2^3 mod N is an equality of integers. Formula (4) is then literally a covariance of a centered kernel weight and a mean-one visible-fraction kernel.

## 3. Low-order characters are already below the target scale

Let X6={chi:chi^6=1}. The unit-group structure at prime powers gives |X6|<=6^omega(N): an odd prime power contributes at most 6, and the 2-primary factor at most 4. Since

    |B3|<=M,  |C_H|<=H^2,

we have the elementary bound

    |Avg_(chi in X6) (|B3|^2-M) C_H|
       <= 6^omega(N) (M^2+M) H^2 / phi(N).                  (5)

At N comparable to T^5, H<=C*T^2 and M=O(T), this is O_epsilon(T^(1+epsilon)). In particular all real/quadratic and all cubic characters belong to this small contribution. They cannot alone account for a failure of a T^(2-eta) target with 0<eta<=1.

Thus, setting

    W(N;B,H)=Avg_(chi^6!=1) (|B3(chi)|^2-M) C_H(chi),

we have R_off(H)=W+O_epsilon(T^(1+epsilon)). For a given 0<eta<=1, an upper bound R_off(H)<<T^(2-eta+epsilon) is equivalent, with the usual subpower reallocation, to the same UPPER bound on this signed high-order covariance. This proves no such saving. Both chi^2 and chi^3 are nonprincipal in W. A method controlling only real characters at the fixed modulus does not address the remaining term.

## 4. Why the elementary moment route stops at TH

Assume also (2T)^2<N in this section. When h^2<N, ordinary fourth character moments over u<=h count literal equal products u1*u2=u3*u4, because both products are less than N. Divisor bounds give

    Avg |sum_(u<=h) chi(u)|^4 <<_epsilon h^2 N^epsilon.

The same holds for B, using its subset of an interval of length O(T). The power maps on the character group have at most N^epsilon fibres. Consequently

    || |B3|^2-M ||_2 << T N^epsilon,
    ||C_H||_2 <= sum_(d<=H) || |U2(chi;H/d)|^2 ||_2
              << H N^epsilon log(2H).

Cauchy-Schwarz in (4) therefore gives only TH*N^epsilon, or T^(3+epsilon) at the full height. No deep multiplicative-energy input was needed for this particular bound. The raw energy also has a diagonal of order TH when M is comparable to T and the modulus is prime, so bounding the raw energy by the desired T^(2-eta) would be false.

The norm scales themselves can be sharp. Under (2T)^3<N, trivial permutations of two elements of B show

    Avg |B3|^4 >=2M^2-M,
    || |B3|^2-M ||_2^2 >=M^2-M.

For prime N=p with 2H^2<p, the squaring map is injective on positive primitive fractions u/v of height at most H modulo p: if two squares agree, one of their cross-difference or cross-sum is divisible by p, and both have size less than p. Positivity excludes the sum. Thus

    ||C_H-1||_2^2 = #{(u,v)<=H : gcd(u,v)=1}-1,

which is comparable to H^2. For example, a union bound over possible common prime divisors leaves at least (2-zeta(2))*H^2 coprime pairs. Hence individual L2 norm estimates cannot supply the missing factor T^(-1). A correlation estimate exploiting their relation is required; this is not an impossibility result for other methods.

The displayed fourth-moment estimate in this height range is elementary: its congruence is an equality of integer products. It supplies no diagonal-free covariance estimate.

## 5. A toric parameter with no joint root-colour ambiguity

In any commutative ring, on the locus of units x,y satisfying y^2=x^3, there is exactly one unit z with

    x=z^2, y=z^3; namely z=y/x.

Indeed z^2=y^2/x^2=x and z^3=y^3/x^3=y. This is valid at primes 2 and 3 too. Applied to

    x=b1/b2 mod N,  y=v/u mod N,

equation (1) becomes this unit cusp, with

    z=b2*v/(b1*u) mod N.                                    (6)

Let F_H be the set of positive primitive fractions v/u, u,v<=H and unit modulo N. Their residue map is injective, because a nonzero cross-difference has absolute value less than H^2<N. Let w_B(x) count ordered kernel pairs with residue b1/b2=x. Then

    R_off(H)=sum_(z in (Z/N)^*, z!=1) w_B(z^2) 1_(F_H)(z^3).  (7)

The z=1 term is exactly M. No other z^2=1 contributes, because then y=z has y^2=1, which forces its positive primitive representative to be 1 by |u^2-v^2|<N.

When x,y are regarded as positive rationals rather than residue classes, the equality y^2=x^3 has no off-diagonal kernel points. If g=gcd(b1,b2), the numerator and denominator of x=(b1/g)/(b2/g) are coprime squarefree integers. A rational square has even valuations, so x^3 being a square forces both to be 1. Hence b1=b2, and the primitive positive y is 1. The target is genuinely modular rational-height cusp points away from all true rational cusp points; it is not an integer-coordinate small-box theorem or a count of rational points ON the curve over Q.

## 6. Exact pair/witness equality for prime moduli

If N=p is prime and 2H^2<p, two witnesses u1/v1,u2/v2 for one kernel pair satisfy

    p | (u1*v2-u2*v1)(u1*v2+u2*v1).

The difference has absolute value below p; the positive sum is below p. Therefore the difference is zero. Reducedness forces the witnesses to coincide. Consequently

    K(H)=R_off(H)                                             (8)

for prime moduli in this range. For general N, the existing CRT square-root bound instead gives K<=R_off<=O_epsilon(N^epsilon)K. Formula (8) is an exact prime simplification, not a pointwise saving.

## 7. Both conditional directions and the precise long complement

Fix c,C>0. All following bounds are uniform over N,T with c*T^5<=N<=C*T^5, and constants may depend on c,C and the fixed parameters stated. Finite initial ranges of T can be absorbed into those constants. Let B contain ALL squarefree integers in [T,2T) coprime to N, and put H_full=floor(4T^2). This is exactly the full-height pair quantity K_N(T) in the [balanced kernel-pair target](POWERFUL_SUM_GRAPHS.md#square-full-kernel-pairs-at-the-balanced-scale). Write K_full=K(H_full), R_full=R_off(H_full), and W_full=W(N;B,H_full).

The existing unit-root and rational-reconstruction argument gives, for every epsilon>0,

    K_full <= R_full <<_epsilon T^epsilon K_full.              (9)

Indeed each fixed pair's primitive fraction is uniquely determined by its unit residue because H_full^2<N for largeT, and its residue solves one unit quadratic congruence moduloN. CRT gives only O_epsilon(N^epsilon) roots. Changing epsilon absorbs N comparable to T^5. This is a multiplicity statement, not a new frequency bound.

By (4)-(5), independently of the short-height theorem,

    R_full = W_full + O_epsilon(T^(1+epsilon)).               (10)

For a real quantity A(N,T), define U_eta(A) to mean the ONE-SIDED uniform assertion

    for every epsilon>0 there exists C_epsilon such that
    A(N,T) <= C_epsilon T^(2-eta+epsilon).                    (11)

This convention is deliberate because W_full is signed. For K_full and R_full, which are nonnegative, it is the usual upper bound.

**Proposition 3 (same exponent in both directions).** For every fixed0<eta<=1,

    U_eta(K_full)  IFF  U_eta(R_full)  IFF  U_eta(W_full).      (12)

Proof. The first equivalence is (9), reallocating the arbitrary epsilon between its subpower factor and (11). If U_eta(R_full) holds, (10) bounds W_full from above by R_full plus O(T^(1+epsilon)). Conversely the same identity bounds R_full by W_full plus that error. Since eta<=1, the error is at most the scale in(11), again with epsilon reallocated. Neither direction requires W_full itself to be nonnegative. This does not say W_full is nonnegative or bounded by its absolute-moment estimate at the desired scale.

Now fix0<delta<=1/4, set H0=floor(T^(5/4-delta)), and define

    K_short=K(H0),
    K_long=#{ordered kernel pairs with SOME witness<=H_full
                                  and NO witness<=H0}.

A pair with both short and long witnesses belongs only to K_short. The definitions therefore give the exact disjoint partition

    K_full=K_short+K_long.                                   (13)

The established [short-witness theorem](ERDOS_943_SHORT_WITNESS_KERNEL_PAIRS.md), Corollary5, gives

    K_short <<_epsilon T^(2-2delta+epsilon).                  (14)

No statement counting all long witnesses or all pairs having a long witness is substituted for K_long.

**Proposition 4 (both quantitative long-complement directions).** For every fixed0<eta<=1:

* U_eta(K_full) implies U_eta(K_long), by containment. Via(12), U_eta(W_full) also implies U_eta(K_long).
* U_eta(K_long) implies U_min(eta,2delta)(K_full), by(13)-(14), and hence U_min(eta,2delta)(W_full) by(12).

The exponent can drop to min(eta,2delta) in the second direction. No same-exponent equivalence between K_long and K_full is asserted when eta>2delta.

**Corollary 5 (three-way existence of a positive saving).** For every one fixed delta in(0,1/4], the following are equivalent:

1. U_eta(K_full) for SOME0<eta<=1;
2. U_eta(K_long) for SOME0<eta<=1;
3. U_eta(W_full) for SOME0<eta<=1.

Proof. Propositions3-4 prove every direction, with a possibly smaller positive eta. If a larger positive exponent is ever supplied, it may first be weakened to eta<=1, so the stated existence criterion loses no positive-saving possibility. This is an equivalence of precise predicates for this family. None of the three predicates is established by the equivalence theorem.

## 8. Logical inputs and limits

The finite identities, the removal of the primitive diagonal, the low-order character estimate, the toric parametrization, prime uniqueness and the moment-norm bounds are proved in Sections1-6. The elementary structure of prime-power unit groups gives the power-map fibre bound: odd prime-power units are cyclic, while for2^e the unit group has at most two cyclic2-primary factors. Thus order-dividing6 characters number at most6^omega(N). This is also why chi-to-chi^2 and chi-to-chi^3 have subpower fibres.

The only reused quantitative inputs for the final bridge are the existing [unit CRT/rational reconstruction argument](ERDOS_943_KERNEL_GRAPH.md) for(9), and the exact [short-witness bound and complement argument](ERDOS_943_SHORT_WITNESS_KERNEL_PAIRS.md) for(13)-(14). The latter's signed graph, normalization and divisor-sum sources retain their existing scope. No new external point-count theorem or primary PDF is required here.

The equivalences do not prove a new bound for the full pair count, its exact long complement, the signed covariance, six-cycles or prescribed sums. They do not replace the full target by a count of unreduced witnesses, by an unsigned raw energy, or by an integer-coordinate box on the cusp. In particular they introduce no implication whose conclusion is an already open full-frequency predicate without a hypothesis proving one of its equivalent quantitative conditions. The [pointwise kernel-pair](POWERFUL_SUM_GRAPHS.md#square-full-kernel-pairs-at-the-balanced-scale) and [six-cycle](ERDOS_943_SIX_CYCLES.md) obligations remain as stated.

