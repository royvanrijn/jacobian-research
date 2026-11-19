# Primitive Image counterexamples in asymptotically all prime balanced degrees

Date: 19 November 2026. This is a number-theoretic corollary of the retained all-order construction, with no additional analytic construction of an Image witness.

## 1. Statement

Work over any characteristic-zero field K. Let H be the set of prime balanced degrees supplied by Theorem HyperbolicOddHeightImage in [the retained construction](HYPERBOLIC_ODD_HEIGHT_IMAGE_WITNESSES.md), Section6. Its parameters are every odd integer a>=1 and integer n>=1, with

    d=(2a+3)n+a+2,  2d-1=(2n+1)(2a+3).                (P1)

At every such degree the theorem supplies an explicit two-pair polynomial F that is radial-primitive and not a nonzero scalar times a proper power, even over coefficient-field extensions. Every positive pure contraction E2(F^m) vanishes. The ONE fixed original marker Z has E2(ZF^(2r)) nonzero for every r>=1.

**Corollary PrimeDegreeDensityOne.** The explicitly constructed set H has relative natural density one among primes:

    #{p<=X: p in H}/pi(X) -> 1,
    #{p<=X: p in H} ~ X/log X.                        (P2)

Consequently the set C of ALL prime degrees admitting radial-primitive, non-proper-power two-pair Image counterexamples also has relative natural density one. This does not assert that every sufficiently large prime belongs to H or C, or give an effective error term for the omitted prime degrees.

**Corollary CubicHeightQuarterDensity.** The a=1 subfamily alone gives every prime p=3 mod10 other than3. For X>=3 its exact prime-degree count is

    pi(X;10,3)-1 ~ X/(4 log X).                       (P3)

Thus this subfamily has relative density1/4 among primes. The quarter is not an upper bound on the full construction: allowing a to vary yields(P2).

Balanced degree p means bidegree(p,p), or ordinary total degree2p in the four original variables. The statement concerns a set of degrees, not a count of polynomials or parameter pairs.

## 2. The quarter-density corollary

Set a=1 in(P1), giving d=5n+3. For d an odd prime, n is even: odd n gives an even integer d>=8. Thus n=2k, k>=1, and d=10k+3. Conversely a prime p=3 mod10 with p>3 has n=(p-3)/5 a positive even integer, so the retained theorem applies. The only removed prime in this residue class is3.

Dirichlet's theorem gives infinitude because gcd(3,10)=1. The stronger fixed-modulus prime number theorem gives pi(X;10,3)~X/(phi(10)log X), and phi(10)=4. Removing3 changes the count by1 and leaves the asymptotic unchanged. Comparing with pi(X)~X/log X proves(P3) and relative density1/4.

## 3. Prime divisors give many other progressions

Let q be any prime with q=1 mod4, so q>=5. Choose

    a=(q-3)/2,  r_q=(q+1)/2.

Here a is odd and at least1, r_q is odd, and gcd(r_q,q)=1. Formula(P1) becomes

    d=q*n+r_q.                                       (P4)

Every odd prime p=r_q modq has even n=(p-r_q)/q. Once p>r_q, this n is at least2, so p is supplied by the construction. The only possible positive prime exception in that residue class is p=r_q (the inadmissible n=0 case).

Equivalently, define P_* by the existence of a PRIME proper divisor q of2p-1 with q=1 mod4. Then P_* is exactly the prime-degree spectrum of the subfamily with2a+3 prime. The inverse parameters are

    a=(q-3)/2, n=((2p-1)/q-1)/2.                     (P5)

Because p is odd and q=1 mod4, the proper cofactor is1 mod4 and at least5; hence n is even and at least2. Conversely(P1) gives that divisor for this subfamily. The prime2 is absent and immaterial to density.

This divisor criterion is sufficient for the full family; it is NOT necessary there. For example41 is supplied by a=3,n=4, with2*41-1=9*9, and has no prime divisor1 mod4. The retained full-family factor criterion allows composite factors1 mod4 as well.

## 4. Finite CRT unions have an exact relative density

Fix a FINITE set S of primes q=1 mod4 and put M=product_(q in S)q. Let U_S be the primes satisfying p=r_q modq for at least one q in S. Finite exceptions from primes dividing M and the n=0 cases have no effect on relative density.

Among the phi(M)=product(q-1) reduced residues modulo M, avoiding the one selected nonzero residue r_q at every q leaves exactly product(q-2) residues, by the Chinese remainder theorem. The prime number theorem in each fixed reduced class modulo this FIXED M therefore gives

    lim_(X->infinity) #{p<=X:p in U_S}/pi(X)
       = 1-product_(q in S)(q-2)/(q-1)
       = 1-product_(q in S)(1-1/(q-1)).               (P6)

Every prime in U_S except finitely many belongs to P_* and thus H. Consequently the liminf of the relative density of H is at least the right side of(P6), for EVERY fixed finite S. The same lower bound holds for P_*.

## 5. The products tend to zero

The fixed-progression prime number theorem with modulus4 gives

    pi(t;4,1) ~ t/(2 log t).

Partial summation implies divergence of the prime reciprocal sum in that class:

    sum_(q<=Y,q=1 mod4) 1/q
       = pi(Y;4,1)/Y + integral_2^Y pi(t;4,1)/t^2 dt
       -> infinity.                                 (P7)

For example, beyond a fixed threshold the integrand is at least1/(4t log t). Thus sum 1/(q-1) diverges as well. Since log(1-v)<=-v for0<v<1,

    product_(q<=Y,q=1 mod4)(1-1/(q-1)) -> 0.           (P8)

Apply(P6) first with S fixed and let X tend to infinity. Only AFTER taking that limit enlarge S through the primes1 mod4. The lower bounds tend to1, while each counting ratio is at most1. This proves an actual relative natural density1 for P_*, then for H and C by inclusion and the same upper bound. Finally pi(X)~X/log X gives(P2). Equivalently the prime degrees omitted by H are o(pi(X))=o(X/log X).

No estimate uniform in M, interchange of limits, sieve bound, or conjecture about simultaneous prime values is used. Infinitude alone of the primes1 mod4 would NOT suffice for(P8); their divergent reciprocal sum is proved explicitly from fixed-modulus PNT.

## 6. Sources and evidence

The essential mathematical input is the already retained HyperbolicOddHeightImage theorem, especially(H18)--(H22): all parameters, all moments, the fixed marker, radial primitivity and proper-power exclusion. This corollary does not replace or strengthen its proof assurance.

The classical analytic import is Ivan Soprounov, *A short proof of the prime number theorem for arithmetic progressions*, Theorem2.1, printedp1, [author-hosted PDF](https://academic.csuohio.edu/soprunov-ivan/wp-content/uploads/sites/93/2023/02/primes.pdf), [retained source](sources/soprounov_pnt_arithmetic_progressions.pdf). Its actual statement and proof were inspected: fixed positive coprime residue/modulus, with count asymptotic X/(phi(q)log X). The q=1 specialization supplies ordinary PNT. The note attributes the Newman/Zagier approach and explicitly imports nonvanishing L(1,chi) in its introduction; this is a classical unconditional theorem, not any repository-specific uniform zero-free hypothesis. No year or priority is inferred from the server upload date.

Dirichlet infinitude can also use the existing retained Williams fixed-progression Mertens input; it is already a consequence of the stronger PNT import used here. CRT residue counting, partial summation and the sequential liminf squeeze are supplied explicitly above.

The [source ledger](prime_image_degree_density_source_ledger.json), [checker](verify_prime_image_degree_density.py) and [receipt](prime_image_degree_density_verification.json) pin the exact construction and primary input. Finite controls compare independent parameter and divisor descriptions, enumerate reduced CRT residues, and include the missing prime3 and the41 counterexample to an overstrong prime-divisor converse. They do not establish an asymptotic by numerical prime counts; Sections4--5 give its proof. All review and replay recorded here are internal. No external human review, formal verification or historical-priority claim is asserted.
