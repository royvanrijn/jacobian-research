# Dimension-dependent prime divisibility and Frobenius carries

16 November 2026. Internal written proofs covering arbitrary dimension,
polynomial phase degree and formal symbol order. No novelty, independent
replay, formal verification or external review is claimed. See the
[checker](verify_prime_carries.py), [receipt](prime_carry_verification.json)
and [source audit](prime_carry_source_audit.json).

The previous [prime-index drift theorem](PRIME_INDEX_CONGRUENCES_AND_ALL_ORDER_DRIFT_REDUCTION.md)
isolated the first-order symbol. Its higher-order replacement contains a
genuine interaction term. We prove the OPTIMAL uniform p-divisibility
exponent, including the precise dependence on spatial dimension, and give
the explicit second-level congruence. Neither proves the general formal
apolar-gap implication or decides either full exceptional MZ case.

## 1. Arithmetic setting and the two Frobenius carries

Fix a prime p and a p-torsion-free Z_(p)-algebra R. Let P in R[x_1,...,x_n]
and lambda in R[[Z_1,...,Z_n]] have order at least r>=1. In particular
lambda(0)=0. Write T=lambda(partial). Every action below is locally finite.

For a polynomial or formal series f=sum_a f_a X^a define

    F_p(f)=sum_a f_a^p X^(p a),
    C_p(f)=(f^p-F_p(f))/p.                                 (1)

The numerator in (1) is divisible by p by the multinomial theorem;
there is no division by p of an arbitrary p-adic coefficient. These are
explicit coefficientwise constructions. F_p is not asserted additive
over R, nor a ring endomorphism of the characteristic-zero coefficient
ring. Modulo p it is the usual Frobenius.

Expand the actual output as FOUR terms:

    T^pP^p = F_p(lambda)(partial)F_p(P)
       +p F_p(lambda)(partial)C_p(P)
       +p C_p(lambda)(partial)F_p(P)
       +p^2 C_p(lambda)(partial)C_p(P).                    (2)

All carries and operators can be truncated before acting, so (2) is an
identity of polynomial outputs over R, even for infinite symbols.

## 2. The optimal universal valuation exponent

**GVCPRIMEVALUATION20261116.** In the setting above every coefficient of
T^pP^p is divisible by p^e, equivalently by (p!)^e in R, where

    e=e(n,r)=min(r, max(2,r-n+3)).                         (3)

The formula includes r=1, giving e=1. For r=2 it gives e=2 in EVERY
dimension. For n<=3 it gives e=r in EVERY operator order. For n>=4,r>=3
the guaranteed exponent can be smaller than r. This is a universal
coefficientwise bound, not a nonvanishing assertion for a specified pair.

Here is the proof for every prime, including p=2. The diagonal operator
F_p(lambda)(partial) has terms partial^(p alpha) with |alpha|>=r.
Each such derivative on an integral polynomial is divisible by
product_i (p alpha_i)!, hence by p^r. The first two terms of (2) are
therefore divisible by p^r and p^(r+1), respectively.

Every monomial partial^gamma in C_p(lambda)(partial) has |gamma|>=rp.
On a Frobenius monomial x^(p beta), the derivative coefficient contains
at least ceil(gamma_i/p) factors of p in coordinate i (or is zero).
Thus its valuation is at least sum_i ceil(gamma_i/p)>=r. The third
term of (2) is divisible by p^(r+1).

On the unrestricted polynomial C_p(P), the same derivative is divisible
by product_i gamma_i!, with valuation at least

    sum_i floor(gamma_i/p)
      >=floor(|gamma|/p)-(n-1)>=r-n+1.

Use zero instead if this lower bound is negative. The fourth term in
(2) is divisible by p^(2+max(0,r-n+1)). Taking the minimum with r
proves (3). No analytic growth, convergence or infinite-ring ACC is used.

For degP<=D with D>=r, only the symbol jet

    |alpha|<=(D-r)p+r                                     (4)

can affect T^pP^p: one higher factor plus p-1 factors of order>=r
differentiates more than Dp times. If D<r the entire output is zero.
Thus a p-integrality hypothesis for a displayed formal symbol is needed
only through (4), not silently beyond every finite jet.

## 3. Sharpness for fixed dimension and order, in all large primes

The exponent (3) cannot be increased as a universal polynomial
coefficient identity, even with a homogeneous finite symbol and homogeneous
phase of degree exactly r. The examples below are FIXED for each n,r;
only the tested prime varies.

If e=r, take T=partial_1^r and P=x_1^r, padding unused variables.
Then T^pP^p=(rp)!, which has valuation EXACTLY r for p>r.

Suppose e<r. Then n>=4,r>=3. Put s=min(n,r+1), so s>=4, and use

    lambda=sum_(i=1)^s u_i Z_i^r,
    P=sum_(i=1)^s x_i^r,                                 (5)

over R=Z_(p)[u_1,...,u_s], with unused variables padded. The phase and
symbol in (5) do not depend on p. For each large prime set

    k=floor((s-1)p/(s r)),
    k_2=...=k_s=k,   k_1=p-(s-1)k.                        (6)

For p>s r(s-1), all k_i are positive and less than p. Moreover

    floor(r k_i/p)=0 for i>=2,
    floor(r k_1/p)=r-s+1.                                (7)

Indeed r k/p is just below (s-1)/s; the error in r k_1/p from
r-s+1+(s-1)/s is nonnegative and less than 1/s. Also p>r makes every
r k_i<p^2, so no additional p^2-factor enters its factorial.

The coefficient of u_1^k_1...u_s^k_s in the CONSTANT output T^pP^p is

    (p!/(product_i k_i!))^2 product_i(r k_i)!.              (8)

There is just one symbol allocation and one phase allocation for that
parameter monomial. Its valuation, by (7), is exactly

    2+r-s+1 = max(2,r-n+3) = e(n,r).                      (9)

There is no cancellation between different u-monomials. One may regard
(5) as one fixed pair over Q(u_1,...,u_s) with the coefficientwise
p-adic Gauss valuation. Equations (6)-(9) show sharpness at EVERY
sufficiently large prime, without a prime search or changing the pair.

For example n=4,r=3 has e=2, so dividing by (p!)^3 need not even leave
p-integral coefficients. In dimensions n<=3, in contrast, (p!)^r
division is always p-integral. These statements concern NORMALIZATION,
not positive GVC or the actual nilpotence-height slope. The sharpness
families (5) are not asserted pure or eventually safe.

## 4. The exact second-level congruence

Assume now r>=2 and let lambda_2 be the ordinary homogeneous quadratic
symbol component (zero if r>2). Then

**GVCQUADCARRY20261116.** In R[x],

    T^pP^p/(p!)^2
      = (lambda_2(partial)P)^p
         +C_p(lambda)(partial)C_p(P) modulo p.            (10)

Thus the normalized residue is NOT determined by the quadratic symbol
acting once on P. The interaction of BOTH carries is essential.

The second and third terms of (2) are divisible by p^3 when
r>=2. In its first term, diagonal symbol components of order>2 are
also divisible by p^3. For |alpha|=2 the Hasse derivative formula gives

    partial^(p alpha)F_p(P)/(p!)^2
       = (partial^alpha P)^p modulo p.                   (11)

To check the constants in (11), factor the derivative as
(product_i(p alpha_i)!)*Hasse_(p alpha). Frobenius makes its Hasse
coefficient (Hasse_alpha P)^p modulo p. The remaining factorial ratio
modulo p is product_i alpha_i!, by the central-binomial congruence
binomial(2p,p)=2 modulo p. Multiplying back agrees with
(partial^alphaP)^p, including the diagonal alpha_i=2,p=2 case where
both sides are zero modulo p. Summing the quadratic terms gives the
first summand in (10).

Finally p^2/(p!)^2=((p-1)!)^-2 is 1 modulo p, by Wilson's theorem
(also for p=2). This gives exactly the carry coefficient written in
(10). These elementary constants can equivalently be checked by the
usual finite-field product and binomial coefficient proofs.

An all-order countercontrol to omitting the carry is already quadratic:

    T=partial_x^2+partial_y^2+partial_z^2,
    P=x^2+y^2+z^2.

Direct radial differentiation gives T^mP^m=(2m+1)! for EVERY m. Hence
for EVERY odd prime p,

    T^pP^p/(p!)^2=(2p+1)binomial(2p,p)=2 modulo p,
    (lambda_2(partial)P)^p=6 modulo p,
    C_p(lambda)(partial)C_p(P)=-4 modulo p.               (12)

The carry is nonzero at every odd prime. These inputs are NOT pure;
they refute the proposed universal extraction identity, not a GVC theorem.

## 5. Higher symbol orders survive the carry even at cubic phase degree

There is also no rule discarding higher symbol orders from (10). In four
variables put

    P=x_1^3+x_2^3+x_3^3+x_4^3,
    lambda=Z_1 Z_2+sum_(i=1)^4 u_i Z_i^3.                 (13)

Here lambda_2(partial)P=0. The CONSTANT output at moment p uses only
the cubic symbol terms: one quadratic factor would leave positive
output degree. Thus its coefficients (8), with r=3,s=4, survive in
the residue after (p!)^2 division. For every sufficiently large prime
the profile (6) has all 3k_i<p, so its normalized coefficient is a
unit modulo p. Consequently (10)'s carry term is genuinely nonzero
and depends on HIGHER symbol orders, despite the zero quadratic action.
At p=7 the explicit profile (2,2,2,1) has normalized coefficient
180^3*6=1 modulo7. No four-variable expansion is required.

Again (13) is a coefficient-identity boundary, not a pure pair or a new
counterexample to eventuality-to-linearity. Imposing pure nullity instead
gives the necessary arithmetic relation

    C_p(lambda)(partial)C_p(P)
       =-(lambda_2(partial)P)^p modulo p                  (14)

at every good prime. The two carries in (14) vary with p. Their vanishing
cannot be inferred from symbol order or a fixed low jet, and the previous
drift-specialization proof cannot turn (14) into a characteristic-zero
identity for lambda_2P by discarding them.

## 6. Scope and the next decisive gate

Equations (3) and (10) cover the ENTIRE residual formal-symbol class once
the relevant jets are p-integral, with no bound on phase degree or
operator order. Arbitrary formal denominators retain the earlier explicit
arithmetic boundary. The formulas constrain pure pairs but do not make
all pure pairs safe, prove quasi-subadditivity of actual heights, or decide
unbounded actual deficit versus a positive deficit fraction.

The retained nonzero-drift dimension classification is unchanged. The
order>=2 case now has a correct replacement for its tempting but false
prime extraction. There is no justified larger prime campaign: the
remaining issue is structural control of the varying carry interaction
in actual pure/eventually safe pairs, or another general apolar-gap proof.

Exploratory low-degree apolar checks during this turn did not establish
general nesting of Ann(P^m), a finite Rees model, or a height inequality.
No such statement is registered. The recent symmetric-powers source
retains its stated subalgebra/quotient and maximum-dimension hypotheses;
bounded catalecticant ranks are not substituted for those hypotheses.

The checker verifies a few exact carry identities, factorial valuations
and sharp parameter-monomial coefficients, plus refusal of nonintegral
reduction. Universal optimality and all-prime statements are written
proofs above, not conclusions from samples. All protected Gaussian/Haar,
ordinary-ideal, cubic-prefix, phase-quotient, drift and radial-transfer
records remain unchanged. Both complete exceptional MZ questions and the
general formal apolar-gap implication stay OPEN;external review is deferred.
