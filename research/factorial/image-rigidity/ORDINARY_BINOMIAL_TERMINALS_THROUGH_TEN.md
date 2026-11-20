# Every ordinary binomial terminal through ten

For integers a>=0,d>=1 define the ACTUAL reversed moment polynomial

    A_(m;a,d)(x)=sum_(k=0)^m binom(m,k)*(am+1)_(dk)*x^(m-k).

**Theorem OrdinaryBinomialTerminalsThroughTen.** For every terminal
2<=N<=10 and every a>=0,d>=1, A_(N-1;a,d) and A_(N;a,d) are coprime
over Q. Equivalently, over EVERY characteristic-zero field, no nonzero
ordinary binomial X^a*(A+B X^d) has both factorial moments N-1,N zero
at any of these terminals. This is a finite TERMINAL prefix with
unbounded exponents and arbitrary field coefficients. No all-terminal
theorem or pairwise coprimality of all ten moment polynomials is claimed.

The proof uses ONLY universal complement arguments and one new direct
finite computation. No earlier terminal-six/fifth certificate, all-order
constant-family finite prefix or adjacent-exponent theorem is a premise.
The computation treats every remaining support after the universal
complements; it is not an extrapolation from sampled exponents.

## 1. The exact high-gap complement already in the repository

Theorem thm:rational-annulus-cutoff of
lower_exponent_initial_and_boundary.tex applies at a terminal N to any
positive rational C_(N-1),C_N satisfying

    (1+1/C_j)^j<2, j=N-1,N.

For a positive integer d_0 put B_j=((d_0+1)j)!/j!. Its other hypothesis is

    B_N^(N-1)>(C_(N-1)*C_N)^(N*(N-1))*B_(N-1)^N.        (F1)

It then proves coprimality of the actual pair for EVERY a>=1,d>=d_0.
The endpoint-normalized root annuli are disjoint at a=1,d=d_0 by (F1).
The same source's lemma lem:monotone-centers proves that the ratio of
their centers increases with BOTH positive real exponents a,d, so this
propagation is uniform in all larger integer exponents.

The following entirely rational choices satisfy all hypotheses:

    N              3    4    5    6    7    8    9    10
    C_(N-1)        3    4    6    7    9   10   12    13
    C_N            4    6    7    9   10   12   13    14
    d_0            7   12   18   24   30   37   44    50.

The checker verifies (C_j+1)^j<2*C_j^j and (F1) directly with integers.
It does not infer the inequality from approximate logarithms or from a
numerical plot. The retained theorem uses the unnormalized forward
moment polynomial M_j(t)=sum binom(j,k)(aj+dk)!t^k. Division by
(aj)! and reversal transform it into A_j; all constant and leading
coefficients are nonzero, so these changes preserve pair coprimality.

The retained sparse theorem d^2>=4(a+d)log(8N) is another sufficient
complement, but is NOT needed by this smaller proof. A proposed weaker
log(100N) criterion would add no mathematical scope here.

## 2. The sharper lower-exponent complement and boundary supports

Corollary ImprovedLinearLowerExponentCriterion of
GAMMA_PRODUCT_SHARP_LOWER_ANNULUS_AND_SHAPE_SEPARATION.md proves that

    a>=d(N-2)                                          (F2)

makes all first N actual positive-Gamma moment polynomials pairwise
coprime. For ordinary shape beta=1 it applies to the actual order shapes
am+1. Equality in (F2) is included. This is a classical Gamma-root
bound/Grace-composition application, not a new root-circle priority claim.

Consequently at positive a the only remaining supports at terminal N
are

    1<=d<d_0(N), 1<=a<d(N-2).                          (F3)

INCLUDING d=1 in this finite set avoids any adjacent-exponent all-order
premise. For a=0, Theorem thm:large in
ALL_EXPONENTS_AND_LARGE_EXPONENT_CUTOFF.tex says the constant-plus-power
moment polynomial of order N is squarefree if 2^d>=4N^2. All d>=9
satisfy this for N<=10, since 512>=400. The actual reversed polynomials
have the exact Appell relation A_N'=N*A_(N-1); hence this squarefreeness
proves consecutive coprimality. Only a=0,1<=d<=8 remains to be tested.
No completed all-exponent constant-family finite prefix is a premise.
Terminal N=2 needs no computation: (F2) holds for every a>=0 because
its threshold is d*(2-2)=0. Every remaining terminal is3 through10.

## 3. The complete remaining domain and good reduction

For EACH N=3,...,10 take(F3) and the eight supports
(a,d)=(0,1),...,(0,8). The eight case counts are respectively

    23,129,450,1089,2154,3968,6587,9759; total24159.       (F4)

For positive a the count is
(N-2)*(d_0-1)*d_0/2-(d_0-1); add eight at a=0.
The largest factorial argument in either actual polynomial is 4400,
attained at N=10,d=49,a=391. At lower order the ACTUAL shape is
a(N-1)+1, while at upper order it is aN+1. A shared-shape discriminant
is not substituted for this pair.

The exact verifier tests EVERY support in (F4). Its default prime is
65537; any exceptional reduction is retried at the verified primes
1000003 and1000033. Each exceeds4400. Therefore every factorial used in
normalizing either polynomial is a field unit; the full degrees N-1,N
and the monic leading coefficients are preserved at every attempt.

For every attempt it obtains the full extended gcd and checks the
polynomial identity U*A_N+V*A_(N-1)=g. It accepts a support only when
g=1. The receipt records all default-prime exceptions, every attempted
prime, complete per-terminal counts and case/polynomial/Bézout streams.
Coefficient vectors are independently reconstructed at EVERY attempt
by a binomial/sliding-rising recurrence, as well as by factorial tables.

A successful full-degree modular gcd implies that the actual integer
resultant has nonzero reduction at that prime. Thus it is nonzero over
Q. Failure at a chosen prime is not treated as a rational common root,
and success at any prescribed prime is not claimed necessary. The
completed deterministic replay, not a sampled subset, supplies (F4).
The actual run accepts all24159 supports. The sole exceptional default
reduction is(N,a,d)=(9,54,8), with a linear gcd modulo65537; it succeeds
modulo1000003. Therefore24160 full polynomial Bezout identities and48320
independent coefficient-vector comparisons were actually checked.

## 4. Completion and genuine field equivalence

For terminals3 through10, Sections1--2 prove all supports outside(F4);
the complete finite replay proves all supports inside it. Terminal2 is
entirely covered by(F2). This proves the stated rational coprimality
without importing any earlier finite terminal computation.

For B!=0 the exact factorial identity is

    Lambda[(A X^a+B X^(a+d))^m]
       =B^m*(am)!*A_(m;a,d)(A/B).                       (F5)

A rational Bezout identity transfers to every characteristic-zero field
and excludes simultaneous vanishing. If one phase coefficient is zero,
each positive moment is a nonzero coefficient power times a factorial.
Conversely a nonconstant rational gcd would have a complex root x and
the phase X^a*(x+X^d) would violate(F5). Thus universal complex or
universal characteristic-zero coefficient validity is equivalent to
the rational pair statement; this is not an iff over one prescribed
non-algebraically-closed field.

Reproduce the finite computation with the Python-FLINT interpreter:

    python verify_ordinary_binomial_terminal_ten.py --repo-root REPO \
      --output review.json

The mandatory ledger pins all proof inputs and the exact domain. Exact
identities, source guards and complete finite replay are supporting
computational proof evidence, not formal verification or external review.
The theorem leaves terminals N>=11 and the unrestricted all-order
ordinary binomial problem open.
