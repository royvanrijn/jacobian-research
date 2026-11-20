# Independent complete finite confirmation for the relative-gap-three family

## The exact finite statement

For positive integers a,d and m>=1 define the actual integer moment polynomial

    P_(m;a,d)(t)=sum_(k=0)^m binom(m,k)*(am+1)_(dk)*t^k.

Here (x)_j is the rising factorial. Its constant coefficient is1 and its
leading coefficient is ((a+d)m)!/(am)!. The lower and upper moments at
N-1,N have DIFFERENT actual shapes a(N-1)+1 and aN+1.

**Theorem IndependentRelativeGapThreeFiniteConfirmation.** Every actual
pair P_(N-1;a,d),P_(N;a,d) is coprime over Q in exactly this finite domain:

    9<=d<=63, 3<=a<=floor(d/3),
    2<=N<max(4000,128(a+d)).

There are532 support pairs and3,694,040 terminal pairs. The largest upper
degree is10751 and the largest factorial argument is903084. This note
establishes that finite statement only. An all-terminal relative-gap
theorem additionally uses the separately retained analytical tail,
sparse overlap and ordinary lower-two theorem. No assertion for all
ordinary binomials or all positive Gamma shapes follows here.

## Coefficients, degree preservation and exact certificates

For f=X^a(A+B X^d), with A!=0,

    L(f^m)=A^m*(am)!*P_(m;a,d)(B/A).

Monomial cases have nonzero factorial moments. The executed source builds
each coefficient using the recurrence

    c_0=1,
    c_(k+1)=c_k*(m-k)/(k+1)
                *product_(j=1)^d(am+dk+j).

Its modular inverse table uses inv(j)=-(p//j)*inv(p%j) modp; its independent
sliding product table stores R(t)=product_(j=0)^(d-1)(t+j), with
R(t+1)=R(t)*(t+d)/t. All used indices are within the allocated tables:
the largest required rising-product start is (a+d)m-d+1, whose product
ends at (a+d)m<=903084. This construction uses no factorial tables.
The preceding polynomial is reused only if it was built at its OWN
previous order for the same support and prime.

The first trial prime is1000033, distinct from the original primary
confirmation's1000003. The ordered fallback bank is
1000037,1000039,1000003. Trial division certifies all four primes.
Every prime exceeds903084, so all factorial factors, denominators and
full leading coefficients are nonzero in the corresponding field.
The exact full degrees N-1,N and constant coefficients1 are checked at
every attempt.

At EVERY attempt the executed source computes q,s,t over the finite field
and multiplies out s*P_N+t*P_(N-1)=q. A successful attempt requires q=1.
This identity certifies modular coprimality independently of trusting only
a returned gcd degree. The identity is also checked at unsuccessful
attempts before trying the next prime.

**Lemma FiniteDegreePreservingTransfer.** If an integer-polynomial pair
retains both full degrees modulo a prime and its modular reductions have
a Bezout identity equal to1, then it is coprime over Q. Indeed the
full-degree Sylvester resultant is nonzero modulo that prime, hence is
a nonzero integer. Equivalently, a primitive common integer factor would
remain a nonconstant common modular factor by Gauss's lemma. This is a
ONE-WAY implication: coprimality over Q does not force success at a
prescribed prime.

## The actual completed execution

The retained executed source is
[independent_verify_relative_gap_three_finite_box.py](independent_verify_relative_gap_three_finite_box.py),
with SHA256
12d554bd8fe2dc2a933d3be47b6ba791c8c1d689126588be1df59425184e5a38.
It is preserved with its exact executed bytes. Its actual completed
[receipt](independent_relative_gap_three_finite_box.json) has SHA256
b5336126eb70bbb25258b215af1a3e3234ee58e6737ddb337a80dc9b4412a865.
The final receipt exists only after the full worker domain completes;
it records complete_finite_box=true,532 supports and3,694,040 terminal
pairs. Benchmark/partial runs do not satisfy this flag and full domain.

The run checked3,694,045 full Bezout identities. There were3,694,035
successes at1000033 and exactly five successes at1000037:

| a | d | N | first-prime gcd degree | success prime |
|---|---|---|---|---|
|5|31|1168|1|1000037|
|9|33|2932|1|1000037|
|10|50|7079|1|1000037|
|5|57|7734|1|1000037|
|4|60|7495|1|1000037|

No attempt beyond the second prime was needed. Each support row stores
its full endpoint/count/degree bounds, successful-prime counts, complete
fallback attempt list and a SHA256 stream of(a,d,N,successful_prime).
The stream is checked over EVERY listed terminal, without sampled
anchors or an early exit. The receipt does not store all millions of
coefficient vectors or Bezout polynomials: those were computed and
verified by the executed source. A digest authenticates coverage and
provenance; it does not itself recompute polynomial arithmetic.

Execution used Python3.12.3, python-flint0.9.0 and FLINT3.6.0, with four
workers, and completed in7528.239744734019 seconds. The source's54 literal
coefficient vectors and actual wrong-shape/common-factor/corrupt-Bezout
controls all passed. This is additional INTERNAL, differently implemented
finite confirmation. It is not external human review or formal verification.

## Portable artifact confirmation

Run the short artifact validator from a checkout:

    .venv/bin/python research/factorial/image-rigidity/verify_relative_gap_three_independent_confirmation.py \
      --output /tmp/relative_gap_independent_review.json

All essential files are local to the packet, so no repository math file
is needed to prove this FINITE statement. The mandatory source ledger
pins this note, the validator, its input, the exact executed source and
the actual completed receipt. The validator reconstructs all532 support
rows and3,694,040 terminal tags, checks every per-prime/attempt count,
freshly replays54 literal vectors and the source's negative controls,
and reconstructs ALL TEN recorded fallback attempts from a separate
factorial-table coefficient builder before checking their full Bezout
identities. It does not rerun the full3,694,040-case job.

To repeat the full executed computation separately, preserving the
retained completed receipt:

    .venv/bin/python research/factorial/image-rigidity/independent_verify_relative_gap_three_finite_box.py \
      --jobs 4 --output /tmp/relative_gap_independent_full_replay.json

The original primary source/receipt and the analytical tail are unchanged.
The finite theorem above can supply a separate proof route for the SAME
finite obligation; it does not supply the analytical premises by itself.

Updated 20 November 2026. No historical-priority, external-review or
complete formal-verification assurance is asserted.
