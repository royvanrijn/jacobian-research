# Strong Factorial for the ordinary relative-gap-three family

## 1. The complete all-terminal statement

**Theorem RelativeGapThreeStrongFactorial.**

Let K be any field of characteristic zero. For all integers a>=2,
d>=3a and N>=2, and all A,B in K not both zero, at least one of
\[
 {\cal L}\!\left([X^a(A+BX^d)]^{N-1}\right),\qquad
 {\cal L}\!\left([X^a(A+BX^d)]^{N}\right)
\]
is nonzero, where the ordinary factorial functional is
\({\cal L}(X^j)=j!\).

The finite obligation in the independently retained
[uniform relative-gap tail and finite-equivalence theorem](RELATIVE_GAP_THREE_LOWER_WEIGHTED_SADDLE_TAIL.md)
is now completed. Its analytical proof and source bytes are unchanged.
This theorem combines that proof with the completed exact finite
calculation described below. It concerns the ordinary factorial
functional and the domain d>=3a; it does not settle arbitrary binomial
supports, all positive Gamma parameters, or unrestricted Strong Factorial.

## 2. The actual polynomials and good reduction

For m>=1 define the INTEGER polynomial
\[
 P_{m;a,d}(t)=\sum_{k=0}^{m}\binom mk
       \frac{(am+dk)!}{(am)!}\,t^k
   =\sum_{k=0}^{m}\binom mk(am+1)_{dk}\,t^k.                       \tag{C1}
\]
It has constant coefficient1 and exact degree m. The lower order m=N-1
and upper order m=N retain their DIFFERENT actual parameters
a(N-1)+1 and aN+1 throughout.

If A!=0 and t=B/A, then
\[
 {\cal L}([X^a(A+BX^d)]^m)=A^m(am)!\,P_{m;a,d}(t).                \tag{C2}
\]
If A=0 or B=0, the binomial is a nonzero monomial, and every positive
factorial moment is nonzero in characteristic zero. Hence it suffices
to prove coprimality of the actual consecutive rational polynomials
in(C1).

**Lemma DegreePreservingGoodReduction.**

Let f,g be integer polynomials and p a prime preserving BOTH degrees.
If their reductions are coprime over F_p, then f,g are coprime over Q.

Indeed, a nonconstant common factor over Q yields a primitive integer
common factor by Gauss's lemma, after removing contents. Since p
preserves both full degrees, its leading coefficient and the leading
coefficients of the corresponding quotient factors do not vanish
modulo p. Its reduction is consequently a nonconstant common factor,
a contradiction. Equivalently the full-degree Sylvester resultant
specializes to a nonzero element modulo p and therefore is nonzero
over Q. This implication is ONE-WAY: rational coprimality does not force
success at a prescribed prime.

For(C1) a prime p>(a+d)m preserves the leading coefficient:
\[
 \operatorname{lc}P_{m;a,d}=\frac{((a+d)m)!}{(am)!}\not\equiv0
       \pmod p.
\]
All factorial arguments are strictly below p. They may therefore be
computed with factorial and inverse-factorial tables. The actual code uses
\[
 c_k=m!\,((am)!)^{-1}(am+dk)!\,(k!)^{-1}((m-k)!)^{-1}\pmod p,      \tag{C3}
\]
which equals the coefficient in(C1). No term or top coefficient is
discarded.

## 3. The completed finite calculation

Put T(a,d)=max(4000,128(a+d)). The completed domain is EXACTLY
\[
 9\le d\le63,\quad 3\le a\le\lfloor d/3\rfloor,\quad
 2\le N<T(a,d).                                                   \tag{C4}
\]
There are532 distinct support pairs(a,d) and3,694,040 terminal pairs.
The largest terminal is10751 and the largest factorial argument is
903084. The two trial-division-certified primes used are
\[
 p_0=1000003,\qquad p_1=1000033.
\]
Both exceed903084, so the preceding good-reduction lemma applies to
EVERY actual pair in(C4).

The exact executed checker is retained as
[verify_relative_gap_three_finite_box.py](verify_relative_gap_three_finite_box.py),
with SHA256
e33fe9e57513b9a0ff74c1e5a5334718e57179450d87cdae2964da27d34ebe32.
It constructs the complete polynomial at every order m, verifies its
degree m, and tests gcd(P_(N-1),P_N) at every listed terminal. The final
[completed receipt](relative_gap_three_finite_box_verification.json)
has SHA256
95bf53dba99f375ea908aedcaa80a36162a9824a0e193840150d81663f6d9036.
It contains all532 completed support rows, with exact counts and
terminal-stream digests recording complete coverage of all3,694,040
tests, without an early exit or degree-dropping exception. Individual
per-terminal coefficient vectors and Bezout witnesses are not stored in
the receipt; they are computed and checked by the executed source.

At five terminals the first prime gives a linear gcd. The checker then
constructs both actual full-degree polynomials over the second prime
and computes and verifies the complete Bezout identity
s P_N+t P_(N-1)=1:

| a | d | N | degree of gcd at p0 | success prime |
|---|---|---|---|---|
|11|41|6142|1|1000033|
|4|48|4491|1|1000033|
|13|61|3147|1|1000033|
|19|61|3843|1|1000033|
|19|62|2254|1|1000033|

All other terminals have monic constant gcd at p0. Thus at each terminal
there is an actual degree-preserving prime at which the full pair is
coprime. The good-reduction lemma proves every rational coprimality
assertion in(C4).

The execution used Python3.12.3, python-flint0.9.0 and FLINT3.6.0, and
completed in4768.533664967981 seconds. Its source and completed receipt
are copied with their EXACT executed bytes. In particular their
launch-time scope string calls the analytical tail a separate
"UNPROVED premise". That string delimits what the FINITE CHECKER itself
establishes. The separate analytical theorem is now proved in the
unchanged retained source cited in Section1. The finite source is not
rewritten to claim it proved that analytical theorem.

The companion completion validator checks the pinned source/receipt,
every support and terminal stream digest, exact degree-preservation
bounds, the literal coefficient formula and all five full Bezout
fallbacks. It does not pretend to rerun all3,694,040 gcd computations.
The full finite run can be reproduced separately with the executed
checker. A second independently implemented runner was still running
when this proof was assembled; no completed independent-run result or
external/formal assurance is asserted.

## 4. Why the finite box closes every terminal in the stated family

Fix a>=2,d>=3a,N>=2.

* If a=2, use the retained
  [ordinary lower-two Strong Factorial theorem](ORDINARY_LOWER_TWO_STRONG_FACTORIAL.md),
  which covers EVERY gap and terminal.
* Suppose a>=3. If d>=64, the retained relative-gap proof, Section7,
  already proves every terminal by overlap of the sparse theorem with
  the uniform weighted tail.
* If a>=3 and d<=63, then9<=d<=63 and a<=floor(d/3).
  For N>=T(a,d), apply the retained uniform weighted tail.
  For2<=N<T(a,d), apply the completed finite calculation(C4).

These cases cover all parameters, including the equality d=3a and
N=T(a,d). This is precisely the reverse implication in the retained
genuine finite-equivalence theorem, now supplied with its indispensable
complete finite source. Conversely the all-terminal statement implies
the finite rational box by restriction over C. No reverse implication
to prescribed-prime modular success is used.

For clarity, the sparse/analytic overlap in the d>=64 branch is an
already proved all-parameter argument: D=a+d<=4d/3, and
d^2/(4D)>=3d/16. At d=64 the exact exponential comparison implies
d^2>=4D log(8N) for2<=N<=128D, with the comparison increasing thereafter.
The analytical theorem covers N>=max(4000,128D); here128D>4000.
This unbounded branch is not extrapolated from finite modular data.

Finally rational Bezout identities for each actual pair(C1) remain
valid over EVERY characteristic-zero extension field K. Formula(C2)
then excludes a common zero in K. This proves the stated theorem,
including arbitrary coefficients in K. The equivalence of universal
characteristic-zero validity with validity over C does not assert that
a failed rational gcd has a root in every specified nonclosed field.

## 5. Reproduction and verification scope

From the repository root, replay the pinned completed artifact and the
five fallback computations:

    .venv/bin/python research/factorial/image-rigidity/verify_relative_gap_three_completion.py \
      --repo-root . --output /tmp/relative_gap_three_completion_review.json

To repeat the ENTIRE3,694,040-case calculation, leaving the retained
completed receipt unchanged:

    .venv/bin/python research/factorial/image-rigidity/verify_relative_gap_three_finite_box.py \
      --jobs 4 --output /tmp/relative_gap_three_finite_full_replay.json

Replay the separate analytical support:

    python3 research/factorial/image-rigidity/verify_relative_gap_three_weighted_tail.py \
      --repo-root . --output /tmp/relative_gap_three_analytic_review.json

The imported finite-field arithmetic consists of python-flint/FLINT
exact polynomial gcd and extended gcd. The degree-preserving transfer,
normalization and all-terminal decomposition are written proofs.
Internal source/code reviews and exact computations do not constitute
external human review or complete formal verification.

Updated 20 November 2026.
