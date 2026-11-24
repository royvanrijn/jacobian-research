# Exact sums of two primitive discriminant −43 norms

Let
\[
 Q(x,y)=x^2+xy+11y^2,\qquad
 P=\{Q(x,y)>0:x,y\in\mathbb Z,\ \gcd(x,y)=1\}.
\]
Let \(C\) contain 1 and all positive integers supported on primes
\(p\) with \((p/43)=1\). Signs and zero coordinates are allowed;
summands are positive numerical values, irrespective of the number of
coordinate representations.

**Theorem M43 (Minus43ExactPrimitiveSumset).** Put
\[
 E=\{4,6,8,10,16,20,38,50,56,74\}.
\]
Then
\[
             C+C=\{n\ge2:2\mid n,\ 43\nmid n\}\setminus E,
                                                               \tag{M43.1}
\]
and
\[
 P+P=\{N\ge2:2\mid N,\ 1849\nmid N\}\setminus\mathcal E,
                                                               \tag{M43.2}
\]
where the complete set of eighteen primitive exceptions is
\[
\begin{split}
 \mathcal E=\{&4,6,8,10,16,20,38,50,172,258,344,430,\\
              &688,860,1634,2150,2408,3182\}.
\end{split}
\]

The proof imports the fully proved analytic finite reduction at
\(D=43\), then checks **every required finite target through
\(2^{32}=4,294,967,296\)**. There are 2,097,542,158 such targets after
the local and explicitly listed core exclusions. The exact segmented
certificate has no missing target. This closes the finite obligation
for −43 left open in the separate refinement note; it supplies no
finite coverage at −67 or −163.

## 1. The analytic and arithmetic interface

The theorem uses [HeegnerCoreAndQuotientOnsets and
HeegnerFiniteClassificationReduction](PRIMITIVE_NORM_HEEGNER_REFINED_BOUNDS.md),
Theorems H1–H2. At \(D=43\) their precise conclusions are:

- Every even 43-unit inert core at least \(2^{32}\) lies in \(C+C\).
- Every number \(ep\ge2^{32}\), with \(e\in E\) and \(p\) a split
  prime, lies in \(C+C\). This includes \(p\mid e\).
- Complete finite coverage through \(2^{32}\), outside \(E\), implies
  the global core identity by multiplicative descent.

The classical premises are specified and pinned in that note: a
maximal discriminant-43 quaternion order, complete ideal-class
quotient, principal class mass \(w=1/7\), ramified fine-point weight
\(\nu=1/154\), and the all-integral-ideal Hecke transfer and Deligne
bound outside the level. The valuation branches, odd-lattice estimate,
Euler frontier, all finite rational inequalities and the infinite
tail are proved there. This packet imports those statements; the
finite scanner neither reproves the spectral theorem nor treats a
hash check as mathematical verification.

The same source proves the necessary arithmetic identities
\[
 P=C\ \dot\cup\ 43C,\qquad
                 43m\in P+P\ \Longleftrightarrow\ m\in C+C.
                                                               \tag{M43.3}
\]
They follow from class number one and primitive ideal orientations.
All values in \(C\) are odd, nonzero squares modulo 43; \(-1\) is
nonsquare modulo 43. Thus \(C+C\) has no multiple of 43, and \(P+P\)
has no odd number or multiple of \(43^2\).

For clarity, the two scalar probabilities here are
\[
 \mathfrak m=(43-1)/12=7/2,\quad
 w=(1/2)/(7/2)=1/7,\quad
 \nu=12/(43^2-1)=1/154.
\]
There is no quaternion class-number-one assumption, no specialized
\(a_2\) cancellation, and no additional analytic conjecture.

## 2. Exact membership on each segment

For an odd positive integer \(b\le H\),
\[
 b\in C\quad\Longleftrightarrow\quad
 (b/43)=1\ \text{and no inert prime }p\le\sqrt H
                                    \text{ divides }b.    \tag{M43.4}
\]
The forward implication is immediate. For the reverse, a remaining
inert prime factor must exceed \(\sqrt H\). There is at most one,
and its exponent is one, since otherwise its product would exceed
\(H\). The character test excludes 43. All other factors are split,
so the remaining inert factor would give character −1, a contradiction.
This is a complete membership test; it does not assume a candidate is
prime or squarefree.

Let a segment contain the even targets from \(L\) through \(H\),
and let \(A=C\cap[1,49999]\). There are exactly 4635 members of
\(A\), determined by (M43.4). Put
\(b_0=\max(1,L-\max A)\), which is odd. Construct the membership
flags of \(C\) for every odd \(b\in[b_0,H-1]\).
The flags start with the 43-periodic character mask. For each odd
inert prime \(p\le\sqrt{H-1}\), delete the positions
\[
                    j\equiv-b_0\,2^{-1}\pmod p.
\]
These are exactly the positions with \(p\mid b_0+2j\).
Primes only through \(\sqrt{2^{32}}=65536\) are needed for the
whole certificate. Equation (M43.4) proves that the resulting flags
are exact for the entire segment, including its endpoints.

Pack the flags into an integer \(B\), with bit \(j\) recording
\(b_0+2j\in C\). For each \(a\in A\), shift \(B\) right by
\((L-a-b_0)/2\), or left by its negative. Bit \(i\) now records
exactly \(L+2i-a\in C\), with a positive second summand in the
stored range. OR these shifted integers. Every set target bit has
an explicit positive pair in \(C+C\); no heuristic density estimate
is involved. Extra bits outside the target segment are masked away.

The required mask contains every even target in \([L,H]\) that is
a 43-unit and not in \(E\). Its construction uses a separate
43-periodic mask and deletes only the listed exception positions.
In particular it includes targets with split factors, not just the
inert cores used in the analytic theorem.

## 3. Complete finite coverage and controls

The checker uses 128 consecutive target segments, each of span
\(2^{25}\):
\[
 L_i=2+i2^{25},\qquad H_i=\min(L_i+2^{25}-2,2^{32}),
                       \quad0\le i<128.                  \tag{M43.5}
\]
It asserts \(L_0=2\), \(L_{i+1}=H_i+2\), and \(H_{127}=2^{32}\).
Every required target bit is covered. The complete total is checked
independently against
\[
 \frac{2^{32}}2-\left\lfloor\frac{2^{32}}{86}\right\rfloor-|E|
                                  =2,097,542,158.         \tag{M43.6}
\]
All ten exceptional targets are even 43-units within the interval, so
this subtraction counts exactly the required domain. The largest
number of short summands needed on any segment was 176. The full
4635-element list remains available; the short cutoff is a verified
finite certificate, not a claimed short-sum theorem at arbitrary size.

Independent implementations within the checker compare all odd
membership flags through \(2^{20}\) against complete primitive
coordinate enumeration. The primitive flags agree with the union
\(C\cup43C\), and deleting ramified values recovers the core flags.
A separate literal trial-factor implementation checks the first 5000
odd flags and additional flags in low and high windows. Both bit-shift
directions are tested against literal membership. Complete small
core and primitive sumsets are computed through 4096, which includes
every claimed exception. Deliberately false summands, omitted targets,
an incorrect exception list and invalid membership bytes are rejected.

The [receipt](primitive_norm_minus43_checks.json) retains every segment's
endpoints, membership and target-mask digests, count of required targets,
number of short summands used, and coverage digest. The verifier
reconstructs the actual masks and tests every required bit; the digests
are supplementary. Only a successful full run emits a complete receipt.
The optional four-segment pilot has a different result type and an
explicitly incomplete scope.

## 4. Descent and the primitive exceptions

Suppose an additional even 43-unit core exception exists, and let
\(n\) be the least. Section 3 forces \(n>2^{32}\). If a split prime
\(p\) divides \(n\), then \(n/p\) is a smaller even 43-unit. If it
is outside \(E\), minimality gives a representation, and multiplication
by \(p\) preserves \(C\). If \(n/p=e\in E\), the analytic
exceptional-quotient theorem applies to \(n=ep\). If no split prime
divides \(n\), all its prime divisors are inert, and the analytic
inert-core theorem applies. Both alternatives contradict the choice
of \(n\).

For the excluded set, literal split-factor membership gives
\[
 C\cap[1,74]=\{1,11,13,17,23,31,41,47,53,59,67\}.
\]
Its small sumset omits precisely the ten members of \(E\) among
the eligible targets up to 74. Positive summands in any representation
of \(e\le74\) would lie in this displayed list, so those exclusions
are unconditional. This proves (M43.1).

For 43-unit targets outside \(E\), the core representations are
primitive by (M43.3). Within \(E\), the mixed representations
\[
                         56=13+43,\qquad74=31+43
\]
fill exactly two holes. The other entries below 44 cannot be mixed
sums, since the least mixed sum is \(1+43=44\). The only other entry
at least 44 is 50, and a mixed representation would require
\(50=7+43\), with \(7\notin C\). Two ramified summands have sum at
least 86. Thus the 43-unit primitive exceptions are exactly
\(E\setminus\{56,74\}\). At targets divisible by 43, the equivalence
in (M43.3) gives exactly \(43E\). Their union is the displayed
\(\mathcal E\), proving (M43.2). \(\square\)

## Evidence and scope

Run the full exact checker from the repository root:

```sh
python3 research/arithmetic/verify_primitive_norm_minus43.py \
  --output research/arithmetic/primitive_norm_minus43_checks.json
```

The analytic proof, analytic checker and analytic receipt are pinned
inside the scanner and the [source ledger](primitive_norm_minus43_source_ledger.json).
Changing a mathematical premise makes the scan command reject its
stale import rather than silently reuse it. The full run covers all
128 segments and records zero missing targets. It uses exact integers,
Python's standard library and active assertions. No floating-point
acceptance is used; timings are informational only.

This is a complete numerical-value classification. It is not a count
asymptotic, an optimized cutoff claim, or a closure of the −67/−163
cases. The ledger distinguishes the author replay, separate internal
reviews and their actual domains. No historical priority, formal
verification, publication, or external mathematical review is claimed.
