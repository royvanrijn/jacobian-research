# MW16 at 3/17: targeted 2-descent preparation

Mathematical status is recorded in `MATH_STATUS.json` under
`EC-SMALL-CONDUCTOR-DESCENT-SHORTCUT-20260906` and
`EC-SMALL-CONDUCTOR-NORM-BATCH-20260906`, and
`EC-SMALL-CONDUCTOR-SPECIAL-PRIMES-20260906`, and
`EC-SMALL-CONDUCTOR-SMALLER-BASE-20260906`. The
[curve certificate](NEW_SMALL_CONDUCTOR_CURVE_2026-09-05.md) still proves
rank **at least 22**, with no unconditional rank upper bound. This study specializes a
known class-group criterion and tests better norm coordinates; it does not
complete a 2-descent or exhibit another point.

## What was already reached

The [relative class-group retry](R17_MW29_RELATIVE_SCLASS_RETRY_2026-09-04.md),
[BNF-free residual machinery](BNF_FREE_RESIDUAL_2SELMER.md), and
[small-field laboratory](R17_SMALL_FIELD_CLASS_QUOTIENT_LAB.md) already
identify the difficulty: certified global generation and productive principal
relations. Factor hints, maximal-order caches, mod-2 linear algebra and
relative/local descent are existing machinery. Recommending those again
would not resolve the obstruction. The earlier 180-second `ellrankinit`
attempt on this particular curve also returned no upper bound.

The new diagnostic isolates setup from relation collection. All arms use
one worker, a 1.5 GiB RSS limit and a 512 MB PARI stack. The field arm has a
30-second cap; the three class-group arms have 30, 30 and 60 seconds.
PARI 2.17.3, scripts, logs, supervisor records and hashes are retained under
`artifacts/local/elliptic-curves/small-conductor-descent-profile-v1/`.

| Arm | Observed field setup | Terminal result |
| --- | ---: | --- |
| Factored discriminant, maximal order only | 1 ms; certification 1 ms | Complete; 0.081 seconds including process overhead |
| Cold `ellrankinit` | 6,543 ms | 30-second timeout |
| Factor-hinted `ellrankinit` | 4 ms | 30-second timeout |
| Cached maximal order, `bnfinit` flag 0 | 0 ms | 60-second timeout |

Termination grace accounts for about one extra second in timeout wall times.
The three class arms each use 3,878 factor-base ideals and remain in relation
collection. No BNF certificate is reached. The old profile's `gp_error` boolean
matches normal debug text containing `***`; timeout outcomes and retained logs,
not that boolean, determine the interpretation.

## Exact stopping target

Write `g = dim_F2 Cl(K)[2]`, where the cubic 2-division field is

```text
K = Q[z] / (z^3 + z^2 - 2919231625641258502793755607986240*z
             + 45440201616242830029801770634418828098464819545088).
```

Its maximal-order discriminant is

```text
128900477062442043600727490102612931938219670661531295245188752203875468
```

The Brumer–Kramer bound, in the notation of
[Klagsbrun–Sherman–Weigandt, Proposition 3.1](https://arxiv.org/html/1606.07178#S3.SS1),
specializes to

```text
22 <= rank E(Q) <= dim_F2 Sel_2(E/Q) <= g + 7.
```

Here `u=2`. The multiplicative primes with even minimal discriminant valuation
are `2,5,13,19,71`, contributing five. The sole additive prime is 17 and has
one prime above it in K, contributing zero. An exact Newton polygon of
`F(z+11)` has vertex endpoints `(0,2),(3,0)` and coefficient valuations
`[2,2,1,0]` at 17. Its slope proves irreducibility and total ramification of
the cubic over Q_17, independently of the number-field calculation.

The global root number is +1 and rational 2-torsion is trivial. The proved
2-parity theorem makes the **2-Selmer dimension even**; see
[Dokchitser–Dokchitser, Theorem 1.4](https://annals.math.princeton.edu/2010/172-1/p11)
and its explicit finite-Selmer formulation in
[Klagsbrun–Sherman–Weigandt, Appendix A](https://arxiv.org/html/1606.07178#A1).
This is not an assertion of unconditional algebraic rank parity.

Consequently `g >= 15` unconditionally, and **any certified `g <= 16` implies
exact rank 22**: the Selmer bound is at most 23, hence at most 22 by parity.
Conversely rank at least 23 would force Selmer dimension at least 24 and
therefore `g >= 17`. A future GRH-dependent class-group upper bound would
yield a GRH-dependent rank statement. No such upper bound currently exists.

## Productive norm coordinates

Constructing a binary cubic from the maximal-order multiplication table and
reducing its positive Hessian gives

```text
F_red(m,n) = -582653777472000*m^3
            -69418898236910689*m^2*n
            +2253317262489553494*m*n^2
            +28630981474137424362*n^3.
```

The maximum coefficient size falls from 165 to 65 bits. This uses a
nonmonic binary form: `polredabs` on the original monic polynomial merely
changes the generator's sign. Binary-cubic reduction is an established
number-field sieve technique, also used in the cited high-rank computations;
the new result here is its exact implementation and measured yield on this
specific field.

The certificate retains the integral basis, multiplication table, determinant-one
coordinate matrix M, and element w. With `a=-582653777472000`, set

```text
(u,v) = M * (m,n)
beta = a*u + v*w
Norm(beta) = a^2 * F_red(m,n).
```

The fixed square factor **a² is included** in every ideal relation. Discriminant
identities alone do not certify maximality: a separate `nfcertify` replay
does that. The pure rational checker verifies the complete polynomial norm
identity and the exact coordinate transport.

Before execution, the pilot fixes all 5,039 coprime pairs with
`-64 <= m <= 64`, `1 <= n <= 64`, and smoothness bound 400,000. Repeated gcd
with the exact primorial removes supported factors; unresolved remainders
are saved without expensive factorization. Each arm is capped at 120 seconds
with one worker and 1.5 GiB; checkpoints occur every eight m values.

| Norm polynomial | Completely smooth values | Observed arithmetic time |
| --- | ---: | ---: |
| Original monic cubic | 0 / 5,039 | 4.504 seconds |
| Reduced maximal-order binary cubic | 18 / 5,039 | 2.378 seconds |

These are different, explicitly specified sets of principal elements, so the
comparison measures this finite choice of coordinates, not a universal speedup.
An exact ideal-factorization audit verifies all 18 principal relations,
including the fixed norm factor. On the displayed finite ideal columns their
parity rows add **18 independent rows modulo the canonical rational-prime
relations**. This means additional relations, not independent class-group
elements. The columns are not proved to generate the class group modulo squares.

## Remaining gate and reuse

Preparation and a productive small pilot now take seconds. Full 2-descent
runtime is still unmeasured. Before a larger run, choose and justify a generating
factor base for `Cl(K)/2`, then feed these exact norm relations into the existing
mod-2 engine. With generation established, enough relations to leave matrix
nullity at most 16 suffice; a complete class-group structure is unnecessary.
If the residual bound is larger, more relations or tighter local Selmer
conditions may still be needed. A larger bound does not prove another point.

An unconditional generating-set certificate remains the principal proof gate.
A GRH-based bound must be explicitly labelled as conditional throughout.
The sieve bound 400,000 in this pilot is only a finite smoothness cutoff, not
a proved generation certificate. The initial pilot launched no larger sieve
or sweep of previous curves. Reuse should proceed only after a further bounded batch shows
sustained independent relation yield on a justified factor base; each new
field also needs its own arithmetic and local correction terms.

## Replay

The aggregate checker replays the original rank/conductor proof, theorem
specialization, root numbers, integral norm identity, every pilot value, and
all principal-ideal identities and parity ranks. It checks historical profile
hashes rather than rerunning the timed-out class jobs.

```bash
sage -python elliptic-curves/cas/check_small_conductor_descent_shortcut.sage --check
```

Inputs and outputs are retained in the
[aggregate certificate](../../artifacts/generated-results/elliptic-curves/small_conductor_descent_shortcut_v1.json)
and its source bindings. The
[evidence archive](../../artifacts/generated-results/elliptic-curves/small_conductor_descent_shortcut_evidence_v1.zip)
contains the replay inputs; the
[isolated replay report](../../artifacts/generated-results/elliptic-curves/small_conductor_descent_shortcut_portable_replay_v1.json)
records verification from a fresh extracted directory. Replay requires Sage
with PARI 2.17.3 and the pinned `/usr/bin/gp` used for root numbers.

## Authorized 512-box continuation

The subsequent user-authorized batch freezes the same field and smoothness
cutoff, and every coprime pair with `-512 <= m <= 512`, `1 <= n <= 512`.
It includes the original box for an exact incremental comparison. The
protocol caps the worker at 180 seconds, scalar replay at 300 seconds,
and each ideal-matrix audit at 180 seconds, with one worker and 1.5 GiB RSS.
Completed chunks of 16 consecutive m values and progress records are retained
under `artifacts/local/elliptic-curves/small-conductor-norm-batch-v1/`.

The worker uses an exact product/remainder tree to compute the primorial
modulo all norms in a chunk, then removes all supported prime powers by gcd.
This is an established batch arithmetic technique. It processes **319,407
pairs in 18.117 seconds**, finding **296 smooth norms**. A separate scalar
implementation recomputes every value and smoothness remainder without the
tree, passing in **137.946 seconds**. The approximately 7.6-fold difference
is for these complete implementations on this fixed input, including their
respective record handling; it is not a calibrated general speedup claim.

All **296 principal ideal relations** verify exactly and contribute independent
mod-2 rows beyond the canonical rational-prime relations. Thus there are
**278 additional rows** beyond the previous 18. The finite smoothness rate is
296/319,407, about 0.093%, compared with 18/5,039, about 0.357%, in the initial box.
The smaller pilot's rate should not be extrapolated to larger coordinates.

This time the auditor enumerates **all 62,003 prime ideals above the 33,860
rational primes at most 400,000**. Each complete prime decomposition and
canonical principal product is checked. Using the existing interval-arithmetic
generation-bound routine gives an upper endpoint of **321,720** for
`12 log(|disc K|)^2`. Therefore this full base generates the ideal class group
under the GRH hypothesis of Bach's theorem; it is not certified to generate
unconditionally. Including ideals of norm greater than the cutoff above these
rational primes only enlarges the base and does not invalidate generation.

| Included box | Verified norm relations | Additional row rank | Full quotient dimension |
| --- | ---: | ---: | ---: |
| 64 | 18 | 18 | 28,125 |
| 128 | 38 | 38 | 28,105 |
| 256 | 98 | 98 | 28,045 |
| 512 | 296 | 296 | 27,847 |

The 33,860 canonical rational-prime rows are independent because their
supports are disjoint and each has an odd exponent. Before adding norm
relations, the quotient dimension is 28,143. The final **27,847** is a
GRH-dependent upper bound for `g`, not its computed value. The corresponding
coarse curve-rank upper bound is 27,854 under GRH, far from the desired 22.
For this particular generating base and criterion, at least **27,831 further
independent rows** would be needed to reach quotient dimension 16.

The larger test therefore establishes productive and faster relation
collection, but it does **not** finish a useful rank bound. A further campaign
should target remaining factor-base directions and justify a smaller base
where possible, with a declared batch and residual-dimension stopping rule.
It should measure independent row gain, since smooth-norm count alone will
eventually overstate progress. No further batch or previous-curve sweep is
launched by this continuation.

The [batch certificate](../../artifacts/generated-results/elliptic-curves/small_conductor_norm_batch_v1.json)
binds every chunk. The
[ideal-matrix certificate](../../artifacts/generated-results/elliptic-curves/small_conductor_norm_batch_relations_v1.json)
contains the full ideal columns, principal factorizations and sparse rows.
The [continuation evidence archive](../../artifacts/generated-results/elliptic-curves/small_conductor_norm_batch_evidence_v1.zip)
includes the preceding proof bundle, the new sources and all batch records.
Its [isolated matrix replay](../../artifacts/generated-results/elliptic-curves/small_conductor_norm_batch_portable_replay_v1.json)
rebuilds the full factor base and verifies every principal relation. The scalar
smoothness run is retained with its hashes and is not repeated by that matrix
replay.

```bash
sage -python elliptic-curves/cas/audit_small_conductor_norm_batch.sage --check
```

The checker reuses `bounds` and `packed_rank` from
`audit_bnf_free_s_class_quotient.py`; it does not create a separate class-group
engine. Neither the GRH hypothesis nor the large residual is suppressed in
the certificate or status entry.

## Targeting 64 previously untouched prime ideals

The next user-authorized continuation targets missing factor-base directions
explicitly. After eliminating the canonical rational-prime rows, the existing
296 norm relations have nonzero entries on 1,449 quotient columns. The target
selector excludes these covered columns and the canonical pivots. It keeps
unramified degree-one ideals over primes at least 1,000, excludes divisors of
the field discriminant, defining-order index and fixed norm coefficient a,
and chooses one eligible ideal per rational prime. This leaves 21,595 eligible
distinct rational primes. The protocol freezes 64 equally spaced indices in
that ordered list, including endpoints, spanning **1,021 through 399,989**.

For a selected ideal `(q, theta-r)`, the exact map
`beta = a*(M00*m+M01*n) + w*(M10*m+M11*n)` gives a linear condition
`m = s*n mod q`. The script verifies the ideal HNF against the corresponding
root, computes s exactly, and constructs the index-q lattice with basis
`(q,0),(s,1)`. Gauss reduction in the positive binary-cubic Hessian metric
retains that lattice while reducing its norm sizes. Every resulting candidate
therefore belongs to the intended ideal. This is the established special-prime
relation-sieving method; see the divisibility construction in
[Klagsbrun–Sherman–Weigandt, section 4.4](https://arxiv.org/html/1606.07178#S4.SS4).

For each of the 64 reduced lattice bases, the fixed box is all coprime
coefficient pairs `-64<=u<=64`, `1<=v<=64`; mapped pairs `(m,n)` are retained
only when primitive. Occurrences across targets remain separately recorded.
The smoothness cutoff remains 400,000. The protocol caps the worker at
180 seconds, independent scalar replay at 300 seconds and each ideal audit
at 180 seconds, using one worker and 1.5 GiB. Every target has its own retained
chunk and a completed-target checkpoint under
`artifacts/local/elliptic-curves/small-conductor-special-primes-v1/`.

The worker processes **322,466 candidate occurrences in 19.253 seconds** and
finds **536 smooth norms**, with at least one for every target. The scalar
replay checks every map, norm and smoothness remainder in **150.649 seconds**.
The target ranges show that the large-prime directions remain productive:

| Target prime range | Targets | Candidate occurrences | Smooth occurrences |
| --- | ---: | ---: | ---: |
| 1,021–90,947 | 16 | 80,608 | 198 |
| 96,847–189,479 | 16 | 80,616 | 130 |
| 196,081–292,661 | 16 | 80,619 | 108 |
| 299,623–399,989 | 16 | 80,623 | 100 |

The exact ideal audit verifies all 536 principal products, including the fixed
norm factor a² and positive valuation at each intended target ideal. All 536
pairs are distinct up to sign, and **all 536 rows are independent modulo the
entire previous relation matrix**. The projected matrix on the 64 selected
columns has rank **64**; these were zero columns in the previous canonically
reduced norm matrix. Thus the new relations independently reach every selected
direction, beyond merely containing a forced rational prime factor.

The cumulative noncanonical row rank is now **832**, and the full quotient
dimension falls from **27,847 to 27,311**. On the same factor base this leaves
**27,295 additional independent rows** to reach the sufficient class-2-rank
target of 16. The corresponding coarse upper bounds are `g<=27311` and
`rank E<=27318`, conditional on the same GRH generation hypothesis. The
unconditional rank upper bound and exact rank remain unknown.

This finite trial supports targeting untouched ideal directions. It does not
establish a sustainable rate for a full campaign: this selection uses one
ideal per rational prime, and later residual directions may be more difficult.
The remaining eligible ideals, including other ideals over the same primes,
are not exhausted. No broader run is launched by this 64-target protocol.

The [targeted batch certificate](../../artifacts/generated-results/elliptic-curves/small_conductor_special_primes_v1.json)
binds the frozen selection, every lattice basis and every retained norm chunk.
The [new relation certificate](../../artifacts/generated-results/elliptic-curves/small_conductor_special_prime_relations_v1.json)
records the exact generators, ideal factorizations, rank gain and projection
rank. The [targeted evidence archive](../../artifacts/generated-results/elliptic-curves/small_conductor_special_primes_evidence_v1.zip)
includes both earlier proof bundles and the new records. Its
[isolated replay report](../../artifacts/generated-results/elliptic-curves/small_conductor_special_primes_portable_replay_v1.json)
checks the deterministic selection and lattice construction, rebuilds the
previous full factor-base matrix, then verifies all new ideal identities and
rank calculations. The completed scalar replay is retained with its hashes;
the isolated matrix check does not repeat the scalar run.

```bash
sage -python elliptic-curves/cas/target_small_conductor_prime_ideals.sage audit-check
```

## Route 1: certify a smaller generating base

The user requested the smaller generating-base route first, with direct
Selmer work as fallback. This route succeeds in reducing the certified
GRH-dependent class-group problem: the cutoff **37,638** is now justified
independently of the earlier timed-out BNF run, and the valid residual dimension
is **3,210**, compared with 27,311 on the old base. Exact rank remains unknown.
No further norm sieve, point search or direct Selmer computation is launched
by this finite generation-and-elimination audit.

The protocol freezes cutoffs `37638,40000,50000,75000,102772`, 256-bit
interval arithmetic, 4,096 series terms, one worker, 1.5 GiB RSS and a
120-second limit per local stage. V1 stops before evaluating any generation
test because the Sage MPFI field has no `catalan_constant` method. Its source,
protocol, log and supervisor failure are retained. V2 uses Arb's rigorously
enclosed Catalan constant converted to an MPFI interval, with unchanged
cutoffs, budgets and mathematical test.

### Explicit-formula certificate

Let `L=log(T)` and `F(x)=max(1-x/L,0)` for nonnegative x. Its even extension
is the normalized convolution square of the characteristic function of
`[-L/2,L/2]`, so its Fourier transform is nonnegative. With field degree
`n=3` and `r1=3`, compute

```text
S = sum_P log(NP) * sum_{m>=1} F(m*log(NP))/(NP)^(m/2)
I = integral_0^infinity (1-F(x))/(2*sinh(x/2)) dx
J = integral_0^infinity F(x)/(2*cosh(x/2)) dx
margin = 2*S - log(abs(disc K)) + 3*(EulerGamma + log(8*pi)) - 3*I + 3*J.
```

The prime sum is finite: only prime powers of norm below T contribute.
All contributing prime ideals are already present in the independently
replayed 400,000 base. A strictly positive margin proves generation by prime
ideals of norm below T under GRH for the nontrivial ideal-class characters,
by [Grenié–Molteni (2017), Theorem 2.4](https://arxiv.org/pdf/1607.02430),
restating the theorem of Belabas–Diaz y Diaz–Friedman.

The archimedean integrals use the exact expansions, with `lambda=k+1/2`,

```text
I = (pi^2/2 - sum_{k>=0} exp(-lambda*L)/lambda^2)/L
J = pi/2 - 4*Catalan/L
    + sum_{k>=0} (-1)^k*exp(-lambda*L)/(L*lambda^2).
```

After N terms, the absolute tail of either exponential series is at most
`exp(-(N+1/2)*L)/((N+1/2)^2*(1-exp(-L)))`. The implementation encloses this
tail and all arithmetic in intervals. At **T=37,638**, the certified margin is

```text
0.01280065990512546338470769612035110997184938609287673150357607764823574164517
  <= margin <=
0.01280065990512546338470769612035110997184938609287673150357607764823574888485.
```

Thus the positivity has been checked independently; the old PARI debug cutoff
is not being promoted to a proof by itself. As a separate comparison,
[Grenié–Molteni (2025), Theorem 1](https://doi.org/10.1090/mcom/4114) gives the
GRH bound `(23/6)*log(abs(disc K))^2` because the cubic's log discriminant is
about 163.7374, below the theorem's threshold 353. Its interval upper endpoint
rounded up is 102,772. The field-specific inequality is substantially stronger
here. The cutoff 37,638 is the smallest in the frozen test list; it is not
claimed to be optimal.

### Keeping only valid relations on the smaller base

The complete old matrix and all 832 principal relations replay first. After
canonical rational-prime elimination, write each norm row as its small-prime
and outside-prime components. Gaussian elimination on the outside components
tracks the exact XOR of original row indices. Only combinations whose outside
component vanishes are admitted to the smaller matrix. Each saved combination
is reconstructed and checked; no outside coordinate is simply discarded.

| Cutoff | Positive generation margin | Initial quotient dimension | Valid supported relation rank | Remaining dimension |
| --- | ---: | ---: | ---: | ---: |
| 37,638 | 0.01280066 | 3,243 | 33 | 3,210 |
| 40,000 | 3.65536207 | 3,427 | 38 | 3,389 |
| 50,000 | 17.85631108 | 4,194 | 54 | 4,140 |
| 75,000 | 47.46493367 | 6,060 | 97 | 5,963 |
| 102,772 | 74.43881665 | 8,120 | 163 | 7,957 |

The 37,638-base outside projection has rank 799, leaving an intersection of
rank `832-799=33` with the small-prime coordinate space. This reduces the
certified GRH-dependent class-2-rank upper bound to **3,210** and the coarse
curve-rank upper bound, including Selmer parity, to **3,216**. On this base,
**3,194 further independent supported relations** would suffice for the
class-2-rank target of 16. Neither this numerical gap nor successful generation
proves that the target is attainable: a larger actual class-2-rank could
require tighter local Selmer conditions to determine the curve rank.

The [smaller-base certificate](../../artifacts/generated-results/elliptic-curves/small_conductor_smaller_base_v2.json)
contains every interval and every surviving relation combination. The
[evidence archive](../../artifacts/generated-results/elliptic-curves/small_conductor_smaller_base_evidence_v2.zip)
includes the prior evidence, both primary-source PDFs and the retained V1
failure. Its [isolated replay report](../../artifacts/generated-results/elliptic-curves/small_conductor_smaller_base_portable_replay_v2.json)
recomputes the inherited matrix, the explicit-formula inequalities and the
relation intersections from a fresh extracted directory.

```bash
sage -python elliptic-curves/cas/certify_small_conductor_smaller_base_v2.sage --check
```

## Continuation: 512 targets directly on the smaller base

The continuation keeps the certified cutoff at **37,638** and starts from
the 33 exactly supported relations above. Its mathematical gate is to reduce
the remaining 3,210-dimensional quotient using relations entirely supported
on this base. The factor-base generation statement remains conditional on
GRH; exact principal-ideal identities and finite matrix ranks do not.

The frozen protocol chooses 512 of 2,428 eligible prime ideals, at equally
spaced indices in the existing column order, with rational primes from 1,009
through 37,633. It excludes canonical pivots, all columns touched by the
33 supported relations, repeated rational primes, ramified or higher-degree
ideals, and divisors of the defining-order index or fixed norm multiplier.
Each target uses an exactly checked Hessian-reduced index-prime lattice.
The fixed box is `-64 <= u <= 64`, `1 <= v <= 64`, with both `(u,v)` and
the mapped `(m,n)` primitive. Each completed target receives a checkpoint.

The one-worker limits are 300 seconds for relation generation, 300 seconds
for scalar replay and 180 seconds for ideal and matrix auditing, each with
1.5 GiB RSS. Smoothness uses only rational primes at most 37,638. The scalar
replay independently checks every candidate against the full primorial;
the ideal audit reconstructs `beta`, verifies `Norm(beta)=a^2*F(m,n)`,
factors its principal ideal and checks the product exactly. It rejects any
factor outside the smaller base. New row rank is measured modulo the 33
existing supported rows, after canonical rational-prime elimination.

The batch completes **2,579,492 candidate occurrences in 35.21 seconds**.
It finds 383 smooth occurrences, comprising 331 pairs up to sign and 52
duplicates. **All 331 distinct relations add independent rows** modulo the
old matrix. There are smooth hits in 262 of the 512 target boxes. After
canonical elimination, the complete new matrix touches 313 of the selected
columns and its projection onto those columns has rank 276; these are
different measurements from either successful boxes or total row gain.
The independently implemented scalar replay passes in approximately
119 seconds, and the exact ideal and matrix audit completes within
4.1 seconds including supervisor overhead.

| Quantity | Before | After |
| --- | ---: | ---: |
| Supported principal-relation rank | 33 | 364 |
| Remaining quotient dimension | 3,210 | 2,879 |
| Further independent rows sufficient to reach dimension 16 | 3,194 | 2,863 |

Thus this finite continuation removes **10.3%** of the previous residual
dimension. It gives the GRH-conditional bound `dim Cl(K)[2] <= 2879`, hence
the coarse GRH-conditional curve-rank bound 2,886 after Selmer parity. The
unconditional lower bound remains 22; exact rank is still unknown. This is
progress in relation collection, not a completed 2-descent. The observed
yield does not certify a runtime for finishing the matrix or guarantee that
dimension 16 can be reached. Further batches would need to account for the
newly covered directions and increasing dependencies.

The [search certificate](../../artifacts/generated-results/elliptic-curves/small_conductor_small_base_targets_v1.json)
pins the complete candidate chunks and scalar result. The
[ideal and matrix certificate](../../artifacts/generated-results/elliptic-curves/small_conductor_small_base_relations_v1.json)
contains all exact principal-ideal factorizations and the verified rank gain.
The [portable evidence archive](../../artifacts/generated-results/elliptic-curves/small_conductor_small_base_targets_evidence_v1.zip)
includes the inherited proofs and every new checkpoint. Its
[isolated replay report](../../artifacts/generated-results/elliptic-curves/small_conductor_small_base_targets_portable_replay_v1.json)
records reconstruction of the certified smaller base, supported old matrix,
target selection and all new ideal relations from an extracted directory.
It hash-checks the retained candidate chunks and locally passed scalar
report without repeating the scalar search.

```bash
sage -python elliptic-curves/cas/target_small_conductor_small_base.sage audit-check
```
