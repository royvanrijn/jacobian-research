# Bad-place support and unramified class-field blocks

The remaining bad primes reveal **joint** arithmetic structure that was absent
from the seven-place audit. Every tested individual local point image is still
filled by the marked generic subgroup. Nevertheless the known exceptional
points enlarge its joint image by two dimensions on the MW16-04 high fibre
and one on the R17 high fibre. On MW16-05 high they enlarge it by zero.

Three discriminant factorizations completed within the frozen limits, including
both fibres of the strongest scale-matched R17 pair. Their strict local kernels
give independent unramified quadratic extensions of the cubic field. The
resulting ordinary class-group 2-rank lower bounds are **10, 8 and 6**.
These are proved lower bounds from rational witnesses, not computed full class
groups, exact Selmer dimensions, or new curve ranks.

This is **incidence** evidence: it locates known rational classes within the
global arithmetic. It does not explain their simultaneous rational solubility
or supply a point-blind predictor.

## The fixed experiment

The [protocol](REMAINING_BAD_PRIMES_PROTOCOL.json) retains the three
[original high/low pairs](ANALYSIS.md), their exact independent point order and
their integral short cubics. It allows six factorization attempts of 20 seconds
and six local workers of 40 seconds. Trial factors through 10,000 are retained
before each completion attempt. A timeout leaves its residual cofactor
unfactored and its all-bad-place coverage UNKNOWN.

The computation adds only certified odd discriminant primes above 13.
Local minimal reduction removes nonminimal-model factors from the bad-place
dictionary. The retained local data at 2,3,5,7,11,13 and infinity are reused.
There are 32 new curve-prime computations; all completed. Factorization
completed on three curves and timed out on the other three. No limit was
increased, parameter added, point searched for, or class group computed.

The [local inputs](../../artifacts/generated-results/elliptic-curves/rank_jump_remaining_bad_prime_inputs_v1.json)
retain prime factorizations, incomplete cofactors, local orders, reduction
data and point squareclasses. The
[support report](../../artifacts/generated-results/elliptic-curves/rank_jump_remaining_bad_primes_v1.json)
contains both the extended dictionary and the filtered bad places plus
2 and infinity.

## Three paired comparisons

Let \(m\) be the marked generic rank and \(q\) the rank of the certified
exceptional witness quotient. Let \(S\) contain 2, infinity and the bad primes
whose factors are certified. Write \(g_S,a_S\) for the generic and full witness
image ranks in \(\prod_{v\in S}E(\mathbb Q_v)/2E(\mathbb Q_v)\).
The last column is the rank of the whole local product, not its globally
realizable subspace.

| Fibre | \(m+q\) | All bad places covered? | \(|S|\) | \(g_S\) | \(a_S\) | \(a_S-g_S\) | Local product dimension |
|---|---:|---|---:|---:|---:|---:|---:|
| MW16-05 \(307/206\) | 16+9 | yes | 14 | 15 | 15 | 0 | 16 |
| MW16-05 \(-3158/1291\) | 16+0 | UNKNOWN | 7 | 8 | 8 | 0 | 8 |
| MW16-04 \(-1647/91\) | 16+9 | UNKNOWN | 13 | 15 | 17 | 2 | 17 |
| MW16-04 \(-2177/2397\) | 16+0 | UNKNOWN | 8 | 8 | 8 | 0 | 8 |
| R17 \(-2300/843\) | 17+7 | yes | 13 | 15 | 16 | 1 | 17 |
| R17 \(-1561/3133\) | 17+0 | yes | 10 | 11 | 11 | 0 | 12 |

The low-gain rows remain censored search controls. Their full ranks are not
proved to equal \(m\). Missing factors on the three incomplete rows prohibit
full-support and class-field conclusions for those rows.

**MW16-05:** even all bad places together expose none of the nine known
exceptional quotient directions beyond the generic image. Each exceptional
class can be corrected by one generic combination to become locally square
at all these places simultaneously. The corresponding relative unramified
character block has dimension nine.

**MW16-04:** no single place exposes a new direction, but the tested joint
map exposes two. This is a real refinement of the earlier small-place result.
Unknown bad primes could increase the joint residual rank. Its seven-dimensional
tested relative kernel is not yet a certified unramified character block.

**Published R17:** both sides have complete bad-place coverage, similar
coefficient sizes and the same retained search policy. The high fibre exposes
one joint residual dimension, while six of its seven exceptional directions
can be corrected into the strict kernel. This is a complete arithmetic
decomposition of the known quotient for this dictionary; it is not a
decomposition of the entire unknown Selmer quotient.

## Joint conditions group directions, without asserting a common construction

The [strict-kernel certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_strict_class_blocks_v1.json)
stores every kernel mask and exact linear functional on the local characters.
All coordinates refer to the already certified independent input basis,
with generic points first and exceptional points \(Q_1,\ldots,Q_q\) next.

For R17 high, the condition for an exceptional word to admit a simultaneous
generic correction is precisely
\[
c_2+c_6=0.
\]
Thus \(Q_2\) and \(Q_6\) have the same nonzero joint residual image.
Their sum can be corrected into the strict kernel, as can
\(Q_1,Q_3,Q_4,Q_5,Q_7\). These are six independent relative classes.

For MW16-04 high, over its tested dictionary the conditions reduce to
\[
c_2+c_4+c_7+c_9=0,\qquad c_3+c_5=0.
\]
These conditions have exact retained witnesses in the local character
coordinates. They are not inferred from proximity of points or from a
height correlation.

The quotient map and its kernel are intrinsic to the specified subgroup and
set of places. Their displayed basis labels and chosen linear functionals
are not canonical arithmetic factors. In particular, a nonzero joint
residual image is an obstruction to **matching a generic combination**,
not an obstruction to rational solubility: all these \(Q_i\) are already
rational points. This must not be called a CT or Sha obstruction.

## From a strict kernel to unramified quadratic extensions

Let \(K=\mathbb Q(\theta)\) be the irreducible cubic algebra of a fibre, and
let \(W\subset K^\times/K^{\times2}\) be its known independent Kummer subspace.
For a set \(S\) containing every bad prime, 2 and infinity, define
\[
W^0_S=\ker\left(W\longrightarrow
             \prod_{v\in S}H^1(\mathbb Q_v,E[2])\right).
\]
Every nonzero \(\beta\in W^0_S\) defines a quadratic extension
\(K(\sqrt\beta)/K\) unramified at all places:

1. Above \(S\), the retained local squareclass is zero, so the extension
   splits. The real signs are positive, so real places stay real.
2. At every remaining rational prime \(p\), the curve has good reduction
   and \(p\ne2\). The multiplication-by-two map on its smooth integral
   model is étale, and rational point Kummer classes are unramified.
   Under the étale cubic Kummer description this says precisely that
   \(K(\sqrt\beta)/K\) is unramified above \(p\).

Independent global squareclasses give a multiquadratic compositum of degree
\(2^{\dim W^0_S}\). The global independence is inherited from the exact
point Kummer certificate and rechecked on the same independent input indices.

By the Hilbert class-field theorem, these are independent characters of
the ordinary ideal class group. They even vanish on the classes of all
prime ideals above the finite part of \(S\), hence factor through the
corresponding \(S\)-ideal class group. See
[Milne, Class Field Theory, Theorem 0.4, Chapter V Example 3.9 and Exercise 3.15](https://www.jmilne.org/math/CourseNotes/CFT.pdf).
In particular,
\[
\dim_{\mathbf F_2}\mathrm{Cl}(K)/2\mathrm{Cl}(K)
\ \ge\ \dim W^0_S.
\]

This is a statement about **class-group characters**. It does not assert that
the half ideals obtained from \((\beta)=\mathfrak a^2\) are independent in
\(\mathrm{Cl}(K)[2]\); unit classes can lie in the kernel of that different
map. No full ideal-class basis or maximal-order class-group calculation is
being smuggled into the certificate.

## Subtract the generic character block

For the generic Kummer subspace \(G\subset W\), elementary kernel accounting
gives
\[
\dim G^0_S=m-g_S,\qquad
\dim W^0_S=m+q-a_S,
\]
and therefore
\[
\boxed{\dim(W^0_S/G^0_S)=q-(a_S-g_S).}
\]
The certificate retains explicit lifts of a basis of this relative kernel.

| Completely covered fibre | Generic unramified characters | Total witnessed unramified characters | Added relative characters | Compositum degree over \(K\) |
|---|---:|---:|---:|---:|
| MW16-05 high | 1 | 10 | 9 | 1024 |
| R17 high | 2 | 8 | 6 | 256 |
| R17 low | 6 | 6 | 0 | 64 |

These are exact dimensions of the **displayed character subspaces** and
lower bounds for the full class-group 2-ranks. They do not compare exact
class-group ranks of the high and low curves. The low R17 control already
has six independent unramified characters from generic points alone,
whereas the high R17 fibre has only two from its generic subspace.

Thus raw cubic class pressure cannot be interpreted as exceptional incidence
without subtracting the marked generic contribution. Conversely, the high
fibres now provide exact relative class-field blocks, rather than merely
sharing an unspecified cubic field.

## What this changes in the mechanism ranking

1. **Incidence: a relative unramified character block is now witnessed on
   actual MW16/MW17 high fibres.** Nine and six exceptional dimensions,
   respectively, lie in this part after generic correction. This is a
   necessary arithmetic location of their rank contribution.
2. **Incidence: joint local support can carry a small remainder.** It carries
   one of seven R17 dimensions, and at least two of nine MW16-04 dimensions.
   A one-prime-at-a-time comparison misses it.
3. **Weak explanation: newly available individual local components.**
   Every tested individual point image is already filled by generic points.
   This explanation is now excluded at every bad place for the completely
   covered high MW16-05 and high R17 fibres.
4. **Still missing: global solubility of point-blind classes.** The present
   unramified characters were extracted from known points. An arbitrary
   class-group character need not satisfy the elliptic local conditions;
   even a Selmer class need not be rational rather than Sha.

The next useful experiment should construct a small **point-blind relative
unramified Selmer subspace** on the completely covered R17 pair, then test
its global obstruction rather than increase a chart budget. A full class-group
campaign is not authorized by this note. It first needs a cheap independent
source of classes or an exact bound that avoids constructing the entire
class group. Computing only generic local deficits is already point-blind,
but this experiment gives no evidence that those deficits predict the large
rational kernel.

For Agent 1, the usable conclusion is a constraint on feature design:
separate generic and relative class information, and distinguish joint
local support from individual component counts. No score or search policy
is proposed or changed.

## Replay

From the repository root:
```sh
python3 elliptic-curves/rank-jump/remaining_bad_primes.py check
python3 elliptic-curves/rank-jump/strict_class_blocks.py check
sage -python elliptic-curves/rank-jump/verify_remaining_bad_primes.py --index 0
```
Use indices 0 through 5 for the six local replays. The arithmetic replay
verifies all certified factors by primality and product, recomputes the new
local signatures and minimal reductions, and checks every retained strict
kernel generator by PARI's separate local-power interface. It also checks
positivity at every real embedding. All six passed, totalling 1,048
independent local-power checks. No integer factorization, point search or
class-group calculation is needed for replay.

Construction used Sage 10.9 and PARI 2.17.3. Original factor timeouts and
per-curve logs are retained under the ignored
`artifacts/local/rank-jump-remaining-bad-primes-v1` directory.
Incomplete rows stay incomplete on replay.
