# Full-extension calibration V2

Result: [independently replayed experiment and certificates](../../artifacts/generated-results/elliptic-curves/adaptive_visibility_cascade_v2.json),
SHA256 `66f3f6639ea7b7ef8a2653d7dbbfa13377a41d0ab719f544d1bd25bf39ece925`.
All 906 certificate-to-transcript bindings and all 1,169 recorded V1
preservation hashes pass. Every landscape's metric, coverage, shortlist and
exact-CVP refinement is replayed; each chart's exact witnesses and mod-2 rank
audit are checked. Every gain also passes mod-3/mod-5 checks, and the terminal
no-gain cloud passes the separate post-run mod-3/mod-5 diagnostic.

**The frozen V2 run reaches 30, not 31.** It certifies
`17 → 18 → 19 → 20 → 21 → 22 → 23 → 24 → 25 → 26 → 27 → 28 → 29 → 30`,
then completes a `30 → 30` no-gain epoch and stops. All **906** charts complete
without timeout. Across 14 landscapes it scores **524,256** cosets, refines
**5,600** with exact rounded-metric CVP, and cancels **2,834** stale chart slots
after gains. The final cloud also has certified finite column rank 30 modulo
both 3 and 5; this post-run diagnostic did not expose a missed extra direction.
These are lower-bound certificates and a bounded experimental stop, not an
upper bound on the curve or proof of dependence of every other returned point.

The target of autonomous rank 31 under this rule is **not achieved**. The
original stop, anchor bank, shortlist sizes, coordinate budget and score have
not been amended. No V2 continuation or residual-labelled repair is folded
into the result.

| Input → output | Masks per anchor | Charts executed | Stale slots cancelled |
| --- | ---: | ---: | ---: |
| 17 → 18 | 1 | 17 | 65 |
| 18 → 19 | 2 | 37 | 77 |
| 19 → 20 | 4 | 1 | 177 |
| 20 → 21 | 8 | 22 | 284 |
| 21 → 22 | 16 | 31 | 275 |
| 22 → 23 | 32 | 81 | 225 |
| 23 → 24 | 64 | 28 | 278 |
| 24 → 25 | 128 | 8 | 298 |
| 25 → 26 | 256 | 41 | 265 |
| 26 → 27 | 512 | 30 | 276 |
| 27 → 28 | 1,024 | 48 | 258 |
| 28 → 29 | 2,048 | 127 | 179 |
| 29 → 30 | 4,096 | 129 | 177 |
| 30 → 30 | 8,192 | 306 | 0 |

The five gains from 25 through 30 use extension encodings 190, 331, 815,
729 and 1702, respectively. All are outside a first-64 coordinate prefix **in
their own displayed bases**. These labels were recorded after selection, not
used as priorities. This demonstrates actual use of formerly truncated parts
of these landscapes; it is not a controlled claim that changing only V1's
cutoff causes all the gains, since V2 also changes rebuilding, tie handling
and CVP refinement. Subgroups of equal rank in V1 and V2 are not identified.
Full mask scoring does not imply exhaustive chart coverage: the retained
anchor bank, heuristic shortlist, choice among minimum representatives, and
finite chart budget remain possible sources of the final miss. No missing
point was used to decide among these explanations.

V2 is a separate experiment; V1's frozen scripts, inputs, charts and
certificates are preserved. V2 starts again with the generic MW17 specialization
on curve 302. Its designer knows the earlier calibration, but its guarded
worker receives no exceptional point, residual label, preferred orbit ID, or
earlier search checkpoint. This is calibration, not an out-of-sample experiment.

Execution is implemented in
[`adaptive_visibility_cascade_v2.sage`](../cas/adaptive_visibility_cascade_v2.sage),
with the exact lattice routine in
[`visibility_lattice_v2.py`](../cas/visibility_lattice_v2.py) and read-only replay in
[`check_visibility_cascade_v2.sage`](../cas/check_visibility_cascade_v2.sage).
The immutable local protocol is
`artifacts/local/elliptic-curves/adaptive-visibility-cascade-v2/protocol.json`,
SHA256 `2d6e1458017b871b063af51b09474f81ae6023ba450a679f3a785b8b54b5de38`.

## Frozen scheduling rule

At every certified subgroup, retain 16 generic anchors in each base shell 8
and 10, ordered by descending specialized canonical-word norm. Each anchor
receives **every** one of the `2^(r-17)` extension masks. A vectorized LLL/Babai
pass scores them all. Per anchor, refine the greatest 16 Babai scores with
exact CVP and retain up to eight untested centres by descending exact norm.
These are interleaved, by norm, with the next 25 untested canonical centres
per shell. Ties use actual finite-reduction coset fingerprints, then lane and
exact centre coordinates—not the integer encoding of the extension mask.
Actual centres are deduplicated up to sign.

Only selected chart finalists receive factor-free quartic construction. Each
chart uses height 125,000 and at most ten seconds. Its exact point witnesses
are replayed and its accumulated point cloud audited immediately. A certified
gain cancels all remaining charts on the old subgroup. Independent mod-2 and
mod-3/mod-5 replay precedes the next landscape. The limits are 4,096 charts,
20 epochs, and a stop at rank at least 31 or the first complete no-gain epoch.
A CVP node-budget failure aborts; it is not interpreted as a poor coset score.
Timeouts are censored searches, not absence certificates.

Every full landscape is retained, including residues, transported residues,
Babai words and integer norms. Exact-CVP records include all minimizing words
and enumeration node counts. Individual chart transcripts, accumulated clouds,
independent certificates and counts of cancelled stale charts are retained.

## What full extension enumeration proves

Suppose the displayed independent generators give a literal direct sum
`M_r = M_17 ⊕ ZQ_18 ⊕ … ⊕ ZQ_r`. The inverse image of any fixed base parity
under the corresponding projection `M_r/2M_r → M_17/2M_17` has exactly
`2^(r-17)` elements. Enumerating every binary extension covers that fibre once:
two masks differ by an element of `2M_r` precisely when their bits agree.
Changing coordinates while transporting this same fibre changes its labels,
not its elements. No saturation claim about the ambient Mordell–Weil group is
needed. Changing the splitting can change what is meant by a fixed base fibre;
the basis-prefix requirement fixes that splitting in this experiment.

For the tie key, reduce actual points into the finite groups `E(F_p)/2E(F_p)`
at the fixed good primes below 1,000. Their direct product gives a homomorphism
on `M_r/2M_r`. The implementation checks that its columns have rank `r`, hence
that the fingerprint is injective on this finite parity space. Its value for
a coset is independent of how that coset is expressed in the displayed basis.
Reversing mask enumeration must reproduce the identical ordered shortlist.

There is a further exact consequence of the fingerprint certificate. If
`2Q ∈ M_r`, its finite mod-2 fingerprints vanish. Injectivity implies
`2Q = 2m` for some `m ∈ M_r`. The separately certified absence of rational
2-torsion gives `Q = m`. Thus each certified subgroup is **2-saturated** in
`E(Q)`. This does not establish saturation at odd primes, equality with the
whole Mordell–Weil group, or an upper rank bound. In particular, a remaining
independent direction cannot be a literal rational half of a point of `M_r`.

This does **not** prove invariance of the entire heuristic under arbitrary
rebasing. LLL and Babai may select different feasible representatives, and
rounding a freshly computed height matrix need not commute with rebasing.
The full set and semantic ties remove the coordinate-prefix bias of V1; they
do not turn a Babai shortlist into a canonical quotient optimizer.

## Replay

The main replay rederives the anchor bank, full extension coverage, shortlist,
exact minima, semantic tie choices, chart witnesses, and per-chart rank audits.
A separate
[`metric replay`](../cas/check_visibility_metric_v2.sage) recomputes every height
decision form and LLL transport from the certified points. The
[`finalizer`](../cas/finalize_visibility_cascade_v2.py) additionally binds each
rank certificate's point cloud to the exact accumulated search transcripts,
checks disjointness of anchor fibres, verifies the V1 preservation manifest,
and packages the compact generated result. Verified replay prefixes are reused
only when their checker, protocol and all checkpoint hashes still match.

From the repository root, after the worker finishes:

```sh
sage -python research/elliptic-curves/cas/check_visibility_cascade_v2.sage --progress research/artifacts/local/elliptic-curves/adaptive-visibility-cascade-v2/replay-progress.json --output research/artifacts/local/elliptic-curves/adaptive-visibility-cascade-v2/replay.json
sage -python research/elliptic-curves/cas/check_visibility_metric_v2.sage
python3 research/elliptic-curves/cas/finalize_visibility_cascade_v2.py
```

The finalizer refuses to overwrite an existing packaged result. Raw checkpoints
are retained under the local V2 directory; the generated certificate hashes
them. Sage 10.9/PARI library 2.17.3 constructs the maps, while the separately
hash-pinned external GP 2.15.4 executes the bounded searches. The immutable
worker protocol predates the first chart; adding independent replay checks
does not amend its selection policy.

## Exact-CVP claim and proof boundary

The decision form is the positive definite integer matrix obtained by rounding
the 384-bit canonical-height matrix at scale `10^6`. It is not the exact
canonical-height pairing. The LLL change of basis is checked unimodular and
this integer form is transported exactly.

For its rational LDL decomposition, write

`G[w] = Σ_i d_i (w_i + Σ_(j>i) μ_(j,i) w_j)^2`, with all `d_i > 0`.

An exactly evaluated Babai word supplies a feasible initial radius. Descending
recursion fixes coordinates with the prescribed parity. Each allowed interval
is computed by rational arithmetic and an integer square root. Positivity
makes every excluded branch exceed the current feasible radius. The finite
closed ellipsoid contains every possible improvement and tie. Exhausting it
therefore proves the returned minimum in the rounded metric and lists all
minimizers. A node-budget exception supplies no certificate. Twenty small
parity problems are independently compared against explicit box enumeration,
including multiple-minimum cases.

Only shortlisted cosets receive this exact computation. A Babai norm is an
upper bound on that coset's minimum; ranking these upper bounds is heuristic.
Consequently full mask scoring does not prove that the best chart or deepest
coset survives the shortlist, and a terminal no-gain epoch cannot prove that
the subgroup contains every rational point.
