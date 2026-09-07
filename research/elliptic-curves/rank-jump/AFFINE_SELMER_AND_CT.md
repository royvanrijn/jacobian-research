# The transporter adds Selmer classes outside the inherited space

The common shift in the [linear-twist transporter](LINEAR_TWIST_SOLUBLE_BLOCKS.md)
exposes a previously untested part of the fixed-field descent. Its affine
coset meets the Selmer group at `u=-1,2,3` and is excluded at
`u=-3,-2,1`. The three nonempty intersections each add one independent
Selmer class. Four additional exact CT pairings show that the enlarged
restricted CT radicals have dimensions **3, 2 and 2**, respectively.

At `u=-1` the new class is represented by an explicit rational point,
already available through the negative-square section construction.
At `u=2,3` the new radical directions remain **point-or-Sha UNKNOWN**.
The full Selmer groups and curve-rank upper bounds are not computed.

## What changed in the descent space

Keep the previous fixed cubic `K=Q(θ)`, the twenty-dimensional
independent anchor space `W`, and the complete local intersections
`W_u=W∩Sel_2(E_u)`. The transporter gives classes in

\[
 \eta_u+W,\qquad \eta_u=[D(u)(1-u\theta)]
\]

using additive notation for squareclasses. This is an affine coset,
not necessarily a subset of W. Its local compatibility was not answered
by the old homogeneous descent on W.

At each of the six nonzero parameters, a newly bad prime has an ideal
where the valuation of `η_u` is odd and all twenty anchor valuations
are even. The certificate retains these exact valuations. Thus
`η_u∉W` globally, and `W+〈η_u〉` has dimension 21.
These are six different enlarged spaces: `η_u` varies with u.

The [local protocol](AFFINE_SELMER_PROTOCOL.json) evaluates this one
additional class against the **existing complete local point images**.
It reuses their retained rational x-coordinate witnesses and checks their
local solubility and independence. It uses all the previously certified
bad places, 2 and infinity, without searching for a new witness.

No place outside that support is missing: `N(η_u)=D(u)^4`, and its
ideal valuations are zero outside D. At other odd good places all classes
in the enlarged space are unramified norm-square classes, hence lie in
the local point image. The earlier support argument already covers W.

## Exact affine intersections

| u | dim W_u | Is eta+W locally admissible? | One anchor correction mask | dim Selmer intersection with W+eta |
|---|---:|---|---:|---:|
| -3 | 17 | No | — | 17 |
| -2 | 13 | No | — | 13 |
| -1 | 18 | Yes | 0 | 19 |
| 1 | 13 | No | — | 13 |
| 2 | 13 | Yes | 591872 | 14 |
| 3 | 15 | Yes | 659456 | 16 |

Masks use zero-based bits in the original twenty-point anchor basis.
For each nonempty row, the complete set of admissible anchor masks is
the displayed particular mask plus W_u. A mask gives a cohomology
class, not a point.

Every individual local affine system is soluble, even on the three
excluded rows. What fails is the existence of **one common global
anchor mask** satisfying them all. The portable certificate gives an
explicit binary sum of equations equal to `0=1`. The participating
places are:

| u | Places in the retained contradiction |
|---|---|
| -3 | 13, 6807347, 5070588247, 2517352081717261588189318068375793 |
| -2 | 61, 463, 1033, 88377350996665273, 830763901348107124747241 |
| 1 | 23, 31, 103, 735074448304715650993556022553346486326391 |

This is an exact compatibility obstruction inside the declared global
class space, not evidence that any single completion has no points on E.
The contradiction is stronger than a failed fixed-basis square test:
on these three fibres **no anchor point whose Kummer class lies in W
can satisfy the nonzero transport square condition**, regardless of
which representative of that class is chosen.

## Four CT entries determine the enlarged ranks

Let B be the old alternating CT matrix on W_u, let R be its radical,
and let `ζ∈η_u+W` be the retained locally admissible class.
The enlarged matrix has the form

\[
 \begin{pmatrix}B&v\\v^t&0\end{pmatrix}.
\]

If `v|_R=0`, then `v∈im(B)`. Changing ζ by an element of W_u
removes the new column, so the enlarged rank equals `rank(B)`.
If `v|_R≠0`, the enlarged rank is `rank(B)+2`.
Hence only pairings with a basis of the old radical are needed to
determine the rank; the other new entries are not needed for this
endpoint.

The [CT protocol](AFFINE_CT_PROTOCOL.json) freezes all four required
entries. Each worker constructs and verifies three norm-cover quartics
using the existing backend, then computes Fisher's pairing with exact
local witnesses. Its fact store is isolated under the rank-jump local
artifact directory. Each worker is capped at 60 seconds; all four
completed.

| u | Old radical mask | CT(zeta, mask) |
|---|---:|---:|
| -1 | 450876 | 0 |
| -1 | 596921 | 0 |
| 2 | 513585 | 0 |
| 3 | 602844 | 0 |

The new classes are represented in the cover construction by
`κ_u times the indicated anchor correction`, where

\[
 \kappa_u=1+u\theta+u^2(A+\theta^2),\qquad
 D(u)(1-u\theta)=\kappa_u(1-u\theta)^2.
\]

This reduces coefficient size without changing the class. Each quartic
has a verified square root identifying its cubic invariant with that
class. The pairing implementation and formula are inherited from
[Fisher's binary-quartic method](https://antsmath.org/ANTSXV/papers/ANTS-XV_fisher.pdf);
the four input classes and computations are new here.

| u | Enlarged Selmer-subspace dimension | Exact enlarged CT rank | Restricted radical dimension |
|---|---:|---:|---:|
| -1 | 19 | 16 | 3 |
| 2 | 14 | 12 | 2 |
| 3 | 16 | 14 | 2 |

In each case the affine coset contains a vector in this enlarged
restricted radical. The two-dimensional radicals at `u=2,3` are not
forced by odd-dimensional alternating parity: their spaces have even
dimension. They still do not prove rational points or survival against
Selmer classes outside the computed subspace.

The exact affine radical masks at `u=2,3` have not yet been identified.
Doing that requires the remaining new CT column entries. The table
certifies dimensions and existence within the restricted radical,
without pretending that those entries were computed.

## The extra rational class at u=-1

The repository already proves

\[
 F_u(-Au-1/u)=-D(u)^2/u^3.
\]

Thus `u=-v²` supplies a rational point, as recorded in the
[reassessment](../notes/RANK_JUMP_REASSESSMENT_2026-09-05.md).
At `u=-1` it is

\[
 P_\eta=(A+1,\ 1+A-B).
\]

The new information is its exact class and independence from W.
With `α=θ-θ²` and `γ=1+θ`,

\[
 (x(P_\eta)-\alpha)\gamma^2=D(-1)\gamma=\eta_{-1}.
\]

The new-prime odd valuation proves `η_-1∉W`, so this is a
nonzero rational Kummer direction outside the inherited space.
Both zero CT entries at `u=-1` also follow theoretically from that
rational point; their explicit arithmetic replay supplies a regression.
No new high-rank curve or new point search is involved.

This sharpens the earlier transport cap. Subtracting `P_eta` from a
transported point removes the common shift. Therefore its anchor class
itself must be in the old CT radical, of dimension two. At `u=-1`
at most **two**, rather than the earlier bound of three, independent
anchor classes can be transported simultaneously.

Across the six nonzero controls the resulting caps on independent
anchor-class transport are

    u:       -3  -2  -1   1   2   3
    cap:      0   0   2   0   2   2.

These are not curve-rank upper bounds. Arbitrarily many points can have
dependent Kummer classes, and other global classes lie outside the
declared construction.

## Reproducibility and research decision

The exact records are:

- [Local inputs](../../artifacts/generated-results/elliptic-curves/rank_jump_affine_selmer_inputs_v1.json)
  and [affine systems](../../artifacts/generated-results/elliptic-curves/rank_jump_affine_selmer_v1.json).
- [Four cover/CT witnesses](../../artifacts/generated-results/elliptic-curves/rank_jump_affine_ct_v1.json).
- [Rank and transport consequences](../../artifacts/generated-results/elliptic-curves/rank_jump_affine_selmer_analysis_v1.json).

    python3 elliptic-curves/rank-jump/affine_selmer.py check
    sage -python elliptic-curves/rank-jump/affine_selmer.py verify
    sage -python elliptic-curves/rank-jump/affine_ct.py verify
    python3 elliptic-curves/rank-jump/affine_selmer_analysis.py check

The local replay uses Sage 10.9 / PARI 2.17.3, recorded in its inputs.
It checks exact real signs, reuses complete local image witnesses and
independently verifies admissible products with local-power tests.
The CT replay checks cover maps, invariant square roots, complete
Hilbert support, local witnesses and the final products. It performs
no elliptic-point or conic-witness search. Missing setup arithmetic
can be reconstructed in the isolated rank-jump fact store.

The conclusions are ranked as follows:

1. **Incidence:** the uniform affine shift genuinely enlarges three
   inherited Selmer subspaces. This is more information than the old
   homogeneous dimensions.
2. **Simultaneous compatibility:** local admissibility at every place
   separately can fail to admit one global anchor correction. The
   exact contradictions exclude whole transporter cosets.
3. **Solubility:** the new class at `u=-1` is rational; the added
   CT-compatible directions at `u=2,3` remain unresolved.
4. **Weak explanation:** counting locally admissible or restricted
   radical classes as points. Neither count predicts the large
   Mordell–Weil jumps on its own.

The next bounded step is to finish the new CT column at **u=2 only**,
against a twelve-dimensional complement to its old radical. Together
with the entry already computed, this identifies an exact affine
radical representative and its two-cover equation. Only then should
one test rational solubility or a higher-descent obstruction for that
specific cover. No parameter sweep, broad point campaign or modification
of Agent 1's selector is warranted by the present results.
