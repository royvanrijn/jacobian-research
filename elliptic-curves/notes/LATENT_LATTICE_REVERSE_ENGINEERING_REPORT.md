# Finite-aware latent Mordell--Weil lattice recovery

## Abstract

This experiment asks whether a generic Mordell--Weil subgroup can be recovered
from several high-rank specializations without knowing corresponding section
labels, the generic rank, or the originating elliptic-surface family.  The
result is a calibrated computational technique with strong positive and
negative controls, followed by two independently frozen target failures.  No
unknown generic lattice was recovered.  The project is parked pending new
record fibres.

The useful methodological observation is that short-vector height shells are
not sufficiently discriminating, but exact relation replay can be.  Candidate
subspaces are proposed from unlabeled additive-relation components, pruned by
finite `F_2`/`F_3` ranks and height-angle compatibility, saturated exactly,
and ranked by held-out rays replayed per added rank.  The numerical data choose
where to look; exact integer linear algebra and elliptic-curve arithmetic say
what was actually recovered.

## Calibrated method

For a displayed independent subgroup of rank `r`, the method performs the
following finite computation.

1. Enumerate primitive unoriented short integer combinations using the
   canonical-height Gram matrix.  Use the full displayed subgroup and an
   adaptive height bound.
2. Attach exact arithmetic data where available: elliptic-curve sums,
   integrality and coordinate complexity, additive relations, saturation
   indices, and reduction quotient codes.
3. Forget the submitted labels and extract maximal or distinctive relation
   components: stars, overlapping stars, dense `2`-cores, biconnected
   components, and sampled dense hyperplanes.
4. Merge components only after finite-rank and height-angle pruning.  After
   every merge, take the exact primitive rational closure and replay the whole
   relation cloud in the proper subspace.
5. Prefer exact held-out rays replayed per added rank to local relation score.
   Scan dimensions rather than fixing rank 17.

The reusable implementation is in
[`latent_lattice/`](../latent_lattice/).  In particular,
[`components.py`](../latent_lattice/components.py) implements component
extraction, exact rational-space deduplication, proper-subspace replay, and
joint component ledgers; [`codes.py`](../latent_lattice/codes.py) implements
the finite quotient signatures.

## Controls

The calibrated pipeline has four results worth retaining.

- Exact graph-walk consensus recovers the published primitive R17
  specialization in each of the four rank-25 through rank-28 control fibres.
- Cross-dimension persistence selects Fermigier's rank-12 primitive closure at
  both independent height bounds on ICARM 245.  It does not reproduce the old
  false rank-17 core.
- In the rank-25 R17 control, an unlabeled relation-component proposal finds a
  primitive rank-16 component without a supplied center and completes it
  exactly to the published primitive rank-17 subgroup.
- The later strict joint component audit has weaker selection performance than
  these existence/one-fibre results: it identifies the particular source
  component on two of four R17 fibres and completes three of four.  This is an
  important limit on any claim that the component selector itself is already a
  universal four-fibre recovery algorithm.

The detailed artifacts and exact/heuristic boundary are catalogued in
[`LATENT_LATTICE_CALIBRATION.md`](LATENT_LATTICE_CALIBRATION.md).

## Frozen target experiments

No parameter was changed after target inspection under either tag.

### ICARM 351, 356, 376, 377, 385

The first tag was `LATENT-LATTICE-WGXLI-FROZEN-2026-09-01-v1`.  The selected
dimensions were `10, FAIL, 10, 13, FAIL`.  No dimension recurred in four
fibres, so the frozen gate prohibited component matching and every later
stage.  The result is recorded in
[`latent_lattice_wgxli_frozen_dimension_v1.json`](../../artifacts/generated-results/elliptic-curves/latent_lattice_wgxli_frozen_dimension_v1.json).

### E29 and ICARM 398--400

The second tag was
`LATENT-LATTICE-E29-398-400-FROZEN-2026-09-01-v1`.  It inherited every
computational parameter from the first freeze and changed only the four target
identities and the corresponding four-of-four hold-out wording.

| fibre | ambient rank | frozen cloud | dimension result |
|---|---:|---:|---:|
| E29 / ICARM 12 | 29 | bound 62, 2,294 rays | fail: dimensions 19 and 20 absent from at least one cutoff |
| ICARM 398 | 30 | bound 20 did not finish in 600 seconds | fail closed: frozen resource bound |
| ICARM 399 | 29 | bound 52, 1,978 rays | 12 |
| ICARM 400 | 28 | bound 46, 1,958 rays | 16 |

Curve 398's failure was isolated only so that it could not suppress the 399
and 400 audits; its timeout was not enlarged.  Thus there was no common
dimension and no authority to run relation-component, primitive-Hermite,
finite-index, abstract-height-lattice, K3, or equation stages.  The consolidated
artifact is
[`latent_lattice_secondary_frozen_dimension_v1.json`](../../artifacts/generated-results/elliptic-curves/latent_lattice_secondary_frozen_dimension_v1.json),
with SHA-256
`e89f6156a22f7a5b762e6c525db7b95b4d625ac1a96d712ca8f4053aa05f62cb`.

## Mathematical status

What is exact within the recorded finite computations includes primitive
integer coordinates, additive relations, rational-span saturation, Smith
factors, elliptic-curve group-law replay, finite quotient codes, and the public
input hashes.  Canonical heights, likelihood scores, component proposals, beam
coverage, and inferred dimension are numerical or statistical.

Neither target failure is a nonexistence theorem for a common primitive
sublattice of dimensions 10 through 20.  In particular, the values 10, 12, 13,
and 16 are selector outputs, not generic-rank theorems.  The experiments prove
only that the two precommitted bounded protocols did not pass their own gates.
They therefore do not meet the original success level C, which required an
exhaustive bounded nonexistence result.

The publishable result is presently methodological: controls demonstrate that
finite-aware unlabeled relation replay can distinguish a genuine rank-12
generic closure from a spurious rank-17 short-vector core and can recover a
non-oracle rank-16 component of R17 before exact completion.  Evidence for an
unknown family would require a new, independently frozen test set.  Until more
record fibres become available, latent-lattice reverse engineering is retired
rather than retuned on either failed target cluster.

