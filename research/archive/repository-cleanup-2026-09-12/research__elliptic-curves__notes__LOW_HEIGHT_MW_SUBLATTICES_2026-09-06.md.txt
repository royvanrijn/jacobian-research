# Low-height MW sublattices and common-cover reconstruction

Two bounded combination searches on ICARM 302 and the masked ICARM 245
control are complete. **Neither new proposal method passes the required
Fermigier reconstruction calibration.** The candidates are explicit integer
sublattices of the certified displayed groups, not generic-family identifications.
The authoritative result is `EC-LOW-HEIGHT-MW-SUBLATTICES-20260906` in
[`MATH_STATUS.json`](../../MATH_STATUS.json).

The [summary certificate](../../artifacts/generated-results/elliptic-curves/low_height_mw_sublattices_summary_v1.json)
binds the selections, calibration results and arithmetic replay. No rank bound
changes. This is a user-authorized reopening for curves 245 and 302, not a
resumption of the retired wgxli or E29/398–400 target campaigns.

## What was searched

Both curves use every displayed point, respectively 20 and 31 independent
generators. PARI computes canonical-height matrices at 80 decimal digits and
unimodular ambient LLL bases. The short-ball bounds are 28 and 70, retaining
1,928 and 3,458 primitive unoriented vectors. These are numerical height-ball
enumerations with a checked 10,000-line cap, not interval-certified boundaries.
Independent enumeration at 100 digits returns exactly the same coordinate sets.

The first pass adds every signed combination of up to three ambient LLL
basis vectors, with coefficients ±1 and one orientation. The sparse shells
have 4,960 and 18,941 lines; their unions with the balls have 6,798 and 22,278.
These are combinations in all public generators; neither the ball nor the
sparse shell is a subset of the printed point list.

Forty-eight growing bases use 16 seeds each for conditional determinant,
unit-relation degree and deterministic random relation-weighted growth. Every
path visits ranks 8–20, giving 624 proposals per curve. Exact saturation and
Hermite normal forms leave 409 and 420 distinct primitive spaces. Two finalists
per rank and score (determinant, ball density and parity collision count), with
duplicates removed, retain 49 and 48 candidates. Saturation here is inside the
displayed free group, not a claim about saturation at every prime in all of
`E(Q)`.

The second pass was frozen separately after the first calibration failure.
For **every** short-ball parity representative `c`, it examines

```text
c + 2z, where z = 0 or a signed support-at-most-two ambient LLL vector.
```

This gives **1,544,328** representative examinations on 245 and **6,649,734**
on 302, including repeated rays. In each coset, the first 128 by height are
available for greedy independent growth through ranks 8–20. There are 25,064
and 44,954 rank-specific proposals. Eight per rank and score (generated
determinant, maximum generator height, height ratio) enter exact scoring:
282/279 shortlisted presentations and 223/273 distinct primitive spaces.
The final selections retain 46/41 candidates. This finite shifted box is not
a complete closest-vector search in a coset.

The same selection code and budgets apply to both curves. Each selection is
written before its separate evaluator reads the generic embedding. These are
retrospective masked controls: the family and previous calibration history
were already known. They are not fresh, previously unseen statistical tests.

## Curve-302 candidates

The table shows the smallest **primitive-closure** determinant found by the
first pass at each rank, the number of its rays in the height-70 ball and the
heights of its reduced basis. Determinants are compared within a rank; these
values are not global minima or a significance test.

| Rank | Determinant | Ball lines | Reduced basis heights |
|---:|---:|---:|---:|
| 8 | 2.45016e12 | 44 | 51.329–65.059 |
| 9 | 5.64524e13 | 64 | 51.329–66.605 |
| 10 | 1.40303e15 | 84 | 52.274–68.821 |
| 11 | 3.18072e16 | 112 | 52.274–68.821 |
| 12 | 5.31221e17 | 152 | 51.329–69.452 |
| 13 | 8.64417e18 | 224 | 51.876–64.657 |
| 14 | 1.27016e20 | 321 | 51.329–67.693 |
| 15 | 1.89856e21 | 454 | 51.329–63.870 |
| 16 | 3.15187e22 | 577 | 53.248–66.555 |
| 17 | 3.62282e23 | 864 | 55.112–64.064 |
| 18 | 5.96955e24 | 991 | 52.274–64.064 |
| 19 | 9.60924e25 | 1123 | 54.195–65.028 |
| 20 | 1.69183e27 | 1225 | 53.510–64.331 |

The selected rank-17 lattice is **exactly equal**, by integer HNF, to the
previous [rank-17 core](RECORD_CURVES_28_29_273_302_HEIGHT_LATTICES.md).
Its minimum is 51.32851052925419, below its displayed reduced-basis diagonals.
This is a successful rediscovery of an existing structural candidate, not
new evidence of generic rank 17.

Every finalist includes its generated basis, primitive closure, reduced
basis and exact index in public-point coordinates:

- [First-pass curve-302 selection](../../artifacts/generated-results/elliptic-curves/low_height_mw_sublattices_v1_302_selection.json).
- [Common-cover curve-302 selection](../../artifacts/generated-results/elliptic-curves/common_cover_mw_sublattices_v1_302_selection.json).
- [Actual generated-subgroup reductions](../../artifacts/generated-results/elliptic-curves/common_cover_generated_lattices_v1.json), retaining their larger determinants and reduced Grams separately from the primitive closures.

## Parity, covers and complements

Fresh finite-reduction certificates have mod-2 column ranks 20 and 31. Thus
coefficient parity exactly distinguishes the displayed subgroups' images in
`E(Q)/2E(Q)`. Equal parity proves that the associated points lift to the same
2-cover class; different parity is separated by the finite certificates.
No explicit quartic equation for each cover is constructed here.

There are no repeated parity classes in either short ball. For distinct
unoriented vectors `v,w` with equal parity, both `(v+w)/2` and `(v-w)/2` are
nonzero lattice vectors. Their heights sum to `(h(v)+h(w))/2`. Consequently
a ball below twice the ambient minimum cannot contain such a pair.

The extended sparse pools do have collisions, but their largest multiplicity
four and collision counts 7,030/27,435 are entirely accounted for by sign
choices on supports two and three. They are not exceptional cover evidence.
The second pass constructs up to 20 independent points in one cover by design.
Any rank-r subgroup generated in one nonzero parity class has index divisible
by `2^(r-1)` in its primitive closure. In particular, merely seeing index
`2^11` at rank 12 is not recognition of Fermigier's actual subgroup.

For each proper-rank finalist, an exact unimodular completion supplies the
**projected quotient** Gram through a Schur complement. It is not an exact
orthogonal integral complement. LLL, scalar normalization and rounding give
only a coarse diagnostic; the table's candidates have rounding RMS roughly
0.28–0.31. No exact root-lattice complement or generic integral Gram is
identified. Those fields remain `UNKNOWN`.

## Calibration and verification

The first pass's best rank-12 intersection with the hidden generic space is
**8/12**. The common-cover pass improves this to **10/12**, evaluated over
all 1,928 rank-12 coset proposals before preselection. Neither ledger contains
the exact rank-12 space, so this is a proposal-recall failure, not just a poor
final ordering. A containing rank-18 or rank-20 space does not pass recovery.

As a separate positive regression, the existing frozen 128-proposal
graph-consensus ledger selects source index 65 and exactly recovers the
primitive Fermigier rank-12 closure. This repeats the earlier successful
control; it does not validate the two new generators. The actual generic
subgroup has index 2,048 in that closure. Neither new pass recovers it.

All **184 finalists** pass exact rank, primitive index, basis-lattice equality,
span support and parity replay. Fresh 100-digit height matrices differ from
the 80-digit inputs by less than `5e-79`. Independent minimum searches agree,
and 24 distinct minimum witnesses are materialized by exact elliptic group
law and checked against directly computed canonical heights. Four focused
regressions pass. The separate common-cover audit verifies actual generated
subgroup reductions and the determinant identity `det(A)=index^2 det(S)`.

The [preverification archive](../../archive/elliptic-curves/low-height-mw-sublattices-preverification-2026-09-06/README.md)
preserves implementation corrections. In particular, bounded `qfminim`'s
second field is not generally the minimum; the active script explicitly
minimizes the returned norms. This correction changes reported minima and
witnesses, not proposal generation or selection. Failed controls were not
retuned after this correction.

## Replay

From the repository root, set `OPENBLAS_NUM_THREADS=1`. Both search scripts
have immutable checkpoints; remove or relocate only this experiment's own
checkpoints to regenerate them from scratch. Every GP invocation is bounded
by 120 seconds. No broader search is started by the verifiers.

```sh
.venv/bin/python elliptic-curves/cas/search_low_height_mw_sublattices.py --curve 245
.venv/bin/python elliptic-curves/cas/search_low_height_mw_sublattices.py --curve 302
.venv/bin/python elliptic-curves/cas/search_low_height_mw_sublattices.py --evaluate
.venv/bin/python elliptic-curves/cas/search_common_cover_mw_sublattices.py --curve 245
.venv/bin/python elliptic-curves/cas/search_common_cover_mw_sublattices.py --curve 302
.venv/bin/python elliptic-curves/cas/search_common_cover_mw_sublattices.py --evaluate
.venv/bin/python elliptic-curves/cas/verify_low_height_mw_sublattices.py --check
.venv/bin/python elliptic-curves/cas/audit_common_cover_generated_lattices.py --check
.venv/bin/python elliptic-curves/cas/summarize_low_height_mw_sublattices.py --check
PYTHONPATH=elliptic-curves:elliptic-curves/cas .venv/bin/python -m unittest discover -s elliptic-curves/tests -p test_low_height_mw_sublattices.py
```
