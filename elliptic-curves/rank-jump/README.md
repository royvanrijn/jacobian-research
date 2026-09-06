# Rank-jump understanding

Retrospective analyses with pinned inputs, independent of the active curve search.

- [Latest: explicit projection fibres and the residual solubility obstruction](EXPLICIT_PROJECTION_FIBRES.md)
- [Ordinary Jacobian local conditions and the CT comparison](JACOBIAN_LOCAL_CONDITIONS_AND_CT.md)
- [The transporter in the torsion-difference and CT extension](TORSION_DIFFERENCE_AND_CT.md)
- [One Jacobian for the Sha block and a soluble splitting event](JACOBIAN_GLUING_AND_SHA_BLOCKS.md)
- [Exact relative full-Selmer dimensions and a common unknown quotient](RELATIVE_FULL_SELMER_THEOREM.md)
- [Full-Selmer ramification block and two unresolved affine covers](FULL_SELMER_RAMIFICATION_BLOCK.md)
- [New affine Selmer classes and their CT obstructions](AFFINE_SELMER_AND_CT.md)
- [An explicit soluble block and a genus-two splitting event](LINEAR_TWIST_SOLUBLE_BLOCKS.md)
- [Shared quadrics and distinct ruling base changes](SHARED_QUADRICS_AND_RULINGS.md)
- [Branch characters and exact specialization collapse](BRANCH_BLOCKS_AND_SPECIALIZATION.md)
- [Two-adic and real completion of the local-support test](DYADIC_REAL_QUOTIENT_SUPPORT.md)
- [Small bad-prime quotient supports and the global kernel](BAD_PRIME_QUOTIENT_SUPPORT.md)
- [Independent halving fields and arithmetic blocks](HALVING_FIELDS_AND_BLOCKS.md)
- [Cubic-field bridge and the norm-solubility gap](CUBIC_BRIDGE_AND_NORM_GAP.md)
- [Simultaneous CT blocks and normalization-invariant interpretation](CT_VARIATION_AND_BLOCKS.md)
- [Exact local collision cuts, reciprocity and 4-division separation](LOCAL_COLLISION_AND_RECIPROCITY.md)
- [Initial findings, paired studies, mechanisms and missing implications](ANALYSIS.md)
- [Per-fibre quotient-block reports](FIBRE_REPORTS.md)
- [Comparison panel (CSV)](../../artifacts/generated-results/elliptic-curves/rank_jump_comparison_panel_v1.csv)
- [Frozen small experiment](EXPERIMENT.json)

Replay from the repository root, using Python 3.12 and its standard library:

```sh
python3 -m unittest discover -s elliptic-curves/rank-jump -p 'test_*.py'
python3 elliptic-curves/rank-jump/retrospective.py check
python3 elliptic-curves/rank-jump/blocks.py check
python3 elliptic-curves/rank-jump/cover_experiment.py check
python3 elliptic-curves/rank-jump/render.py --check
python3 elliptic-curves/rank-jump/local_collision.py check
python3 elliptic-curves/rank-jump/reciprocity.py check
python3 elliptic-curves/rank-jump/four_division.py check
python3 elliptic-curves/rank-jump/ct_variation.py check
python3 elliptic-curves/rank-jump/cubic_bridge.py check
python3 elliptic-curves/rank-jump/halving_fields.py check
python3 elliptic-curves/rank-jump/bad_prime_support.py check
python3 elliptic-curves/rank-jump/dyadic_real_support.py check
python3 elliptic-curves/rank-jump/branch_blocks.py check
python3 elliptic-curves/rank-jump/quadric_rulings.py check
python3 elliptic-curves/rank-jump/linear_twist_blocks.py check
python3 elliptic-curves/rank-jump/affine_selmer.py check
python3 elliptic-curves/rank-jump/affine_selmer_analysis.py check
python3 elliptic-curves/rank-jump/ramification_block.py check
python3 elliptic-curves/rank-jump/u2_affine_radical.py check
python3 elliptic-curves/rank-jump/selmer_comparison.py check
python3 elliptic-curves/rank-jump/jacobian_sha_blocks.py check
python3 elliptic-curves/rank-jump/torsion_difference.py check
python3 elliptic-curves/rank-jump/jacobian_local_conditions.py check
```

Optional independent Hilbert-symbol check: `sage -python -m unittest discover -s elliptic-curves/rank-jump -p test_local_collision.py`.

Explicit Jacobian fibre and CT replay: `sage -python elliptic-curves/rank-jump/projection_pairings.py verify`.

The portable `rank_jump_*inputs_v1.json` files retain the exact projections
needed for replay. `capture` and `build` refuse to overwrite outputs; they
are construction commands, not replay requirements. Capture additionally
needs the pinned git history and hash-matching historical transcripts or
bundles. No command runs a point search, parameter sweep or full descent.
`MATH_STATUS.json`, generated `STATUS.md`, search policies and live outputs
are unchanged.
