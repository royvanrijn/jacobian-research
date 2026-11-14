# Symmetry groups of Keller maps

The [theory and proofs](THEORY.md) distinguish a permissible permutation action, a cover, a determinant-one chart, and an actual polynomial self-map of affine space.

The affine obstruction is stronger than exclusion of regular actions: a point stabilizer must normally generate the geometric monodromy group. Consequently no nontrivial nilpotent group, and no natural even-degree dihedral action, occurs. For PSL2(11), finite inertia must fix a sheet. The [explicit F11 audit](F11.md) replaces the modular genus-one presentation by a known genus-zero (5,5,2) rational cover, while exposing an arithmetic constant-field enlargement.

[Exact checker](verify.py) · [coefficient model](model.py) · [receipt](receipt.json) · [primary sources](sources/manifest.json) · [mathematical status](../../MATH_STATUS.json).

Run from the repository root:

    .venv/bin/python research/jacobian-followup/monodromy-theory/verify.py

The universal obstruction has a written proof using purity and characteristic-zero simple connectedness of affine space. The checker verifies finite-group consequences, the F11 polynomial identities, its boundary ledger and arithmetic discriminant. It is not a formal proof or an independent review. The full action spectrum and absolute PSL2(11) realization remain open.
