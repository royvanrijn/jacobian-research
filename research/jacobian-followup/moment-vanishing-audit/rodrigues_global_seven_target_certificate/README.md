# Global Rodrigues gate: seven exact membership tests

Read `PROOF.md` for the complete statement, qualifications and proof. The global length question remains open; this packet gives an equivalent, smaller certificate to attempt.

## Replay

Requirements: Python 3 and SymPy.

Run from this directory:

```sh
python verify_gate.py
```

This runs the complete invariant-coordinate adapter once, then checks the adapted target powers, the exact cubic relation, and the four-parameter support certificate. Typical runtime in the producing environment was under one minute. Every arithmetic step uses exact rational or integer arithmetic.

The original seventeen-dimensional local certificate is included unchanged in `rodrigues_local_length17_certificate/`. Its retained `certificate_data.json` is checked against its manifest before use. The new adapter imports the included Taylor scripts and uses the actual local multiplication matrices.

New output JSON files are regenerated deterministically. The runner checks that they match the included files byte-for-byte. Logs contain timing information and are not compared byte-for-byte.

## Contents

- `PROOF.md`: both equivalent seven-condition global gates, the augmented global ideal identity, and the twelve-component grading.
- `global_gate_targets.json`: machine-readable target definitions, powers, cubic relation coefficients, degrees and weights modulo twelve.
- `GLOBAL_TARGET_NILPOTENCE_ADAPTER.md`: details of normalization and the individual target exponent profile.
- `global_target_nilindices.py` and `.json`: full cyclic adapter and actual normalized local vectors.
- `global_adapted_targets.py` and `.json`: exact powers of the adapted target functions.
- `global_cubic_relation.py` and `.json`: the exact cubic replacement, including its primitive integer coefficient vector.
- `odd_slice_support_certificate/`: complete four-parameter reduction, exact moment polynomials, rational Bezout identities, Jacobian, and independent replay.
- `rodrigues_local_length17_certificate/`: the retained local proof and data.

## What is established

The included calculations and arguments prove the first-six support result on the selected fibre, its first-nine reduced scheme closure, all local powers used in the gate, and the exact local cubic relation. They do not produce a membership witness for any of the seven requested polynomials in the unrestricted global ideal `I14`.

Use the full `I14`, including moments thirteen and fourteen. The first-six and first-nine statements concern the locus obtained after adjoining the seven selected coordinate equations.
