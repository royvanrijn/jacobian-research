# Eight term Keller collision exclusion

Every normalized three-variable Keller collision with raw coordinate degrees
at most seven requires **at least nine nonlinear monomial occurrences**,
counted with coordinate multiplicity.

All 17,482 eight-term support orbits are excluded over Q in Singular and
SymPy. Independent support enumeration and derivative-level coefficient
compilation verify the complete cover. Another 2,326 smaller orbits are
retained as regressions.

- [Proof and scope](PROOF.md)
- [Coordinate gates and a sharp dimension-four special form](COORDINATE_GATES.md)
- [Complete manifest](results/manifest.json)
- [Exact replay validation](validation.json)
- [Independent coverage](coverage-validation.json)
- [Structural classification](structural-validation.json)
- [Scoped metadata checks and navigation limitation](integration.json)

From the repository root:

```sh
.venv/bin/python research/jacobian-followup/sparse-eight/verify.py --replay-algebra --replay-coverage
.venv/bin/python research/jacobian-followup/sparse-eight/verify_coordinate_gates.py
```

Nine-term attainment remains open. Global degree minima and unrestricted
support or special-form dimension minima are separate questions.

