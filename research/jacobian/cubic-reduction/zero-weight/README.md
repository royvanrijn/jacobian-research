# KCUB39 exact controls

The [canonical proof](../ZERO_WEIGHT_CLASSIFICATION.md) classifies every
three-variable Keller map with exactly one zero diagonal weight in all
degrees and proves sharp collision dimension four for that form. The
bounded files here are controls, not the proof of the unbounded assertion.

- `p1q1.json`, `p1q2.json`, `p1q3.json`, `p1q4.json`, `p1q5.json`, `p1q6.json`,
  `p2q3.json`, `p2q5.json`, `p3q4.json`: complete degree-seven native modules,
  determinants, coefficient ideals and principal faces.
- [checkpoint.json](checkpoint.json): completed producer provenance;
  240 atoms, 438 coefficient generators, 162 removed pure-x atoms,
  78 remaining atoms and twelve positive principal faces.
- [validation.json](validation.json): independent Fraction reconstruction,
  all faces, 27 rank-one and 16 tail leading controls, nine generic shear
  inverses, same-sign block identities and the full stabilized old upper map.
- [integration-validation.json](integration-validation.json): source/import
  pins, current metadata, work reviews and navigation-check boundary.

Reproduce with `python zero_weight.py`, then `python verify_zero_weight.py`
from the parent directory. No solver query, coefficient grid or timeout
increase is involved. Old111 timeouts and unused primary evidence are
unchanged. The full111 raw<=12 frontier remains open; global positive bounds
and whole-profile counts do not change.
