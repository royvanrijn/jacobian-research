# KCUB43 complete affine degree-eight shear family

The [canonical proof](../QUADRATIC_TARGET_AFFINE_DEGREE8.md) gives a uniform
lower16 for every affine source/output frame of every nonzero shear
H+(0,L,M)*H1^2. The existing19-term control gives upper19;the sharp nonlinear
minimum within[16,19] is not established. No member improves raw13.

- [input.json](input.json): complete translated homogeneous matrices,
  degrees2..8,483 shear entries,all base entries and square rows;all three
  centre coordinates and both shear parameters remain formal.
- [rigidity.json](rigidity.json): coupled top-degree and cubic equality
  obstructions with a free third source row.
- [checkpoint.json](checkpoint.json): completed producer provenance.
- [validation.json](validation.json): independent full translation,all rank
  separators/zero branches/equality coefficients and matching19-term control.
- [integration-validation.json](integration-validation.json): current ledger,
  parent receipts,source pins,reviews and navigation boundary.

Reproduce from the parent with `python quadratic_target_affine8.py` and
`python verify_quadratic_target_affine8.py`. No solver or grid is used.
Universal rank/UFD/squarefreeness proofs are written,not formalized or
externally reviewed. KCUB41/42 imports retain their stated boundaries.
Global positive bounds and unrestricted whole-profile counts are unchanged.
