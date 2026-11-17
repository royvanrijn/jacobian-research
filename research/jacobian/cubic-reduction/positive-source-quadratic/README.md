# KCUB44 positive-axis source-quadratic field obstruction

The [canonical proof](../POSITIVE_SOURCE_QUADRATIC_FIELD_OBSTRUCTION.md)
proves,for every nonzero q(b,c) of degree two,that the target field of
H composed E_q intersects source degree<=8 in span(1,H1 composed E_q).
Thus no uncapped rational target repair can be Keller through degree eight.
The zero quadratic is retained separately. Endpoint degree is bounded[9,10],
with no degree-nine existence claim or global improvement.

- [input.json](input.json), [checkpoint.json](checkpoint.json): complete
  three-parameter source shear/inverse/composition,augmented weights and
  five formal Taylor obstruction heads with distinct a exponents.
- [base-field9/input.json](base-field9/input.json),
  [base-field9/checkpoint.json](base-field9/checkpoint.json): all220 original
  source columns,213 independent rows,107 good exact target points.
  This original-field degree9 gap is not a transformed-field degree9 theorem.
- [discovery.json](discovery.json): seven fixed-coefficient rank pilots,
  explicitly unused as universal proof.
- [validation.json](validation.json), [integration-validation.json](integration-validation.json):
  independent full model,rank minor,obstruction heads,zero branches,
  degree10 control and current metadata/provenance checks.

Reproduce from the parent with `python raw13_target_field9.py`,
`python positive_source_quadratic.py` and `python verify_positive_source_quadratic.py`.
Universal coefficient-descent,formal induction and weighted initial-space
arguments are written,not formal or externally reviewed. No polynomial
ideal solver is used. All global positive bounds and profile counts remain.
