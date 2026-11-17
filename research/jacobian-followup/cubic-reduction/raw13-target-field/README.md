# KCUB42 full target-field controls

The [canonical proof](../RAW13_TARGET_FUNCTION_FIELD_RIGIDITY.md) proves
uncapped rational-target rigidity through source degree seven and classifies
the complete first nonlinear layer through degree eight. All degree-eight
Keller targets normalize to Y+(0,lambda,mu)*Y1^2. With the source fixed,
a nonzero shear needs19 terms;affine source frames of that family remain open.

- [input.json](input.json), [checkpoint.json](checkpoint.json): complete
  degree<=7 source basis120,rank116,58 exact target specialization points.
- [degree8/input.json](degree8/input.json),
  [degree8/checkpoint.json](degree8/checkpoint.json): complete basis165,
  rank158,79 points. These are exact minor certificates,not candidate grids.
- [known-degree8.json](known-degree8.json): exact Q(H),S(H),T(H) controls.
- [target-locus/input.json](target-locus/input.json): full nine-parameter
  tangent-identity determinant model with27 generators.
- [target-locus/certificate.log](target-locus/certificate.log): seven exact
  linear memberships,17 multiplier terms,maxdegree1;all multiplied independently.
- [discovery.json](discovery.json): real leading cancellations that invalidate
  a naive leading-form shortcut;no SAGBI completeness result is assumed.
- [validation.json](validation.json), [integration-validation.json](integration-validation.json):
  independent cubic/recovery identities,companion/adjugate reconstruction,
  both full modular minors,target model and current repository integration.

Reproduce with `python raw13_target_field.py`, `python raw13_target_field8.py`
and `python verify_raw13_target_field.py` from the parent. The retained
Singular scripts/logs suffice for exact replay of the small locus proof.
Unknown target degrees and denominators are unbounded. Universal arguments
are written,not formal or externally reviewed. Global minima remain open.
