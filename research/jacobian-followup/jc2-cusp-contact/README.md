# F2 cusp contact and adapted-quintic filtration

14 November 2026. This follows the [cusp-rigidity adoption](../jc2-cusp-filtration-rigidity/README.md)
under the explicit ambitious F2 research goal.

- **JC2CFSB1:** on actual contact507 squarefree candidates, the first two
  source bands vanish and the next root correction has one parameter. The
  double-carrier first parameter survives its nonlinear quadratic row.
- **JC2CFGA1:** affine target normalization reduces eight contact values
  to five normal forms:507,520,585,715,910.
- **JC2CF7F1:** the four remaining carrier invariants give a finite flat
  degree-seven map on the normalized rational-quintic coefficient family.
- **JC2CFMV1:** a mixed-volume and mandatory boundary-intersection audit
  bounds the geometric degree by2019, or1347 on the contact507 branch.
- **JC2CFJ1:** the zero-seven-centre target-transport scheme is exactly the
  adapted-quintic Euler scheme. This is a conditional attachment test;
  generic carrier centres do not force the global attachment.
- **JC2CFN1:** combine the two retained Newton faces to classify the complete
  highest contact spectrum:507,520,585,650,715,780,845,910. The first
  branch has a forced nonzero coefficient and Wronskian degree51.
- **JC2CFW1:** the whole-source Keller restriction has derivative ideal (z).
  A Wronskian identity improves the F2 contact bound from456 to457,
  excluding the contact-degree-456 stratum. Centering at the actual
  source-origin image gives bound459, with exact tangent branches461/462.
- **JC2CFQ1:** classify every preserving polynomial derivation of the
  complete weight-(3,5,7) deformed-quintic filtration. On the entire
  retained four-parameter rational quintic family, its eligibility ideal
  is a length-seven scheme supported only at the monomial E8 endpoint.
  Every other row has no polynomial Euler field and no polynomial
  filtered conjugacy to any binomial-cusp filtration.
- **JC2CFC1:** classify the first five high-contact coefficient strata
  on the commuting control slice. These Jacobian-zero controls explain
  why early resonant freedoms cannot establish the desired upper bound.

[Actual source coefficient obstruction](SOURCE_BANDS.md) ·
[Affine gauge and degree-seven interface](CARRIER_QUOTIENT.md) ·
[Contact507 degree atlas and jet compatibility](BRANCH507.md) ·
[Complete proofs and exact remaining gates](PROOF.md) ·
[Exact checker](verify.py) · [Replay receipt](receipt.json).

Run from the repository root with SymPy:

    .venv/bin/python research/jacobian-followup/jc2-cusp-contact/verify.py

Additional exact branch checks:

    .venv/bin/python research/jacobian-followup/jc2-cusp-contact/verify_branch507.py

    .venv/bin/python research/jacobian-followup/jc2-cusp-contact/verify_carrier_quotient.py
    .venv/bin/python research/jacobian-followup/jc2-cusp-contact/verify_source_bands.py

The large historic F2 circuit artifact is absent and was not reconstructed.
The results above use retained source statements and new small exact
calculations. They do not certify a full nonlinear-circuit point, exclude
all F2 candidates, settle the E8 global endpoint, or prove JC2/HC4.
The full nonlinear coefficient gate remains open; the exact leading strata
and the complete adapted-quintic Euler test are established here.

[MATH_STATUS.json](../../MATH_STATUS.json) is authoritative. Universal proofs
are written; finite checks support the stated algebraic scopes. No formal
verification, external review or novelty priority is claimed.
