# Strong quantum residue under polynomial Weyl shears

The [QFL1 unrestricted formal lifting theorem](../dixmier-formal-lifting/PROOF.md)
now gives a polynomial right inverse for every fiber defect. Explicit
corrections of Z orders 9 and 8 cancel the retained order-five defect.
The remaining general question is finite or arithmetic realization;
QRT1--QRT3 below remain restricted-category obstructions.

The [QRT3 uniform positional theorem](UNIFORM_POSITIONAL.md) now covers
positional polynomial Hamiltonians of arbitrary degree. The quartic test
passes independently; boundary-weight decomposition then proves period
one on an unbounded coefficient category with fixed leading-Z bounds.
The generic rank minors and every exceptional integer weight have an
independent exact replay. A first mixed-polarization test shows why a
native extension must include the transverse momentum and all six relations.

The next gate is now closed by [QRT2](COMBINED_COCYCLE.md): a new 16-term
strong cocycle annihilates the combined native and transported correction
spaces and has period one on every compatible lower lift. It allows all
finite X-Laurent poles and covers the whole positional cubic shear class.
The finite certificate has a complete independent Fraction replay.
[Combined-space validation](combined_validation.json) records the exact
coverage, source checks and completion audit.

The [proof](PROOF.md) transports the retained order-five strong residue
through arbitrary positional polynomial Weyl shears, including the
nonlinear example p_x -> p_x+lambda X^2. Exact quantum transport preserves
its nonzero value -47547660815739/190658 on the transported correction
category.

The support comparison finds a sharp limitation: under that example, the
native allowed correction X^12 Q^6 Z pairs to -1296 lambda. Keeping the
old numerical degree bounds is therefore not a valid transport of this
cocycle. Adding the corresponding single new direction does not solve
all order-five equations at the selected lower lift.

The [symbolic replay](verify_transport.py) checks all 614 current columns,
all 945 lower-family coefficients and 35 gauge controls. A
[separate scalar replay](check_independent.py) verifies the period,
support failure and six ambient shear relations. The all-degree assertion
uses the written automorphism proof, with the source finite-Laurent
strength inherited from the retained residue theorem.

```sh
.venv/bin/python research/jacobian-followup/dixmier-residue-transport/verify_transport.py
.venv/bin/python research/jacobian-followup/dixmier-residue-transport/check_independent.py
```

[Inputs](inputs.json) · [receipt](receipt.json) ·
[independent coverage](independent.json) ·
[mathematical authority](../../MATH_STATUS.json).

QRT1 is the original orbit-covariant restricted result. QRT2 supplies the
combined-space extension, and QRT3 removes the positional Hamiltonian degree
bound. QFL1 settles unrestricted formal fiber lifting. Finite realization,
compatible fourth images and prescribed arithmetic lifts remain separate.
No DC2 solution or novelty priority is claimed.

[Uniform validation](uniform_validation.json) records the quartic and
all-degree certificate coverage and the pre-existing navigation blockers.

[Scoped validation](validation.json) checks the new sources, hashes and
reviewed metadata. Repository-wide rendering remains blocked by the
pre-existing missing KDSQ6 consumer and METHOD-EC-CENTRE-POLICY artifact;
generated navigation was not rewritten around those gaps.

<!-- status-consumer: QRT1 eadd3731ecd8c9bd -->

<!-- status-consumer: QRT2 c4022a1fc6968780 -->

<!-- status-consumer: QRT3 d003fd8d2ba9bcc3 -->

<!-- status-consumer: QFL1 9163a4d2d52475ed -->
