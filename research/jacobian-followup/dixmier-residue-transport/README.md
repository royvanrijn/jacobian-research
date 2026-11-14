# Strong quantum residue under polynomial Weyl shears

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

QRT1 is an orbit-covariant restricted result. A cocycle surviving the sum
of transported and native correction spaces remains open. No unrestricted
quantization obstruction, DC2 solution or novelty priority is claimed.

[Scoped validation](validation.json) checks the new sources, hashes and
reviewed metadata. Repository-wide rendering remains blocked by the
pre-existing missing KDSQ6 consumer and METHOD-EC-CENTRE-POLICY artifact;
generated navigation was not rewritten around those gaps.

<!-- status-consumer: QRT1 eadd3731ecd8c9bd -->
