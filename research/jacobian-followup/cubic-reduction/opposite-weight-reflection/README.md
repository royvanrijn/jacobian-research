# Opposite-weight model and replay artifacts

The accepted result is [KCUB33](../OPPOSITE_WEIGHT_CLASSIFICATION.md).
Run [the independent checker](../verify_opposite_weight_reflection.py);
its `validation.json` is the accepted exact-arithmetic receipt.

`checkpoint.json` and the fourteen `m*q*.json` inputs retain complete
source bases, boundary reductions and invariant dictionaries. The
`universal-identities.json` input retains the three symbolic-weight
coefficients, seven support masks and nine univariate leading cases.
These are independently reconstructed, rather than trusted as primary
symbolic answers.

The five `*-unit.json` certificates, their Singular scripts/logs and
`lift-checkpoint.json` are primary-only discovery evidence. They are
explicitly unused in KCUB33: the short diagonal-coefficient proof removes
the need to multiply these large certificates. The earlier `weighted212/`
screen and its producer are also retained unused. No timeout or generic
primary result is promoted to an accepted proof premise.

The initial producer's empty-support Python integer/substitution error
was corrected before its complete checkpoint. The independent replay
uses a separate rational sparse-polynomial implementation. The theorem
still has the written, unformalized field/module/degree trust boundary
specified in the canonical proof. It does not cover m=1 or improve a
positive counterexample bound.
