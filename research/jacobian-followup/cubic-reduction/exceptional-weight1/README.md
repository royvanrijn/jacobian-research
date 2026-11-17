# Exceptional positive-unit weighted family

The accepted all-r>=4 result is [KCUB34](../EXCEPTIONAL_UNIT_WEIGHT_CLASSIFICATION.md).
Run [the independent checker](../verify_exceptional_weight1.py). Its
`validation.json` is the exact rational replay receipt.

`checkpoint.json` pins the four complete primary models r4,r5,r6,r7.
Each `r*.json` retains its full source basis, source and quotient outputs,
all native coefficient equations and universal determinant. The primary
`classification_status` remains a model-only status; the separate accepted
classification is proved by KCUB34, not by a solver answer in that model.
`structural-identities.json` retains the universal univariate/graph
identities, normalized r5 residual and its five-term unit, and both
constant-first quadratic-tail exclusion identities.

The preliminary nine-equation r5 Groebner basis[1] is primary discovery
only and is unused as a proof premise. The retained E0,E1 forms reconstruct
all nine equations. The accepted proof instead uses three denominator-free
linear eliminations and a directly multiplied five-term univariate unit.

The first independent replay stopped before any algebra at an export
format mismatch: the producers initially wrote [exponent,coefficient]
pairs, while the common decoder expects named fields. Both producers were
changed to the repository format and their full models regenerated; the
independent checker then passed. No failed run or primary verdict is
interpreted as a mathematical exclusion. Earlier r2/r3 UNKNOWN/timeouts
are outside this new calculation and remain unchanged.

The all-r modules, rational derivative and root arguments and complete
normalizations are written, not formalized or independently externally
reviewed. No external plane theorem is imported. There is no improved
positive witness or global simplicity claim.
