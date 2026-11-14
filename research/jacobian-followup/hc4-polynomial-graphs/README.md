# All-degree Meng--Yang polynomial-graph exclusion

For the scaled Meng--Yang potential and every polynomial graph
r=R(x,y,p,q), a nonzero constant Hessian determinant forces the pulled-back
gradient to be a polynomial automorphism. Consequently **no graph in this
standard chart can both have nonzero constant Hessian determinant and
retain a gradient collision, in any degree**.

The [proof](PROOF.md) combines two results:

- **HC4MYGAP1:** an elementary fixed-slice and polynomial-degree argument
  excludes every constant two-dimensional affine fibre plane. It also gives
  an explicit nonrational algebraic formal graph with constant determinant,
  with exact all-order nontermination and marked-point ramification.
- **HC4MYGALL1:** Tang's [reproduced structural alternative](../tang-hc4-audit/STRUCTURAL_REVIEW.md)
  then forces invertibility for any polynomial graph with unit Hessian
  determinant. This uses no plane degree bound and does not assume JC2.

The result replaces the previous degree-48 marked-graph frontier. It does
not exclude every unmarked constant-determinant graph or arbitrary nonlinear
reparameterizations, other graph charts, mixed Schur pivots or coisotropic
reductions. OP-HC4-SCHUR remains open in those mechanisms.

[Exact checker](verify.py) and [actual replay receipt](receipt.json):

```sh
timeout 60s .venv/bin/python research/jacobian-followup/hc4-polynomial-graphs/verify.py
```

The all-degree arguments are written proofs; finite normal-series checks
are regressions. No formal verification, external review or novelty
priority is claimed. Mathematical status is authoritative in
[MATH_STATUS.json](../../MATH_STATUS.json).
