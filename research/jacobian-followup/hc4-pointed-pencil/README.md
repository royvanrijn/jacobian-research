# Pointed Keller pencils: weighted boundary budgets and a genus control

For the singleton-zero relative-linear branch isolated by
[HC4MYOR1/2](../hc4-relative-row-orbit/README.md), blowing up the origin
gives a smooth pencil with an exceptional section. The section accounts
for exactly one sheet; it does not account for the other sheets.

The [written proof](PROOF.md) records:

- **HC4MYPP1:** on the finite normalization of the target blowup, the
  omitted radial divisors satisfy sum(e_i*f_i)=d-1. Here e_i is radial
  ramification and f_i is degree over the exceptional projective line.
  A generic pencil curve also satisfies an exact Euler deficit budget.
- **HC4MYPP2:** even elementary matrices in the same oriented row orbit
  give every nonnegative generic pencil genus while preserving the same
  zero ideal and the same nonconstant Jacobian polynomial. These are
  controls, not Keller counterexamples.

The published Nguyen Van Chau pencil theorem requires irreducibility and
rationality of every pencil curve. Those hypotheses have not been proved
for our unrestricted Keller rows. The new budgets specify checks for a
complete boundary calculation; they do not supply that calculation.

[Checker](verify.py) · [replay](receipt.json) · [primary-source pins](sources/manifest.json).

```sh
timeout 60s .venv/bin/python research/jacobian-followup/hc4-pointed-pencil/verify.py
```

Geometric arguments are written proofs, not certified by the symbolic
checker. No formal verification, external review or novelty priority is
claimed. The singleton-zero problem, OP-HC4-SCHUR and JC2 remain open;
[MATH_STATUS.json](../../MATH_STATUS.json) is authoritative.
