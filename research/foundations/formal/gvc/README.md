# Generalized Vanishing in Lean

This project supplies formal proofs for the [current GVC paper](../../../../results/generalized-vanishing-two-variables/README.md); the [deposited edition](../../../../results/generalized-vanishing-two-variables/published/README.md) is preserved separately. Lean and Mathlib versions are pinned locally.

```sh
lake build
```

The ternary counterexample is formalized from its literal polynomial and differential symbol through the pure and mixed identities. Coefficient base change and unused-variable padding prove failure in every characteristic-zero dimension at least three. The complete winding--profile--radial counterexample family is also formalized under its stated hypotheses.

The binary theorem is **partially formalized**. Hall localization and supporting algebraic interfaces are present, but the complete moving-envelope and binary proof are not an end-to-end Lean theorem. The full dimensional classification must therefore not be described as completely machine-checked.

Read [GVC.lean](GVC.lean) and the modules in [GVC/](GVC/) for the exact exported statements and assumptions.
