# Conditional Lean supplement

The formal sources cited in the PDF are accessible at the [versioned GMC(2) project](https://github.com/royvanrijn/jacobian-research/tree/c32b9b4d273cf7988e85184a0d5eae9c488727fd/research/foundations/formal/gmc2). This source pin is a public repository edition, not a new deposit. Its README, `lean-toolchain`, `lakefile.toml` and `lake-manifest.json` identify the build and dependencies.

Use Lean `v4.32.1` and Mathlib `v4.32.1`, pinned in the manifest to commit `520045ab14e26149ee970e2e617ca04b09bde5d6`. From the project directory:

```sh
lake build GMC2
```

The exported endpoint is `GMC2.gaussianMomentsConjecture_two_variables_bivariate` in `GMC2/BivariateGaussian.lean`. It uses ordinary bivariate polynomials with the circular Wick expectation and asserts eventual mixed vanishing. The explicit cutoff and optimality are written results in the paper, not exported endpoints.

To reproduce its transitive axiom report, create `Audit.lean` in the project containing:

```lean
import GMC2
#print axioms GMC2.gaussianMomentsConjecture_two_variables_bivariate
```

Then run `lake env lean Audit.lean`. The recorded audit lists the mathematical inputs `GMC2.duistermaat_van_der_kallen` and `GMC2.exists_goodReduction`, besides Lean's foundational choice, quotient and propositional-extensionality axioms. The former supplies the constant-term theorem; the latter supplies good finite-characteristic reductions of finite-type characteristic-zero domains, preserving a selected nonzero coefficient. Neither mathematical input is proved inside this package.

The source interfaces and endpoint were inspected during this revision. No fresh Lean build or axiom audit was performed, and this supplement does not upgrade the historical audit to a new replay. Independent mathematical review remains pending.
