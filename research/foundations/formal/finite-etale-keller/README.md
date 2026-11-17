# Finite etale Keller fibres in Lean

This project formalizes the characteristic-zero construction in [Common Arithmetic Fibres](../../../../unpublished/common-arithmetic-fibers/README.md). The paper certificate proves automatic admissible translation, the arbitrary-degree polynomial map, its Jacobian-one normalization and degree bound, the literal full fibre and naturality, and the passage to abstract finite etale algebras.

The publication entry point is [PaperCertificate.lean](FiniteEtaleKeller/PaperCertificate.lean). It imports the focused theorem chain rather than every companion development.

```sh
lake build FiniteEtaleKeller.PaperCertificate
```

The toolchain and Mathlib versions are pinned in this project. `lake build FiniteEtaleKeller` also builds the broader companion modules.

Formal boundaries remain explicit:

- The broader supplied-presentation proposition in characteristic different from two is not formalized end to end.
- Selberg--Delange counting, the analytic first-prime-moment input and the classical rank-two exclusion remain outside the focused certificate.
- Geometric bridges for monodromy, base change and stable atomicity are not supplied by their abstract algebraic interfaces.
- Companion HC4 and Image finite-sum modules prove their stated algebraic portions; they do not establish the unrestricted conjectures.

The source modules are the authority for individual formal statements and assumptions.
