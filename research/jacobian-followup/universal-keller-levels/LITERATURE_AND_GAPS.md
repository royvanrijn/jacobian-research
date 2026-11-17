# Literature transfer: finite residues and affine finite-field lines

14 November 2026. The extension uses established foundational mathematics with explicit Keller constructions.

| Primary source | Exact use |
|---|---|
| [Stacks Project, Proposition10.35.19](https://stacks.math.columbia.edu/tag/00GB) and [Lemma10.35.20](https://stacks.math.columbia.edu/tag/00GC) | Finite-type models over localized integers have closed points with finite residue fields. |
| [Stacks Project, Section32.10](https://stacks.math.columbia.edu/tag/01ZL) | Morphisms and inverse identities of finitely presented affine schemes descend to a finite coefficient model. |
| [Stacks Project, Lemma49.3.1](https://stacks.math.columbia.edu/tag/0BJF) | A finite algebra over a field is etale exactly when its trace pairing is nondegenerate. |
| [Jensen--Ledet--Yui, Generic Polynomials, Section2.1 and AppendixA.1](https://library.slmath.org/books/Book45/files/book45.pdf) | Cubic generic-polynomial context and additive Hilbert90; multiplicative Hilbert90 and generic additive-polynomial theory are classical inputs to the affine-line realization argument. |

The finite-field counts in UABLEVEL1 come from a written four-stratum decomposition, not a search. Spreading a hypothetical affine-space isomorphism gives a finite-field contradiction. For active coefficients the nonzero-level polynomial has the same Euler-characteristic value as affine space; the point-count polynomial retains the missing information. Those levels are also smooth, rational and factorial with only constant units, so all of these simpler source tests pass. This avoids assuming that a formal polynomial in the Lefschetz class is automatically nonzero over an arbitrary field.

The Frobenius map and its affine-group Galois structure belong to classical additive-polynomial theory. We give direct proofs of the complete-fibre criterion and monodromy. The explicit source and target group actions give effective polynomial gluing on affine bundles, with a specified marked root fibre. We make no claim that generic affine-group polynomials, Hilbert90, or Frobenius perturbations themselves are new.

## What changes in the research programme

- The B=constant source gate is closed for independently promoted coefficients. In characteristic zero none of the three moving target coordinates has a stably affine source level.
- The exact characteristic-three exception produces a universal cubic plane map over every field, including finite fields. Together with the univariate discriminant obstruction it gives Keller target dimension two in the identical positive-characteristic version of the programme definition.
- Canonical marked descent is possible over the affine finite-field symmetry group using an affine source bundle without an origin. At q=3 this covers BS3 itself.
- The fixed-space obstruction extends to all characteristics for linear representations, including wild actions. The affine source without an origin explains exactly why the new construction is outside that obstruction.

## Boundaries retained

The source-level theorem does not cover relations among seed coordinates, nonzero fixed seed specializations, arbitrary nonlinear coefficient sections, or another completion. The prime-power plane map is universal only for algebras carrying an affine F_q-line structure; for q>=4 it is not universal for all rank-q etale algebras. The canonical object lives on bundles; it does not choose their trivializations as fixed affine spaces.

The characteristic-three result does not improve any characteristic-zero plane bound or solve JC(2). In characteristic zero finite affine-linear group actions have fixed points by averaging, so the modular source-without-origin mechanism does not give the same escape.

See the [proof](PROOF.md) and the [unified theory map](../universal-keller-receiver/THEORY_MAP.md). Independent mathematical review remains pending.
