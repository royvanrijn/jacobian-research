# Lean formalization: finite étale Keller fibers

This project formalizes *Finite Étale Algebras as Keller Fibers* in stages.
It uses Lean `v4.33.0-rc1` and Mathlib at the matching release candidate.

## Completed in stage 1

The current files contain no `sorry` and certify:

- the displayed denominator-free quintic map has Jacobian determinant `-722`;
- its determinant-`-2` quadratic-gauge normalization and determinant-one output normalization;
- the exact output scaling `diag(1,19,19)` and target conversion `-2 -> -38`;
- the inverse-polynomial identity for
  `(T^3 - 19)(T^2 + T + 1)`;
- an explicit Bézout identity for the quintic and its derivative;
- the resulting derivative inverse in `Q[T]/(P_5)`.

## Planned stages

1. **Universal marked-line identities.** Formalize the determinant cancellation
   from the source chart and the marked-line map for a general seed.
2. **Scheme reconstruction.** Construct the two coordinate-ring homomorphisms
   and prove the full fiber algebra is `K[S]/(E)` when `E` is squarefree.
3. **Polynomial realization.** Formalize translation
   `G(S)=P(a+S)-P(a)` and the prescribed fiber theorem.
4. **Finite étale realization.** Add monogenicity over infinite fields.
5. **Rank classification.** Add the historical degree-two Galois theorem or
   retain it as an explicitly named external hypothesis until formalized.

## Build

```bash
cd formal/finite-etale-keller
lake build
```

The repository CI runs this build independently of the existing external Lean
certificate for the foundational three-dimensional map.
