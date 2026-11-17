# Lean formalization: support saturation

This package kernel-checks the module-theoretic core of
[`../../verified/SUPPORT_SATURATION_PRINCIPLE.md`](../../verified/SUPPORT_SATURATION_PRINCIPLE.md).
It is pinned to Lean `v4.32.1` and Mathlib `v4.32.1`.

<!-- status-consumer: SST1F 838e558b5fcb9d81 -->

Build it with:

```sh
lake build
```

## Checked statements

- `Ideal.primaryComponent M I`, Mathlib's union of the annihilator modules
  `0 :_M I^n`, vanishes when `I` contains an `M`-regular element.
- Avoidance of every associated prime produces such a regular element by
  finite prime avoidance.
- Vanishing of the primary component is equivalent to avoidance of every
  associated prime and to existence of a regular element in the ideal.
- If every associated prime is minimal over `ann(M)` and `I` avoids those
  minimal primes, then `Ideal.primaryComponent M I = bot`.
- A nonzero element annihilated by `I` forces a nonzero primary component;
  this checks the mechanism behind the counterexample `R ⊕ R/I`.
- For a presentation `F -> F/N`, vanishing of the primary component proves
  the elementwise saturation statement: if `I^n * x` lies in `N`, then
  `x` already lies in `N`.
- Conversely, elementwise saturation of `N` forces the quotient primary
  component to vanish.

The last statement is the associated-prime form of the repository theorem
“`M` is `S1` and `I` has positive relative height on `M`.” Mathlib does not
currently provide the module-local-depth interface needed to state that
geometric wrapper directly, so the wrapper remains an ordinary proof from
the standard equivalence between `S1` and absence of embedded associated
primes for finite modules over Noetherian rings.

## Deliberate formalization boundary

The conductor theorem is a second layer: it glues compatible primitives
through the exact conductor square and thereby puts a defect in the primary
component. The derived theorem additionally uses the local-cohomology
spectral sequence, and the Rees corollary applies the module theorem to a
separately verified strict Rees presentation. These layers are not claimed
by this Lean package.  The checked application routing, including the exact
boundary between this formal core and the cubic, plane-JC, and restricted-Weyl
inputs, is recorded in
[`../../verified/SUPPORT_SATURATION_PATHS.json`](../../verified/SUPPORT_SATURATION_PATHS.json).

There are no `sorry`, `admit`, or explicit `axiom` declarations.
