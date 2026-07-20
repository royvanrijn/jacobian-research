# Audited downstream consequences

This note audits the claims in Zihan Zhang's 20 July 2026 expository page.

## Directly verified premise

The exact scripts verify the constant nonzero Jacobian and explicit collision.
Therefore `JC(3)` is false, and padding by identity coordinates makes `JC(n)`
false for every `n>=3`.

## Consequences verified by contraposition of published implications

- **Mathieu conjecture for `SU(3)`: false.** The cited fixed-dimension
  implication is `Mathieu(SU(N)) => JC(N)`. Taking `N=3` and contraposing gives
  the conclusion. The logical step is valid; the fixed-dimensional implication
  is literature-dependent and does not construct explicit witness functions.
- **Gaussian Moments Conjecture: false in some finite dimension.** Derksen, van
  den Essen, and Zhao prove that its validity in every dimension implies the
  all-dimensional Jacobian conjecture. The counterexample negates the latter,
  so at least one finite-dimensional GMC statement fails. It does not identify
  the least dimension or an explicit pair of polynomials.
- **Zhao's Vanishing Conjecture: false in some finite dimension.** Zhao's paper
  states equivalence between the all-dimensional Jacobian conjecture and the
  quartic homogeneous vanishing conjecture, and identifies its hypothesis with
  Hessian nilpotency. Thus an existential quartic Hessian-nilpotent counterexample
  follows. No small explicit quartic follows automatically from the 3D formula.
- **Image Conjecture: false in some finite dimension.** The published implication
  is `Image Conjecture => Vanishing Conjecture`; since the latter fails in some
  dimension globally, the all-dimensional Image Conjecture cannot hold. This is
  again existential.

These conclusions are logically sound given the cited implication/equivalence
theorems. Only the first identifies a failing dimension (`SU(3)`); the others do
not. None of the implication theorems automatically supplies compact explicit
witnesses from the displayed map.

