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
- **Zhao's Vanishing Conjecture: explicitly false in dimension 190.**  The C15
  symmetric-gradient construction gives the 2012-term quartic HN polynomial
  `R` stored in `results/c15_consequence_witnesses.json`.  Its gradient map has
  a transported collision, so Zhao's inversion formula proves
  `Delta^m R^(m+1)` is not eventually zero.
- **Special Image Conjecture: explicitly false in dimension 190.**  For this
  `R`, take `f=(sum zeta_i^2)R` and `g=R`.  Every `f^m` lies in the operator
  image, while `g f^m` fails to do so infinitely often.
- **Dixmier conjecture: explicitly false for the third Weyl algebra.**  The
  normalized C01 map gives `Psi(x_i)=F_i` and
  `Psi(d_i)=sum_j(JF^(-1))_(ji)d_j`.  The exact Weyl relations hold, while an
  inverse would force a polynomial left inverse to the colliding map.

The first two bullets remain nonconstructive.  The final three are explicit;
their formulas, hypothesis audit, and remaining minimality questions are in
[C15_INDEPENDENT_AUDIT.md](C15_INDEPENDENT_AUDIT.md).
