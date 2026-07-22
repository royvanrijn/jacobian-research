# Audited downstream consequences

> **Archived scope.**  This was the 20 July consequence audit.  Christopher
> D. Long's 20--21 July preprints subsequently supplied direct explicit
> three-variable GMC, `(xz)`, and `SU(2)` witnesses.  The current consequence
> graph and the distinction between the nonexplicit `GMC(158)` route and the
> direct `GMC(n>=3)` witness are maintained in
> [External consequences and provenance](../../extended-geometry/EXTERNAL_CONSEQUENCES_AND_PROVENANCE.md).

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
- **Zhao's Vanishing Conjecture: explicitly false in dimension 190.**  The Stable Normal-Form Consequences
  symmetric-gradient construction gives the 2012-term quartic HN polynomial
  `R` stored in `artifacts/generated-results/stable_normal_form_consequence_witnesses.json`.  Its gradient map has
  a transported collision, so Zhao's inversion formula proves
  `Delta^m R^(m+1)` is not eventually zero.
- **Special Image Conjecture: explicitly false in dimension 190.**  For this
  `R`, take `f=(sum zeta_i^2)R` and `g=R`.  Every `f^m` lies in the operator
  image, while `g f^m` fails to do so infinitely often.
- **Dixmier conjecture: explicitly false for the third Weyl algebra.**  The
  normalized foundational Keller map gives `Psi(x_i)=F_i` and
  `Psi(d_i)=sum_j(JF^(-1))_(ji)d_j`.  The exact Weyl relations hold, while an
  inverse would force a polynomial left inverse to the colliding map.

The first two bullets describe only the original route and remain
nonconstructive as written.  Long's later direct witnesses are explicit and
strictly stronger on Gaussian dimension.  The final three bullets here are explicit;
their formulas, hypothesis audit, and remaining minimality questions are in
[STABLE_NORMAL_FORM_AUDIT.md](../../extended-geometry/STABLE_NORMAL_FORM_AUDIT.md).
