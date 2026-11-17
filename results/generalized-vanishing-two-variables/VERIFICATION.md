# Verification matrix

This manuscript is a conceptual proof, not a computer-assisted theorem.
The matrix records where each load-bearing step is proved and what it uses.

| Claim | Paper location | Input |
|---|---|---|
| Degree at most the lowest positive operator order | Proposition 2.3 | Split-symbol polarization and Duistermaat--van der Kallen |
| Initial strict horizontal gap | Lemma 3.1 | Complete polarization, Duistermaat--van der Kallen, Hall's marriage theorem |
| Factorial valuation at a good prime | Lemma 3.2 | Legendre's formula at an unramified prime ideal |
| Shifted unequal-weight faces cannot overlap | Proposition 4.1 | Number-field specialization, Frobenius, Lemma 3.2, and unique-minimum valuation |
| Unequal-weight pure-zero faces are horizontally disjoint | Corollary 4.3 | Proposition 4.1 |
| A common unequal threshold is terminal for every fixed multiplier | Proposition 5.1 | Equal-face coordinate gap and bounded weight defect |
| The global envelopes reach a common threshold | Section 6 | Finite piecewise-linear envelopes and the initial Hall gap |
| Finite support certificate over the original field | Section 6.1, Binary support certificate | Unique Hall direction, strict refinement at the threshold, and the mediant weight |
| Sharp operator-independent cutoff | Section 6.1, Optimal operator-independent cutoff | Positive separator margin and the explicit monomial sharpness family |
| Unrestricted GVC in two variables | Theorem 1.1 | The preceding rows |
| Long's adjacent-coefficient mechanism | Section 7 | Long's Gaussian moment identity in the cited public preprint |
| Homogeneous GVC(3) counterexample | Theorem 8.1 | Quadric identity, exact spherical coefficients, and the Reynolds--apolar identity |
| Full winding--admissible-profile--radial failure family | Theorem 9.1 | The same phase extraction with arbitrary endpoint contact and admissible profile |
| `GVC(n)` and homogeneous `GVC(n)` hold exactly for `n <= 2` | Theorem 10.2 | Theorems 1.1 and 8.1, de Bondt's split-symbol theorem, and Lemma 10.1 |

No bounded search, generated artifact, or symbolic certificate is used in
the proof.  The Ferrers support statement in Section 12 is a consequence of
horizontal separation on unequal-weight walls with operator weight below phase weight, with the initial endpoint coefficients retained.  It does not prohibit overlap when the operator face has larger weight.  The dependency-free GVC(3) checker is an independent
bounded replay, not the all-order proof.

The Lean audit in [`research/foundations/formal/gvc`](../../research/foundations/formal/gvc) checks the
literal ternary polynomials and their degree-twelve homogeneity, the cusp
identity, the algebraic beta evaluation of the endpoint moment, the complete
binomial coefficient ladder for an arbitrary rational endpoint profile, the
literal multivariate profile family and its degree/order formula, and
positivity of the exact mixed scalar.  It also proves the coefficientwise
equal-degree apolar contraction, composition of differential symbols,
the full concrete Reynolds/Laurent phase extraction, coefficientwise
algebraic change of height variable, and identification with the endpoint
kernel.  The resulting bridge is constructed in Lean, yielding the pure
identity, exact mixed scalar, and literal ternary counterexample.  Coefficient
base change and unused-variable padding then prove GVC failure over every
characteristic-zero field in every finite dimension at least three.  Thus
Theorem 8.1 and the negative half of Theorem 10.2 are fully Lean-verified.
For Theorem 9.1 it also proves the literal arbitrary-profile quadric
restriction under the stated declared-degree condition, general even-phase
coefficient extraction, coefficientwise height substitution, endpoint
factorization, and the shifted primitive identity.  The resulting
`verifiedProfileFamilyBridge` yields the pure identity, complete exact
multiplier ladder, and GVC counterexample over `Q`; injective coefficient
base change proves the same failure over every characteristic-zero field
from exactly the theorem's nonzero-moment hypothesis.  Thus Theorem 9.1 is
fully Lean-verified.
Lean also proves abstractly that binary GVC implies unary GVC.  On the binary
side it proves the rational `p`-adic lower and exact factorial valuations
from Lemma 3.2 and the negative-final-slope plus intermediate-value envelope
step.  For Lemma 3.1 it now proves the complete finite Hall-deficiency core:
Mathlib's Hall theorem produces a deficient set, binary linear algebra puts
that set in one direction class, and exact counting gives at least
`d - e + 1` annihilated polynomial factors.  The split-factor interface
records the literal symbol and polynomial products, and Lean proves that the
full direction power divides the symbol while the matching perpendicular
power divides the polynomial.  These are the coordinate-free versions of
the displayed normal forms.  It also proves the
complete finite-support deduction in Proposition
5.1: finite maximization turns an ordering of the equal faces into an
integral coordinate cut, which refines the common threshold to a strict
positive weight whose linearly growing power gap eventually dominates every
fixed multiplier.  Empty equality faces are handled separately.

The audit still exposes, but does not yet prove, the translated
polarization/Duistermaat--van der Kallen implication that rules out a
split-factor matching, the literal field-extension and coordinate-change
and exact-quotient packaging of Lemma 3.1, or the
shifted-ray/no-reversal/global-envelope
obligations for the binary theorem.  In particular, the transfer of the
rational valuation statement to an unramified prime ideal and the resulting
`delta = 0` equal-face ordering remain part of the written number-field
proof.  A spherical integral
realization of the coefficientwise Reynolds functional is not checked, but
the formal differential identity uses the algebraic top contraction directly
and does not require that realization.
Therefore Theorem 1.1 and the complete biconditional in Theorem 10.2 should
not yet be described as fully Lean-verified.

Build the current manuscript from the repository root with:

```bash
latexmk -cd -pdf -interaction=nonstopmode -halt-on-error \
  -jobname=generalized-vanishing-two-variables results/generalized-vanishing-two-variables/main.tex
cd research/foundations/formal/gvc
lake build
```

The Lean command requires the project's pinned toolchain and dependencies.
It was not replayed during the 17 November 2026 editorial review. The older
standalone Python checker commands are not present in the current checkout;
the all-order written argument and retained Lean sources are the proof
support. The current [verification record](verification.json) records source
and PDF identities and the limited scope of this review.

## Prior release-review corrections, 17 November 2026

The Ferrers consequence now explicitly assumes `W_A < W_B` and retains the initial endpoint terms. The example `lambda = X^2 + Y^4`, `P = x` demonstrates why the unrestricted paragraph was false: all pure contractions vanish, but at weight `(2,1)` its projected faces overlap with operator weight four above phase weight two. This local consequence is corrected without extending the theorem's hypotheses or the formalization's scope.

The shifted-ray proof now uses containment of specialized Newton sets, preserving only the two selected endpoint coefficients. The horizontal-separation proof supplies the rationality argument for the vertical shift, including the two singleton case. The Reynolds identity identifies the one-dimensional invariant dual space before concluding proportionality. The abstract specifies an admissible profile, and the Feng--Sun reference uses the inspected published record, [Czechoslovak Mathematical Journal 69 (2019), 1061--1068](https://articles.math.cas.cz/10.21136/CMJ.2019.0049-18).

The formal sources are publicly accessible at the [versioned GVC supplement](https://github.com/royvanrijn/jacobian-research/tree/c32b9b4d273cf7988e85184a0d5eae9c488727fd/research/foundations/formal/gvc), with the toolchain, dependency pins and build instructions. This citation identifies the existing source artifact; no new deposit, Lean replay, transitive axiom audit or external review is claimed.
