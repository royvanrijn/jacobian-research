# Factorial moments and rigidity

<!-- Generated from research/MATH_STATUS.json; run make render-navigation. -->

Factorial moment functionals, composition rigidity and Strong Factorial families.

**Authored assessment:** Credited binary audits and the stated universal families are complete internally; unrestricted Factorial and Strong Factorial remain open.

| Current result | Evidence / publication | Graph status | Main source |
|---|---|---|---|
| Binary Factorial and composition rigidity reproductions | Internal result | proved | [Source](image-rigidity/LONG_EIT_AND_FC2_AUDIT.tex) |
| Restricted all-dimensional families | Internal result | proved | [Source](image-rigidity/SIMPLEX_PENCILS_AND_AFFINE_LIMITS.tex) |
| All ordinary constant-plus-power SFC | Internal result | proved | [Source](image-rigidity/ALL_UNIVARIATE_TWO_TERM_SFC.tex) |
| Other universal two-term SFC families | Internal result | proved | [Source](image-rigidity/UNIFORM_PADE_PARABOLA_SFC.tex) |
| Current general-gap algebraic criterion | Partial result | proved | [Source](image-rigidity/HIGHER_GAP_ADJOINT_RANK.tex) |
| Ordinary initial order five, Gamma initial order four and support reductions | Internal result | stale | [Source](image-rigidity/lower_exponent_initial_and_boundary.tex) |
| Arithmetic terminal criteria and a resultant stability obstruction | Internal result | proved | [Source](image-rigidity/ARITHMETIC_FRONTIER_TERMINAL_ORDERS.tex) |
| Gap-three whole-matrix certificates and uniform initial minors | Internal result | proved | [Source](image-rigidity/gap_three_hurwitz.tex) |
| All-gap discriminants, cyclic Hurwitz minors and auxiliary factorizations | Internal result | proved | [Source](image-rigidity/gap_three_hurwitz_structure.tex) |
| All prime-gap lower-one and degree-four binomial SFC | Internal result | proved | [Source](image-rigidity/ARITHMETIC_FRONTIER_TERMINAL_ORDERS.tex) |

**Computed obligations:** 22 unresolved; 0 resolved within their exact scope.


Scope limits:

- Unrestricted FC and SFC remain open. Universal special-family theorems do not supply an arbitrary-support reduction.
- Every prime-gap lower-one family X(A+BX^p) now satisfies ordinary SFC at every index, and every ordinary binomial of degree at most four is covered. The proof is an all-prime/all-index cyclic Fitting argument, not a finite prefix or a proof of the all-real-shape H3 total-nonnegativity target.
- The pencil and coupled-family descriptions retain the scale, Gamma-shape, block and ordinary coefficient-ring normalization restrictions stated in their proofs.
- Higher-gap rank and backward-pivot results are proved criteria. Uniform actual-family rank/nonzero-pivot assertions remain open; finite controls do not discharge them. Priority and independent review are separate obligations.
- New gap-three certificates cover entire infinite matrices for n=2..16 and all beta>=0; initial minors and pivots have separate unbounded-index scopes. Prime-power terminal families and arbitrary-support initial pairs do not close all univariate binomial SFC.
- Effective large-gap and large-lower-exponent bounds reduce each fixed ordinary terminal order to an explicit finite support rectangle; all rectangles as N varies are still unchecked. Fixed-shape auxiliary positive networks do not yet transfer through the signed cyclic Wronskian map.
- All lower exponents and gaps are separated through terminal five for ordinary factorials and through terminal four for every positive Gamma shape. Arithmetic criteria, all-j root-scale monotonicity and the diagonal derivative embedding do not discharge the remaining all-N>=6 rectangles or uniform gap-three H3 total nonnegativity.
- A uniform all-gap identity now normalizes the moment discriminant as an (n-1)-square cyclic Hurwitz determinant. Its nonnegativity is proved for every positive shape, giving the last three gap-three special minors nonnegative at all indices. Strictness is now proved at prime gaps for positive rational p-integral shapes; all-real-shape strictness and the remaining intermediate minors are open.
- The all-order formal Gamma-power cumulant fraction and its positive tree coefficients are recorded separately from fixed-shape Hankel nonnegativity. The latter retains an explicit PSZ theorem/application verification input; transferring either structure to consecutive binomial moments at two different shapes remains an open mathematical gap.

Tracked obligations (computed states):

- **mathematical / open:** [Higher-dimensional FC and unrestricted SFC.](README.md) — FACTORIAL-GENERAL-FC
- **mathematical / open:** [All lower-one higher-gap Gamma pairs at every terminal order; lower exponents greater than one require a separate argument.](README.md) — FACTORIAL-LOWER-ONE-HIGHER-GAP
- **verification / open:** [Independent mathematical and priority review.](README.md) — GAP-FACTORIAL-3
- **mathematical / open:** [Uniform total nonnegativity of the actual gap-three H3 family, for example via all-size staircase-minor positivity or a parameter-compatible positive factorization.](image-rigidity/gap_three_hurwitz.tex) — FACTORIAL-GAP3-HURWITZ-TN
- **mathematical / stale:** [All ordinary univariate binomial SFC, equivalently uniform consecutive coprimality in the finite support rectangles for every terminal N>=6.](image-rigidity/lower_exponent_initial_and_boundary.tex) — FACTORIAL-ALL-UNIVARIATE-BINOMIAL-SFC
- **mathematical / open:** [Prove consecutive coprimality of the actual diagonal derivatives in the fixed lower-one base sequence; generic seed simplicity and stability do not imply this.](image-rigidity/lower_exponent_initial_and_boundary.tex) — FACTORIAL-BINOMIAL-DIAGONAL-PAIRS
- **mathematical / open:** [Strict positivity of the actual terminal Hurwitz determinant at every lower-one gap and index; its uniform nonnegativity is proved.](image-rigidity/gap_three_hurwitz_structure.tex) — FACTORIAL-LOWER-ONE-TERMINAL-DETERMINANT
- **mathematical / open:** [Exclude zero in the final independent special minor for the actual coupled gap-three family; full H3 total nonnegativity remains a separate route.](image-rigidity/gap_three_hurwitz_structure.tex) — FACTORIAL-GAP3-TERMINAL-MINOR-STRICT

[Exact scopes and verification limits](../STATUS.md). Publication, formal verification and independent review are separate statuses. [Dependency routes](../DEPENDENCIES.md) and [scope audit](../SCOPES.md) expose remaining gaps.
