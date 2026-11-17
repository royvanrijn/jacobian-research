# Verification matrix

This file records the proof layer for every load-bearing claim in
*Quantitative Hasse-Principle Failures in the Fibers of a Fixed Keller Map*.
“Lean” means a theorem in the pinned [finite-etale Keller project](../../research/foundations/formal/finite-etale-keller/)
with no `sorry`, `admit`, or project-specific axiom. Exact symbolic and
integer checks are separate audits; bounded enumeration is never used as
proof of an asymptotic.

The isolated publication entry point is
`FiniteEtaleKeller/FixedHassePaperCertificate.lean`.

The rational extension has a written proof and a fresh finite audit. The
retained Lean statements below cover their original sufficient integer
family and algebraic assertions; they do not formalize the rational
converse, general rational height or rational asymptotic.

| Claim | Paper proof | Lean certificate | Symbolic audit / external input |
|---|---|---|---|
| The displayed map is exactly `Φ = L ∘ F₀ ∘ A` | Definition (2.5) | `FixedHasseMap.lean`; definition `paperMap` | `verify_infinite_hasse_keller_fibers.py` expands the source and target scalings |
| `det DΦ = 1` | Chain rule from the determinant `-2` base map and the two diagonal transformations | `jacobianDet_paperMap` | Exact symbolic Jacobian check |
| The exact displayed map has geometric degree five | Function-field reconstruction for `F₀`; source and target polynomial automorphisms preserve degree | `baseMap_geometricDegree` is the complete function-field computation; `paperMap_normalization_inverse` formally recovers `baseMap` through the displayed automorphisms; `paperMapGeometricDegree` records the degree transported along that normalization. The elementary invariance of function-field degree under polynomial automorphisms is not separately re-encoded as a second scaled-field `AlgEquiv` | The base inverse polynomial has irreducible generic degree five |
| The moving inverse equation is `((S-1/2)^3-a)(S^2+3/4)` | Centered identity (2.17)–(2.18) | `inversePolynomial_eq_centered` and `inversePolynomial_eq_translate` | Exact symbolic factorization |
| The literal paper fiber is represented by `Q[X]/((X^3-a)(X^2+X+1))`, over every commutative rational test algebra | Reconstruction in Section 2, including derivative-unit argument | `paperFiberRepresentingEquiv` | Quotient reconstruction in `verify_infinite_hasse_keller_fibers.py` |
| For `a ≠ 0,1` the quotient is reduced, finite étale, and rank five | Discriminants and resultant in (2.20) | `polynomial_separable`, `quotient_etale`, and `quotient_rank` | Exact discriminant and resultant checks |
| An admissible parameter gives no rational point | Neither the cubic nor cyclotomic factor has a rational root | `paperFiberPoint_rat_isEmpty` | Rational-root audit |
| An admissible parameter gives a real point | The cubic factor has a real root | `paperFiberPoint_real_nonempty` | Elementary continuity argument |
| An admissible parameter gives a point over every `Q_p` | Place-by-place proof: exceptional primes `2,3`, cyclotomic factor for `p = 1 mod 3`, and bijective cubing for `p = 2 mod 3` | `polynomial_has_padic_root` and `paperFiberPoint_hasse_certificate` | Exact checks of the exceptional substitutions |
| The target has primitive coordinates `[9:-9:32a:24a+3]` and height `32a` | Elementary gcd and maximum calculation | `targetProjectiveContent_eq_one` and `targetProjectiveHeight_eq` | Integer checker |
| Different parameters give different targets | The second target coordinate is linear with nonzero slope | `rationalTarget_injective` | Direct rational arithmetic |
| The congruence-and-support core is multiplicatively closed | Elementary congruence and prime-support argument | `HasseCoreCondition.mul` and `HasseCoreCondition.pow` | Integer checker |
| The complete admissible set is not claimed to be a semigroup | Perfect cubes are removed after counting the multiplicative core | `AdmissibleHasseParameter` keeps the noncube condition separate | Example: `19` and `19²` are admissible, while `19³` is excluded |
| Every prime `ℓ = 1 mod 9` is admissible | Prime support and valuation obstruction to being a rational cube | `prime_not_rational_cube` and `primeParameter_certificate` | Prime examples in the symbolic checker |
| Full-family asymptotic `G₃(1) X/(3√π√log X)` | Euler product, Selberg–Delange, and explicit character orthogonality on `{1,4,7}` | Not formalized | External analytic input: character-twisted Selberg–Delange |
| Clean-family asymptotic `G₉(1) X/(Γ(1/6)(log X)^(5/6))` | Euler product with exponent `1/6` and Selberg–Delange | Not formalized | External analytic input: Selberg–Delange |
| Degree five is minimal | Low-degree factorization cases plus the Berend–Bilu fixed-point criterion / Chebotarev | Not formalized | External arithmetic input: Berend–Bilu |
| Counts through `a = 10^6` | Not used as proof | Not applicable | Historical dependency-free enumeration with a recorded JSON digest; the checker and artifact are absent from the current checkout |

## Rational extension: proof and current finite audit

| Claim | Current proof | Assurance |
|---|---|---|
| Every Hasse-failing rational parameter on the entire line | `thm:rational-classification` in `sections/03-local-global.tex`, including the exceptional fibres at `a=0,1` | Written local proof; not formally verified |
| Exact rational projective content `gcd(n,8) gcd(m,3)` and height polygon | `eq:rational-height` and `eq:rational-area` in `sections/04-rational-counting.tex` | Elementary written argument, direct gcd and independent denominator clearing on 20,086 reduced signed pairs |
| Coprime-pair density and coupled cubic-character condition | `eq:rational-sd` through `eq:rational-polygon-count` | Inspected de la Bretèche–Tenenbaum Theorem 1.2; convolution, uniform tail bound and polygon approximation supplied in the paper |
| Full rational count `C B²/log B`, with `C=11/(288√3 π²) ∏(1+p⁻³)/(1−p⁻³)` over inert primes | `thm:rational-count`; valuation summation, content weights and global-cube removal | Written analytic proof; finite factor algebra and enumeration are separate audits |

Run the standard-library checker from the repository root:

```sh
python3 unpublished/fixed-map-hasse-failures/verify_rational_line.py --maximum-height 10000 --output unpublished/fixed-map-hasse-failures/rational_line_verification.json
```

The retained JSON identifies the actual checker bytes and Python version.
It checks all 20,086 reduced signed pairs with numerator absolute value and
positive denominator at most 128, comparing heights with independent
denominator clearing and local conditions with unit-cube residues at all
primes at most 127. Complete height enumeration yields 617, 4377 and 41415
targets at bounds 1000, 3000 and 10000. A separate naive enumeration agrees
through height 80. Exact arithmetic also checks both sign areas and the
weights `72/7` and `18/13`. These are internal finite checks only.

Two contextual agents independently reconstructed the local/height
calculations and constant during this revision. Their checks are internal
AI review; agreement does not substitute for the written proof or external
human review. Review corrected an initial misuse of the root-scheme
description at `a=0`: the repeated zero root is absent from the actual
fibre, whose quadratic component has no real point.

The manuscript was rebuilt with the command below. No new Lean build is
claimed for this revision.

## Internal artifact review

An internal PDF reader reviewed the 16-page edition with SHA-256
`1487360cd37327e4dddcdda3777695c12b3001b07f3fe7a2daa8db752c8146ae`.
It received only the PDF, but subsequently saw two agents' short completion
summaries while checking collaboration availability. Its review is therefore
not labelled strictly blind or independent. It found no concrete error in
the rational classification/count, identified the Selberg–Delange interface
as the main imported step, and confirmed that the limits on targets versus
algebra classes and on formal verification are readable from the PDF.

The reader flagged the inherited fixed-map assertion's implication for the
Jacobian conjecture. The introduction now states that consequence explicitly
and distinguishes it from the rational root-scheme arithmetic. This is an
exposition correction, not a new external validation of the geometric claim.
The revised PDF was rebuilt and all 16 pages were rendered and inspected;
the final LaTeX log has no unresolved references or overfull boxes.

External human review, priority review and a new Lean replay remain
unperformed. Internal AI artifact review does not upgrade any of these
assurance flags.

An additional **isolated internal AI blind reader** received only the PDF
and did not inspect repository sources, prior messages or other agents.
It completely read the edition with hash `1487360cd37327e4dddcdda3777695c12b3001b07f3fe7a2daa8db752c8146ae`,
then checked the revised introduction in edition
`758e9f96e3ef8bbd6d3c04267c4b9dee3c39b8e4d0cb08366820fd66ffd6da86`.
It identified the fixed-map/rational-classification/exact-height count as
the contribution; the split/inert factor coverage and two sparse integer
counts as the proof mechanism; and the coprime-pair polygon limit plus
uniform valuation tails as the hardest steps. It found no definite error
in the rational argument, independently reproduced `N(1000)=617` from the
PDF, and identified the limits on line targets, algebra classes, integrality
and formal coverage correctly. Its PDF-only remit could not audit the
external Selberg–Delange source; the separate source audit above does so.

That reader found a precision gap in the inherited rank-versus-degree
argument: Zariski shrinking alone need not produce distinct finite étale
branches. Section 6 now uses the complex inverse function theorem and
intersects a common Euclidean neighborhood with the generic Zariski open.
The same isolated reader rechecked the repaired passage and found the
objection resolved. The final 16-page PDF has SHA-256
`6c4464de50e7f0a635d0928fa1e34ed2d11a919df6d1d0646b7f1314b9b8b8ce`.
All pages were rendered for layout inspection. The final build, exact
rational audit, Python compilation, focused `git diff --check`, and all
59 tooling tests passed. An earlier `make check` also passed in full.
The last global run stopped on an unrelated concurrent source-pin change
in `research/moments/vanishing/GAUSSIAN_AFFINE_ISOTROPIC_CHART.md` (two
declarations). That document was not changed by this arithmetic task.
All 39 rational-manuscript pin declarations across 13 files were separately
checked against their current bytes. These remain internal checks.

## Current build commands

From the repository root, build the focused Lean certificate with:

```sh
cd research/foundations/formal/finite-etale-keller
lake build FiniteEtaleKeller.FixedHassePaperCertificate
```

Build the manuscript from the repository root with:

```sh
latexmk -cd -pdf -jobname=paper -interaction=nonstopmode -halt-on-error \
  unpublished/fixed-map-hasse-failures/main.tex
```

The historical Python audits named in the matrix,
`verify_infinite_hasse_keller_fibers.py` and
`verify_multiplicative_hasse_artifact.py`, and their pinned count artifact
are not present in the current checkout. Their recorded checks are not
current replay instructions. The current rational audit and manuscript
build are recorded above; no fresh replay of the historical symbolic
checkers or Lean build is claimed. The finite count is not used as proof
of the asymptotic.
