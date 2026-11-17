# Verification matrix

This file records the proof layer for every load-bearing claim in
*Quantitative Hasse-Principle Failures in the Fibers of a Fixed Keller Map*.
“Lean” means a theorem in the pinned `formal/finite-etale-keller` project
with no `sorry`, `admit`, or project-specific axiom. Exact symbolic and
integer checks are independent audits; bounded enumeration is never used as
proof of an asymptotic.

The isolated publication entry point is
`FiniteEtaleKeller/FixedHassePaperCertificate.lean`.

| Claim | Paper proof | Lean certificate | Independent audit / external input |
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
| Counts through `a = 10^6` | Not used as proof | Not applicable | Dependency-free enumeration reproduces the pinned JSON and SHA-256 digest |

## Reproduction commands

From the repository root:

```bash
.venv/bin/python scripts/verify_infinite_hasse_keller_fibers.py
.venv/bin/python scripts/verify_multiplicative_hasse_artifact.py
cd formal/finite-etale-keller
lake build FiniteEtaleKeller.FixedHassePaperCertificate
```

Build the manuscript with:

```bash
latexmk -cd -pdf -interaction=nonstopmode -halt-on-error \
  papers/fixed-map-hasse-failures/main.tex
python3 scripts/check_latex_log.py \
  papers/fixed-map-hasse-failures/main.log
```
