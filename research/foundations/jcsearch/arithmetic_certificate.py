"""Portable proof objects for local-to-global arithmetic Keller compilation.

The compiler side may use SymPy, but the emitted JSON contains only integers,
canonical rational strings, sparse-polynomial metadata, and SHA-256 hashes.
The independent verifiers deliberately do not import this module.
"""
from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

import sympy as sp

from jcsearch.keller_fiber import (
    KellerFiberCompilation,
    compile_polynomial_to_keller_fiber,
)
from jcsearch.local_global import (
    LocalAlgebraStabilityCertificate,
    automatic_local_stability_certificate,
    rational_valuation,
    synthesize_monic_polynomial,
)


SCHEMA = "proof-carrying-arithmetic-compilation/v1"
SPEC_SCHEMA = "proof-carrying-arithmetic-specification/v1"
MAP_SERIALIZATION = "expanded-sparse-rational-json/v1"


def fraction_text(value: Any) -> str:
    """Return a unique, language-neutral representation of a rational."""
    rational = Fraction(value)
    if rational.denominator == 1:
        return str(rational.numerator)
    return f"{rational.numerator}/{rational.denominator}"


def parse_fraction(value: str | int) -> Fraction:
    """Parse the restricted rational syntax used by specifications."""
    if not isinstance(value, (str, int)):
        raise ValueError("rationals must be encoded as strings or integers")
    return Fraction(value)


def sympy_rational(value: Fraction) -> sp.Rational:
    return sp.Rational(value.numerator, value.denominator)


def polynomial_from_coefficients(
    coefficients: Sequence[str | int],
    variable: sp.Symbol,
) -> sp.Poly:
    """Build a polynomial from coefficients in ascending degree order."""
    expression = sum(
        sympy_rational(parse_fraction(coefficient)) * variable**degree
        for degree, coefficient in enumerate(coefficients)
    )
    return sp.Poly(expression, variable, domain=sp.QQ)


def polynomial_coefficients(polynomial: sp.Poly) -> list[str]:
    """Serialize all coefficients, including zeroes, in ascending order."""
    variable = polynomial.gens[0]
    return [
        fraction_text(polynomial.coeff_monomial(variable**degree))
        for degree in range(polynomial.degree() + 1)
    ]


def canonical_expanded_map(
    mapping: Sequence[sp.Expr],
    variables: Sequence[sp.Symbol],
) -> list[list[list[Any]]]:
    """Return the canonical sparse representation hashed by both verifiers."""
    coordinates: list[list[list[Any]]] = []
    for component in mapping:
        polynomial = sp.Poly(sp.expand(component), *variables, domain=sp.QQ)
        terms = [
            [list(exponents), fraction_text(coefficient)]
            for exponents, coefficient in polynomial.terms()
        ]
        terms.sort(key=lambda term: tuple(term[0]))
        coordinates.append(terms)
    return coordinates


def expanded_map_sha256(
    mapping: Sequence[sp.Expr],
    variables: Sequence[sp.Symbol],
) -> tuple[str, list[int]]:
    """Hash the fully expanded rational map in the canonical JSON encoding."""
    expanded = canonical_expanded_map(mapping, variables)
    payload = json.dumps(
        expanded, ensure_ascii=True, separators=(",", ":"), sort_keys=False
    ).encode("ascii")
    return hashlib.sha256(payload).hexdigest(), [len(part) for part in expanded]


def _real_isolating_intervals(polynomial: sp.Poly) -> tuple[int, list[dict[str, Any]]]:
    """Compute disjoint rational intervals, one for each real root."""
    intervals = sp.intervals(polynomial, eps=sp.Rational(1, 10_000))
    serialized = []
    for (left, right), multiplicity in intervals:
        serialized.append(
            {
                "left": fraction_text(left),
                "right": fraction_text(right),
                "multiplicity": int(multiplicity),
            }
        )
    return sum(item["multiplicity"] for item in serialized), serialized


def _stability_record(
    certificate: LocalAlgebraStabilityCertificate,
    claimed_precision: int,
) -> dict[str, Any]:
    return {
        "prime": certificate.prime,
        "model_coefficients_ascending": polynomial_coefficients(
            certificate.polynomial
        ),
        "discriminant": fraction_text(certificate.discriminant),
        "discriminant_valuation": certificate.discriminant_valuation,
        "universal_precision": certificate.coefficient_precision,
        "claimed_precision": claimed_precision,
        "claim": "same finite etale Q_p-algebra",
        "theorem": "monic coefficient congruence modulo p^(2*v_p(discriminant)+1)",
    }


def _normalize_action_certificate(
    action: Mapping[str, Any],
    variable: sp.Symbol,
) -> dict[str, Any]:
    """Normalize polynomial fields while leaving descriptive labels intact."""
    normalized = dict(action)
    normalized["prime"] = int(action["prime"])
    normalized_factors = []
    for factor in action["factors"]:
        factor_record = dict(factor)
        polynomial = polynomial_from_coefficients(
            factor["coefficients_ascending"], variable
        )
        factor_record["coefficients_ascending"] = polynomial_coefficients(polynomial)
        normalized_factors.append(factor_record)
    normalized["factors"] = normalized_factors
    return normalized


def build_certificate(specification: Mapping[str, Any]) -> tuple[dict[str, Any], sp.Poly]:
    """Compile one JSON local-field specification into a portable certificate."""
    if specification.get("schema") != SPEC_SCHEMA:
        raise ValueError(f"expected specification schema {SPEC_SCHEMA!r}")
    variable = sp.Symbol(str(specification.get("variable", "T")))
    inverse_variable = sp.Symbol("S")
    sources = sp.symbols("x y z")
    degree = int(specification["degree"])

    local_models: dict[int, tuple[int | None, sp.Poly]] = {}
    stability_data: list[
        tuple[LocalAlgebraStabilityCertificate, int]
    ] = []
    for record in specification["local_models"]:
        prime = int(record["prime"])
        polynomial = polynomial_from_coefficients(
            record["coefficients_ascending"], variable
        )
        automatic = automatic_local_stability_certificate(
            polynomial, variable, prime
        )
        raw_precision = record.get("precision", "automatic")
        precision = (
            automatic.coefficient_precision
            if raw_precision == "automatic"
            else int(raw_precision)
        )
        if precision < automatic.coefficient_precision:
            raise ValueError(
                f"precision at {prime} is below the universal stability radius"
            )
        local_models[prime] = (precision, polynomial)
        stability_data.append((automatic, precision))

    coefficient_intervals = [
        (parse_fraction(left), parse_fraction(right))
        for left, right in specification["real_coefficient_intervals"]
    ]
    if len(coefficient_intervals) != degree:
        raise ValueError("the real coefficient box must have one entry per degree")
    synthesis = synthesize_monic_polynomial(
        local_models, coefficient_intervals, variable
    )
    polynomial = synthesis.polynomial

    translation_value = specification.get("translation", "automatic")
    stable_parameter = specification.get("stable_parameter")
    if stable_parameter is not None and (
        isinstance(stable_parameter, bool)
        or not isinstance(stable_parameter, int)
        or stable_parameter < 0
    ):
        raise ValueError("stable_parameter must be a nonnegative integer")
    compilation = compile_polynomial_to_keller_fiber(
        polynomial,
        variable,
        translation=(
            None
            if translation_value == "automatic"
            else parse_fraction(translation_value)
        ),
        inverse_variable=inverse_variable,
        source_variables=sources,
        stable_parameter=stable_parameter,
    )

    irreducibility_prime = int(specification["irreducibility_witness_prime"])
    denominator = sp.ilcm(
        *[int(sp.denom(coefficient)) for coefficient in polynomial.all_coeffs()]
    )
    if denominator % irreducibility_prime == 0:
        raise ValueError("irreducibility witness prime divides a denominator")

    root_count, root_intervals = _real_isolating_intervals(polynomial)
    map_hash, term_counts = expanded_map_sha256(
        compilation.determinant_one_map, sources
    )
    crt = synthesis.certificate
    local_moduli = dict(crt.local_moduli)
    local_records = []
    for automatic, claimed_precision in sorted(
        stability_data, key=lambda item: item[0].prime
    ):
        record = _stability_record(automatic, claimed_precision)
        record["coefficient_modulus"] = local_moduli[automatic.prime]
        local_records.append(record)

    global_coefficients = polynomial_coefficients(polynomial)
    translated_polynomial = sp.Poly(
        compilation.polynomial.as_expr().subs(
            variable, compilation.translation + inverse_variable
        ),
        inverse_variable,
        domain=sp.QQ,
    )
    inverse_polynomial = sp.Poly(
        compilation.inverse_polynomial,
        inverse_variable,
        domain=sp.QQ,
    )
    seed_polynomial = sp.Poly(
        compilation.seed, inverse_variable, domain=sp.QQ
    )
    g1 = seed_polynomial.coeff_monomial(inverse_variable)

    certificate = {
        "schema": SCHEMA,
        "name": str(specification["name"]),
        "provenance": {
            "generator": "scripts/compile_arithmetic_keller_certificate.py",
            "input_schema": SPEC_SCHEMA,
            "arithmetic_backend": "SymPy 1.14 exact QQ arithmetic",
            "regeneration_command": (
                ".venv/bin/python "
                "scripts/compile_arithmetic_keller_certificate.py"
            ),
        },
        "base_field": "Q",
        "variable": str(variable),
        "degree": degree,
        "local_models_and_precision_claims": local_records,
        "coefficient_crt": {
            "base_denominator": crt.base_denominator,
            "local_moduli": [
                {"prime": prime, "modulus": modulus}
                for prime, modulus in crt.local_moduli
            ],
            "crt_modulus": crt.crt_modulus,
            "coefficient_residues_ascending": list(crt.coefficient_residues),
            "multiplier": crt.multiplier,
            "common_denominator": crt.common_denominator,
            "real_coefficient_intervals": [
                [fraction_text(left), fraction_text(right)]
                for left, right in coefficient_intervals
            ],
        },
        "global_polynomial": {
            "monic_coefficients_ascending": global_coefficients,
            "primitive_integer_coefficients_ascending": [
                int(coefficient * denominator)
                for coefficient in [
                    Fraction(value) for value in global_coefficients
                ]
            ],
            "common_denominator": int(denominator),
        },
        "real_roots": {
            "count": root_count,
            "isolating_intervals": root_intervals,
        },
        "irreducibility_witness": {
            "prime": irreducibility_prime,
            "method": "irreducible reduction of a primitive polynomial",
        },
        "selected_translation": fraction_text(compilation.translation),
        "keller_map": {
            "normalization": "diag(1,-1/2,1)",
            "jacobian_determinant": "1",
            "coordinate_degrees": list(compilation.coordinate_degrees),
            "geometric_degree": compilation.geometric_degree,
            "expanded_map_hash": {
                "algorithm": "sha256",
                "serialization": MAP_SERIALIZATION,
                "digest": map_hash,
                "term_counts": term_counts,
            },
        },
        "inverse_polynomial_identity": {
            "seed_coefficients_ascending": polynomial_coefficients(seed_polynomial),
            "linear_coefficient": fraction_text(g1),
            "formula": "seed - linear_coefficient * target[2] / 2",
            "inverse_coefficients_ascending": polynomial_coefficients(
                inverse_polynomial
            ),
            "translated_input_coefficients_ascending": polynomial_coefficients(
                translated_polynomial
            ),
        },
        "target": [fraction_text(value) for value in compilation.target],
        "local_action_certificates": [
            _normalize_action_certificate(action, variable)
            for action in specification.get("local_action_certificates", [])
        ],
    }
    stable_certificate = compilation.stable_multiplicity
    if stable_certificate is None:
        certificate["lean_instantiation"] = {
            "module": str(specification["lean_module"]),
            "theorems": [
                "compiledAutomaticPageOne",
                "compiledMap_jacobianDet",
                "compiled_inversePolynomial",
                "compiledFiberEquiv",
            ],
        }
    else:
        certificate["keller_map"]["stable_multiplicity"] = {
            "family_parameter": stable_certificate.family_parameter,
            "gauge_exponent": stable_certificate.gauge_exponent,
            "separation_invariant": stable_certificate.separation_invariant,
            "separation_value": stable_certificate.separation_value,
            "fitting_support": [
                list(point) for point in stable_certificate.fitting_support
            ],
            "boundary_prime_count": stable_certificate.boundary_prime_count,
            "boundary_ramification_index": (
                stable_certificate.boundary_ramification_index
            ),
        }
        if "lean_module" in specification:
            certificate["lean_instantiation"] = {
                "module": str(specification["lean_module"]),
                "theorems": [
                    "compiledMap_jacobianDet",
                    "compiled_inversePolynomial",
                    "compiledFiberEquiv",
                ],
            }
    return certificate, polynomial


def write_certificate(certificate: Mapping[str, Any], path: Path) -> None:
    """Write deterministic, human-readable JSON."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(certificate, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _lean_rational(value: Any) -> str:
    rational = Fraction(value)
    if rational.denominator == 1:
        return f"({rational.numerator} : ℚ)"
    return (
        f"(({rational.numerator} : ℚ) / "
        f"({rational.denominator} : ℚ))"
    )


def _lean_polynomial_expression(coefficients: Sequence[Any]) -> str:
    terms = []
    for degree, raw_coefficient in enumerate(coefficients):
        coefficient = Fraction(raw_coefficient)
        if not coefficient:
            continue
        scalar = _lean_rational(coefficient)
        if degree == 0:
            term = f"Polynomial.C {scalar}"
        elif degree == 1:
            term = f"Polynomial.C {scalar} * Polynomial.X"
        else:
            term = f"Polynomial.C {scalar} * Polynomial.X ^ {degree}"
        terms.append(term)
    return "\n    + ".join(terms) if terms else "0"


def _lean_integer_polynomial_expression(coefficients: Sequence[int]) -> str:
    """Render integer coefficients as polynomial numerals understood by `ring`."""
    terms = []
    for degree, raw_coefficient in enumerate(coefficients):
        coefficient = int(raw_coefficient)
        if not coefficient:
            continue
        scalar = f"({coefficient} : Polynomial ℚ)"
        if degree == 0:
            term = scalar
        elif degree == 1:
            term = f"{scalar} * Polynomial.X"
        else:
            term = f"{scalar} * Polynomial.X ^ {degree}"
        terms.append(term)
    return "\n    + ".join(terms) if terms else "0"


def render_lean_instantiation(
    certificate: Mapping[str, Any],
    polynomial: sp.Poly,
) -> str:
    """Generate a concrete Lean specialization of the formal realization layer."""
    coefficients = [
        Fraction(value)
        for value in certificate["global_polynomial"][
            "monic_coefficients_ascending"
        ]
    ]
    variable = polynomial.gens[0]
    primitive_denominator = int(
        certificate["global_polynomial"]["common_denominator"]
    )
    primitive_coefficients = [
        int(value)
        for value in certificate["global_polynomial"][
            "primitive_integer_coefficients_ascending"
        ]
    ]
    primitive_polynomial = polynomial_from_coefficients(
        primitive_coefficients, variable
    )
    bezout_u, bezout_v, gcd = sp.gcdex(
        primitive_polynomial, primitive_polynomial.diff()
    )
    if gcd != 1:
        raise ValueError("the compiled polynomial is not squarefree")
    bezout_denominator = 1
    for coefficient in (*bezout_u.all_coeffs(), *bezout_v.all_coeffs()):
        bezout_denominator = sp.ilcm(
            bezout_denominator, int(sp.denom(coefficient))
        )
    u_coefficients = [
        int(
            bezout_denominator
            * bezout_u.coeff_monomial(variable**degree)
        )
        for degree in range(max(0, bezout_u.degree()) + 1)
    ]
    v_coefficients = [
        int(
            bezout_denominator
            * bezout_v.coeff_monomial(variable**degree)
        )
        for degree in range(max(0, bezout_v.degree()) + 1)
    ]
    translation = _lean_rational(certificate["selected_translation"])
    target_c = _lean_rational(certificate["target"][2])
    primitive_expression = _lean_integer_polynomial_expression(
        primitive_coefficients
    )
    u_expression = _lean_integer_polynomial_expression(u_coefficients)
    v_expression = _lean_integer_polynomial_expression(v_coefficients)
    lean_namespace = str(certificate["lean_instantiation"]["module"])
    scale_coefficient = _lean_rational(
        Fraction(primitive_denominator, bezout_denominator)
    )
    polynomial_coefficient = _lean_rational(
        Fraction(1, primitive_denominator)
    )
    bezout_coefficient = _lean_rational(Fraction(1, bezout_denominator))
    if int(certificate["degree"]) == 3:
        degree_tactic = (
            "compute_degree! <;> simp [Polynomial.coeff_one]"
        )
        bezout_numeral_tactic = "simp only [Polynomial.C_ofNat]"
        cubic_nonzero_tactic = (
            "norm_num [compiledPolynomial, primitivePolynomial, "
            "Polynomial.coeff_X,\n    Polynomial.coeff_one]"
        )
    else:
        degree_tactic = "compute_degree!"
        bezout_numeral_tactic = (
            "rw [Polynomial.C_ofNat 2, Polynomial.C_ofNat 3, "
            "Polynomial.C_ofNat 4]"
        )
        cubic_nonzero_tactic = (
            "norm_num [compiledPolynomial, primitivePolynomial, "
            "Polynomial.coeff_X]"
        )
    stable_record = certificate["keller_map"].get("stable_multiplicity")
    if stable_record is None:
        extra_import = ""
        compiled_map_expression = "generalGaugeJacobianOneMap compiledSeed"
        seed_degree_block = ""
        determinant_body = """\
  apply jacobianDet_generalGaugeJacobianOneMap
  · exact realizationSeed_linear_ne_zero compiledPolynomial selectedTranslation
      selectedTranslation_linear_ne_zero
  · exact realizationSeed_cubic_ne_zero compiledPolynomial selectedTranslation
      selectedTranslation_cubic_ne_zero"""
        fiber_codomain = """\
      GeneralGaugeJacobianOneFiberPoint
        compiledSeed 1 compiledTargetC A"""
        fiber_body = """\
  rw [compiledTargetC_eq]
  exact realizationJacobianOneFiberRepresentingEquiv
    compiledPolynomial selectedTranslation compiledPolynomial_squarefree
    selectedTranslation_linear_ne_zero selectedTranslation_cubic_ne_zero"""
        final_theorems = """\
theorem compiledAutomaticPageOne :
    AutomaticPageOneCertificate compiledPolynomial
      compiledPolynomial_squarefree (by
        rw [compiledPolynomial_natDegree]
        omega) :=
  automaticRealization_pageOne compiledPolynomial
    compiledPolynomial_squarefree (by
      rw [compiledPolynomial_natDegree]
      omega)

#print axioms compiledMap_jacobianDet
#print axioms compiled_inversePolynomial
#print axioms compiledFiberEquiv
#print axioms compiledAutomaticPageOne"""
    else:
        exponent = int(stable_record["gauge_exponent"])
        extra_import = "\nimport FiniteEtaleKeller.StableGaugeFiber"
        if int(certificate["degree"]) == 3:
            compiled_map_expression = (
                "cubicLiftGaugeJacobianOneMap compiledSeed "
                f"{exponent}"
            )
            seed_degree_block = """\

theorem compiledSeed_natDegree : compiledSeed.natDegree = 3 := by
  unfold compiledSeed
  rw [realizationSeed_natDegree compiledPolynomial selectedTranslation]
  · exact compiledPolynomial_natDegree
  · rw [compiledPolynomial_natDegree]
    omega
"""
            determinant_body = f"""\
  apply jacobianDet_cubicLiftGaugeJacobianOneMap
  · exact compiledSeed_natDegree
  · norm_num
  · exact realizationSeed_linear_ne_zero compiledPolynomial selectedTranslation
      selectedTranslation_linear_ne_zero
  · exact realizationSeed_cubic_ne_zero compiledPolynomial selectedTranslation
      selectedTranslation_cubic_ne_zero"""
            fiber_body = f"""\
  rw [compiledTargetC_eq]
  simpa [compiledMap, compiledSeed] using
    (cubicLiftGaugeRealizationFiberRepresentingEquiv (A := A)
      compiledPolynomial selectedTranslation {exponent} (by norm_num)
      compiledPolynomial_squarefree selectedTranslation_linear_ne_zero
      selectedTranslation_cubic_ne_zero)"""
        else:
            compiled_map_expression = (
                "powerShiftedGaugeJacobianOneMap compiledSeed "
                f"{exponent}"
            )
            seed_degree_block = ""
            determinant_body = """\
  apply jacobianDet_powerShiftedGaugeJacobianOneMap
  · exact realizationSeed_linear_ne_zero compiledPolynomial selectedTranslation
      selectedTranslation_linear_ne_zero
  · exact realizationSeed_cubic_ne_zero compiledPolynomial selectedTranslation
      selectedTranslation_cubic_ne_zero"""
            fiber_body = f"""\
  rw [compiledTargetC_eq]
  simpa [compiledMap, compiledSeed] using
    (powerShiftedGaugeRealizationFiberRepresentingEquiv (A := A)
      compiledPolynomial selectedTranslation {exponent}
      compiledPolynomial_squarefree selectedTranslation_linear_ne_zero
      selectedTranslation_cubic_ne_zero)"""
        fiber_codomain = """\
      StableGaugeFiberPoint compiledMap compiledTargetC A"""
        final_theorems = """\
#print axioms compiledMap_jacobianDet
#print axioms compiled_inversePolynomial
#print axioms compiledFiberEquiv"""
    return f"""\
/- This file is generated by scripts/compile_arithmetic_keller_certificate.py. -/
import FiniteEtaleKeller.PageOneTheorem{extra_import}

/-!
# Generated arithmetic Keller specialization

This file contains only the algebra-to-Keller layer.  The local-field, CRT,
real-isolation, and finite-field claims live in the portable JSON certificate
and are replayed by the two independent arithmetic verifiers.
-/

noncomputable section

open Polynomial

namespace {lean_namespace}

def primitivePolynomial : Polynomial ℚ :=
  {primitive_expression}

def compiledPolynomial : Polynomial ℚ :=
  Polynomial.C {_lean_rational(Fraction(1, primitive_denominator))} *
    primitivePolynomial

theorem compiledPolynomial_natDegree :
    compiledPolynomial.natDegree = {certificate["degree"]} := by
  unfold compiledPolynomial primitivePolynomial
  {degree_tactic}

def bezoutUInt : Polynomial ℚ :=
  {u_expression}

def bezoutVInt : Polynomial ℚ :=
  {v_expression}

theorem primitive_bezout :
    bezoutUInt * primitivePolynomial +
      bezoutVInt * primitivePolynomial.derivative =
        ({bezout_denominator} : Polynomial ℚ) := by
  unfold bezoutUInt bezoutVInt primitivePolynomial
  simp
  {bezout_numeral_tactic}
  ring

def bezoutU : Polynomial ℚ :=
  Polynomial.C {scale_coefficient} *
    bezoutUInt

def bezoutV : Polynomial ℚ :=
  Polynomial.C {scale_coefficient} *
    bezoutVInt

theorem compiled_bezout :
    bezoutU * compiledPolynomial +
      bezoutV * compiledPolynomial.derivative = 1 := by
  calc
    bezoutU * compiledPolynomial +
        bezoutV * compiledPolynomial.derivative =
      Polynomial.C {bezout_coefficient} *
        (bezoutUInt * primitivePolynomial +
          bezoutVInt * primitivePolynomial.derivative) := by
            unfold bezoutU bezoutV compiledPolynomial
            simp only [Polynomial.derivative_mul, Polynomial.derivative_C,
              zero_mul, zero_add]
            calc
              Polynomial.C {scale_coefficient} * bezoutUInt *
                    (Polynomial.C {polynomial_coefficient} *
                      primitivePolynomial) +
                  Polynomial.C {scale_coefficient} * bezoutVInt *
                    (Polynomial.C {polynomial_coefficient} *
                      primitivePolynomial.derivative) =
                (Polynomial.C {scale_coefficient} *
                    Polynomial.C {polynomial_coefficient}) *
                  (bezoutUInt * primitivePolynomial +
                    bezoutVInt * primitivePolynomial.derivative) := by ring
              _ = Polynomial.C {bezout_coefficient} *
                  (bezoutUInt * primitivePolynomial +
                    bezoutVInt * primitivePolynomial.derivative) := by
                    rw [← Polynomial.C_mul]
                    norm_num
    _ = Polynomial.C {bezout_coefficient} *
        ({bezout_denominator} : Polynomial ℚ) := by
          rw [primitive_bezout]
    _ = 1 := by
      change Polynomial.C {bezout_coefficient} *
        Polynomial.C ({bezout_denominator} : ℚ) = 1
      rw [← Polynomial.C_mul]
      norm_num

theorem compiledPolynomial_squarefree : Squarefree compiledPolynomial := by
  exact ((Polynomial.separable_def' compiledPolynomial).2
    ⟨bezoutU, bezoutV, compiled_bezout⟩).squarefree

def selectedTranslation : ℚ := {translation}

theorem selectedTranslation_linear_ne_zero :
    compiledPolynomial.derivative.eval selectedTranslation ≠ 0 := by
  norm_num [compiledPolynomial, primitivePolynomial, selectedTranslation]

theorem selectedTranslation_cubic_ne_zero :
    (Polynomial.hasseDeriv 3 compiledPolynomial).eval selectedTranslation ≠ 0 := by
  rw [show selectedTranslation = 0 by rfl]
  rw [← Polynomial.coeff_zero_eq_eval_zero]
  rw [Polynomial.hasseDeriv_coeff]
  {cubic_nonzero_tactic}

def compiledSeed : Polynomial ℚ :=
  realizationSeed compiledPolynomial selectedTranslation{seed_degree_block}

def compiledMap : Fin 3 → GaugePolynomial ℚ :=
  {compiled_map_expression}

def compiledTargetC : ℚ := {target_c}

theorem compiledTargetC_eq :
    compiledTargetC =
      realizationTargetC compiledPolynomial selectedTranslation
        (compiledPolynomial.derivative.eval selectedTranslation) := by
  norm_num [compiledTargetC, compiledPolynomial, selectedTranslation,
    primitivePolynomial, realizationTargetC]

theorem compiledMap_jacobianDet : jacobianDet compiledMap = 1 := by
{determinant_body}

theorem compiled_inversePolynomial :
    generalGaugeInversePolynomial compiledSeed 1 0 compiledTargetC =
      translatePolynomial compiledPolynomial selectedTranslation := by
  rw [compiledTargetC_eq]
  exact generalGaugeInversePolynomial_realization
    compiledPolynomial selectedTranslation selectedTranslation_linear_ne_zero

def compiledFiberEquiv
    (A : Type*) [CommRing A] [Algebra ℚ A] :
    (AdjoinRoot compiledPolynomial →ₐ[ℚ] A) ≃
{fiber_codomain} := by
{fiber_body}

{final_theorems}

end {lean_namespace}
"""


def write_lean_instantiation(
    certificate: Mapping[str, Any],
    polynomial: sp.Poly,
    path: Path,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        render_lean_instantiation(certificate, polynomial),
        encoding="utf-8",
    )
