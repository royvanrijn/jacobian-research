#!/usr/bin/env sage
"""Build the pre-solver T-arithmetic ledger for the rank-seven catalogue.

For each distinct catalogue (T,NS), this script first passes from the literal
ternary Gram matrix to the primitive integral quadratic form with the same
rational orthogonal group.  It then computes rational isotropy, the even
Clifford algebra, and the integral even-Clifford order discriminant.

The output is deliberately fail-closed.  An integral U-splitting certifies a
split Eichler order and hence Gamma_0(N), including exact cusps and genus.
For a general split order, or a division order not independently certified to
be Eichler, the order is still identified exactly but Gamma_0(N) / X_0^D(N)
labels are retained only as typed candidates.  The determinant-948 H3 row is
the pinned positive control and imports the exact (D,N)=(6,79), Atkin--Lehner,
and CM-anchor identifications already replayed by the dedicated scripts.

status: ACTIVE_SEARCH_INFRASTRUCTURE
claim: exact T-similarity normalization, rational isotropy, Clifford algebra,
  Clifford-order reduced discriminant, certified Gamma_0(N) curves for the
  integral-U split rows, and exact H3 arithmetic positive-control metadata.
output: artifacts/generated-results/elkies-k3-rank7-t-arithmetic-v1.json
"""

import argparse
import hashlib
import json
from collections import Counter
from math import gcd as python_gcd, isqrt
from pathlib import Path

from sage.all import (
    CliffordAlgebra,
    Gamma0,
    QQ,
    QuaternionAlgebra,
    ZZ,
    QuadraticForm,
    factor,
    gcd,
    kronecker_symbol,
    lcm,
    matrix,
    pari,
    prod,
    vector,
    xgcd,
)


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CATALOGUE = (
    ROOT / "artifacts/generated-results/elkies-k3-rank7-auxiliary-catalogue-v1.json"
)
DEFAULT_OUTPUT = (
    ROOT / "artifacts/generated-results/elkies-k3-rank7-t-arithmetic-v1.json"
)
H3_SURFACE_ID = "K3-8188cdcda8c57b2d"
H3_T = [[-2, 0, 1], [0, 4, 0], [1, 0, 118]]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path):
    return str(path.resolve().relative_to(ROOT))


def rows(value):
    return [list(map(int, row)) for row in value.rows()]


def factor_primes(value):
    return [int(prime) for prime, _ in factor(ZZ(value))]


def primitive_vector(values):
    denominators = [QQ(value).denominator() for value in values]
    scale = lcm(denominators)
    integral = vector(ZZ, [ZZ(scale * QQ(value)) for value in values])
    content = gcd([abs(value) for value in integral])
    integral = integral / content
    for value in integral:
        if value:
            return -integral if value < 0 else integral
    raise ArithmeticError("zero isotropic vector")


def primitive_integral_hessian(literal):
    """Return a primitive integral quadratic Hessian with O(H)=O(literal)."""
    content = gcd([abs(value) for value in literal.list()])
    primitive_rational = literal / content
    assert all(QQ(value).denominator() == 1 for value in primitive_rational.list())
    primitive = primitive_rational.change_ring(ZZ)
    quadratic_scale = 1 if all(primitive[i, i] % 2 == 0 for i in range(3)) else 2
    hessian = quadratic_scale * primitive
    assert all(hessian[i, i] % 2 == 0 for i in range(3))
    coefficients = [
        hessian[0, 0] // 2,
        hessian[1, 1] // 2,
        hessian[2, 2] // 2,
        hessian[1, 2],
        hessian[0, 2],
        hessian[0, 1],
    ]
    assert gcd([abs(value) for value in coefficients]) == 1
    return content, primitive, quadratic_scale, hessian


def bezout_vector(covector):
    a, b, c = map(ZZ, covector)
    first_gcd, first_a, first_b = xgcd(a, b)
    total_gcd, combine, last = xgcd(first_gcd, c)
    assert total_gcd == 1
    result = vector(ZZ, [combine * first_a, combine * first_b, last])
    assert result * vector(ZZ, covector) == 1
    return result


def integral_u_split(hessian, isotropic):
    covector = hessian * isotropic
    divisibility = gcd([abs(value) for value in covector])
    if divisibility != 1:
        return {
            "status": "UNKNOWN_NO_DIVISIBILITY_ONE_ISOTROPIC_VECTOR_FROM_EXACT_WITNESS",
            "isotropic_divisibility": int(divisibility),
            "basis_change_rows": None,
        }
    companion = bezout_vector(covector)
    square = ZZ(companion * hessian * companion)
    assert square % 2 == 0
    companion -= (square // 2) * isotropic
    assert isotropic * hessian * companion == 1
    assert companion * hessian * companion == 0
    orthogonality = matrix(ZZ, [hessian * isotropic, hessian * companion])
    complement_basis = orthogonality.right_kernel_matrix()
    assert complement_basis.nrows() == 1
    complement = complement_basis.row(0)
    complement /= gcd([abs(value) for value in complement])
    complement_square = ZZ(complement * hessian * complement)
    assert complement_square > 0 and complement_square % 2 == 0
    change = matrix(ZZ, [isotropic, complement, companion])
    assert abs(change.det()) == 1
    level = complement_square // 2
    standard = matrix(ZZ, [[0, 0, 1], [0, 2 * level, 0], [1, 0, 0]])
    assert change * hessian * change.transpose() == standard
    return {
        "status": "PASS_EXACT_INTEGRAL_U_SPLIT_EICHLER_ORDER",
        "isotropic_divisibility": 1,
        "basis_change_rows": rows(change),
        "standard_hessian": rows(standard),
        "eichler_level": int(level),
    }


def clifford_data(hessian):
    quadratic_form = QuadraticForm(QQ, hessian)
    diagonal = quadratic_form.rational_diagonal_form()
    diagonal_values = [QQ(diagonal[i, i]) for i in range(3)]
    hilbert_a = -diagonal_values[0] * diagonal_values[1]
    hilbert_b = -diagonal_values[0] * diagonal_values[2]
    quaternion = QuaternionAlgebra(QQ, hilbert_a, hilbert_b)
    ramified = tuple(map(int, quaternion.ramified_primes()))
    algebra_discriminant = int(prod(ramified))

    clifford = CliffordAlgebra(quadratic_form)
    keys = list(clifford.basis().keys())
    order_basis = [
        clifford.one(),
        clifford.gen(0) * clifford.gen(1),
        clifford.gen(0) * clifford.gen(2),
        clifford.gen(1) * clifford.gen(2),
    ]

    def even_coordinates(value):
        coefficients = value.monomial_coefficients()
        return vector(
            QQ,
            [coefficients.get(keys[index], 0) for index in (0, 4, 5, 6)],
        )

    def left_matrix(value):
        return matrix(
            QQ, [even_coordinates(value * basis) for basis in order_basis]
        ).transpose()

    def reduced_trace(value):
        return left_matrix(value).trace() / 2

    trace_pairing = matrix(
        QQ,
        4,
        4,
        lambda left, right: reduced_trace(order_basis[left] * order_basis[right]),
    )
    assert all(QQ(value).denominator() == 1 for value in trace_pairing.list())
    trace_pairing = trace_pairing.change_ring(ZZ)
    determinant = abs(trace_pairing.det())
    reduced_discriminant = determinant.sqrt()
    assert reduced_discriminant**2 == determinant
    assert reduced_discriminant % algebra_discriminant == 0
    local_level_index = reduced_discriminant // algebra_discriminant
    return {
        "rational_diagonal_quadratic_coefficients": [
            str(value) for value in diagonal_values
        ],
        "hilbert_symbol": [str(hilbert_a), str(hilbert_b)],
        "finite_ramified_primes": list(ramified),
        "quaternion_discriminant": algebra_discriminant,
        "integral_even_clifford_order": {
            "basis": ["1", "e0*e1", "e0*e2", "e1*e2"],
            "reduced_trace_pairing": rows(trace_pairing),
            "reduced_discriminant": int(reduced_discriminant),
            "local_level_index": int(local_level_index),
        },
    }


def modular_curve(level):
    group = Gamma0(level)
    cusps = []
    widths = []
    for cusp in group.cusps():
        width = int(group.cusp_width(cusp))
        widths.append(width)
        cusps.append(
            {
                "representative": str(cusp),
                "width": width,
            }
        )
    index = int(group.index())
    genus = int(group.genus())
    nu2 = int(group.nu2())
    nu3 = int(group.nu3())
    assert sum(widths) == index
    if genus == 0:
        if level == 1:
            hauptmodul = {
                "status": "PASS_EXACT_REGISTRY_HAUPTMODUL",
                "formula": "j(tau)",
                "normalization": "classical j-invariant",
            }
        elif level in (2, 3, 4, 5, 7, 9, 13, 25):
            exponent = 24 // (level - 1)
            hauptmodul = {
                "status": "PASS_EXACT_REGISTRY_HAUPTMODUL",
                "formula": f"(eta(tau)/eta({level}*tau))^{exponent}",
                "normalization": "q^(-1)+O(1) at Infinity",
            }
        else:
            hauptmodul = {
                "status": "SEARCH_TRIGGERED_GENUS_ZERO_NO_LOCAL_FORMULA_REGISTRY_HIT",
                "formula": None,
                "normalization": None,
            }
    else:
        hauptmodul = {
            "status": "NOT_APPLICABLE_POSITIVE_GENUS",
            "formula": None,
            "normalization": None,
        }
    torsion_free = nu2 == 0 and nu3 == 0
    if torsion_free and genus == 0 and index == 12:
        surface_match = "SEMISTABLE_RATIONAL_MODULAR_ELLIPTIC_SURFACE_NUMERICAL_MATCH"
    elif torsion_free and genus == 0 and index == 24:
        surface_match = "SEMISTABLE_MODULAR_ELLIPTIC_K3_NUMERICAL_MATCH"
    else:
        surface_match = "NO_DIRECT_SEMISTABLE_RATIONAL_OR_K3_INDEX_MATCH"
    return {
        "status": "PASS_EXACT_GAMMA0_CURVE",
        "group": f"Gamma_0({level})",
        "level": int(level),
        "index_in_PSL2Z": index,
        "genus": genus,
        "elliptic_orbits_order_2": nu2,
        "elliptic_orbits_order_3": nu3,
        "cusp_count": len(cusps),
        "cusps": cusps,
        "cusp_width_multiset": sorted(widths),
        "hauptmodul": hauptmodul,
        "explicit_modular_elliptic_surface_comparison": {
            "status": surface_match,
            "torsion_free_orbifold": torsion_free,
            "semistable_euler_number_if_torsion_free": index if torsion_free else None,
            "fibre_multiplicities_from_cusp_widths": sorted(widths),
            "equation_registry_match": "UNKNOWN_NOT_IDENTIFIED_BY_CUSP_DATA_ALONE",
        },
    }


def is_squarefree(value):
    return all(exponent == 1 for _, exponent in factor(ZZ(value)))


def is_fundamental_discriminant(value):
    if value >= 0:
        return False
    if value % 4 == 1:
        return is_squarefree(abs(value))
    if value % 4 == 0:
        quotient = value // 4
        return is_squarefree(abs(quotient)) and quotient % 4 in (2, 3)
    return False


def fundamental_discriminant_and_conductor(discriminant):
    absolute = abs(discriminant)
    for conductor in range(isqrt(absolute), 0, -1):
        if absolute % (conductor * conductor):
            continue
        field_discriminant = discriminant // (conductor * conductor)
        if is_fundamental_discriminant(field_discriminant):
            return field_discriminant, conductor
    raise ArithmeticError(f"no fundamental decomposition for {discriminant}")


def class_number_order(discriminant):
    maximum_a = isqrt(abs(discriminant) // 3) + 2
    count = 0
    for aa in range(1, maximum_a + 1):
        for bb in range(-aa, aa + 1):
            numerator = bb * bb - discriminant
            denominator = 4 * aa
            if numerator % denominator:
                continue
            cc = numerator // denominator
            if aa > cc or python_gcd(python_gcd(aa, abs(bb)), cc) != 1:
                continue
            if (abs(bb) == aa or aa == cc) and bb < 0:
                continue
            count += 1
    return count


def eichler_symbol(discriminant, prime):
    field_discriminant, conductor = fundamental_discriminant_and_conductor(
        discriminant
    )
    if conductor % prime == 0:
        return 1
    return int(kronecker_symbol(field_discriminant, prime))


def cm_screen(discriminant, level, maximum):
    if gcd(discriminant, level) != 1 or not is_squarefree(level):
        return {
            "status": "UNKNOWN_GENERAL_NONCOPRIME_OR_NONSQUAREFREE_LEVEL_FACTORS_NOT_IMPLEMENTED",
            "discriminant_bound": maximum,
            "small_embedding_orbits": [],
        }
    ramified = factor_primes(discriminant)
    level_primes = factor_primes(level)
    candidates = []
    for absolute in range(3, maximum + 1):
        delta = -absolute
        if delta % 4 not in (0, 1):
            continue
        field_discriminant, conductor = fundamental_discriminant_and_conductor(delta)
        if python_gcd(conductor, discriminant) != 1:
            continue
        if any(kronecker_symbol(field_discriminant, prime) == 1 for prime in ramified):
            continue
        class_number = class_number_order(delta)
        embedding_count = class_number
        for prime in ramified:
            embedding_count *= 1 - eichler_symbol(delta, prime)
        for prime in level_primes:
            embedding_count *= 1 + eichler_symbol(delta, prime)
        if not 0 < embedding_count <= 4:
            continue
        candidates.append(
            {
                "order_discriminant": delta,
                "field_discriminant": field_discriminant,
                "conductor": conductor,
                "ring_class_number": class_number,
                "optimal_embedding_classes_before_Atkin_Lehner": embedding_count,
                "rationality_status": "UNKNOWN_EMBEDDING_COUNT_IS_NOT_A_FIELD_OF_DEFINITION_CERTIFICATE",
            }
        )
    candidates.sort(
        key=lambda row: (
            row["optimal_embedding_classes_before_Atkin_Lehner"],
            abs(row["order_discriminant"]),
        )
    )
    return {
        "status": "PASS_EXACT_SQUAREFREE_EICHLER_OPTIMAL_EMBEDDING_SCREEN",
        "discriminant_bound": maximum,
        "small_embedding_orbits": candidates[:24],
        "truncated_to": 24,
    }


def shimura_genus(discriminant, level):
    assert gcd(discriminant, level) == 1
    ramified = factor_primes(discriminant)
    level_primes = factor_primes(level)
    mu = ZZ(discriminant * level)
    for prime in ramified:
        mu = mu * (prime - 1) / prime
    for prime in level_primes:
        mu = mu * (prime + 1) / prime
    assert mu in ZZ
    e2 = prod(1 - kronecker_symbol(-4, prime) for prime in ramified) * prod(
        1 + kronecker_symbol(-4, prime) for prime in level_primes
    )
    e3 = prod(1 - kronecker_symbol(-3, prime) for prime in ramified) * prod(
        1 + kronecker_symbol(-3, prime) for prime in level_primes
    )
    genus = 1 + QQ(mu) / 12 - QQ(e2) / 4 - QQ(e3) / 3
    assert genus in ZZ and genus >= 0
    return {
        "mu": int(mu),
        "elliptic_orbits_order_2": int(e2),
        "elliptic_orbits_order_3": int(e3),
        "cusp_count": 0,
        "genus": int(genus),
    }


def atkin_lehner_labels(discriminant, level):
    prime_powers = []
    for prime, exponent in factor(ZZ(discriminant * level)):
        prime_powers.append(int(prime**exponent))
    labels = [1]
    for prime_power in prime_powers:
        labels += [value * prime_power for value in tuple(labels)]
    return sorted(labels)


def h3_override(row, cm):
    assert row["surface_id"] == H3_SURFACE_ID
    assert row["literal_transcendental_gram"] == H3_T
    assert row["clifford"]["quaternion_discriminant"] == 6
    assert row["clifford"]["integral_even_clifford_order"]["local_level_index"] == 79
    assert row["rational_isotropy"]["status"] == "PASS_EXACT_ANISOTROPIC_LOCAL_OBSTRUCTION"
    anchors = [
        {
            "order_discriminant": -3,
            "curve_coordinates": "t=Infinity; in s=1/t chart (s,v)=(0,+/-4)",
            "atkin_lehner": "fixed by w_3",
            "singular_k3_transcendental_gram": [[2, 1], [1, 2]],
            "transported_primitive_k3_vector": [81, 95, -52],
        },
        {
            "order_discriminant": -24,
            "curve_coordinates": "(t,u)=(+/-2,+/-32)",
            "atkin_lehner": "one free four-point visible Klein orbit",
            "singular_k3_transcendental_gram": [[4, 0], [0, 6]],
            "transported_primitive_k3_vector": [70, 86, -3],
        },
    ]
    return {
        "status": "PASS_EXACT_H3_EICHLER_AND_ARITHMETIC_SOURCE_IDENTIFICATION",
        "order_type": {
            "status": "PASS_EXACT_EICHLER_ORDER_POSITIVE_CONTROL",
            "quaternion_discriminant": 6,
            "eichler_level": 79,
            "certificate_scripts": [
                "elkies-k3/scripts/construct_exact_gross_lattice.sage",
                "elkies-k3/scripts/map_clifford_to_k3_T.sage",
            ],
        },
        "base_curve": {
            "status": "PASS_EXACT_SHIMURA_CURVE",
            "label": "X_0^6(79)",
            **shimura_genus(6, 79),
        },
        "atkin_lehner": {
            "status": "PASS_EXACT_SELECTED_H3_QUOTIENT",
            "available_labels": [f"w_{value}" for value in atkin_lehner_labels(6, 79)],
            "selected_subgroup": ["1", "w_474"],
            "quotient_genus": 2,
            "quotient_model": "u^2 = 16*t^6 - 19*t^4 + 88*t^2 - 48",
            "further_kumar_quotient": "t-line (genus 0)",
            "certificate_script": "elkies-k3/scripts/deconstruct_x0679_quotients.sage",
        },
        "cm_screen": cm,
        "explicit_cm_anchors": {
            "status": "PASS_EXACT_TWO_H3_CM_ANCHORS",
            "anchors": anchors,
            "certificate_scripts": [
                "elkies-k3/scripts/transport_cm_delta3_to_k3.sage",
                "elkies-k3/scripts/deconstruct_x0679_quotients.sage",
            ],
        },
        "explicit_source_comparison": {
            "status": "PASS_EXACT_HISTORICAL_H3_SOURCE_MATCH",
            "non_cm_rational_point": "(t,u)=(14/13, 64*251/13^3)",
            "source": "Elkies X_0^6(79)/w_474 genus-two route",
        },
    }


def build_surface(surface, cm_maximum):
    literal = matrix(ZZ, surface["surface_key"]["transcendental_gram"])
    content, primitive, quadratic_scale, hessian = primitive_integral_hessian(literal)
    clifford = clifford_data(hessian)
    isotropy_result = pari(primitive).qfsolve()
    if isotropy_result.type() == "t_INT":
        rational_isotropy = {
            "status": "PASS_EXACT_ANISOTROPIC_LOCAL_OBSTRUCTION",
            "isotropic": False,
            "pari_qfsolve_obstruction_prime": int(isotropy_result),
            "primitive_isotropic_vector": None,
        }
        split = None
    else:
        isotropic = primitive_vector(isotropy_result)
        assert isotropic * primitive * isotropic == 0
        split = integral_u_split(hessian, isotropic)
        rational_isotropy = {
            "status": "PASS_EXACT_RATIONAL_ISOTROPIC_WITNESS",
            "isotropic": True,
            "pari_qfsolve_obstruction_prime": None,
            "primitive_isotropic_vector": list(map(int, isotropic)),
            "integral_u_split": split,
        }
    row = {
        "surface_id": surface["surface_id"],
        "determinant": int(surface["determinant"]),
        "literal_transcendental_gram": rows(literal),
        "similarity_normalization": {
            "literal_content": int(content),
            "primitive_bilinear_gram": rows(primitive),
            "quadratic_integrality_scale": int(quadratic_scale),
            "primitive_integral_quadratic_hessian": rows(hessian),
        },
        "rational_isotropy": rational_isotropy,
        "clifford": clifford,
    }
    algebra_discriminant = clifford["quaternion_discriminant"]
    level = clifford["integral_even_clifford_order"]["local_level_index"]
    if rational_isotropy["isotropic"]:
        assert algebra_discriminant == 1
        if split["status"] == "PASS_EXACT_INTEGRAL_U_SPLIT_EICHLER_ORDER":
            assert split["eichler_level"] == level
            arithmetic = {
                "status": "PASS_EXACT_SPLIT_EICHLER_MODULAR_CURVE",
                "order_type": {
                    "status": "PASS_EXACT_EICHLER_ORDER_FROM_INTEGRAL_U_SPLIT",
                    "quaternion_discriminant": 1,
                    "eichler_level": int(level),
                },
                "base_curve": modular_curve(level),
                "polarization_quotient": {
                    "status": "UNKNOWN_STABLE_ORTHOGONAL_DISCRIMINANT_KERNEL_NOT_COMPUTED",
                    "warning": "Gamma_0(N) is the norm-one Clifford-order curve; the final marked-NS quotient may be a finite cover or quotient.",
                },
            }
        else:
            arithmetic = {
                "status": "PARTIAL_EXACT_GENERAL_SPLIT_CLIFFORD_ORDER",
                "order_type": {
                    "status": "PASS_EXACT_SPLIT_QUATERNION_ORDER_EICHLER_LABEL_NOT_CERTIFIED",
                    "quaternion_discriminant": 1,
                    "reduced_discriminant_level_candidate": int(level),
                },
                "base_curve": {
                    "status": "UNKNOWN_GENERAL_CONGRUENCE_FUCHSIAN_SIGNATURE_NOT_COMPUTED",
                    "group": "projective norm-one group of the displayed split even-Clifford order",
                    "genus": None,
                    "cusps": None,
                },
                "hauptmodul": {
                    "status": "BLOCKED_UNTIL_GENERAL_SPLIT_ORDER_SIGNATURE_IS_COMPUTED",
                    "formula": None,
                },
                "explicit_modular_elliptic_surface_comparison": {
                    "status": "BLOCKED_UNTIL_CUSP_WIDTHS_AND_GENUS_ARE_COMPUTED"
                },
            }
    else:
        assert algebra_discriminant > 1
        coprime = gcd(algebra_discriminant, level) == 1
        if coprime:
            curve = {
                "status": "CONDITIONAL_X0DN_IF_CLIFFORD_ORDER_IS_EICHLER",
                "label": f"X_0^{algebra_discriminant}({level})",
                **shimura_genus(algebra_discriminant, level),
            }
            atkin = {
                "status": "CONDITIONAL_ATKIN_LEHNER_GROUP_IF_EICHLER",
                "available_labels": [
                    f"w_{value}"
                    for value in atkin_lehner_labels(algebra_discriminant, level)
                ],
                "selected_polarization_quotient": None,
                "quotient_genus": None,
            }
        else:
            curve = {
                "status": "UNKNOWN_NONMAXIMAL_AT_RAMIFIED_PRIME_NOT_AN_EICHLER_ORDER",
                "label": None,
                "genus": None,
                "cusp_count": 0,
            }
            atkin = {
                "status": "UNKNOWN_GENERAL_ORDER_NORMALIZER_NOT_COMPUTED",
                "available_labels": [],
                "selected_polarization_quotient": None,
                "quotient_genus": None,
            }
        cm = cm_screen(algebra_discriminant, level, cm_maximum)
        arithmetic = {
            "status": "PARTIAL_EXACT_ANISOTROPIC_CLIFFORD_ORDER",
            "order_type": {
                "status": (
                    "EICHLER_LEVEL_CANDIDATE_FROM_REDUCED_DISCRIMINANT_LOCAL_TYPE_NOT_CERTIFIED"
                    if coprime
                    else "PASS_EXACT_NON_EICHLER_AT_RAMIFIED_PRIME"
                ),
                "quaternion_discriminant": int(algebra_discriminant),
                "eichler_level_candidate": int(level) if coprime else None,
                "local_level_index": int(level),
            },
            "base_curve": curve,
            "atkin_lehner": atkin,
            "cm_screen": cm,
            "explicit_cm_anchors": {
                "status": "SEARCH_TRIGGERED_NO_EXPLICIT_CM_TRANSPORT_REGISTRY_HIT",
                "anchors": [],
                "order_discriminant_candidates": [
                    candidate["order_discriminant"]
                    for candidate in cm.get("small_embedding_orbits", [])[:8]
                ],
            },
        }
        if surface["surface_id"] == H3_SURFACE_ID:
            arithmetic = h3_override(row, cm)
    row["arithmetic_source"] = arithmetic
    identification_complete = arithmetic["status"].startswith("PASS_EXACT")
    row["pre_solver_gate"] = {
        "status": (
            "PASS_T_ARITHMETIC_CURVE_IDENTIFIED_SOLVER_GATE_OPEN"
            if identification_complete
            else "BLOCKED_T_ARITHMETIC_CURVE_IDENTIFICATION_INCOMPLETE"
        ),
        "arithmetic_attempt_recorded": True,
        "equation_solver_may_launch": identification_complete,
        "identification_complete": identification_complete,
        "policy": "No equation target may be emitted without a hash-matched, certified modular/Shimura curve row; typed arithmetic unknowns are retained, never imputed, and keep the solver gate closed.",
    }
    return row


def build(catalogue, catalogue_path, cm_maximum):
    assert catalogue["schema"] == "elkies-k3.rank7-auxiliary-catalogue.v1"
    surfaces = [build_surface(surface, cm_maximum) for surface in catalogue["surfaces"]]
    surfaces.sort(key=lambda row: row["surface_id"])
    isotropic = [row for row in surfaces if row["rational_isotropy"]["isotropic"]]
    anisotropic = [row for row in surfaces if not row["rational_isotropy"]["isotropic"]]
    gamma0 = [
        row
        for row in isotropic
        if row["arithmetic_source"]["status"]
        == "PASS_EXACT_SPLIT_EICHLER_MODULAR_CURVE"
    ]
    genus_zero = [
        row
        for row in gamma0
        if row["arithmetic_source"]["base_curve"]["genus"] == 0
    ]
    hauptmodul = [
        row
        for row in genus_zero
        if row["arithmetic_source"]["base_curve"]["hauptmodul"]["status"]
        == "PASS_EXACT_REGISTRY_HAUPTMODUL"
    ]
    candidate_eichler = [
        row
        for row in anisotropic
        if row["arithmetic_source"]["order_type"]["status"]
        == "EICHLER_LEVEL_CANDIDATE_FROM_REDUCED_DISCRIMINANT_LOCAL_TYPE_NOT_CERTIFIED"
    ]
    non_eichler = [
        row
        for row in anisotropic
        if row["arithmetic_source"]["order_type"]["status"]
        == "PASS_EXACT_NON_EICHLER_AT_RAMIFIED_PRIME"
    ]
    h3 = next(row for row in surfaces if row["surface_id"] == H3_SURFACE_ID)
    assert len(surfaces) == 827
    assert len(isotropic) == 550
    assert len(anisotropic) == 277
    assert len(gamma0) == 313
    assert len(candidate_eichler) == 171
    assert len(non_eichler) == 105
    assert sum(
        row["pre_solver_gate"]["equation_solver_may_launch"] for row in surfaces
    ) == 314
    assert h3["arithmetic_source"]["status"] == (
        "PASS_EXACT_H3_EICHLER_AND_ARITHMETIC_SOURCE_IDENTIFICATION"
    )
    return {
        "schema": "elkies-k3.rank7-t-arithmetic.v1",
        "status": "PASS_EXACT_PRE_SOLVER_T_ARITHMETIC_LEDGER_WITH_TYPED_OPEN_CURVE_IDENTIFICATIONS",
        "input": {
            "catalogue": relative(catalogue_path),
            "catalogue_sha256": digest(catalogue_path),
        },
        "parameters": {
            "cm_order_discriminant_absolute_bound": cm_maximum,
        },
        "proof_scope": {
            "proved": (
                "Every imported (T,NS) has an exact primitive-similarity normalization, "
                "rational isotropy decision, even Clifford quaternion algebra, and integral "
                "Clifford-order reduced discriminant. Divisibility-one isotropic witnesses "
                "give exact integral U splittings and Gamma_0(N) cusp/genus data. The H3 "
                "positive control retains its separately certified X_0^6(79) arithmetic."
            ),
            "not_proved": (
                "A reduced-discriminant quotient is not by itself an Eichler local-type "
                "certificate. General split-order signatures, general non-Eichler Shimura "
                "curves, polarization discriminant-kernel quotients, Atkin--Lehner fixed-point "
                "genera, rationality of CM orbits, and new explicit CM/K3 anchors remain typed "
                "unknown unless a dedicated certificate is attached."
            ),
        },
        "accounting": {
            "surfaces": len(surfaces),
            "rationally_isotropic": len(isotropic),
            "anisotropic": len(anisotropic),
            "certified_integral_U_split_Gamma0": len(gamma0),
            "general_split_orders_signature_open": len(isotropic) - len(gamma0),
            "certified_Gamma0_genus_zero": len(genus_zero),
            "registry_Hauptmodul_hits": len(hauptmodul),
            "anisotropic_H3_certified_Eichler": 1,
            "anisotropic_Eichler_level_candidates_local_type_open": len(candidate_eichler),
            "anisotropic_non_Eichler_at_ramified_prime": len(non_eichler),
            "pre_solver_arithmetic_attempt_rows": len(surfaces),
            "pre_solver_gate_open_rows": sum(
                row["pre_solver_gate"]["equation_solver_may_launch"]
                for row in surfaces
            ),
            "isotropy_obstruction_prime_distribution": dict(
                sorted(
                    Counter(
                        str(row["rational_isotropy"]["pari_qfsolve_obstruction_prime"])
                        for row in anisotropic
                    ).items()
                )
            ),
        },
        "surfaces": surfaces,
        "literature_boundary": [
            {
                "citation": "Elkies, Shimura curve computations via K3 surfaces of Neron--Severi rank at least 19",
                "arxiv": "0802.1301",
                "use": "arithmetic source identification and CM rank-20 anchors before equation construction",
            },
            {
                "citation": "Voight, Quaternion Algebras",
                "url": "https://jvoight.github.io/quat-book.pdf",
                "use": "ternary quadratic forms, Clifford orders, Eichler orders, and arithmetic Fuchsian groups",
            },
            {
                "citation": "Padurariu--Saia, Shimura curve Atkin--Lehner quotients of genus at most two",
                "arxiv": "2509.25368",
                "use": "independent low-genus quotient equations for the H3 positive control",
            },
        ],
        "reproduce": (
            "sage -python elkies-k3/scripts/build_rank7_t_arithmetic.sage "
            f"--cm-disc-max {cm_maximum}"
        ),
    }


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--catalogue", type=Path, default=DEFAULT_CATALOGUE)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--cm-disc-max", type=int, default=200)
parser.add_argument("--check", action="store_true")
arguments = parser.parse_args()

catalogue_path = arguments.catalogue.resolve()
payload = build(
    json.loads(catalogue_path.read_text()),
    catalogue_path,
    arguments.cm_disc_max,
)
encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
output = arguments.output.resolve()
if arguments.check:
    if not output.exists() or output.read_text() != encoded:
        raise SystemExit("rank-seven T-arithmetic artifact is stale")
else:
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(encoded)

accounting = payload["accounting"]
print(
    "RANK7TARITH|surfaces={}|isotropic={}|anisotropic={}|gamma0={}|"
    "split_open={}|h3=1|status=PASS".format(
        accounting["surfaces"],
        accounting["rationally_isotropic"],
        accounting["anisotropic"],
        accounting["certified_integral_U_split_Gamma0"],
        accounting["general_split_orders_signature_open"],
    )
)
