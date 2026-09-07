#!/usr/bin/env sage-python
"""Certify hidden rational-surface base changes in repeated-fibre sources.

The audit has two exact promoted-model targets:

* the rational Golay-chart ``3I6+6I1`` K3, and
* the marked NS0031 ``I2+I8+I8+6I1`` model over ``GF(7)``.

It first performs a complete j-stabilizer census on the exhaustive normalized
fibre charts over GF(5) and GF(7).  It then proves coefficient-level quadratic
descent for the two promoted models, reconstructs the rational elliptic
quotients, and checks the induced deck action on the available sections.
Repeated fibre multiplicities alone are never used as a descent certificate.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing, QQ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
GEN = ROOT / "artifacts/generated-results"
DEFAULT_OUTPUT = GEN / "elkies-k3-repeated-fibre-rational-base-change-audit-v1.json"

GOLAY_FIBRES = {
    5: GEN / "elkies-k3-golay-det720-3a5-source-ansatz-mod5-v1.json",
    7: GEN / "elkies-k3-golay-det720-3a5-source-ansatz-mod7-v1.json",
}
A1_2A7_FIBRES = {
    5: GEN / "elkies-k3-lattice-foundry-ns0005-a1-2a7-source-ansatz-mod5-v1.json",
    7: GEN / "elkies-k3-lattice-foundry-ns0005-a1-2a7-source-ansatz-mod7-v1.json",
}
GOLAY_PAIRS_7 = GEN / "elkies-k3-golay-det720-3a5-pole0-pairs-mod7-nonsquare-v1.json"
NS0031_MARKING_7 = GEN / "elkies-k3-lattice-foundry-ns0031-a1-2a7-marking-mod7-v1.json"
GOLAY_QQ = GEN / "elkies-k3-golay-det720-3a5-source-qq-v1.json"
GOLAY_SATURATION = GEN / "elkies-k3-golay-det720-3a5-saturation-rejection-v1.json"
GOLAY_RATIONAL_SCANS = sorted(GEN.glob("elkies-k3-golay-det720-3a5-rational-parameter-scan*-v1.json"))


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def serialize_polynomial(poly):
    return [str(value) for value in poly.list()]


def j_map(A, B):
    field = A.parent().fraction_field()
    return field(4 * A**3) / field(4 * A**3 + 27 * B**2)


def support_maps(t):
    """The five nonidentity maps permuting the normalized set {0,1,infinity}."""

    return {
        "1-t": 1 - t,
        "1/t": 1 / t,
        "1/(1-t)": 1 / (1 - t),
        "(t-1)/t": (t - 1) / t,
        "t/(t-1)": t / (t - 1),
    }


def symmetry_type(symmetries):
    symmetries = tuple(symmetries)
    if not symmetries:
        return "trivial"
    if len(symmetries) == 1:
        return f"C2:{symmetries[0]}"
    if set(symmetries) == {"1/(1-t)", "(t-1)/t"}:
        return "C3"
    if len(symmetries) == 5:
        return "S3"
    raise ArithmeticError(f"unexpected normalized support stabilizer {symmetries}")


def model_j_symmetries(example, prime, candidates):
    field = GF(prime)
    ring = PolynomialRing(field, "t")
    t = ring.gen()
    function_field = ring.fraction_field()
    A = ring(example["A_coefficients_low_to_high"])
    B = ring(example["B_coefficients_low_to_high"])
    j = j_map(A, B)
    return [
        name
        for name, image in candidates(t).items()
        if j == j(function_field(image))
    ]


def golay_candidates(t):
    return support_maps(t)


def a1_2a7_candidates(t):
    # Fibre types force 0:I2 to be fixed and allow only 1:I8 <-> infinity:I8.
    return {"t/(t-1)": t / (t - 1)}


def census(path, candidates):
    payload = json.loads(path.read_text())
    prime = int(payload["prime"])
    histogram = Counter()
    symmetric_indices = []
    for index, example in enumerate(payload["examples"]):
        symmetries = model_j_symmetries(example, prime, candidates)
        kind = symmetry_type(symmetries)
        histogram[kind] += 1
        if symmetries:
            symmetric_indices.append(index)
    return {
        "input": relative(path),
        "prime": prime,
        "model_count": len(payload["examples"]),
        "nontrivial_j_stabilizer_count": len(symmetric_indices),
        "stabilizer_histogram": dict(sorted(histogram.items())),
        "symmetric_example_indices": symmetric_indices,
    }


def fit_descended_polynomial(target, quotient_coordinate, prefactor, degree):
    """Solve target=prefactor*bar(quotient_coordinate) coefficientwise."""

    ring = target.parent()
    field = ring.base_ring()
    basis = [ring(prefactor * quotient_coordinate**index) for index in range(degree + 1)]
    row_count = max(poly.degree() for poly in basis) + 1
    system = matrix(
        field,
        [[poly[row] for poly in basis] for row in range(row_count)],
    )
    coefficients = system.solve_right(vector(field, [target[row] for row in range(row_count)]))
    reconstructed = sum(
        (coefficients[index] * basis[index] for index in range(degree + 1)),
        ring.zero(),
    )
    if reconstructed != target:
        raise ArithmeticError("coefficient-level quotient reconstruction failed")
    quotient_ring = PolynomialRing(field, "u")
    return quotient_ring(list(coefficients))


def quotient_record(A, B, quotient_coordinate, a_prefactor, b_prefactor, fibre_profile):
    a = fit_descended_polynomial(A, quotient_coordinate, a_prefactor, 4)
    b = fit_descended_polynomial(B, quotient_coordinate, b_prefactor, 6)
    discriminant = 4 * a**3 + 27 * b**2
    numerator = 4 * a**3
    if numerator.gcd(discriminant) != 1:
        raise ArithmeticError("quotient j-map is not reduced")
    return {
        "equation": "Y^2=X^3+a(u)X+b(u)",
        "a_coefficients_low_to_high": serialize_polynomial(a),
        "b_coefficients_low_to_high": serialize_polynomial(b),
        "a_degree": int(a.degree()),
        "b_degree": int(b.degree()),
        "discriminant_factorization": str(discriminant.factor()),
        "discriminant_degree": int(discriminant.degree()),
        "j_map_degree": int(max(numerator.degree(), discriminant.degree())),
        "fibre_profile": fibre_profile,
        "euler_number": 12,
        "rational_elliptic_surface_degree_bounds": bool(a.degree() <= 4 and b.degree() <= 6),
    }


def exact_golay_quotient(model, saturation):
    ring = PolynomialRing(QQ, "t")
    t = ring.gen()
    function_field = ring.fraction_field()
    A = ring(model["weierstrass_model"]["A_coefficients_low_to_high"])
    B = ring(model["weierstrass_model"]["B_coefficients_low_to_high"])
    discriminant = 4 * A**3 + 27 * B**2
    involution = 1 / t
    j = j_map(A, B)
    candidates = support_maps(t)
    symmetries = [name for name, image in candidates.items() if j == j(function_field(image))]
    if symmetries != ["1/t"]:
        raise ArithmeticError("canonical rational 3I6 j-stabilizer changed")
    if function_field(A(involution)) * t**8 != A:
        raise ArithmeticError("A does not have the required reciprocal covariance")
    if function_field(B(involution)) * t**12 != B:
        raise ArithmeticError("B does not have the required reciprocal covariance")
    if function_field(discriminant(involution)) * t**24 != discriminant:
        raise ArithmeticError("the full discriminant support is not reciprocal")

    u = t + 1 / t
    quotient = quotient_record(A, B, u, t**4, t**6, "I6+I3+3I1")
    quotient_ring = PolynomialRing(QQ, "u")
    u_bar = quotient_ring.gen()
    a = quotient_ring(quotient["a_coefficients_low_to_high"])
    b = quotient_ring(quotient["b_coefficients_low_to_high"])
    d_bar = 4 * a**3 + 27 * b**2
    if d_bar.valuation(u_bar - 2) != 3 or d_bar(-2) == 0:
        raise ArithmeticError("the quotient branch-fibre profile changed")

    curve = EllipticCurve(function_field, [0, 0, 0, function_field(A), function_field(B)])
    P_record = model["marked_sections"][0]
    P = curve(
        function_field(ring(P_record["X_coefficients_low_to_high"])),
        function_field(ring(P_record["Y_coefficients_low_to_high"])),
    )
    torsion_record = saturation["torsion_section"]
    T = curve(
        function_field(ring(torsion_record["X_coefficients_low_to_high"])),
        function_field(ring(torsion_record["Y_coefficients_low_to_high"])),
    )
    half_record = saturation["half_of_displayed_Q"]
    H = curve(
        function_field(ring(half_record["X_coefficients_low_to_high"])),
        function_field(ring(half_record["Y_coefficients_low_to_high"])),
    )

    def deck(point):
        return curve(
            function_field(t**4 * point[0](involution)),
            function_field(t**6 * point[1](involution)),
        )

    if deck(T) != T or deck(H) != H:
        raise ArithmeticError("known quotient sections are not deck-invariant")
    if P + deck(P) != 2 * T:
        raise ArithmeticError("trace of P under the deck involution changed")
    if deck(P - T) != -(P - T):
        raise ArithmeticError("P-T is not anti-invariant")

    # Explicit quotient coordinates of the invariant torsion and free sections.
    Hx = 3 * u_bar**2 + QQ(36) / 5 * u_bar + QQ(2988) / 25
    Hy = quotient_ring(QQ(41472) / 25)
    Tx = 3 * u_bar**2 + QQ(36) / 5 * u_bar + QQ(108) / 25
    Ty = QQ(1728) / 5 * (u_bar - 2)
    quotient_curve = EllipticCurve(
        quotient_ring.fraction_field(),
        [0, 0, 0, quotient_ring.fraction_field()(a), quotient_ring.fraction_field()(b)],
    )
    H_bar = quotient_curve(Hx, Hy)
    T_bar = quotient_curve(Tx, Ty)
    if 3 * T_bar != quotient_curve(0):
        raise ArithmeticError("descended torsion section does not have order three")
    if 2 * H_bar == quotient_curve(0):
        raise ArithmeticError("descended free section unexpectedly became 2-torsion")

    quotient["free_mordell_weil_rank_from_shioda_tate"] = 1
    quotient["torsion_section_exact_order"] = 3
    quotient["free_generator_height_from_degree_two_pullback"] = "1/2"
    return {
        "source": "rational s6=10 Golay-chart model",
        "base_field": "QQ",
        "source_fibre_profile": "3I6+6I1",
        "weighted_reducible_support_automorphism_group": "S3",
        "full_j_and_discriminant_stabilizer": ["identity", "1/t"],
        "cubic_support_maps_fail_j_invariance": ["1/(1-t)", "(t-1)/t"],
        "deck_involution": "t -> 1/t",
        "quotient_coordinate": "u=t+1/t",
        "branch_points_on_source": ["t=1", "t=-1"],
        "branch_values_on_quotient": ["u=2", "u=-2"],
        "weierstrass_scaling": "x=t^2*X, y=t^3*Y",
        "j_map_decomposition": "j_K3(t)=j_RES(t+1/t)",
        "source_j_map_degree": 24,
        "quotient": quotient,
        "sections": {
            "invariant_torsion": "T, exact order 3",
            "invariant_free_section": "H with 2H=Q",
            "anti_invariant_free_section": "P-T",
            "deck_trace_identity": "P+sigma(P)=2T",
            "rank_explanation": (
                "The free rank splits as invariant rank 1 from the rational surface "
                "plus anti-invariant rank 1 from its quadratic twist."
            ),
            "descended_H": {
                "x": str(Hx),
                "y": str(Hy),
                "height": "1/2",
            },
            "descended_T": {
                "x": str(Tx),
                "y": str(Ty),
                "exact_order": 3,
            },
        },
    }


def ns0031_quotient(fibres, marking):
    prime = 7
    field = GF(prime)
    ring = PolynomialRing(field, "t")
    t = ring.gen()
    function_field = ring.fraction_field()
    example_index = 157
    example = fibres["examples"][example_index]
    marked_model = next(row for row in marking["models"] if row["example_index"] == example_index)
    if marking["accounting"]["marked_mw2_pairs"] != 2:
        raise ArithmeticError("NS0031 marked-pair count changed")
    A = ring(example["A_coefficients_low_to_high"])
    B = ring(example["B_coefficients_low_to_high"])
    discriminant = 4 * A**3 + 27 * B**2
    involution = function_field(t) / (t - 1)
    j = j_map(A, B)
    if j != j(involution):
        raise ArithmeticError("NS0031 j-map lost its involution")
    if function_field(A(involution)) * (t - 1) ** 8 != A:
        raise ArithmeticError("NS0031 A covariance failed")
    if function_field(B(involution)) * (t - 1) ** 12 != B:
        raise ArithmeticError("NS0031 B covariance failed")
    if function_field(discriminant(involution)) * (t - 1) ** 24 != discriminant:
        raise ArithmeticError("NS0031 full discriminant support is not invariant")

    u = t**2 / (t - 1)
    quotient = quotient_record(A, B, u, (t - 1) ** 4, (t - 1) ** 6, "I8+4I1")
    quotient_ring = PolynomialRing(field, "u")
    u_bar = quotient_ring.gen()
    a = quotient_ring(quotient["a_coefficients_low_to_high"])
    b = quotient_ring(quotient["b_coefficients_low_to_high"])
    d_bar = 4 * a**3 + 27 * b**2
    if d_bar.valuation(u_bar) != 1 or d_bar(4) == 0:
        raise ArithmeticError("NS0031 quotient branch-fibre profile changed")
    quotient["free_mordell_weil_rank_from_shioda_tate"] = 1

    zero_record = marked_model["pole_zero_sections"][1]
    one_record = marked_model["pole_one_sections"][0]
    if one_record["marked_pairs"][0]["pole_zero_index"] != 1:
        raise ArithmeticError("selected NS0031 marked pair changed")
    curve = EllipticCurve(function_field, [0, 0, 0, function_field(A), function_field(B)])
    P = curve(
        function_field(ring(zero_record["X_coefficients_low_to_high"])),
        function_field(ring(zero_record["Y_coefficients_low_to_high"])),
    )
    denominator = ring(one_record["C_coefficients_low_to_high"])
    Q = curve(
        function_field(ring(one_record["X_numerator_coefficients_low_to_high"])) / denominator**2,
        function_field(ring(one_record["Y_numerator_coefficients_low_to_high"])) / denominator**3,
    )

    def deck(point):
        return curve(
            function_field((t - 1) ** 4 * point[0](involution)),
            function_field((t - 1) ** 6 * point[1](involution)),
        )

    anti = 3 * P - 2 * Q
    trace = P + deck(P)
    if deck(anti) != -anti:
        raise ArithmeticError("the NS0031 marked pair lost its anti-invariant direction")
    if trace == curve(0) or deck(trace) != trace:
        raise ArithmeticError("the NS0031 invariant trace section failed")

    return {
        "source": "NS0031 model 157 with a complete marked MW2 pair",
        "base_field": "GF(7)",
        "source_fibre_profile": "I2+I8+I8+6I1",
        "weighted_reducible_support_automorphism_group": "C2",
        "full_j_and_discriminant_stabilizer": ["identity", "t/(t-1)"],
        "cubic_support_maps": "none compatible with the weighted reducible supports",
        "deck_involution": "t -> t/(t-1)",
        "quotient_coordinate": "u=t^2/(t-1)",
        "branch_points_on_source": ["t=0", "t=2"],
        "branch_values_on_quotient": ["u=0", "u=4"],
        "weierstrass_scaling": "x=(t-1)^2*X, y=(t-1)^3*Y",
        "j_map_decomposition": "j_K3(t)=j_RES(t^2/(t-1))",
        "source_j_map_degree": 24,
        "quotient": quotient,
        "sections": {
            "selected_marked_basis": ["P=pole-zero section 1", "Q=pole-one section 0"],
            "invariant_nonzero_section": "P+sigma(P)",
            "anti_invariant_nonzero_section": "3P-2Q",
            "rank_explanation": (
                "The marked pair is compatible with one rational-surface direction "
                "and one quadratic-twist direction, but the full GF(7)(t) MW group "
                "and any characteristic-zero lift are not certified here."
            ),
        },
    }


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--check", action="store_true")
arguments = parser.parse_args()
output_path = arguments.output.resolve()

input_paths = (
    list(GOLAY_FIBRES.values())
    + list(A1_2A7_FIBRES.values())
    + [GOLAY_PAIRS_7, NS0031_MARKING_7, GOLAY_QQ, GOLAY_SATURATION]
    + GOLAY_RATIONAL_SCANS
)
for path in input_paths:
    if not path.exists():
        raise FileNotFoundError(path)

golay_censuses = [census(GOLAY_FIBRES[prime], golay_candidates) for prime in (5, 7)]
a1_2a7_censuses = [census(A1_2A7_FIBRES[prime], a1_2a7_candidates) for prime in (5, 7)]
expected_golay = {
    5: {
        "C2:1-t": 8,
        "C2:1/t": 8,
        "C2:t/(t-1)": 8,
        "C3": 4,
        "S3": 3,
        "trivial": 42,
    },
    7: {
        "C2:1-t": 22,
        "C2:1/t": 22,
        "C2:t/(t-1)": 22,
        "C3": 16,
        "S3": 5,
        "trivial": 150,
    },
}
for row in golay_censuses:
    if row["stabilizer_histogram"] != expected_golay[row["prime"]]:
        raise ArithmeticError("Golay j-stabilizer census changed")
expected_a1_2a7 = {
    5: {"C2:t/(t-1)": 15, "trivial": 56},
    7: {"C2:t/(t-1)": 43, "trivial": 228},
}
for row in a1_2a7_censuses:
    if row["stabilizer_histogram"] != expected_a1_2a7[row["prime"]]:
        raise ArithmeticError("A1+2A7 j-stabilizer census changed")

golay_fibres_7 = json.loads(GOLAY_FIBRES[7].read_text())
golay_pairs_7 = json.loads(GOLAY_PAIRS_7.read_text())
promoted_golay_models = []
for row in golay_pairs_7["models"]:
    symmetries = model_j_symmetries(
        golay_fibres_7["examples"][row["example_index"]], 7, golay_candidates
    )
    promoted_golay_models.append(
        {
            "example_index": int(row["example_index"]),
            "marked_mw2_pairs": len(row["marked_mw2_pairs"]),
            "j_stabilizer": ["identity"] + symmetries,
            "stabilizer_type": symmetry_type(symmetries),
        }
    )
if sum(row["marked_mw2_pairs"] for row in promoted_golay_models) != 24:
    raise ArithmeticError("promoted Golay pair accounting changed")
if any(
    row["stabilizer_type"].startswith("C3") or row["stabilizer_type"] == "S3"
    for row in promoted_golay_models
    if row["marked_mw2_pairs"]
):
    raise ArithmeticError("a promoted Golay marked pair unexpectedly has cubic symmetry")

rational_models = []
for scan_path in GOLAY_RATIONAL_SCANS:
    scan = json.loads(scan_path.read_text())
    for row in scan.get("exact_rational_points", []):
        ring = PolynomialRing(QQ, "t")
        t = ring.gen()
        function_field = ring.fraction_field()
        coordinates = row["coordinates"]
        A = ring(coordinates[0:9])
        B = ring(coordinates[9:22])
        j = j_map(A, B)
        symmetries = [
            name
            for name, image in support_maps(t).items()
            if j == j(function_field(image))
        ]
        rational_models.append(
            {
                "input": relative(scan_path),
                "parameter": row["parameter"],
                "j_stabilizer": ["identity"] + symmetries,
                "stabilizer_type": symmetry_type(symmetries),
            }
        )
if len(rational_models) != 3 or any(row["stabilizer_type"].startswith("C3") for row in rational_models):
    raise ArithmeticError("bounded rational Golay base-change census changed")

golay_model = json.loads(GOLAY_QQ.read_text())
golay_saturation = json.loads(GOLAY_SATURATION.read_text())
ns0031_fibres = json.loads(A1_2A7_FIBRES[7].read_text())
ns0031_marking = json.loads(NS0031_MARKING_7.read_text())
exact_golay = exact_golay_quotient(golay_model, golay_saturation)
exact_ns0031 = ns0031_quotient(ns0031_fibres, ns0031_marking)

payload = {
    "schema": "elkies-k3.repeated-fibre-rational-base-change-audit.v1",
    "status": "PASS_EXACT_QUADRATIC_RATIONAL_SURFACE_DESCENT_FOR_PROMOTED_REPEATED_FIBRE_MODELS",
    "reproduce": (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/certify_repeated_fibre_rational_base_changes.sage"
    ),
    "inputs": {relative(path): digest(path) for path in input_paths},
    "finite_field_j_stabilizer_censuses": {
        "3I6": golay_censuses,
        "I2+2I8": a1_2a7_censuses,
    },
    "promoted_finite_field_models": {
        "Golay_3I6_GF7_models_with_both_section_types": promoted_golay_models,
        "NS0031_I2_2I8_GF7": {
            "example_index": 157,
            "marked_mw2_pairs": 2,
            "j_stabilizer": ["identity", "t/(t-1)"],
            "stabilizer_type": "C2:t/(t-1)",
        },
    },
    "bounded_exact_rational_Golay_models": rational_models,
    "exact_promoted_model_certificates": [exact_golay, exact_ns0031],
    "conclusion": {
        "quadratic_base_change_detected": True,
        "cubic_base_change_detected_for_a_promoted_marked_source": False,
        "interpretation": (
            "Both promoted repeated-fibre patterns arise from quadratic pullback of "
            "a rational elliptic surface.  The rational 3I6 model descends to "
            "I6+I3+3I1, and the marked NS0031 model over GF(7) descends to I8+4I1."
        ),
    },
    "proof_boundary": {
        "proved": (
            "All displayed j-invariance, full-discriminant covariance, weighted "
            "Weierstrass descent, quotient equations, quotient fibre profiles, and "
            "section deck identities are exact over the stated fields.  The finite-field "
            "censuses exhaust the stored normalized fibre charts."
        ),
        "not_proved": (
            "The NS0031 result is over GF(7), not a characteristic-zero or rational "
            "source construction.  The audit does not prove that every model in either "
            "fibre stratum is a base change, that the finite-field MW groups are saturated, "
            "or that a cubic construction is impossible outside the normalized weighted "
            "support automorphism groups.  The rational 3I6 quotient concerns the already "
            "rejected determinant-20 specialization, not the determinant-720 target K3."
        ),
    },
}

encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
if arguments.check:
    if not output_path.exists() or output_path.read_text() != encoded:
        raise SystemExit(f"stale artifact: {output_path}")
else:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(encoded)

print(
    "REPEATEDFIBREBASECHANGE|"
    "G720=quadratic_RES_I6_I3_3I1|"
    "NS0031_mod7=quadratic_RES_I8_4I1|"
    "promoted_cubic=none|status=PASS",
    flush=True,
)
print(f"OUTPUT|{output_path}", flush=True)
