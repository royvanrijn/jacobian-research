#!/usr/bin/env sage -python
"""Exhaust polynomial sections on the exact q8 child modulo 103.

status: EXPERIMENT
claim: bounded finite-field enumeration split into four-variable blocks
input: exact A11-to-2A5 q8 equation and equation marking
output: artifacts/local/elkies-k3/q24-2a5-zero-pole-sections-p103.json

For X=x0+...+x4*T^4 and Y=y0+...+y6*T^6, the leading pair
(x4,y6) lies on the fibre at infinity.  At p=103 that fibre has 96 affine
points and none has y6=0.  Fixing one leading point lets coefficients 11 down
to 6 of the Weierstrass identity determine y5,...,y0 recursively.  Each
remaining calculation is therefore only six equations in x3,...,x0.

The modular points are construction aids.  Component incidence and pair
intersection gates do not identify or lift a characteristic-zero section.
"""

import hashlib
import json
import time
from pathlib import Path

from sage.all import GF, PolynomialRing, PowerSeriesRing


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
EXACT_EQUATION = LOCAL / "q24-a11-to-2a5-q8-resolved-rr-qq.json"
EXACT_MARKING = LOCAL / "q24-a11-to-2a5-q8-equation-marking-qq.json"
WORD = LOCAL / "q24-2a5-q6o1307-horizontal-word.json"
OUTPUT = LOCAL / "q24-2a5-zero-pole-sections-p103.json"

PRIME = 103
F = GF(PRIME)

# Reducing the 76 MB exact certificate is slower than this bounded search.
# These coefficients are pinned to the exact input hashes checked below.
A_COEFFICIENTS = [34, 39, 41, 40, 86, 72, 82, 7, 61]
B_COEFFICIENTS = [55, 85, 90, 96, 0, 80, 84, 24, 24, 42, 6, 67, 63]
AFFINE_X_COEFFICIENTS = [28, 36, 58, 25, 72]
AFFINE_Y_COEFFICIENTS = [89, 6, 89, 46, 43, 50, 65]
EXPECTED_HASHES = {
    EXACT_EQUATION: "9e7b3c2fab5a1b8d95820f4e9bcde127908588c3efbe65678669efcd49d2ac3b",
    EXACT_MARKING: "b9d46929f37de95c68910f6246ee07b32d7bdb9476b79e6629fda0cc24cb23f9",
    WORD: "7103f196450fa7ef691927259ed56739f68829f4a78ff93eff9e87584ca6f3ec",
}


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


for path, expected in EXPECTED_HASHES.items():
    actual = sha256(path)
    if actual != expected:
        raise RuntimeError(f"input hash changed for {path}: {actual} != {expected}")

RT = PolynomialRing(F, "T")
T = RT.gen()
A = RT(A_COEFFICIENTS)
B = RT(B_COEFFICIENTS)
Delta = -F(16) * (4 * A**3 + 27 * B**2)

RX = PolynomialRing(F, "x")
x = RX.gen()
i6_fibres = []
for factor, exponent in Delta.factor():
    if int(exponent) != 6:
        continue
    assert factor.degree() == 1
    root = F(-factor[0] / factor[1])
    cubic = x**3 + A(root) * x + B(root)
    repeated = cubic.gcd(cubic.derivative())
    assert repeated.degree() == 1
    node_x = F(-repeated[0] / repeated[1])
    i6_fibres.append((root, node_x))
i6_fibres.sort(key=lambda item: int(item[0]))
assert [(int(root), int(node)) for root, node in i6_fibres] == [(68, 90), (89, 65)]

affine = (RT(AFFINE_X_COEFFICIENTS), RT(AFFINE_Y_COEFFICIENTS))
assert affine[1]**2 == affine[0]**3 + A * affine[0] + B

# Fix the split-I6 orientations by the exact old-A11 affine section.  It meets
# old component 1 over T=68 and old component 0 over T=89; both are oriented
# component 4 in the Tate six-cycle.  This pins more than node/depth data.
I6_ORIENTATION = {
    68: {
        "rho_constant": 95,
        "labels_by_oriented_component": {
            0: "old_A11_component_8",
            1: "old_A11_component_7",
            2: "old_A11_component_6",
            3: "old_A11_component_2",
            4: "old_A11_component_1",
            5: "second_I6_affine_component",
        },
    },
    89: {
        "rho_constant": 35,
        "labels_by_oriented_component": {
            0: "old_A11_component_10",
            1: "old_A11_component_5",
            2: "old_A11_component_4",
            3: "old_A11_component_3",
            4: "old_A11_component_0",
            5: "first_I6_affine_component",
        },
    },
}


def shifted_series(poly, root, series_ring):
    s = series_ring.gen()
    return sum(series_ring(poly[index]) * (F(root) + s)**index for index in range(poly.degree() + 1))


def oriented_i6_component(section, root, node_x):
    section_x, section_y = section
    if section_x(F(root)) != F(node_x) or section_y(F(root)) != 0:
        return 0
    series_ring = PowerSeriesRing(F, "s", default_prec=10)
    A_local = shifted_series(A, root, series_ring)
    B_local = shifted_series(B, root, series_ring)
    X_local = shifted_series(section_x, root, series_ring)
    Y_local = shifted_series(section_y, root, series_ring)
    center = series_ring(F(node_x))
    for unused in range(5):
        center = (center + (-A_local / 3) / center) / 2
    assert (center**3 + A_local * center + B_local).valuation() == 6
    rho = series_ring(F(I6_ORIENTATION[int(root)]["rho_constant"]))
    for unused in range(5):
        rho = (rho + (X_local + 2 * center) / rho) / 2
    first = int((Y_local - rho * (X_local - center)).valuation())
    second = int((Y_local + rho * (X_local - center)).valuation())
    assert first + second == 6 and 1 <= first <= 5
    return first


assert [
    oriented_i6_component(affine, root, node) for root, node in i6_fibres
] == [4, 4]


def padded_coefficients(poly, length):
    return [int(poly[index]) if index <= poly.degree() else 0 for index in range(length)]


def resolved_pair_intersection(left, right):
    """Compute P.Q=(P-Q).O from the cancelled group-law denominator."""
    x_left, y_left = left
    x_right, y_right = right
    difference = x_left - x_right
    if difference == 0:
        return None
    numerator = (y_left + y_right)**2 - difference**2 * (x_left + x_right)
    first = difference.gcd(numerator)
    numerator_once = numerator // first
    cancellation = difference.gcd(numerator_once)
    reduced_denominator = difference // cancellation
    reduced_numerator = numerator // cancellation**2
    finite = max(0, int(reduced_denominator.degree()))
    excess = int(reduced_numerator.degree()) - 2 * finite - 4
    if excess > 0:
        assert excess % 2 == 0
    return finite + (excess // 2 if excess > 0 else 0)


leading_points = []
for x4 in F:
    infinity_rhs = x4**3 + A[8] * x4 + B[12]
    for y6 in F:
        if y6**2 == infinity_rhs:
            leading_points.append((x4, y6))
assert len(leading_points) == 96
assert all(y6 for _, y6 in leading_points)

started = time.monotonic()
sections = {}
blocks = []
for block_index, (x4, y6) in enumerate(leading_points):
    block_started = time.monotonic()
    S = PolynomialRing(F, names=("x3", "x2", "x1", "x0"), order="degrevlex")
    x3, x2, x1, x0 = S.gens()
    ST = PolynomialRing(S, "T")
    TT = ST.gen()
    X = S(x4) * TT**4 + x3 * TT**3 + x2 * TT**2 + x1 * TT + x0
    AS = ST([S(value) for value in A])
    BS = ST([S(value) for value in B])
    rhs = X**3 + AS * X + BS

    y_coefficients = {6: S(y6)}
    for degree in range(11, 5, -1):
        index = degree - 6
        known = sum(
            y_coefficients[left] * y_coefficients[degree - left]
            for left in y_coefficients
            if degree - left in y_coefficients
        )
        y_coefficients[index] = (rhs[degree] - known) / (2 * S(y6))
    Y = sum(y_coefficients[index] * TT**index for index in range(7))
    residual = Y**2 - rhs
    equations = [S(residual[degree]) for degree in range(6)]
    assert all(equations)
    ideal = S.ideal(equations)
    dimension = int(ideal.dimension())
    assert dimension in (-1, 0)
    quotient_dimension = int(ideal.vector_space_dimension())
    solutions = ideal.variety()

    for solution in solutions:
        values = {generator: F(solution[generator]) for generator in S.gens()}
        X_value = RT([values[x0], values[x1], values[x2], values[x3], x4])
        Y_value = RT([
            F(y_coefficients[index].subs(values)) for index in range(7)
        ])
        assert Y_value**2 == X_value**3 + A * X_value + B
        key = (
            tuple(padded_coefficients(X_value, 5)),
            tuple(padded_coefficients(Y_value, 7)),
        )
        sections[key] = (X_value, Y_value)

    blocks.append({
        "index": block_index,
        "leading_point": [int(x4), int(y6)],
        "quotient_dimension_over_algebraic_closure": quotient_dimension,
        "F103_rational_solutions": len(solutions),
        "elapsed_seconds": round(time.monotonic() - block_started, 6),
    })
    print(
        f"A5A5P103BLOCK|index={block_index}|lead={int(x4)},{int(y6)}|"
        f"vdim={quotient_dimension}|rational={len(solutions)}|status=PASS_SMALL_BLOCK",
        flush=True,
    )

records = []
for index, (key, section) in enumerate(sorted(sections.items())):
    X_value, Y_value = section
    incidence = []
    for root, node_x in i6_fibres:
        value = [int(X_value(root)), int(Y_value(root))]
        hits_node = value == [int(node_x), 0]
        oriented_component = oriented_i6_component(section, root, node_x)
        incidence.append({
            "base_root": int(root),
            "nodal_x": int(node_x),
            "section_value": value,
            "identity_component": not hits_node,
            "oriented_component_in_I6_cycle": oriented_component,
            "physical_component_label": I6_ORIENTATION[int(root)][
                "labels_by_oriented_component"
            ][oriented_component],
        })
    records.append({
        "index": index,
        "X_coefficients_low_to_high": list(key[0]),
        "Y_coefficients_low_to_high": list(key[1]),
        "I6_incidence": incidence,
        "identity_at_both_I6": all(item["identity_component"] for item in incidence),
        "intersection_with_known_old_A11_affine_section": resolved_pair_intersection(section, affine),
        "is_known_old_A11_affine_section": section == affine,
    })

identity_identity = [record for record in records if record["identity_at_both_I6"]]
carrier_candidates = [
    record for record in identity_identity
    if record["intersection_with_known_old_A11_affine_section"] == 3
]
known_affine_records = [record for record in records if record["is_known_old_A11_affine_section"]]
assert len(known_affine_records) == 1
assert (len(records), len(identity_identity), len(carrier_candidates)) == (130, 18, 6)

payload = {
    "schema": "elkies-k3.q24-2a5-zero-pole-sections-p103.v1",
    "status": "PASS_BOUNDED_MOD103_ZERO_POLE_SECTION_ENUMERATION",
    "software": "SageMath 10.9 (conda-forge pinned repository environment)",
    "prime": PRIME,
    "method": {
        "leading_fibre_affine_points": len(leading_points),
        "leading_points_with_y6_zero": 0,
        "blocks": len(blocks),
        "variables_per_block": 4,
        "equations_per_block": 6,
        "upper_coefficients_solved_recursively": ["y5", "y4", "y3", "y2", "y1", "y0"],
        "single_large_Groebner_basis_used": False,
        "blocks_detail": blocks,
    },
    "surface_mod_103": {
        "A_coefficients_low_to_high": A_COEFFICIENTS,
        "B_coefficients_low_to_high": B_COEFFICIENTS,
        "I6_base_root_and_node_x": [
            [int(root), int(node)] for root, node in i6_fibres
        ],
        "I6_orientation": I6_ORIENTATION,
    },
    "counts": {
        "distinct_F103_polynomial_sections": len(records),
        "identity_at_both_I6": len(identity_identity),
        "identity_at_both_I6_and_affine_intersection_3": len(carrier_candidates),
    },
    "known_affine_section_index": known_affine_records[0]["index"],
    "P1229_mod103_candidates": [record["index"] for record in carrier_candidates],
    "sections": records,
    "elapsed_seconds": round(time.monotonic() - started, 6),
    "proof_boundary": (
        "This is exhaustive only for polynomial sections over F_103. The identity-component "
        "and affine-intersection gates are necessary fingerprints for P1229, but they do not "
        "identify or lift a characteristic-zero section. Exact QQ reconstruction and literal "
        "substitution remain open."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in EXPECTED_HASHES],
        "sha256": {
            str(path.relative_to(ROOT)): sha256(path) for path in EXPECTED_HASHES
        },
        "pinned_reduction_note": (
            "The p=103 equation and affine-section coefficients were reduced from the "
            "hash-pinned exact inputs to avoid reparsing multi-megabyte integers."
        ),
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "A5A5P103|"
    f"sections={len(records)}|identity_identity={len(identity_identity)}|"
    f"P1229_candidates={len(carrier_candidates)}|"
    f"candidate_indices={tuple(record['index'] for record in carrier_candidates)}|"
    f"elapsed={payload['elapsed_seconds']}|status={payload['status']}|output={OUTPUT}",
    flush=True,
)
