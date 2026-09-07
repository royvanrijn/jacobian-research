#!/usr/bin/env sage -python
"""Point the reconstructed q8/orbit376 quartic using known exact parent sections.

The unpointed RR compiler emits a quartic

    W^2 = Q(t,U)

over QQ(U), where t is the old q4/orbit164 base and U is the new q8 base.
Any exact parent section C whose restriction U|_C is Mobius gives a rational
point (t(U),W(U)) on this quartic. Translating the quartic at that point
attaches an exact rational origin to the invariant 4A1 Jacobian.

This script tests the exact polynomial B0,...,B7 basis, the certified C8
opposite section, and the two transported primitive one-node sections when
available. It does not assume that P1229 is a section of the parent
fibration. The selected origin is recorded by its actual equation curve;
reframing it to the preferred P1229-zero marked child is a separate NS gate.
"""

import hashlib
import json
import time
from math import comb
from pathlib import Path

from sage.all import PolynomialRing, QQ, ZZ


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
MODEL = LOCAL / "q4o164-compact-weierstrass-qq.json"
HORIZONTAL = LOCAL / "q4o164-q8o376-horizontal-crt-qq.json"
RR = LOCAL / "q4o164-q8o376-rr-p2-qq.json"
BASIS = LOCAL / "q4o164-integral-basis-qq.json"
C8 = LOCAL / "q4o164-c8-equation-marking-qq.json"
TRANSPORTED = LOCAL / "q4o1584-degree1-sections-to-q4o164-qq.json"
OUTPUT = LOCAL / "q4o164-q8o376-known-section-pointing-qq.json"

started = time.monotonic()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def coefficient_bits(value):
    value = QQ(value)
    return max(abs(ZZ(value.numerator())).nbits(), value.denominator().nbits())


def rational_bits(value):
    return max(
        coefficient_bits(coefficient)
        for polynomial in (value.numerator(), value.denominator())
        for coefficient in polynomial.list()
    )


def rational_record(value, ring):
    value = ring.fraction_field()(value)
    numerator = ring(value.numerator())
    denominator = ring(value.denominator())
    common = numerator.gcd(denominator)
    numerator //= common
    denominator //= common
    scale = denominator.leading_coefficient()
    numerator /= scale
    denominator /= scale
    return {
        "numerator_coefficients_low_to_high": [str(value) for value in numerator.list()],
        "denominator_coefficients_low_to_high": [str(value) for value in denominator.list()],
        "degrees_numerator_denominator": [
            int(numerator.degree()), int(denominator.degree())
        ],
        "maximum_rational_bits": max(
            coefficient_bits(value)
            for value in list(numerator) + list(denominator)
        ),
    }


for path in (MODEL, HORIZONTAL, RR, BASIS):
    if not path.exists():
        raise SystemExit(f"missing prerequisite: {path}")

model = json.loads(MODEL.read_text())
horizontal = json.loads(HORIZONTAL.read_text())
rr = json.loads(RR.read_text())
basis = json.loads(BASIS.read_text())
assert model["status"] == "PASS_EXACT_QQ_Q4O164_COMPACT_WEIERSTRASS_NORMALIZATION"
assert horizontal["status"] == "PASS_EXACT_QQ_Q4O164_Q8O376_HORIZONTAL_CRT"
assert rr["status"] == "PASS_EXACT_QQ_Q4O164_Q8O376_UNPOINTED_RR_AND_4A1_JACOBIAN"
assert basis["status"] in {
    "PASS_EXACT_QQ_Q4O164_PAIR_NODE_SECTION_SUBGROUP_RANK8",
    "PASS_EXACT_QQ_Q4O164_INTEGRAL_SECTION_BASIS_RANK9",
}

RT = PolynomialRing(QQ, "t")
t = RT.gen()
KT = RT.fraction_field()
RU = PolynomialRing(QQ, "U")
U = RU.gen()
KU = RU.fraction_field()

A = RT([QQ(value) for value in model["compact_model"]["A_coefficients_low_to_high"]])
B = RT([QQ(value) for value in model["compact_model"]["B_coefficients_low_to_high"]])
A_child = RU([QQ(value) for value in rr["child"]["A_coefficients_low_to_high"]])
B_child = RU([QQ(value) for value in rr["child"]["B_coefficients_low_to_high"]])

section = horizontal["section"]


def read_rational(record, ring):
    numerator = ring([QQ(value) for value in record["numerator_coefficients_low_to_high"]])
    denominator = ring([QQ(value) for value in record["denominator_coefficients_low_to_high"]])
    return ring.fraction_field()(numerator) / ring.fraction_field()(denominator)


xH = read_rational(section["x"], RT)
yH = read_rational(section["y"], RT)
assert yH**2 == xH**3 + KT(A) * xH + KT(B)
xH_denominator = RT(xH.denominator())
if not xH_denominator or xH_denominator.leading_coefficient() != 1:
    raise ArithmeticError("horizontal x denominator is not monic")
Z = RT.one()
for factor, multiplicity in xH_denominator.factor():
    if int(multiplicity) % 2:
        raise ArithmeticError("horizontal x denominator is not a square")
    Z *= factor.monic() ** (int(multiplicity) // 2)
if Z**2 != xH_denominator or Z.degree() != 4:
    raise ArithmeticError("failed to recover the degree-four horizontal projective Z")

basis_rows = rr["resolved_RR"]["basis"]
if len(basis_rows) != 2:
    raise ArithmeticError("unpointed RR artifact does not contain a two-plane")
nonconstant = basis_rows[1]
AA = RT([QQ(value) for value in nonconstant["AA_coefficients_low_to_high"]])
BB = RT([QQ(value) for value in nonconstant["BB_coefficients_low_to_high"]])
square_root = RT([
    QQ(value)
    for value in rr["resolved_RR"][
        "finite_vertical_square_root_coefficients_low_to_high"
    ]
])
if not BB or BB.degree() > 2:
    raise ArithmeticError("stored nonconstant RR direction has invalid BB")
if not square_root:
    square_root = RT.one()

quartic_coefficients = [
    RU([QQ(value) for value in row])
    for row in rr["quartic"]["coefficients_in_old_t_low_to_high"]
]
if len(quartic_coefficients) != 5:
    raise ArithmeticError("stored q8 quartic is not degree four in old t")


def evaluate_in_u(polynomial, argument):
    target = argument.parent()
    return sum(
        target(coefficient) * argument**index
        for index, coefficient in enumerate(polynomial)
    )


def substitute_rational(value, argument, source_ring):
    value = source_ring.fraction_field()(value)
    target = argument.parent()
    numerator = sum(
        target(coefficient) * argument**index
        for index, coefficient in enumerate(source_ring(value.numerator()))
    )
    denominator = sum(
        target(coefficient) * argument**index
        for index, coefficient in enumerate(source_ring(value.denominator()))
    )
    return numerator / denominator


candidates = []


def append_candidate(name, x_value, y_value, provenance):
    x_value = KT(x_value)
    y_value = KT(y_value)
    if y_value**2 != x_value**3 + KT(A) * x_value + KT(B):
        raise ArithmeticError(f"{name} misses the compact q4/orbit164 model")
    for sign in (1, -1):
        y_signed = sign * y_value
        if x_value == xH and y_signed in (yH, -yH):
            continue
        chord = (y_signed - yH) / (x_value - xH)
        new_base = KT(AA) / KT(Z**2) + KT(BB) / KT(Z) * chord
        numerator = RT(new_base.numerator())
        denominator = RT(new_base.denominator())
        common = numerator.gcd(denominator)
        numerator //= common
        denominator //= common
        degree = max(int(numerator.degree()), int(denominator.degree()))

        raw_ordinate = 2 * x_value + xH - chord**2
        ordinate = raw_ordinate * KT(BB**2) / (KT(Z) * KT(square_root))
        quartic_on_curve = sum(
            evaluate_in_u(coefficient, new_base) * KT(t)**index
            for index, coefficient in enumerate(quartic_coefficients)
        )
        if ordinate**2 != quartic_on_curve:
            raise ArithmeticError(f"{name} sign {sign} misses the stored q8 quartic")

        record = {
            "name": name,
            "sign": sign,
            "provenance": provenance,
            "new_base_degree": degree,
            "new_base": rational_record(new_base, RT),
            "quartic_ordinate_on_parent_curve": rational_record(ordinate, RT),
            "exact_parent_section_identity": True,
            "exact_surface_to_quartic_identity": True,
        }
        if degree == 1:
            if numerator.degree() > 1 or denominator.degree() > 1:
                raise ArithmeticError("degree-one reduction did not remain linear")
            n0, n1 = QQ(numerator[0]), QQ(numerator[1])
            d0, d1 = QQ(denominator[0]), QQ(denominator[1])
            old_base = KU(
                (KU(n0) - KU(U) * KU(d0))
                / (KU(U) * KU(d1) - KU(n1))
            )
            if substitute_rational(new_base, old_base, RT) != KU(U):
                raise ArithmeticError("Mobius inversion failed")
            W0 = substitute_rational(ordinate, old_base, RT)
            quartic_KU = [KU(value) for value in quartic_coefficients]
            if W0**2 != sum(
                coefficient * old_base**index
                for index, coefficient in enumerate(quartic_KU)
            ):
                raise ArithmeticError("inverted degree-one point misses the quartic")
            if not W0:
                raise ArithmeticError("degree-one quartic point has zero ordinate")

            translated = [
                sum(
                    quartic_KU[index] * QQ(comb(index, order))
                    * old_base**(index - order)
                    for index in range(order, 5)
                )
                for order in range(5)
            ]
            e, d, c, b, a = translated
            if e != W0**2:
                raise ArithmeticError("quartic translation constant changed")
            a1 = d / W0
            a2 = c - d**2 / (4 * W0**2)
            a3 = 2 * W0 * b
            a4 = -4 * W0**2 * a
            a6 = a2 * a4
            b2 = a1**2 + 4 * a2
            b4 = 2 * a4 + a1 * a3
            b6 = a3**2 + 4 * a6
            c4 = b2**2 - 24 * b4
            c6 = -b2**3 + 36 * b2 * b4 - 216 * b6
            if KU(81 * (-c4 / 48)) != KU(A_child):
                raise ArithmeticError("pointed quartic c4 misses invariant child")
            if KU(729 * (-c6 / 864)) != KU(B_child):
                raise ArithmeticError("pointed quartic c6 misses invariant child")
            record.update({
                "old_base_as_function_of_new_base": rational_record(old_base, RU),
                "quartic_ordinate_as_function_of_new_base": rational_record(W0, RU),
                "pointed_generalized_coefficients": {
                    "a1": rational_record(a1, RU),
                    "a2": rational_record(a2, RU),
                    "a3": rational_record(a3, RU),
                    "a4": rational_record(a4, RU),
                    "a6": rational_record(a6, RU),
                },
                "exact_pointed_child_invariant_identity": True,
                "maximum_pointing_rational_bits": max(
                    rational_bits(old_base), rational_bits(W0),
                    rational_bits(a1), rational_bits(a2), rational_bits(a3),
                    rational_bits(a4), rational_bits(a6),
                ),
            })
        candidates.append(record)


# Exact polynomial B-basis sections on the compact model.
for index, record in enumerate(basis["resolved_hensel"]["sections"]):
    x_value = RT([QQ(value) for value in record["x_coefficients_low_to_high"]])
    y_value = RT([QQ(value) for value in record["y_coefficients_low_to_high"]])
    append_candidate(
        f"B{index}",
        x_value,
        y_value,
        "q4o164-integral-basis-qq.json resolved Hensel section",
    )

# Certified C8-opposite section, stored before the compact rescaling.
if C8.exists():
    c8 = json.loads(C8.read_text())
    if c8.get("status") == "PASS_EXACT_QQ_Q4O164_C8_EQUATION_MARKING":
        old = c8["opposite_constant_support_section"]
        x_old = read_rational(old["x"], RT)
        y_old = read_rational(old["y"], RT)
        base_scale = QQ(model["exact_coordinate_change"]["c"])
        xy_scale = QQ(model["exact_coordinate_change"]["s"])
        x_compact = (
            substitute_rational(x_old, KT(base_scale * t), RT) / xy_scale**2
        )
        y_compact = (
            substitute_rational(y_old, KT(base_scale * t), RT) / xy_scale**3
        )
        append_candidate(
            "C8opposite",
            x_compact,
            y_compact,
            "certified C8 opposite section transported to compact coordinates",
        )

# Two exact primitive one-node sections already stored in compact coordinates.
if TRANSPORTED.exists():
    transported = json.loads(TRANSPORTED.read_text())
    if transported.get("status") == "PASS_EXACT_QQ_Q4O164_TWO_PRIMITIVE_ONE_NODE_SECTIONS":
        for index, row in enumerate(transported["degree_one_sections"]):
            append_candidate(
                f"transported_one_node_{index}",
                read_rational(row["compact_child_x"], RT),
                read_rational(row["compact_child_y"], RT),
                "degree-one q4/orbit1584 curve transported through the pointed q4/orbit164 quartic",
            )

degree_one = [row for row in candidates if row["new_base_degree"] == 1]
degree_one.sort(key=lambda row: (
    row["maximum_pointing_rational_bits"], row["name"], row["sign"]
))
selected = degree_one[0] if degree_one else None
status = (
    "PASS_EXACT_QQ_Q4O164_Q8O376_KNOWN_SECTION_POINTING"
    if selected
    else "NO_KNOWN_SECTION_Q4O164_Q8O376_DEGREE_ONE_POINTING"
)

input_paths = [MODEL, HORIZONTAL, RR, BASIS] + [
    path for path in (C8, TRANSPORTED) if path.exists()
]
payload = {
    "schema": "elkies-k3.q4o164-q8o376-known-section-pointing-qq.v1",
    "status": status,
    "candidate_count": len(candidates),
    "restriction_degrees": [
        {
            "name": row["name"],
            "sign": row["sign"],
            "new_base_degree": row["new_base_degree"],
        }
        for row in candidates
    ],
    "degree_one_count": len(degree_one),
    "degree_one_candidates": degree_one,
    "selected_origin": selected,
    "child": {
        "A_coefficients_low_to_high": [str(value) for value in A_child.list()],
        "B_coefficients_low_to_high": [str(value) for value in B_child.list()],
        "ADE": "4A1",
    },
    "method": {
        "large_Groebner_required": False,
        "new_section_ansatz_required": False,
        "construction": (
            "exact restriction of the reconstructed q8 RR pencil to known parent "
            "sections, Mobius inversion, and pointed binary-quartic invariants"
        ),
        "runtime_seconds": time.monotonic() - started,
    },
    "proof_boundary": (
        "A PASS attaches an exact rational origin to the reconstructed 4A1 Jacobian "
        "using the displayed equation curve. It does not identify that curve with the "
        "preferred P1229 zero, reframe the child NS lattice, or certify the final "
        "q8/orbit376 marked equation edge. A NO_KNOWN result is exhaustive only for "
        "the exact section list loaded by this script."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in input_paths],
        "sha256": {
            str(path.relative_to(ROOT)): sha256(path) for path in input_paths
        },
    },
}
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "Q8O376POINT|tested={}|degree1={}|selected={}|status={}|output={}".format(
        len(candidates),
        len(degree_one),
        None if selected is None else f"{selected['name']}:{selected['sign']}",
        status,
        OUTPUT,
    ),
    flush=True,
)
