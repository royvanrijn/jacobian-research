#!/usr/bin/env sage -python
"""Identify the two transported q4/orbit164 degree-one sections in marked MW9.

The q4/o164 exact equation side contains an exact rank-eight section basis
B0,...,B7.  Its height Gram and the selected embedding (the same embedding
used by the q8/o376 Abel-trace horizontal certificate) place those eight
sections inside the C8-pointed marked MW9 lattice.

For each transported one-node section P we compute, using exact function-field
group law only,

    h(P), h(P+B_i), i=0,...,7.

Polarization gives all eight pairings <P,B_i>.  Enumerating the tiny shell of
4*MW9 with norm 4*h(P) then identifies the unique marked MW tail compatible
with those pairings.  The tail is compared directly with the named inherited
curve classes in the physical q4/o164 marking, in particular P1229.

No resolved I4 chart, Groebner basis, or nonlinear section solve is used.
"""

import hashlib
import json
import time
from pathlib import Path

from sage.all import EllipticCurve, PolynomialRing, QQ, ZZ, matrix, pari, vector


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
GENERATED = ROOT / "artifacts/generated-results"
MODEL = LOCAL / "q4o164-compact-weierstrass-qq.json"
BASIS = LOCAL / "q4o164-integral-basis-qq.json"
HEIGHT_AUDIT = LOCAL / "q4o164-integral-basis-height-gram-audit-qq.json"
TRANSPORTED = LOCAL / "q4o1584-degree1-sections-to-q4o164-qq.json"
MARKING = GENERATED / "elkies-k3-h3-q4o208-q4o1584-q4o164-old_a11_component_8-marking.json"
HORIZONTAL = LOCAL / "q4o164-q8o376-horizontal-crt-qq.json"
OUTPUT = LOCAL / "q4o164-transported-degree1-marked-classes-qq.json"
INPUTS = (MODEL, BASIS, HEIGHT_AUDIT, TRANSPORTED, MARKING, HORIZONTAL)

started = time.monotonic()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_matrix(path):
    return matrix(ZZ, [
        [ZZ(value) for value in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ])


def read_rational(record, ring):
    numerator = ring([QQ(value) for value in record["numerator_coefficients_low_to_high"]])
    denominator = ring([QQ(value) for value in record["denominator_coefficients_low_to_high"]])
    return ring.fraction_field()(numerator) / ring.fraction_field()(denominator)


for path in INPUTS:
    if not path.exists():
        raise SystemExit(f"missing prerequisite: {path}")

model = json.loads(MODEL.read_text())
basis = json.loads(BASIS.read_text())
audit = json.loads(HEIGHT_AUDIT.read_text())
transported = json.loads(TRANSPORTED.read_text())
marking = json.loads(MARKING.read_text())
horizontal = json.loads(HORIZONTAL.read_text())

assert model["status"] == "PASS_EXACT_QQ_Q4O164_COMPACT_WEIERSTRASS_NORMALIZATION"
assert basis["status"] == "PASS_EXACT_QQ_Q4O164_PAIR_NODE_SECTION_SUBGROUP_RANK8"
assert audit["status"] == "PASS_EXACT_QQ_Q4O164_FOURFOLD_HEIGHT_GRAM_AND_C8_MARKED_EMBEDDING_CENSUS"
assert transported["status"] == "PASS_EXACT_QQ_Q4O164_TWO_PRIMITIVE_ONE_NODE_SECTIONS"
assert marking["status"] == "PASS_EXACT_Q4O164_PHYSICAL_EFFECTIVE_ZERO_MARKING"
assert marking["zero"] == "old_A11_component_8"
assert horizontal["status"] == "PASS_EXACT_QQ_Q4O164_Q8O376_HORIZONTAL_CRT"

# The exact horizontal was reconstructed only after the modular Abel trace
# selected embedding 15.  Lock the same embedding here; this is not a new
# choice of marking.
embedding_index = 15
embedding_record = next(
    record
    for record in audit["marked_embedding_enumeration"]["embeddings"]
    if int(record["embedding_index"]) == embedding_index
)
if not embedding_record["compatible_with_first_seven_stored_profiles_up_to_fibre_symmetry"]:
    raise ArithmeticError("selected embedding is not component-compatible")
if not embedding_record["contains_q8_residual"]:
    raise ArithmeticError("selected embedding does not contain the certified q8 residual")
B_marked = matrix(ZZ, embedding_record["rows_B0_through_B7_in_marked_MW9"])
assert B_marked.nrows() == 8 and B_marked.ncols() == 9 and B_marked.rank() == 8

frame = load_matrix(ROOT / marking["frame_output"])
assert frame.nrows() == frame.ncols() == 17
root_gram = frame[:8, :8]
coupling = frame[:8, 8:]
marked_height = frame[8:, 8:].change_ring(QQ) - (
    coupling.transpose() * root_gram.inverse() * coupling
)
assert marked_height.nrows() == marked_height.ncols() == 9
scaled_height = (4 * marked_height).change_ring(ZZ)
assert scaled_height.is_positive_definite()

R = PolynomialRing(QQ, "t")
t = R.gen()
K = R.fraction_field()
A = R([QQ(value) for value in model["compact_model"]["A_coefficients_low_to_high"]])
B = R([QQ(value) for value in model["compact_model"]["B_coefficients_low_to_high"]])
E = EllipticCurve(K, [0, 0, 0, K(A), K(B)])


def basis_point(record):
    x = R([QQ(value) for value in record["x_coefficients_low_to_high"]])
    y = R([QQ(value) for value in record["y_coefficients_low_to_high"]])
    return E(K(x), K(y))


def transported_point(row):
    return E(
        read_rational(row["compact_child_x"], R),
        read_rational(row["compact_child_y"], R),
    )


def fourfold_height(point):
    fourfold = 4 * point
    x = fourfold[0]
    y = fourfold[1]
    x_num, x_den = R(x.numerator()), R(x.denominator())
    y_num, y_den = R(y.numerator()), R(y.denominator())
    pole_degree = max(x_den.degree(), x_num.degree() - 4)
    if pole_degree < 0 or pole_degree % 2:
        raise ArithmeticError(f"invalid fourfold pole degree {pole_degree}")
    if x_den**3 != y_den**2:
        raise ArithmeticError("fourfold x/y denominators have incompatible powers")
    expected = (
        4 + pole_degree,
        pole_degree,
        6 + 3 * pole_degree // 2,
        3 * pole_degree // 2,
    )
    actual = (x_num.degree(), x_den.degree(), y_num.degree(), y_den.degree())
    if actual != expected:
        raise ArithmeticError(f"fourfold degree fingerprint changed: {actual} != {expected}")
    return QQ(4 + pole_degree) / 16, int(pole_degree), list(map(int, actual))


B_points = [basis_point(record) for record in basis["resolved_hensel"]["sections"]]
assert len(B_points) == 8
B_heights = [fourfold_height(point)[0] for point in B_points]
for index, value in enumerate(B_heights):
    expected = QQ(audit["fourfold_height_audit"]["height_gram"][index][index])
    if value != expected:
        raise ArithmeticError(f"B{index} height mismatch: {value} != {expected}")

explicit = {
    name: vector(ZZ, values)
    for name, values in marking["equation_explicit_curves_in_child"].items()
}
# A section has old-fibre degree one, i.e. coefficient 1 against F=e0 in the
# U block.  Its final nine coordinates are exactly its marked MW tail.
named_section_tails = {}
for name, curve in explicit.items():
    if len(curve) != 19:
        continue
    if curve[1] == 1:
        named_section_tails.setdefault(tuple(map(int, curve[-9:])), []).append(name)


def shell_vectors(target_norm):
    target_norm = ZZ(target_norm)
    if target_norm <= 0:
        return []
    raw = matrix(ZZ, pari(scaled_height).qfminim(target_norm)[2]).columns()
    candidates = [vector(ZZ, column) for column in raw]
    candidates += [-candidate for candidate in candidates]
    unique = {}
    for candidate in candidates:
        if candidate * scaled_height * candidate != target_norm:
            continue
        unique[tuple(map(int, candidate))] = candidate
    return list(unique.values())


results = []
for index, row in enumerate(transported["degree_one_sections"]):
    P = transported_point(row)
    hP, pole_degree, degree_fingerprint = fourfold_height(P)
    pairings = []
    sum_heights = []
    for b_index, (Bpoint, hB) in enumerate(zip(B_points, B_heights)):
        hsum, unused_pole, unused_degree = fourfold_height(P + Bpoint)
        pairing = (hsum - hP - hB) / 2
        pairings.append(pairing)
        sum_heights.append(hsum)

    scaled_norm = 4 * hP
    if scaled_norm.denominator() != 1:
        raise ArithmeticError(f"transported section {index} has nonintegral scaled norm {scaled_norm}")
    shell = shell_vectors(ZZ(scaled_norm))
    compatible = []
    for candidate in shell:
        if all(
            candidate * marked_height * B_marked.row(b_index) == pairings[b_index]
            for b_index in range(8)
        ):
            compatible.append(candidate)
    if len(compatible) != 1:
        raise ArithmeticError(
            f"transported section {index} has {len(compatible)} marked tails; "
            f"shell={len(shell)}, h={hP}, pairings={pairings}"
        )
    tail = compatible[0]
    direct_names = sorted(named_section_tails.get(tuple(map(int, tail)), []))
    negative_names = sorted(named_section_tails.get(tuple(map(int, -tail)), []))
    results.append({
        "transported_index": index,
        "pair_index": row.get("pair_index"),
        "sign": row.get("sign"),
        "canonical_height": str(hP),
        "fourfold_pole_degree": pole_degree,
        "fourfold_degree_fingerprint": degree_fingerprint,
        "height_of_P_plus_B0_through_B7": [str(value) for value in sum_heights],
        "pairing_with_B0_through_B7": [str(value) for value in pairings],
        "enumerated_exact_norm_shell_size": len(shell),
        "marked_MW9_tail": list(map(int, tail)),
        "direct_named_section_matches": direct_names,
        "negative_named_section_matches": negative_names,
        "matches_P1229": "P1229" in direct_names,
        "negative_matches_P1229": "P1229" in negative_names,
    })
    print(
        "Q4O164D1MARK|index={}|height={}|shell={}|tail={}|names={}|neg_names={}".format(
            index, hP, len(shell), list(map(int, tail)), direct_names, negative_names
        ),
        flush=True,
    )

p1229_matches = [record for record in results if record["matches_P1229"]]
negative_p1229_matches = [record for record in results if record["negative_matches_P1229"]]
status = (
    "PASS_EXACT_QQ_Q4O164_TRANSPORTED_DEGREE1_P1229_IDENTIFIED"
    if len(p1229_matches) == 1
    else "PASS_EXACT_QQ_Q4O164_TRANSPORTED_DEGREE1_MARKED_CLASSES"
)

payload = {
    "schema": "elkies-k3.q4o164-transported-degree1-marked-classes-qq.v1",
    "status": status,
    "selected_embedding_index": embedding_index,
    "marked_height_determinant": str(marked_height.det()),
    "transported_section_count": len(results),
    "sections": results,
    "P1229_direct_match_count": len(p1229_matches),
    "P1229_negative_match_count": len(negative_p1229_matches),
    "P1229_transport_index": (
        p1229_matches[0]["transported_index"] if len(p1229_matches) == 1 else None
    ),
    "method": {
        "large_Groebner_required": False,
        "resolved_component_chart_required": False,
        "construction": (
            "exact fourfold pole heights, polarization against B0..B7, and finite shell "
            "enumeration in the selected C8-pointed marked MW9 lattice"
        ),
        "runtime_seconds": time.monotonic() - started,
    },
    "proof_boundary": (
        "Each transported exact q4/o164 section is assigned a unique marked MW9 tail by "
        "its exact height and eight exact pairings with the certified equation basis. "
        "A direct named match identifies the same inherited curve class in the physical "
        "marking. This identifies curve classes; it does not yet compile or point q8/o376."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS],
        "sha256": {str(path.relative_to(ROOT)): sha256(path) for path in INPUTS},
    },
}
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    f"Q4O164D1MARK|P1229={len(p1229_matches)}|status={status}|output={OUTPUT}",
    flush=True,
)
