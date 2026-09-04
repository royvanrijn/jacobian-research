#!/usr/bin/env sage-python
"""Generate fresh target-free A1/MW16 specializations for a blind search."""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
from importlib.machinery import SourceFileLoader
import json
from math import gcd
from pathlib import Path

from sage.all import EllipticCurve, PolynomialRing, QQ, ZZ, vector


ROOT = Path(__file__).resolve().parents[2]
TEMPLATE = ROOT / "elliptic-curves/data/a1_mw16_family_template_v1.json"
MODEL = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json"
CHORD = ROOT / "elkies-k3/scripts/construct_elkies_2026_bisections.sage"
BLIND = ROOT / "elliptic-curves/data/a1_mw16_blind_parameter_experiment_v1.json"
KEY = ROOT / "artifacts/local/elliptic-curves/a1-mw16-parameter-experiment/unblinding-v1.json"
SALT = "a1-mw16-target-free-parameter-experiment-20260904-v1"


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def qtext(value) -> str:
    value = QQ(value)
    return str(value.numerator()) if value.denominator() == 1 else f"{value.numerator()}/{value.denominator()}"


def poly_from_record(record, ring):
    numerator = ring([QQ(value) for value in record["numerator_coefficients_low_to_high"]])
    denominator = ring([QQ(value) for value in record["denominator_coefficients_low_to_high"]])
    return ring.fraction_field()(numerator) / ring.fraction_field()(denominator)


def child_geometry(trace, old_a, old_ring, chord):
    frame = chord.trace_chord_frame(trace[0], trace[1], old_ring)
    h, nx, ny, m0 = (frame[key] for key in ("h", "Nx", "Ny", "M0"))
    parameter_ring = PolynomialRing(QQ, "lambda")
    parameter = parameter_ring.gen()
    bivariate = PolynomialRing(parameter_ring, "t")
    hh, nnx, nny, mm0 = map(bivariate, (h, nx, ny, m0))
    slope = mm0 + parameter * hh**2
    numerator = slope**4 - 6*slope**2*nnx - 8*slope*nny - 3*nnx**2 - 4*bivariate(old_a)*hh**4
    quartic, remainder = numerator.quo_rem(hh**6)
    if remainder or quartic.degree() != 4:
        raise ArithmeticError("canonical class did not produce a quartic pencil")
    coefficients = [parameter_ring(quartic[index]) for index in range(5)]
    e, d, c, b, a = coefficients
    invariant_i = 12*a*e - 3*b*d + c**2
    invariant_j = 72*a*c*e + 9*b*c*d - 27*a*d**2 - 27*b**2*e - 2*c**3
    return h, nx, m0, quartic, -27*invariant_i, -27*invariant_j


def specialize(parameter, selected, new_zero, trace, old_basis, old_curve, old_a, old_ring, h, nx, m0, quartic, child_a, child_b):
    fixed_m = m0 + parameter*h**2
    fixed_quartic = old_ring([QQ(quartic[index](parameter)) for index in range(5)])
    sum_x = old_ring((fixed_m**2 - nx) // h**2)
    quartic_points = []
    for section_vector in (new_zero,) + selected:
        point = sum((coefficient*basis_point for coefficient, basis_point in zip(section_vector, old_basis) if coefficient), old_curve(0))
        base_map = (((point[1] + trace[1]) / (point[0] - trace[0])) * h - m0) / h**2
        old_parameter = old_curve.base_field()(base_map).numerator() - parameter * old_curve.base_field()(base_map).denominator()
        roots = old_ring(old_parameter).roots(QQ)
        if len(roots) != 1 or roots[0][1] != 1:
            raise ArithmeticError("section base map did not invert uniquely")
        old_value = QQ(roots[0][0])
        x_value = QQ(point[0](old_value))
        y_value = QQ(point[1](old_value))
        w_value = (2*x_value - QQ(sum_x(old_value))) / QQ(h(old_value))
        if w_value**2 != fixed_quartic(old_value):
            raise ArithmeticError("specialized section missed the quartic")
        quartic_points.append((old_value, w_value))

    t0, w0 = quartic_points[0]
    if not w0:
        raise ArithmeticError("pointed quartic origin is singular")
    shift_ring = PolynomialRing(QQ, "z")
    shifted = shift_ring(fixed_quartic(t0 + shift_ring.gen()))
    ee, dd, cc, bb, aa = [QQ(shifted[index]) for index in range(5)]
    if ee != w0**2:
        raise ArithmeticError("pointed quartic constant changed")
    a1g = dd/w0
    a2g = cc - dd**2/(4*w0**2)
    a3g = 2*w0*bb
    a4g = -4*w0**2*aa
    a6g = a2g*a4g
    b2g = a1g**2 + 4*a2g
    specialized_a = QQ(child_a(parameter))
    specialized_b = QQ(child_b(parameter))
    curve = EllipticCurve(QQ, [specialized_a, specialized_b])
    if curve.discriminant() == 0:
        raise ArithmeticError("singular specialization")
    points = []
    for old_value, w_value in quartic_points[1:]:
        z = old_value - t0
        if not z:
            raise ArithmeticError("generic section met the new zero")
        xg = (2*w0*(w_value+w0) + dd*z)/z**2
        yg = (4*w0**2*(w_value+w0) + 2*w0*dd*z + (2*w0*cc-dd**2/(2*w0))*z**2)/z**3
        points.append(curve(9*(xg+b2g/12), 27*(yg+(a1g*xg+a3g)/2)))
    if len(points) != 16:
        raise ArithmeticError("specialized MW16 basis is incomplete")
    return curve, tuple(points)


def parameter_stream():
    counter = 0
    seen = set()
    while True:
        raw = sha256(f"{SALT}:{counter}".encode()).digest()
        counter += 1
        numerator = int.from_bytes(raw[:8], "big") % 1023 - 511
        denominator = int.from_bytes(raw[8:16], "big") % 511 + 1
        if gcd(abs(numerator), denominator) != 1:
            continue
        value = QQ(numerator) / denominator
        if value in seen:
            continue
        seen.add(value)
        yield value, counter - 1


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--count", type=int, default=32)
    parser.add_argument("--blind-output", type=Path, default=BLIND)
    parser.add_argument("--key-output", type=Path, default=KEY)
    args = parser.parse_args()
    if not 1 <= args.count <= 256:
        raise SystemExit("--count must be between 1 and 256")
    template = json.loads(TEMPLATE.read_text())
    if template.get("status") != "PASS_TARGET_FREE_A1_MW16_FAMILY_TEMPLATE":
        raise ArithmeticError("family template status changed")
    model = json.loads(MODEL.read_text())
    chord = SourceFileLoader("a1_mw16_parameter_chord", str(CHORD)).load_module()
    old_ring = PolynomialRing(QQ, "t")
    old_field = old_ring.fraction_field()
    old_a = old_ring(model["weierstrass_model"]["A_coefficients_low_to_high"])
    old_b = old_ring(model["weierstrass_model"]["B_coefficients_low_to_high"])
    old_curve = EllipticCurve(old_field, [old_a, old_b])
    old_basis = tuple(old_curve(poly_from_record(row["X"], old_ring), poly_from_record(row["Y"], old_ring)) for row in model["sections"]["records"])
    marking = template["source_marking"]
    trace_vector = vector(ZZ, marking["trace_section_basis_w"])
    trace = sum((coefficient*point for coefficient, point in zip(trace_vector, old_basis) if coefficient), old_curve(0))
    selected = tuple(vector(ZZ, row) for row in marking["generic_source_section_basis_coordinates"])
    new_zero = vector(ZZ, marking["new_zero_source_section_basis_coordinates"])
    h, nx, m0, quartic, child_a, child_b = child_geometry(trace, old_a, old_ring, chord)
    if [qtext(child_a[index]) for index in range(child_a.degree()+1)] != template["pencil"]["A_coefficients_low_to_high"]:
        raise ArithmeticError("template A polynomial changed")
    if [qtext(child_b[index]) for index in range(child_b.degree()+1)] != template["pencil"]["B_coefficients_low_to_high"]:
        raise ArithmeticError("template B polynomial changed")

    rows = []
    key_rows = []
    j_values = set()
    stream = parameter_stream()
    while len(rows) < args.count:
        parameter, counter = next(stream)
        try:
            curve, points = specialize(parameter, selected, new_zero, trace, old_basis, old_curve, old_a, old_ring, h, nx, m0, quartic, child_a, child_b)
        except (ArithmeticError, ZeroDivisionError):
            continue
        j_value = curve.j_invariant()
        if j_value in j_values:
            continue
        j_values.add(j_value)
        row_hash = sha256(f"{SALT}:row:{counter}:{qtext(parameter)}".encode()).hexdigest()
        row_id = f"a1-{row_hash[:16]}"
        rows.append({
            "row_id": row_id,
            "short_model": [qtext(value) for value in curve.a_invariants()],
            "specialized_generic_points": [{"x": qtext(point[0]), "y": qtext(point[1])} for point in points],
            "generic_height_gram": template["generic_height_gram"],
        })
        key_rows.append({"row_id": row_id, "parameter": qtext(parameter), "sampler_counter": counter})
        print(f"A1MW16PREP|row={len(rows)}/{args.count}|id={row_id}|status=PASS", flush=True)

    blind = {
        "schema": "elliptic-curves.a1-mw16-blind-parameter-experiment.v1",
        "status": "PASS_TARGET_FREE_BLIND_INPUT",
        "family_id": template["family_id"],
        "row_count": len(rows),
        "rows": rows,
        "blindness": {
            "parameter_in_worker_rows": False,
            "known_target_data_loaded": False,
            "public_points_loaded": False,
            "public_ranks_loaded": False,
        },
        "inputs": {str(path.relative_to(ROOT)): digest(path) for path in (TEMPLATE, MODEL, CHORD, Path(__file__))},
    }
    key = {
        "schema": "elliptic-curves.a1-mw16-parameter-unblinding.v1",
        "status": "SEALED_UNTIL_SEARCH_MERGE",
        "salt": SALT,
        "rows": key_rows,
    }
    args.blind_output.parent.mkdir(parents=True, exist_ok=True)
    args.key_output.parent.mkdir(parents=True, exist_ok=True)
    args.blind_output.write_text(json.dumps(blind, indent=2, sort_keys=True)+"\n")
    args.key_output.write_text(json.dumps(key, indent=2, sort_keys=True)+"\n")
    print(f"A1MW16PREP|rows={len(rows)}|blind_sha256={digest(args.blind_output)}|status=PASS_TARGET_FREE_BLIND_INPUT")


if __name__ == "__main__":
    main()
