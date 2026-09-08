#!/usr/bin/env sage -python
"""Second bounded direct determinant slice in E302's maximal resolvent fibre."""

from __future__ import annotations

import argparse
import importlib.util
import itertools
import json
from pathlib import Path
from importlib.machinery import SourceFileLoader

from sage.all import QQ, ZZ, PolynomialRing, matrix, vector


HERE = Path(__file__).resolve()
ROOT = HERE.parents[2]
ART = ROOT / "artifacts" / "generated-results" / "elliptic-curves"
PROTOCOL = HERE.with_name("CURVE302_FIXED_RESOLVENT_ORBIT_PROTOCOL_V2.json")
OUTPUT = ART / "curve302_fixed_resolvent_orbits_v2.json"
WORK = ROOT / "artifacts" / "local" / "elliptic-curves" / "curve302-fixed-resolvent-orbits-v2"
V1 = HERE.with_name("curve302_fixed_resolvent_orbits.sage")

spec = importlib.util.spec_from_loader("curve302_fixed_resolvent_orbits_v1", SourceFileLoader("curve302_fixed_resolvent_orbits_v1", str(V1)))
base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)
base.PROTOCOL = PROTOCOL


def read(path):
    return json.loads(path.read_text())


def put_new(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")


def coeff(poly, x, y, degree_x):
    return QQ(poly.monomial_coefficient(x**degree_x * y**(1 - degree_x)))


def degree_three_coefficients(poly, x, y):
    return [QQ(poly.monomial_coefficient(x**(3 - i) * y**i)) for i in range(4)]


def linear_quotient(denominator, numerator, x, y):
    """Return T=t_x*x+t_y*y iff denominator*T=numerator exactly."""
    d = [QQ(denominator.monomial_coefficient(x ** (2 - i) * y**i)) for i in range(3)]
    n = degree_three_coefficients(numerator, x, y)
    equations = matrix(QQ, [[d[0], 0], [d[1], d[0]], [d[2], d[1]], [0, d[2]]])
    try:
        answer = equations.solve_right(vector(QQ, n))
    except ValueError:
        return None
    trial = answer[0] * x + answer[1] * y
    return trial if denominator * trial == numerator else None


def base_pencil(form, px, py):
    """Swap coordinates of the split chart before its direct perturbation."""
    A, B = base.matrix_pair_for_affine_slice(form, px, py)
    swap = matrix(ZZ, [[0, 1, 0], [1, 0, 0], [0, 0, 1]])
    return swap.transpose() * A * swap, swap.transpose() * B * swap


def direct_perturbations(form, parameter_bound, perturbation_bound):
    R = PolynomialRing(QQ, names=("x", "y"))
    x, y = R.gens()
    F = sum(QQ(value) * x ** (3 - i) * y**i for i, value in enumerate(form))
    coefficients = range(-perturbation_bound, perturbation_bound + 1)
    for px, py in itertools.product(range(-parameter_bound, parameter_bound + 1), repeat=2):
        A, B = base_pencil(form, px, py)
        M = x * A - y * B
        L, m, n, r, s = M[0, 0], M[0, 1], M[0, 2], M[1, 1], M[1, 2]
        for rx, ry, sx, sy in itertools.product(coefficients, repeat=4):
            Rdelta, Sdelta = rx * x + ry * y, sx * x + sy * y
            rr, ss = r + Rdelta, s + Sdelta
            denominator = L * rr - m * m
            numerator = F + L * ss * ss - 2 * m * n * ss + n * n * rr
            if denominator == 0:
                continue
            tt = linear_quotient(denominator, numerator, x, y)
            if tt is None:
                continue
            P = matrix(QQ, [[L, m, n], [m, rr, ss], [n, ss, tt]])
            Anew = matrix(QQ, 3, 3, lambda i, j: coeff(P[i, j], x, y, 1))
            Bnew = matrix(QQ, 3, 3, lambda i, j: -coeff(P[i, j], x, y, 0))
            yield (px, py, rx, ry, sx, sy), Anew, Bnew


def compute():
    protocol = read(PROTOCOL)
    preflight, raw, nf = base.maximal_cubic_preflight()
    assert preflight["status"] == "PASS"
    form = [ZZ(value) for value in preflight["binary_cubic_descending"]]
    settings = protocol["bounds"]
    seen, candidates, accepted, attempted = set(), [], [], 0
    for parameters, A, B in direct_perturbations(form, settings["base_affine_parameter_bound"], settings["perturbation_coefficient_bound"]):
        attempted += 1
        coefficients = base.pair_coefficients(A, B)
        key = tuple(coefficients)
        if key in seen:
            continue
        seen.add(key)
        if not base.primitive(coefficients):
            continue
        determinant = base.determinant_coefficients(A, B)
        assert determinant == form
        table, multiplication = base.quartic_multiplication(coefficients)
        probe = base.algebra_probe(multiplication, settings["quartic_generator_l1_bound"])
        row = {"slice_parameters": {"p_x": parameters[0], "p_y": parameters[1], "r_x": parameters[2], "r_y": parameters[3], "s_x": parameters[4], "s_y": parameters[5]},
               "pair_coefficients": [str(value) for value in coefficients],
               "determinant_form_descending": [str(value) for value in determinant],
               "primitive": True, "quartic_ring_associative_integral": True, "algebra_probe": probe}
        if probe["status"] != "FIELD_GENERATOR":
            row["status"] = "REJECTED_ZERO_DIVISOR" if probe["status"] == "ZERO_DIVISOR" else "UNKNOWN_NO_FIELD_WITNESS"
            candidates.append(row)
            continue
        Q = PolynomialRing(QQ, "X")
        polynomial = Q([QQ(value) for value in probe["polynomial_ascending"]])
        structural = base.structural_field_certificate(polynomial, ZZ(preflight["field_discriminant"]), settings["finite_field_witness_prime_bound"])
        row["quartic_generator_polynomial_ascending"] = probe["polynomial_ascending"]
        row["structural"] = structural
        if structural["status"] != "PASS":
            row["status"] = structural["status"]
            candidates.append(row)
            continue
        extraction = base.extract_squareclass(polynomial, raw)
        row["squareclass_extraction"] = {key: ([str(value) for value in value] if key == "alpha_raw_power_basis" and isinstance(value, list) else value) for key, value in extraction.items()}
        if extraction["status"] != "PASS":
            row["status"] = extraction["status"]
            candidates.append(row)
            continue
        strict = base.strict_class_certificate(extraction["alpha_raw_power_basis"], raw, nf, settings["finite_field_witness_prime_bound"])
        row["strict"] = strict
        if strict["status"] != "PASS":
            row["status"] = strict["status"]
            candidates.append(row)
            continue
        row["status"] = "ACCEPTED_STRICT_INDEPENDENT"
        candidates.append(row)
        accepted.append(row)
        row["cover_handoff"] = base.handoff_to_sealed_cover(row, raw, protocol)
        break
    status = "PASS_NEW_STRICT_CHARACTER_HANDED_OFF" if accepted else "BOUNDED_SLICE_NO_ACCEPTED_CHARACTER"
    return {"schema": "elliptic-curves.curve302-fixed-resolvent-orbits.v2", "status": status,
            "protocol": base.relative(PROTOCOL), "protocol_sha256": base.digest(PROTOCOL),
            "bindings": {base.relative(path): base.digest(path) for path in (PROTOCOL, V1, base.MAXIMUM_FORM, base.STRICT, base.FILTRATION, base.CURVE, base.HIDDEN_SOLVER, HERE)},
            "maximal_cubic_preflight": preflight,
            "bounded_search": {"normalization": protocol["normalization"], "declared_determinantal_attempts": settings["determinantal_attempts"],
                                "solved_determinantal_pairs_before_deduplication": attempted, "distinct_primitive_pairs": len(candidates),
                                "raw_pair_candidates": candidates, "accepted_count": len(accepted), "completeness": "NOT_CLAIMED"},
            "next_action": "SEALED_2COVER_POINTSQI" if accepted else "STOP_BOUNDED_SLICE_NO_CLASS_CREATED",
            "boundary": protocol["boundary"]}


def capture():
    if OUTPUT.exists():
        raise FileExistsError(OUTPUT)
    WORK.mkdir(parents=True, exist_ok=True)
    put_new(OUTPUT, compute())
    print(read(OUTPUT)["status"], flush=True)


def check():
    stored = read(OUTPUT)
    assert stored["bindings"] == {base.relative(path): base.digest(path) for path in (PROTOCOL, V1, base.MAXIMUM_FORM, base.STRICT, base.FILTRATION, base.CURVE, base.HIDDEN_SOLVER, HERE)}
    if stored["status"] != "PASS_NEW_STRICT_CHARACTER_HANDED_OFF":
        assert compute() == stored
    print("PASS direct fixed-resolvent determinant slice replay", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("capture", "check"))
    args = parser.parse_args()
    (capture if args.mode == "capture" else check)()
