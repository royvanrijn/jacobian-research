#!/usr/bin/env sage -python
"""Construct the first deterministic NS-lattice bisection on the 302 parent.

This is the equation-side pilot for the complete degree-two quotient in
``enumerate_curve302_parent_degree2_multisections.py``.  It selects the
minimum-``l1`` exported rational orbit (then its mask), so t=0, curve302's
exceptional points, and all splitting outcomes are absent from selection.

For a minimum vector w of norm ten, use the explicit section P_{-w}.  Since

  P_{-w}.O = 3,
  (2O+4F+phi(w)) + P_{-w} = 3O+9F,

the unique element of H^0(3O+9F) through P_{-w} has a residual bisection.
In a Weierstrass chart its equation is

  f0(t) + f1(t)*x + f2(t)*y = 0,

with degree bounds (9,5,3).  The script verifies the interpolation identity,
removes the known trace section, proves the residual quadratic is generically
irreducible, and factors it at t=0.  It does not infer anything about the
remaining 40,916 lattice bisection orbits.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import QQ, ZZ, PolynomialRing, gcd, lcm, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
LATTICE = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "curve302_parent_degree2_multisection_lattice_v1.json"
)
ORBITS = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "curve302_parent_degree2_multisection_orbits_v1.tsv"
)
PARENT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "curve302_recovered_mw17_parent_v1.json"
)
LOADER = ROOT / "elliptic-curves/cas/load_curve302_recovered_parent.sage"
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "curve302_parent_cheapest_lattice_bisection_v1.json"
)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def qtext(value) -> str:
    value = QQ(value)
    if value.denominator() == 1:
        return str(value.numerator())
    return f"{value.numerator()}/{value.denominator()}"


def polynomial_record(polynomial) -> list[str]:
    polynomial = polynomial.parent()(polynomial)
    if not polynomial:
        return ["0"]
    return [qtext(polynomial[index]) for index in range(polynomial.degree() + 1)]


def function_record(function) -> dict[str, list[str]]:
    return {
        "numerator_coefficients_low_to_high": polynomial_record(function.numerator()),
        "denominator_coefficients_low_to_high": polynomial_record(function.denominator()),
    }


def selected_orbit() -> dict:
    lattice = json.loads(LATTICE.read_text())
    if lattice["status"] != "PASS_COMPLETE_DEGREE2_TRANSLATION_QUOTIENT":
        raise ArithmeticError("the degree-two lattice certificate is unavailable")
    if digest(ORBITS) != lattice["orbits_tsv_sha256"]:
        raise ArithmeticError("the lattice orbit table does not match its certificate")
    rows = []
    for line in ORBITS.read_text().splitlines()[1:]:
        fields = line.split("\t")
        if fields[1] != "rational":
            continue
        word = tuple(int(value) for value in fields[4].split())
        if len(word) != 17 or int(fields[2]) != 10:
            raise ArithmeticError("unexpected rational degree-two orbit row")
        rows.append(
            {
                "orbit_mask": int(fields[0]),
                "minimum_norm": int(fields[2]),
                "parent_MW17_w": word,
            }
        )
    if len(rows) != lattice["rational_bisections"]["translation_orbits"]:
        raise ArithmeticError("incomplete rational orbit table")
    return min(
        rows,
        key=lambda row: (
            sum(abs(value) for value in row["parent_MW17_w"]),
            max(abs(value) for value in row["parent_MW17_w"]),
            row["orbit_mask"],
        ),
    )


def primitive_kernel_relation(point, ring):
    """Return f0,f1,f2 in H0(3O+9F) vanishing on ``point`` exactly."""

    field = point.curve().base_ring()
    t = ring.gen()
    x_coordinate, y_coordinate = point[0], point[1]
    x_denominator = ring(x_coordinate.denominator())
    if not x_denominator.is_square():
        raise ArithmeticError("trace section has a nonsquare x denominator")
    h = ring(x_denominator.sqrt())
    if h.leading_coefficient() < 0:
        h = -h
    numerator_x = ring(x_coordinate * h**2)
    numerator_y = ring(y_coordinate * h**3)
    if (
        h.degree(),
        numerator_x.degree(),
        numerator_y.degree(),
    ) != (3, 10, 15):
        raise ArithmeticError("the selected norm-ten trace has the wrong pole fingerprint")
    # 1, x and y have fibrewise pole orders 0,2,3 and vertical weights 0,4,6.
    # After multiplication by h^3, all twenty columns have degree at most18.
    bounds_and_terms = ((9, h**3), (5, numerator_x * h), (3, numerator_y))
    columns = []
    for bound, term in bounds_and_terms:
        for power in range(bound + 1):
            polynomial = ring(term * t**power)
            columns.append([polynomial[index] for index in range(19)])
    interpolation = matrix(QQ, 19, 20, lambda row, column: columns[column][row])
    if interpolation.rank() != 19 or interpolation.right_kernel().dimension() != 1:
        raise ArithmeticError("H0(3O+9F) interpolation is not a line")
    raw = interpolation.right_kernel().basis()[0]
    raw *= lcm(value.denominator() for value in raw)
    raw = vector(ZZ, raw)
    common = gcd(list(raw))
    if not common:
        raise ArithmeticError("zero interpolation relation")
    raw = vector(ZZ, [value // common for value in raw])
    first = next(value for value in raw if value)
    if first < 0:
        raw = -raw
    pieces = []
    offset = 0
    for bound in (9, 5, 3):
        pieces.append(ring(list(raw[offset : offset + bound + 1])))
        offset += bound + 1
    f0, f1, f2 = pieces
    if ring.gcd(ring.gcd(f0, f1), f2).degree() != 0:
        raise ArithmeticError("the linear-system member has a vertical fibre factor")
    if f0 * h**3 + f1 * numerator_x * h + f2 * numerator_y:
        raise ArithmeticError("interpolation relation does not contain the trace")
    return {
        "f0": f0,
        "f1": f1,
        "f2": f2,
        "h": h,
        "numerator_x": numerator_x,
        "numerator_y": numerator_y,
        "rank": interpolation.rank(),
        "kernel_coefficients": raw,
    }


def residual_quadratic(curve, point, relation):
    """Eliminate y and divide by the known trace section over QQ(t)."""

    field = curve.base_ring()
    x_ring = PolynomialRing(field, "x")
    x = x_ring.gen()
    f0, f1, f2 = (field(relation[key]) for key in ("f0", "f1", "f2"))
    u = f0 + f1 * x
    a1, a2, a3, a4, a6 = curve.a_invariants()
    elimination = (
        u**2
        - a1 * x * u * f2
        - a3 * u * f2
        - f2**2 * (x**3 + a2 * x**2 + a4 * x + a6)
    )
    residual, remainder = elimination.quo_rem(x - field(point[0]))
    if remainder or residual.degree() != 2:
        raise ArithmeticError("trace removal did not leave a quadratic cover")
    factors = [
        int(factor.degree())
        for factor, exponent in residual.factor()
        for unused in range(int(exponent))
    ]
    if factors != [2]:
        raise ArithmeticError("the residual bisection is generically reducible")
    return residual


def build() -> dict:
    import runpy

    selected = selected_orbit()
    parent = json.loads(PARENT.read_text())
    loader = runpy.run_path(str(LOADER))
    curve, basis, t0 = loader["load_curve302_recovered_parent"]()
    if t0 != 0 or len(basis) != 17:
        raise ArithmeticError("the literal curve302 specialization changed")
    w = vector(ZZ, selected["parent_MW17_w"])
    gram = matrix(ZZ, parent["generic_height_gram"])
    if int(w * gram * w) != 10:
        raise ArithmeticError("selected vector lost its norm-ten attachment")
    trace = -sum((coefficient * section for coefficient, section in zip(w, basis)), curve(0))
    if trace.is_zero():
        raise ArithmeticError("the selected trace vanished")
    ring = curve.base_ring().ring()
    relation = primitive_kernel_relation(trace, ring)
    residual = residual_quadratic(curve, trace, relation)
    zero_ring = PolynomialRing(QQ, "x")
    zero_residual = zero_ring([QQ(coefficient(t0)) for coefficient in residual.list()])
    if zero_residual.degree() != 2 or not zero_residual.discriminant():
        raise ArithmeticError("the t=0 residual is not an etale quadratic")
    zero_factor_degrees = [
        int(factor.degree())
        for factor, exponent in zero_residual.factor()
        for unused in range(int(exponent))
    ]
    split = zero_factor_degrees == [1, 1]
    if split:
        # This pilot is selection-blind but no longer quotient-blind after a
        # split.  A future splitter must attach every root to D/sp(M); do not
        # serialize a guessed quotient coordinate here.
        raise ArithmeticError("unexpected split: install the quotient-attachment gate first")
    if zero_factor_degrees != [2]:
        raise ArithmeticError("unexpected t=0 factor pattern")
    return {
        "schema": "elliptic-curves.curve302-parent-cheapest-lattice-bisection.v1",
        "status": "PASS_EXACT_NONSPLIT_CHEAPEST_LATTICE_BISECTION",
        "scope": (
            "One deterministic minimum-l1 orbit of the complete degree-two lattice "
            "quotient. The selection uses only exported lattice coordinates. Its residual "
            "bisection is constructed and generically irreducible over QQ(t), but is "
            "nonsplit at t=0. This is one exact negative control, not a statement about "
            "the other 40916 rational bisection orbits."
        ),
        "inputs": {
            relative(path): digest(path)
            for path in (LATTICE, ORBITS, PARENT, LOADER, Path(__file__))
        },
        "selection": {
            "rule": "lexicographically minimize (parent-MW17 l1, parent-MW17 linfinity, orbit mask)",
            "orbit_mask": selected["orbit_mask"],
            "minimum_norm": selected["minimum_norm"],
            "bisection_phi_w_parent_MW17_coordinates": list(map(int, w)),
            "trace_section_phi_coordinates": list(map(int, -w)),
            "selection_uses_t0_or_exceptional_points": False,
        },
        "RR_construction": {
            "identity": "(2O+4F+phi(w))+P_(-w)=3O+9F",
            "ambient": "H0(3O+9F)=f0(t)+f1(t)x+f2(t)y, degree bounds (9,5,3)",
            "interpolation_rank": relation["rank"],
            "trace_pole_data": {
                "h_degree": int(relation["h"].degree()),
                "numerator_x_degree": int(relation["numerator_x"].degree()),
                "numerator_y_degree": int(relation["numerator_y"].degree()),
            },
            "f0_coefficients_low_to_high": polynomial_record(relation["f0"]),
            "f1_coefficients_low_to_high": polynomial_record(relation["f1"]),
            "f2_coefficients_low_to_high": polynomial_record(relation["f2"]),
            "primitive_kernel_maximum_coefficient_bits": max(
                abs(int(value)).bit_length() for value in relation["kernel_coefficients"]
            ),
            "generic_residual_x_coefficients_low_to_high": [
                function_record(coefficient) for coefficient in residual.list()
            ],
            "generic_residual_factor_degrees": [2],
        },
        "zero_fibre": {
            "parameter": "0",
            "residual_x_polynomial_coefficients_low_to_high": polynomial_record(zero_residual),
            "factor_degrees": zero_factor_degrees,
            "discriminant": qtext(zero_residual.discriminant()),
            "split_over_Q": split,
            "rational_bisection_fibre_points": [],
            "quotient_attachment": "VACUOUS_NONSPLIT",
        },
        "boundary": (
            "No missing displayed quotient direction is represented in this nonsplit fibre. "
            "The result supplies the uniform RR/elimination interface and one exact control; "
            "it does not rank, discard, or sample the remaining lattice quotient."
        ),
        "reproducing_command": (
            "sage -python elliptic-curves/cas/"
            "construct_curve302_parent_cheapest_lattice_bisection.sage --check"
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    payload = build()
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    if arguments.check:
        if json.loads(arguments.output.read_text()) != payload:
            raise ArithmeticError("stored cheapest-bisection certificate does not replay")
    else:
        if arguments.output.exists():
            raise FileExistsError("refusing to overwrite the generated certificate")
        arguments.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        "CURVE302D2EQUATION|orbit={}|generic_degree=2|t0_split=0|"
        "status=PASS_EXACT_NONSPLIT_CHEAPEST_LATTICE_BISECTION".format(
            payload["selection"]["orbit_mask"]
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
