#!/usr/bin/env sage-python
"""Construct exact rootless R17 bisections from ranked height-ten traces.

For a height-ten section

    tau = (Nx/h^2, Ny/h^3),   deg(h)=3,

the residual line through ``-tau`` has slope ``m=M/h``.  Its unique integral
degree bound is obtained by the linear congruence

    M*Nx + Ny == 0 (mod h^2),   deg(M)<6.

The elliptic equation then forces the chord-discriminant numerator

    M^4 - 6*M^2*Nx - 8*M*Ny - 3*Nx^2 - 4*A*h^4

to be divisible by ``h^6``.  The quotient is quadratic.  This script replays
that construction exactly for a prefix of the equation-priority table and
writes input accepted directly by ``hash_bisection_extensions.py``.
"""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path

from sage.all import EllipticCurve, PolynomialRing, QQ, ZZ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
MODEL = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_model.json"
SECTIONS = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_sections.json"
PINNED = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
DEFAULT_PRIORITY = ROOT / "artifacts/generated-results/elkies-2026-bisection-equation-priority.tsv"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-2026-equation-bisections.json"
DEFAULT_ORBITS = ROOT / "artifacts/generated-results/elkies-2026-equation-bisections-orbits.tsv"


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def display_path(path: Path) -> str:
    try:
        return relative(path)
    except ValueError:
        return str(path.resolve())


def load_matrix(path: Path):
    return matrix(
        ZZ,
        [
            [ZZ(value) for value in line.split()]
            for line in path.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ],
    )


def parse_vector(text: str):
    result = vector(ZZ, [ZZ(value) for value in text.split()])
    if len(result) != 17:
        raise ValueError("expected 17 coordinates")
    return result


def rational_text(value) -> str:
    value = QQ(value)
    return str(value.numerator()) if value.denominator() == 1 else f"{value.numerator()}/{value.denominator()}"


def polynomial_coefficients(polynomial) -> list[str]:
    return [rational_text(polynomial[index]) for index in range(polynomial.degree() + 1)]


def reconstruct_basis(R, A, B, section_data):
    points = []
    for expected_index, record in enumerate(section_data["sections"]):
        assert int(record["basis_index"]) == expected_index
        x_coordinate = R([QQ(value) for value in record["x_coefficients_low_to_high"]])
        if expected_index == 0:
            y_coordinate = R([QQ(value) for value in record["y_coefficients_low_to_high"]])
        else:
            chord = record["chord"]
            reference_x, reference_y = points[int(chord["reference_basis_index"])]
            slope = R([QQ(value) for value in chord["slope_coefficients_low_to_high"]])
            y_coordinate = reference_y + slope * (x_coordinate - reference_x)
        assert y_coordinate**2 == x_coordinate**3 + A * x_coordinate + B
        points.append((x_coordinate, y_coordinate))
    return points


def construct_record(row, *, R, K, E, basis, A, B, Delta, pinned):
    published_vector = parse_vector(row["published_basis_w"])
    pinned_vector = parse_vector(row["pinned_rank17_w"])
    assert pinned_vector * pinned * pinned_vector == 10
    tau = sum((coefficient * point for coefficient, point in zip(published_vector, basis)), E(0))
    if tau.is_zero():
        raise ArithmeticError("height-ten trace unexpectedly vanished")
    X, Y = tau[0], tau[1]
    x_denominator = R(X.denominator())
    if not x_denominator.is_square():
        raise ArithmeticError("trace x denominator is not a square")
    h = R(x_denominator.sqrt())
    h /= h.leading_coefficient()
    assert X.denominator() == h**2
    assert Y.denominator() == h**3
    assert h.degree() == 3
    Nx = R(X * h**2)
    Ny = R(Y * h**3)
    modulus = h**2
    if Nx.gcd(modulus) != 1:
        raise ArithmeticError("trace numerator is not invertible modulo h^2")
    M = R((-Ny * Nx.inverse_mod(modulus)) % modulus)
    assert M.degree() < 6
    assert (M * Nx + Ny) % modulus == 0
    assert (M**2 - Nx) % modulus == 0

    discriminant_numerator = (
        M**4 - 6 * M**2 * Nx - 8 * M * Ny - 3 * Nx**2 - 4 * A * h**4
    )
    quotient, remainder = discriminant_numerator.quo_rem(h**6)
    if remainder:
        raise ArithmeticError("chord discriminant is not divisible by h^6")
    q = R(quotient)
    if q.degree() != 2 or q.gcd(q.derivative()).degree() != 0:
        raise ArithmeticError("residual cover is not a squarefree quadratic")
    branch_fibres_smooth = q.gcd(Delta).degree() == 0

    sum_x = R((M**2 - Nx) // h**2)
    product_x = K(((M * Nx + Ny) ** 2 - B * h**6) / (h**4 * Nx))
    if product_x.denominator() != 1:
        raise ArithmeticError("residual quadratic has a nonintegral constant term")
    product_x = R(product_x)
    assert sum_x**2 - 4 * product_x == h**2 * q

    # Verify the two residual points over u^2=q coefficientwise, without
    # relying on a numerical or specialization check.
    x0 = sum_x / 2
    x1 = h / 2
    n = K(-(Ny + M * Nx) / h**3)
    y0 = K(M / h) * x0 + n
    y1 = K(M / h) * x1
    if K(y0).denominator() != 1 or K(y1).denominator() != 1:
        raise ArithmeticError("residual bisection coordinates are not integral")
    x0, x1, y0, y1 = map(R, (x0, x1, y0, y1))
    assert y0**2 + y1**2 * q == x0**3 + 3 * x0 * x1**2 * q + A * x0 + B
    assert 2 * y0 * y1 == 3 * x0**2 * x1 + x1**3 * q + A * x1

    label = f"orbit-{int(row['orbit_mask'], 0):05x}"
    return {
        "label": label,
        "lattice_orbit_mask": int(row["orbit_mask"], 0),
        "pinned_rank17_w": [int(entry) for entry in pinned_vector],
        "published_basis_w": [int(entry) for entry in published_vector],
        "priority_rank": int(row["priority_rank"]),
        "equation_rank": int(row["equation_rank"]),
        "disjoint_degree_in_priority_pool": int(row["disjoint_degree_in_pool"]),
        "equation_complexity": {
            key: int(row[key])
            for key in (
                "group_addition_upper_bound", "support_count", "dependency_count",
                "coordinate_input_bits", "maximum_absolute_coefficient", "coefficient_l1",
            )
        },
        "trace_section": {
            "h_coefficients": polynomial_coefficients(h),
            "Nx_coefficients": polynomial_coefficients(Nx),
            "Ny_coefficients": polynomial_coefficients(Ny),
        },
        "residual_chord": {
            "slope": "M(t)/h(t)",
            "M_coefficients": polynomial_coefficients(M),
            "linear_congruence": "M*Nx+Ny == 0 mod h^2",
            "discriminant_identity": "sum_x^2-4*product_x=h^2*q",
            "q_coefficients": polynomial_coefficients(q),
            "branch_fibres_smooth": branch_fibres_smooth,
        },
        "quadratic_cover": {
            "leading_coefficients": ["1"],
            "linear_coefficients": polynomial_coefficients(-sum_x),
            "constant_coefficients": polynomial_coefficients(product_x),
        },
        "lifted_section": {
            "cover": "u^2=q(t)",
            "x0_coefficients": polynomial_coefficients(x0),
            "x1_coefficients": polynomial_coefficients(x1),
            "y0_coefficients": polynomial_coefficients(y0),
            "y1_coefficients": polynomial_coefficients(y1),
            "galois_trace": "P(t,u)+P(t,-u)=tau",
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--priority-table", type=Path, default=DEFAULT_PRIORITY)
    parser.add_argument("--limit", type=int, default=100)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--orbits-output", type=Path, default=DEFAULT_ORBITS)
    arguments = parser.parse_args()
    if arguments.limit <= 0:
        parser.error("--limit must be positive")

    model = json.loads(MODEL.read_text())
    section_data = json.loads(SECTIONS.read_text())
    pinned = load_matrix(PINNED)
    R = PolynomialRing(QQ, "t")
    K = R.fraction_field()
    A = R([QQ(value) for value in model["A_coefficients_low_to_high"]])
    B = R([QQ(value) for value in model["B_coefficients_low_to_high"]])
    Delta = R(-16 * (4 * A**3 + 27 * B**2))
    points = reconstruct_basis(R, A, B, section_data)
    E = EllipticCurve(K, [A, B])
    basis = [E(K(x_coordinate), K(y_coordinate)) for x_coordinate, y_coordinate in points]
    with arguments.priority_table.open(newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    selected = rows[: min(arguments.limit, len(rows))]
    if not selected:
        raise ValueError("priority table contains no rows")
    records = [
        construct_record(
            row, R=R, K=K, E=E, basis=basis, A=A, B=B, Delta=Delta, pinned=pinned
        )
        for row in selected
    ]

    arguments.orbits_output.parent.mkdir(parents=True, exist_ok=True)
    with arguments.orbits_output.open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(["orbit_mask", "hex", "pinned_rank17_w"])
        for record in records:
            orbit = record["lattice_orbit_mask"]
            writer.writerow([
                orbit, f"0x{orbit:05x}",
                " ".join(str(value) for value in record["pinned_rank17_w"]),
            ])

    payload = {
        "schema": "elkies-k3.bisection-extension-input.v1",
        "base_parameter": "t",
        "invariant_mw_rank": 17,
        "required_lattice_orbits": {
            "table": display_path(arguments.orbits_output),
            "sha256": digest(arguments.orbits_output),
        },
        "construction": {
            "schema": "elkies-k3.elkies-2026-direct-rootless-bisection-batch.v1",
            "status": "PASS_EXACT_ELKIES_2026_ROOTLESS_BISECTION_BATCH",
            "method": (
                "For tau=(Nx/h^2,Ny/h^3), solve the unique linear congruence "
                "M*Nx+Ny=0 mod h^2 with deg(M)<6, then divide the exact chord "
                "discriminant numerator by h^6."
            ),
            "inputs": {
                relative(path): digest(path)
                for path in (MODEL, SECTIONS, PINNED, arguments.priority_table)
            },
            "record_count": len(records),
            "smooth_branch_fibre_record_count": sum(
                record["residual_chord"]["branch_fibres_smooth"] for record in records
            ),
            "proof_boundary": (
                "Every record is an exact equation-level rational bisection and exact "
                "quadratic extension. A rank-19 conclusion still requires two distinct "
                "orbit records with equal squareclass and an exact anti-invariant height "
                "matrix; no such conclusion is encoded in this input."
            ),
        },
        "bisections": records,
    }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        "ELKIES2026BISECTIONS|records={}|status={}|output={}".format(
            len(records), payload["construction"]["status"], arguments.output
        )
    )


if __name__ == "__main__":
    main()
