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
that construction exactly for an interval of the equation-priority table and
writes input accepted directly by ``hash_bisection_extensions.py``.  If a
trace has a pole over infinity, the same calculation is performed after
``t -> 1/t`` and the resulting relation is transported back exactly.
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


def reciprocal_with_bound(polynomial, bound, R):
    """Return ``t^bound*polynomial(1/t)`` in the current polynomial ring."""

    polynomial = R(polynomial)
    if polynomial.degree() > bound:
        raise ArithmeticError("reciprocal degree exceeds its geometric bound")
    t = R.gen()
    return R(sum(polynomial[index] * t ** (bound - index) for index in range(polynomial.degree() + 1)))


def invert_rational(function, weight, R, K):
    """Return ``t^weight*function(1/t)`` without symbolic substitution."""

    function = K(function)
    numerator = R(function.numerator())
    denominator = R(function.denominator())
    t = R.gen()
    reverse_numerator = reciprocal_with_bound(numerator, numerator.degree(), R)
    reverse_denominator = reciprocal_with_bound(denominator, denominator.degree(), R)
    exponent = weight - numerator.degree() + denominator.degree()
    return K(t**exponent * reverse_numerator / reverse_denominator)


def trace_chord_frame(X, Y, R):
    """Return the denominator frame and least regular residual-chord slope.

    The congruence is independent of the genus of the residual bisection.  A
    rational bisection uses the least representative ``M0``.  For a genus-one
    bisection pencil the regular representatives are ``M0 + lambda*h^2``.
    """
    x_denominator = R(X.denominator())
    if not x_denominator.is_square():
        raise ArithmeticError("trace x denominator is not a square")
    h = R(x_denominator.sqrt())
    h /= h.leading_coefficient()
    assert X.denominator() == h**2
    assert Y.denominator() == h**3
    Nx = R(X * h**2)
    Ny = R(Y * h**3)
    modulus = h**2
    if Nx.gcd(modulus) != 1:
        raise ArithmeticError("trace numerator is not invertible modulo h^2")
    M0 = R((-Ny * Nx.inverse_mod(modulus)) % modulus)
    assert M0.degree() < 2 * h.degree()
    assert (M0 * Nx + Ny) % modulus == 0
    assert (M0**2 - Nx) % modulus == 0
    assert M0.gcd(h).degree() == 0
    return {"h": h, "Nx": Nx, "Ny": Ny, "M0": M0}


def chord_data_from_slope_numerator(h, Nx, Ny, M, A, B, Delta, R, K, expected_q_degree):
    """Compile and verify a residual chord from a regular slope numerator."""

    modulus = h**2
    M = R(M)
    if (M * Nx + Ny) % modulus:
        raise ArithmeticError("slope numerator is not regular modulo h^2")

    discriminant_numerator = (
        M**4 - 6 * M**2 * Nx - 8 * M * Ny - 3 * Nx**2 - 4 * A * h**4
    )
    quotient, remainder = discriminant_numerator.quo_rem(h**6)
    if remainder:
        raise ArithmeticError("chord discriminant is not divisible by h^6")
    q = R(quotient)
    if q.degree() != expected_q_degree or q.gcd(q.derivative()).degree() != 0:
        raise ArithmeticError(
            f"residual cover is not squarefree of degree {expected_q_degree}"
        )

    sum_x = R((M**2 - Nx) // h**2)
    product_x = K(((M * Nx + Ny) ** 2 - B * h**6) / (h**4 * Nx))
    if product_x.denominator() != 1:
        raise ArithmeticError("residual quadratic has a nonintegral constant term")
    product_x = R(product_x)
    assert sum_x**2 - 4 * product_x == h**2 * q

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
    return {
        "h": h, "Nx": Nx, "Ny": Ny, "M": M, "q": q,
        "sum_x": sum_x, "product_x": product_x,
        "x0": x0, "x1": x1, "y0": y0, "y1": y1,
        "branch_fibres_smooth": q.gcd(Delta).degree() == 0,
    }


def local_chord_data(X, Y, A, B, Delta, R, K):
    """Construct the rational chord where all three trace poles are finite."""

    frame = trace_chord_frame(X, Y, R)
    if frame["h"].degree() != 3:
        raise ArithmeticError(
            f"trace denominator has degree {frame['h'].degree()}, not three"
        )
    return chord_data_from_slope_numerator(
        frame["h"], frame["Nx"], frame["Ny"], frame["M0"],
        A, B, Delta, R, K, expected_q_degree=2,
    )


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
    chart = "finite"
    if h.degree() == 3:
        data = local_chord_data(X, Y, A, B, Delta, R, K)
    else:
        # A degree drop means that one or more intersections with the zero
        # section lie over infinity.  Move infinity to zero, carry out the
        # identical resolved Riemann--Roch construction, then transport the
        # bisection and its quadratic cover back to the published t-chart.
        chart = "inverted_at_infinity"
        A_inverse = reciprocal_with_bound(A, 8, R)
        B_inverse = reciprocal_with_bound(B, 12, R)
        Delta_inverse = reciprocal_with_bound(Delta, 24, R)
        X_inverse = invert_rational(X, 4, R, K)
        Y_inverse = invert_rational(Y, 6, R, K)
        inverse = local_chord_data(
            X_inverse, Y_inverse, A_inverse, B_inverse, Delta_inverse, R, K
        )
        data = {
            "h": reciprocal_with_bound(inverse["h"], 3, R),
            "M": reciprocal_with_bound(inverse["M"], 5, R),
            "q": reciprocal_with_bound(inverse["q"], 2, R),
            "sum_x": reciprocal_with_bound(inverse["sum_x"], 4, R),
            "product_x": reciprocal_with_bound(inverse["product_x"], 8, R),
            "x0": reciprocal_with_bound(inverse["x0"], 4, R),
            "x1": reciprocal_with_bound(inverse["x1"], 3, R),
            "y0": reciprocal_with_bound(inverse["y0"], 6, R),
            "y1": reciprocal_with_bound(inverse["y1"], 5, R),
        }
        data["Nx"] = R(X * data["h"]**2)
        data["Ny"] = R(Y * data["h"]**3)
        data["branch_fibres_smooth"] = data["q"].gcd(Delta).degree() == 0
        assert data["sum_x"]**2 - 4 * data["product_x"] == data["h"]**2 * data["q"]
        assert data["x0"] == data["sum_x"] / 2
        assert data["x1"] == data["h"] / 2
        assert data["y0"]**2 + data["y1"]**2 * data["q"] == (
            data["x0"]**3 + 3 * data["x0"] * data["x1"]**2 * data["q"]
            + A * data["x0"] + B
        )
        assert 2 * data["y0"] * data["y1"] == (
            3 * data["x0"]**2 * data["x1"] + data["x1"]**3 * data["q"]
            + A * data["x1"]
        )
    h, Nx, Ny, M, q = (data[key] for key in ("h", "Nx", "Ny", "M", "q"))
    sum_x, product_x = data["sum_x"], data["product_x"]
    x0, x1, y0, y1 = (data[key] for key in ("x0", "x1", "y0", "y1"))
    branch_fibres_smooth = data["branch_fibres_smooth"]

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
            "construction_chart": chart,
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
    parser.add_argument(
        "--start", type=int, default=0,
        help="zero-based first row of the priority table to construct",
    )
    parser.add_argument("--limit", type=int, default=100)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--orbits-output", type=Path, default=DEFAULT_ORBITS)
    arguments = parser.parse_args()
    if arguments.start < 0:
        parser.error("--start must be nonnegative")
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
    assert Delta.degree() == 24
    assert Delta.gcd(Delta.derivative()).degree() == 0
    points = reconstruct_basis(R, A, B, section_data)
    E = EllipticCurve(K, [A, B])
    basis = [E(K(x_coordinate), K(y_coordinate)) for x_coordinate, y_coordinate in points]
    with arguments.priority_table.open(newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    selected = rows[arguments.start : arguments.start + arguments.limit]
    if not selected:
        raise ValueError("selected priority-table interval contains no rows")
    records = []
    for row in selected:
        try:
            records.append(
                construct_record(
                    row, R=R, K=K, E=E, basis=basis, A=A, B=B, Delta=Delta, pinned=pinned
                )
            )
        except Exception as error:
            raise RuntimeError(
                "failed priority_rank={} equation_rank={} orbit={}".format(
                    row["priority_rank"], row["equation_rank"], row["orbit_hex"]
                )
            ) from error

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
                "discriminant numerator by h^6; use the reciprocal base chart "
                "when a pole of tau lies over infinity."
            ),
            "inputs": {
                relative(path): digest(path)
                for path in (MODEL, SECTIONS, PINNED, arguments.priority_table)
            },
            "record_count": len(records),
            "smooth_branch_fibre_record_count": sum(
                record["residual_chord"]["branch_fibres_smooth"] for record in records
            ),
            "construction_chart_counts": {
                chart: sum(
                    record["residual_chord"]["construction_chart"] == chart
                    for record in records
                )
                for chart in ("finite", "inverted_at_infinity")
            },
            "proof_boundary": (
                "Every record is an exact equation-level rational bisection and exact "
                "quadratic extension. A rank-19 conclusion still requires two distinct "
                "orbit records with equal squareclass and an exact anti-invariant height "
                "matrix; no such conclusion is encoded in this input."
            ),
        },
        "individual_base_change_certificate": {
            "status": "PASS_EACH_COVER_GENERIC_RANK_AT_LEAST_18",
            "record_count": len(records),
            "invariant_mw_rank": 17,
            "source_fibres": "24I1",
            "source_discriminant_degree": int(Delta.degree()),
            "source_discriminant_squarefree": True,
            "branch_fibres_smooth_for_every_record": all(
                record["residual_chord"]["branch_fibres_smooth"] for record in records
            ),
            "base_change_chi": 4,
            "base_change_rootless": True,
            "lift_self_intersection": -4,
            "conjugate_lift_intersection": 2,
            "anti_invariant_height": 12,
            "height_formula": "2*(P.sigma(P)-P.P)=2*(2-(-4))=12",
            "generic_mw_rank_lower_bound": 18,
            "proof": (
                "The squarefree degree-24 discriminant gives 24I1 fibres. Each squarefree "
                "quadratic q is coprime to that discriminant, so the degree-two pullback "
                "branches only at two smooth fibres and remains rootless with chi=4. "
                "The conjugate lifted sections meet transversely at exactly those two "
                "branch points; their difference therefore has height 12 and is non-torsion."
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
