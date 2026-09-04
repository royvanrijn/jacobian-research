#!/usr/bin/env sage
"""Eliminate the first jets of a literal five-fibre rootless-K3 interpolation.

The input consists of five canonical short fibres and seventeen labelled
points on each fibre.  After the weighted PGL2 and Weierstrass gauges

    t_351=0, t_356=1, t_385=-1, u_351=1,

the script exhausts the remaining ordered pair (t_376,t_377) over P1(GF(p)),
including the two boundary charts in which exactly one residual node is at
infinity.
For each pair it eliminates the 34 free ordinate-interpolation coefficients
from the differentiated section identities and asks msolve whether the
resulting first-jet system has a solution over the algebraic closure.

An empty finite-field chart is a necessary-condition obstruction only.  It
does not exclude a rational model whose parameters or scalings have bad or
colliding reduction at the selected prime, and it assumes literal point
labels and signs rather than merely a nearby Mordell--Weil basis.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
import hashlib
import json
from multiprocessing import get_context
from pathlib import Path
import shutil
import subprocess
import sys
import time

from sage.all import GF, PolynomialRing, matrix, prod, vector


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "icarm_wgxli_rank17_lineage_v1.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "icarm_wgxli_rank17_first_jet_mod17_v1.json"
)
DEFAULT_WORK = ROOT / "artifacts/local/elliptic-curves/wgxli-r17-first-jet"
FIBRE_IDS = (351, 356, 376, 377, 385)
FIXED_NODE_IDS = {351: 0, 356: 1, 385: -1}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--work-dir", type=Path, default=DEFAULT_WORK)
    parser.add_argument("--prime", type=int, default=17)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument(
        "--jobs",
        type=int,
        default=1,
        help="number of independent parameter charts solved concurrently",
    )
    parser.add_argument("--pair-timeout", type=float, default=10.0)
    parser.add_argument(
        "--reuse-unit-outputs",
        action="store_true",
        help=(
            "reuse an existing msolve output only when it is exactly the unit ideal; "
            "all missing, partial, timeout, or nonunit outputs are solved again"
        ),
    )
    parser.add_argument("--max-pairs", type=int)
    parser.add_argument("--check", action="store_true")
    return parser.parse_args()


def sha256_bytes(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def sha256_file(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def mod_fraction(field, value: str | int):
    rational = Fraction(value)
    if rational.denominator % int(field.characteristic()) == 0:
        raise ZeroDivisionError(f"coordinate denominator vanishes modulo {field.characteristic()}")
    return field(rational.numerator) / field(rational.denominator)


def msolve_version() -> str:
    completed = subprocess.run(
        ["msolve", "-h"], text=True, stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT, check=False,
    )
    first = completed.stdout.splitlines()[0] if completed.stdout.splitlines() else "msolve"
    return first.strip()


def load_fibres(path: Path, prime: int):
    raw = path.read_bytes()
    payload = json.loads(raw)
    rows = payload["rootless_k3_interpolation_input"]["fibres"]
    by_id = {int(row["curve_id"]): row for row in rows}
    if tuple(by_id) != FIBRE_IDS:
        raise AssertionError("unexpected five-fibre order")
    admissible = payload["rootless_k3_interpolation_input"][
        "admissible_primes_below_300"
    ]
    if prime < 300 and prime not in admissible:
        raise ValueError(
            f"prime {prime} is not admissible for all five fibres and 85 coordinates"
        )
    return raw, by_id


def node_coordinates(field, node):
    """Return the fixed homogeneous representative of a point of P1."""
    return (field(1), field(0)) if node is None else (field(node), field(1))


def evaluation_row(field, degree: int, node):
    a_value, b_value = node_coordinates(field, node)
    return vector(
        field,
        [a_value**index * b_value ** (degree - index) for index in range(degree + 1)],
    )


def derivative_row(field, degree: int, node):
    """Differentiate in t on [t:1], and in s on [1:s] at infinity."""
    if node is None:
        return vector(
            field,
            [field(index == degree - 1) for index in range(degree + 1)],
        )
    t_value = field(node)
    return vector(
        field,
        [
            field(0) if index == 0 else field(index) * t_value ** (index - 1)
            for index in range(degree + 1)
        ],
    )


def interpolate_binary_form(ring, field, nodes, values, degree: int):
    """Choose one homogeneous degree-d form having the displayed values."""
    evaluation = matrix(field, [evaluation_row(field, degree, node) for node in nodes])
    if evaluation.rank() != 5:
        raise ValueError("the projective interpolation nodes must be distinct")
    pivots = evaluation.pivots()
    if len(pivots) != 5:
        raise AssertionError("binary-form evaluation did not have five pivots")
    square = evaluation.matrix_from_columns(pivots)
    coefficients_at_pivots = square.change_ring(ring).solve_right(vector(ring, values))
    coefficients = [ring.zero()] * (degree + 1)
    for index, coefficient in zip(pivots, coefficients_at_pivots):
        coefficients[index] = coefficient
    return vector(ring, coefficients)


def node_polynomial_derivative(field, nodes, index: int):
    """Derivative of prod(T-tZ), including the Z factor for infinity."""
    node = nodes[index]
    a_value, b_value = node_coordinates(field, node)
    answer = field(1)
    for other, other_node in enumerate(nodes):
        if other == index:
            continue
        if other_node is None:
            answer *= b_value
        else:
            answer *= a_value - field(other_node) * b_value
    return answer


def build_pair_system(
    fibres, prime: int, node_376, node_377, require_generic_equation_count=True
):
    field = GF(prime)
    names = (
        ("w356", "w376", "w377", "w385")
        + tuple(f"alpha{index}" for index in range(5))
        + tuple(f"beta{index}" for index in range(5))
        + ("inverse",)
    )
    ring = PolynomialRing(field, names=names, order="degrevlex")
    variables = ring.gens_dict()
    weights = [
        ring(1),
        variables["w356"],
        variables["w376"],
        variables["w377"],
        variables["w385"],
    ]
    alphas = [variables[f"alpha{index}"] for index in range(5)]
    betas = [variables[f"beta{index}"] for index in range(5)]
    inverse = variables["inverse"]
    nodes = [0, 1, node_376, node_377, -1]
    coordinates = [node_coordinates(field, node) for node in nodes]
    if len(set(coordinates)) != 5:
        raise ValueError("the normalized nodes must be distinct")
    node_derivatives = [
        ring(node_polynomial_derivative(field, nodes, index)) for index in range(5)
    ]

    family_A_values = [
        ring(mod_fraction(field, fibres[curve_id]["short_model"][0])) * weight**4
        for curve_id, weight in zip(FIBRE_IDS, weights)
    ]
    equations = []
    for position in range(17):
        x_values = []
        y_values = []
        for curve_id, weight in zip(FIBRE_IDS, weights):
            x_value, y_value = fibres[curve_id]["short_points_first_17"][position]
            x_values.append(ring(mod_fraction(field, x_value)) * weight**2)
            y_values.append(ring(mod_fraction(field, y_value)) * weight**3)
        x_interpolant = interpolate_binary_form(ring, field, nodes, x_values, 4)
        y_interpolant = interpolate_binary_form(ring, field, nodes, y_values, 6)

        ordinate_coefficients = []
        derivative_right_sides = []
        for index, node in enumerate(nodes):
            x_value = x_values[index]
            y_value = y_values[index]
            x_derivative = ring(
                derivative_row(field, 4, node) * x_interpolant
            )
            y_derivative = ring(
                derivative_row(field, 6, node) * y_interpolant
            )
            ordinate_coefficients.append(2 * y_value * node_derivatives[index])
            known_part = (
                2 * y_value * y_derivative
                - 3 * x_value**2 * x_derivative
                - family_A_values[index] * x_derivative
            )
            derivative_right_sides.append(
                x_value * alphas[index] + betas[index] - known_part
            )

        # The quotient
        #   (x*alpha+beta-known)/(2*y*L')
        # must be affine in the base parameter because every sextic ordinate
        # differs from its five-point quartic interpolant by L(t)*(r+s*t).
        for index in range(2, 5):
            a_value, b_value = coordinates[index]
            a_value = ring(a_value)
            b_value = ring(b_value)
            equations.append(
                ordinate_coefficients[0]
                * ordinate_coefficients[1]
                * derivative_right_sides[index]
                - ordinate_coefficients[index]
                * (
                    (b_value - a_value)
                    * ordinate_coefficients[1]
                    * derivative_right_sides[0]
                    + a_value
                    * ordinate_coefficients[0]
                    * derivative_right_sides[1]
                )
            )

    # A binary octic has nine coefficients.  Its five values and five local
    # first derivatives therefore satisfy the unique left-kernel relation of
    # this 10-by-9 evaluation matrix.  This is the projective replacement for
    # killing the t^9 coefficient of the affine Hermite interpolant.
    value_rows = [evaluation_row(field, 8, node) for node in nodes]
    derivative_rows = [derivative_row(field, 8, node) for node in nodes]
    hermite_matrix = matrix(field, value_rows + derivative_rows)
    hermite_relations = hermite_matrix.left_kernel().basis()
    if len(hermite_relations) != 1:
        raise AssertionError("binary-octic first jets did not have one relation")
    relation = hermite_relations[0]
    equations.append(
        sum(ring(relation[index]) * family_A_values[index] for index in range(5))
        + sum(ring(relation[index + 5]) * alphas[index] for index in range(5))
    )
    equations.append(inverse * prod(weights[1:]) - 1)
    equations = [equation for equation in equations if equation]
    # There are nominally 53 equations.  At p=17 the binary-octic first-jet
    # relation vanishes on the literal input; retain the resulting 52 nonzero
    # equations.
    if require_generic_equation_count and len(equations) not in (52, 53):
        raise AssertionError(
            f"expected 52 or 53 nonzero first-jet equations, got {len(equations)}"
        )
    return names, equations


def render_msolve(path: Path, names, equations, prime: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    # The header lines are not comma-terminated; polynomial lines are.
    with path.open("w") as handle:
        handle.write(",".join(names) + "\n")
        handle.write(str(prime) + "\n")
        for index, equation in enumerate(equations):
            handle.write(str(equation).replace("**", "^"))
            handle.write(",\n" if index + 1 < len(equations) else "\n")


def solve_chart(
    msolve_input: Path,
    msolve_output: Path,
    threads: int,
    pair_timeout: float,
) -> dict[str, object]:
    """Solve one already-rendered chart and return its deterministic result."""
    pair_started = time.monotonic()
    try:
        completed = subprocess.run(
            [
                "msolve", "-t", str(threads), "-v", "0",
                "-f", str(msolve_input), "-o", str(msolve_output),
            ],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            timeout=pair_timeout,
            check=False,
        )
        output_text = msolve_output.read_text(errors="replace").strip()
        no_solution = output_text in ("[-1]", "[-1]:")
        status = (
            "NO_GEOMETRIC_SOLUTION"
            if no_solution
            else "SOLUTION_OR_NONZERO_DIMENSION"
        )
        returncode = completed.returncode
    except subprocess.TimeoutExpired:
        output_text = ""
        status = "TIMEOUT"
        returncode = None
    return {
        "status": status,
        "msolve_returncode": returncode,
        "output_prefix": output_text[:200],
        "runtime_seconds": time.monotonic() - pair_started,
    }


_WORKER_FIBRES = None
_WORKER_PRIME = None
_WORKER_WORK = None
_WORKER_THREADS = None
_WORKER_TIMEOUT = None
_WORKER_REUSE_UNIT_OUTPUTS = None
_WORKER_TRUSTED_INPUT_HASHES = None


def initialize_chart_worker(
    fibres, prime, work, threads, pair_timeout, reuse_unit_outputs,
    trusted_input_hashes,
) -> None:
    """Install read-only run data once in each chart worker."""
    global _WORKER_FIBRES, _WORKER_PRIME, _WORKER_WORK
    global _WORKER_THREADS, _WORKER_TIMEOUT, _WORKER_REUSE_UNIT_OUTPUTS
    global _WORKER_TRUSTED_INPUT_HASHES
    _WORKER_FIBRES = fibres
    _WORKER_PRIME = prime
    _WORKER_WORK = work
    _WORKER_THREADS = threads
    _WORKER_TIMEOUT = pair_timeout
    _WORKER_REUSE_UNIT_OUTPUTS = reuse_unit_outputs
    _WORKER_TRUSTED_INPUT_HASHES = trusted_input_hashes


def build_and_solve_chart(pair) -> dict[str, object]:
    """Construct and solve one chart in an independent worker process."""
    node_376, node_377 = pair
    names, equations = build_pair_system(
        _WORKER_FIBRES, _WORKER_PRIME, node_376, node_377
    )
    left_text = "inf" if node_376 is None else str(node_376)
    right_text = "inf" if node_377 is None else str(node_377)
    stem = f"a{left_text}-b{right_text}"
    msolve_input = _WORKER_WORK / f"{stem}.ms"
    msolve_output = _WORKER_WORK / f"{stem}.solve"
    render_msolve(msolve_input, names, equations, _WORKER_PRIME)
    input_hash = sha256_file(msolve_input)
    cached_output = (
        msolve_output.read_text(errors="replace").strip()
        if (
            _WORKER_REUSE_UNIT_OUTPUTS
            and _WORKER_TRUSTED_INPUT_HASHES.get((node_376, node_377)) == input_hash
            and msolve_output.exists()
        )
        else ""
    )
    if cached_output in ("[-1]", "[-1]:"):
        result = {
            "status": "NO_GEOMETRIC_SOLUTION",
            "msolve_returncode": 0,
            "output_prefix": cached_output,
            "runtime_seconds": 0.0,
        }
    else:
        result = solve_chart(
            msolve_input, msolve_output, _WORKER_THREADS, _WORKER_TIMEOUT
        )
    result["nodes_376_377"] = [node_376, node_377]
    result["msolve_input_sha256"] = input_hash
    result["equation_count"] = len(equations)
    return result


def positive_control(prime: int) -> dict[str, object]:
    """Check the eliminated identities on the exact published R17 model."""
    field = GF(prime)
    model_path = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_model.json"
    sections_path = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_sections.json"
    if not model_path.exists() or not sections_path.exists():
        return {"available": False}
    model = json.loads(model_path.read_text())
    sections = json.loads(sections_path.read_text())
    try:
        A = vector(field, [
            mod_fraction(field, value)
            for value in model["A_coefficients_low_to_high"]
        ])
        B = vector(field, [
            mod_fraction(field, value)
            for value in model["B_coefficients_low_to_high"]
        ])
    except ZeroDivisionError:
        return {
            "available": False,
            "reason": "published-R17 control model has a coefficient denominator at this prime",
        }
    control_polynomial_ring = PolynomialRing(field, "q")
    reconstructed_sections = []
    for row in sections["sections"]:
        try:
            X_polynomial = control_polynomial_ring([
                mod_fraction(field, value)
                for value in row["x_coefficients_low_to_high"]
            ])
            if not reconstructed_sections:
                Y_polynomial = control_polynomial_ring([
                    mod_fraction(field, value)
                    for value in row["y_coefficients_low_to_high"]
                ])
            else:
                chord = row["chord"]
                reference_x, reference_y = reconstructed_sections[
                    int(chord["reference_basis_index"])
                ]
                slope = control_polynomial_ring([
                    mod_fraction(field, value)
                    for value in chord["slope_coefficients_low_to_high"]
                ])
                Y_polynomial = reference_y + slope * (X_polynomial - reference_x)
        except ZeroDivisionError:
            return {
                "available": False,
                "reason": "published-R17 control section has a coefficient denominator at this prime",
            }
        reconstructed_sections.append((X_polynomial, Y_polynomial))
    section_coefficients = [
        (
            vector(field, list(X) + [field(0)] * (5 - len(list(X)))),
            vector(field, list(Y) + [field(0)] * (7 - len(list(Y)))),
        )
        for X, Y in reconstructed_sections
    ]

    allowed_control_nodes = [
        value for value in range(prime) if field(value) not in (0, 1, -1)
    ]
    control_chart_candidates = (
        (
            (left, right)
            for left in allowed_control_nodes
            for right in allowed_control_nodes
            if left != right
        ),
        ((None, right) for right in allowed_control_nodes),
        ((left, None) for left in allowed_control_nodes),
    )
    chart_records = []
    for candidates in control_chart_candidates:
        selected = None
        degenerate_fallback = None
        for node_376, node_377 in candidates:
            nodes = [0, 1, node_376, node_377, -1]
            synthetic = {}
            for curve_id, node in zip(FIBRE_IDS, nodes):
                a_value = evaluation_row(field, 8, node) * A
                b_value = evaluation_row(field, 12, node) * B
                points = []
                for X, Y in section_coefficients:
                    x_value = evaluation_row(field, 4, node) * X
                    y_value = evaluation_row(field, 6, node) * Y
                    if y_value**2 != x_value**3 + a_value * x_value + b_value:
                        raise AssertionError("published-R17 positive-control section failed")
                    points.append((str(int(x_value)), str(int(y_value))))
                synthetic[curve_id] = {
                    "short_model": [str(int(a_value)), str(int(b_value))],
                    "short_points_first_17": points,
                }
            # At p=17 the fixed control node q=1 is singular and several
            # reduced equations vanish.  It is still a valid degenerate
            # algebraic witness for the eliminated identities.  The theorem
            # uses the nonsingular 53- and 67-controls; do not silently drop
            # the deliberately recorded mod-17 fallback.
            names, equations = build_pair_system(
                synthetic,
                prime,
                node_376,
                node_377,
                require_generic_equation_count=False,
            )
            assignments = {name: field(1) for name in names}
            for index, node in enumerate(nodes):
                assignments[f"alpha{index}"] = derivative_row(field, 8, node) * A
                assignments[f"beta{index}"] = derivative_row(field, 12, node) * B
            ring = equations[0].parent()
            substitution = {ring(name): value for name, value in assignments.items()}
            if any(equation.subs(substitution) != 0 for equation in equations):
                raise AssertionError("published-R17 eliminated first-jet control failed")
            record = {
                "nodes_376_377": [
                    "infinity" if node_376 is None else node_376,
                    "infinity" if node_377 is None else node_377,
                ],
                "nonzero_equation_count": len(equations),
                "status": (
                    "PASS_ELIMINATED_SYSTEM_WITNESS"
                    if len(equations) in (52, 53)
                    else "PASS_DEGENERATE_ELIMINATED_SYSTEM_WITNESS"
                ),
            }
            if degenerate_fallback is None or len(equations) > degenerate_fallback[
                "nonzero_equation_count"
            ]:
                degenerate_fallback = record
            if len(equations) in (52, 53):
                selected = record
                break
        if selected is None:
            selected = degenerate_fallback
        if selected is None:
            raise AssertionError("no published-R17 control chart found")
        chart_records.append(selected)
    return {
        "available": True,
        "model": str(model_path.relative_to(ROOT)),
        "sections": str(sections_path.relative_to(ROOT)),
        "verified_section_count": len(section_coefficients),
        "normalized_chart_controls": chart_records,
        "status": "PASS_EXACT_PROJECTIVE_FIRST_JET_SYSTEM",
    }


def main() -> None:
    arguments = parse_args()
    if shutil.which("msolve") is None:
        raise SystemExit("msolve is required")
    started = time.monotonic()
    input_raw, fibres = load_fibres(arguments.input, arguments.prime)
    field = GF(arguments.prime)
    allowed = [
        value for value in range(arguments.prime) if field(value) not in (0, 1, -1)
    ] + [None]
    pairs = [(left, right) for left in allowed for right in allowed if left != right]
    if arguments.max_pairs is not None:
        pairs = pairs[: arguments.max_pairs]
    work = arguments.work_dir / f"p{arguments.prime}"
    work.mkdir(parents=True, exist_ok=True)
    trusted_input_hashes = {}
    if arguments.reuse_unit_outputs:
        if not arguments.output.exists():
            raise SystemExit(
                "--reuse-unit-outputs requires the previous output artifact"
            )
        previous = json.loads(arguments.output.read_text())
        trusted_input_hashes = {
            tuple(row["nodes_376_377"]): row["msolve_input_sha256"]
            for row in previous.get("pair_records", [])
        }
    records = []
    solution_pairs = []
    timeout_pairs = []
    equation_count = None
    if arguments.jobs < 1 or arguments.threads < 1:
        raise ValueError("--jobs and --threads must both be positive")
    with ProcessPoolExecutor(
        max_workers=arguments.jobs,
        mp_context=get_context("fork"),
        initializer=initialize_chart_worker,
        initargs=(
            fibres,
            arguments.prime,
            work,
            arguments.threads,
            arguments.pair_timeout,
            arguments.reuse_unit_outputs,
            trusted_input_hashes,
        ),
    ) as executor:
        # executor.map preserves canonical pair order, so artifacts are
        # independent of scheduling after runtimes are removed.
        results = executor.map(build_and_solve_chart, pairs, chunksize=1)
        for pair_index, result in enumerate(results, 1):
            node_376, node_377 = result["nodes_376_377"]
            if equation_count is None:
                equation_count = result["equation_count"]
            elif equation_count != result["equation_count"]:
                raise AssertionError("the per-chart nonzero equation count changed")
            status = result["status"]
            if status == "TIMEOUT":
                timeout_pairs.append([node_376, node_377])
            elif status != "NO_GEOMETRIC_SOLUTION":
                solution_pairs.append([node_376, node_377])
            records.append(
                {
                    "nodes_376_377": [node_376, node_377],
                    "status": status,
                    "msolve_returncode": result["msolve_returncode"],
                    "msolve_input_sha256": result["msolve_input_sha256"],
                    "output_prefix": result["output_prefix"],
                    "runtime_seconds": result["runtime_seconds"],
                }
            )
            if pair_index % 25 == 0 or pair_index == len(pairs):
                print(
                    f"WGXLIFIRSTJET|prime={arguments.prime}|"
                    f"tested={pair_index}/{len(pairs)}|"
                    f"solutions={len(solution_pairs)}|"
                    f"timeouts={len(timeout_pairs)}",
                    flush=True,
                )

    complete = arguments.max_pairs is None and len(pairs) == (arguments.prime - 2) * (
        arguments.prime - 3
    )
    status = (
        "PASS_COMPLETE_DISTINCT_CHART_EMPTY"
        if complete and not solution_pairs and not timeout_pairs
        else "PASS_BOUNDED_FIRST_JET_AUDIT"
    )
    payload = {
        "schema": "icarm.wgxli-rank17-first-jet-elimination.v1",
        "status": status,
        "prime": arguments.prime,
        "input": {
            "path": str(arguments.input.resolve().relative_to(ROOT)),
            "sha256": sha256_bytes(input_raw),
        },
        "gauge": {
            "t_351": 0,
            "t_356": 1,
            "t_385": -1,
            "inverse_weierstrass_scale_351": 1,
            "enumerated_ordered_pair": ["t_376", "t_377"],
        },
        "elimination": {
            "original_continuous_unknown_count": 52,
            "free_ordinate_coefficients_eliminated": 34,
            "surface_A_B_coefficients_eliminated_before_first_jet": 12,
            "per_chart_variables": 15,
            "per_chart_variables_description": (
                "four inverse fibre scalings, five A first derivatives, five B first "
                "derivatives, and one saturation inverse"
            ),
            "per_chart_nonzero_equations": equation_count,
            "equations_description": (
                "three affine-jet compatibility equations for each of seventeen "
                "labelled sections, up to one nonzero deg(A)<=8 Hermite equation, "
                "and one nonzero-scaling saturation equation"
            ),
        },
        "chart": {
            "required_distinct_nodes": [0, 1, -1, "t_376", "t_377"],
            "parameter_space": "ordered pairs in P1(F_p) minus the three fixed nodes",
            "normalized_chart_types": [
                "both residual nodes finite",
                "t_376 at infinity",
                "t_377 at infinity",
            ],
            "ordered_pair_count": len(pairs),
            "complete": complete,
            "solution_pair_count": len(solution_pairs),
            "solution_pairs": solution_pairs,
            "timeout_pair_count": len(timeout_pairs),
            "timeout_pairs": timeout_pairs,
        },
        "positive_control": positive_control(arguments.prime),
        "solver": {
            "name": "msolve",
            "banner": msolve_version(),
            "threads": arguments.threads,
            "chart_jobs": arguments.jobs,
            "pair_timeout_seconds": arguments.pair_timeout,
            "meaning_of_empty_output": (
                "A reduced Groebner basis [1] proves no solution over the algebraic "
                "closure for that fixed parameter pair."
            ),
        },
        "pair_records": records,
        "proof_boundary": (
            "This is an exact necessary-condition computation in the distinct normalized "
            "mod-p chart for literal point labels and signs. It does not exclude a QQ "
            "model with bad/colliding reduction at this prime, a relabelled or integrally "
            "changed Mordell-Weil basis, or degree bounds outside the rootless-K3 ansatz."
        ),
        "runtime_seconds": time.monotonic() - started,
    }
    if arguments.reuse_unit_outputs:
        payload["solver"]["reuse_policy"] = (
            "reuse only an existing exact unit-ideal output [-1] or [-1]: whose "
            "rendered input hash matches the previous artifact; rerun every other chart"
        )
    if arguments.check:
        if not arguments.output.exists():
            raise SystemExit(f"missing output artifact: {arguments.output}")
        expected = json.loads(arguments.output.read_text())
        for row in expected.get("pair_records", []):
            row.pop("runtime_seconds", None)
        expected.pop("runtime_seconds", None)
        for row in payload["pair_records"]:
            row.pop("runtime_seconds", None)
        payload.pop("runtime_seconds", None)
        if expected != payload:
            raise SystemExit("stale first-jet elimination artifact")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
        arguments.output.write_text(rendered)
        print(f"WGXLIFIRSTJET|output={arguments.output}|sha256={sha256_bytes(rendered.encode())}")
    print(f"WGXLIFIRSTJET|status={status}")


if __name__ == "__main__":
    main()
