#!/usr/bin/env sage
"""Eliminate the first jets of a literal five-fibre rootless-K3 interpolation.

The input consists of five canonical short fibres and seventeen labelled
points on each fibre.  After the weighted PGL2 and Weierstrass gauges

    t_351=0, t_356=1, t_385=-1, u_351=1,

the script exhausts the remaining ordered pair (t_376,t_377) over GF(p).
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

from sage.all import GF, PolynomialRing, prod


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


def interpolate(polynomial_ring, nodes, values):
    t = polynomial_ring.gen()
    return sum(
        polynomial_ring(
            values[index]
            / prod(nodes[index] - nodes[other] for other in range(5) if other != index)
        )
        * prod(t - nodes[other] for other in range(5) if other != index)
        for index in range(5)
    )


def build_pair_system(fibres, prime: int, node_376: int, node_377: int):
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
    polynomial_ring = PolynomialRing(ring, "t")
    t = polynomial_ring.gen()
    nodes = list(map(field, (0, 1, node_376, node_377, -1)))
    if len(set(nodes)) != 5:
        raise ValueError("the normalized nodes must be distinct")
    node_polynomial = prod(t - node for node in nodes)
    node_derivatives = [ring(node_polynomial.derivative()(node)) for node in nodes]

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
        x_interpolant = interpolate(polynomial_ring, nodes, x_values)
        y_interpolant = interpolate(polynomial_ring, nodes, y_values)
        x_derivative = x_interpolant.derivative()
        y_derivative = y_interpolant.derivative()

        ordinate_coefficients = []
        derivative_right_sides = []
        for index, node in enumerate(nodes):
            x_value = x_values[index]
            y_value = y_values[index]
            ordinate_coefficients.append(2 * y_value * node_derivatives[index])
            known_part = (
                2 * y_value * ring(y_derivative(node))
                - 3 * x_value**2 * ring(x_derivative(node))
                - family_A_values[index] * ring(x_derivative(node))
            )
            derivative_right_sides.append(
                x_value * alphas[index] + betas[index] - known_part
            )

        # The quotient
        #   (x*alpha+beta-known)/(2*y*L')
        # must be affine in the base parameter because every sextic ordinate
        # differs from its five-point quartic interpolant by L(t)*(r+s*t).
        for index in range(2, 5):
            node = ring(nodes[index])
            equations.append(
                ordinate_coefficients[0]
                * ordinate_coefficients[1]
                * derivative_right_sides[index]
                - ordinate_coefficients[index]
                * (
                    (1 - node)
                    * ordinate_coefficients[1]
                    * derivative_right_sides[0]
                    + node
                    * ordinate_coefficients[0]
                    * derivative_right_sides[1]
                )
            )

    # Five values and five first derivatives have a unique degree-at-most-nine
    # Hermite interpolant.  The K3 bound deg(A)<=8 kills its leading term.
    hermite_A = polynomial_ring.zero()
    for index, node in enumerate(nodes):
        lagrange = prod(
            (t - nodes[other]) / (node - nodes[other])
            for other in range(5)
            if other != index
        )
        lagrange_derivative = ring(lagrange.derivative()(node))
        hermite_A += (
            family_A_values[index]
            * (1 - 2 * lagrange_derivative * (t - node))
            * lagrange**2
            + alphas[index] * (t - node) * lagrange**2
        )
    equations.append(ring(hermite_A[9]))
    equations.append(inverse * prod(weights[1:]) - 1)
    equations = [equation for equation in equations if equation]
    # There are nominally 53 equations.  At p=17 the Hermite leading-term
    # equation vanishes identically because the reduced published A-values do;
    # retain the resulting 52 nonzero equations.
    if len(equations) not in (52, 53):
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


def initialize_chart_worker(fibres, prime, work, threads, pair_timeout) -> None:
    """Install read-only run data once in each chart worker."""
    global _WORKER_FIBRES, _WORKER_PRIME, _WORKER_WORK
    global _WORKER_THREADS, _WORKER_TIMEOUT
    _WORKER_FIBRES = fibres
    _WORKER_PRIME = prime
    _WORKER_WORK = work
    _WORKER_THREADS = threads
    _WORKER_TIMEOUT = pair_timeout


def build_and_solve_chart(pair) -> dict[str, object]:
    """Construct and solve one chart in an independent worker process."""
    node_376, node_377 = pair
    names, equations = build_pair_system(
        _WORKER_FIBRES, _WORKER_PRIME, node_376, node_377
    )
    stem = f"a{node_376}-b{node_377}"
    msolve_input = _WORKER_WORK / f"{stem}.ms"
    msolve_output = _WORKER_WORK / f"{stem}.solve"
    render_msolve(msolve_input, names, equations, _WORKER_PRIME)
    result = solve_chart(
        msolve_input, msolve_output, _WORKER_THREADS, _WORKER_TIMEOUT
    )
    result["nodes_376_377"] = [node_376, node_377]
    result["msolve_input_sha256"] = sha256_file(msolve_input)
    result["equation_count"] = len(equations)
    return result


def positive_control(prime: int) -> dict[str, object]:
    """Check the eliminated identities on the exact published R17 model."""
    field = GF(prime)
    model_path = ROOT / "artifacts/local/elkies-k3/q12o5867-smooth-rr-qq.json"
    sections_path = (
        ROOT / "artifacts/local/elkies-k3/q12o5867-rootless-selected-basis-qq.json"
    )
    if not model_path.exists() or not sections_path.exists():
        return {"available": False}
    model = json.loads(model_path.read_text())
    sections = json.loads(sections_path.read_text())
    try:
        A = [
            mod_fraction(field, value)
            for value in model["child"]["minimal_A_coefficients_low_to_high"]
        ]
        B = [
            mod_fraction(field, value)
            for value in model["child"]["minimal_B_coefficients_low_to_high"]
        ]
    except ZeroDivisionError:
        return {
            "available": False,
            "reason": "published-R17 control model has a coefficient denominator at this prime",
        }
    nodes = list(map(field, (0, 1, 2, 3, -1)))

    def evaluate(coefficients, value):
        return sum(coefficient * value**index for index, coefficient in enumerate(coefficients))

    def derivative_value(coefficients, value):
        return sum(
            index * coefficient * value ** (index - 1)
            for index, coefficient in enumerate(coefficients)
            if index
        )

    # Literal differentiation is enough to validate the eliminated affine-jet
    # identity, including the signs and the deg(A)<=8 Hermite constraint.
    verified = 0
    for row in sections["sections"]:
        try:
            X = [
                mod_fraction(field, value)
                for value in row["section"]["x_coefficients_low_to_high"]
            ]
            Y = [
                mod_fraction(field, value)
                for value in row["section"]["y_coefficients_low_to_high"]
            ]
        except ZeroDivisionError:
            return {
                "available": False,
                "reason": "published-R17 control section has a coefficient denominator at this prime",
            }
        for node in nodes:
            x_value = evaluate(X, node)
            y_value = evaluate(Y, node)
            if y_value**2 != x_value**3 + evaluate(A, node) * x_value + evaluate(B, node):
                raise AssertionError("published-R17 positive-control section failed")
            if (
                2 * y_value * derivative_value(Y, node)
                != 3 * x_value**2 * derivative_value(X, node)
                + derivative_value(A, node) * x_value
                + evaluate(A, node) * derivative_value(X, node)
                + derivative_value(B, node)
            ):
                raise AssertionError("published-R17 differentiated identity failed")
        verified += 1
    return {
        "available": True,
        "model": str(model_path.relative_to(ROOT)),
        "sections": str(sections_path.relative_to(ROOT)),
        "nodes": [0, 1, 2, 3, -1],
        "verified_section_count": verified,
        "status": "PASS_EXACT_DIFFERENTIATED_IDENTITIES",
    }


def main() -> None:
    arguments = parse_args()
    if shutil.which("msolve") is None:
        raise SystemExit("msolve is required")
    started = time.monotonic()
    input_raw, fibres = load_fibres(arguments.input, arguments.prime)
    field = GF(arguments.prime)
    allowed = [value for value in range(arguments.prime) if field(value) not in (0, 1, -1)]
    pairs = [(left, right) for left in allowed for right in allowed if left != right]
    if arguments.max_pairs is not None:
        pairs = pairs[: arguments.max_pairs]
    work = arguments.work_dir / f"p{arguments.prime}"
    work.mkdir(parents=True, exist_ok=True)
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

    complete = arguments.max_pairs is None and len(pairs) == (arguments.prime - 3) * (
        arguments.prime - 4
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
