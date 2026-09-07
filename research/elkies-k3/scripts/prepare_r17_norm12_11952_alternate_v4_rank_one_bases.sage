#!/usr/bin/env sage-python
"""Prepare the seventeen exact rank-one native alternate-Q80 V4 bases.

status: ACTIVE_SEARCH
claim: primitive free generators and exact maps to the alternate parameter u
inputs: exact V4 shortlist and exact PARI rank-one screen
outputs: artifacts/generated-results/elkies-k3-r17-norm12-11952-v4-rank-one-bases-v1.json

The rank screen is performed on the Jacobian of ``w^2=q_i(u)q_j(u)``.  The
actual V4 base is a pointed quartic obtained after parametrizing the first
conic; its Jacobian is 2-isogenous to that third quotient.  This script builds
that isogeny, saturates the image on the actual paired base, and stores the
primitive rank-one generator together with the exact pointed-quartic inverse
map to ``(u,s,t)``.  Each base runs in an isolated subprocess with a hard
timeout.  Failure remains UNKNOWN and is not silently omitted.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from hashlib import sha256
import json
from pathlib import Path
import shlex
import shutil
import subprocess
import sys

from sage.all import EllipticCurve, PolynomialRing, QQ, ZZ, pari
from sage.env import SAGE_VERSION


ROOT = Path(__file__).resolve().parents[2]
SHORTLIST = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-11952-v4-pair-shortlist-64-v1.json"
RANK_SCREEN = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-11952-v4-base-rank-screen-v1.json"
OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-11952-v4-rank-one-bases-v1.json"
CHECKPOINTS = ROOT / "artifacts/local/elkies-k3/r17-norm12-11952-v4-rank-one-bases"
SCHEMA = "elkies-k3.r17-norm12-11952-v4-rank-one-bases.v1"
WORKER_PREFIX = "ALTV4BASEWORKER|"


def digest(path: Path) -> str:
    result = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            result.update(block)
    return result.hexdigest()


def display_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path.resolve())


def rational_text(value) -> str:
    value = QQ(value)
    return str(value.numerator()) if value.denominator() == 1 else f"{value.numerator()}/{value.denominator()}"


def point_text(point) -> list[str]:
    if point.is_zero():
        return ["0"]
    return [rational_text(point[0]), rational_text(point[1])]


def pointed_curve(coefficients, parameter, ordinate, variable_name):
    ring = PolynomialRing(QQ, variable_name)
    variable = ring.gen()
    quartic = ring([QQ(value) for value in coefficients])
    shifted = ring(quartic(variable + parameter))
    e, d, c, b, a = (QQ(shifted[index]) for index in range(5))
    if ordinate == 0 or ordinate**2 != e:
        raise ArithmeticError("invalid pointed-quartic base point")
    curve = EllipticCurve(
        QQ,
        [
            d / ordinate,
            c - d**2 / (4 * ordinate**2),
            2 * ordinate * b,
            -4 * ordinate**2 * a,
            a * (d**2 - 4 * ordinate**2 * c),
        ],
    )
    opposite_x = d**2 / (4 * ordinate**2) - c
    return ring, quartic, curve, opposite_x, (c, d, ordinate, parameter)


def inverse_pointed(point, constants):
    c, d, ordinate, parameter = constants
    if point.is_zero():
        return parameter, ordinate
    if point[1] == 0:
        return None
    x_value, y_value = point[:2]
    local_parameter = (4 * ordinate**2 * (x_value + c) - d**2) / (2 * ordinate * y_value)
    if local_parameter == 0:
        return parameter, -ordinate
    quartic_ordinate = (x_value * local_parameter**2 - d * local_parameter) / (2 * ordinate) - ordinate
    return parameter + local_parameter, quartic_ordinate


def worker(shortlist_path: Path, screen_path: Path, shortlist_rank: int) -> None:
    shortlist = json.loads(shortlist_path.read_text())
    screen = json.loads(screen_path.read_text())
    pair = next(row for row in shortlist["pairs"] if int(row["shortlist_rank"]) == shortlist_rank)
    rank_row = next(row for row in screen["results"] if int(row["shortlist_rank"]) == shortlist_rank)
    if (rank_row.get("rank_lower_bound"), rank_row.get("rank_upper_bound")) != (1, 1):
        raise ArithmeticError("worker target is not certified rank one")

    q_left_coefficients, q_right_coefficients = pair["branch_quadratics_coefficients_low_to_high"]
    u0 = QQ(pair["v4_base_point"]["u"])
    s0 = QQ(pair["v4_base_point"]["left_square_root"])
    t0 = QQ(pair["v4_base_point"]["right_square_root"])
    w0 = QQ(pair["v4_base_point"]["product_square_root"])

    # Produce a rational point on the third quotient, then map it to the exact
    # integral model used by the rank screen and saturate it there.
    _ring_w, product_quartic, third_pointed, opposite_x, third_constants = pointed_curve(
        pair["product_quartic_coefficients_low_to_high"], u0, w0, "z"
    )
    visible = next(point for point in third_pointed.lift_x(opposite_x, all=True) if point[1])
    third_integral = EllipticCurve(QQ, [QQ(value) for value in pair["base_jacobian_integral_a1_a2_a3_a4_a6"]])
    third_to_integral = third_pointed.isomorphism_to(third_integral)
    visible_integral = third_to_integral(visible)
    third_basis, third_index, third_regulator = third_integral.saturation([visible_integral])
    if len(third_basis) != 1 or third_basis[0].has_finite_order():
        raise ArithmeticError("third-quotient saturation did not return one nontorsion point")

    # Parametrize s^2=q_left(u) through (u0,s0).  The remaining equation is a
    # pointed quartic in the line slope r and is the actual V4 base.
    u_ring = PolynomialRing(QQ, "u")
    r_ring = PolynomialRing(QQ, "r")
    r_field = r_ring.fraction_field()
    r = r_ring.gen()
    q_left = u_ring([QQ(value) for value in q_left_coefficients])
    q_right = u_ring([QQ(value) for value in q_right_coefficients])
    leading = QQ(q_left[2])
    derivative = QQ(q_left.derivative()(u0))
    conic_denominator = 1 - leading * r**2
    u_of_r = r_field(u0 + (derivative * r**2 - 2 * s0 * r) / conic_denominator)
    q_right_at_u = sum(QQ(q_right[index]) * u_of_r**index for index in range(3))
    paired_quartic = r_ring(conic_denominator**2 * q_right_at_u)
    if paired_quartic[0] != t0**2 or not paired_quartic.is_squarefree():
        raise ArithmeticError("paired V4 quartic construction failed")
    _ring_r, rebuilt_paired_quartic, paired_pointed, _paired_opposite_x, paired_constants = pointed_curve(
        [paired_quartic[index] for index in range(5)], QQ(0), t0, "r"
    )
    if rebuilt_paired_quartic != paired_quartic:
        raise ArithmeticError("paired pointed quartic changed on rebuild")

    invariant_i = QQ(pair["binary_quartic_invariants"]["I"])
    invariant_j = QQ(pair["binary_quartic_invariants"]["J"])
    third_raw = EllipticCurve(QQ, [0, 0, 0, -27 * invariant_i, -27 * invariant_j])
    integral_to_raw = third_integral.isomorphism_to(third_raw)
    isogenies = []
    for root, _multiplicity in third_raw.division_polynomial(2).roots(QQ):
        try:
            isogenies.append(third_raw.isogeny(third_raw(root, 0), codomain=paired_pointed))
        except (ArithmeticError, ValueError):
            continue
    if len(isogenies) != 1:
        raise ArithmeticError("the expected third-quotient to V4-base 2-isogeny is not unique")
    isogeny = isogenies[0]
    image = isogeny(integral_to_raw(third_basis[0]))
    if image.is_zero() or image.has_finite_order():
        raise ArithmeticError("the 2-isogeny killed the free generator")

    paired_integral = paired_pointed.global_integral_model()
    pointed_to_integral = paired_pointed.isomorphism_to(paired_integral)
    image_integral = pointed_to_integral(image)
    # The cokernel on free Mordell--Weil groups of a degree-two isogeny is
    # killed by two.  Since the source generator is fully saturated, only
    # 2-saturation is needed on the image to obtain the primitive paired-base
    # generator; no odd prime can divide this index.
    paired_basis, paired_index, paired_regulator = paired_integral.saturation(
        [image_integral], min_prime=2, max_prime=2
    )
    if len(paired_basis) != 1 or paired_basis[0].has_finite_order():
        raise ArithmeticError("paired-base saturation did not return one nontorsion point")
    generator_integral = paired_basis[0]
    if generator_integral[1] < 0:
        generator_integral = -generator_integral
    integral_to_pointed = ~pointed_to_integral
    generator_pointed = integral_to_pointed(generator_integral)
    inverse = inverse_pointed(generator_pointed, paired_constants)
    if inverse is None:
        # At Y=0 the usual inverse has the indeterminate form 0/0.  Recover
        # the finite quartic preimage from the forward X-map instead.
        x_value = QQ(generator_pointed[0])
        inverse_polynomial = r_ring(
            ((x_value * r**2 - paired_constants[1] * r) / (2 * t0) - t0) ** 2
            - paired_quartic
        )
        candidates = []
        for root, _multiplicity in inverse_polynomial.roots(QQ):
            if root == 0:
                continue
            ordinate_value = (x_value * root**2 - paired_constants[1] * root) / (2 * t0) - t0
            if ordinate_value**2 == paired_quartic(root):
                candidates.append((root, ordinate_value))
        if len(candidates) != 1:
            raise ArithmeticError(f"exceptional pointed inverse has {len(candidates)} finite candidates")
        r_value, paired_ordinate = candidates[0]
        denominator_value = conic_denominator(r_value)
        if denominator_value == 0:
            raise ArithmeticError("exceptional generator maps to the conic infinity chart")
        u_value = u_of_r(r_value)
        s_value = s0 + (u_value - u0) / r_value
        t_value = paired_ordinate / denominator_value
    else:
        r_value, paired_ordinate = inverse
        if paired_ordinate**2 != paired_quartic(r_value):
            raise ArithmeticError("primitive generator misses the paired quartic")
        denominator_value = conic_denominator(r_value)
        if denominator_value == 0:
            raise ArithmeticError("primitive generator maps to the conic infinity chart")
        u_value = u_of_r(r_value)
        s_value = s0 if r_value == 0 else s0 + (u_value - u0) / r_value
        t_value = paired_ordinate / denominator_value
    if s_value**2 != q_left(u_value) or t_value**2 != q_right(u_value):
        raise ArithmeticError("primitive generator misses the V4 cover")

    result = {
        "shortlist_rank": shortlist_rank,
        "pair_key": pair["pair_key"],
        "labels": pair["labels"],
        "lattice_orbit_masks": pair["lattice_orbit_masks"],
        "rank_screen_interval": [1, 1],
        "third_quotient": {
            "integral_a1_a2_a3_a4_a6": [rational_text(value) for value in third_integral.a_invariants()],
            "visible_point": point_text(visible_integral),
            "saturated_generator": point_text(third_basis[0]),
            "visible_subgroup_saturation_index": int(third_index),
            "saturated_regulator_approx": str(third_regulator),
        },
        "paired_v4_base": {
            "pointed_a1_a2_a3_a4_a6": [rational_text(value) for value in paired_pointed.a_invariants()],
            "integral_a1_a2_a3_a4_a6": [rational_text(value) for value in paired_integral.a_invariants()],
            "primitive_generator_integral": point_text(generator_integral),
            "primitive_generator_pointed": point_text(generator_pointed),
            "isogeny_image_saturation_index": int(paired_index),
            "primitive_generator_height_approx": str(generator_integral.height()),
            "saturated_regulator_approx": str(paired_regulator),
            "torsion_order": int(paired_integral.torsion_order()),
        },
        "map_to_v4_cover": {
            "origin": {"u0": rational_text(u0), "s0": rational_text(s0), "t0": rational_text(t0)},
            "left_q_coefficients_low_to_high": [rational_text(value) for value in q_left],
            "right_q_coefficients_low_to_high": [rational_text(value) for value in q_right],
            "paired_quartic_coefficients_low_to_high": [rational_text(paired_quartic[index]) for index in range(5)],
            "pointed_inverse_constants": {
                "c": rational_text(paired_constants[0]),
                "d": rational_text(paired_constants[1]),
                "v0": rational_text(paired_constants[2]),
            },
            "formulas": {
                "r": "(4*v0^2*(X+c)-d^2)/(2*v0*Y)",
                "paired_ordinate": "(X*r^2-d*r)/(2*v0)-v0",
                "conic_denominator": "1-q_left[2]*r^2",
                "u": "u0+(q_left'(u0)*r^2-2*s0*r)/(1-q_left[2]*r^2)",
                "s": "s0+(u-u0)/r",
                "t": "paired_ordinate/conic_denominator",
            },
            "primitive_generator_image": {
                "r": "infinity" if r_value is None else rational_text(r_value),
                "u": rational_text(u_value),
                "s": rational_text(s_value),
                "t": rational_text(t_value),
            },
            "identities_verified": True,
        },
    }
    print(WORKER_PREFIX + json.dumps(result, sort_keys=True), flush=True)


def verify(path: Path, shortlist_path: Path, screen_path: Path) -> None:
    document = json.loads(path.read_text())
    if document.get("schema") != SCHEMA or len(document.get("bases", [])) != 17:
        raise ValueError("unexpected rank-one V4 base artifact")
    for name, expected in document["inputs"].items():
        if digest(ROOT / name) != expected:
            raise ArithmeticError(f"rank-one base input changed: {name}")
    shortlist = json.loads(shortlist_path.read_text())
    pairs = {row["pair_key"]: row for row in shortlist["pairs"]}
    for base in document["bases"]:
        pair = pairs[base["pair_key"]]
        data = base["map_to_v4_cover"]
        u_value = QQ(data["primitive_generator_image"]["u"])
        s_value = QQ(data["primitive_generator_image"]["s"])
        t_value = QQ(data["primitive_generator_image"]["t"])
        q_left = [QQ(value) for value in pair["branch_quadratics_coefficients_low_to_high"][0]]
        q_right = [QQ(value) for value in pair["branch_quadratics_coefficients_low_to_high"][1]]
        evaluate = lambda coefficients: sum(value * u_value**index for index, value in enumerate(coefficients))
        if s_value**2 != evaluate(q_left) or t_value**2 != evaluate(q_right):
            raise ArithmeticError(f"stored primitive map failed for {base['pair_key']}")
    print(f"ALTV4BASECHECK|bases=17|output={display_path(path)}", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shortlist", type=Path, default=SHORTLIST)
    parser.add_argument("--rank-screen", type=Path, default=RANK_SCREEN)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--timeout", type=float, default=180.0)
    parser.add_argument("--jobs", type=int, default=2)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--worker-shortlist-rank", type=int, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.worker_shortlist_rank is not None:
        worker(args.shortlist, args.rank_screen, args.worker_shortlist_rank)
        return
    if args.check:
        verify(args.output, args.shortlist, args.rank_screen)
        return
    if args.timeout <= 0 or args.jobs < 1:
        parser.error("timeout and jobs must be positive")

    shortlist = json.loads(args.shortlist.read_text())
    screen = json.loads(args.rank_screen.read_text())
    exact_ranks = sorted(
        int(row["shortlist_rank"])
        for row in screen["results"]
        if row.get("rank_lower_bound") == 1 and row.get("rank_upper_bound") == 1
    )
    if len(exact_ranks) != 17:
        raise ArithmeticError(f"expected seventeen exact rank-one bases, found {len(exact_ranks)}")
    pair_keys = {int(row["shortlist_rank"]): row["pair_key"] for row in shortlist["pairs"]}
    sage = shutil.which("sage")
    if sage is None:
        raise FileNotFoundError("the Sage launcher is required")
    CHECKPOINTS.mkdir(parents=True, exist_ok=True)

    def run(shortlist_rank: int) -> dict[str, object]:
        checkpoint = CHECKPOINTS / f"pair-{shortlist_rank:03d}.json"
        run_key = {
            "script_sha256": digest(Path(__file__).resolve()),
            "shortlist_sha256": digest(args.shortlist),
            "rank_screen_sha256": digest(args.rank_screen),
            "pair_key": pair_keys[shortlist_rank],
            "timeout_seconds": args.timeout,
        }
        if checkpoint.exists():
            old = json.loads(checkpoint.read_text())
            if old.get("run_key") == run_key:
                return old
        command = [
            sage,
            "-python",
            str(Path(__file__).resolve()),
            "--shortlist",
            str(args.shortlist),
            "--rank-screen",
            str(args.rank_screen),
            "--worker-shortlist-rank",
            str(shortlist_rank),
        ]
        try:
            completed = subprocess.run(command, cwd=ROOT, text=True, capture_output=True, timeout=args.timeout)
        except subprocess.TimeoutExpired:
            result = {"run_key": run_key, "status": "timeout", "pair_key": pair_keys[shortlist_rank], "shortlist_rank": shortlist_rank}
        else:
            marker = next(
                (line[len(WORKER_PREFIX):] for line in completed.stdout.splitlines() if line.startswith(WORKER_PREFIX)),
                None,
            )
            if completed.returncode or marker is None:
                result = {
                    "run_key": run_key,
                    "status": "error",
                    "pair_key": pair_keys[shortlist_rank],
                    "shortlist_rank": shortlist_rank,
                    "returncode": completed.returncode,
                    "stderr_tail": completed.stderr[-2000:],
                }
            else:
                result = {"run_key": run_key, "status": "completed", "base": json.loads(marker)}
        checkpoint.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
        return result

    results = []
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = {executor.submit(run, rank): rank for rank in exact_ranks}
        for future in as_completed(futures):
            result = future.result()
            results.append(result)
            print(
                f"ALTV4BASE|rank={result.get('shortlist_rank', result.get('base', {}).get('shortlist_rank'))}|"
                f"status={result['status']}",
                flush=True,
            )
    completed = [row["base"] for row in results if row["status"] == "completed"]
    completed.sort(key=lambda row: int(row["shortlist_rank"]))
    output = {
        "schema": SCHEMA,
        "status": "PASS_EXACT_RANK_ONE_V4_BASE_GENERATORS_AND_MAPS" if len(completed) == 17 else "INCOMPLETE_RANK_ONE_V4_BASE_PREPARATION",
        "inputs": {
            display_path(path): digest(path)
            for path in (Path(__file__).resolve(), args.shortlist, args.rank_screen)
        },
        "limits": {"worker_timeout_seconds": args.timeout, "concurrent_workers": args.jobs},
        "summary": {
            "rank_one_targets": 17,
            "completed": len(completed),
            "timeouts": sum(row["status"] == "timeout" for row in results),
            "errors": sum(row["status"] == "error" for row in results),
        },
        "bases": completed,
        "incomplete": [row for row in results if row["status"] != "completed"],
        "software_assumptions": {"sage": SAGE_VERSION, "pari": ".".join(map(str, pari.version()))},
        "reproducing_command": shlex.join(sys.argv),
        "proof_boundary": (
            "The stored base ranks are the exact [1,1] intervals from the pinned PARI screen. "
            "Each completed row constructs the actual V4 base through an exact 2-isogeny, "
            "fully saturates the third-quotient point, then uses the degree-two isogeny "
            "cokernel bound and exact 2-saturation to obtain a primitive paired-base generator, and "
            "checks its exact map to both quadratic covers. It makes no specialization-rank claim."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(f"ALTV4BASES|completed={len(completed)}/17|status={output['status']}|output={display_path(args.output)}", flush=True)


if __name__ == "__main__":
    main()
