#!/usr/bin/env python3
"""Build exact relative-descent inputs for four low-conductor near misses.

The output pins one exact mod-2-independent Mordell--Weil subgroup on each
curve, together with the 2-division cubic and conductor metadata.  It does not
compute a Selmer group.  For the split-infinity Mestre fibres, PARI's bounded
``ellsaturation`` is used only to discover a better exact point basis; exact
membership and a full-rank finite-quotient certificate are checked separately.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from math import isqrt
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Any, Sequence


ROOT = Path(__file__).resolve().parents[2]
CAS = Path(__file__).resolve().parent
sys.path.insert(0, str(CAS))

from mod2_reduction_independence import (  # noqa: E402
    combined_mod2_rank,
    find_mod2_reduction_certificate,
)
from search_mestre_dsquare_four import (  # noqa: E402
    FAMILIES,
    base_parameter,
    known_jacobian_points,
)
from search_mestre_root_tuple_scale import (  # noqa: E402
    point_digest,
    quartic_point_to_jacobian,
)
from structural_search import two_division_cubic  # noqa: E402


Q = Fraction
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves/"
    "conductor_first_near_miss_descent_targets_v1.json"
)
ICARM245 = (
    ROOT
    / "archive/elliptic-curves/artifacts/snapshots/pre-cleanup-2026-08-24/"
    "icarm_curve245_rank20_v1.json"
)
FERMIGIER = (
    ROOT
    / "artifacts/generated-results/elliptic-curves/"
    "elliptic_curve_candidate_fermigier_mestre_v1_u28917_20.json"
)
MESTRE = (
    ROOT
    / "archive/elliptic-curves/artifacts/snapshots/pre-cleanup-2026-08-24/"
    "elliptic_mestre_dsquare_rank19_frontiers.json"
)
EXPECTED_HASHES = {
    ICARM245: "487d6e072ed7a2508d7ab12663910b3028c8b23362039c3e8b93a278809a2cbd",
    FERMIGIER: "0a45bd473d1eba34e1b548eb71c29b0e12bfcc90d74059e5be4aefd7e236d149",
    MESTRE: "e78613cc35ad523242a6d3af529a4b59bece136a8f2b7880bead3ef7094144be",
}
SATURATION_PRIME_BOUND = 3
GP_STACK_BYTES = 2_000_000_000
GP_TIMEOUT_SECONDS = 60.0
CERTIFICATE_PRIME_BOUND = 2_000


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def checked_json(path: Path) -> dict[str, Any]:
    actual = sha256_file(path)
    expected = EXPECTED_HASHES[path]
    if actual != expected:
        raise AssertionError(f"changed input {path}: {actual} != {expected}")
    return json.loads(path.read_text())


def rational_text(value: Fraction | int | str) -> str:
    value = Q(value)
    return str(value.numerator) if value.denominator == 1 else str(value)


def point_record(point: tuple[Fraction, Fraction]) -> list[str]:
    return [rational_text(point[0]), rational_text(point[1])]


def on_curve(
    model: Sequence[Fraction], point: tuple[Fraction, Fraction]
) -> bool:
    a1, a2, a3, a4, a6 = map(Q, model)
    x, y = map(Q, point)
    return y * y + a1 * x * y + a3 * y == x**3 + a2 * x * x + a4 * x + a6


def discriminant(model: Sequence[Fraction]) -> Fraction:
    a1, a2, a3, a4, a6 = map(Q, model)
    b2 = a1 * a1 + 4 * a2
    b4 = a1 * a3 + 2 * a4
    b6 = a3 * a3 + 4 * a6
    b8 = a1 * a1 * a6 + 4 * a2 * a6 - a1 * a3 * a4 + a2 * a3 * a3 - a4 * a4
    return -b2 * b2 * b8 - 8 * b4**3 - 27 * b6**2 + 9 * b2 * b4 * b6


def short_model_and_points(
    model: Sequence[Fraction], points: Sequence[tuple[Fraction, Fraction]]
) -> tuple[
    tuple[Fraction, ...],
    tuple[tuple[Fraction, Fraction], ...],
    dict[str, str],
]:
    model = tuple(map(Q, model))
    if not any(model[:3]):
        return model, tuple(points), {"X": "x", "Y": "y"}
    a1, a2, a3, a4, a6 = model
    b2 = a1 * a1 + 4 * a2
    b4 = a1 * a3 + 2 * a4
    b6 = a3 * a3 + 4 * a6
    c4 = b2 * b2 - 24 * b4
    c6 = -b2**3 + 36 * b2 * b4 - 216 * b6
    short = (Q(0), Q(0), Q(0), -27 * c4, -54 * c6)
    transformed = tuple(
        (36 * x + 3 * b2, 108 * (2 * y + a1 * x + a3))
        for x, y in points
    )
    if any(not on_curve(short, point) for point in transformed):
        raise ArithmeticError("generalized-to-short point transformation failed")
    return short, transformed, {
        "X": f"36*x+({rational_text(3 * b2)})",
        "Y": f"108*(2*y+({rational_text(a1)})*x+({rational_text(a3)}))",
    }


def coefficient_mod(value: Fraction, prime: int) -> int:
    value = Q(value)
    if value.denominator % prime == 0:
        raise ZeroDivisionError
    return value.numerator * pow(value.denominator, -1, prime) % prime


def point_count_mod_prime(model: Sequence[Fraction], prime: int) -> int:
    coefficients = [coefficient_mod(Q(value), prime) for value in model]
    a1, a2, a3, a4, a6 = coefficients
    total = prime + 1
    character_sum = 0
    for x in range(prime):
        right = (x**3 + a2 * x * x + a4 * x + a6) % prime
        dy = ((a1 * x + a3) ** 2 + 4 * right) % prime
        if dy:
            character_sum += 1 if pow(dy, (prime - 1) // 2, prime) == 1 else -1
    return total + character_sum


def small_primes(limit: int) -> tuple[int, ...]:
    result = []
    for value in range(3, limit + 1, 2):
        if all(value % divisor for divisor in range(2, isqrt(value) + 1)):
            result.append(value)
    return tuple(result)


def no_rational_two_torsion_witness(model: Sequence[Fraction]) -> dict[str, int]:
    delta = discriminant(model)
    for prime in small_primes(1_000):
        try:
            delta_mod = coefficient_mod(delta, prime)
            order = point_count_mod_prime(model, prime)
        except ZeroDivisionError:
            continue
        if delta_mod and order % 2:
            return {"prime": prime, "group_order": order}
    raise ArithmeticError("no odd good-reduction order found through 1000")


def gp_rational(value: Fraction) -> str:
    value = Q(value)
    return f"({value.numerator}/{value.denominator})"


def pari_small_prime_saturation(
    model: Sequence[Fraction], points: Sequence[tuple[Fraction, Fraction]]
) -> tuple[tuple[Fraction, Fraction], ...]:
    gp = shutil.which("gp")
    if gp is None:
        raise FileNotFoundError("PARI/GP executable 'gp' is required")
    model_text = ",".join(gp_rational(Q(value)) for value in model)
    lines = [
        f"default(parisizemax,{GP_STACK_BYTES});",
        f"E=ellinit([{model_text}]);",
        "P=[];",
    ]
    for x, y in points:
        lines.append(
            f"P=concat(P,[[{gp_rational(x)},{gp_rational(y)}]]);"
        )
    lines.extend(
        (
            f"S=ellsaturation(E,P,{SATURATION_PRIME_BOUND});",
            'print("COUNT\\t",#S);',
            'for(i=1,#S,print("POINT\\t",S[i][1],"\\t",S[i][2]));',
            "quit",
        )
    )
    completed = subprocess.run(
        [gp, "-q", "-s", str(GP_STACK_BYTES), "-f"],
        input="\n".join(lines) + "\n",
        text=True,
        capture_output=True,
        check=True,
        timeout=GP_TIMEOUT_SECONDS,
    )
    if any(
        "***" in line and "Warning: new maximum stack size" not in line
        for line in (completed.stdout + completed.stderr).splitlines()
    ):
        raise RuntimeError(completed.stdout + completed.stderr)
    count = None
    saturated = []
    for line in completed.stdout.splitlines():
        fields = line.split("\t")
        if fields[0] == "COUNT":
            count = int(fields[1])
        elif fields[0] == "POINT":
            saturated.append((Q(fields[1]), Q(fields[2])))
    if count != len(points) or len(saturated) != count:
        raise ArithmeticError("PARI saturation output changed rank or shape")
    return tuple(saturated)


def mod2_certificate(
    model: Sequence[Fraction], points: Sequence[tuple[Fraction, Fraction]]
) -> dict[str, Any]:
    short_model, short_points, transformation = short_model_and_points(model, points)
    signatures = find_mod2_reduction_certificate(
        short_model, short_points, prime_bound=CERTIFICATE_PRIME_BOUND
    )
    rank = combined_mod2_rank(signatures, len(points))
    if rank != len(points):
        raise ArithmeticError(f"mod-2 Kummer image rank {rank} != {len(points)}")
    return {
        "method": "exact images in products of E(F_p)/2E(F_p)",
        "short_model": [rational_text(value) for value in short_model],
        "point_transformation": transformation,
        "combined_binary_rank": rank,
        "prime_bound": CERTIFICATE_PRIME_BOUND,
        "primes": [row.prime for row in signatures],
        "rows": [
            {
                "prime": row.prime,
                "group_order": row.group_order,
                "quotient_dimension": row.quotient_dimension,
                "matrix_rows": [list(bits) for bits in row.rows],
            }
            for row in signatures
        ],
    }


def target_record(
    *,
    target_id: str,
    label: str,
    model: Sequence[Fraction],
    global_minimal_model: Sequence[Fraction],
    points: Sequence[tuple[Fraction, Fraction]],
    known_rank: int,
    conductor: str,
    log_conductor: str,
    root_number: int,
    source_artifact: Path,
    basis_provenance: dict[str, Any],
) -> dict[str, Any]:
    model = tuple(map(Q, model))
    points = tuple((Q(x), Q(y)) for x, y in points)
    if len(points) != known_rank or len(set(points)) != known_rank:
        raise AssertionError(f"{target_id} basis count changed")
    if any(not on_curve(model, point) for point in points):
        raise ArithmeticError(f"{target_id} basis contains an off-curve point")
    certificate = mod2_certificate(model, points)
    torsion_witness = no_rational_two_torsion_witness(model)
    return {
        "id": target_id,
        "label": label,
        "descent_model": [rational_text(value) for value in model],
        "global_minimal_model": [rational_text(value) for value in global_minimal_model],
        "conductor": str(conductor),
        "log_conductor": str(log_conductor),
        "root_number": int(root_number),
        "certified_known_rank": known_rank,
        "rank21_directions_needed": 21 - known_rank,
        "rational_two_torsion_dimension": 0,
        "no_rational_two_torsion_witness": torsion_witness,
        "two_division_cubic_coefficients_ascending": list(two_division_cubic(model)),
        "known_basis": [point_record(point) for point in points],
        "known_basis_sha256": point_digest(points),
        "known_basis_mod2_certificate": certificate,
        "known_kummer_image_dimension": certificate["combined_binary_rank"],
        "basis_provenance": basis_provenance,
        "source_artifact": {
            "path": str(source_artifact.relative_to(ROOT)),
            "sha256": EXPECTED_HASHES[source_artifact],
        },
        "next_exact_gate": (
            "compute the complete 2-Selmer dimension; subtract the certified "
            "known Kummer image dimension, then materialize only a basis of "
            "the residual quotient as minimized locally tested covers"
        ),
    }


def build_icarm245(data: dict[str, Any]) -> dict[str, Any]:
    model = tuple(Q(value) for value in data["curve"]["global_minimal_model"])
    points = tuple((Q(row[0]), Q(row[1])) for row in data["points"])
    return target_record(
        target_id="icarm-245",
        label="ICARM 245",
        model=model,
        global_minimal_model=model,
        points=points,
        known_rank=20,
        conductor=data["curve"]["conductor"],
        log_conductor=data["curve"]["log_conductor_numeric_80_digits"],
        root_number=data["curve"]["root_number"],
        source_artifact=ICARM245,
        basis_provenance={
            "method": "public basis with exact repository finite-quotient replay",
            "claim": "twenty exact independent points with full mod-2 Kummer image",
        },
    )


def build_fermigier(data: dict[str, Any]) -> dict[str, Any]:
    model = tuple(
        Q(value)
        for value in data["models"]["legacy_normalized_short_jacobian"]["coefficients"]
    )
    minimal = tuple(
        Q(value) for value in data["models"]["global_minimal"]["coefficients"]
    )
    saturation = data["bounded_saturation_status"]
    points = tuple(
        (Q(record["x"]), Q(record["y"]))
        for record in saturation["returned_legacy_basis"]
    )
    arithmetic = data["global_arithmetic"]
    return target_record(
        target_id="fermigier-u28917-20",
        label="Fermigier--Mestre u=28917/20",
        model=model,
        global_minimal_model=minimal,
        points=points,
        known_rank=20,
        conductor=arithmetic["conductor"],
        log_conductor=arithmetic["log_conductor"],
        root_number=arithmetic["root_number"],
        source_artifact=FERMIGIER,
        basis_provenance={
            "method": "PARI bounded small-prime saturation candidate, followed by exact independent certification",
            "prime_bound_strict_upper_limit": saturation["prime_bound_strict_upper_limit"],
            "scope_warning": saturation["scope_warning"],
            "claim": "twenty exact independent points with full mod-2 Kummer image; no global saturation claim",
        },
    )


def mestre_seed_basis(record: dict[str, Any]) -> tuple[tuple[Fraction, Fraction], ...]:
    family = FAMILIES[int(record["family_index"])]
    parameter_u = Q(record["u"])
    parameter_t = base_parameter(family, parameter_u)
    by_x = {point[0]: point for point in known_jacobian_points(family, parameter_u)}
    for row in record["point_search"]["searched_points"]:
        point = quartic_point_to_jacobian(
            family.construction,
            parameter_t,
            (Q(row["quartic_x"]), Q(row["quartic_y_positive"])),
        )
        by_x.setdefault(point[0], point)
    pool = tuple(by_x.values())
    indices = tuple(
        int(value) - 1
        for value in record["exact_rank_certificate"][
            "independent_subset_indices_one_based"
        ]
    )
    seed = tuple(pool[index] for index in indices)
    expected = record["exact_rank_certificate"]["independent_subset_sha256"]
    if point_digest(seed) != expected:
        raise AssertionError(f"{record['label']} seed basis hash changed")
    return seed


def build_mestre(record: dict[str, Any]) -> dict[str, Any]:
    model = tuple(Q(value) for value in record["short_weierstrass_coefficients"])
    seed = mestre_seed_basis(record)
    points = pari_small_prime_saturation(model, seed)
    return target_record(
        target_id=record["label"].replace("_", "-"),
        label=record["label"],
        model=model,
        global_minimal_model=tuple(Q(value) for value in record["minimal_model"]),
        points=points,
        known_rank=19,
        conductor=record["conductor"],
        log_conductor=record["log_conductor"],
        root_number=record["root_number"],
        source_artifact=MESTRE,
        basis_provenance={
            "method": "PARI ellsaturation discovery from the certified mod-3 basis, followed by exact mod-2 finite-quotient certification",
            "saturation_prime_bound": SATURATION_PRIME_BOUND,
            "seed_basis_sha256": record["exact_rank_certificate"][
                "independent_subset_sha256"
            ],
            "claim": "nineteen exact independent points with full mod-2 Kummer image; no global saturation claim",
        },
    )


def build_manifest() -> dict[str, Any]:
    icarm = checked_json(ICARM245)
    fermigier = checked_json(FERMIGIER)
    mestre = checked_json(MESTRE)
    targets = [
        build_icarm245(icarm),
        build_fermigier(fermigier),
        *(build_mestre(record) for record in mestre["frontiers"]),
    ]
    return {
        "schema": "elliptic-curves.conductor-first-near-miss-descent-targets.v1",
        "claim_level": "exact descent inputs and known-subgroup Kummer certificates",
        "programme_order": [
            "complete 2-Selmer",
            "quotient by the pinned known Kummer image",
            "close the fibre if the residual quotient is zero",
            "otherwise construct, minimize, and locally test residual covers",
            "use 4-descent or independent 3-descent only on surviving classes",
        ],
        "targets": targets,
        "limitations": [
            "this artifact does not compute a complete Selmer group or rank upper bound",
            "bounded PARI saturation is used only to discover exact points and is not promoted to a global saturation theorem",
            "a positive residual Selmer quotient may consist partly or entirely of Tate--Shafarevich classes",
        ],
        "reproducing_command": (
            ".venv/bin/python elliptic-curves/cas/"
            "build_conductor_first_near_miss_targets.py --check"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    payload = build_manifest()
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.exists() or args.output.read_text() != rendered:
            raise SystemExit(f"stale or missing artifact: {args.output}")
        print(f"PASS {args.output}")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
        print(f"WROTE {args.output}")
    for target in payload["targets"]:
        print(
            f"{target['id']}: rank>={target['certified_known_rank']} "
            f"known_mod2={target['known_kummer_image_dimension']} "
            f"logN={target['log_conductor']}"
        )


if __name__ == "__main__":
    main()
