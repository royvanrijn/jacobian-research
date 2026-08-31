#!/usr/bin/env python3
"""Pin the public wgxli lineages and export the five-fibre rank-17 target.

This is a source-comparative and numerical construction fingerprint.  It does
not prove that any cluster is a common family, that its generic rank is 17, or
that it is a rootless elliptic K3.  The public database snapshot and the five
target curve records are hash-pinned as retrieved on 2026-09-01.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from math import isqrt, sqrt
from pathlib import Path
import subprocess
from urllib.request import urlopen


DATABASE_SOURCE = (
    "https://elliptic-rank.icarm.cloud/database.json",
    "18699517c2969c8c3a250ae612d5caae9fb23c379fe054ba3c7fdf2ec2a83e50",
)
DATABASE_COUNT = 474
SUBMITTER = "wgxli"
SUBMITTER_CURVE_IDS = (
    351,
    356,
    363,
    364,
    376,
    377,
    378,
    385,
    389,
    390,
    391,
    393,
    395,
)
TARGET_CURVE_IDS = (351, 356, 376, 377, 385)
TARGET_SOURCES = {
    351: "02c0de1801d0c925dd6e42204f8461e99595926e95221e91da0c09466a6f67fd",
    356: "58afbc62dbb6e01b47266c90edcf0e09bb003bb6a558333422b332e42546e89e",
    376: "3ef328db3226392a86d6eb9563f90f89b1e061c1042f1568d85db744ccacd4e2",
    377: "b22e30ccadb8243a8b0c02c18761251e0e94835eee3e4e250731e94e1cd8d62a",
    385: "633538058f79df2ac75871ab8e8e892c776244590577e3fec4505ec639e58bee",
}
EXPECTED_COMPONENTS_AT_RESIDUAL_POINT_TWO = (
    (351, 356, 376, 377, 385),
    (363, 364, 378),
    (389, 390, 391),
    (393,),
    (395,),
)
EXPECTED_TARGET_FITS = {
    (351, 356): (1.4208782482875446, 0.11220111822209557, 0.9748839794656168),
    (351, 376): (0.8807498725869669, 0.1811110135333225, 0.9228352639666967),
    (351, 377): (0.8765727158608414, 0.15875125232239204, 0.9414994738265168),
    (351, 385): (1.5869071464473365, 0.14730984771630737, 0.9511000489324858),
    (356, 376): (0.6213797830800714, 0.12518078764546747, 0.9570467079699343),
    (356, 377): (0.6150947428286059, 0.13644362992241268, 0.9497468962459801),
    (356, 385): (1.1158381159097792, 0.10498027500108868, 0.9695273331633893),
    (376, 377): (0.9779111787271342, 0.1642202062367851, 0.9276831403025061),
    (376, 385): (1.7734918634396495, 0.14149560851021403, 0.9453988926225649),
    (377, 385): (1.7820511589831465, 0.16596006132173968, 0.925706239233385),
}
EXPECTED_MINIMUM_CROSS_COMPONENT_RESIDUAL = (356, 390, 0.348332561670562)


def fetch_bytes(url: str, expected_hash: str) -> bytes:
    with urlopen(url, timeout=60) as response:
        raw = response.read()
    observed_hash = hashlib.sha256(raw).hexdigest()
    if observed_hash != expected_hash:
        raise AssertionError(f"public source changed: {observed_hash} != {expected_hash}")
    return raw


def load_sources() -> tuple[dict[int, dict[str, object]], dict[int, dict[str, object]]]:
    database_raw = fetch_bytes(*DATABASE_SOURCE)
    database = json.loads(database_raw)
    if database.get("count") != DATABASE_COUNT:
        raise AssertionError("the hash-pinned database count is inconsistent")
    submitter_records = {
        int(record["id"]): record
        for record in database["curves"]
        if record.get("submitter") == SUBMITTER
    }
    if tuple(sorted(submitter_records)) != SUBMITTER_CURVE_IDS:
        raise AssertionError("the hash-pinned submitter inventory is inconsistent")

    target_records = {}
    for curve_id, expected_hash in TARGET_SOURCES.items():
        url = f"https://elliptic-rank.icarm.cloud/curve/{curve_id}.json"
        target = json.loads(fetch_bytes(url, expected_hash))
        database_record = submitter_records[curve_id]
        for key in ("id", "ainvs", "points", "rank_lower_bound", "submitter"):
            if target[key] != database_record[key]:
                raise AssertionError(f"curve {curve_id} disagrees with database at {key}")
        target_records[curve_id] = target
    return submitter_records, target_records


def rational_text(value: str | int) -> str:
    value_q = Fraction(value)
    if value_q.denominator == 1:
        return str(value_q.numerator)
    return f"({value_q.numerator}/{value_q.denominator})"


def denominator_roots(record: dict[str, object], count: int = 17) -> tuple[int, ...]:
    roots = []
    for x_value, _y_value in record["points"][:count]:
        denominator = Fraction(x_value).denominator
        root = isqrt(denominator)
        if root * root != denominator:
            raise AssertionError("an elliptic x-coordinate denominator is not a square")
        roots.append(root)
    return tuple(roots)


def height_matrices(
    records: dict[int, dict[str, object]], count: int = 17
) -> tuple[dict[int, list[float]], str]:
    program = ["default(realprecision,80);"]
    for curve_id, record in sorted(records.items()):
        coefficients = ",".join(str(value) for value in record["ainvs"])
        points = ",".join(
            f"[{rational_text(x_value)},{rational_text(y_value)}]"
            for x_value, y_value in record["points"][:count]
        )
        program.extend(
            (
                f"E=ellinit([{coefficients}]);P=[{points}];H=ellheightmatrix(E,P);",
                f'print("BEGIN|{curve_id}");',
                (
                    f'for(i=1,{count},for(j=1,{count},if(j>1,print1("|"));'
                    "print1(H[i,j]));print());"
                ),
            )
        )
    program.append('print("PARI|",version());')
    completed = subprocess.run(
        ["gp", "-q"],
        input="\n".join(program) + "\n",
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=120,
        check=True,
    )
    if completed.stderr.strip():
        raise RuntimeError(f"PARI/GP stderr: {completed.stderr.strip()}")

    matrices: dict[int, list[float]] = {}
    current: int | None = None
    pari_version = ""
    for line in completed.stdout.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith("BEGIN|"):
            current = int(line.split("|", 1)[1])
            matrices[current] = []
        elif line.startswith("PARI|"):
            pari_version = line.removeprefix("PARI|")
        else:
            if current is None:
                raise AssertionError("unexpected PARI/GP matrix output")
            matrices[current].extend(float(value) for value in line.split("|"))
    if sorted(matrices) != sorted(records) or not pari_version:
        raise AssertionError("incomplete PARI/GP output")
    if any(len(matrix) != count * count for matrix in matrices.values()):
        raise AssertionError("PARI/GP returned a matrix of the wrong size")
    return matrices, pari_version


def dot(left: list[float], right: list[float]) -> float:
    return sum(x_value * y_value for x_value, y_value in zip(left, right))


def norm(values: list[float]) -> float:
    return sqrt(dot(values, values))


def correlation(left: list[float], right: list[float]) -> float:
    left_mean = sum(left) / len(left)
    right_mean = sum(right) / len(right)
    left_centered = [value - left_mean for value in left]
    right_centered = [value - right_mean for value in right]
    return dot(left_centered, right_centered) / (
        norm(left_centered) * norm(right_centered)
    )


def fit(left: list[float], right: list[float]) -> tuple[float, float, float]:
    scale = dot(left, right) / dot(left, left)
    residual = [y_value - scale * x_value for x_value, y_value in zip(left, right)]
    return scale, norm(residual) / norm(right), correlation(left, right)


def residual_components(
    matrices: dict[int, list[float]], cutoff: float
) -> tuple[tuple[int, ...], ...]:
    parents = {curve_id: curve_id for curve_id in matrices}

    def find(curve_id: int) -> int:
        while parents[curve_id] != curve_id:
            parents[curve_id] = parents[parents[curve_id]]
            curve_id = parents[curve_id]
        return curve_id

    def union(left: int, right: int) -> None:
        left_root = find(left)
        right_root = find(right)
        if left_root != right_root:
            parents[right_root] = left_root

    curve_ids = sorted(matrices)
    for left_index, left in enumerate(curve_ids):
        for right in curve_ids[left_index + 1 :]:
            if fit(matrices[left], matrices[right])[1] <= cutoff:
                union(left, right)
    components: dict[int, list[int]] = {}
    for curve_id in curve_ids:
        components.setdefault(find(curve_id), []).append(curve_id)
    return tuple(sorted((tuple(value) for value in components.values())))


def minimum_cross_component_residual(
    matrices: dict[int, list[float]], components: tuple[tuple[int, ...], ...]
) -> tuple[int, int, float]:
    component_index = {
        curve_id: index
        for index, component in enumerate(components)
        for curve_id in component
    }
    candidates = []
    curve_ids = sorted(matrices)
    for left_index, left in enumerate(curve_ids):
        for right in curve_ids[left_index + 1 :]:
            if component_index[left] != component_index[right]:
                candidates.append((fit(matrices[left], matrices[right])[1], left, right))
    residual, left, right = min(candidates)
    return left, right, residual


def invariants(ainvs: list[str]) -> tuple[int, int, int, int]:
    a1, a2, a3, a4, a6 = (int(value) for value in ainvs)
    b2 = a1 * a1 + 4 * a2
    b4 = a1 * a3 + 2 * a4
    b6 = a3 * a3 + 4 * a6
    c4 = b2 * b2 - 24 * b4
    c6 = -b2 * b2 * b2 + 36 * b2 * b4 - 216 * b6
    return b2, c4, c6, c4 * c4 * c4 - c6 * c6


def canonical_short_record(record: dict[str, object]) -> dict[str, object]:
    a1, _a2, a3, _a4, _a6 = (int(value) for value in record["ainvs"])
    b2, c4, c6, invariant_difference = invariants(record["ainvs"])
    points = []
    short_a = -27 * c4
    short_b = -54 * c6
    for x_value, y_value in record["points"][:17]:
        x_q = Fraction(x_value)
        y_q = Fraction(y_value)
        short_x = 36 * x_q + 3 * b2
        short_y = 216 * y_q + 108 * (a1 * x_q + a3)
        if short_y * short_y != short_x * short_x * short_x + short_a * short_x + short_b:
            raise AssertionError("canonical short-model point transport failed")
        points.append((str(short_x), str(short_y)))
    return {
        "curve_id": record["id"],
        "rank_lower_bound": record["rank_lower_bound"],
        "created_at": record["created_at"],
        "source_ainvs": record["ainvs"],
        "short_model": [str(short_a), str(short_b)],
        "short_points_first_17": points,
        "transport": {
            "X": f"36*x+{3 * b2}",
            "Y": f"216*y+108*({a1}*x+{a3})",
        },
        "1728_discriminant": str(invariant_difference),
    }


def is_prime(value: int) -> bool:
    if value < 2:
        return False
    return all(value % divisor for divisor in range(2, isqrt(value) + 1))


def admissible_interpolation_primes(
    fibres: list[dict[str, object]], limit: int = 300
) -> list[int]:
    admissible = []
    for prime in range(5, limit):
        if not is_prime(prime):
            continue
        usable = True
        for fibre in fibres:
            short_a, short_b = (int(value) for value in fibre["short_model"])
            discriminant = -16 * (4 * short_a**3 + 27 * short_b**2)
            if discriminant % prime == 0:
                usable = False
                break
            if any(
                Fraction(value).denominator % prime == 0
                for point in fibre["short_points_first_17"]
                for value in point
            ):
                usable = False
                break
        if usable:
            admissible.append(prime)
    return admissible


def build_payload() -> dict[str, object]:
    submitter_records, target_records = load_sources()
    matrices, pari_version = height_matrices(submitter_records)
    components = residual_components(matrices, 0.2)
    if components != EXPECTED_COMPONENTS_AT_RESIDUAL_POINT_TWO:
        raise AssertionError(f"the bounded residual clustering changed: {components}")
    cross_gap = minimum_cross_component_residual(matrices, components)
    expected_left, expected_right, expected_residual = (
        EXPECTED_MINIMUM_CROSS_COMPONENT_RESIDUAL
    )
    if cross_gap[:2] != (expected_left, expected_right) or abs(
        cross_gap[2] - expected_residual
    ) > 1e-12:
        raise AssertionError(f"the cross-component residual gap changed: {cross_gap}")

    roots = {
        curve_id: denominator_roots(record)
        for curve_id, record in target_records.items()
    }
    target_fits = {}
    for pair, expected in EXPECTED_TARGET_FITS.items():
        observed = fit(matrices[pair[0]], matrices[pair[1]])
        for label, observed_value, expected_value in zip(
            ("scale", "residual", "correlation"), observed, expected
        ):
            if abs(observed_value - expected_value) > 1e-12:
                raise AssertionError(
                    f"target {pair} {label} changed: {observed_value} != {expected_value}"
                )
        target_fits[f"{pair[0]}/{pair[1]}"] = {
            "scale": observed[0],
            "relative_frobenius_residual": observed[1],
            "pearson_correlation_all_289_entries": observed[2],
            "ordered_denominator_root_matches": [
                [index, left]
                for index, (left, right) in enumerate(
                    zip(roots[pair[0]], roots[pair[1]]), 1
                )
                if left == right
            ],
        }

    interpolation_records = [
        canonical_short_record(target_records[curve_id])
        for curve_id in TARGET_CURVE_IDS
    ]
    return {
        "status": "bounded numerical lineage fingerprint; not family recognition",
        "snapshot": {
            "date": "2026-09-01",
            "database_url": DATABASE_SOURCE[0],
            "database_sha256": DATABASE_SOURCE[1],
            "database_curve_count": DATABASE_COUNT,
            "submitter": SUBMITTER,
            "submitter_curve_ids": list(SUBMITTER_CURVE_IDS),
            "target_curve_json_sha256": {
                str(key): value for key, value in TARGET_SOURCES.items()
            },
        },
        "residual_graph": {
            "definition": (
                "edge when the ordered first-17 canonical-height matrices have "
                "best scalar-fit relative Frobenius residual <= 0.2"
            ),
            "components": [list(component) for component in components],
            "minimum_cross_component_residual": {
                "pair": [cross_gap[0], cross_gap[1]],
                "relative_frobenius_residual": cross_gap[2],
            },
            "pari_version": pari_version,
        },
        "target_cluster": {
            "curve_ids": list(TARGET_CURVE_IDS),
            "ordered_first_17_x_denominator_roots": {
                str(curve_id): list(roots[curve_id]) for curve_id in TARGET_CURVE_IDS
            },
            "pairwise_fits": target_fits,
            "persistent_nontrivial_anchors": {"5": 71, "17": 41},
        },
        "rootless_k3_interpolation_input": {
            "hypothesis_only": True,
            "degree_bounds": {"A": 8, "B": 12, "section_x": 4, "section_y": 6},
            "short_model_convention": "Y^2=X^3-27*c4*X-54*c6",
            "unknown_fibre_scaling": (
                "For a proposed family fibre use X_k=u_k^2*x(t_k), "
                "Y_k=u_k^3*y(t_k), A_k=u_k^4*A(t_k), B_k=u_k^6*B(t_k)."
            ),
            "admissible_primes_below_300": admissible_interpolation_primes(
                interpolation_records
            ),
            "fibres": interpolation_records,
        },
        "interpretation": (
            "Curves 376, 377, and 385 are three additional public candidates for "
            "the labelled 351/356 lineage.  Exact common-family membership, base "
            "parameters, generic rank, and a rootless-K3 model remain unproved."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--write-artifact",
        type=Path,
        help="write the complete deterministic JSON payload to this path",
    )
    args = parser.parse_args()
    payload = build_payload()
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.write_artifact:
        args.write_artifact.parent.mkdir(parents=True, exist_ok=True)
        args.write_artifact.write_text(rendered)
        print(f"wrote {args.write_artifact}")
        print(f"sha256 {hashlib.sha256(rendered.encode()).hexdigest()}")
    else:
        compact = dict(payload)
        compact["rootless_k3_interpolation_input"] = {
            "hypothesis_only": True,
            "degree_bounds": payload["rootless_k3_interpolation_input"]["degree_bounds"],
            "fibre_ids": list(TARGET_CURVE_IDS),
            "full_payload": "use --write-artifact PATH",
        }
        print(json.dumps(compact, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
