#!/usr/bin/env sage -python
"""Bound the signed-permutation rebasing of the five wgxli rank-17 records.

This script is a numerical proposal and exact-fingerprint stage, not a family
recognition theorem.  It exhausts every relative diagonal sign for each pair
of canonical-height Grams, constructs a finite permutation graph from stable
denominator anchors, integrality, exact identity-component data, and
sign-invariant height rows, and enumerates every perfect matching in that
declared graph.  Only retained signed permutations are inputs to the exact
point/group-law and modular first-jet stages.  The default source is the later
committed sufficient projection of all thirteen public records;
``--live-pinned-source`` retains the original raw database/hash audit.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from math import isqrt, sqrt
from pathlib import Path
import subprocess
import time
from urllib.request import urlopen

import numpy
from sage.all import EllipticCurve, ZZ


ROOT = Path(__file__).resolve().parents[2]
LINEAGE = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "icarm_wgxli_rank17_lineage_v1.json"
)
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "icarm_wgxli_rank17_signed_permutation_rebasing_v1.json"
)
DATABASE_URL = "https://elliptic-rank.icarm.cloud/database.json"
DATABASE_SHA256 = "18699517c2969c8c3a250ae612d5caae9fb23c379fe054ba3c7fdf2ec2a83e50"
PUBLIC_FIBRE_PROJECTION = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-icarm-public-fibres-v1.json"
)
PUBLIC_FIBRE_PROJECTION_SHA256 = (
    "9a2675ab48cc37111d1f4050bd1797fc84c98b7839668d292d11406efe7a9eaa"
)
TARGETS = (351, 356, 376, 377, 385)
REFERENCE = 356
EXPECTED_COMPONENTS = (
    (351, 356, 376, 377, 385),
    (363, 364, 378),
    (389, 390, 391),
    (393,),
    (395,),
)
SIGN_NEAR_RELATIVE_BAND = 0.10
PERMUTATION_NEAR_RELATIVE_BAND = 0.10
HEIGHT_ROW_ABSOLUTE_SLACK = 0.025
HEIGHT_ROW_RELATIVE_SLACK = 0.20
LOCAL_COMPONENT_SCORE_SLACK = 0.20


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=LINEAGE)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    parser.add_argument(
        "--live-pinned-source",
        action="store_true",
        help="require the original 474-row public database byte hash",
    )
    return parser.parse_args()


def sha256_bytes(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def load_database(*, live_pinned_source):
    if live_pinned_source:
        with urlopen(DATABASE_URL, timeout=60) as response:
            raw = response.read()
        source_hash = sha256_bytes(raw)
        if source_hash != DATABASE_SHA256:
            raise AssertionError("the pinned ICARM database snapshot changed")
        payload = json.loads(raw)
        records = {
            int(record["id"]): record
            for record in payload["curves"]
            if record.get("submitter") == "wgxli"
        }
        return source_hash, records

    raw = PUBLIC_FIBRE_PROJECTION.read_bytes()
    observed = sha256_bytes(raw)
    if observed != PUBLIC_FIBRE_PROJECTION_SHA256:
        raise AssertionError(
            f"committed public-fibre projection changed: {observed} != "
            f"{PUBLIC_FIBRE_PROJECTION_SHA256}"
        )
    payload = json.loads(raw)
    by_id = {int(record["id"]): record for record in payload.get("records", [])}
    records = {}
    for curve_id in sorted({value for component in EXPECTED_COMPONENTS for value in component}):
        record = by_id.get(curve_id)
        if record is None:
            raise AssertionError(f"committed public-fibre projection omitted {curve_id}")
        if len(record.get("points", [])) < 17:
            raise AssertionError(f"committed public-fibre projection lacks 17 points for {curve_id}")
        records[curve_id] = {
            **record,
            "rank_lower_bound": record["snapshot_rank_lower_bound"],
            "submitter": "wgxli",
        }
    return DATABASE_SHA256, records


def rational_text(value):
    value = Fraction(value)
    return (
        str(value.numerator)
        if value.denominator == 1
        else f"({value.numerator}/{value.denominator})"
    )


def height_matrices(records):
    program = ["default(realprecision,80);"]
    for curve_id, record in sorted(records.items()):
        coefficients = ",".join(str(value) for value in record["ainvs"])
        points = ",".join(
            f"[{rational_text(x_value)},{rational_text(y_value)}]"
            for x_value, y_value in record["points"][:17]
        )
        program.extend((
            f"E=ellinit([{coefficients}]);P=[{points}];H=ellheightmatrix(E,P);",
            f'print("BEGIN|{curve_id}");',
            "for(i=1,17,for(j=1,17,if(j>1,print1(\"|\"));print1(H[i,j]));print());",
        ))
    program.append('print("PARI|",version());')
    completed = subprocess.run(
        ["gp", "-q"],
        input="\n".join(program) + "\n",
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=180,
        check=True,
    )
    if completed.stderr.strip():
        raise RuntimeError(completed.stderr.strip())
    matrices = {}
    current = None
    pari_version = None
    for line in completed.stdout.splitlines():
        line = line.strip()
        if line.startswith("BEGIN|"):
            current = int(line.split("|", 1)[1])
            matrices[current] = []
        elif line.startswith("PARI|"):
            pari_version = line.split("|", 1)[1]
        elif line:
            matrices[current].extend(map(float, line.split("|")))
    matrices = {
        key: numpy.array(values, dtype=float).reshape(17, 17)
        for key, values in matrices.items()
    }
    if pari_version is None or any(matrix.shape != (17, 17) for matrix in matrices.values()):
        raise AssertionError("incomplete PARI height output")
    return matrices, pari_version


def all_signs():
    signs = numpy.ones((1 << 16, 17), dtype=numpy.int8)
    bits = numpy.arange(1 << 16, dtype=numpy.uint32)[:, None]
    signs[:, 1:] = 1 - 2 * ((bits >> numpy.arange(16, dtype=numpy.uint32)) & 1)
    return signs


SIGNS = all_signs()


def sign_string(signs) -> str:
    return "".join("+" if value == 1 else "-" for value in signs)


def exhaustive_sign_fit(left, right, retain_count=32):
    coefficient = left * right
    dots = numpy.einsum("ni,ij,nj->n", SIGNS, coefficient, SIGNS, optimize=True)
    left_norm_squared = float(numpy.sum(left * left))
    right_norm_squared = float(numpy.sum(right * right))
    scales = dots / left_norm_squared
    residual_squared = 1.0 - dots * dots / (left_norm_squared * right_norm_squared)
    residual_squared = numpy.maximum(residual_squared, 0.0)
    residual_squared = numpy.where(scales > 0, residual_squared, 1.0)
    order = numpy.argsort(residual_squared)
    retained = []
    for index in order[:retain_count]:
        retained.append({
            "relative_signs": sign_string(SIGNS[index]),
            "scale_left_to_right": float(scales[index]),
            "relative_frobenius_residual": sqrt(float(residual_squared[index])),
        })
    return retained


def fit_residual(left, right):
    scale = float(numpy.sum(left * right) / numpy.sum(left * left))
    residual = float(numpy.linalg.norm(right - scale * left) / numpy.linalg.norm(right))
    return scale, residual


def residual_components(matrices, cutoff=0.2):
    parents = {curve_id: curve_id for curve_id in matrices}

    def find(curve_id):
        while parents[curve_id] != curve_id:
            parents[curve_id] = parents[parents[curve_id]]
            curve_id = parents[curve_id]
        return curve_id

    def union(left, right):
        left = find(left)
        right = find(right)
        if left != right:
            parents[right] = left

    pair_records = {}
    curve_ids = sorted(matrices)
    for left_index, left_id in enumerate(curve_ids):
        for right_id in curve_ids[left_index + 1 :]:
            best = exhaustive_sign_fit(matrices[left_id], matrices[right_id], 2)
            pair_records[f"{left_id}/{right_id}"] = best[0]
            if best[0]["relative_frobenius_residual"] <= cutoff:
                union(left_id, right_id)
    grouped = {}
    for curve_id in curve_ids:
        grouped.setdefault(find(curve_id), []).append(curve_id)
    return tuple(sorted(tuple(value) for value in grouped.values())), pair_records


def valuation(value, prime):
    value = Fraction(value)
    if not value:
        return 10**9
    numerator = abs(value.numerator)
    denominator = value.denominator
    answer = 0
    while numerator % prime == 0:
        numerator //= prime
        answer += 1
    while denominator % prime == 0:
        denominator //= prime
        answer -= 1
    return answer


def singular_reduction(point, coefficients, prime):
    x_value, y_value = map(Fraction, point)
    if valuation(x_value, prime) < 0 or valuation(y_value, prime) < 0:
        return False
    a1, a2, a3, a4, _a6 = map(Fraction, coefficients)
    derivative_x = a1 * y_value - 3 * x_value**2 - 2 * a2 * x_value - a4
    derivative_y = 2 * y_value + a1 * x_value + a3
    return valuation(derivative_x, prime) > 0 and valuation(derivative_y, prime) > 0


def local_component_fingerprints(record):
    curve = EllipticCurve(list(map(ZZ, record["ainvs"])))
    places = []
    for prime_text in record["bad_primes"]:
        prime = ZZ(prime_text)
        data = curve.local_data(prime)
        if data.minimal_model().a_invariants() != curve.a_invariants():
            raise AssertionError(
                f"curve {record['id']} is not minimal at the fingerprint prime {prime}"
            )
        tamagawa = int(data.tamagawa_number())
        if tamagawa > 1:
            places.append((int(prime), str(data.kodaira_symbol()), tamagawa))
    point_records = []
    for index, point in enumerate(record["points"][:17], 1):
        statuses = [
            int(singular_reduction(point, record["ainvs"], prime))
            for prime, _symbol, _tamagawa in places
        ]
        weighted = sum(
            status * (tamagawa - 1) / tamagawa
            for status, (_prime, _symbol, tamagawa) in zip(statuses, places)
        )
        point_records.append({
            "point_index": index,
            "identity_component_bits": statuses,
            "nonidentity_place_count": sum(statuses),
            "normalized_nonidentity_fraction": (
                sum(statuses) / len(places) if places else 0.0
            ),
            "normalized_tamagawa_weight": weighted / len(places) if places else 0.0,
        })
    return {
        "places": [
            {"prime": prime, "kodaira_symbol": symbol, "tamagawa_number": tamagawa}
            for prime, symbol, tamagawa in places
        ],
        "points": point_records,
        "meaning": (
            "Exact identity/nonidentity component status on the supplied global minimal "
            "model; oriented component labels are deliberately not inferred."
        ),
    }


def height_row_features(gram, local_records):
    diagonal = numpy.diag(gram)
    correlations = numpy.abs(gram) / numpy.sqrt(diagonal[:, None] * diagonal[None, :])
    answer = []
    for index in range(17):
        height_part = numpy.concatenate((
            [2.0 * diagonal[index] / numpy.mean(diagonal)],
            numpy.sort(numpy.delete(correlations[index], index)),
        ))
        local = local_records[index]
        answer.append(numpy.concatenate((
            height_part,
            [
                0.1 * local["normalized_nonidentity_fraction"],
                0.1 * local["normalized_tamagawa_weight"],
            ],
        )))
    return numpy.array(answer)


def enumerate_matchings(candidate_sets, cap=100000):
    order = sorted(range(17), key=lambda index: (len(candidate_sets[index]), index))
    answers = []
    permutation = [-1] * 17

    def recurse(depth, used):
        if len(answers) >= cap:
            raise ArithmeticError("declared permutation graph exceeded its safety cap")
        if depth == 17:
            answers.append(tuple(permutation))
            return
        source = order[depth]
        for target in candidate_sets[source]:
            if target in used:
                continue
            permutation[source] = target
            used.add(target)
            recurse(depth + 1, used)
            used.remove(target)

    recurse(0, set())
    return answers


def permutation_search(lineage, matrices, local_data):
    roots = {
        int(curve_id): values
        for curve_id, values in lineage["target_cluster"][
            "ordered_first_17_x_denominator_roots"
        ].items()
    }
    stable_positions = []
    stable_records = []
    for index in range(17):
        values = [roots[curve_id][index] for curve_id in TARGETS]
        counts = {value: values.count(value) for value in set(values)}
        mode = min(counts, key=lambda value: (-counts[value], value))
        if counts[mode] >= 4:
            stable_positions.append(index)
            stable_records.append({
                "label": index + 1,
                "modal_denominator_root": mode,
                "support": counts[mode],
            })

    reference_features = height_row_features(
        matrices[REFERENCE], local_data[REFERENCE]["points"]
    )
    records = {}
    retained_per_fibre = {}
    for curve_id in TARGETS:
        if curve_id == REFERENCE:
            identity = tuple(range(17))
            records[str(curve_id)] = {
                "candidate_sets": [[index + 1] for index in identity],
                "perfect_matching_count": 1,
                "matching_records": [{
                    "permutation_new_label_to_old_point": [index + 1 for index in identity],
                    "signless_residual_against_356": 0.0,
                }],
            }
            retained_per_fibre[curve_id] = [(identity, "+" * 17, 0.0)]
            continue

        target_features = height_row_features(
            matrices[curve_id], local_data[curve_id]["points"]
        )
        costs = numpy.linalg.norm(
            reference_features[:, None, :] - target_features[None, :, :], axis=2
        )
        candidate_sets = []
        for index in range(17):
            if index in stable_positions:
                candidates = [index]
            else:
                baseline = costs[index, index]
                candidates = []
                for target in range(17):
                    if target in stable_positions:
                        continue
                    local_left = local_data[REFERENCE]["points"][index]
                    local_right = local_data[curve_id]["points"][target]
                    local_gap = abs(
                        local_left["normalized_nonidentity_fraction"]
                        - local_right["normalized_nonidentity_fraction"]
                    )
                    if (
                        costs[index, target] <= baseline + HEIGHT_ROW_ABSOLUTE_SLACK
                        and costs[index, target]
                        <= (1 + HEIGHT_ROW_RELATIVE_SLACK) * baseline
                        and local_gap <= LOCAL_COMPONENT_SCORE_SLACK
                    ):
                        candidates.append(target)
                # The displayed label is the control member of every ambiguity
                # class; it is never removed by a noisy numerical fingerprint.
                if index not in candidates:
                    candidates.append(index)
            candidate_sets.append(sorted(set(candidates)))

        matchings = enumerate_matchings(candidate_sets)
        identity = tuple(range(17))
        if identity not in matchings:
            raise AssertionError("the displayed-order control left the permutation graph")
        matching_records = []
        for permutation in matchings:
            permuted = matrices[curve_id][list(permutation)][:, list(permutation)]
            _scale, residual = fit_residual(matrices[REFERENCE], permuted)
            matching_records.append({
                "permutation_new_label_to_old_point": [value + 1 for value in permutation],
                "signless_residual_against_356": residual,
            })
        matching_records.sort(
            key=lambda row: (
                row["signless_residual_against_356"],
                row["permutation_new_label_to_old_point"],
            )
        )
        identity_residual = next(
            row["signless_residual_against_356"]
            for row in matching_records
            if row["permutation_new_label_to_old_point"] == list(range(1, 18))
        )
        retained = []
        for row in matching_records:
            if row["signless_residual_against_356"] > (
                (1 + PERMUTATION_NEAR_RELATIVE_BAND) * identity_residual
            ):
                continue
            permutation = tuple(
                value - 1 for value in row["permutation_new_label_to_old_point"]
            )
            permuted = matrices[curve_id][list(permutation)][:, list(permutation)]
            best_sign = exhaustive_sign_fit(matrices[REFERENCE], permuted, 2)[0]
            retained.append((
                permutation,
                best_sign["relative_signs"],
                best_sign["relative_frobenius_residual"],
            ))
        retained_per_fibre[curve_id] = retained
        records[str(curve_id)] = {
            "candidate_sets_new_label_to_allowed_old_points": [
                [value + 1 for value in candidates] for candidates in candidate_sets
            ],
            "perfect_matching_count": len(matchings),
            "matching_records": matching_records,
            "retained_within_10_percent_of_displayed_residual": [
                {
                    "permutation_new_label_to_old_point": [value + 1 for value in permutation],
                    "relative_signs": signs,
                    "residual_against_356": residual,
                }
                for permutation, signs, residual in retained
            ],
        }

    return {
        "stable_denominator_anchors": stable_records,
        "reference_curve": REFERENCE,
        "height_row_absolute_slack": HEIGHT_ROW_ABSOLUTE_SLACK,
        "height_row_relative_slack": HEIGHT_ROW_RELATIVE_SLACK,
        "local_component_fraction_slack": LOCAL_COMPONENT_SCORE_SLACK,
        "per_fibre": records,
        "retained_per_fibre": retained_per_fibre,
    }


def main():
    arguments = parse_args()
    started = time.monotonic()
    lineage_raw = arguments.input.read_bytes()
    lineage = json.loads(lineage_raw)
    database_source_hash, submitter_records = load_database(
        live_pinned_source=arguments.live_pinned_source
    )
    matrices, pari_version = height_matrices(submitter_records)

    components, all_signed_pair_fits = residual_components(matrices)
    if components != EXPECTED_COMPONENTS:
        raise AssertionError(f"signed negative-control components changed: {components}")

    local_data = {
        curve_id: local_component_fingerprints(submitter_records[curve_id])
        for curve_id in TARGETS
    }

    sign_pairs = {}
    pair_minimum_sum = 0.0
    minimum_nontrivial_increment = float("inf")
    for left_index, left_id in enumerate(TARGETS):
        for right_id in TARGETS[left_index + 1 :]:
            retained = exhaustive_sign_fit(matrices[left_id], matrices[right_id], 32)
            best = retained[0]
            second = retained[1]
            if best["relative_signs"] != "+" * 17:
                raise AssertionError("a target pair changed its optimum relative sign")
            best_squared = best["relative_frobenius_residual"] ** 2
            second_squared = second["relative_frobenius_residual"] ** 2
            pair_minimum_sum += best_squared
            minimum_nontrivial_increment = min(
                minimum_nontrivial_increment, second_squared - best_squared
            )
            cutoff = (1 + SIGN_NEAR_RELATIVE_BAND) * best[
                "relative_frobenius_residual"
            ]
            sign_pairs[f"{left_id}/{right_id}"] = {
                "best": best,
                "near_optimal_relative_signs": [
                    row for row in retained
                    if row["relative_frobenius_residual"] <= cutoff
                ],
                "next_distinct": second,
            }

    joint_cutoff = (1 + SIGN_NEAR_RELATIVE_BAND) ** 2 * pair_minimum_sum
    only_displayed_jointly_near = (
        pair_minimum_sum + minimum_nontrivial_increment > joint_cutoff
    )
    if not only_displayed_jointly_near:
        raise AssertionError("the pairwise gap no longer certifies the joint near-optimum")

    permutation = permutation_search(lineage, matrices, local_data)
    retained = permutation.pop("retained_per_fibre")
    if any(len(records) != 1 for records in retained.values()):
        raise AssertionError("more than one signed permutation survived the declared band")
    signed_permutation_candidates = []
    candidate = {}
    for curve_id in TARGETS:
        permutation_value, signs, residual = retained[curve_id][0]
        candidate[str(curve_id)] = {
            "permutation_new_label_to_old_point": [value + 1 for value in permutation_value],
            "signs": signs,
            "residual_against_356": residual,
        }
    signed_permutation_candidates.append(candidate)

    exact_point_replay = {}
    fibres_by_id = {
        int(fibre["curve_id"]): fibre
        for fibre in lineage["rootless_k3_interpolation_input"]["fibres"]
    }
    for curve_id in TARGETS:
        fibre = fibres_by_id[curve_id]
        A_value, B_value = map(Fraction, fibre["short_model"])
        transform = candidate[str(curve_id)]
        transformed_points = []
        for old_index, sign in zip(
            transform["permutation_new_label_to_old_point"], transform["signs"]
        ):
            x_value, y_value = map(
                Fraction, fibre["short_points_first_17"][old_index - 1]
            )
            if sign == "-":
                y_value = -y_value
            if y_value**2 != x_value**3 + A_value * x_value + B_value:
                raise AssertionError("exact signed-permuted point missed its short curve")
            transformed_points.append([str(x_value), str(y_value)])
        exact_point_replay[str(curve_id)] = {
            "verified_point_count": len(transformed_points),
            "transformed_points_sha256": sha256_bytes(
                json.dumps(transformed_points, separators=(",", ":")).encode()
            ),
        }

    payload = {
        "schema": "icarm.wgxli-rank17-signed-permutation-rebasing.v1",
        "status": "PASS_BOUNDED_SIGN_AND_PERMUTATION_ALIGNMENT",
        "inputs": {
            str(arguments.input.relative_to(ROOT)): sha256_bytes(lineage_raw),
            DATABASE_URL: database_source_hash,
        },
        "software": {"pari_gp": pari_version, "numpy": numpy.__version__},
        "sign_alignment": {
            "relative_signs_enumerated_per_pair": 1 << 16,
            "global_sign_gauge": (
                "the reference fibre is fixed to +...+; a common diagonal sign on all "
                "five fibres is independent global section inversion, and D and -D act "
                "identically in each pairwise Gram fit"
            ),
            "joint_objective": "sum over ten pairs of squared scale-fitted relative Frobenius residuals",
            "joint_minimum": pair_minimum_sum,
            "near_relative_band": SIGN_NEAR_RELATIVE_BAND,
            "minimum_increment_for_any_nontrivial_pair_sign": minimum_nontrivial_increment,
            "only_displayed_signs_jointly_near_optimal": only_displayed_jointly_near,
            "pair_records": sign_pairs,
        },
        "exact_local_component_fingerprints": {
            str(curve_id): record for curve_id, record in local_data.items()
        },
        "bounded_permutation_search": permutation,
        "retained_signed_permutation_candidates": signed_permutation_candidates,
        "exact_signed_permutation_group_law_replay": {
            "operation": "permutation followed by elliptic inversion (x,y)->(x,-y)",
            "curves": exact_point_replay,
            "all_85_transformed_points_verified": True,
        },
        "negative_controls": {
            "signed_residual_cutoff": 0.2,
            "components": [list(component) for component in components],
            "required_separate_components": [[363, 364, 378], [389, 390, 391]],
            "all_pair_best_signed_fits": all_signed_pair_fits,
            "status": "PASS_CONTROLS_REMAIN_SEPARATE",
        },
        "proof_boundary": (
            "The sign enumeration is exhaustive and its joint optimum follows because all "
            "ten pairwise minima are simultaneously attained. The permutation enumeration is "
            "complete only inside the explicitly declared fingerprint graph and numerical "
            "slacks. Gram alignment proposes rebasing; it neither constructs exact elliptic "
            "points nor proves common-family membership."
        ),
        "runtime_seconds": time.monotonic() - started,
    }
    if arguments.check:
        expected = json.loads(arguments.output.read_text())
        expected.pop("runtime_seconds", None)
        payload.pop("runtime_seconds", None)
        if expected != payload:
            raise SystemExit("stale signed-permutation rebasing artifact")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
        arguments.output.write_text(rendered)
        print(f"WGXLIREBASE|output={arguments.output}|sha256={sha256_bytes(rendered.encode())}")
    print(
        "WGXLIREBASE|pair_signs=65536x10|signed_candidates=1|"
        "negative_controls=PASS|status=PASS_BOUNDED_SIGN_AND_PERMUTATION_ALIGNMENT"
    )


if __name__ == "__main__":
    main()
