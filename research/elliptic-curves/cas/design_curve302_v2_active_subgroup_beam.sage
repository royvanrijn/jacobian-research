#!/usr/bin/env sage-python
"""Geometry-only, basis-invariant active-subgroup beam design for V2.

No point search is imported or run.  The only prospective inputs are the V2
generic shell landscape, the current certified subgroup, and a hypothetical
batch of independently certified newly discovered points.  The fixed batch in
this calibration is made from three already-V2-discovered points only to test
the transition rule offline; its identities are not used as score labels.

The score sees every extension parity of every retained generic anchor, exact
CVP minima, their multiplicity, and exact reduced-chart coefficient sizes.
It has no dependency on the 302 exceptional-direction diagnostic.  A separate
evaluator may read that diagnostic after this program has sealed its result.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.machinery
import json
import math
import random
from fractions import Fraction as F
from pathlib import Path

import numpy as np
from sage.all import RealField, ZZ, identity_matrix, matrix, pari, vector


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
V2 = ROOT / "artifacts/local/elliptic-curves/adaptive-visibility-cascade-v2"
CALIBRATION = V2 / "calibration302"
BASE_SELECTION = CALIBRATION / "epoch-07/selection.json"  # certified rank 24 before its gain
BATCH_AUDITS = (
    CALIBRATION / "epoch-07/mod2-007.json",  # 24 -> 25
    CALIBRATION / "epoch-08/mod2-040.json",  # 25 -> 26
    CALIBRATION / "epoch-09/mod2-029.json",  # 26 -> 27
)
V2_PROTOCOL = V2 / "protocol.json"
V2_GENERIC = V2 / "generic.json"
V2_METRIC_REPLAY = V2 / "metric-replay.json"
OUTPUT = ROOT / "artifacts/generated-results/elliptic-curves/curve302_v2_active_subgroup_beam_v1.json"
ANCHORS_PER_SHELL = 4
BEAM_WIDTH = 4
NODE_LIMIT = 2_000_000
SCALE = 1_000_000
FORBIDDEN_SELECTION_TOKENS = (
    "exceptional_subgroup_landscape",
    "residual_visibility",
    "residual-strict",
    "curve302_visibility",
)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def read(path: Path) -> dict:
    return json.loads(path.read_text())


def load(name: str):
    return importlib.machinery.SourceFileLoader("v2_beam_" + name.replace(".", "_"), str(CAS / name)).load_module()


def qpoint(row) -> tuple[F, F]:
    return tuple(F(value) for value in row)


def canonical_point_key(point) -> tuple[str, str]:
    x, y = point
    return str(x), str(abs(y))


def key_digest(value) -> str:
    return hashlib.sha256(json.dumps(value, separators=(",", ":"), sort_keys=True).encode()).hexdigest()


def qquantile(values: list[float], numerator: int, denominator: int) -> float:
    if not values:
        raise ArithmeticError("empty score distribution")
    ordered = sorted(values)
    index = (len(ordered) - 1) * numerator // denominator
    return ordered[index]


def chart_bits(mapping) -> int:
    coefficients = [F(value) for value in mapping["discriminant_quartic"]]
    return max(max(abs(value.numerator).bit_length(), value.denominator.bit_length()) for value in coefficients)


def rounded_gram(geometry, model, basis):
    height, asymmetry = geometry.canonical_height_gram(model, basis)
    gram = matrix(ZZ, geometry.rounded_gram(height, SCALE))
    if not gram.is_positive_definite():
        raise ArithmeticError("rounded active-subgroup metric is not positive definite")
    return gram, str(asymmetry)


def logarithmic_geometric_mean(gram) -> float:
    field = RealField(256)
    return math.exp(float(field(abs(gram.det())).log() / gram.nrows()))


def finite_fingerprint_columns(model, basis):
    """Canonical finite-reduction labels of actual mod-two cosets.

    The reduction tables are evaluated on actual rational points.  Changing an
    input basis changes the coordinate word but not the final fingerprint of
    the represented point.
    """

    mod2 = load("audit_recorded_point_mod2_rank_v3.py")
    from mod2_reduction_independence import _primes_up_to
    from research_runtime.finite_reduction import ReductionCache
    from research_runtime.memory_store import MemoryFactStore
    cache = ReductionCache(MemoryFactStore())
    columns = [0] * len(basis)
    pivots = {}
    offset = 0
    for prime in _primes_up_to(1000):
        if prime == 2:
            continue
        try:
            signature = mod2.signature(cache, model, basis, prime)
        except ValueError:
            continue
        for row in signature.rows:
            mod2.insert(pivots, row)
            for index, bit in enumerate(row):
                columns[index] |= int(bit) << offset
            offset += 1
    if len(pivots) != len(basis):
        raise ArithmeticError("finite fingerprints are not injective on this active subgroup")
    return columns


def fingerprint(word, columns) -> int:
    answer = 0
    for coefficient, column in zip(word, columns):
        if int(coefficient) & 1:
            answer ^= int(column)
    return answer


def transform_word(word, transport):
    result = vector(ZZ, word) * transport
    if any(value.denominator() != 1 for value in result):
        raise ArithmeticError("unimodular rebase gave a nonintegral representative word")
    return [int(value) for value in result]


def prepare_exact_cvp(gram, exact_parity):
    change = matrix(ZZ, pari(gram).qflllgram()).transpose()
    if abs(change.det()) != 1:
        raise ArithmeticError("LLL transport is not unimodular")
    inverse = change.inverse()
    reduced = change * gram * change.transpose()
    return change, inverse, exact_parity(reduced.rows())


def coset_minimum(model, basis, columns, exact, change, inverse, source_word, group, mapper, chart_cache, map_charts):
    """Exact shortest representatives of one actual degree-two parity class."""

    reduced_residue = vector(ZZ, source_word) * inverse
    parity = [int(value) & 1 for value in reduced_residue]
    seeds, _ = exact.babai(np.asarray([parity], dtype=np.int64))
    proof = exact.solve(parity, seeds[0], NODE_LIMIT)
    centres = []
    bits = []
    for reduced_word in proof["minima"]:
        word = [int(value) for value in vector(ZZ, reduced_word) * change]
        if [value & 1 for value in word] != [value & 1 for value in source_word]:
            raise ArithmeticError("exact CVP returned the wrong original parity")
        point = group.linear_combination(model, basis, word)
        if point is None:
            raise ArithmeticError("CVP returned the point at infinity as a centre")
        key = canonical_point_key(point)
        centres.append(key)
        if key not in chart_cache:
            if not map_charts:
                raise ArithmeticError("rebased CVP produced an unseen centre")
            chart_cache[key] = chart_bits(mapper.mapping(model, basis, {"representative": word}))
        bits.append(chart_cache[key])
    return {
        "coset_fingerprint": str(fingerprint(source_word, columns)),
        "exact_norm": int(proof["norm"]),
        "multiplicity": len(proof["minima"]),
        "cvp_nodes": int(proof["nodes"]),
        "chart_bits_min": min(bits),
        "centre_key_digest": key_digest(sorted(set(centres))),
    }


def score_state(model, basis, anchors, word_transport, geometry, exact_parity, group, mapper, chart_cache, thresholds=None, map_charts=True, gram_override=None):
    """Score one active subgroup only from V2 generic-orbit and own-point data."""

    rank = len(basis)
    # Entrywise rounding is not covariant: recomputing a decimal height matrix
    # in a new basis and then rounding can change an integral entry by one.  A
    # certified state therefore freezes its rounded form once; every rebase
    # transports that same finite decision metric by U*G*U^T.
    if gram_override is None:
        gram, asymmetry = rounded_gram(geometry, model, basis)
        metric_origin = "fresh 384-bit pairing, rounded entrywise at 1e6"
    else:
        gram = matrix(ZZ, gram_override)
        if not gram.is_positive_definite():
            raise ArithmeticError("transported rounded metric is not positive definite")
        asymmetry, metric_origin = "not recomputed under rebase", "exact unimodular transport of the frozen rounded metric"
    change, inverse, exact = prepare_exact_cvp(gram, exact_parity)
    columns = finite_fingerprint_columns(model, basis)
    records = []
    for anchor in anchors:
        anchor_word = list(anchor["representative"]) + [0] * (rank - 24)
        extension_dimension = rank - 24
        for mask in range(1 << extension_dimension):
            source_word = anchor_word[:]
            for index in range(extension_dimension):
                source_word[24 + index] = (mask >> index) & 1
            source_word = transform_word(source_word, word_transport)
            record = coset_minimum(
                model, basis, columns, exact, change, inverse, source_word,
                group, mapper, chart_cache, map_charts,
            )
            record["generic_shell"] = int(anchor["shell"])
            record["anchor_fingerprint"] = str(anchor["fingerprint"])
            records.append(record)
    normalizer = logarithmic_geometric_mean(gram)
    normalized = [row["exact_norm"] / normalizer for row in records]
    chart_values = [row["chart_bits_min"] for row in records]
    multiplicities = [row["multiplicity"] for row in records]
    if thresholds is None:
        thresholds = {
            "lower_quartile": qquantile(normalized, 1, 4),
            "median": qquantile(normalized, 1, 2),
            "chart_median": qquantile(chart_values, 1, 2),
        }
    low_one = sum(value <= thresholds["lower_quartile"] for value in normalized)
    low_two = sum(value <= thresholds["median"] for value in normalized)
    low_chart = sum(value <= thresholds["chart_median"] for value in chart_values)
    # The ordering is deliberately lexicographic: no uncalibrated numerical
    # trade-off turns a chart-bit reduction into an arbitrary number of CVP
    # shell hits.  Every component is an invariant of actual parity classes.
    score_key = (
        -low_one / len(records),
        -low_two / len(records),
        sum(normalized) / len(normalized),
        -sum(math.log2(value) for value in multiplicities) / len(multiplicities),
        qquantile(chart_values, 1, 2),
        qquantile(chart_values, 9, 10),
        tuple(sorted(row["coset_fingerprint"] for row in records)),
    )
    signature = [
        (row["coset_fingerprint"], row["exact_norm"], row["multiplicity"], row["chart_bits_min"], row["centre_key_digest"])
        for row in sorted(records, key=lambda row: int(row["coset_fingerprint"]))
    ]
    return {
        "rank": rank,
        "coset_count": len(records),
        "metric": {
            "rounded_gram_sha256": key_digest([list(map(int, row)) for row in gram.rows()]),
            "height_asymmetry": asymmetry,
            "metric_origin": metric_origin,
            "geometric_mean_norm": f"{normalizer:.12g}",
        },
        "rounded_gram": [list(map(int, row)) for row in gram.rows()],
        "profile": {
            "low_shell_thresholds": {key: f"{value:.12g}" for key, value in thresholds.items()},
            "low_shell_counts": {"lower_quartile": low_one, "median": low_two},
            "low_shell_rates": {"lower_quartile": f"{low_one/len(records):.12g}", "median": f"{low_two/len(records):.12g}"},
            "exact_CVP_normalized_minimum_quantiles": {
                "q10": f"{qquantile(normalized, 1, 10):.12g}", "median": f"{qquantile(normalized, 1, 2):.12g}", "q90": f"{qquantile(normalized, 9, 10):.12g}",
                "mean": f"{sum(normalized)/len(normalized):.12g}",
            },
            "short_representative_multiplicity": {
                "mean_log2": f"{sum(math.log2(value) for value in multiplicities)/len(multiplicities):.12g}",
                "maximum": max(multiplicities), "total": sum(multiplicities),
            },
            "reduced_chart_complexity_bits": {
                "at_or_below_baseline_median": low_chart,
                "median": qquantile(chart_values, 1, 2), "q90": qquantile(chart_values, 9, 10), "maximum": max(chart_values),
            },
        },
        "selection_key": [
            f"{score_key[0]:.12g}", f"{score_key[1]:.12g}", f"{score_key[2]:.12g}", f"{score_key[3]:.12g}", score_key[4], score_key[5],
            key_digest(score_key[6]),
        ],
        "comparison_key": score_key,
        "invariant_record_signature": key_digest(signature),
        "records": records,
    }


def random_unimodular(rank: int, seed: int):
    rng = random.Random(seed)
    transform = identity_matrix(ZZ, rank)
    for _ in range(4 * rank + 5):
        left, right = rng.sample(range(rank), 2)
        choice = rng.randrange(3)
        if choice == 0:
            transform.swap_rows(left, right)
        elif choice == 1:
            transform.rescale_row(left, -1)
        else:
            transform.add_multiple_of_row(left, right, rng.choice((-2, -1, 1, 2)))
    if abs(transform.det()) != 1:
        raise ArithmeticError("random elementary transform stopped being unimodular")
    return transform


def rebase_basis(model, basis, transform, group):
    rebased = []
    for row in transform.rows():
        point = group.linear_combination(model, basis, [int(value) for value in row])
        if point is None:
            raise ArithmeticError("unimodular basis row evaluated to infinity")
        rebased.append(point)
    return tuple(rebased)


def same_score(left, right):
    return left["invariant_record_signature"] == right["invariant_record_signature"] and left["selection_key"] == right["selection_key"]


def exact_batch():
    selection = read(BASE_SELECTION)
    if selection["rank"] != 24 or len(selection["basis"]) != 24:
        raise ArithmeticError("chosen V2 base selection is not rank 24")
    basis = tuple(qpoint(point) for point in selection["basis"])
    batch = []
    current = basis
    for audit_path in BATCH_AUDITS:
        audit = read(audit_path)
        if audit.get("status") != "COMPLETE_DECLARED_FINITE_AUDIT":
            raise ArithmeticError("V2 batch source lacks an exact finite rank certificate")
        points = tuple(qpoint(point) for point in audit["independent_points"])
        if points[:len(current)] != current or len(points) != len(current) + 1:
            raise ArithmeticError("V2 batch inputs are not a one-direction prefix chain")
        batch.append(points[-1])
        current = points
    return selection, basis, tuple(batch)


def selected_anchors(selection):
    groups = {8: [], 10: []}
    for row in selection["anchors"]:
        shell = int(row["anchor"]["shell"])
        if shell not in groups:
            raise ArithmeticError("unexpected generic shell in V2 landscape")
        groups[shell].append(row["anchor"])
    anchors = []
    for shell in (8, 10):
        group = sorted(groups[shell], key=lambda row: int(row["fingerprint"]))
        if len(group) < ANCHORS_PER_SHELL:
            raise ArithmeticError("V2 generic landscape has too few anchors")
        anchors += group[:ANCHORS_PER_SHELL]
    return anchors


def state_key(batch, mask):
    return tuple(sorted(canonical_point_key(batch[index]) for index in range(len(batch)) if (mask >> index) & 1))


def choose_beam(branches):
    """Keep diversity across all intermediate subgroup ranks before filling."""

    ordered = sorted(branches, key=lambda row: (row["comparison_key"], row["state_key"]))
    chosen = []
    for increase in sorted({row["increment_rank"] for row in ordered}):
        candidate = next(row for row in ordered if row["increment_rank"] == increase)
        chosen.append(candidate)
    for candidate in ordered:
        if len(chosen) == BEAM_WIDTH:
            break
        if candidate not in chosen:
            chosen.append(candidate)
    if len(chosen) > BEAM_WIDTH:
        raise ArithmeticError("beam diversity rule exceeded its fixed width")
    return chosen, ordered


def build() -> dict:
    protocol, selection = read(V2_PROTOCOL), read(BASE_SELECTION)
    if not protocol.get("calibration_only") or not protocol.get("oracle_boundary"):
        raise ArithmeticError("V2 generic-input boundary is absent")
    source_paths = (Path(__file__), V2_PROTOCOL, V2_GENERIC, V2_METRIC_REPLAY, BASE_SELECTION, *BATCH_AUDITS,
                    CAS / "visibility_lattice_v2.py", CAS / "factor_free_pari_mapping.sage", CAS / "prospective_half_lattice_v3.sage")
    if any(token in rel(path) for token in FORBIDDEN_SELECTION_TOKENS for path in source_paths):
        raise ArithmeticError("selection input list contains a forbidden diagnostic artifact")
    inputs = {rel(path): sha(path) for path in source_paths}
    geometry = load("prospective_half_lattice_v3.sage")
    exact_module = load("visibility_lattice_v2.py")
    group = load("half_lattice_pointed_sieve.py")
    mapper = load("factor_free_pari_mapping.sage")
    mapper.pari.allocatemem(256_000_000, silent=True)
    base_selection, base, batch = exact_batch()
    model = tuple(F(value) for value in read(BATCH_AUDITS[0])["curve"])
    anchors = selected_anchors(selection)
    chart_cache = {}
    identity_cache = identity_matrix(ZZ, 24)
    base_score = score_state(
        model, base, anchors, identity_cache, geometry, exact_module.ExactParity, group, mapper, chart_cache, thresholds=None,
    )
    thresholds = {key: float(value) for key, value in base_score["profile"]["low_shell_thresholds"].items()}
    branches = []
    for mask in range(1, 1 << len(batch)):
        indices = [index for index in range(len(batch)) if (mask >> index) & 1]
        basis = (*base, *(batch[index] for index in indices))
        score = score_state(
            model, basis, anchors, identity_matrix(ZZ, len(basis)), geometry, exact_module.ExactParity, group, mapper, chart_cache, thresholds,
        )
        score.update({
            "branch_mask": mask,
            "increment_rank": len(indices),
            "state_key": state_key(batch, mask),
            "member_point_keys": [canonical_point_key(batch[index]) for index in indices],
        })
        branches.append(score)
        print("V2 BEAM BRANCH", mask, "rank", len(basis), "cosets", score["coset_count"], flush=True)
    beam, ordered = choose_beam(branches)
    # Rebase every branch once, then the retained beam three times.  The exact
    # CVP/centre records and chart-bit-derived score must agree bit for bit.
    invariance = []
    requested = [(row, 1000 + row["branch_mask"]) for row in branches]
    for row in beam:
        requested += [(row, 2000 + row["branch_mask"]), (row, 3000 + row["branch_mask"])]
    for row, seed in requested:
        indices = [index for index in range(len(batch)) if (row["branch_mask"] >> index) & 1]
        basis = (*base, *(batch[index] for index in indices))
        transform = random_unimodular(len(basis), seed)
        rebased = rebase_basis(model, basis, transform, group)
        recomputed = score_state(
            model, rebased, anchors, transform.inverse(), geometry, exact_module.ExactParity, group, mapper, chart_cache, thresholds,
            map_charts=False, gram_override=transform * matrix(ZZ, row["rounded_gram"]) * transform.transpose(),
        )
        if not same_score(row, recomputed):
            raise ArithmeticError("basis rebase changed an intrinsic branch score")
        invariance.append({
            "branch_mask": row["branch_mask"], "seed": seed, "rank": len(basis), "determinant": int(transform.det()),
            "score_signature": recomputed["invariant_record_signature"], "status": "PASS_EXACT_CVP_AND_CHART_KEY_INVARIANCE",
        })
        print("V2 BEAM REBASE", row["branch_mask"], seed, flush=True)
    def render(row):
        return {key: value for key, value in row.items() if key not in {"comparison_key", "records", "rounded_gram"}}
    return {
        "schema": "elliptic-curves.curve302-v2-active-subgroup-beam.v1",
        "status": "PASS_GEOMETRY_ONLY_BASIS_INVARIANT_BEAM_DESIGN",
        "inputs": inputs,
        "prospective_boundary": {
            "allowed": "V2 frozen generic shell/orbit landscape; current certified subgroup; newly certified independent points from the current wave; exact height/CVP/chart computations on those data.",
            "forbidden": "exceptional directions, residual labels, diagnostic score values, catalogue ranks, point-search transcripts, and all point-search execution.",
            "selection_read_paths": sorted(inputs),
        },
        "synthetic_batch_calibration": {
            "meaning": "Three consecutively V2-certified points are coalesced retrospectively into one hypothetical rank-24 wave solely to exercise every nonzero intermediate subgroup. Their point identities are opaque state keys to the score.",
            "base_rank": 24, "batch_rank": 3, "nonzero_intermediate_subgroups": 7,
            "batch_point_keys": [canonical_point_key(point) for point in batch],
        },
        "anchor_policy": {
            "source": "V2 epoch-07 frozen generic shell landscape", "shells": [8, 10],
            "anchors_per_shell": ANCHORS_PER_SHELL,
            "selection": "Within each shell, retain the fixed generic classes with lexicographically least finite-reduction fingerprints. These are actual parity classes, not coordinate-prefix IDs.",
            "anchor_fingerprints": [str(row["fingerprint"]) for row in anchors],
        },
        "score_definition": {
            "cosets": "For every retained generic anchor and every binary extension in the candidate subgroup, certify the exact minimum in that M/2M parity coset by rational LDL ellipsoid enumeration.",
            "lexicographic_priority": [
                "larger proportion below the baseline lower-quartile normalized exact-CVP shell",
                "larger proportion below the baseline median normalized exact-CVP shell",
                "smaller mean determinant-normalized exact-CVP minimum",
                "larger mean log2 multiplicity of exact shortest representatives",
                "smaller median then 90th-percentile exact reduced-chart coefficient bits",
                "finite-reduction coset fingerprint tuple only as deterministic tie key",
            ],
            "normalization": "Each squared CVP norm is divided by det(G)^(1/r), so equal raw norms at different active ranks are not compared directly.",
            "finite_metric_transport": "At a certified subgroup, round the 384-bit height Gram once. Under any rebasing, use the exact congruent integer form U*G*U^T; do not independently re-round decimal height entries.",
            "beam": "Keep the best branch at each nonzero increment rank, then fill by the same lexicographic score to a hard maximum of four states.",
        },
        "baseline_profile": {key: value for key, value in base_score.items() if key not in {"comparison_key", "records"}},
        "all_branch_scores_in_selection_order": [render(row) for row in ordered],
        "beam": [render(row) for row in beam],
        "basis_invariance": {
            "random_unimodular_rebases": len(invariance),
            "coverage": "all seven candidate states once; each retained beam state two additional independent rebases",
            "checks": invariance,
            "criterion": "exact-CVP coset norms, multiplicities, actual-centre identities, derived reduced-chart bits, and lexicographic score agree exactly after transporting the same subgroup and cosets.",
        },
        "chart_cache_centre_count": len(chart_cache),
        "point_searches_run": 0,
        "reproducing_command": "sage -python elliptic-curves/cas/design_curve302_v2_active_subgroup_beam.sage --check",
        "boundary": "This freezes a branching policy design and an offline V2 batch calibration. It does not select a future search centre, execute a chart, recover a point, or prove any rank statement.",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.build == args.check:
        parser.error("choose exactly one of --build or --check")
    payload = build()
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.build:
        if OUTPUT.exists():
            raise FileExistsError("preserve the immutable beam-design artifact")
        OUTPUT.write_text(rendered)
    elif not OUTPUT.is_file() or OUTPUT.read_text() != rendered:
        raise ArithmeticError("beam-design replay differs")
    print("V2_ACTIVE_SUBGROUP_BEAM|states=7|beam<=4|point_searches=0|status=PASS", flush=True)


if __name__ == "__main__":
    main()
