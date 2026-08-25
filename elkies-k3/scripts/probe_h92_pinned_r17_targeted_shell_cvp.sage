#!/usr/bin/env sage -python
"""Probe fixed-q pinned-R17 bisection shells by target-directed exact CVP.

This is a bounded candidate generator, not a complete shell enumeration.  For
``D=(q/2,2,w)`` and a marked isotropic target ``T=(A,B,t)``, minimizing
``D.T`` at fixed ``w^2=2q`` amounts to maximizing ``w.M.t``.  The latter is
probed by enumerating close lattice vectors around rational points on the ray
``sqrt(2q/(t.M.t))*t``.  Every retained candidate is then checked exactly for
norm, primitivity, and all-section nefness.
"""

import argparse
import hashlib
import json
from itertools import islice
from pathlib import Path

from sage.all import *
from sage.modules.free_quadratic_module_integer_symmetric import IntegralLattice


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
PINNED_FRAME = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
MARKING = GENERATED / "elkies-k3-h3-pinned-r17-equation-marking.json"
Q4_RANKING = GENERATED / "elkies-k3-h3-pinned-r17-q4-degree2-targeted-ranking.json"
DEFAULT_OUTPUT = GENERATED / "elkies-k3-h3-pinned-r17-targeted-shell-cvp.json"
U2 = matrix(ZZ, ((0, 1), (1, 0)))


def load_matrix(path):
    return matrix(
        ZZ,
        [[ZZ(value) for value in line.split()] for line in path.read_text().splitlines()
         if line.strip() and not line.lstrip().startswith("#")],
    )


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--target", default="orbit12")
parser.add_argument("--frame", type=Path, default=PINNED_FRAME)
parser.add_argument("--marking", type=Path, default=MARKING)
parser.add_argument("--root-rank", type=int, default=0)
parser.add_argument("--q", type=int, action="append", dest="q_values")
parser.add_argument("--degree", type=int, default=2, choices=(1, 2, 3, 4))
parser.add_argument("--scale-min", type=int, default=20,
                    help="minimum percent of the real Lagrange scale")
parser.add_argument("--scale-max", type=int, default=180,
                    help="maximum percent of the real Lagrange scale")
parser.add_argument("--scale-step", type=int, default=2)
parser.add_argument("--close-count", type=int, default=128)
parser.add_argument("--retain", type=int, default=200)
parser.add_argument("--analyze-child-top", type=int, default=0)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
args = parser.parse_args()
q_values = tuple(args.q_values or (4, 6, 8))
assert all(q > 0 and q % args.degree == 0 for q in q_values)
assert 0 < args.scale_min <= args.scale_max and args.scale_step > 0
assert args.close_count > 0 and args.retain > 0
assert args.analyze_child_top >= 0
OUTPUT = args.output.resolve()
FRAME = args.frame.resolve()
MARKING_INPUT = args.marking.resolve()

pinned = load_matrix(FRAME)
g = block_diagonal_matrix(U2, -pinned)
marking = json.loads(MARKING_INPUT.read_text())
q4_ranking = json.loads(Q4_RANKING.read_text())
assert q4_ranking["status"] == "PASS_EXACT_PINNED_R17_Q4_TARGETED_RANKING"
assert 0 <= args.root_rank <= 17
target_field = (
    "target_fibres_in_root_adapted_hub"
    if "target_fibres_in_root_adapted_hub" in marking
    else "target_fibres_in_child"
)
targets = {
    name: vector(ZZ, value)
    for name, value in marking[target_field].items()
}
assert args.target in targets
target = targets[args.target]
target_tail = vector(ZZ, target[2:])
target_tail_norm = ZZ(target_tail * pinned * target_tail)
assert target_tail_norm == 2 * target[0] * target[1] and target_tail_norm > 0
lattice = IntegralLattice(pinned)
real = RealField(200)
rational_denominator = ZZ(10) ** 40


def rational_approximation(value):
    return QQ(ZZ((value * rational_denominator).round())) / rational_denominator


def bezout_vector_for_pairing(ns, fibre):
    current = ZZ(0)
    result = [ZZ(0)] * ns.nrows()
    for index, value in enumerate(ns * fibre):
        if value == 0:
            continue
        divisor, left, right = xgcd(current, ZZ(value))
        result = [left * entry for entry in result]
        result[index] += right
        current = divisor
    assert abs(current) == 1
    if current == -1:
        result = [-entry for entry in result]
    return vector(ZZ, result)


def child_root_data(fibre):
    mate = bezout_vector_for_pairing(g, fibre)
    mate -= ZZ(mate * g * mate // 2) * fibre
    kernel = matrix(ZZ, [list(fibre * g), list(mate * g)]).right_kernel_matrix()
    child = -(kernel * g * kernel.transpose())
    transition = matrix(ZZ, [list(fibre), list(mate)] + list(kernel.rows()))
    assert abs(transition.det()) == 1
    minimum = pari(child).qfminim(2)
    root_count = ZZ(minimum[0])
    if root_count == 0:
        return [0, 0, 1]
    half = [vector(ZZ, column) for column in matrix(ZZ, minimum[2]).columns()]
    roots = half + [-root for root in half]
    basis = matrix(ZZ, [list(root) for root in roots]).row_module().basis_matrix()
    root_gram = basis * child * basis.transpose()
    return [int(basis.rank()), int(root_count), int(abs(root_gram.det()))]


def components(cartan):
    unseen = set(range(cartan.nrows()))
    result = []
    while unseen:
        todo = [unseen.pop()]
        component = []
        while todo:
            node = todo.pop()
            component.append(node)
            for other in tuple(unseen):
                if cartan[node, other]:
                    unseen.remove(other)
                    todo.append(other)
        result.append(tuple(sorted(component)))
    return tuple(result)


def highest_roots(cartan):
    if cartan.nrows() == 0:
        return ()
    half = matrix(ZZ, pari(cartan).qfminim(2)[2]).transpose().rows()
    roots = tuple(half) + tuple(-item for item in half)
    answer = []
    for component in components(cartan):
        candidates = [
            root for root in roots
            if all(value >= 0 for value in root)
            and all(index in component or root[index] == 0
                    for index in range(cartan.nrows()))
        ]
        answer.append(max(candidates, key=lambda root: sum(root)))
    return tuple(answer)


cartan = pinned[:args.root_rank, :args.root_rank]
tops = highest_roots(cartan)


def negative_horizontal_walls(w, degree):
    """Exactly enumerate every negative old-horizontal (-2)-curve wall."""
    walls = []
    w_norm = ZZ(w * pinned * w)
    for old_degree in range(1, int(degree) + 1):
        m = ZZ(old_degree)
        cross = -degree * m * pinned * w.column()
        augmented = block_matrix(ZZ, [
            [degree**2 * pinned, cross],
            [cross.transpose(), matrix(ZZ, [[m**2 * w_norm + 1]])],
        ])
        # For z=1 the augmented norm is ||degree*x-m*w||^2+1.
        # Negative intersection is equivalent to the first term < 2*degree^2.
        result = pari(augmented).qfminim(2 * degree**2 - 1)
        half = matrix(ZZ, result[2]).transpose().rows()
        normalized = set()
        for candidate in half:
            if abs(candidate[-1]) != 1:
                continue
            value = candidate if candidate[-1] == 1 else -candidate
            normalized.add(tuple(value))
        for value in normalized:
            x = vector(ZZ, value[:-1])
            x_norm = ZZ(x * pinned * x)
            if (x_norm - 2) % (2 * m):
                continue
            k = ZZ((x_norm - 2) // (2 * m))
            intersection = ZZ((w_norm // (2 * degree)) * m + degree * k - w * pinned * x)
            y = degree * x - m * w
            assert y * pinned * y == 2 * degree * m * intersection + 2 * degree**2
            if intersection < 0:
                walls.append({
                    "old_fibre_degree": int(m),
                    "curve": [int(k), int(m)] + list(map(int, x)),
                    "intersection": int(intersection),
                    "congruence_vector_norm": int(y * pinned * y),
                })
    return sorted(walls, key=lambda item: (item["intersection"], item["curve"]))


def candidate_record(q_value, w, scale_percent):
    degree = ZZ(args.degree)
    a = ZZ(q_value // degree)
    fibre = vector(ZZ, [a, degree] + list(w))
    assert fibre * g * fibre == 0
    if gcd(tuple(g * fibre)) != 1:
        return None
    labels = w * pinned[:, :args.root_rank]
    affine = [ZZ(degree - top * labels) for top in tops]
    if min(tuple(labels) + tuple(affine) + (ZZ(0),)) < 0:
        return None
    center = vector(QQ, w) / degree
    closest = vector(ZZ, next(lattice.enumerate_close_vectors(center)))
    distance = (closest - center) * pinned * (closest - center)
    minimum_section = degree * (distance - 2) / 2
    if minimum_section < 0:
        return None
    horizontal_walls = negative_horizontal_walls(w, degree)
    if horizontal_walls:
        return None
    degrees = {
        name: int(fibre * g * value) for name, value in targets.items()
    }
    return {
        "candidate_id": {
            "q": int(q_value),
            "old_fibre_degree": int(degree),
            "discovery_scale_percent": int(scale_percent),
        },
        "fibre_in_pinned_R17": list(map(int, fibre)),
        "P_dot_O": int(a - degree),
        "minimum_section_distance": str(distance),
        "minimum_section_intersection": str(minimum_section),
        "exact_negative_horizontal_walls": horizontal_walls,
        "exact_horizontal_nef_gate": True,
        "component_pairings": list(map(int, labels)),
        "affine_component_pairings": list(map(int, affine)),
        "marked_target_degrees": degrees,
        "coordinate_growth_max": int(max(abs(value) for value in fibre)),
    }


searches = {}
for q_value in q_values:
    optimum = sqrt(real(2 * q_value) / real(target_tail_norm))
    seen = set()
    accepted = {}
    enumerated = 0
    for scale_percent in range(args.scale_min, args.scale_max + 1, args.scale_step):
        scalar = rational_approximation(optimum * real(scale_percent) / 100)
        center = scalar * target_tail
        for close in islice(lattice.enumerate_close_vectors(center), args.close_count):
            enumerated += 1
            w = vector(ZZ, close)
            key = tuple(w)
            if key in seen:
                continue
            seen.add(key)
            if w * pinned * w != 2 * q_value:
                continue
            record = candidate_record(q_value, w, scale_percent)
            if record is not None:
                accepted[key] = record
    ranking = sorted(
        accepted.values(),
        key=lambda item: (
            item["marked_target_degrees"][args.target],
            item["coordinate_growth_max"],
            tuple(item["fibre_in_pinned_R17"]),
        ),
    )
    for record in ranking[:args.analyze_child_top]:
        record["child_root_data"] = child_root_data(
            vector(ZZ, record["fibre_in_pinned_R17"])
        )
        record["child_mw_rank"] = 17 - record["child_root_data"][0]
    searches[str(q_value)] = {
        "q": int(q_value),
        "old_fibre_degree": args.degree,
        "lagrange_scale_200bit": str(optimum),
        "close_vectors_enumerated_with_repetitions": enumerated,
        "distinct_vectors_seen": len(seen),
        "primitive_nef_fixed_norm_candidates": len(accepted),
        "rankings": ranking[:args.retain],
    }
    best = ranking[0] if ranking else None
    print(
        "R17CVP|target={}|q={}|seen={}|accepted={}|best_degree={}|best_fibre={}".format(
            args.target, q_value, len(seen), len(accepted),
            None if best is None else best["marked_target_degrees"][args.target],
            None if best is None else ",".join(map(str, best["fibre_in_pinned_R17"])),
        ),
        flush=True,
    )

calibration_applicable = (
    args.root_rank == 0
    and FRAME == PINNED_FRAME.resolve()
    and args.degree == 2
    and args.target in q4_ranking["rankings_top_200"]
)
q4_exact = (
    q4_ranking["rankings_top_200"][args.target][0]
    if calibration_applicable else None
)
q4_discovered = searches.get("4", {}).get("rankings", [])
q4_calibration = {
    "applicable": calibration_applicable,
    "performed": calibration_applicable and "4" in searches,
    "exact_best_degree": (
        q4_exact["marked_target_degrees"][args.target] if q4_exact else None
    ),
    "discovered_best_degree": (
        q4_discovered[0]["marked_target_degrees"][args.target]
        if q4_discovered else None
    ),
}
q4_calibration["reproduces_exact_best_degree"] = (
    q4_calibration["performed"]
    and q4_calibration["discovered_best_degree"] == q4_calibration["exact_best_degree"]
)

inputs = tuple(dict.fromkeys((FRAME, MARKING_INPUT, Q4_RANKING)))
payload = {
    "schema": "elkies-k3.h3-root-adapted-targeted-shell-cvp.v1",
    "status": "PASS_EXACT_CANDIDATE_DISCOVERY_BOUNDED_CVP",
    "source_hub": marking.get("source_hub", marking.get("hub", "marked_child")),
    "frame": str(FRAME.relative_to(ROOT)),
    "root_rank": args.root_rank,
    "target": args.target,
    "search_parameters": {
        "q_values": list(q_values),
        "old_fibre_degree": args.degree,
        "scale_percent_min": args.scale_min,
        "scale_percent_max": args.scale_max,
        "scale_percent_step": args.scale_step,
        "close_vectors_per_scale": args.close_count,
        "retained_per_q": args.retain,
        "child_root_data_analyzed_per_q": args.analyze_child_top,
        "rational_approximation_denominator": str(rational_denominator),
    },
    "q4_exact_calibration": q4_calibration,
    "searches": searches,
    "proof_boundary": (
        "Every retained vector has exact fixed norm and yields an exact primitive "
        "isotropic class. Simple and affine component pairings are nonnegative. "
        "All old-fibre sections are checked by exact closest-vector "
        "enumeration. Every possible negative old-horizontal (-2)-curve of degree "
        "at most the candidate degree is excluded by an exact augmented-lattice "
        "short-vector enumeration. The ray/scale sample is bounded and does not "
        "prove target optimality or shell exhaustiveness. Child roots, chambers, and "
        "full transports require the separate candidate certifier."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in inputs],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in inputs
        },
    },
}
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "R17CVP|q4_calibrated={}|status={}|output={}".format(
        int(q4_calibration["reproduces_exact_best_degree"]), payload["status"], OUTPUT
    ),
    flush=True,
)
