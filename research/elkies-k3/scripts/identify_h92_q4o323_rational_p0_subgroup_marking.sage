#!/usr/bin/env sage-python
"""Embed the exact rational q323 P.O=0 height graph in the marked MW lattice.

The observed data are the exact QQ height Gram and the three component
profiles of ten polynomial sections.  Exhaust all global component
orientations and graph embeddings into the 258 marked lattice P.O=0 classes,
then test whether the marked q207 MW tail lies in the rational span and in the
integral subgroup.  This is finite exact lattice arithmetic only.
"""

import json
from itertools import product
from pathlib import Path

from sage.all import QQ, ZZ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
SUBGROUP = LOCAL / "q4o323-rational-p0-subgroup-qq.json"
ANCHOR = LOCAL / "q4o323-p0-shell-anchor-domains-mod61.json"
MARKING = LOCAL / "q4o323-reflected-fixed-suffix-component2-marking.json"
OUTPUT = LOCAL / "q4o323-rational-p0-subgroup-marking.json"

subgroup = json.loads(SUBGROUP.read_text())
anchor = json.loads(ANCHOR.read_text())
marking = json.loads(MARKING.read_text())
assert subgroup["status"] == "PASS_EXACT_QQ_Q4O323_RATIONAL_P0_SUBGROUP_HEIGHT_GRAM"
assert anchor["status"] == "PASS_MOD61_Q4O323_REGULAR_P0_SHELL_ANCHOR_DOMAINS"
assert marking["status"] == "PASS_EXACT_Q4O323_REFLECTED_FIXED_SUFFIX_MARKING"

frame = matrix(ZZ, [
    list(map(int, line.split()))
    for line in (ROOT/marking["frame_output"]).read_text().splitlines()
    if line.strip() and not line.lstrip().startswith("#")
])
root = frame[:7, :7]
coupling = frame[:7, 7:]
tail = frame[7:, 7:]
height = tail-coupling.transpose()*root.inverse()*coupling

lattice = anchor["lattice"]["sections"]
lattice_tails = [vector(ZZ, record["NS_coordinates"][-10:]) for record in lattice]
lattice_profiles = [tuple(record["component_profile_I4_I3_I3"]) for record in lattice]
lattice_heights = [value*height*value for value in lattice_tails]

representatives = list(map(int, subgroup["representative_shell_indices"]))
raw_profiles = [
    tuple(subgroup["raw_component_profiles_I3_I3_I4"][str(index)])
    for index in representatives
]
observed = matrix(QQ, subgroup["height_gram"])
assert observed.rank() == subgroup["height_gram_rank"] == 8

q207 = marking["fixed_suffix_transport"]["q207_component_reduction"]["equation_preflight"]
q207_tail = vector(ZZ, q207["horizontal_section"][-10:])
assert q207_tail*height*q207_tail == QQ(q207["horizontal_height"])


def mapped_profile(raw, finite_swap, signs):
    finite = (raw[1], raw[0]) if finite_swap else (raw[0], raw[1])
    return (
        signs[0]*raw[2] % 4,
        signs[1]*finite[0] % 3,
        signs[2]*finite[1] % 3,
    )


solutions = []
orientation_summaries = []
for finite_swap in (False, True):
    for signs in product((1, -1), repeat=3):
        pools = []
        for position, raw in enumerate(raw_profiles):
            profile = mapped_profile(raw, finite_swap, signs)
            pools.append([
                index for index in range(len(lattice))
                if lattice_profiles[index] == profile
                and lattice_heights[index] == observed[position, position]
            ])
        order = sorted(range(len(representatives)), key=lambda position: len(pools[position]))
        assigned = {}
        used = set()
        count_before = len(solutions)

        def search(depth):
            if depth == len(order):
                rows = [lattice_tails[assigned[position]] for position in range(len(representatives))]
                generator_matrix = matrix(ZZ, rows)
                rational_member = matrix(QQ, rows).stack(matrix(QQ, [q207_tail])).rank() == 8
                integral_member = bool(q207_tail in generator_matrix.row_module())
                solutions.append({
                    "finite_I3_swap": finite_swap,
                    "signs_I4_I3_I3": list(signs),
                    "shell_to_lattice": [
                        {
                            "shell_index": representatives[position],
                            "lattice_index": assigned[position],
                        }
                        for position in range(len(representatives))
                    ],
                    "generator_rank": int(generator_matrix.rank()),
                    "q207_in_rational_span": rational_member,
                    "q207_in_integral_subgroup": integral_member,
                })
                return
            position = order[depth]
            for candidate in pools[position]:
                if candidate in used:
                    continue
                if any(
                    lattice_tails[candidate]*height*lattice_tails[assigned[other]]
                    != observed[position, other]
                    for other in assigned
                ):
                    continue
                assigned[position] = candidate
                used.add(candidate)
                search(depth+1)
                used.remove(candidate)
                del assigned[position]

        search(0)
        orientation_summaries.append({
            "finite_I3_swap": finite_swap,
            "signs_I4_I3_I3": list(signs),
            "pool_sizes": list(map(len, pools)),
            "embedding_count": len(solutions)-count_before,
        })

rational_hits = [record for record in solutions if record["q207_in_rational_span"]]
integral_hits = [record for record in solutions if record["q207_in_integral_subgroup"]]
payload = {
    "schema": "elkies-k3.h92-q4o323-rational-p0-subgroup-marking.v1",
    "status": "PASS_EXACT_Q4O323_RATIONAL_P0_SUBGROUP_MARKING_AUDIT",
    "embedding_count": len(solutions),
    "q207_rational_span_embedding_count": len(rational_hits),
    "q207_integral_subgroup_embedding_count": len(integral_hits),
    "orientation_summaries": orientation_summaries,
    "embeddings": solutions,
    "q207_MW_tail": list(map(int, q207_tail)),
    "method": {
        "exact_height_graph_embedding": True,
        "global_component_orientations_tested": 16,
        "large_Groebner_required": False,
        "equation_elimination_required": False,
    },
    "proof_boundary": (
        "This identifies every lattice embedding compatible with the exact rational "
        "height graph and modular component profiles. Any residual embeddings require "
        "one more exact section/component orientation anchor before a unique group-law word."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
print(
    "Q4O323P0MARK|embeddings={}|q207_Qspan={}|q207_Zspan={}|output={}".format(
        len(solutions), len(rational_hits), len(integral_hits), OUTPUT,
    ), flush=True,
)
