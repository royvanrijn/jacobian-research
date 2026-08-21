#!/usr/bin/env sage
"""Bounded q4 shell-window search for an easier fifth q80 neighbor.

The pinned fifth step has horizontal class ``O+(-R)`` on the compact fourth
child, but its geometric zero pulls back to old-fiber degree 18.  This search
enumerates an indexed window in the norm-eight shell of the generic
``4A1/MW13`` frame, retains q4 neighbors with root rank at most one,
transports each class to the CM24 fourth child, and scores its
chamber-reduced pullback to the explicit
third-child Weierstrass surface.  It also records whether the horizontal MW
class is a sum of two of the five known degree-one sections.

This is a finite discovery calculation, not a nefness or equation certificate.
"""

import argparse
import json
import sys
from pathlib import Path

from sage.all import QQ, ZZ, matrix, pari, vector


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
search_parser = argparse.ArgumentParser(description=__doc__)
search_parser.add_argument("--start-half", type=int, default=13_500)
search_parser.add_argument("--stop-half", type=int, default=14_500)
search_parser.add_argument("--a", type=int, default=2)
search_parser.add_argument("--b", type=int, default=2)
search_parser.add_argument(
    "--mw-coordinates",
    default=None,
    help="optional comma-separated exact horizontal MW coordinates",
)
search_parser.add_argument(
    "--mw-up-to-sign",
    action="store_true",
    help="accept the negative of --mw-coordinates as well",
)
search_parser.add_argument(
    "--target-cm-roots",
    default=None,
    help="optional comma-separated (rank,count,determinant) filter",
)
search_parser.add_argument("--stop-after-target", action="store_true")
search_parser.add_argument(
    "--output",
    type=Path,
    default=None,
    help="optional local JSON output path (useful for disjoint windows)",
)
search_arguments = search_parser.parse_args()
search_start_half = search_arguments.start_half
search_stop_half = search_arguments.stop_half
search_a = ZZ(search_arguments.a)
search_b = ZZ(search_arguments.b)
if search_a*search_b != 4:
    search_parser.error("this norm-eight q4 shell requires --a * --b = 4")
search_mw_coordinates = (
    None
    if search_arguments.mw_coordinates is None
    else tuple(QQ(value) for value in search_arguments.mw_coordinates.split(","))
)
search_target_cm_roots = (
    None
    if search_arguments.target_cm_roots is None
    else tuple(int(value) for value in search_arguments.target_cm_roots.split(","))
)
if search_target_cm_roots is not None and len(search_target_cm_roots) != 3:
    search_parser.error("--target-cm-roots requires rank,count,determinant")
# The loaded exact readiness checker has its own command-line parser.  Hide
# this wrapper's indexed-window options while loading it, then restore argv.
saved_argv = list(sys.argv)
try:
    sys.argv = [sys.argv[0]]
    load(str(HERE / "analyze_q80_fifth_q4_cm24_readiness.sage"))
finally:
    sys.argv = saved_argv


def child_root_data(frame):
    minimum = pari(frame).qfminim(2)
    count = ZZ(minimum[0])
    if not count:
        return (0, 0, 1)
    # A reduced rank-one root system is A1 and has exactly two roots.
    # Any larger count already fails the only filter used below, so avoid an
    # expensive HNF/rank computation for rejected children.
    if count == 2:
        return (1, 2, 2)
    return (2, int(count), 0)


known_sections = tuple(
    {
        "x_coefficients": list(map(int, row[0])),
        "old_mw_coordinates": list(map(int, row[1])),
        "new_mw_coordinates": vector(ZZ, row[-1]),
    }
    for row in candidate_rows
)
pair_sums = {}
for left in range(len(known_sections)):
    for right in range(left, len(known_sections)):
        total = tuple(
            known_sections[left]["new_mw_coordinates"]
            + known_sections[right]["new_mw_coordinates"]
        )
        pair_sums.setdefault(total, []).append((left, right))

fourth_child_frame = generic_fourth_frame
shell_result = pari(fourth_child_frame).qfminim(8)
shell_half = matrix(ZZ, shell_result[2]).transpose()
norm_eight_half = tuple(
    vector(ZZ, row)
    for row in shell_half.rows()
    if row*fourth_child_frame*row == 8
)
print(
    "Q80FIFTHQ4SEARCH|"
    f"stage=shell|norm8_half={len(norm_eight_half)}|"
    f"window={search_start_half}:{search_stop_half or len(norm_eight_half)}",
    flush=True,
)

rows = []
primitive_count = 0
low_root_count = 0
window = norm_eight_half[search_start_half:search_stop_half]
mw_coordinate_map = (
    fourth_embedding
    * special_fourth
    * new_projected_basis.transpose()
    * new_optimal_height.inverse()
)
pair_candidates = []
for representative in window:
    for sign in (1, -1):
        candidate_v = sign*representative
        coordinate_key = tuple(candidate_v*mw_coordinate_map)
        coordinate_matches = (
            search_mw_coordinates is not None
            and (
                coordinate_key == search_mw_coordinates
                or (
                    search_arguments.mw_up_to_sign
                    and coordinate_key
                    == tuple(-value for value in search_mw_coordinates)
                )
            )
        )
        if (
            (search_mw_coordinates is None and coordinate_key in pair_sums)
            or coordinate_matches
        ):
            pair_candidates.append((candidate_v, coordinate_key))
print(
    "Q80FIFTHQ4SEARCH|"
    f"stage=pair_prefilter|oriented={2*len(window)}|"
    f"matched={len(pair_candidates)}",
    flush=True,
)
target_hits = 0
for candidate_v, coordinate_key in pair_candidates:
        try:
            child, candidate_transport = neighbor(
                fourth_child_frame, ZZ(4), search_a, search_b, candidate_v
            )
        except AssertionError:
            continue
        primitive_count += 1
        if primitive_count % 250 == 0:
            print(
                "Q80FIFTHQ4SEARCH|"
                f"stage=classify|primitive={primitive_count}|"
                f"root_rank_le1={low_root_count}",
                flush=True,
            )
        roots = child_root_data(child)
        if roots[0] > 1:
            continue
        low_root_count += 1
        special_child, _ = enhance_neighbor(
            candidate_transport, fourth_embedding, special_fourth
        )
        special_roots = tuple(map(int, root_invariants(special_child)[:3]))
        if (
            search_target_cm_roots is not None
            and special_roots != search_target_cm_roots
        ):
            continue
        target_hits += 1

        special_candidate = vector(
            ZZ,
            [search_a, search_b] + list(candidate_v*fourth_embedding),
        )
        reduced_candidate, reflections = chamber_reduce(
            special_candidate, new_curves, new_ns
        )
        assert reduced_candidate*new_ns*reduced_candidate == 0
        old_pullback = vector(ZZ, reduced_candidate*special_fourth_basis)
        frame_part = vector(QQ, reduced_candidate[2:])
        projection = new_project_mw(frame_part)
        coordinates = new_projected_basis.solve_left(projection)
        assert coordinates*new_projected_basis == projection
        assert tuple(coordinates) == coordinate_key

        rows.append(
            {
                "v": list(map(int, candidate_v)),
                "generic_child_roots": list(roots),
                "generic_child_mw_rank": int(17-roots[0]),
                "cm24_child_roots": list(special_roots),
                "cm24_child_mw_rank": 18-special_roots[0],
                "cm24_reflection_count": len(reflections),
                "cm24_reduced": list(map(int, reduced_candidate)),
                "cm24_D_dot_O": int(
                    intersection(reduced_candidate, new_zero, new_ns)
                ),
                "cm24_mw_coordinates": [str(value) for value in coordinates],
                "known_section_pair_representations": [
                    list(pair) for pair in pair_sums.get(coordinate_key, ())
                ],
                "old_pullback": list(map(int, old_pullback)),
                "old_fiber_degree": int(
                    intersection(old_pullback, old_fiber, old_special_ns)
                ),
                "old_zero_pairing": int(
                    intersection(old_pullback, old_zero, old_special_ns)
                ),
            }
        )
        print(
            "Q80FIFTHQ4SEARCH|"
            f"stage=hit|hit={low_root_count}|roots={roots}|"
            f"cm24_roots={special_roots}|"
            f"old_degree={rows[-1]['old_fiber_degree']}|"
            f"section_pairs={tuple(map(tuple, rows[-1]['known_section_pair_representations']))}|"
            f"v={tuple(map(int, candidate_v))}",
            flush=True,
        )
        if search_arguments.stop_after_target:
            break

rows.sort(
    key=lambda row: (
        0 if row["known_section_pair_representations"] else 1,
        row["old_fiber_degree"],
        row["old_zero_pairing"],
        row["generic_child_roots"],
        row["v"],
    )
)
payload = {
    "schema": "q80-fifth-q4-low-degree-neighbor-search-gf73-v1",
    "status": "bounded_indexed_norm_eight_shell_window_experiment",
    "prime": 73,
    "neighbor_factorization": [int(search_a), int(search_b)],
    "norm_eight_half_shell": len(norm_eight_half),
    "half_shell_window": [
        search_start_half,
        search_stop_half or len(norm_eight_half),
    ],
    "oriented_pair_prefilter_matches": len(pair_candidates),
    "oriented_primitive_candidates": primitive_count,
    "generic_root_rank_at_most_one": low_root_count,
    "known_degree_one_sections": [
        {
            key: value
            for key, value in section.items()
            if key != "new_mw_coordinates"
        }
        | {
            "new_mw_coordinates": list(
                map(int, section["new_mw_coordinates"])
            )
        }
        for section in known_sections
    ],
    "candidates": rows,
    "rank_claim": None,
    "reproduce": (
        "sage elkies-k3/scripts/"
        "search_q80_fifth_q4_low_degree_neighbors.sage "
        f"--start-half {search_start_half} --a {search_a} --b {search_b}"
        + (
            f" --stop-half {search_stop_half}"
            if search_stop_half is not None else ""
        )
    ),
}
output = search_arguments.output or (
    ROOT / "artifacts/local/"
    "q80-fifth-q4-low-degree-neighbor-search-gf73-v1.json"
)
output.parent.mkdir(parents=True, exist_ok=True)
output.write_text(json.dumps(payload, indent=2, sort_keys=True, default=int)+"\n")

print(
    "Q80FIFTHQ4SEARCH|"
    f"norm8_half={len(norm_eight_half)}|"
    f"window={search_start_half}:{search_stop_half or len(norm_eight_half)}|"
    f"pair_matched={len(pair_candidates)}|"
    f"primitive={primitive_count}|"
    f"root_rank_le1={low_root_count}|target_hits={target_hits}|"
    f"retained={len(rows)}|"
    f"output={output}|status=PASS_BOUNDED_SHELL",
    flush=True,
)
for index, row in enumerate(rows[:20]):
    print(
        "Q80FIFTHQ4SEARCH|"
        f"rank={index+1}|child_roots={tuple(row['generic_child_roots'])}|"
        f"cm24_child_roots={tuple(row['cm24_child_roots'])}|"
        f"mw={row['generic_child_mw_rank']}|"
        f"old_degree={row['old_fiber_degree']}|"
        f"old_O={row['old_zero_pairing']}|"
        f"mw_coordinates={tuple(row['cm24_mw_coordinates'])}|"
        f"section_pairs={tuple(map(tuple, row['known_section_pair_representations']))}|"
        f"v={tuple(row['v'])}",
        flush=True,
    )
