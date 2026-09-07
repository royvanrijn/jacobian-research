#!/usr/bin/env sage
"""Run the immediate source-first gate for the determinant-720 Golay frame.

This deterministic bounded Kneser-neighbour beam searches the positive-
definite genus of the exact Golay-octad rank-17 frame for a primitive-root
companion of root rank at least 15 (MW rank at most two for rho=19).

The emitted frame, ADE data, root primitivity, MW height quotient, genus
identity, and p-neighbour replay are exact.  The p-neighbour path is discovery
provenance, not an elliptic-neighbour corridor.  Because the discriminant
group is noncyclic, a separate full discriminant-gluing certificate is still
required before the source is promoted to a primitive embedding of the same
auxiliary in a named Niemeier lattice.

status: ACTIVE_SEARCH
claim: exact bounded source-first genus search; no classification or Niemeier
  embedding claim.
inputs: artifacts/generated-results/elkies-k3-golay-octad-rank17-det720.json
outputs: artifacts/generated-results/elkies-k3-golay-octad-det720-source-hunt.json
"""

import argparse
import hashlib
import json
import random as pyrandom
from collections import Counter
from pathlib import Path

from sage.all import Genus, ZZ, matrix, pari, set_random_seed, vector


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DEFAULT_INPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-golay-octad-rank17-det720.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-golay-octad-det720-source-hunt.json"
)

# Reuse the exact genus/root/frame kernel of the established foundry hunter,
# stopping before its CLI.  This avoids a second implementation of ADE and
# primitive root-plus-MW extraction.
_shared_path = HERE / "hunt_lattice_foundry_rootful_source.sage"
_shared_source = _shared_path.read_text().split(
    "parser = argparse.ArgumentParser", 1
)[0]
_shared = {"__file__": str(_shared_path)}
exec(compile(_shared_source, str(_shared_path), "exec"), _shared)

form_from_gram = _shared["form_from_gram"]
reduced_gram = _shared["reduced_gram"]
reduced_key = _shared["reduced_key"]
root_rank_and_count = _shared["root_rank_and_count"]
minimize_child_frame = _shared["minimize_child_frame"]
ade_name = _shared["ade_name"]
gram_digest = _shared["gram_digest"]
rows = _shared["rows"]
rational_rows = _shared["rational_rows"]
replay_path = _shared["replay_path"]


def search(arguments, target_payload):
    target = matrix(ZZ, target_payload["frame"]["gram"])
    determinant = int(target.det())
    assert target.nrows() == 17 and determinant == 720
    assert Genus(target).discriminant_form().normal_form().invariants() == (2, 6, 60)

    pyrandom.seed(arguments.seed)
    set_random_seed(arguments.seed)
    requested_primes = [ZZ(value) for value in arguments.primes.split(",")]
    primes = [prime for prime in requested_primes if determinant % prime]
    if not primes:
        raise ValueError("every requested p-neighbour prime divides 720")

    start_rank, start_count = root_rank_and_count(target)
    frontier = [
        {
            "root_rank": start_rank,
            "root_count": start_count,
            "form": form_from_gram(target),
            "path": [],
            "path_digests": [],
        }
    ]
    retained = list(frontier)
    seen = {reduced_key(target)}
    best = frontier[0]
    accounting = []

    print(
        f"GOLAYSOURCE|stage=start|det={determinant}|root_rank={start_rank}"
        f"|target_root_rank={arguments.target_root_rank}",
        flush=True,
    )
    for generation in range(1, arguments.generations + 1):
        candidates = []
        failures = 0
        for parent in frontier:
            for unused_sample in range(arguments.samples_per_parent):
                prime = pyrandom.choice(primes)
                try:
                    witness = parent["form"].find_primitive_p_divisible_vector__random(
                        prime
                    )
                    child = parent["form"].find_p_neighbor_from_vec(prime, witness)
                    gram = child.Hessian_matrix()
                    if gram.det() != determinant or any(
                        gram[index, index] % 2 for index in range(17)
                    ):
                        failures += 1
                        continue
                    key = reduced_key(gram)
                    if key in seen:
                        continue
                    seen.add(key)
                    root_rank, root_count = root_rank_and_count(gram)
                except Exception:
                    failures += 1
                    continue
                candidates.append(
                    {
                        "root_rank": root_rank,
                        "root_count": root_count,
                        "form": child,
                        "path": parent["path"]
                        + [
                            {
                                "prime": int(prime),
                                "primitive_p_divisible_vector": list(
                                    map(int, witness)
                                ),
                            }
                        ],
                        "path_digests": parent["path_digests"]
                        + [gram_digest(gram)],
                    }
                )

        candidates.sort(
            key=lambda item: (-item["root_rank"], -item["root_count"])
        )
        next_frontier = []
        per_profile = Counter()
        for candidate in candidates:
            profile = (candidate["root_rank"], candidate["root_count"])
            if per_profile[profile] >= arguments.per_profile:
                continue
            per_profile[profile] += 1
            next_frontier.append(candidate)
            if len(next_frontier) == arguments.beam:
                break
        if not next_frontier:
            raise RuntimeError("Golay source beam died before producing a candidate")
        frontier = next_frontier
        retained.extend(frontier)
        if (
            frontier[0]["root_rank"],
            frontier[0]["root_count"],
        ) > (best["root_rank"], best["root_count"]):
            best = frontier[0]
        accounting.append(
            {
                "generation": generation,
                "new_candidates": len(candidates),
                "failures": failures,
                "visited_reduced_keys": len(seen),
                "best_root_rank": best["root_rank"],
                "best_signed_root_count": best["root_count"],
                "frontier_profiles": [
                    [item["root_rank"], item["root_count"]]
                    for item in frontier
                ],
            }
        )
        print(
            f"GOLAYSOURCE|stage=generation|generation={generation}"
            f"|new={len(candidates)}|seen={len(seen)}"
            f"|best_root_rank={best['root_rank']}|best_roots={best['root_count']}",
            flush=True,
        )
        if best["root_rank"] >= arguments.target_root_rank:
            break

    certifiable = None
    for candidate in sorted(
        retained,
        key=lambda item: (
            -item["root_rank"],
            -item["root_count"],
            len(item["path"]),
        ),
    ):
        minimized = minimize_child_frame(candidate["form"].Hessian_matrix())
        if minimized["root_lattice_primitive"] and minimized["mw_height"] is not None:
            certifiable = candidate, minimized
            break
    if certifiable is None:
        raise RuntimeError("bounded Golay source hunt found no primitive-root frame")

    winner, minimized = certifiable
    if (
        winner["root_rank"] < arguments.target_root_rank
        and not arguments.allow_below_target
    ):
        raise RuntimeError(
            f"best primitive-root source has rank {winner['root_rank']}, below "
            f"the requested {arguments.target_root_rank}"
        )

    source_raw = winner["form"].Hessian_matrix()
    source, source_change = reduced_gram(source_raw)
    source_minimized = minimize_child_frame(source)
    source_ade, components = ade_name(source)
    source_rank, source_count = root_rank_and_count(source)
    assert Genus(source) == Genus(target)
    assert source_rank == winner["root_rank"]
    assert source_minimized["root_lattice_primitive"]
    assert source_minimized["mw_height"] is not None

    replayed, replayed_digests = replay_path(
        target,
        winner["path"],
        winner["path_digests"],
        determinant,
    )
    assert replayed.Hessian_matrix() == source_raw

    return {
        "schema": "elkies-k3.golay-octad-source-hunt.v1",
        "status": (
            f"PASS_EXACT_BOUNDED_DET720_ROOT_RANK_{source_rank}_GENUS_SOURCE"
            "_NIEMEIER_GLUE_OPEN"
        ),
        "proof_scope": {
            "proved": (
                "The emitted rank-17 source frame is in the exact genus of the "
                "Golay target, has the displayed complete ADE root system, a "
                "primitive root lattice, and the displayed exact MW height "
                "quotient. Every Kneser p-neighbour edge replays exactly."
            ),
            "not_proved": (
                "The bounded beam is not a genus classification. The source "
                "has not yet been glued through the noncyclic discriminant "
                "module to the pinned auxiliary in a named Niemeier lattice. "
                "The p-neighbour path is not an elliptic-neighbour corridor, "
                "and no equation or arithmetic descent is claimed."
            ),
        },
        "input": {
            "artifact": str(arguments.input.resolve().relative_to(ROOT)),
            "sha256": hashlib.sha256(arguments.input.read_bytes()).hexdigest(),
        },
        "search": {
            "seed": arguments.seed,
            "generations_bound": arguments.generations,
            "generations_used": len(accounting),
            "beam": arguments.beam,
            "samples_per_parent": arguments.samples_per_parent,
            "per_profile": arguments.per_profile,
            "primes": list(map(int, primes)),
            "target_root_rank": arguments.target_root_rank,
            "target_reached_by_emitted_primitive_source": (
                source_rank >= arguments.target_root_rank
            ),
            "allow_below_target": arguments.allow_below_target,
            "visited_reduced_keys": len(seen),
            "generation_accounting": accounting,
        },
        "target": {
            "gram": rows(target),
            "gram_sha256": gram_digest(target),
            "rank": 17,
            "root_rank": start_rank,
            "mw_rank_for_rho_19": 17 - start_rank,
            "determinant": determinant,
        },
        "source": {
            "gram": rows(source),
            "gram_sha256": gram_digest(source),
            "raw_search_gram_sha256": gram_digest(source_raw),
            "reduced_basis_rows_in_raw_basis": rows(source_change),
            "root_type": source_ade,
            "root_components": components,
            "root_rank": source_rank,
            "signed_root_count": source_count,
            "mw_rank_for_rho_19": 17 - source_rank,
            "root_lattice_primitive": True,
            "root_smith_invariants": list(
                map(int, source_minimized["root_smith_invariants"])
            ),
            "root_adapted_gram": rows(source_minimized["frame"]),
            "root_adapted_basis_rows_in_source_basis": rows(
                source_minimized["basis"]
            ),
            "mw_height_gram": rational_rows(source_minimized["mw_height"]),
            "mw_regulator": str(abs(source_minimized["mw_height"].det())),
            "genus_equals_target": True,
            "discriminant_form_equals_target": True,
            "determinant": determinant,
        },
        "kneser_p_neighbor_provenance": {
            "edge_count": len(winner["path"]),
            "edges": winner["path"],
            "child_raw_gram_sha256_by_edge": winner["path_digests"],
            "exact_replay_passed": replayed_digests == winner["path_digests"],
            "warning": (
                "These positive-definite p-neighbours discover a genus mate; "
                "they are not elliptic-neighbour edges on a K3 surface."
            ),
        },
    }


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--generations", type=int, default=18)
parser.add_argument("--beam", type=int, default=20)
parser.add_argument("--samples-per-parent", type=int, default=100)
parser.add_argument("--per-profile", type=int, default=3)
parser.add_argument("--primes", default="7,11,13,17,19,23")
parser.add_argument("--seed", type=int, default=20260901)
parser.add_argument("--target-root-rank", type=int, default=15)
parser.add_argument("--allow-below-target", action="store_true")
parser.add_argument("--check", action="store_true")
arguments = parser.parse_args()

target_payload = json.loads(arguments.input.read_text())
result = search(arguments, target_payload)
rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
if arguments.check:
    assert arguments.output.read_text() == rendered
else:
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(rendered)

print(
    f"GOLAYSOURCE|stage=done|source={result['source']['root_type']}"
    f"|root_rank={result['source']['root_rank']}"
    f"|mw_rank={result['source']['mw_rank_for_rho_19']}"
    f"|edges={result['kneser_p_neighbor_provenance']['edge_count']}"
    "|status=PASS_EXACT_GENUS_SOURCE_NIEMEIER_GLUE_OPEN"
)
