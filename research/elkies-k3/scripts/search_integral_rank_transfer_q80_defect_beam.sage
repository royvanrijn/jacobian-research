#!/usr/bin/env sage-python
"""Transport and annihilate Q80 masked witnesses through directed neighbors."""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import random as python_random
import runpy

from sage.all import ZZ, pari, set_random_seed, vector


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
DIRECTED_SCRIPT = ROOT / "elkies-k3/scripts/search_integral_rank_transfer_q80_defect_neighbors.sage"
OUTPUT = GENERATED / "elkies-k3-integral-rank-transfer-q80-defect-beam-v1.json"

RANDOM_SEED = 161803
PRIMES = (7, 11, 13, 17, 19, 29, 31, 37, 41, 43)
GENERATIONS = 5
INITIAL_SAMPLES = 10000
SAMPLES_PER_PARENT = 500
BEAM_WIDTH = 20


def relative(path):
    return str(Path(path).resolve().relative_to(ROOT))


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def line_key(witness, prime):
    values = [int(value % prime) for value in witness]
    pivot = next(value for value in values if value)
    inverse = pow(pivot, -1, prime)
    return tuple((value * inverse) % prime for value in values)


def distinct_isometry_elite(candidates, width, base):
    """Keep the lowest-defect candidates in distinct integral isometry classes."""

    selected = []
    reduced = []
    comparisons = 0
    for candidate in sorted(candidates, key=lambda row: row[:2]):
        reduced_gram = base["lll_reduce"](candidate[2].Hessian_matrix())
        duplicate = False
        for known in reduced:
            comparisons += 1
            if pari(reduced_gram).qfisom(pari(known)):
                duplicate = True
                break
        if duplicate:
            continue
        selected.append(candidate)
        reduced.append(reduced_gram)
        if len(selected) == width:
            break
    return selected, comparisons


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--generations", type=int, default=GENERATIONS)
    parser.add_argument("--initial-samples", type=int, default=INITIAL_SAMPLES)
    parser.add_argument("--samples-per-parent", type=int, default=SAMPLES_PER_PARENT)
    parser.add_argument("--beam-width", type=int, default=BEAM_WIDTH)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    arguments = parser.parse_args()
    if min(
        arguments.generations,
        arguments.initial_samples,
        arguments.samples_per_parent,
        arguments.beam_width,
    ) <= 0:
        parser.error("all search bounds must be positive")

    directed = runpy.run_path(str(DIRECTED_SCRIPT))
    base = runpy.run_path(str(directed["BASE_SCRIPT"]))
    search = runpy.run_path(str(directed["SEARCH_SCRIPT"]))
    control = runpy.run_path(str(directed["CONTROL_SCRIPT"]))
    core = runpy.run_path(str(directed["CORE_SCRIPT"]))
    reverse = runpy.run_path(str(directed["REVERSE_SCRIPT"]))
    prepared, bridge, initial_gram, prefix_states, _, _ = directed["initial_q80"](
        base, search, control, core, reverse
    )
    initial_form = base["quadratic_form"](initial_gram)
    frontier = [(initial_form, [])]
    set_random_seed(RANDOM_SEED)
    python_random.seed(RANDOM_SEED)
    seen_neighbors = set()
    generation_records = []
    hit = None

    for generation in range(1, arguments.generations + 1):
        parent_count = len(frontier)
        candidates = []
        root_histogram = Counter()
        defect_histogram = Counter()
        directed_lines = 0
        rootless_neighbors = 0
        old_witness_count = 0
        samples_per_parent = (
            arguments.initial_samples
            if generation == 1
            else arguments.samples_per_parent
        )
        for parent_index, (parent, path) in enumerate(frontier, start=1):
            parent_gram = parent.Hessian_matrix()
            parent_mask, _, witnesses = directed["masked_witness_data"](
                parent_gram,
                bridge,
                prepared["order"],
                base,
                reverse,
            )
            assert parent_mask["occupied_forbidden_cells"] > 0
            old_witness_count += len(witnesses)
            inverse = parent_gram.inverse()
            seen_lines = set()
            for _ in range(samples_per_parent):
                prime = python_random.choice(PRIMES)
                raw_witness = parent.find_primitive_p_divisible_vector__random(prime)
                key = (prime, line_key(raw_witness, prime))
                if key in seen_lines:
                    continue
                seen_lines.add(key)
                residue_pairings = [
                    int(witness * raw_witness) % prime for witness in witnesses
                ]
                if any(value == 0 for value in residue_pairings):
                    continue
                directed_lines += 1
                transform = parent.find_p_neighbor_from_vec(
                    prime, raw_witness, return_matrix=True
                )
                neighbor = parent.find_p_neighbor_from_vec(prime, raw_witness)
                neighbor_gram = neighbor.Hessian_matrix()
                neighbor_key = hashlib.sha256(
                    repr(tuple(neighbor_gram.list())).encode()
                ).hexdigest()
                if neighbor_key in seen_neighbors:
                    continue
                seen_neighbors.add(neighbor_key)

                ambient_basis = transform.transpose()
                for witness in witnesses:
                    dual_vector = witness * inverse
                    pairings = ambient_basis * parent_gram * dual_vector.column()
                    assert not all(value in ZZ for value in pairings)

                roots = int(pari(neighbor_gram).qfminim(2)[0])
                root_histogram[roots] += 1
                if roots:
                    continue
                rootless_neighbors += 1
                masks, _, _ = base["mask_profile"](
                    neighbor_gram, [bridge], reverse
                )
                defect = masks[0]["occupied_forbidden_cells"]
                defect_histogram[defect] += 1
                new_path = path + [
                    {"prime": prime, "witness": list(map(int, raw_witness))}
                ]
                if defect == 0:
                    hit = {
                        "form": neighbor,
                        "path": new_path,
                        "mask": masks[0],
                    }
                    break
                candidates.append(
                    (defect, python_random.random(), neighbor, new_path)
                )
            if hit is not None:
                break
            print(
                f"generation={generation} parent={parent_index}/{len(frontier)} "
                f"directed={directed_lines} rootless={rootless_neighbors} "
                f"defects={dict(defect_histogram)}",
                flush=True,
            )
        if hit is not None:
            generation_records.append(
                {
                    "generation": generation,
                    "parent_count": parent_count,
                    "samples_per_parent": samples_per_parent,
                    "directed_lines": directed_lines,
                    "rootless_neighbors": rootless_neighbors,
                    "root_count_histogram": {
                        str(key): root_histogram[key] for key in sorted(root_histogram)
                    },
                    "mask_defect_histogram": {
                        str(key): defect_histogram[key]
                        for key in sorted(defect_histogram)
                    },
                    "hit": True,
                }
            )
            break
        selected, comparisons = distinct_isometry_elite(
            candidates, arguments.beam_width, base
        )
        frontier = [(row[2], row[3]) for row in selected]
        generation_records.append(
            {
                "generation": generation,
                "parent_count": parent_count,
                "samples_per_parent": samples_per_parent,
                "mean_old_physical_witnesses_per_parent": str(
                    old_witness_count / max(1, parent_count)
                ),
                "directed_lines": directed_lines,
                "rootless_neighbors": rootless_neighbors,
                "root_count_histogram": {
                    str(key): root_histogram[key] for key in sorted(root_histogram)
                },
                "mask_defect_histogram": {
                    str(key): defect_histogram[key]
                    for key in sorted(defect_histogram)
                },
                "candidate_count": len(candidates),
                "isometry_comparisons": comparisons,
                "distinct_isometry_frontier": len(frontier),
                "best_defect": selected[0][0] if selected else None,
                "hit": False,
            }
        )
        print(
            f"generation={generation} selected={len(frontier)} "
            f"best={selected[0][0] if selected else None}",
            flush=True,
        )
        if not frontier:
            break

    if hit is not None:
        reduced = base["lll_reduce"](hit["form"].Hessian_matrix())
        masks, _, _ = base["mask_profile"](reduced, [bridge], reverse)
        assert masks[0]["zero_mask_accepts"]
        completion = control["completion"](
            reduced, prepared, masks[0], base, core
        )
        hit_record = {
            "directed_path_after_near_miss": hit["path"],
            "full_path_from_canonical_seed": prefix_states + [
                {
                    "step": len(prefix_states) + index,
                    **edge,
                }
                for index, edge in enumerate(hit["path"], start=1)
            ],
            "reduced_core_gram": [list(map(int, row)) for row in reduced.rows()],
            "completion": completion,
        }
    else:
        hit_record = None

    payload = {
        "schema": "elkies-k3.integral-rank-transfer-q80-defect-beam.v1",
        "status": (
            "PASS_MULTISTEP_DEFECT_DIRECTED_Q80_COMPLETION"
            if hit is not None
            else "PASS_BOUNDED_MULTISTEP_DEFECT_DIRECTED_Q80_MISS"
        ),
        "inputs": {
            relative(DIRECTED_SCRIPT): digest(DIRECTED_SCRIPT),
            relative(directed["BRIDGES"]): digest(directed["BRIDGES"]),
            relative(directed["THETA"]): digest(directed["THETA"]),
            relative(directed["BASE_SCRIPT"]): digest(directed["BASE_SCRIPT"]),
            relative(directed["SEARCH_SCRIPT"]): digest(directed["SEARCH_SCRIPT"]),
            relative(directed["CONTROL_SCRIPT"]): digest(directed["CONTROL_SCRIPT"]),
            relative(directed["CORE_SCRIPT"]): digest(directed["CORE_SCRIPT"]),
            relative(directed["REVERSE_SCRIPT"]): digest(directed["REVERSE_SCRIPT"]),
        },
        "search_rule": {
            "random_seed": RANDOM_SEED,
            "primes": list(PRIMES),
            "generations": arguments.generations,
            "initial_samples": arguments.initial_samples,
            "samples_per_parent": arguments.samples_per_parent,
            "beam_width": arguments.beam_width,
            "mandatory_filter": (
                "the isotropic line is nonorthogonal to every current physical "
                "masked witness"
            ),
            "frontier_rule": (
                "lowest defect, at most one representative per integral "
                "isometry class"
            ),
        },
        "generations": generation_records,
        "unique_constructed_neighbors": len(seen_neighbors),
        "zero_mask_hit": hit is not None,
        "hit": hit_record,
        "proof_boundary": (
            "Every constructed edge kills every physical witness responsible "
            "for its parent defect. Any child defect is therefore newly "
            "generated. The beam is finite and sampled; a miss is not a "
            "genus-wide obstruction or a proof that defect descent is impossible."
        ),
        "reproduce": (
            "sage -python elkies-k3/scripts/"
            "search_integral_rank_transfer_q80_defect_beam.sage --check"
        ),
    }
    output = arguments.output if arguments.output.is_absolute() else ROOT / arguments.output
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.check:
        if not output.exists():
            raise SystemExit(f"missing artifact: {output}")
        if output.read_text() != encoded:
            raise SystemExit(f"stale artifact: {output}")
        print(payload["status"])
        return
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(encoded)
    print(relative(output))


if __name__ == "__main__":
    main()
