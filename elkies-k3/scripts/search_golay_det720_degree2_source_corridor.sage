#!/usr/bin/env sage-python
"""Exhaust all section-nonnegative degree-two fibres on the Golay MW17 target.

The complete degree-two spectrum has 64,355 genus-one translation cosets of
minimum norm eight.  On the rootless fibration these give the natural cheapest
candidate elliptic-neighbour fibres.  This script constructs the exact child
frame for every such coset.  It also proves that the only remaining eligible
cosets are 160 classes of minimum norm twelve, constructs those children, and
tests every primitive child for integral isometry to the marked
``G720-S0128`` 3A5/MW2 source.

The census is complete for degree-two nef candidates on a rootless frame:
section nonnegativity requires coset minimum at least eight, and isotropic
integrality requires that minimum to be zero modulo four.  It does not search
degrees three or more or multi-edge routes.
"""

from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
import multiprocessing
import os
from pathlib import Path

from sage.all import QQ, ZZ, block_diagonal_matrix, gcd, matrix, pari, vector
from sage.modules.free_quadratic_module_integer_symmetric import IntegralLattice


ROOT = Path(__file__).resolve().parents[2]
TARGET = ROOT / "artifacts/generated-results/elkies-k3-golay-octad-rank17-det720.json"
SOURCES = (
    ROOT
    / "artifacts/generated-results/elkies-k3-golay-octad-det720-prescribed-root-sources-v1.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-golay-det720-degree2-direct-3a5-corridor-v1.json"
)
SOURCE_ID = "G720-S0128"
U = matrix(ZZ, ((0, 1), (1, 0)))

engine_path = ROOT / "elkies-k3/scripts/exact_neighbor_engine.sage"
engine = {"__file__": str(engine_path), "__name__": "golay_degree2_engine"}
exec(compile(engine_path.read_text(), str(engine_path), "exec"), engine)
primitive_hyperbolic_split = engine["primitive_hyperbolic_split"]
roots_and_data = engine["roots_and_data"]

NS = None
SOURCE_GRAM = None


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def rows(value) -> list[list[int]]:
    return [[int(entry) for entry in row] for row in value.rows()]


def integral_isometry(left, right):
    raw = pari(left).qfisom(pari(right))
    if raw == 0:
        return None
    candidate = matrix(ZZ, raw)
    for value in (candidate, candidate.transpose()):
        if value * left * value.transpose() == right:
            return value
        if value * right * value.transpose() == left:
            return value.inverse().change_ring(ZZ)
    raise ArithmeticError("PARI returned an unrecognized isometry orientation")


def classify_chunk(task):
    chunk_index, representatives = task
    histogram = Counter()
    divisibility_histogram = Counter()
    candidates = 0
    hits = []
    minimum_histogram = Counter()
    for mask, w_entries, minimum_norm in representatives:
        w = vector(ZZ, w_entries)
        if int(w * (-NS[2:, 2:]) * w) != int(minimum_norm):
            raise ArithmeticError("stored coset representative norm changed")
        minimum_histogram[int(minimum_norm)] += 1
        divisor = vector(ZZ, [int(minimum_norm) // 4, 2] + list(w))
        divisibility = int(gcd(list(NS * divisor)))
        divisibility_histogram[divisibility] += 1
        if divisibility != 1:
            continue
        split = primitive_hyperbolic_split(NS, divisor)
        child = split["child_frame"]
        root_rank, signed_roots, root_determinant = roots_and_data(child)[2]
        signature = f"r{root_rank}-n{signed_roots}-d{root_determinant}"
        histogram[signature] += 1
        if (int(root_rank), int(signed_roots), int(root_determinant)) != (
            15,
            90,
            216,
        ):
            continue
        candidates += 1
        isometry = integral_isometry(child, SOURCE_GRAM)
        if isometry is None:
            continue
        hits.append(
            {
                "residue_mask": int(mask),
                "minimum_norm_representative": list(map(int, w)),
                "fibre_in_U_plus_reduced_target_frame": list(map(int, divisor)),
                "child_frame": rows(child),
                "target_fibre_to_child_transport": rows(split["transport"]),
                "child_to_marked_source_isometry": rows(isometry),
            }
        )
    return {
        "chunk_index": chunk_index,
        "classes": len(representatives),
        "root_signature_histogram": dict(histogram),
        "minimum_norm_histogram": dict(minimum_histogram),
        "fibre_divisibility_histogram": dict(divisibility_histogram),
        "three_a5_signature_candidates": candidates,
        "marked_source_hits": hits,
    }


def build_representatives(target_gram):
    reduction = target_gram.LLL_gram().transpose()
    reduced = reduction * target_gram * reduction.transpose()
    if abs(reduction.det()) != 1:
        raise ArithmeticError("target LLL change is not unimodular")
    pari.allocatemem(4 * 1024**3)
    shell = pari(reduced).qfminim(10)
    columns = matrix(ZZ, shell[2].sage()).columns()
    masks_by_norm = {4: set(), 6: set(), 8: set(), 10: set()}
    representative_by_norm_and_mask = {}
    for column in columns:
        value = vector(ZZ, column)
        norm = int(value * reduced * value)
        if norm not in masks_by_norm:
            continue
        mask = sum((int(entry) % 2) << index for index, entry in enumerate(value))
        masks_by_norm[norm].add(mask)
        representative_by_norm_and_mask.setdefault(
            (norm, mask), list(map(int, value))
        )
    seen = set().union(*masks_by_norm.values())
    missing = sorted(set(range(2**17)) - seen)
    lattice = IntegralLattice(reduced)
    minimum_by_mask = {0: (0, [0] * 17)}
    prior = set()
    for norm in sorted(masks_by_norm):
        for mask in sorted(masks_by_norm[norm] - prior):
            minimum_by_mask[mask] = (
                norm,
                representative_by_norm_and_mask[(norm, mask)],
            )
        prior.update(masks_by_norm[norm])
    for mask in missing:
        if mask == 0:
            continue
        residue = vector(ZZ, [(mask >> index) & 1 for index in range(17)])
        closest = vector(
            ZZ,
            next(lattice.enumerate_close_vectors(vector(QQ, residue) / 2)),
        )
        representative = residue - 2 * closest
        norm = int(representative * reduced * representative)
        if norm <= 10:
            raise ArithmeticError("shell-complement CVP contradicts qfminim(10)")
        minimum_by_mask[mask] = (norm, list(map(int, representative)))
    if len(minimum_by_mask) != 2**17:
        raise ArithmeticError("degree-two coset minimum coverage is incomplete")
    minimum_histogram = Counter(norm for norm, unused in minimum_by_mask.values())
    if set(minimum_histogram) != {0, 4, 6, 8, 10, 12}:
        raise ArithmeticError("unexpected degree-two covering spectrum")
    eligible = [
        (mask, representative, norm)
        for mask, (norm, representative) in sorted(minimum_by_mask.items())
        if norm >= 8 and norm % 4 == 0
    ]
    return (
        reduction,
        reduced,
        int(shell[0]),
        masks_by_norm,
        minimum_histogram,
        len(missing),
        eligible,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--workers", type=int, default=min(16, os.cpu_count() or 1))
    parser.add_argument("--chunk-size", type=int, default=512)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    if arguments.workers <= 0 or arguments.chunk_size <= 0:
        parser.error("workers and chunk size must be positive")

    target_payload = json.loads(TARGET.read_text())
    source_payload = json.loads(SOURCES.read_text())
    target_gram = matrix(ZZ, target_payload["frame"]["gram"])
    source_row = next(
        row for row in source_payload["sources"] if row["source_id"] == SOURCE_ID
    )
    source_gram = matrix(ZZ, source_row["source"]["gram"])
    (
        reduction,
        reduced,
        signed_shell_count,
        masks_by_norm,
        coset_minimum_histogram,
        shell_complement_classes,
        representatives,
    ) = (
        build_representatives(target_gram)
    )
    if len(representatives) != 64_515:
        raise ArithmeticError("degree-two eligible-class count changed")

    global NS, SOURCE_GRAM
    NS = block_diagonal_matrix(U, -reduced)
    SOURCE_GRAM = source_gram
    tasks = [
        (index, representatives[start : start + arguments.chunk_size])
        for index, start in enumerate(
            range(0, len(representatives), arguments.chunk_size)
        )
    ]
    results = []
    with ProcessPoolExecutor(
        max_workers=arguments.workers,
        mp_context=multiprocessing.get_context("fork"),
    ) as executor:
        futures = [executor.submit(classify_chunk, task) for task in tasks]
        for future in as_completed(futures):
            results.append(future.result())
    results.sort(key=lambda row: row["chunk_index"])

    histogram = Counter()
    divisibility_histogram = Counter()
    classified_minimum_histogram = Counter()
    hits = []
    signature_candidates = 0
    for result in results:
        histogram.update(result["root_signature_histogram"])
        divisibility_histogram.update(result["fibre_divisibility_histogram"])
        classified_minimum_histogram.update(result["minimum_norm_histogram"])
        signature_candidates += result["three_a5_signature_candidates"]
        hits.extend(result["marked_source_hits"])
    if sum(divisibility_histogram.values()) != len(representatives):
        raise ArithmeticError("corridor census accounting mismatch")
    primitive_classes = int(divisibility_histogram[1])
    if sum(histogram.values()) != primitive_classes:
        raise ArithmeticError("primitive child accounting mismatch")

    output = {
        "schema": "elkies-k3.golay-det720-degree2-direct-3a5-corridor.v1",
        "status": (
            "PASS_COMPLETE_SECTION_NONNEGATIVE_DEGREE2_DIRECT_SOURCE_HIT"
            if hits
            else "PASS_COMPLETE_SECTION_NONNEGATIVE_DEGREE2_DIRECT_SOURCE_EMPTY"
        ),
        "inputs": {relative(TARGET): digest(TARGET), relative(SOURCES): digest(SOURCES)},
        "source": {
            "source_id": SOURCE_ID,
            "source_gram_sha256": source_row["source"]["gram_sha256"],
            "root_type": source_row["source"]["root_type"],
            "mw_rank": int(source_row["source"]["mw_rank_for_rho_19"]),
        },
        "target": {
            "frame_id": "G720-F001",
            "target_basis_to_reduced_basis": rows(reduction),
            "reduced_target_frame": rows(reduced),
        },
        "degree_two_census": {
            "pari_signed_vectors_through_norm_ten": signed_shell_count,
            "norm_four_residue_masks": len(masks_by_norm[4]),
            "norm_six_residue_masks": len(masks_by_norm[6]),
            "norm_eight_residue_masks": len(masks_by_norm[8]),
            "norm_ten_residue_masks": len(masks_by_norm[10]),
            "shell_complement_cosets_checked_by_exact_closest_vectors": shell_complement_classes,
            "all_coset_minimum_norm_histogram": {
                str(key): value for key, value in sorted(coset_minimum_histogram.items())
            },
            "eligible_genus_one_degree_two_fibre_classes": len(representatives),
            "eligible_minimum_norm_histogram": {
                str(key): value for key, value in sorted(classified_minimum_histogram.items())
            },
            "classes_tested": sum(row["classes"] for row in results),
            "fibre_divisibility_histogram": {
                str(key): value for key, value in sorted(divisibility_histogram.items())
            },
            "primitive_elliptic_fibre_classes": primitive_classes,
            "root_signature_histogram": dict(sorted(histogram.items())),
            "three_a5_signature_candidates": signature_candidates,
            "marked_source_isometry_hits": len(hits),
        },
        "hits": hits,
        "proof_boundary": {
            "proved": (
                "Every target M/2M coset minimum is determined through norm ten by "
                "complete PARI shells and on the 161-class shell complement by exact "
                "closest vectors. Every section-nonnegative integral genus-one class "
                "is represented, its fibre divisibility is computed, every divisibility-one "
                "elliptic fibre is split, and every child with the 3A5 root signature "
                "is tested for an integral isometry "
                "to G720-S0128."
            ),
            "not_proved": (
                "The search excludes higher-degree edges and multi-edge corridors. A "
                "direct lattice hit would still require a "
                "separate equation-level neighbour construction; an empty result is not "
                "a lower bound on the unrestricted corridor length."
            ),
        },
        "reproduce": (
            "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
            "elkies-k3/scripts/search_golay_det720_degree2_source_corridor.sage "
            f"--workers {arguments.workers} --chunk-size {arguments.chunk_size}"
        ),
    }
    serialized = json.dumps(output, indent=2, sort_keys=True) + "\n"
    output_path = arguments.output.resolve()
    if arguments.check:
        if output_path.read_text() != serialized:
            raise SystemExit("Golay determinant-720 degree-two corridor artifact is stale")
    else:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(serialized)
    print(
        "GOLAY720D2CORRIDOR|classes={}|3A5_signature={}|hits={}|status=PASS".format(
            len(representatives), signature_candidates, len(hits)
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
