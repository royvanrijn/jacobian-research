#!/usr/bin/env sage-python
"""Find a unimodular old-curve section basis for a norm-12 two-neighbor.

This is a lattice-only planning utility.  It enumerates the published R17
sections of height at most eight that become degree-one curves for the chosen
neighbor and tests every committed rational bisection as a possible index-two
glue section.  A successful output consists of sixteen old sections and one
glue curve whose child-frame coordinate matrix has determinant one.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from sage.all import QQ, ZZ, block_diagonal_matrix, matrix, pari, vector


ROOT = Path(__file__).resolve().parents[2]
PINNED = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
SPLITTING = ROOT / "artifacts/generated-results/elkies-k3-r17-genus-one-bisection-splitting-search-v1.json"
BISECTIONS = ROOT / "artifacts/generated-results/elkies-2026-equation-bisections-full.json"


def load_matrix(path: Path):
    return matrix(
        ZZ,
        [
            [ZZ(value) for value in line.split()]
            for line in path.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ],
    )


def find_record(payload, label):
    return next(record for record in payload["construction"]["records"] if record["label"] == label)


def improve_basis(glue_row, candidates):
    rows = [glue_row]
    chosen = []
    for index, (_, child_row) in enumerate(candidates):
        trial = matrix(QQ, rows + [child_row])
        if trial.rank() > len(rows):
            rows.append(child_row)
            chosen.append(index)
        if len(rows) == 17:
            break
    if len(rows) != 17:
        return None
    basis = matrix(ZZ, rows)
    while abs(basis.det()) != 1:
        determinant = abs(basis.det())
        inverse = basis.inverse()
        best = None
        for candidate_index, (_, child_row) in enumerate(candidates):
            if candidate_index in chosen:
                continue
            coordinates = vector(QQ, child_row) * inverse
            for row_index in range(1, 17):
                new_determinant = abs(coordinates[row_index] * determinant)
                if (
                    not new_determinant
                    or new_determinant.denominator() != 1
                    or new_determinant >= determinant
                ):
                    continue
                key = (ZZ(new_determinant), candidate_index, row_index)
                if best is None or key < best[0]:
                    best = (key, candidate_index, row_index)
        if best is None:
            return None
        _, candidate_index, row_index = best
        chosen[row_index - 1] = candidate_index
        rows[row_index] = candidates[candidate_index][1]
        basis = matrix(ZZ, rows)
    return [candidates[index][0] for index in chosen], basis


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-label", required=True)
    args = parser.parse_args()

    pinned = load_matrix(PINNED)
    splitting = json.loads(SPLITTING.read_text())
    source = find_record(splitting, args.source_label)
    ns = block_diagonal_matrix(matrix(ZZ, [[0, 1], [1, 0]]), -pinned)
    fibre = vector(ZZ, [3, 2] + list(source["pinned_rank17_w"]))
    old_zero = vector(ZZ, [-1, 1] + [0] * 17)
    mate = fibre + old_zero
    complement = matrix(ZZ, [list(fibre * ns), list(mate * ns)]).right_kernel_matrix()
    transport = matrix(ZZ, [list(fibre), list(mate)] + [list(row) for row in complement.rows()])
    if abs(transport.det()) != 1:
        raise ArithmeticError("neighbor transport is not unimodular")
    transport_inverse = transport.inverse()

    short = matrix(ZZ, pari(pinned).qfminim(8)[2])
    old_vectors = set()
    for column in short.columns():
        entries = tuple(map(int, column))
        old_vectors.add(entries)
        old_vectors.add(tuple(-value for value in entries))
    candidates = []
    for entries in old_vectors:
        old_mw = vector(ZZ, entries)
        height = ZZ(old_mw * pinned * old_mw)
        old_class = vector(ZZ, [(height - 2) // 2, 1] + list(old_mw))
        child_class = old_class * transport_inverse
        if child_class[1] != 1 or any(value not in ZZ for value in child_class):
            continue
        child_row = vector(ZZ, child_class[2:])
        candidates.append((old_mw, child_row))
    candidates.sort(
        key=lambda item: (
            ZZ(item[0] * pinned * item[0]),
            sum(abs(value) for value in item[0]),
            max(abs(value) for value in item[0]),
            tuple(item[0]),
        )
    )

    full_bisections = json.loads(BISECTIONS.read_text())
    glue_candidates = []
    for record in full_bisections["bisections"]:
        glue_class = vector(ZZ, [2, 2] + list(record["pinned_rank17_w"]))
        child_class = glue_class * transport_inverse
        if child_class[1] == 1 and all(value in ZZ for value in child_class):
            glue_candidates.append((record["label"], vector(ZZ, child_class[2:])))
    glue_candidates.sort(key=lambda item: item[0])

    for glue_label, glue_row in glue_candidates:
        result = improve_basis(glue_row, candidates)
        if result is None:
            continue
        selected, basis = result
        payload = {
            "source_label": args.source_label,
            "old_degree_one_candidate_count": len(candidates),
            "rational_bisection_glue_candidate_count": len(glue_candidates),
            "glue_label": glue_label,
            "selected_old_vectors": [list(map(int, row)) for row in selected],
            "child_coordinate_determinant": int(basis.det()),
        }
        print(json.dumps(payload, indent=2, sort_keys=True))
        return
    raise ArithmeticError("no unimodular old-section plus rational-bisection basis found")


if __name__ == "__main__":
    main()
