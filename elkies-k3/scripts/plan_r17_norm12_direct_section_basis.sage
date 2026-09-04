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
import hashlib
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


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def reduce_basis(fixed_rows, candidates):
    rows = list(fixed_rows)
    fixed_count = len(rows)
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
            for row_index in range(fixed_count, 17):
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
            break
        _, candidate_index, row_index = best
        chosen[row_index - fixed_count] = candidate_index
        rows[row_index] = candidates[candidate_index][1]
        basis = matrix(ZZ, rows)
    return [candidates[index][0] for index in chosen], basis


def improve_basis(fixed_rows, candidates):
    result = reduce_basis(fixed_rows, candidates)
    if result is None or abs(result[1].det()) != 1:
        return None
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-label", required=True)
    parser.add_argument("--height-bound", type=int, default=8)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.height_bound < 4 or args.height_bound > 12 or args.height_bound % 2:
        parser.error("--height-bound must be an even integer from 4 through 12")

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

    def emit(payload):
        payload.update(
            {
                "schema": "elkies-k3.r17-norm12-direct-section-basis-plan.v1",
                "inputs": {
                    relative(path): digest(path)
                    for path in (PINNED, SPLITTING, BISECTIONS)
                },
                "proof_boundary": (
                    (
                        "For a degree-one old section of height h, <w,v>=h and "
                        "(v-w)^2=12-h. Rootlessness excludes h=10 and positive "
                        "definiteness forces v=w at h=12, so the height-eight "
                        "Fincke-Pohst shell plus v=w is complete through the "
                        "Cauchy bound. "
                        if args.height_bound == 12
                        else f"The old-section enumeration is bounded at height {args.height_bound}. "
                    )
                    + "The rational-bisection scan covers every record in the "
                    "committed 39120-class equation atlas."
                ),
            }
        )
        serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
        if args.output is not None:
            if args.check:
                if not args.output.exists() or args.output.read_text() != serialized:
                    raise ArithmeticError("stored basis-plan artifact differs from replay")
            else:
                args.output.parent.mkdir(parents=True, exist_ok=True)
                args.output.write_text(serialized)
        elif args.check:
            parser.error("--check requires --output")
        print(serialized, end="")

    # For a degree-one old section of height h, <w,v>=h.  With w^2=12,
    # (v-w)^2=12-h.  The published R17 lattice is rootless, so h=10 is
    # impossible; h=12 forces v=w.  Thus the height-at-most-eight enumeration
    # plus that single equality case is the complete search through height 12.
    short = matrix(ZZ, pari(pinned).qfminim(min(args.height_bound, 8))[2])
    old_vectors = set()
    for column in short.columns():
        entries = tuple(map(int, column))
        old_vectors.add(entries)
        old_vectors.add(tuple(-value for value in entries))
    if args.height_bound >= 12:
        old_vectors.add(tuple(map(int, source["pinned_rank17_w"])))
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

    old_only = improve_basis([], candidates)
    if old_only is not None:
        selected, basis = old_only
        payload = {
            "status": "PASS_EXACT_UNIMODULAR_OLD_SECTION_MARKING",
            "source_label": args.source_label,
            "height_bound": args.height_bound,
            "old_degree_one_candidate_count": len(candidates),
            "rational_bisection_glue_candidate_count": len(glue_candidates),
            "glue_label": None,
            "selected_old_vectors": [list(map(int, row)) for row in selected],
            "child_coordinate_determinant": int(basis.det()),
        }
        emit(payload)
        return

    for glue_label, glue_row in glue_candidates:
        result = improve_basis([glue_row], candidates)
        if result is None:
            continue
        selected, basis = result
        payload = {
            "status": "PASS_EXACT_UNIMODULAR_OLD_SECTION_PLUS_BISECTION_MARKING",
            "source_label": args.source_label,
            "height_bound": args.height_bound,
            "old_degree_one_candidate_count": len(candidates),
            "rational_bisection_glue_candidate_count": len(glue_candidates),
            "glue_label": glue_label,
            "selected_old_vectors": [list(map(int, row)) for row in selected],
            "child_coordinate_determinant": int(basis.det()),
        }
        emit(payload)
        return
    candidate_matrix = matrix(ZZ, [list(row) for _, row in candidates])
    candidate_module = candidate_matrix.row_module()
    reduced = reduce_basis([], candidates)
    if reduced is None:
        raise ArithmeticError("old degree-one candidates do not span rank 17")
    independent_old_vectors, independent_matrix = reduced
    diagnostic = {
        "source_label": args.source_label,
        "height_bound": args.height_bound,
        "old_degree_one_candidate_count": len(candidates),
        "old_degree_one_span_rank": int(candidate_matrix.rank()),
        "old_degree_one_lattice_saturation_index": int(
            candidate_module.index_in(candidate_module.saturation())
        ),
        "selected_independent_old_vectors": [
            list(map(int, row)) for row in independent_old_vectors
        ],
        "selected_child_coordinate_determinant": int(independent_matrix.det()),
        "rational_bisection_glue_candidate_count": len(glue_candidates),
        "status": "NO_UNIMODULAR_MARKING_IN_SEARCHED_CURVES",
    }
    emit(diagnostic)


if __name__ == "__main__":
    main()
