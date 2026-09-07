#!/usr/bin/env sage-python
"""Classify all 43 exact norm-twelve R17 genus-one bisection frames.

For every stored trace vector ``w`` the divisor ``D=(3,2,w)`` is a primitive
isotropic class in ``NS=U+R17(-1)``.  The exact bisection construction supplies
an irreducible smooth genus-one representative, so ``D`` is nef and ``|D|`` is
an elliptic pencil.  Since the old zero has intersection one with ``D``, it is
also a zero for the new pencil.  This checker splits off

    <D, O_old + D> = U,

classifies the orthogonal frame, and compares it with the two mass-complete
rootless determinant-948 J2 controls.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

from sage.all import ZZ, block_diagonal_matrix, matrix, pari, vector
from sage.env import SAGE_VERSION


ROOT = Path(__file__).resolve().parents[2]
PINNED = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
SPLITTING = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-genus-one-bisection-splitting-search-v1.json"
)
ALTERNATE = (
    ROOT / "artifacts/generated-results/q80-alternate-fifth-q6-rootless-transport.json"
)
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-isotropic-frame-classification-v1.json"
)


def load_matrix(path: Path):
    return matrix(
        ZZ,
        [
            [ZZ(value) for value in line.split()]
            for line in path.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ],
    )


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def rows(value):
    return [list(map(int, row)) for row in value.rows()]


def matrix_digest(value) -> str:
    payload = json.dumps(rows(value), separators=(",", ":"))
    return hashlib.sha256(payload.encode()).hexdigest()


def qfisometric(left, right) -> bool:
    return pari(left).qfisom(pari(right)) != 0


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    published = load_matrix(PINNED)
    alternate_payload = json.loads(ALTERNATE.read_text())
    alternate = matrix(ZZ, alternate_payload["rootless_frame"])
    splitting = json.loads(SPLITTING.read_text())
    records = [
        record
        for record in splitting["construction"]["records"]
        if int(record["trace_norm"]) == 12
    ]
    if len(records) != 43:
        raise ArithmeticError(f"expected 43 norm-twelve records, found {len(records)}")
    if not splitting["construction"][
        "all_branch_quartics_irreducible_squarefree_and_smoothly_branched"
    ]:
        raise ArithmeticError("the exact genus-one source certificate lost its smoothness gate")

    hyperbolic = matrix(ZZ, [[0, 1], [1, 0]])
    ns = block_diagonal_matrix(hyperbolic, -published)
    old_fibre = vector(ZZ, [1, 0] + [0] * 17)
    old_zero = vector(ZZ, [-1, 1] + [0] * 17)
    old_mate = old_fibre + old_zero

    classifications = []
    for record in records:
        if not all(
            record[key]
            for key in (
                "branch_polynomial_irreducible_over_Q",
                "branch_polynomial_squarefree",
                "branch_polynomial_coprime_to_surface_discriminant",
                "branch_polynomial_coprime_to_trace_denominator",
            )
        ):
            raise ArithmeticError(f"{record['label']} lost an exact branch gate")
        if record["member_selection"] != "unique regular M0 member":
            raise ArithmeticError(f"{record['label']} is not a regular norm-twelve member")

        w = vector(ZZ, record["pinned_rank17_w"])
        if w * published * w != 12:
            raise ArithmeticError(f"{record['label']} no longer has norm twelve")
        fibre = vector(ZZ, [3, 2] + list(w))
        new_mate = fibre + old_zero
        cross = matrix(
            ZZ,
            [
                [old_fibre * ns * fibre, old_fibre * ns * new_mate],
                [old_mate * ns * fibre, old_mate * ns * new_mate],
            ],
        )
        if fibre * ns * fibre != 0:
            raise ArithmeticError(f"{record['label']} is not isotropic")
        if fibre * ns * old_fibre != 2 or fibre * ns * old_zero != 1:
            raise ArithmeticError(f"{record['label']} lost the (degree, zero-cost)=(2,1) profile")
        if new_mate * ns * new_mate != 0 or fibre * ns * new_mate != 1:
            raise ArithmeticError(f"{record['label']} does not split the shared-zero U")
        if cross != matrix(ZZ, [[2, 3], [3, 2]]):
            raise ArithmeticError(f"{record['label']} has an unexpected relative-U matrix")

        complement = matrix(
            ZZ, [list(fibre * ns), list(new_mate * ns)]
        ).right_kernel_matrix()
        transport = matrix(
            ZZ,
            [list(fibre), list(new_mate)] + [list(row) for row in complement.rows()],
        )
        if abs(transport.det()) != 1:
            raise ArithmeticError(f"{record['label']} does not give a primitive U splitting")
        frame = -(complement * ns * complement.transpose())
        if not frame.is_positive_definite() or frame.det() != 948:
            raise ArithmeticError(f"{record['label']} has the wrong frame genus")
        root_count = int(pari(frame).qfminim(2)[0])
        if root_count:
            raise ArithmeticError(f"{record['label']} unexpectedly has {root_count} roots")
        is_published = qfisometric(frame, published)
        is_alternate = qfisometric(frame, alternate)
        if is_published == is_alternate:
            raise ArithmeticError(
                f"{record['label']} must match exactly one rootless J2 control"
            )
        classification = "published-R17" if is_published else "alternate-Q80"
        classifications.append(
            {
                "label": record["label"],
                "lattice_orbit_mask": int(record["lattice_orbit_mask"]),
                "trace_vector": list(map(int, w)),
                "isotropic_fibre": list(map(int, fibre)),
                "old_fibre_degree": 2,
                "old_zero_degree": 1,
                "shared_zero": True,
                "relative_u_cross_matrix": rows(cross),
                "transport_determinant": int(transport.det()),
                "frame_gram_sha256": matrix_digest(frame),
                "frame_class": classification,
                "equation_complexity": record["equation_complexity"],
            }
        )

    counts = Counter(row["frame_class"] for row in classifications)
    if counts != Counter({"published-R17": 33, "alternate-Q80": 10}):
        raise ArithmeticError(f"unexpected frame-class distribution {dict(counts)}")
    alternate_witnesses = sorted(
        (row for row in classifications if row["frame_class"] == "alternate-Q80"),
        key=lambda row: (
            row["equation_complexity"]["coefficient_l1"],
            row["equation_complexity"]["coordinate_input_bits"],
            row["label"],
        ),
    )

    bridge_gram = matrix(ZZ, [[2, 3], [3, 2]]).transpose() * hyperbolic * matrix(
        ZZ, [[2, 3], [3, 2]]
    ) - hyperbolic
    if bridge_gram != matrix(ZZ, [[12, 12], [12, 12]]) or bridge_gram.rank() != 1:
        raise ArithmeticError("shared-zero bridge regression failed")

    result = {
        "schema": "elkies-k3.r17-norm12-isotropic-frame-classification.v1",
        "status": "PASS_EXACT_MINIMAL_J2_ACCESSIBILITY",
        "theorem": {
            "source_frame_class": "published-R17",
            "target_frame_class": "alternate-Q80",
            "elliptic_incidence_distance": 2,
            "lower_bound": (
                "Distinct J2 frame classes cannot have fibre intersection 0 or 1; "
                "intersection 1 spans a U and an Eichler transvection identifies the complements."
            ),
            "attainment": (
                "Ten exact irreducible smooth genus-one bisections D=(3,2,w) are nef "
                "degree-two fibres with the old zero as a section, and their complements "
                "are integrally isometric to the alternate Q80 frame."
            ),
            "minimum_directed_zero_section_cost": 1,
            "shared_zero_attained": True,
            "relative_u_cross_matrix": [[2, 3], [3, 2]],
            "relative_bridge_gram": rows(bridge_gram),
            "relative_bridge_rank": 1,
        },
        "classification": {
            "norm_twelve_record_count": len(classifications),
            "frame_class_counts": dict(sorted(counts.items())),
            "alternate_witness_count": len(alternate_witnesses),
            "cheapest_alternate_witness": alternate_witnesses[0]["label"],
            "cheapest_alternate_equation_complexity": alternate_witnesses[0][
                "equation_complexity"
            ],
            "records": classifications,
        },
        "proof_boundary": (
            "The exact source artifact supplies the 43 irreducible smooth genus-one curves; "
            "this replay certifies their primitive shared-zero U splittings, rootless frames, "
            "and integral J2 classifications. It proves the minimum one-edge fibre-intersection "
            "and zero-section costs. It does not yet compile the two-dimensional |D| pencil, "
            "a Weierstrass equation, or a J1 surface-automorphism classification."
        ),
        "software_assumptions": {
            "sage_version": SAGE_VERSION,
            "pari_version": ".".join(map(str, pari.version())),
            "required_features": ["PARI qfminim", "PARI qfisom"],
        },
        "inputs": {
            relative(path): digest(path) for path in (PINNED, SPLITTING, ALTERNATE)
        },
        "reproducing_command": (
            "sage -python elkies-k3/scripts/classify_r17_norm12_isotropic_frames.sage"
        ),
    }
    serialized = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.exists() or args.output.read_text() != serialized:
            raise ArithmeticError("stored norm-twelve frame classification differs from replay")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(serialized)
    print(
        "R17NORM12FRAMES|records={}|published={}|alternate={}|distance=2|"
        "zero_cost=1|cheapest={}|output={}".format(
            len(classifications),
            counts["published-R17"],
            counts["alternate-Q80"],
            alternate_witnesses[0]["label"],
            relative(args.output),
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
