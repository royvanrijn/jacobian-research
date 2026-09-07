#!/usr/bin/env sage-python
"""Classify the elliptic frame orthogonal to the ``0x103b2`` bisection class.

The norm-twelve trace representative ``w`` defines the primitive isotropic
class ``D=(3,2,w)`` in ``NS=U+R17(-1)``.  This replay splits off an integral
hyperbolic plane containing ``D``, enumerates all roots in the resulting
positive frame, and compares that frame with the two certified rootless J2
classes of determinant 948.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import ZZ, block_diagonal_matrix, gcd, matrix, pari, vector, xgcd
from sage.env import SAGE_VERSION


ROOT = Path(__file__).resolve().parents[2]
PINNED = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
SPLITTING = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-genus-one-bisection-splitting-search-v1.json"
)
ALTERNATE = ROOT / "artifacts/generated-results/q80-alternate-fifth-q6-rootless-transport.json"
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-103b2-isotropic-frame-v1.json"
)
TARGET_LABEL = "norm12-orbit-103b2"


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


def bezout_vector(pairings):
    current = ZZ(0)
    coefficients = [ZZ(0)] * len(pairings)
    for index, pairing in enumerate(pairings):
        pairing = ZZ(pairing)
        if not pairing:
            continue
        divisor, left, right = xgcd(current, pairing)
        coefficients = [left * coefficient for coefficient in coefficients]
        coefficients[index] += right
        current = divisor
    if current == -1:
        coefficients = [-coefficient for coefficient in coefficients]
        current = -current
    if current != 1:
        raise ArithmeticError(f"isotropic class has divisibility {current}, not one")
    return vector(ZZ, coefficients)


def rows(value):
    return [list(map(int, row)) for row in value.rows()]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    pinned = load_matrix(PINNED)
    splitting = json.loads(SPLITTING.read_text())
    target = next(
        record
        for record in splitting["construction"]["records"]
        if record["label"] == TARGET_LABEL
    )
    w = vector(ZZ, target["pinned_rank17_w"])
    if w * pinned * w != 12:
        raise ArithmeticError("the stored 0x103b2 representative no longer has norm twelve")

    hyperbolic = matrix(ZZ, [[0, 1], [1, 0]])
    ns = block_diagonal_matrix(hyperbolic, -pinned)
    old_fibre = vector(ZZ, [1, 0] + [0] * 17)
    old_zero = vector(ZZ, [-1, 1] + [0] * 17)
    fibre = vector(ZZ, [3, 2] + list(w))
    if fibre * ns * fibre != 0 or fibre * ns * old_fibre != 2:
        raise ArithmeticError("D=(3,2,w) is not the expected degree-two isotropic class")

    divisibility = gcd([abs(ZZ(value)) for value in ns * fibre])
    if divisibility != 1:
        raise ArithmeticError(f"D has divisibility {divisibility}, so it does not split an integral U")
    mate = bezout_vector(list(ns * fibre))
    mate -= ZZ(mate * ns * mate) // 2 * fibre
    if fibre * ns * mate != 1 or mate * ns * mate != 0:
        raise ArithmeticError("failed to construct an isotropic mate")

    complement = matrix(
        ZZ, [list(fibre * ns), list(mate * ns)]
    ).right_kernel_matrix()
    transport = matrix(
        ZZ,
        [list(fibre), list(mate)] + [list(row) for row in complement.rows()],
    )
    if abs(transport.det()) != 1:
        raise ArithmeticError("the new U plus frame basis is not unimodular")
    frame = -(complement * ns * complement.transpose())
    if not frame.is_positive_definite() or frame.det() != pinned.det():
        raise ArithmeticError("the orthogonal frame has the wrong signature or determinant")

    root_query = pari(frame).qfminim(2)
    root_count = int(root_query[0])
    root_rank = (
        matrix(ZZ, root_query[2]).rank()
        if root_count
        else 0
    )
    published_isometry = pari(frame).qfisom(pari(pinned))

    alternate_data = json.loads(ALTERNATE.read_text())
    alternate = matrix(ZZ, alternate_data["rootless_frame"])
    alternate_isometry = pari(frame).qfisom(pari(alternate))
    if root_count != 0:
        identification = "rootful"
    elif published_isometry != 0:
        identification = "published-R17-J2-class"
    elif alternate_isometry != 0:
        identification = "alternate-Q80-rootless-J2-class"
    else:
        identification = "unidentified-rootless-J2-class"

    result = {
        "schema": "elkies-k3.r17-norm12-103b2-isotropic-frame.v1",
        "status": "PASS_EXACT_ISOTROPIC_FRAME_CLASSIFICATION",
        "source": {
            "cover_label": TARGET_LABEL,
            "trace_coset": "0x103b2",
            "trace_representative_in_pinned_r17": list(map(int, w)),
            "trace_norm": int(w * pinned * w),
            "isotropic_class_D_in_pinned_ns": list(map(int, fibre)),
            "D_squared": int(fibre * ns * fibre),
            "D_old_fibre_degree": int(fibre * ns * old_fibre),
            "D_old_zero_pairing": int(fibre * ns * old_zero),
            "D_divisibility": int(divisibility),
        },
        "new_marking": {
            "isotropic_mate_in_pinned_ns": list(map(int, mate)),
            "zero_in_pinned_ns": list(map(int, mate - fibre)),
            "basis_to_pinned_ns": rows(transport),
            "transport_determinant": int(transport.det()),
        },
        "orthogonal_positive_frame": {
            "rank": int(frame.rank()),
            "determinant": int(frame.det()),
            "norm_two_root_count": root_count,
            "root_lattice_rank": int(root_rank),
            "gram": rows(frame),
            "gram_sha256": hashlib.sha256(str(rows(frame)).encode()).hexdigest(),
            "integrally_isometric_to_published_r17": published_isometry != 0,
            "integrally_isometric_to_alternate_q80_rootless_frame": alternate_isometry != 0,
            "identification": identification,
            "mw_rank_if_rho_19": 17 - int(root_rank),
        },
        "comparison": {
            "published_isometry_matrix": None if published_isometry == 0 else rows(matrix(ZZ, published_isometry)),
            "alternate_q80_isometry_matrix": None if alternate_isometry == 0 else rows(matrix(ZZ, alternate_isometry)),
        },
        "claim_boundary": (
            "This is an exact lattice/J2 classification of the primitive isotropic class D. "
            "It does not construct the second fibration equation over QQ, identify the pointed "
            "quartic with a fibre of that fibration, or prove an upper bound for its Jacobian rank."
        ),
        "software_assumptions": {
            "sage_version": SAGE_VERSION,
            "pari_version": ".".join(map(str, pari.version())),
            "required_features": ["PARI qfminim", "PARI qfisom"],
        },
        "inputs": {relative(path): digest(path) for path in (PINNED, SPLITTING, ALTERNATE)},
        "reproducing_command": (
            "sage -python elkies-k3/scripts/classify_r17_103b2_isotropic_frame.sage"
        ),
    }
    serialized = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.exists() or args.output.read_text() != serialized:
            raise ArithmeticError("stored isotropic-frame certificate differs from replay")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(serialized)
    print(
        f"R17103B2FRAME|roots={root_count}|det={frame.det()}|"
        f"published={int(published_isometry != 0)}|alternate_q80={int(alternate_isometry != 0)}|"
        f"identification={identification}|output={relative(args.output)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
