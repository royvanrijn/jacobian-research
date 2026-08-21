#!/usr/bin/env sage
"""Verify the equation-friendly alternate q4,q6 ending of the q80 path.

The first four transitions are the pinned q4,q4,q12,q12 chain.  The fifth
q4 class was selected by a bounded low-degree search, and the q6 class below
was found by a bounded small-q continuation.  This verifier rebuilds both
child frames exactly, proves that the q6 child is rootless, distinguishes it
from the previously pinned rank-17 frame, and checks the complete NS transport
from the original q80 frame.  The chosen q4 class is not yet identified with
the separate pair14 marked equation: their CM root data currently differ.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import (
    ZZ,
    block_diagonal_matrix,
    matrix,
    pari,
    vector,
)


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--write-artifact", action="store_true")
arguments = parser.parse_args()
load(str(HERE / "analyze_q80_fifth_q4_chamber.sage"))


def load_matrix(path):
    return matrix(ZZ, [
        [ZZ(value) for value in line.split()]
        for line in Path(path).read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ])


def exact_root_data(frame):
    result = pari(frame).qfminim(2)
    count = ZZ(result[0])
    if not count:
        return (0, 0, 1)
    roots = matrix(ZZ, result[2]).transpose()
    basis = roots.row_module().basis_matrix()
    gram = basis*frame*basis.transpose()
    return (int(basis.rank()), int(count), int(abs(gram.det())))


alternate_q4_v = vector(ZZ, (
    -9, 8, -11, 10, -4, 0, 5, 1, -6,
    6, 1, -2, -1, -1, 1, 2, 0,
))
alternate_child, alternate_transport = neighbor(
    fourth_child_frame, ZZ(4), ZZ(2), ZZ(2), alternate_q4_v
)
assert exact_root_data(alternate_child) == (1, 2, 2)

q6_v = vector(ZZ, (
    0, -2, 4, 2, -1, 2, 1, -1, 1,
    0, 1, -1, 1, 0, 0, 0, 0,
))
rootless_child, q6_transport = neighbor(
    alternate_child, ZZ(6), ZZ(2), ZZ(3), q6_v
)
assert exact_root_data(rootless_child) == (0, 0, 1)
assert rootless_child.det() == 948

target = load_matrix(ROOT / "elkies-k3/data/lattice/rank17_gram.txt")
pinned_isometry = pari(rootless_child).qfisom(pari(target))
assert pinned_isometry == 0

q80_to_fourth = (
    fourth_child_transport
    * third_child_transport
    * second_transport
    * first_transport
)
composite = q6_transport*alternate_transport*q80_to_fourth
q80_ns = block_diagonal_matrix(U, -start)
rootless_ns = block_diagonal_matrix(U, -rootless_child)
assert composite*q80_ns*composite.transpose() == rootless_ns
assert abs(composite.det()) == 1


def canonical_matrix_text(value):
    return "\n".join(
        " ".join(map(str, row)) for row in value.rows()
    )+"\n"


frame_sha256 = hashlib.sha256(
    canonical_matrix_text(rootless_child).encode()
).hexdigest()
transport_sha256 = hashlib.sha256(
    canonical_matrix_text(composite).encode()
).hexdigest()

print(
    "Q80ALTFIFTHROOTLESS|"
    f"q4_v={tuple(alternate_q4_v)}|q4_child=A1/MW16|"
    f"q6_v={tuple(q6_v)}|q6_child=rootless/MW17|"
    f"q4_transport_det={alternate_transport.det()}|"
    f"q6_transport_det={q6_transport.det()}|"
    "pinned_rootless_isometric=0|"
    f"composite_det={composite.det()}|"
    f"frame_sha256={frame_sha256}|transport_sha256={transport_sha256}|"
    "status=PASS",
    flush=True,
)
print("Q80ALTFIFTHROOTLESS|rootless_frame=", flush=True)
print(rootless_child, flush=True)
print("Q80ALTFIFTHROOTLESS|rootless_to_q80=", flush=True)
print(composite, flush=True)

if arguments.write_artifact:
    output = {
        "schema": "q80-alternate-fifth-q6-rootless-transport-v1",
        "q4": {
            "q": 4,
            "a": 2,
            "b": 2,
            "v": list(map(int, alternate_q4_v)),
            "child_root_data": [1, 2, 2],
            "child_MW_rank": 16,
        },
        "q6": {
            "q": 6,
            "a": 2,
            "b": 3,
            "v": list(map(int, q6_v)),
            "child_root_data": [0, 0, 1],
            "child_MW_rank": 17,
        },
        "rootless_frame": [list(map(int, row)) for row in rootless_child.rows()],
        "rootless_frame_sha256": frame_sha256,
        "rootless_to_q80_ns_transport": [
            list(map(int, row)) for row in composite.rows()
        ],
        "rootless_to_q80_ns_transport_sha256": transport_sha256,
        "transport_determinant": int(composite.det()),
        "isometric_to_previously_pinned_rootless_frame": False,
        "reproduce": (
            "sage elkies-k3/scripts/"
            "verify_q80_alternate_fifth_q6_rootless.sage --write-artifact"
        ),
    }
    output_path = (
        ROOT / "artifacts/generated-results/"
        "q80-alternate-fifth-q6-rootless-transport.json"
    )
    encoded = json.dumps(output, indent=2, sort_keys=True, default=int)+"\n"
    output_path.write_text(encoded)
    artifact_sha256 = hashlib.sha256(encoded.encode()).hexdigest()
    print(
        "Q80ALTFIFTHROOTLESS|"
        f"artifact={output_path}|sha256={artifact_sha256}|"
        "status=PASS_ARTIFACT_WRITE",
        flush=True,
    )
