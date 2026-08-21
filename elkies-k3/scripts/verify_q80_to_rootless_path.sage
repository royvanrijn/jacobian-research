#!/usr/bin/env sage
"""Replay the exact q80 E6+D5+A3/MW3 to rootless/MW17 path.

The six neighbor witnesses were found by bounded searches, but this verifier
does not rerun those searches.  It reconstructs every child frame by exact
integral arithmetic, checks the complete root invariants, composes the six
Neron--Severi basis transports, and identifies the terminal frame with the
pinned determinant-948 rootless lattice by an explicit unimodular isometry.
"""

import csv
import hashlib
from pathlib import Path

from sage.all import (
    QQ,
    ZZ,
    block_diagonal_matrix,
    gcd,
    identity_matrix,
    matrix,
    pari,
    vector,
    xgcd,
)


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "elkies-k3" / "data" / "fibrations"
U = matrix(ZZ, [[0, 1], [1, 0]])


def load_matrix(path):
    return matrix(ZZ, [
        [ZZ(value) for value in line.split()]
        for line in Path(path).read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ])


def bezout_vector(pairings):
    current = ZZ(0)
    coefficients = [ZZ(0)] * len(pairings)
    for index, pairing in enumerate(pairings):
        if pairing == 0:
            continue
        new_gcd, left, right = xgcd(current, ZZ(pairing))
        coefficients = [left * value for value in coefficients]
        coefficients[index] += right
        current = new_gcd
    assert abs(current) == 1
    return vector(ZZ, coefficients if current == 1 else [-x for x in coefficients])


def neighbor(parent, qnorm, a, b, coordinates):
    ns = block_diagonal_matrix(U, -parent)
    fiber = vector(ZZ, [a, b] + list(coordinates))
    assert a * b == qnorm
    assert coordinates * parent * coordinates == 2 * qnorm
    assert fiber * ns * fiber == 0
    assert gcd([abs(ZZ(value)) for value in ns * fiber]) == 1
    mate = bezout_vector(list(ns * fiber))
    assert fiber * ns * mate == 1
    mate -= ZZ(mate * ns * mate) // 2 * fiber
    assert mate * ns * mate == 0 and fiber * ns * mate == 1
    complement = matrix(
        ZZ, [list(fiber * ns), list(mate * ns)]
    ).right_kernel_matrix()
    child = -(complement * ns * complement.transpose())
    transport = matrix(ZZ, [list(fiber), list(mate)] + complement.rows())
    assert abs(transport.det()) == 1
    assert transport * ns * transport.transpose() == block_diagonal_matrix(U, -child)
    return child, transport


def root_data(frame):
    result = pari(frame).qfminim(2)
    count = ZZ(result[0])
    if count == 0:
        return 0, 0, 1
    roots = matrix(ZZ, result[2]).transpose()
    basis = roots.row_module().basis_matrix()
    root_gram = basis * frame * basis.transpose()
    return basis.rank(), count, abs(root_gram.det())


def marked_mw_height(frame, ns_coordinates):
    """Project a marked NS direction to the current MW quotient."""
    frame_vector = vector(QQ, ns_coordinates[2:])
    root_result = pari(frame).qfminim(2)
    if ZZ(root_result[0]) == 0:
        projection = frame_vector
    else:
        roots = matrix(ZZ, root_result[2]).transpose()
        root_basis = roots.row_module().basis_matrix()
        root_gram = root_basis*frame*root_basis.transpose()
        root_coordinates = (
            frame_vector*frame*root_basis.transpose()*root_gram.inverse()
        )
        projection = frame_vector-root_coordinates*root_basis
    return projection*frame*projection


start = load_matrix(DATA / "kumar_q80_e6_d5_a3_mw3_frame.txt")
target = load_matrix(ROOT / "elkies-k3" / "data" / "lattice" / "rank17_gram.txt")
isometry = load_matrix(DATA / "kumar_q80_rootless_frame_isometry.txt")
pinned_transport = load_matrix(
    DATA / "kumar_q80_rootless_target_to_q80_ns_transport.txt"
)
assert start.det() == target.det() == 948
assert root_data(start) == (14, 124, 48)

# These are the old Kumar height-4 and level-79 frame directions transported
# into the initial q80 NS basis.  They are printed independently by
# classify_kumar_cm_frame_extensions.sage --print-markings-in-q80.  Their
# initial MW projections recover the certified q80 heights 4 and 120.
markings = {
    "height4": vector(ZZ, (
        -8, -8, 32, 48, 64, 96, 80, 64, 96, 144, 192, 288, 64, 56,
        -154, -67, 16, 30, -16,
    )),
    "Q79": vector(ZZ, (
        -123, -123, 492, 738, 984, 1476, 1230, 984, 1476, 2214, 2952,
        4428, 984, 861, -2368, -1027, 251, 461, -247,
    )),
}
assert marked_mw_height(start, markings["height4"]) == 4
assert marked_mw_height(start, markings["Q79"]) == 120
print(
    "Q80ROOTLESSMARKING|step=0|height4_MW_height=4|Q79_MW_height=120",
    flush=True,
)

with (DATA / "kumar_q80_to_rootless_path.tsv").open() as handle:
    steps = list(csv.DictReader(handle, delimiter="\t"))
assert len(steps) == 6

frame = start
composite = identity_matrix(ZZ, 19)
for row in steps:
    coordinates = vector(ZZ, map(ZZ, row["v"].split(",")))
    qnorm, a, b = map(ZZ, (row["q"], row["a"], row["b"]))
    frame, transition = neighbor(frame, qnorm, a, b, coordinates)
    composite = transition * composite
    markings = {
        name: vector(ZZ, marked*transition.inverse())
        for name, marked in markings.items()
    }
    actual = root_data(frame)
    expected = tuple(map(ZZ, (row["root_rank"], row["roots"], row["rootdet"])))
    assert actual == expected
    assert 17 - actual[0] == ZZ(row["MW"])
    print(
        f"Q80ROOTLESSPATH|step={row['step']}|q={qnorm}|ab={a},{b}|"
        f"ADE={row['ADE']}|root_rank={actual[0]}|roots={actual[1]}|"
        f"rootdet={actual[2]}|MW={row['MW']}|transport_det={transition.det()}",
        flush=True,
    )
    print(
        "Q80ROOTLESSPATH|step={}|q80_fiber={}".format(
            row["step"], tuple(map(ZZ, composite[0]))
        ),
        flush=True,
    )
    print(
        f"Q80ROOTLESSMARKING|step={row['step']}|"
        f"height4_MW_height={marked_mw_height(frame, markings['height4'])}|"
        f"Q79_MW_height={marked_mw_height(frame, markings['Q79'])}",
        flush=True,
    )

initial_ns = block_diagonal_matrix(U, -start)
terminal_ns = block_diagonal_matrix(U, -frame)
assert composite * initial_ns * composite.transpose() == terminal_ns
assert abs(isometry.det()) == 1
assert isometry.transpose() * target * isometry == frame

# qfisom is an independent exact existence check; its chosen matrix can vary
# across PARI versions, so the pinned matrix above is the transport authority.
pari_isometry = matrix(ZZ, pari(frame).qfisom(pari(target)))
assert abs(pari_isometry.det()) == 1
assert pari_isometry.transpose() * target * pari_isometry == frame

target_to_terminal = block_diagonal_matrix(
    identity_matrix(ZZ, 2), isometry.inverse().transpose().change_ring(ZZ)
)
target_to_q80 = target_to_terminal * composite
target_ns = block_diagonal_matrix(U, -target)
assert abs(target_to_q80.det()) == 1
assert target_to_q80 * initial_ns * target_to_q80.transpose() == target_ns
assert target_to_q80 == pinned_transport

transport_text = "\n".join(
    " ".join(map(str, target_to_q80[row]))
    for row in range(target_to_q80.nrows())
) + "\n"
transport_sha256 = hashlib.sha256(transport_text.encode()).hexdigest()

print(
    "Q80ROOTLESSPATH|terminal=rootless|MW=17|det=948|"
    f"composite_transport_det={target_to_q80.det()}|integral_isometry=1|"
    f"transport_sha256={transport_sha256}|status=PASS",
    flush=True,
)
print("Q80ROOTLESSPATH|target_to_q80_transport=", flush=True)
print(target_to_q80, flush=True)
