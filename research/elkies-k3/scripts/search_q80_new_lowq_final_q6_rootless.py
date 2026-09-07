#!/usr/bin/env python3
from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path

SAGE_CODE = r"""
from __future__ import print_function

import csv
import sys
from math import isqrt
from collections import defaultdict
from pathlib import Path
from time import perf_counter

from sage.all import (
    QQ, ZZ, QuadraticForm, block_diagonal_matrix, ceil, floor, gcd, identity_matrix,
    lcm, matrix, vector, xgcd
)

ROOT = Path(sys.argv[1]).resolve()
CERT = Path(sys.argv[2]).resolve()
DATA = ROOT / "elkies-k3" / "data" / "fibrations"
U = matrix(ZZ, [[0, 1], [1, 0]])


def load_matrix(path):
    return matrix(
        ZZ,
        [[ZZ(x) for x in line.split()]
         for line in Path(path).read_text().splitlines()
         if line.strip() and not line.lstrip().startswith("#")],
    )


def bezout_vector(pairings):
    current = ZZ(0)
    coeffs = [ZZ(0)] * len(pairings)
    for i, pairing in enumerate(pairings):
        if not pairing:
            continue
        g, left, right = xgcd(current, ZZ(pairing))
        coeffs = [left * c for c in coeffs]
        coeffs[i] += right
        current = g
    assert abs(current) == 1
    return vector(ZZ, coeffs if current == 1 else [-c for c in coeffs])


def neighbor(parent, qnorm, a, b, coordinates):
    ns = block_diagonal_matrix(U, -parent)
    coordinates = vector(ZZ, coordinates)
    F = vector(ZZ, [a, b] + list(coordinates))
    assert a * b == qnorm
    assert coordinates * parent * coordinates == 2 * qnorm
    assert F * ns * F == 0

    # For (a,b)=(2,3), primitivity is automatic since NS.F contains 2 and 3.
    assert gcd([abs(ZZ(x)) for x in ns * F]) == 1

    mate = bezout_vector(list(ns * F))
    mate -= ZZ(mate * ns * mate) // 2 * F
    complement = matrix(
        ZZ, [list(F * ns), list(mate * ns)]
    ).right_kernel_matrix()

    child = -(complement * ns * complement.transpose())
    transport = matrix(
        ZZ, [list(F), list(mate)] + complement.rows()
    )
    assert abs(transport.det()) == 1
    return child, transport


def roots_of_norm_two(gram):
    half = QuadraticForm(ZZ, gram).short_vector_list_up_to_length(
        2, up_to_sign_flag=True
    )[1]
    half = [vector(ZZ, row) for row in half]
    return half + [-row for row in half]


def lex_positive(row):
    return next(value > 0 for value in row if value)


def deterministic_simple_roots(gram):
    roots = roots_of_norm_two(gram)
    positive = [r for r in roots if lex_positive(r)]
    positive_set = {tuple(r) for r in positive}
    simple = []
    for r in positive:
        if not any(
            tuple(r - s) in positive_set
            for s in positive if s != r
        ):
            simple.append(r)
    simple = matrix(ZZ, [list(r) for r in simple])
    assert simple.nrows() == simple.rank()
    return simple, positive


def connected_components(cartan):
    unseen = set(range(cartan.nrows()))
    answer = []
    while unseen:
        first = min(unseen)
        unseen.remove(first)
        todo = [first]
        component = []
        while todo:
            i = todo.pop()
            component.append(i)
            for j in list(unseen):
                if cartan[i, j]:
                    unseen.remove(j)
                    todo.append(j)
        answer.append(tuple(sorted(component)))
    return tuple(sorted(answer, key=lambda c: (len(c), c)))


def ade_component_name(rank, det):
    rank = int(rank)
    det = int(abs(det))
    if det == rank + 1:
        return f"A{rank}"
    if rank >= 4 and det == 4:
        return f"D{rank}"
    if (rank, det) == (6, 3):
        return "E6"
    if (rank, det) == (7, 2):
        return "E7"
    if (rank, det) == (8, 1):
        return "E8"
    raise ArithmeticError(
        f"Unrecognized ADE component rank={rank} det={det}"
    )


def root_signature(gram):
    simple, positive = deterministic_simple_roots(gram)
    rank = simple.nrows()
    if rank == 0:
        return "rootless", 0, 0, ZZ(1)
    cartan = simple * gram * simple.transpose()
    components = connected_components(cartan)
    parts = []
    for component in components:
        sub = cartan.matrix_from_rows_and_columns(component, component)
        parts.append(
            (len(component),
             ade_component_name(len(component), sub.det()))
        )
    parts.sort(key=lambda x: (-x[0], x[1]))
    return (
        "+".join(name for _, name in parts),
        rank,
        2 * len(positive),
        abs(cartan.det()),
    )


def is_rootless(gram):
    # Exact and cheaper than building a full ADE signature for every child.
    shell = QuadraticForm(ZZ, gram).short_vector_list_up_to_length(
        2, up_to_sign_flag=True
    )[1]
    return len(shell) == 0


def dominant_weights(cartan, max_norm):
    inv = cartan.inverse()
    bounds = []
    for i in range(cartan.nrows()):
        diag = QQ(inv[i, i])
        b = ZZ(0)
        while QQ((b + 1) ** 2) * diag <= max_norm:
            b += 1
        bounds.append(int(b))

    out = []

    def rec(prefix, index):
        if index == cartan.nrows():
            p = vector(ZZ, prefix)
            norm = QQ(p * inv * p)
            if norm <= max_norm:
                out.append((tuple(p), norm))
            return
        for value in range(bounds[index] + 1):
            rec(prefix + [value], index + 1)

    rec([], 0)
    out.sort(key=lambda item: (item[1], item[0]))
    return tuple(out)



def quotient_data(frame):
    simple, _ = deterministic_simple_roots(frame)
    assert simple.nrows() == 1
    root_gram = simple * frame * simple.transpose()
    assert root_gram == matrix(ZZ, [[2]])

    projection = (
        identity_matrix(QQ, 17)
        - frame * simple.transpose() * root_gram.inverse() * simple
    )
    den = lcm(x.denominator() for x in projection.list())
    scaled_projection = (den * projection).change_ring(ZZ)

    projected_integer = scaled_projection.row_module().basis_matrix()
    mw_rank = 17 - simple.nrows()
    assert mw_rank == 16
    assert projected_integer.nrows() == mw_rank

    mw_basis = projected_integer / den
    mw_height = mw_basis * frame * mw_basis.transpose()

    hden = lcm(x.denominator() for x in mw_height.list())
    scaled_height = (2 * hden * mw_height).change_ring(ZZ)
    transform = scaled_height.LLL_gram().transpose()
    assert abs(transform.det()) == 1
    mw_basis = transform * mw_basis
    mw_height = mw_basis * frame * mw_basis.transpose()

    hden = lcm(x.denominator() for x in mw_height.list())
    mw_form = QuadraticForm(
        ZZ, (2 * hden * mw_height).change_ring(ZZ)
    )

    linear_map = scaled_projection.transpose()
    diagonal, left, right = linear_map.smith_form()
    assert left * linear_map * right == diagonal

    def integral_preimage(projected):
        target = vector(QQ, projected) * den
        assert all(x in ZZ for x in target)
        target = vector(ZZ, target)
        transformed = left * target
        smith = vector(ZZ, [0] * 17)
        for i in range(17):
            d = diagonal[i, i]
            if d:
                assert transformed[i] % d == 0
                smith[i] = transformed[i] // d
            else:
                assert transformed[i] == 0
        pre = right * smith
        assert pre * scaled_projection == target
        assert vector(QQ, pre) * projection == vector(QQ, projected)
        return vector(ZZ, pre)

    mw_pre = matrix(
        ZZ, [list(integral_preimage(row)) for row in mw_basis.rows()]
    )
    assert mw_pre * projection == mw_basis

    return {
        "simple": simple,
        "root_inverse": root_gram.inverse(),
        "projection": projection,
        "mw_basis": mw_basis,
        "mw_height": mw_height,
        "mw_form": mw_form,
        "height_den": hden,
        "mw_pre": mw_pre,
    }


def pari_columns(pari_matrix):
    # qfminim returns a PARI matrix whose columns are the stored vectors.
    # Iteration over the matrix gives its columns in cypari/Sage.
    for column in pari_matrix:
        yield vector(ZZ, [ZZ(x) for x in column])


def search_pairing_shell(frame, data, pairing, store_cap, started):
    # For A1 with dominant root pairing p, root correction norm is p^2/2.
    # Full target norm is 12, so required MW height is 12-p^2/2.
    p = ZZ(pairing)
    target_height = QQ(12) - (QQ(p * p) / QQ(2))
    assert target_height >= 0

    hden = ZZ(data["height_den"])
    # QuadraticForm(mw_form) has Q(z)=hden*MWheight.
    target_q = target_height * hden
    assert target_q in ZZ
    target_q = ZZ(target_q)

    # PARI sees the even Gram matrix, hence norm = 2*Q(z).
    pari_bound = ZZ(2 * target_q)

    print(
        f"Q80A1Q6SHELL|p={p}|target_MW_height={target_height}|"
        f"pari_bound={pari_bound}|store_cap_pairs={store_cap}|"
        "status=START_CAPPED_QFMINIM",
        flush=True,
    )

    # Critical memory fix: explicit m prevents PARI from materializing every
    # vector in a potentially enormous radius ball.  qfminim still enumerates
    # the full ball and reports its exact total count.
    result = data["mw_form"].__pari__().qfminim(
        pari_bound, ZZ(store_cap)
    )
    total_signed = ZZ(result[0])
    max_norm_seen = result[1]
    stored_matrix = result[2]

    stored = list(pari_columns(stored_matrix))
    total_pairs = total_signed // 2
    exhaustive_ball = total_pairs <= store_cap

    print(
        f"Q80A1Q6SHELL|p={p}|qfminim_signed_count={total_signed}|"
        f"pair_count_ball={total_pairs}|stored_pairs={len(stored)}|"
        f"exhaustive_ball={int(exhaustive_ball)}|"
        f"max_norm_seen={max_norm_seen}|"
        f"elapsed={perf_counter()-started:.1f}s|"
        "status=PASS_CAPPED_QFMINIM",
        flush=True,
    )

    simple = data["simple"]
    root_inverse = data["root_inverse"]
    projection = data["projection"]
    mw_basis = data["mw_basis"]
    mw_height = data["mw_height"]
    mw_pre = data["mw_pre"]

    tested_shell = 0
    compatible = 0

    # qfminim stores vectors in no guaranteed order.  We therefore filter
    # exactly by the requested shell norm rather than assuming its tail.
    for z in stored:
        q_value = data["mw_form"](z)
        if q_value != target_q:
            continue

        tested_shell += 1
        mw = vector(QQ, z) * mw_basis
        assert QQ(mw * frame * mw) == target_height

        x0 = z * mw_pre
        assert vector(QQ, x0) * projection == vector(QQ, mw)
        p0 = x0 * frame * simple.transpose()
        assert len(p0) == 1 and p0[0] in ZZ
        p0 = vector(ZZ, p0)

        target_pairing = vector(ZZ, [p])
        delta = (
            (vector(QQ, target_pairing) - vector(QQ, p0))
            * root_inverse
        )
        if not all(x in ZZ for x in delta):
            continue

        root_part = (
            vector(QQ, target_pairing) * root_inverse * simple
        )
        vq = vector(QQ, mw) + root_part
        if not all(x in ZZ for x in vq):
            continue

        v = vector(ZZ, vq)
        assert v * frame * v == 12
        compatible += 1

        child, _ = neighbor(
            frame, ZZ(6), ZZ(2), ZZ(3), v
        )
        if is_rootless(child):
            assert child.det() == 948
            assert root_signature(child)[:2] == ("rootless", 0)
            return {
                "found": True,
                "p": p,
                "v": tuple(v),
                "child": child,
                "total_pairs_ball": total_pairs,
                "stored_pairs": len(stored),
                "tested_shell": tested_shell,
                "compatible": compatible,
                "exhaustive_ball": exhaustive_ball,
            }

    return {
        "found": False,
        "p": p,
        "total_pairs_ball": total_pairs,
        "stored_pairs": len(stored),
        "tested_shell": tested_shell,
        "compatible": compatible,
        "exhaustive_ball": exhaustive_ball,
    }



def exact_ldl(gram):
    # Exact rational LDL^T decomposition gram = L*D*L^T.
    gram = matrix(QQ, gram)
    n = gram.nrows()
    assert gram.ncols() == n
    L = identity_matrix(QQ, n)
    D = [QQ(0)] * n

    for j in range(n):
        d = gram[j, j]
        for k in range(j):
            d -= L[j, k] * L[j, k] * D[k]
        assert d > 0
        D[j] = d

        for i in range(j + 1, n):
            value = gram[i, j]
            for k in range(j):
                value -= L[i, k] * L[j, k] * D[k]
            L[i, j] = value / d

    Dmat = matrix.diagonal(D)
    assert L * Dmat * L.transpose() == gram
    return L, tuple(D)


def ceil_sqrt_rational(value):
    value = QQ(value)
    assert value >= 0
    numerator = ZZ(value.numerator())
    denominator = ZZ(value.denominator())
    base = ZZ(isqrt(int(numerator // denominator)))
    while base * base * denominator < numerator:
        base += 1
    while base > 0 and (base - 1) * (base - 1) * denominator >= numerator:
        base -= 1
    return base


def exact_shell_stream(gram, target_norm):
    # Stream one representative for every +/- integral vector z satisfying
    # z^T gram z == target_norm.  All pruning is exact over QQ.
    gram = matrix(QQ, gram)
    target_norm = QQ(target_norm)
    assert target_norm > 0

    L, D = exact_ldl(gram)
    n = gram.nrows()
    z = [ZZ(0)] * n
    emitted = 0

    def recurse(i, partial):
        nonlocal emitted
        if i < 0:
            if partial != target_norm:
                return
            zv = vector(ZZ, z)
            assert QQ(zv * gram * zv) == target_norm
            # Global sign quotient: first nonzero coordinate positive.
            first = next((entry for entry in zv if entry), None)
            assert first is not None
            if first < 0:
                return
            emitted += 1
            yield zv
            return

        remaining = target_norm - partial
        if remaining < 0:
            return

        # (L^T z)_i = z_i + sum_{j>i} L[j,i] z_j.
        center_shift = QQ(0)
        for j in range(i + 1, n):
            center_shift += L[j, i] * z[j]

        allowed_square = remaining / D[i]
        radius = ceil_sqrt_rational(allowed_square)
        center = -center_shift
        low = ZZ(floor(center)) - radius
        high = ZZ(ceil(center)) + radius

        for zi in range(int(low), int(high) + 1):
            yi = QQ(zi) + center_shift
            contribution = D[i] * yi * yi
            if contribution > remaining:
                continue
            z[i] = ZZ(zi)
            yield from recurse(i - 1, partial + contribution)

        z[i] = ZZ(0)

    yield from recurse(n - 1, QQ(0))


def pari_signed_ball_count(mw_form, pari_norm_bound):
    # m=0 stores no vector pairs but still enumerates and reports exact count.
    result = mw_form.__pari__().qfminim(ZZ(pari_norm_bound), ZZ(0))
    assert len(result) == 3
    return ZZ(result[0])


def search_exact_pairing_shell(frame, data, pairing, cert, started):
    p = ZZ(pairing)
    target_height = QQ(12) - (QQ(p * p) / QQ(2))
    assert target_height >= 0

    hden = ZZ(data["height_den"])
    target_q = target_height * hden
    assert target_q in ZZ
    target_q = ZZ(target_q)

    # mw_form has Gram 2*hden*H, while QuadraticForm evaluates Q=xGx/2.
    # Our streamer works with the Gram norm xGx = 2*target_q.
    target_gram_norm = ZZ(2 * target_q)

    simple = data["simple"]
    root_inverse = data["root_inverse"]
    projection = data["projection"]
    mw_basis = data["mw_basis"]
    mw_pre = data["mw_pre"]
    gram = matrix(ZZ, data["mw_form"].matrix())

    print(
        f"Q80A1STREAM|p={p}|target_MW_height={target_height}|"
        f"target_Q={target_q}|target_gram_norm={target_gram_norm}|"
        "status=START_EXACT_STREAM",
        flush=True,
    )

    shell_pairs = 0
    compatible = 0

    for z in exact_shell_stream(gram, target_gram_norm):
        shell_pairs += 1

        mw = vector(QQ, z) * mw_basis
        assert QQ(mw * frame * mw) == target_height

        x0 = z * mw_pre
        assert vector(QQ, x0) * projection == vector(QQ, mw)
        p0 = x0 * frame * simple.transpose()
        assert len(p0) == 1 and p0[0] in ZZ
        p0 = vector(ZZ, p0)

        target_pairing = vector(ZZ, [p])
        delta = (
            (vector(QQ, target_pairing) - vector(QQ, p0))
            * root_inverse
        )
        if not all(entry in ZZ for entry in delta):
            continue

        root_part = (
            vector(QQ, target_pairing) * root_inverse * simple
        )
        vq = vector(QQ, mw) + root_part
        if not all(entry in ZZ for entry in vq):
            continue

        v = vector(ZZ, vq)
        assert v * frame * v == 12
        compatible += 1

        child, _ = neighbor(
            frame, ZZ(6), ZZ(2), ZZ(3), v
        )
        if is_rootless(child):
            assert child.det() == 948
            assert root_signature(child)[:2] == ("rootless", 0)

            cert.parent.mkdir(parents=True, exist_ok=True)
            cert.write_text(
                "source=A1/MW16 candidate1\n"
                "q=6\n"
                "a=2\n"
                "b=3\n"
                f"A1_pairing={p}\n"
                "v=" + ",".join(map(str, tuple(v))) + "\n"
                "child=rootless/MW17\n"
                "det=948\n"
            )

            return {
                "found": True,
                "p": p,
                "v": tuple(v),
                "shell_pairs_seen": shell_pairs,
                "compatible": compatible,
            }

        if compatible <= 5 or compatible % 5000 == 0:
            print(
                f"Q80A1STREAM|p={p}|shell_pairs_seen={shell_pairs}|"
                f"gluing_compatible={compatible}|"
                f"elapsed={perf_counter()-started:.1f}s",
                flush=True,
            )

    return {
        "found": False,
        "p": p,
        "shell_pairs_seen": shell_pairs,
        "compatible": compatible,
    }


# Reconstruct the new low-q path through candidate-1 A1/MW16.
with (DATA / "kumar_q80_to_rootless_path.tsv").open() as handle:
    steps = list(csv.DictReader(handle, delimiter="\t"))
assert len(steps) >= 6
start = load_matrix(DATA / "kumar_q80_e6_d5_a3_mw3_frame.txt")

# Reconstruct the canonical A1/MW16 endpoint explicitly.
canon_a1 = start
for canonical_step in steps[:5]:
    canonical_v = vector(
        ZZ, map(ZZ, canonical_step["v"].split(","))
    )
    canon_a1, _ = neighbor(
        canon_a1,
        ZZ(canonical_step["q"]),
        ZZ(canonical_step["a"]),
        ZZ(canonical_step["b"]),
        canonical_v,
    )
canon_q6_step = steps[5]
canon_q6_v = vector(
    ZZ, map(ZZ, canon_q6_step["v"].split(","))
)
assert root_signature(canon_a1)[:2] == ("A1", 1)
assert canon_q6_v * canon_a1 * canon_q6_v == 12
assert (
    ZZ(canon_q6_step["q"]),
    ZZ(canon_q6_step["a"]),
    ZZ(canon_q6_step["b"]),
) == (6, 2, 3)

first = steps[0]
first_frame, _ = neighbor(
    start,
    ZZ(first["q"]), ZZ(first["a"]), ZZ(first["b"]),
    vector(ZZ, map(ZZ, first["v"].split(","))),
)
second = steps[1]
second_frame, _ = neighbor(
    first_frame,
    ZZ(second["q"]), ZZ(second["a"]), ZZ(second["b"]),
    vector(ZZ, map(ZZ, second["v"].split(","))),
)

escape_v = vector(
    ZZ,
    (-5,-3,6,6,-8,-4,2,4,-1,8,-16,-1,0,3,5,-2,-2),
)
q6_frame, _ = neighbor(
    second_frame, ZZ(6), ZZ(2), ZZ(3), escape_v
)

orbit424_v = vector(
    ZZ,
    (32,48,-21,28,8,-52,-34,0,18,5,-23,43,9,-18,16,-6,-6),
)
a6a4_frame, _ = neighbor(
    q6_frame, ZZ(4), ZZ(2), ZZ(2), orbit424_v
)

orbit1222_v = vector(
    ZZ,
    (10,53,-192,-114,29,-256,-170,-12,-14,74,-32,-14,-6,-26,-58,84,-28),
)
a6a3_frame, _ = neighbor(
    a6a4_frame, ZZ(4), ZZ(2), ZZ(2), orbit1222_v
)

v7774 = vector(
    ZZ,
    (85,2699,1257,7718,3756,-41,3077,-4614,-6615,6032,2584,-1678,121,-736,-913,1,1165),
)
frame7774, _ = neighbor(
    a6a3_frame, ZZ(6), ZZ(2), ZZ(3), v7774
)

v1938 = vector(
    ZZ,
    (-94,-1003,5298,4977,-1431,-1440,100,1,-1632,1893,1634,-1264,-4175,2248,-3111,1561,2842),
)
frame1938, _ = neighbor(
    frame7774, ZZ(4), ZZ(2), ZZ(2), v1938
)

v6855 = vector(
    ZZ,
    (30693,-339,-2534,45446,10413,16390,-11527,5970,-18424,4193,21146,11296,25035,17925,-6032,4304,7717),
)
frame6855, _ = neighbor(
    frame1938, ZZ(4), ZZ(2), ZZ(2), v6855
)

a1_v1 = vector(
    ZZ,
    (21,671,-20182,-10366,27727,30558,5582,20831,-10195,-19691,6086,10389,20928,18651,16123,15473,-11496),
)
frame_a1, _ = neighbor(
    frame6855, ZZ(4), ZZ(2), ZZ(2), a1_v1
)

source_sig = root_signature(frame_a1)
assert source_sig[:2] == ("A1", 1), source_sig
assert frame_a1.det() == 948

print(
    f"Q80A1Q6|source=candidate1|signature={source_sig}|MW=16|"
    "status=PASS_A1_SOURCE",
    flush=True,
)


# Determine the A1 pairing used by the canonical successful q6 move.
canon_simple, _ = deterministic_simple_roots(canon_a1)
assert canon_simple.nrows() == 1
canon_pairings = canon_q6_v * canon_a1 * canon_simple.transpose()
assert len(canon_pairings) == 1
canon_pairing_raw = ZZ(canon_pairings[0])
canon_pairing = abs(canon_pairing_raw)
assert 0 <= canon_pairing <= 4

print(
    f"Q80A1Q6|canonical_success_A1_pairing={canon_pairing}|"
    f"raw_pairing={canon_pairing_raw}|"
    "status=PASS_CANONICAL_PAIRING_GUIDE",
    flush=True,
)

data = quotient_data(frame_a1)
assert data["simple"].nrows() == 1
assert data["height_den"] == 2
assert data["mw_basis"].nrows() == 16
assert data["mw_basis"].ncols() == 17

print(
    f"Q80A1STREAM|candidate1_MW_rank=16|height_den={data['height_den']}|"
    f"canonical_success_A1_pairing={canon_pairing}|"
    "status=PASS_STRUCTURAL_SETUP",
    flush=True,
)

# ------------------------------------------------------------------
# Self-certify the exact streamer on the small p=4 shell.
# p=4 => MW height 4 => Q=8 => Gram norm 16.
# PARI's N is signed, while exact_shell_stream yields one +/- representative.
# ------------------------------------------------------------------
mw_gram = matrix(ZZ, data["mw_form"].matrix())
pari_le_16 = pari_signed_ball_count(data["mw_form"], ZZ(16))
pari_le_14 = pari_signed_ball_count(data["mw_form"], ZZ(14))
pari_exact_signed = pari_le_16 - pari_le_14
assert pari_exact_signed >= 0 and pari_exact_signed % 2 == 0
pari_exact_pairs = pari_exact_signed // 2

stream_p4 = list(exact_shell_stream(mw_gram, ZZ(16)))
assert len(stream_p4) == pari_exact_pairs, (
    len(stream_p4), pari_exact_pairs
)
assert len({tuple(z) for z in stream_p4}) == len(stream_p4)
assert all(
    ZZ(vector(ZZ, z) * mw_gram * vector(ZZ, z)) == 16
    for z in stream_p4
)

print(
    f"Q80A1STREAM|selfcheck_p4|pari_ball16_signed={pari_le_16}|"
    f"pari_ball14_signed={pari_le_14}|"
    f"exact_shell_signed={pari_exact_signed}|"
    f"exact_shell_pairs={pari_exact_pairs}|stream_pairs={len(stream_p4)}|"
    "status=PASS_EXACT_STREAMER_VS_PARI",
    flush=True,
)

started = perf_counter()
found = None

# p=4 is already in memory from the self-check and is cheap, so search it
# first.  If it fails, use the same exact streamer on canonical-success p=1.
def search_precomputed_p4():
    p = ZZ(4)
    target_height = QQ(4)
    simple = data["simple"]
    root_inverse = data["root_inverse"]
    projection = data["projection"]
    mw_basis = data["mw_basis"]
    mw_pre = data["mw_pre"]
    compatible = 0

    for index, z in enumerate(stream_p4, 1):
        mw = vector(QQ, z) * mw_basis
        assert QQ(mw * frame_a1 * mw) == target_height
        x0 = z * mw_pre
        p0 = x0 * frame_a1 * simple.transpose()
        assert len(p0) == 1 and p0[0] in ZZ
        delta = (
            (vector(QQ, [p]) - vector(QQ, p0))
            * root_inverse
        )
        if not all(entry in ZZ for entry in delta):
            continue

        root_part = vector(QQ, [p]) * root_inverse * simple
        vq = vector(QQ, mw) + root_part
        if not all(entry in ZZ for entry in vq):
            continue
        v = vector(ZZ, vq)
        assert v * frame_a1 * v == 12
        compatible += 1

        child, _ = neighbor(
            frame_a1, ZZ(6), ZZ(2), ZZ(3), v
        )
        if is_rootless(child):
            assert child.det() == 948
            assert root_signature(child)[:2] == ("rootless", 0)
            CERT.parent.mkdir(parents=True, exist_ok=True)
            CERT.write_text(
                "source=A1/MW16 candidate1\nq=6\na=2\nb=3\n"
                "A1_pairing=4\n"
                "v=" + ",".join(map(str, tuple(v))) + "\n"
                "child=rootless/MW17\ndet=948\n"
            )
            return {
                "found": True, "p": p, "v": tuple(v),
                "shell_pairs_seen": index, "compatible": compatible,
            }

    return {
        "found": False, "p": p,
        "shell_pairs_seen": len(stream_p4), "compatible": compatible,
    }


p4_result = search_precomputed_p4()
print(
    f"Q80A1STREAM|p=4|shell_pairs={p4_result['shell_pairs_seen']}|"
    f"gluing_compatible={p4_result['compatible']}|"
    f"found_rootless={int(p4_result['found'])}|"
    f"elapsed={perf_counter()-started:.1f}s|"
    "status=PASS_EXHAUSTIVE_P4_SHELL",
    flush=True,
)

if p4_result["found"]:
    found = p4_result
else:
    # Canonical successful terminal move has p=1, so test that exact shell next.
    assert canon_pairing == 1
    p1_result = search_exact_pairing_shell(
        frame_a1, data, ZZ(1), CERT, started
    )
    print(
        f"Q80A1STREAM|p=1|shell_pairs={p1_result['shell_pairs_seen']}|"
        f"gluing_compatible={p1_result['compatible']}|"
        f"found_rootless={int(p1_result['found'])}|"
        f"elapsed={perf_counter()-started:.1f}s|"
        "status=PASS_EXACT_P1_SHELL_SEARCH",
        flush=True,
    )
    if p1_result["found"]:
        found = p1_result

if found is not None:
    print(
        f"Q80A1Q6|pairing={found['p']}|q=6|presentation=2,3|"
        f"v={found['v']}|endpoint=rootless/MW17|det=948|"
        f"certificate={CERT}|elapsed={perf_counter()-started:.1f}s|"
        "status=PASS_NEW_A1_Q6_ROOTLESS",
        flush=True,
    )
    print(
        "Q80A1Q6|path="
        "A6+A3/MW8-q6-A4+A2+A1/MW10-"
        "q4-A3+A2/MW12-q4-4A1/MW13-"
        "q4-A1/MW16-q6-rootless/MW17|"
        "status=PASS_COMPLETE_NEW_LOWQ_ROOTLESS_PATH",
        flush=True,
    )
else:
    print(
        f"Q80A1Q6|p4_rootless=0|p1_rootless=0|"
        f"elapsed={perf_counter()-started:.1f}s|"
        "status=PASS_EXHAUSTIVE_P4_P1_NO_ROOTLESS_CANDIDATE1",
        flush=True,
    )

"""


def find_repo() -> Path:
    candidates = [
        Path.cwd(),
        Path.home() / "Documents" / "jacobian-research",
    ]
    for candidate in candidates:
        if (
            candidate
            / "elkies-k3"
            / "data"
            / "fibrations"
        ).is_dir():
            return candidate
    raise SystemExit(
        "Could not locate jacobian-research; "
        "run from repo or keep it in ~/Documents/jacobian-research"
    )


def main():
    repo = find_repo()
    sage = shutil.which("sage") or "/usr/local/bin/sage"
    if shutil.which("sage") is None and not Path(sage).exists():
        raise SystemExit("sage not found")

    cert = Path.home() / "Downloads" / "q80_new_lowq_rootless_candidate1.txt"

    print(f"repo={repo}", flush=True)
    print(f"sage={sage}", flush=True)
    print(f"certificate={cert}", flush=True)

    with tempfile.TemporaryDirectory() as td:
        script = Path(td) / "q80_a1_candidate1_q6_rootless.sage"
        script.write_text(SAGE_CODE)
        subprocess.run(
            [sage, str(script), str(repo), str(cert)],
            check=True,
        )


if __name__ == "__main__":
    main()
