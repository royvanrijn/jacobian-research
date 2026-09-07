#!/usr/bin/env python3
from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path

SAGE_CODE = r"""
from __future__ import print_function

import csv
import json
import sys
from collections import Counter, defaultdict
from itertools import product
from pathlib import Path
from time import perf_counter

from sage.all import (
    QQ, ZZ, QuadraticForm, block_diagonal_matrix, ceil, floor, gcd,
    identity_matrix, lcm, matrix, vector, xgcd
)

ROOT = Path(sys.argv[1]).resolve()
OUT = Path(sys.argv[2]).resolve()
CHECKPOINT = OUT.with_name(OUT.stem + "_checkpoint.jsonl")
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
    assert gcd([abs(ZZ(x)) for x in ns * F]) == 1
    mate = bezout_vector(list(ns * F))
    mate -= ZZ(mate * ns * mate) // 2 * F
    complement = matrix(
        ZZ, [list(F * ns), list(mate * ns)]
    ).right_kernel_matrix()
    child = -(complement * ns * complement.transpose())
    T = matrix(ZZ, [list(F), list(mate)] + complement.rows())
    assert abs(T.det()) == 1
    return child, T


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


def highest_roots(gram, simple, positive):
    components = connected_components(
        simple * gram * simple.transpose()
    )
    inverse_simple = simple.pseudoinverse()
    result = []
    for component in components:
        candidates = []
        for root in positive:
            coordinates = vector(QQ, root) * inverse_simple
            if not all(x in ZZ and x >= 0 for x in coordinates):
                continue
            support = tuple(
                i for i, x in enumerate(coordinates) if x
            )
            if support and all(i in component for i in support):
                candidates.append(
                    (sum(coordinates), root, coordinates)
                )
        _, root, coordinates = max(
            candidates, key=lambda item: item[0]
        )
        result.append(
            (component, vector(ZZ, root), vector(ZZ, coordinates))
        )
    return result


def intersection(left, right, ns):
    return ZZ(vector(ZZ, left) * ns * vector(ZZ, right))


def chamber_reduce(divisor, curves, ns):
    divisor = vector(ZZ, divisor)
    sequence = []
    while True:
        for name, curve in curves:
            pairing = intersection(divisor, curve, ns)
            if pairing < 0:
                divisor += pairing * curve
                sequence.append((name, pairing))
                assert divisor * ns * divisor == 0
                break
        else:
            return divisor, tuple(sequence)


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
        sub = cartan.matrix_from_rows_and_columns(
            component, component
        )
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


def build_chamber_context(frame):
    ns = block_diagonal_matrix(U, -frame)
    simple, positive = deterministic_simple_roots(frame)
    cartan = simple * frame * simple.transpose()
    components = connected_components(cartan)

    fiber = vector(ZZ, [1, 0] + [0] * frame.nrows())
    zero = vector(ZZ, [-1, 1] + [0] * frame.nrows())

    curves = [("O", zero)]
    for index, root in enumerate(simple.rows(), 1):
        curves.append(
            (f"R{index}", vector(ZZ, [0, 0] + list(root)))
        )
    for component_index, (_, root, _) in enumerate(
        highest_roots(frame, simple, positive), 1
    ):
        curves.append(
            (
                f"Theta0_{component_index}",
                fiber - vector(ZZ, [0, 0] + list(root)),
            )
        )

    return {
        "frame": frame,
        "ns": ns,
        "simple": simple,
        "cartan": cartan,
        "components": components,
        "fiber": fiber,
        "zero": zero,
        "curves": tuple(curves),
        "cvp_cache": {},
    }


def exact_component_nearest_shifts(context, component, target):
    subgram = context["cartan"].matrix_from_rows_and_columns(
        component, component
    )
    target = vector(QQ, target)

    integer_part = vector(
        ZZ, [ZZ(floor(value)) for value in target]
    )
    fractional = vector(QQ, target - integer_part)
    key = (tuple(component), tuple(fractional))

    cached = context["cvp_cache"].get(key)
    if cached is None:
        z0 = vector(
            ZZ, [ZZ(floor(-value)) for value in fractional]
        )
        x0 = fractional + vector(QQ, z0)
        best = QQ(x0 * subgram * x0)

        inverse = subgram.inverse()
        ranges = []
        for i in range(subgram.nrows()):
            radius2 = QQ(best * inverse[i, i])
            radius = ZZ(0)
            while QQ(radius * radius) < radius2:
                radius += 1
            lo = ZZ(ceil(-fractional[i] - radius))
            hi = ZZ(floor(-fractional[i] + radius))
            ranges.append(range(int(lo), int(hi) + 1))

        winners = []
        for entries in product(*ranges):
            z = vector(ZZ, entries)
            x = fractional + vector(QQ, z)
            value = QQ(x * subgram * x)
            if value < best:
                best = value
                winners = [tuple(z)]
            elif value == best:
                winners.append(tuple(z))

        cached = (best, tuple(sorted(set(winners))))
        context["cvp_cache"][key] = cached

    _, fractional_winners = cached
    return tuple(
        vector(ZZ, winner) - integer_part
        for winner in fractional_winners
    )


def exact_nearest_root_shifts(context, root_coordinates):
    component_winners = []
    for component in context["components"]:
        target = vector(
            QQ, [root_coordinates[i] for i in component]
        )
        winners = exact_component_nearest_shifts(
            context, component, target
        )
        component_winners.append((component, winners))

    shifts = []
    for choices in product(
        *(item[1] for item in component_winners)
    ):
        full = [ZZ(0)] * context["simple"].nrows()
        for (component, _), choice in zip(
            component_winners, choices
        ):
            for index, value in zip(component, choice):
                full[index] = ZZ(value)
        shifts.append(vector(ZZ, full))
    return tuple(shifts)


def effective_shortest_section(
    context, reduced, root_coordinates, mw_projection
):
    frame = context["frame"]
    ns = context["ns"]
    simple = context["simple"]
    fiber = context["fiber"]
    zero = context["zero"]

    shifts = exact_nearest_root_shifts(
        context, root_coordinates
    )
    accepted = []
    shortest_norm = None

    for root_shift in shifts:
        short_lift = (
            vector(ZZ, reduced[2:]) + root_shift * simple
        )
        short_norm = ZZ(short_lift * frame * short_lift)
        if shortest_norm is None:
            shortest_norm = short_norm
        else:
            assert short_norm == shortest_norm

        assert short_norm % 2 == 0
        section_pole = short_norm // 2 - 2
        section = vector(
            ZZ,
            [section_pole + 1, 1] + list(short_lift),
        )
        assert section * ns * section == -2
        assert intersection(section, fiber, ns) == 1
        assert intersection(section, zero, ns) == section_pole

        pairings = tuple(
            (name, intersection(section, curve, ns))
            for name, curve in context["curves"]
        )
        if all(value >= 0 for _, value in pairings):
            accepted.append(
                (
                    section, short_lift, short_norm,
                    section_pole, root_shift, pairings,
                )
            )

    if len(accepted) != 1:
        raise ArithmeticError(
            "Expected exactly one effective shortest section; "
            f"found {len(accepted)} for MW={tuple(mw_projection)}"
        )
    return accepted[0], len(shifts)


def score_geometry(context, a, b, v):
    frame = context["frame"]
    ns = context["ns"]
    simple = context["simple"]
    cartan = context["cartan"]
    fiber = context["fiber"]
    zero = context["zero"]

    raw = vector(ZZ, [a, b] + list(v))
    assert raw * ns * raw == 0
    reduced, sequence = chamber_reduce(
        raw, context["curves"], ns
    )

    degree = intersection(reduced, fiber, ns)
    d_o = intersection(reduced, zero, ns)

    frame_part = vector(QQ, reduced[2:])
    root_coordinates = (
        frame_part
        * frame
        * simple.transpose()
        * cartan.inverse()
    )
    root_projection = root_coordinates * simple
    mw_projection = frame_part - root_projection
    mw_height = QQ(
        mw_projection * frame * mw_projection
    )

    selected, nearest_count = effective_shortest_section(
        context, reduced, root_coordinates, mw_projection
    )
    (
        section, short_lift, short_norm, section_pole,
        root_shift, section_pairings,
    ) = selected

    residual = reduced - (degree - 1) * zero - section
    assert intersection(residual, fiber, ns) == 0
    assert residual[1] == 0

    root_coefficients = (
        vector(QQ, residual[2:]) * simple.pseudoinverse()
    )
    assert (
        root_coefficients * simple
        == vector(QQ, residual[2:])
    )
    vertical_integral = all(
        value in ZZ for value in root_coefficients
    )

    support_fibres = 0
    support_components = 0
    for component in context["components"]:
        nz = sum(
            1 for index in component
            if root_coefficients[index] != 0
        )
        support_components += nz
        if nz:
            support_fibres += 1

    vertical_l1 = QQ(
        sum(abs(value) for value in root_coefficients)
    )
    vertical_max = QQ(max(
        [abs(value) for value in root_coefficients]
        + [QQ(0)]
    ))
    fiber_twist = ZZ(residual[0])

    return {
        "reduced": tuple(reduced),
        "degree": int(degree),
        "d_o": int(d_o),
        "mw_height": mw_height,
        "short_norm": int(short_norm),
        "section_pole": int(section_pole),
        "fiber_twist": int(fiber_twist),
        "vertical_integral": bool(vertical_integral),
        "support_fibres": int(support_fibres),
        "support_components": int(support_components),
        "vertical_l1": vertical_l1,
        "vertical_max": vertical_max,
        "reflection_count": len(sequence),
        "root_coefficients": tuple(root_coefficients),
        "nearest_section_lifts": int(nearest_count),
        "section_root_shift": tuple(root_shift),
        "section_nonzero_pairings": tuple(
            (name, value)
            for name, value in section_pairings if value
        ),
    }


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


def quotient_data(frame, max_norm):
    simple, _ = deterministic_simple_roots(frame)
    root_gram = simple * frame * simple.transpose()
    components = connected_components(root_gram)

    projection = (
        identity_matrix(QQ, 17)
        - frame
        * simple.transpose()
        * root_gram.inverse()
        * simple
    )
    den = lcm(
        value.denominator() for value in projection.list()
    )
    scaled_projection = (den * projection).change_ring(ZZ)

    projected_integer = (
        scaled_projection.row_module().basis_matrix()
    )
    mw_rank = 17 - simple.nrows()
    assert projected_integer.nrows() == mw_rank
    mw_basis = projected_integer / den
    mw_height = mw_basis * frame * mw_basis.transpose()

    hden = lcm(
        value.denominator() for value in mw_height.list()
    )
    scaled_height = (
        2 * hden * mw_height
    ).change_ring(ZZ)
    transform = scaled_height.LLL_gram().transpose()
    assert abs(transform.det()) == 1
    mw_basis = transform * mw_basis
    mw_height = mw_basis * frame * mw_basis.transpose()

    hden = lcm(
        value.denominator() for value in mw_height.list()
    )
    mw_form = QuadraticForm(
        ZZ, (2 * hden * mw_height).change_ring(ZZ)
    )
    shells = mw_form.short_vector_list_up_to_length(
        int(max_norm * hden + 1),
        up_to_sign_flag=True,
    )

    coords = {tuple([0] * mw_rank)}
    for shell in shells[1:]:
        for row in shell:
            coords.add(tuple(vector(ZZ, row)))
    coords = tuple(sorted(
        coords,
        key=lambda z: (
            QQ(vector(ZZ, z) * mw_height * vector(ZZ, z)),
            z,
        ),
    ))

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
        assert (
            vector(QQ, pre) * projection
            == vector(QQ, projected)
        )
        return vector(ZZ, pre)

    mw_preimage_basis = matrix(
        ZZ,
        [list(integral_preimage(row))
         for row in mw_basis.rows()],
    )
    assert mw_preimage_basis * projection == mw_basis

    component_weights = []
    for component in components:
        sub = root_gram.matrix_from_rows_and_columns(
            component, component
        )
        component_weights.append(
            dominant_weights(sub, max_norm)
        )

    root_weights_by_norm = defaultdict(list)

    def combine(component_index, entries, total_norm):
        if component_index == len(components):
            p = [ZZ(0)] * simple.nrows()
            for component, values in zip(
                components, entries
            ):
                for i, value in zip(component, values):
                    p[i] = ZZ(value)
            root_weights_by_norm[
                QQ(total_norm)
            ].append(vector(ZZ, p))
            return

        for values, norm in component_weights[component_index]:
            new_norm = total_norm + norm
            if new_norm <= max_norm:
                combine(
                    component_index + 1,
                    entries + [values],
                    new_norm,
                )

    combine(0, [], QQ(0))

    return {
        "simple": simple,
        "root_gram": root_gram,
        "root_inverse": root_gram.inverse(),
        "components": components,
        "projection": projection,
        "mw_basis": mw_basis,
        "mw_height": mw_height,
        "mw_coords": coords,
        "mw_preimage_basis": mw_preimage_basis,
        "root_weights_by_norm": root_weights_by_norm,
        "component_counts": tuple(
            len(x) for x in component_weights
        ),
    }


def dominant_shell_vectors(frame, qnorm):
    N = ZZ(2 * qnorm)
    data = quotient_data(frame, N)

    simple = data["simple"]
    projection = data["projection"]
    mw_basis = data["mw_basis"]
    mw_height = data["mw_height"]
    mw_pre = data["mw_preimage_basis"]
    root_inverse = data["root_inverse"]
    roots_by_norm = data["root_weights_by_norm"]

    seen = set()
    vectors = []

    for z_tuple in data["mw_coords"]:
        z = vector(ZZ, z_tuple)
        mw = z * mw_basis
        mw_norm = QQ(mw * frame * mw)
        needed = QQ(N) - mw_norm
        if needed < 0:
            continue

        x0 = z * mw_pre
        assert (
            vector(QQ, x0) * projection
            == vector(QQ, mw)
        )
        p0 = x0 * frame * simple.transpose()
        assert all(x in ZZ for x in p0)
        p0 = vector(ZZ, p0)

        for p in roots_by_norm.get(needed, ()):
            delta = (
                (vector(QQ, p) - vector(QQ, p0))
                * root_inverse
            )
            if not all(x in ZZ for x in delta):
                continue

            root_part = (
                vector(QQ, p) * root_inverse * simple
            )
            vq = vector(QQ, mw) + root_part
            assert all(x in ZZ for x in vq)
            v = vector(ZZ, vq)
            assert v * frame * v == N

            key = tuple(v)
            assert key not in seen
            seen.add(key)
            vectors.append(v)

    return tuple(vectors), data


def checkpoint_payload(row):
    return {
        key: (
            str(value)
            if key in {
                "mw_height", "vertical_l1", "vertical_max"
            }
            else list(value)
            if key in {
                "v", "root_coefficients",
                "section_root_shift", "reduced"
            }
            else [
                [str(name), int(pairing)]
                for name, pairing in value
            ]
            if key == "section_nonzero_pairings"
            else value
        )
        for key, value in row.items()
    }


# Reconstruct retained Q80 corridor through A6+A3/MW8.
with (DATA / "kumar_q80_to_rootless_path.tsv").open() as handle:
    steps = list(csv.DictReader(handle, delimiter="\t"))

start = load_matrix(
    DATA / "kumar_q80_e6_d5_a3_mw3_frame.txt"
)
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

assert root_signature(a6a3_frame)[:2] == ("A6+A3", 9)

v7774 = vector(
    ZZ,
    (85,2699,1257,7718,3756,-41,3077,-4614,-6615,6032,2584,-1678,121,-736,-913,1,1165),
)
assert v7774 * a6a3_frame * v7774 == 12
frame7774, _ = neighbor(
    a6a3_frame, ZZ(6), ZZ(2), ZZ(3), v7774
)

source_sig = root_signature(frame7774)
assert source_sig[:2] == ("A4+A2+A1", 7), source_sig
assert frame7774.det() == 948

print(
    f"Q80R5SCORE|source=7774|ADE={source_sig[0]}|"
    f"root_rank={source_sig[1]}|MW={17-source_sig[1]}|"
    f"det={frame7774.det()}|status=PASS_SOURCE",
    flush=True,
)

# Enumerate the exact q4 Weyl/sign shell once.
started = perf_counter()
vectors, qdata = dominant_shell_vectors(frame7774, ZZ(4))
assert len(vectors) == 7815, len(vectors)

print(
    f"Q80R5SCORE|q4_shell={len(vectors)}|"
    f"mw_vectors={len(qdata['mw_coords'])}|"
    f"dominant_component_counts={qdata['component_counts']}|"
    "status=PASS_Q4_SHELL_RECONSTRUCTION",
    flush=True,
)

context = build_chamber_context(frame7774)
ns7774 = block_diagonal_matrix(U, -frame7774)

CHECKPOINT.parent.mkdir(parents=True, exist_ok=True)
CHECKPOINT.write_text("")

rank_counter = Counter()
ade5_counter = Counter()
divisibility_counter = Counter()
rows = []
eligible_id = 0

for shell_id, v in enumerate(vectors, 1):
    F = vector(ZZ, [2, 2] + list(v))
    assert F * ns7774 * F == 0
    divisibility = gcd(
        [abs(ZZ(x)) for x in ns7774 * F]
    )
    divisibility_counter[int(divisibility)] += 1
    if divisibility != 1:
        continue

    eligible_id += 1
    child, _ = neighbor(
        frame7774, ZZ(4), ZZ(2), ZZ(2), v
    )
    ade, root_rank, root_count, root_det = root_signature(child)
    rank_counter[root_rank] += 1

    if root_rank != 5:
        continue

    ade5_counter[ade] += 1
    geometry = score_geometry(context, 2, 2, v)

    row = {
        "shell_id": int(shell_id),
        "eligible_id": int(eligible_id),
        "ade": ade,
        "root_rank": int(root_rank),
        "mw_rank": int(17 - root_rank),
        "root_count": int(root_count),
        "root_det": int(root_det),
        "d_f": geometry["degree"],
        "d_o": geometry["d_o"],
        "mw_height": geometry["mw_height"],
        "short_norm": geometry["short_norm"],
        "section_p_o": geometry["section_pole"],
        "fiber_twist": geometry["fiber_twist"],
        "vertical_integral": int(
            geometry["vertical_integral"]
        ),
        "support_fibres": geometry["support_fibres"],
        "support_components": geometry[
            "support_components"
        ],
        "vertical_l1": geometry["vertical_l1"],
        "vertical_max": geometry["vertical_max"],
        "reflection_count": geometry[
            "reflection_count"
        ],
        "nearest_section_lifts": geometry[
            "nearest_section_lifts"
        ],
        "section_root_shift": geometry[
            "section_root_shift"
        ],
        "section_nonzero_pairings": geometry[
            "section_nonzero_pairings"
        ],
        "root_coefficients": geometry[
            "root_coefficients"
        ],
        "reduced": geometry["reduced"],
        "v": tuple(v),
    }
    rows.append(row)

    with CHECKPOINT.open("a") as handle:
        handle.write(
            json.dumps(
                checkpoint_payload(row),
                sort_keys=True,
            )
            + "\n"
        )
        handle.flush()

# Exact regressions from the already completed low-q scan.
assert eligible_id == 7812, eligible_id
assert dict(divisibility_counter) == {1: 7812, 2: 3}
assert dict(rank_counter) == {
    5: 118,
    6: 894,
    7: 2348,
    8: 2437,
    9: 1551,
    10: 429,
    11: 35,
}, dict(rank_counter)
assert len(rows) == 118, len(rows)

expected_ade5 = {
    "A2+A1+A1+A1": 53,
    "A2+A2+A1": 2,
    "A3+A1+A1": 57,
    "A3+A2": 4,
    "A4+A1": 2,
}
assert dict(ade5_counter) == expected_ade5, dict(ade5_counter)


def score_key(row):
    return (
        row["d_f"],
        row["section_p_o"],
        0 if row["vertical_integral"] else 1,
        row["support_fibres"],
        row["support_components"],
        row["vertical_l1"],
        row["vertical_max"],
        row["mw_height"],
        row["short_norm"],
        row["reflection_count"],
        len(row["ade"].split("+")),
        row["ade"],
        row["shell_id"],
    )


rows.sort(key=score_key)

degree_counter = Counter(
    row["d_f"] for row in rows
)
degree_pole_counter = Counter(
    (row["d_f"], row["section_p_o"])
    for row in rows
)
cheap_gate = [
    row for row in rows
    if row["d_f"] == 2
    and row["section_p_o"] <= 1
    and row["vertical_integral"]
]
moderate_gate = [
    row for row in rows
    if row["d_f"] == 2
    and row["section_p_o"] <= 3
    and row["vertical_integral"]
]

fields = [
    "rank", "shell_id", "eligible_id", "ade",
    "root_rank", "mw_rank", "root_count", "root_det",
    "d_f", "d_o", "mw_height", "short_norm",
    "section_p_o", "fiber_twist", "vertical_integral",
    "support_fibres", "support_components",
    "vertical_l1", "vertical_max", "reflection_count",
    "nearest_section_lifts", "section_root_shift",
    "section_nonzero_pairings", "root_coefficients",
    "reduced", "v",
]

# Persist before any human-readable ranking output.
OUT.parent.mkdir(parents=True, exist_ok=True)
with OUT.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle, fieldnames=fields, delimiter="\t"
    )
    writer.writeheader()
    for rank, row in enumerate(rows, 1):
        record = dict(row)
        record["rank"] = rank
        for key in (
            "mw_height", "vertical_l1", "vertical_max"
        ):
            record[key] = str(record[key])
        for key in (
            "section_root_shift", "root_coefficients",
            "reduced", "v",
        ):
            record[key] = ",".join(
                map(str, record[key])
            )
        record["section_nonzero_pairings"] = ";".join(
            f"{name}:{value}"
            for name, value
            in record["section_nonzero_pairings"]
        )
        writer.writerow(record)

checkpoint_lines = [
    line for line in CHECKPOINT.read_text().splitlines()
    if line.strip()
]
assert len(checkpoint_lines) == 118

with OUT.open() as handle:
    persisted = list(
        csv.DictReader(handle, delimiter="\t")
    )
assert len(persisted) == 118

print(
    f"Q80R5SCORE|persisted_rank5=118|tsv={OUT}|"
    f"checkpoint={CHECKPOINT}|"
    "status=PASS_DURABLE_RESULT",
    flush=True,
)
print(
    f"Q80R5SCORE|rank5_ADE="
    f"{tuple(sorted(ade5_counter.items()))}",
    flush=True,
)
print(
    f"Q80R5SCORE|degree_distribution="
    f"{tuple(sorted(degree_counter.items()))}",
    flush=True,
)
print(
    f"Q80R5SCORE|degree_pole_distribution="
    f"{tuple(sorted(degree_pole_counter.items()))}",
    flush=True,
)
print(
    f"Q80R5SCORE|cheap_gate={len(cheap_gate)}|"
    "criterion=D.F=2,P.O<=1,integral",
    flush=True,
)
print(
    f"Q80R5SCORE|moderate_gate={len(moderate_gate)}|"
    "criterion=D.F=2,P.O<=3,integral",
    flush=True,
)

for rank, row in enumerate(rows[:20], 1):
    print(
        f"Q80R5TOP|rank={rank}|shell_id={row['shell_id']}|"
        f"ADE={row['ade']}|D.F={row['d_f']}|"
        f"D.O={row['d_o']}|P.O={row['section_p_o']}|"
        f"MWheight={row['mw_height']}|"
        f"shortnorm={row['short_norm']}|"
        f"twist={row['fiber_twist']}|"
        f"vfibres={row['support_fibres']}|"
        f"vcomponents={row['support_components']}|"
        f"L1={row['vertical_l1']}|"
        f"max={row['vertical_max']}|"
        f"nearest={row['nearest_section_lifts']}|"
        f"refl={row['reflection_count']}|"
        f"v={row['v']}",
        flush=True,
    )

best = rows[0]
print(
    f"Q80R5BEST|shell_id={best['shell_id']}|"
    f"ADE={best['ade']}|D.F={best['d_f']}|"
    f"D.O={best['d_o']}|P.O={best['section_p_o']}|"
    f"MWheight={best['mw_height']}|"
    f"shortnorm={best['short_norm']}|"
    f"twist={best['fiber_twist']}|"
    f"vfibres={best['support_fibres']}|"
    f"vcomponents={best['support_components']}|"
    f"L1={best['vertical_l1']}|"
    f"max={best['vertical_max']}|"
    f"nearest={best['nearest_section_lifts']}|"
    f"refl={best['reflection_count']}|"
    f"rootcoeff={best['root_coefficients']}|"
    f"v={best['v']}",
    flush=True,
)
print(
    f"Q80R5SCORE|elapsed={perf_counter()-started:.1f}s|"
    "status=PASS_7774_Q4_RANK5_GEOMETRY_SCORE",
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
    if (
        shutil.which("sage") is None
        and not Path(sage).exists()
    ):
        raise SystemExit("sage not found")

    out = (
        Path.home()
        / "Downloads"
        / "q80_7774_q4_rank5_scores.tsv"
    )

    print(f"repo={repo}", flush=True)
    print(f"sage={sage}", flush=True)
    print(f"out={out}", flush=True)

    with tempfile.TemporaryDirectory() as td:
        script = Path(td) / "q80_7774_q4_rank5_score.sage"
        script.write_text(SAGE_CODE)
        subprocess.run(
            [sage, str(script), str(repo), str(out)],
            check=True,
        )


if __name__ == "__main__":
    main()
