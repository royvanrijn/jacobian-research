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



# ----------------------------------------------------------------------
# Complete retained new low-q path, now including the newly found rootless q6.
# ----------------------------------------------------------------------
with (DATA / "kumar_q80_to_rootless_path.tsv").open() as handle:
    canonical_steps = list(csv.DictReader(handle, delimiter="\t"))

start = load_matrix(DATA / "kumar_q80_e6_d5_a3_mw3_frame.txt")

# Common certified q4,q4 prefix to D7+D5/MW5.
first = canonical_steps[0]
frame1, _ = neighbor(
    start,
    ZZ(first["q"]), ZZ(first["a"]), ZZ(first["b"]),
    vector(ZZ, map(ZZ, first["v"].split(","))),
)
second = canonical_steps[1]
frame2, _ = neighbor(
    frame1,
    ZZ(second["q"]), ZZ(second["a"]), ZZ(second["b"]),
    vector(ZZ, map(ZZ, second["v"].split(","))),
)

assert root_signature(frame2)[:2] == ("D7+D5", 12)
assert frame2.det() == 948

steps = [
    {
        "name": "escape",
        "q": 6, "a": 2, "b": 3,
        "v": (-5,-3,6,6,-8,-4,2,4,-1,8,-16,-1,0,3,5,-2,-2),
        "child_ade": "D7+D4", "child_rank": 11,
    },
    {
        "name": "orbit424",
        "q": 4, "a": 2, "b": 2,
        "v": (32,48,-21,28,8,-52,-34,0,18,5,-23,43,9,-18,16,-6,-6),
        "child_ade": "A6+A4", "child_rank": 10,
    },
    {
        "name": "orbit1222",
        "q": 4, "a": 2, "b": 2,
        "v": (10,53,-192,-114,29,-256,-170,-12,-14,74,-32,-14,-6,-26,-58,84,-28),
        "child_ade": "A6+A3", "child_rank": 9,
    },
    {
        "name": "q6_7774",
        "q": 6, "a": 2, "b": 3,
        "v": (85,2699,1257,7718,3756,-41,3077,-4614,-6615,6032,2584,-1678,121,-736,-913,1,1165),
        "child_ade": "A4+A2+A1", "child_rank": 7,
    },
    {
        "name": "q4_1938",
        "q": 4, "a": 2, "b": 2,
        "v": (-94,-1003,5298,4977,-1431,-1440,100,1,-1632,1893,1634,-1264,-4175,2248,-3111,1561,2842),
        "child_ade": "A3+A2", "child_rank": 5,
    },
    {
        "name": "q4_6855",
        "q": 4, "a": 2, "b": 2,
        "v": (30693,-339,-2534,45446,10413,16390,-11527,5970,-18424,4193,21146,11296,25035,17925,-6032,4304,7717),
        "child_ade": "A1+A1+A1+A1", "child_rank": 4,
    },
    {
        "name": "q4_a1_candidate1",
        "q": 4, "a": 2, "b": 2,
        "v": (21,671,-20182,-10366,27727,30558,5582,20831,-10195,-19691,6086,10389,20928,18651,16123,15473,-11496),
        "child_ade": "A1", "child_rank": 1,
    },
    {
        "name": "q6_rootless",
        "q": 6, "a": 2, "b": 3,
        "v": (-44717,-282065,63356,564493,-98198,249323,239104,-1054,-22328,-389456,-231271,-641746,-570362,-123785,227276,-186445,89497),
        "child_ade": "rootless", "child_rank": 0,
    },
]

rows = []
parent = frame2
parent_ade, parent_rank, _, _ = root_signature(parent)

for index, step in enumerate(steps, 1):
    v = vector(ZZ, step["v"])
    q = ZZ(step["q"])
    a = ZZ(step["a"])
    b = ZZ(step["b"])

    assert v * parent * v == 2*q

    context = build_chamber_context(parent)
    geometry = score_geometry(context, a, b, v)

    child, _ = neighbor(parent, q, a, b, v)
    child_ade, child_rank, child_roots, child_rootdet = root_signature(child)

    assert child.det() == 948
    assert child_ade == step["child_ade"], (step["name"], child_ade)
    assert child_rank == step["child_rank"], (step["name"], child_rank)

    row = {
        "step": index,
        "name": step["name"],
        "source_ade": parent_ade,
        "source_root_rank": parent_rank,
        "source_mw": 17-parent_rank,
        "q": int(q),
        "a": int(a),
        "b": int(b),
        "child_ade": child_ade,
        "child_root_rank": child_rank,
        "child_mw": 17-child_rank,
        "child_roots": int(child_roots),
        "child_rootdet": int(child_rootdet),
        "d_f": int(geometry["degree"]),
        "d_o": int(geometry["d_o"]),
        "p_o": int(geometry["section_pole"]),
        "mw_height": geometry["mw_height"],
        "short_norm": int(geometry["short_norm"]),
        "fiber_twist": int(geometry["fiber_twist"]),
        "vertical_integral": int(geometry["vertical_integral"]),
        "support_fibres": int(geometry["support_fibres"]),
        "support_components": int(geometry["support_components"]),
        "vertical_l1": geometry["vertical_l1"],
        "vertical_max": geometry["vertical_max"],
        "nearest_section_lifts": int(geometry["nearest_section_lifts"]),
        "reflection_count": int(geometry["reflection_count"]),
        "root_coefficients": tuple(geometry["root_coefficients"]),
        "v": tuple(v),
    }
    rows.append(row)

    print(
        f"Q80NEWPATH|step={index}|name={step['name']}|"
        f"{parent_ade}/MW{17-parent_rank}--q{q}({a},{b})-->"
        f"{child_ade}/MW{17-child_rank}|"
        f"D.F={row['d_f']}|D.O={row['d_o']}|P.O={row['p_o']}|"
        f"MWheight={row['mw_height']}|shortnorm={row['short_norm']}|"
        f"twist={row['fiber_twist']}|"
        f"vfibres={row['support_fibres']}|"
        f"vcomponents={row['support_components']}|"
        f"L1={row['vertical_l1']}|max={row['vertical_max']}|"
        f"nearest={row['nearest_section_lifts']}|"
        f"vertical_integral={row['vertical_integral']}|"
        "status=PASS_STEP",
        flush=True,
    )

    parent = child
    parent_ade = child_ade
    parent_rank = child_rank

assert parent_ade == "rootless"
assert parent_rank == 0
assert parent.det() == 948

# Existing corrected regressions for the two earlier q4 moves.
by_name = {row["name"]: row for row in rows}
assert by_name["orbit424"]["d_f"] == 2
assert by_name["orbit424"]["p_o"] == 1
assert by_name["orbit424"]["support_fibres"] == 1
assert by_name["orbit424"]["support_components"] == 3
assert by_name["orbit424"]["vertical_l1"] == 3

assert by_name["orbit1222"]["d_f"] == 2
assert by_name["orbit1222"]["p_o"] == 1
assert by_name["orbit1222"]["support_fibres"] == 1
assert by_name["orbit1222"]["support_components"] == 2
assert by_name["orbit1222"]["vertical_l1"] == 2

# Persist the complete exact path/geometry table before summary reporting.
OUT = Path(sys.argv[2]).resolve()
OUT.parent.mkdir(parents=True, exist_ok=True)
fields = [
    "step", "name", "source_ade", "source_root_rank", "source_mw",
    "q", "a", "b", "child_ade", "child_root_rank", "child_mw",
    "child_roots", "child_rootdet",
    "d_f", "d_o", "p_o", "mw_height", "short_norm",
    "fiber_twist", "vertical_integral",
    "support_fibres", "support_components",
    "vertical_l1", "vertical_max", "nearest_section_lifts",
    "reflection_count", "root_coefficients", "v",
]
with OUT.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    for row in rows:
        record = dict(row)
        record["mw_height"] = str(record["mw_height"])
        record["vertical_l1"] = str(record["vertical_l1"])
        record["vertical_max"] = str(record["vertical_max"])
        record["root_coefficients"] = ",".join(
            map(str, record["root_coefficients"])
        )
        record["v"] = ",".join(map(str, record["v"]))
        writer.writerow(record)

with OUT.open() as handle:
    persisted = list(csv.DictReader(handle, delimiter="\t"))
assert len(persisted) == len(steps) == 8

degree_distribution = {}
pole_distribution = {}
for row in rows:
    degree_distribution[row["d_f"]] = degree_distribution.get(row["d_f"], 0) + 1
    key = (row["d_f"], row["p_o"])
    pole_distribution[key] = pole_distribution.get(key, 0) + 1

all_degree_two = all(row["d_f"] == 2 for row in rows)
cheap = [
    row["name"] for row in rows
    if row["d_f"] == 2 and row["p_o"] <= 1 and row["vertical_integral"]
]
moderate = [
    row["name"] for row in rows
    if row["d_f"] == 2 and row["p_o"] <= 2 and row["vertical_integral"]
]

print(
    f"Q80NEWPATH|persisted={OUT}|steps={len(rows)}|"
    "status=PASS_DURABLE_COMPLETE_PATH",
    flush=True,
)
print(
    f"Q80NEWPATH|degree_distribution={tuple(sorted(degree_distribution.items()))}|"
    f"degree_pole_distribution={tuple(sorted(pole_distribution.items()))}|"
    f"all_degree_two={int(all_degree_two)}",
    flush=True,
)
print(
    f"Q80NEWPATH|cheap_P.O_le_1={tuple(cheap)}|"
    f"moderate_P.O_le_2={tuple(moderate)}",
    flush=True,
)
print(
    "Q80NEWPATH|path="
    "D7+D5/MW5-q6-D7+D4/MW6-"
    "q4-A6+A4/MW7-q4-A6+A3/MW8-"
    "q6-A4+A2+A1/MW10-q4-A3+A2/MW12-"
    "q4-4A1/MW13-q4-A1/MW16-q6-rootless/MW17|"
    "status=PASS_COMPLETE_NEW_Q80_ROOTLESS_PATH_WITH_GEOMETRY",
    flush=True,
)
"""


def find_repo() -> Path:
    candidates = [
        Path.cwd(),
        Path.home() / "Documents" / "jacobian-research",
    ]
    for candidate in candidates:
        if (candidate / "elkies-k3" / "data" / "fibrations").is_dir():
            return candidate
    raise SystemExit(
        "Could not locate jacobian-research; run from repo or keep it in "
        "~/Documents/jacobian-research"
    )


def main():
    repo = find_repo()
    sage = shutil.which("sage") or "/usr/local/bin/sage"
    if shutil.which("sage") is None and not Path(sage).exists():
        raise SystemExit("sage not found")

    out = Path.home() / "Downloads" / "q80_new_lowq_rootless_geometry.tsv"
    print(f"repo={repo}", flush=True)
    print(f"sage={sage}", flush=True)
    print(f"out={out}", flush=True)

    with tempfile.TemporaryDirectory() as td:
        script = Path(td) / "q80_new_lowq_rootless_verify_geometry.sage"
        script.write_text(SAGE_CODE)
        subprocess.run(
            [sage, str(script), str(repo), str(out)],
            check=True,
        )


if __name__ == "__main__":
    main()
