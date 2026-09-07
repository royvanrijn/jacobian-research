#!/usr/bin/env sage -python
"""
Exhaustively search the ORIGINAL H3 MW rank-2 section lattice for a cheaper
q8/D13 bridge to q24.

Why this search matters
-----------------------
The raw-q6 search proved the best q6-section bridge below degree 206 has
q8 degree 46.  But the original H3 fibration already has two exact rational
sections P1 and P2.  P2 alone has large q8 degree, yet another MW combination
may have a much smaller *effective* section representative after its E7 root
correction.

This script keeps the original H3 presentation throughout:
    F0 = source H3 fibre,
    O0 = source zero,
    roots = source E7+E8 simple roots,
    basis = (-P1, P2) as actual source section classes.

For every MW word n=(a,b), the unique integral effective section is recovered
from its Shioda vector plus one of the two possible III* multiplicity-one
component patterns (affine or the unique non-affine multiplicity-one E7
component).  II* contributes only its affine component.

For each component pattern, q8 degree is an exact positive-definite quadratic
    d8(n,p) = n^T H0 n + ell*n + c_p.
Completing the square yields a finite box containing EVERY source section with
d8 < --threshold.  The search is therefore exhaustive below the threshold.

We test for anchored D13 coordinates z satisfying
    z[1] = -1, z[3] = 1,
because then q24=(2,-1,-1,1) differs from z only by the already explicit
G1=(1,0,0,0) and G3=(0,0,1,0).

Run:
  sage -python ~/Downloads/search_h92_source_mw_for_q24_bridge.sage
"""

import argparse
import json
from math import isqrt
from pathlib import Path

from sage.all import (
    QQ, ZZ, block_diagonal_matrix, identity_matrix, lcm, matrix, pari, vector
)


E7_CARTAN = matrix(ZZ, [
    [2, 0, 0, -1, 0, 0, 0],
    [0, 2, 0, 0, -1, 0, -1],
    [0, 0, 2, -1, 0, -1, -1],
    [-1, 0, -1, 2, 0, 0, 0],
    [0, -1, 0, 0, 2, 0, 0],
    [0, 0, -1, 0, 0, 2, 0],
    [0, -1, -1, 0, 0, 0, 2],
])

# Source simple-root order is the canonical frame's E7_1..E7_7,
# E8_1..E8_8.  This E8 Cartan is read directly from the frame below rather
# than pinned separately.


def locate_repo(explicit=None):
    candidates = []
    if explicit:
        candidates.append(Path(explicit).expanduser())
    cwd = Path.cwd().resolve()
    candidates += [cwd, *cwd.parents]
    home = Path.home()
    candidates += [
        home / "Documents" / "jacobian-research",
        home / "jacobian-research",
        home / "src" / "jacobian-research",
        home / "git" / "jacobian-research",
        home / "projects" / "jacobian-research",
    ]
    seen = set()
    for candidate in candidates:
        try:
            candidate = candidate.resolve()
        except Exception:
            continue
        if candidate in seen:
            continue
        seen.add(candidate)
        if (candidate / "elkies-k3/scripts").is_dir():
            return candidate
    raise SystemExit("Could not locate jacobian-research; pass --repo PATH")


def load_gram(path):
    return matrix(ZZ, [
        [ZZ(v) for v in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ])


def child_frame_with_zero(ns, fibre, zero):
    assert fibre * ns * fibre == 0
    assert zero * ns * zero == -2
    assert zero * ns * fibre == 1
    mate = zero + fibre
    orth = matrix(
        ZZ, [list(fibre * ns), list(mate * ns)]
    ).right_kernel_matrix()
    basis = matrix(
        ZZ, [list(fibre), list(mate)] + [list(row) for row in orth.rows()]
    )
    assert abs(basis.det()) == 1
    child = -(orth * ns * orth.transpose())
    U2 = matrix(ZZ, ((0, 1), (1, 0)))
    assert basis * ns * basis.transpose() == block_diagonal_matrix(U2, -child)
    return child, basis


def roots_and_data(gram):
    result = pari(gram).qfminim(2)
    count = ZZ(result[0])
    if not count:
        return (), matrix(ZZ, 0, gram.nrows()), (0, 0, 1)
    half = [vector(ZZ, c) for c in matrix(ZZ, result[2]).columns()]
    roots = tuple(half + [-r for r in half])
    rb = matrix(ZZ, [list(r) for r in roots]).row_module().basis_matrix()
    rg = rb * gram * rb.transpose()
    return roots, rb, (rb.rank(), count, abs(ZZ(rg.det())))


def deterministic_simple_roots(gram):
    roots, unused, data = roots_and_data(gram)
    rank = data[0]
    regular = None
    for shift in range(1, 1000):
        candidate = vector(ZZ, [
            (i + 1)**2 + shift*(i + 1) + 1
            for i in range(gram.nrows())
        ])
        if all(candidate * root != 0 for root in roots):
            regular = candidate
            break
    assert regular is not None
    positive = [r for r in roots if regular * r > 0]
    pset = {tuple(r) for r in positive}
    simple = [
        r for r in positive
        if not any(tuple(r-left) in pset for left in positive)
    ]
    simple = matrix(ZZ, [list(r) for r in simple])
    assert simple.nrows() == simple.rank() == rank
    return simple, simple * gram * simple.transpose()


def d13_root_adaptation(child):
    unused, root_basis, invariants = roots_and_data(child)
    assert invariants == (13, 312, 4), invariants
    simple, cartan = deterministic_simple_roots(child)
    assert cartan.det() == 4

    smith, left, right = root_basis.smith_form()
    assert smith == left * root_basis * right
    assert tuple(abs(smith[i, i]) for i in range(13)) == (1,) * 13
    completion = right.inverse()
    initial = simple.stack(completion[13:])
    assert abs(initial.det()) == 1

    adapted = initial * child * initial.transpose()
    root = adapted[:13, :13]
    coupling = adapted[:13, 13:]
    tail = adapted[13:, 13:]
    H = tail - coupling.transpose() * root.inverse() * coupling
    scale = ZZ(1)
    for value in H.list():
        scale = lcm(scale, ZZ(QQ(value).denominator()))
    lll = matrix(ZZ, pari((scale*H).change_ring(ZZ)).qflllgram())
    assert abs(lll.det()) == 1

    change = block_diagonal_matrix(identity_matrix(ZZ, 13), lll.transpose())
    basis = change * initial
    adapted = basis * child * basis.transpose()
    root = adapted[:13, :13]
    coupling = adapted[:13, 13:]
    tail = adapted[13:, 13:]
    H = tail - coupling.transpose() * root.inverse() * coupling
    assert H.det() == 237
    return basis, adapted, H


def ceil_sqrt_rational(value):
    value = QQ(value)
    assert value >= 0
    num = ZZ(value.numerator())
    den = ZZ(value.denominator())
    q = (num + den - 1)//den
    k = ZZ(isqrt(int(q)))
    if k*k < q:
        k += 1
    assert QQ(k*k) >= value
    return k


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo", type=Path)
parser.add_argument(
    "--threshold", type=int, default=53,
    help="exhaustively enumerate every source H3 section with q8 degree below this",
)
parser.add_argument("--top", type=int, default=30)
parser.add_argument("--output", type=Path)
args = parser.parse_args()

ROOT = locate_repo(args.repo)
LOCAL = ROOT / "artifacts/local/elkies-k3"
FRAME = ROOT / "elkies-k3/data/fibrations/kumar_e7e8_mw2_frame_3.txt"
TARGET = LOCAL / "q8-target-component-nef.json"
BRANCH = LOCAL / "q8-d13-branch-anchor.json"
G3FILE = LOCAL / "q8-d13-g3-from-e77-bisection.json"
OUTPUT = (
    args.output.resolve()
    if args.output
    else LOCAL / "source-mw-q24-bridge-search.json"
)

for path in (FRAME, TARGET, BRANCH, G3FILE):
    if not path.exists():
        raise SystemExit(f"Missing prerequisite: {path}")

frame = load_gram(FRAME)
target = json.loads(TARGET.read_text())
branch = json.loads(BRANCH.read_text())
g3data = json.loads(G3FILE.read_text())
assert target["status"] == "PASS_EXACT_Q6_CHILD_Q8_PHYSICAL_ROOT_TARGET"
assert target["normalization"]["representative"] == "component-nef"
assert branch["status"] == "PASS_EXACT_D13_BRANCH_ANCHOR"
assert g3data["status"] == "PASS_EXACT_D13_G3_FROM_E77_BISECTION"

U2 = matrix(ZZ, ((0, 1), (1, 0)))
ns = block_diagonal_matrix(U2, -frame)
F0 = vector(ZZ, [1, 0] + [0]*17)
O0 = vector(ZZ, [-1, 1] + [0]*17)

# Canonical source E7+E8 simple roots.
source_roots = matrix(
    ZZ,
    [
        [0, 0] + [ZZ(i == node) for i in range(17)]
        for node in range(15)
    ],
)
assert all(root * ns * F0 == 0 for root in source_roots.rows())
assert all(root * ns * O0 == 0 for root in source_roots.rows())

Groot = source_roots * ns * source_roots.transpose()
source_e7_cartan = -Groot[:7, :7]
source_e8_cartan = -Groot[7:, 7:]
assert source_e7_cartan.det() == 2
assert source_e8_cartan.det() == 1
assert Groot[:7, 7:] == matrix(ZZ, 7, 8)
assert Groot[7:, :7] == matrix(ZZ, 8, 7)
assert abs(Groot.det()) == 2  # E7 det 2 times E8 det 1.

# Source-frame ordering is NOT the formal E7 resolution numbering.  Pin its
# highest root using the coefficients already used by the passing q6
# preflight, and verify them intrinsically in this ordered Cartan matrix.
highest_e7 = vector(ZZ, (2, 2, 3, 4, 3, 2, 1))
affine_e7 = F0 - highest_e7 * source_roots[:7, :]
assert affine_e7 * ns * affine_e7 == -2
assert affine_e7 * ns * F0 == 0
assert affine_e7 * ns * O0 == 1
assert all(affine_e7 * ns * source_roots[7:, :].row(i) == 0 for i in range(8))

# Actual source sections.  The second stored coordinate is the H3-frame P2
# class, not the reconstructed affine sign convention used in local Hensel
# files.
minus_p1 = vector(ZZ, [5, 1] + [
    -2, -3, -4, -6, -5, -4, -3,
    0, 0, 0, 0, 0, 0, 0, 0,
    1, 0,
])
p2 = vector(ZZ, [22, 1] + [0]*16 + [1])

for name, P in (("minus_P1", minus_p1), ("P2", p2)):
    assert len(P) == 19
    assert P * ns * P == -2, (name, P*ns*P)
    assert P * ns * F0 == 1, (name, P*ns*F0)

projection = identity_matrix(QQ, 19) - ns*source_roots.transpose()*Groot.inverse()*source_roots

def shioda(P):
    horizontal = P - O0 - (P*ns*O0 + 2)*F0
    assert horizontal * ns * F0 == 0
    assert horizontal * ns * O0 == 0
    return vector(QQ, horizontal) * projection

basis_sections = (minus_p1, p2)
phis = tuple(shioda(P) for P in basis_sections)
H0 = matrix(QQ, [[-left*ns*right for right in phis] for left in phis])

# Trust the pinned NS quotient basis, not a sign inference from the human
# point labels.  verify_h92_section_descent.sage independently certifies this
# exact H3 quotient Gram from the stored frame.
assert H0 == matrix(QQ, [[QQ(21)/2, 3], [3, 46]]), H0

# III* has exactly two multiplicity-one components: affine and source E7_7.
# The coefficient vector above has a unique entry 1, at source node 7.
# II* has only the affine multiplicity-one component.
multiplicity_one_e7 = [i for i, value in enumerate(highest_e7) if value == 1]
assert multiplicity_one_e7 == [6]
patterns = [
    vector(ZZ, [0]*15),  # affine III* and affine II*
    vector(ZZ, [0]*6 + [1] + [0]*8),  # source E7_7
]


def candidate_for_pattern(n, p):
    n = vector(ZZ, n)
    p = vector(ZZ, p)
    phi = QQ(n[0])*phis[0] + QQ(n[1])*phis[1]
    height = QQ(n*H0*n)
    correction = -QQ(p*Groot.inverse()*p)
    pole = (height + correction - 4)/2
    if pole not in ZZ or pole < 0:
        return None
    root_part = vector(QQ, p)*Groot.inverse()*source_roots
    P = (
        vector(QQ, O0)
        + (QQ(pole)+2)*vector(QQ, F0)
        + phi + root_part
    )
    if not all(value in ZZ for value in P):
        return None
    P = vector(ZZ, P)
    if P*ns*P != -2 or P*ns*F0 != 1:
        return None
    if vector(ZZ, P*ns*source_roots.transpose()) != p:
        return None
    return P, ZZ(pole), correction, p


def section_from_word(n):
    solutions = []
    for p in patterns:
        candidate = candidate_for_pattern(n, p)
        if candidate is not None:
            solutions.append(candidate)
    assert len(solutions) == 1, (tuple(n), len(solutions))
    return solutions[0]


# Exact basis regressions.
assert section_from_word((1, 0))[0] == minus_p1
assert section_from_word((0, 1))[0] == p2

print(
    "H3MWQ24_SETUP|basis=(section_class_A,section_class_B)|height={}|basis_reconstruction=PASS".format(H0),
    flush=True,
)

# Anchored q8/D13 reader.
F8 = vector(ZZ, target["selected_q8"]["source_h3_ns_vector"])
E6 = matrix(
    ZZ, target["selected_q8"]["E6"]["simple_root_vectors_in_source_h3_ns"]
)
E8 = matrix(
    ZZ, target["selected_q8"]["E8"]["simple_root_vectors_in_source_h3_ns"]
)
O8 = vector(ZZ, E8.row(0))
G1_curve = vector(ZZ, E6.row(0))

child8, Bzero = child_frame_with_zero(ns, F8, O8)
A13, adapted, H13 = d13_root_adaptation(child8)
Badapt = block_diagonal_matrix(identity_matrix(ZZ, 2), A13) * Bzero
Gadapt = block_diagonal_matrix(U2, -adapted)
assert Badapt*ns*Badapt.transpose() == Gadapt

def d13_coords(P):
    value = vector(QQ, P)*Badapt.inverse()
    assert all(entry in ZZ for entry in value)
    return vector(ZZ, value)

G1 = vector(ZZ, d13_coords(G1_curve)[-4:])
assert G1 == vector(ZZ, (1, 0, 0, 0))
G3 = vector(ZZ, (0, 0, 1, 0))
q24 = vector(ZZ, (2, -1, -1, 1))
assert G3*H13*G3 == QQ(11)/4
assert q24*H13*q24 == 52

# Reproduce the already observed basis section q8 profiles.
for word, expected_degree, expected_z in (
    ((1, 0), 19, (2, 1, 4, 0)),
    ((0, 1), 350, (22, 0, 66, -1)),
):
    P = section_from_word(word)[0]
    assert P*ns*F8 == expected_degree
    assert vector(ZZ, d13_coords(P)[-4:]) == vector(ZZ, expected_z)

print(
    "H3MWQ24_ANCHOR|minus_P1_degree=19|P2_degree=350|"
    "profiles=PASS|q24=2,-1,-1,1",
    flush=True,
)

# Exact q8-degree quadratics.
#
# For a source section:
#   P = O + ((height+correction)/2)*F0 + phi + root_part.
# Therefore the quadratic height contribution to P.F8 is scaled by
# (F0.F8)/2.  Do not reuse the raw-q6 special case F6.F8=2.
fiber_intersection = ZZ(F0*ns*F8)
assert fiber_intersection > 0
height_scale = QQ(fiber_intersection)/2
linear = vector(QQ, [phi*ns*F8 for phi in phis])
Qdegree = height_scale * H0
Qinv = Qdegree.inverse()
center = -QQ(1)/2 * linear * Qinv
threshold = QQ(args.threshold)

pattern_quadratics = []
global_lower = [None]*2
global_upper = [None]*2

for pattern_index, p in enumerate(patterns):
    correction = -QQ(p*Groot.inverse()*p)
    root_part = vector(QQ, p)*Groot.inverse()*source_roots
    constant = (
        QQ(O0*ns*F8)
        + height_scale*correction
        + QQ(root_part*ns*F8)
    )
    minimum = constant - QQ(1)/4*(linear*Qinv*linear)
    budget = threshold - minimum
    if budget <= 0:
        continue
    bounds = []
    for i in range(2):
        radius = ceil_sqrt_rational(budget*Qinv[i, i])
        lo = ZZ(center[i].floor()) - radius - 1
        hi = ZZ(center[i].ceil()) + radius + 1
        bounds.append((lo, hi))
        global_lower[i] = lo if global_lower[i] is None else min(global_lower[i], lo)
        global_upper[i] = hi if global_upper[i] is None else max(global_upper[i], hi)
    pattern_quadratics.append({
        "pattern_index": pattern_index,
        "pattern": p,
        "correction": correction,
        "constant": constant,
        "minimum": minimum,
        "bounds": bounds,
    })

assert pattern_quadratics
assert all(v is not None for v in global_lower + global_upper)

# Check closed formulas on the two source basis sections.
for word in ((1, 0), (0, 1)):
    n = vector(ZZ, word)
    P, pole, correction, p = section_from_word(n)
    root_part = vector(QQ, p)*Groot.inverse()*source_roots
    formula = (
        height_scale*(QQ(n*H0*n) + correction)
        + linear*n + QQ(O0*ns*F8)
        + QQ(root_part*ns*F8)
    )
    assert formula == P*ns*F8

print(
    "H3MWQ24_BOUND|threshold={}|F0dotF8={}|height_scale={}|center={}|box={}|patterns={}".format(
        args.threshold,
        fiber_intersection,
        height_scale,
        ",".join(map(str, center)),
        ";".join(
            f"{global_lower[i]}..{global_upper[i]}" for i in range(2)
        ),
        len(pattern_quadratics),
    ),
    flush=True,
)

best_fourth = []
best_compatible = []
exact_q24 = []
all_low = []
tested = 0

for a in range(int(global_lower[0]), int(global_upper[0])+1):
    for b in range(int(global_lower[1]), int(global_upper[1])+1):
        if a == b == 0:
            continue
        n = vector(ZZ, (a, b))
        P, pole, correction, p = section_from_word(n)
        tested += 1
        degree8 = int(P*ns*F8)
        if degree8 >= args.threshold:
            continue

        root_part = vector(QQ, p)*Groot.inverse()*source_roots
        formula = (
            height_scale*(QQ(n*H0*n) + correction)
            + linear*n + QQ(O0*ns*F8)
            + QQ(root_part*ns*F8)
        )
        assert formula == degree8
        assert degree8 >= 0

        z = vector(ZZ, d13_coords(P)[-4:])
        record = {
            "source_mw_word_minusP1_P2": [a, b],
            "source_height": str(n*H0*n),
            "source_P_dot_O": int(P*ns*O0),
            "source_component_correction": str(correction),
            "q8_degree": degree8,
            "q8_zero_intersection": int(P*ns*O8),
            "d13_mw": list(map(int, z)),
            "source_h3_ns": list(map(int, P)),
        }
        all_low.append(record)

        if z[3] != 0:
            best_fourth.append(record)

        if z[1] == -1 and z[3] == 1:
            residual = q24-z
            assert residual[1] == 0 and residual[3] == 0
            record["q24_correction_G1"] = int(residual[0])
            record["q24_correction_G3"] = int(residual[2])
            best_compatible.append(record)

        if z == q24:
            exact_q24.append(record)

key = lambda r: (r["q8_degree"], sum(abs(v) for v in r["source_mw_word_minusP1_P2"]))
all_low.sort(key=key)
best_fourth.sort(key=key)
best_compatible.sort(key=key)
exact_q24.sort(key=key)

print(
    f"H3MWQ24_SEARCH|tested_box={tested}|below_threshold={len(all_low)}|"
    f"fourth_candidates={len(best_fourth)}|"
    f"q24_compatible={len(best_compatible)}|exact_q24={len(exact_q24)}",
    flush=True,
)

def emit(prefix, records):
    for rank, record in enumerate(records[:args.top], 1):
        extra = ""
        if "q24_correction_G1" in record:
            extra = (
                f"q24_add={record['q24_correction_G1']}*G1+"
                f"{record['q24_correction_G3']}*G3|"
            )
        print(
            f"{prefix}|rank={rank}|"
            f"word={','.join(map(str,record['source_mw_word_minusP1_P2']))}|"
            f"q8_degree={record['q8_degree']}|"
            f"d13_mw={','.join(map(str,record['d13_mw']))}|"
            f"source_height={record['source_height']}|"
            f"source_PdotO={record['source_P_dot_O']}|"
            + extra
            + "status=CANDIDATE",
            flush=True,
        )

emit("H3MWQ24_LOW", all_low)
emit("H3MWQ24_FOURTH", best_fourth)
emit("H3MWQ24_COMPAT", best_compatible)
emit("H3MWQ24_EXACT", exact_q24)

if best_compatible:
    status = "PASS_FOUND_SOURCE_Q24_BRIDGE"
else:
    status = "PASS_NO_SOURCE_Q24_BRIDGE_BELOW_THRESHOLD"

print(
    "H3MWQ24_RESULT|status={}|best_compatible_degree={}|"
    "threshold={}|exhaustive_below_threshold=1|q6_reference_degree=46".format(
        status,
        best_compatible[0]["q8_degree"] if best_compatible else "none",
        args.threshold,
    ),
    flush=True,
)

payload = {
    "schema": "elkies-k3.h92-source-mw-q24-bridge-search.v1",
    "status": status,
    "basis": {
        "coordinates": ["section_class_A (historically -P1)", "section_class_B (historically P2)"],
        "height_gram": [[str(v) for v in row] for row in H0.rows()],
        "basis_reconstruction_exact": True,
    },
    "search": {
        "threshold": args.threshold,
        "tested_box": tested,
        "below_threshold": len(all_low),
        "source_fibre_intersection_with_q8": int(fiber_intersection),
        "height_quadratic_scale": str(height_scale),
        "completed_square_center": [str(v) for v in center],
        "global_coordinate_bounds": [
            [int(global_lower[i]), int(global_upper[i])] for i in range(2)
        ],
        "pattern_quadratics": [
            {
                "pattern_index": item["pattern_index"],
                "pairing_pattern": list(map(int, item["pattern"])),
                "component_correction": str(item["correction"]),
                "constant": str(item["constant"]),
                "real_minimum": str(item["minimum"]),
                "coordinate_bounds": [
                    [int(a), int(b)] for a, b in item["bounds"]
                ],
            }
            for item in pattern_quadratics
        ],
        "exhaustive_for_q8_degree_below_threshold": True,
    },
    "target": {
        "q24_D13_mw": [2, -1, -1, 1],
        "G1": [1, 0, 0, 0],
        "G3": [0, 0, 1, 0],
        "q6_best_known_bridge_degree": 46,
    },
    "all_low_degree": all_low[:args.top],
    "best_fourth_direction": best_fourth[:args.top],
    "best_q24_compatible": best_compatible[:args.top],
    "exact_q24": exact_q24[:args.top],
    "boundary": (
        "This is an exact exhaustive lattice/effective-section search in the "
        "original H3 fibration below the stated q8-degree threshold. It does "
        "not evaluate any newly found source MW combination by group law."
    ),
}

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(f"OUTPUT|{OUTPUT}")
