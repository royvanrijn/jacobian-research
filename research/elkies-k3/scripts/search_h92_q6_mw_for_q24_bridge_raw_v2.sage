#!/usr/bin/env sage -python
"""
Exhaustively search the ACTUAL raw-q6 MW section lattice for a cheap bridge to
q24 in the anchored q8/D13 fibration.

This version fixes a presentation-mixing bug in the previous search.

The explicit q6 equation / physical-root certificate uses:
    F6 = the raw q6 isotropic class [3,2,w],
    O6 = the old H3 zero transported BACK through the 22 Weyl reflections,
    E6+E8 roots = the exact vectors stored in q8-target-component-nef.json.

The previous script incorrectly used Weyl(F6) together with those raw-presentation
roots.  Here F6, O6, roots, MW Shioda vectors, and searched section classes all
live in the same raw-q6 presentation.

No duplicate qfminim root reconstruction is performed: the exact E6/E8 root
vectors are consumed from the already-passing physical-target artifact.

The q6 MW basis is reconstructed from the last three rows of the stored
root_mw_basis_in_child in elkies-k3-h3-q6-q8-orbits.json, avoiding copied
constants.  For every MW word n, the unique integral section is recovered from
its Shioda vector plus one of the three IV* multiplicity-one component classes.

For each fixed IV* component class p, q8 degree is exactly

    d8(n,p) = n^T H6 n + ell*n + c_p.

Completing this positive-definite quadratic gives a finite box containing EVERY
q6 section with d8 < --threshold.  The enumeration is therefore exhaustive
below the threshold, not heuristic.

We search for anchored D13 Abel-Jacobi coordinates z satisfying
    z[1] = -1 and z[3] = 1,
because q24=(2,-1,-1,1) then differs from z only by the already-explicit
G1=(1,0,0,0) and G3=(0,0,1,0).

Run:
  sage -python ~/Downloads/search_h92_q6_mw_for_q24_bridge_raw_v2.sage
"""

import argparse
import json
from math import isqrt
from pathlib import Path

from sage.all import (
    QQ, ZZ, block_diagonal_matrix, identity_matrix, lcm, matrix, pari,
    vector, xgcd,
)

Q6_REFLECTIONS = (
    1, 2, 4, 3, 5, 4, 2, 6, 5, 4, 3, 1,
    7, 6, 5, 4, 2, 3, 4, 5, 6, 7,
)

E6_CARTAN = matrix(ZZ, [
    [2, -1, 0, 0, 0, 0],
    [-1, 2, -1, 0, 0, 0],
    [0, -1, 2, -1, 0, -1],
    [0, 0, -1, 2, -1, 0],
    [0, 0, 0, -1, 2, 0],
    [0, 0, -1, 0, 0, 2],
])
E8_CARTAN = matrix(ZZ, [
    [2, 0, -1, 0, 0, 0, 0, 0],
    [0, 2, 0, -1, 0, 0, 0, 0],
    [-1, 0, 2, -1, 0, 0, 0, 0],
    [0, -1, -1, 2, -1, 0, 0, 0],
    [0, 0, 0, -1, 2, -1, 0, 0],
    [0, 0, 0, 0, -1, 2, -1, 0],
    [0, 0, 0, 0, 0, -1, 2, -1],
    [0, 0, 0, 0, 0, 0, -1, 2],
])
EXPECTED_Q6_HEIGHT = matrix(QQ, [
    [QQ(8)/3, QQ(1)/3, -1],
    [QQ(1)/3, QQ(8)/3, 1],
    [-1, 1, 46],
])


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


def reflect(row, gram, root):
    # Always make a genuine copy before mutation-sensitive operations.
    row = vector(ZZ, list(row))
    root = vector(ZZ, list(root))
    assert root * gram * root == -2
    return row + (row * gram * root) * root


def isotropic_mate(ns, fibre):
    current = ZZ(0)
    data = [ZZ(0)] * ns.nrows()
    for index, value in enumerate(ns * fibre):
        if not value:
            continue
        divisor, left, right = xgcd(current, ZZ(value))
        data = [left * entry for entry in data]
        data[index] += right
        current = divisor
    assert abs(current) == 1
    if current == -1:
        data = [-entry for entry in data]
    mate = vector(ZZ, data)
    mate -= (mate * ns * mate // 2) * fibre
    assert mate * ns * mate == 0 and mate * ns * fibre == 1
    return mate


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


def highest_root(cartan):
    half = matrix(ZZ, pari(cartan).qfminim(2)[2]).transpose().rows()
    roots = [vector(ZZ, r) for r in half]
    roots += [-vector(ZZ, r) for r in half]
    positive = [r for r in roots if all(v >= 0 for v in r)]
    assert positive
    return max(positive, key=lambda r: sum(r))


def ceil_sqrt_rational(value):
    value = QQ(value)
    assert value >= 0
    num = ZZ(value.numerator())
    den = ZZ(value.denominator())
    q = (num + den - 1) // den
    k = ZZ(isqrt(int(q)))
    if k*k < q:
        k += 1
    assert QQ(k*k) >= value
    return k


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo", type=Path)
parser.add_argument(
    "--threshold", type=int, default=206,
    help="exhaustively enumerate every raw-q6 section with q8 degree below this",
)
parser.add_argument("--top", type=int, default=20)
parser.add_argument("--output", type=Path)
args = parser.parse_args()

ROOT = locate_repo(args.repo)
GEN = ROOT / "artifacts/generated-results"
LOCAL = ROOT / "artifacts/local/elkies-k3"
FRAME = ROOT / "elkies-k3/data/fibrations/kumar_e7e8_mw2_frame_3.txt"
ORBITS = GEN / "elkies-k3-h3-q6-q8-orbits.json"
TARGET = LOCAL / "q8-target-component-nef.json"
BRANCH = LOCAL / "q8-d13-branch-anchor.json"
G3FILE = LOCAL / "q8-d13-g3-from-e77-bisection.json"
OUTPUT = (
    args.output.resolve()
    if args.output
    else LOCAL / "q6-raw-mw-q24-bridge-search.json"
)

for path in (FRAME, ORBITS, TARGET, BRANCH, G3FILE):
    if not path.exists():
        raise SystemExit(f"Missing prerequisite: {path}")

frame = load_gram(FRAME)
orbits = json.loads(ORBITS.read_text())
target = json.loads(TARGET.read_text())
branch = json.loads(BRANCH.read_text())
g3data = json.loads(G3FILE.read_text())
assert orbits["status"] == "PASS_H3_Q6_CHILD_Q8_WEYL_CLASSIFICATION"
assert target["status"] == "PASS_EXACT_Q6_CHILD_Q8_PHYSICAL_ROOT_TARGET"
assert target["normalization"]["representative"] == "component-nef"
assert branch["status"] == "PASS_EXACT_D13_BRANCH_ANCHOR"
assert g3data["status"] == "PASS_EXACT_D13_G3_FROM_E77_BISECTION"

U2 = matrix(ZZ, ((0, 1), (1, 0)))
ns = block_diagonal_matrix(U2, -frame)
source_O = vector(ZZ, [-1, 1] + [0]*17)
source_simple = tuple(
    vector(ZZ, [0, 0] + [ZZ(i == node) for i in range(17)])
    for node in range(15)
)
reflection_roots = tuple(source_simple[node-1] for node in Q6_REFLECTIONS)

# ---------------------------------------------------------------------------
# 1. EXACT raw-q6 presentation used by the explicit child/physical target.
# ---------------------------------------------------------------------------

F6 = vector(ZZ, [3, 2] + [
    0, 0, -1, -1, -1, -1, -1,
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
])
O6 = vector(ZZ, list(source_O))
for root in reversed(reflection_roots):
    O6 = reflect(O6, ns, root)

assert F6 * ns * F6 == 0
assert O6 * ns * O6 == -2
assert O6 * ns * F6 == 1

# Consume physical roots from the already-passing target artifact.  Do not
# reconstruct another qfminim basis here.
e6_roots = matrix(
    ZZ, target["selected_q8"]["E6"]["simple_root_vectors_in_source_h3_ns"]
)
e8_roots = matrix(
    ZZ, target["selected_q8"]["E8"]["simple_root_vectors_in_source_h3_ns"]
)
assert -e6_roots * ns * e6_roots.transpose() == E6_CARTAN
assert -e8_roots * ns * e8_roots.transpose() == E8_CARTAN
assert e6_roots * ns * e8_roots.transpose() == matrix(ZZ, 6, 8)
assert all(root * ns * F6 == 0 for root in e6_roots.rows())
assert all(root * ns * F6 == 0 for root in e8_roots.rows())
assert all(root * ns * O6 == 0 for root in e6_roots.rows())
assert all(root * ns * O6 == 0 for root in e8_roots.rows())

Rroots = e6_roots.stack(e8_roots)
Groot = Rroots * ns * Rroots.transpose()
assert Groot.det() == 3

# ---------------------------------------------------------------------------
# 2. Saturated raw-q6 MW Shioda basis from the stored q6 child certificate.
# ---------------------------------------------------------------------------

raw_mate = isotropic_mate(ns, F6)
raw_orth = matrix(
    ZZ, [list(F6*ns), list(raw_mate*ns)]
).right_kernel_matrix()
raw_transport = matrix(
    ZZ, [list(F6), list(raw_mate)] + [list(row) for row in raw_orth.rows()]
)
assert abs(raw_transport.det()) == 1
raw_child = -(raw_orth * ns * raw_orth.transpose())
assert raw_child.det() == 948

root_mw_basis = matrix(ZZ, orbits["q6"]["root_mw_basis_in_child"])
assert abs(root_mw_basis.det()) == 1
assert (
    root_mw_basis * raw_child * root_mw_basis.transpose()
    == matrix(ZZ, orbits["q6"]["root_adapted_gram"])
)
mw_lifts = root_mw_basis[14:17, :]
assert mw_lifts.nrows() == 3 and mw_lifts.ncols() == 17

# A lift gives the MW direction in the raw child complement.  Root shifts do
# not affect the Shioda class, so these convenient (-2) section classes are
# sufficient to pin the saturated MW basis.
raw_basis_sections = []
for lift in mw_lifts.rows():
    perpendicular = vector(ZZ, lift) * raw_orth
    start = O6 + perpendicular
    delta = ZZ(-2 - start*ns*start)
    assert delta % 2 == 0
    section = start + (delta//2) * F6
    assert section * ns * section == -2
    assert section * ns * F6 == 1
    raw_basis_sections.append(section)

projection = identity_matrix(QQ, 19) - ns*Rroots.transpose()*Groot.inverse()*Rroots

def shioda(P):
    horizontal = P - O6 - (P*ns*O6 + 2)*F6
    assert horizontal * ns * F6 == 0
    assert horizontal * ns * O6 == 0
    return vector(QQ, horizontal) * projection

phis = [shioda(P) for P in raw_basis_sections]
H6 = matrix(QQ, [[-left*ns*right for right in phis] for left in phis])
assert H6 == EXPECTED_Q6_HEIGHT

# IV* has three multiplicity-one components: affine + two outer simple roots.
# II* contributes only its affine multiplicity-one component.
e6_high = highest_root(E6_CARTAN)
outer_indices = [i for i, value in enumerate(e6_high) if value == 1]
assert len(outer_indices) == 2, (e6_high, outer_indices)
pairing_patterns = [vector(ZZ, [0]*14)]
for index in outer_indices:
    data = [0]*14
    data[index] = 1
    pairing_patterns.append(vector(ZZ, data))


def candidate_for_pattern(n, p):
    n = vector(ZZ, n)
    p = vector(ZZ, p)
    phi = sum((QQ(n[i])*phis[i] for i in range(3)), vector(QQ, [0]*19))
    height = QQ(n*H6*n)
    correction = -QQ(p*Groot.inverse()*p)
    pole = (height + correction - 4)/2
    if pole not in ZZ or pole < 0:
        return None
    root_part = vector(QQ, p) * Groot.inverse() * Rroots
    P = (
        vector(QQ, O6)
        + (QQ(pole)+2)*vector(QQ, F6)
        + phi + root_part
    )
    if not all(value in ZZ for value in P):
        return None
    P = vector(ZZ, P)
    if P*ns*P != -2 or P*ns*F6 != 1:
        return None
    if vector(ZZ, P*ns*Rroots.transpose()) != p:
        return None
    return P, ZZ(pole), correction, p


def section_from_word(n):
    solutions = []
    for p in pairing_patterns:
        candidate = candidate_for_pattern(n, p)
        if candidate is not None:
            solutions.append(candidate)
    assert len(solutions) == 1, (tuple(n), len(solutions))
    return solutions[0]

# Validate the section reconstruction on the known explicit old_E7_7 curve.
def inverse_q6_weyl(curve):
    result = vector(ZZ, list(curve))
    for root in reversed(reflection_roots):
        result = reflect(result, ns, root)
    return result

old_e77 = inverse_q6_weyl(source_simple[6])
assert old_e77 * ns * old_e77 == -2
assert old_e77 * ns * F6 == 1
phi_e77 = shioda(old_e77)
pairings_e77 = vector(QQ, [-phi_e77*ns*phi for phi in phis])
word_e77 = pairings_e77 * H6.inverse()
assert all(value in ZZ for value in word_e77)
word_e77 = vector(ZZ, word_e77)
reconstructed_e77 = section_from_word(word_e77)[0]
assert reconstructed_e77 == old_e77

print(
    "Q6RAWSEARCH_SETUP|presentation=raw-q6|physical_roots=artifact|"
    f"E77_word={','.join(map(str,word_e77))}|"
    "E77_reconstruction=PASS|q6_height=PASS",
    flush=True,
)

# ---------------------------------------------------------------------------
# 3. Anchored D13 coordinate reader and known equation-level regressions.
# ---------------------------------------------------------------------------

F8 = vector(ZZ, target["selected_q8"]["source_h3_ns_vector"])
O8 = vector(ZZ, e8_roots.row(0))
G1_curve = vector(ZZ, e6_roots.row(0))
assert F8 * ns * F6 == 2
assert O8 * ns * F8 == 1
assert G1_curve * ns * F8 == 1

child8, Bzero = child_frame_with_zero(ns, F8, O8)
A13, adapted, H13 = d13_root_adaptation(child8)
Badapt = block_diagonal_matrix(identity_matrix(ZZ, 2), A13) * Bzero
Gadapt = block_diagonal_matrix(U2, -adapted)
assert Badapt * ns * Badapt.transpose() == Gadapt


def d13_coords(P):
    value = vector(QQ, P) * Badapt.inverse()
    assert all(entry in ZZ for entry in value)
    return vector(ZZ, value)

G1 = vector(ZZ, d13_coords(G1_curve)[-4:])
assert G1 == vector(ZZ, (1, 0, 0, 0))
G3 = vector(ZZ, (0, 0, 1, 0))
q24 = vector(ZZ, (2, -1, -1, 1))
assert G1 * H13 * G1 == QQ(3)/4
assert G3 * H13 * G3 == QQ(11)/4
assert q24 * H13 * q24 == 52

z_e77 = vector(ZZ, d13_coords(old_e77)[-4:])
assert old_e77 * ns * F8 == 2
assert z_e77 == vector(ZZ, (1, 0, 1, 0))
assert g3data["lattice_certificate"]["AJ_old_E7_7"] == [1, 0, 1, 0]
assert g3data["lattice_certificate"]["G3"] == [0, 0, 1, 0]

print(
    "Q6RAWSEARCH_ANCHOR|G1=1,0,0,0|G3=0,0,1,0|"
    "E77_AJ=1,0,1,0|q24=2,-1,-1,1|status=PASS",
    flush=True,
)

# ---------------------------------------------------------------------------
# 4. Exact positive-definite q8-degree quadratics and exhaustive bound.
# ---------------------------------------------------------------------------

linear = vector(QQ, [phi*ns*F8 for phi in phis])
Hinv = H6.inverse()
center = -QQ(1)/2 * linear * Hinv
threshold = QQ(args.threshold)

pattern_quadratics = []
global_lower = [None]*3
global_upper = [None]*3

for pattern_index, p in enumerate(pairing_patterns):
    p = vector(ZZ, p)
    correction = -QQ(p*Groot.inverse()*p)
    root_part = vector(QQ, p) * Groot.inverse() * Rroots
    constant = QQ(O6*ns*F8) + correction + QQ(root_part*ns*F8)
    minimum = constant - QQ(1)/4*(linear*Hinv*linear)
    budget = threshold - minimum
    if budget <= 0:
        continue
    bounds = []
    for i in range(3):
        radius = ceil_sqrt_rational(budget*Hinv[i, i])
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
assert all(value is not None for value in global_lower + global_upper)

# Closed-form q8-degree regression on E7_7.
P_e77, pole_e77, corr_e77, pattern_e77 = section_from_word(word_e77)
root_e77 = vector(QQ, pattern_e77)*Groot.inverse()*Rroots
formula_e77 = (
    QQ(word_e77*H6*word_e77)
    + linear*word_e77
    + QQ(O6*ns*F8) + corr_e77 + QQ(root_e77*ns*F8)
)
assert formula_e77 == P_e77*ns*F8 == 2

print(
    "Q6RAWSEARCH_BOUND|threshold={}|center={}|box={}|patterns={}|"
    "E77_degree_formula=2".format(
        args.threshold,
        ",".join(map(str, center)),
        ";".join(f"{global_lower[i]}..{global_upper[i]}" for i in range(3)),
        len(pattern_quadratics),
    ),
    flush=True,
)

# ---------------------------------------------------------------------------
# 5. Exhaustive enumeration below threshold.
# ---------------------------------------------------------------------------

best_fourth = []
best_compatible = []
exact_q24 = []
tested_box = 0
below_threshold = 0

for n1 in range(int(global_lower[0]), int(global_upper[0])+1):
    for n2 in range(int(global_lower[1]), int(global_upper[1])+1):
        for n3 in range(int(global_lower[2]), int(global_upper[2])+1):
            if n1 == n2 == n3 == 0:
                continue
            n = vector(ZZ, (n1, n2, n3))
            P, pole, correction, p = section_from_word(n)
            tested_box += 1
            degree8 = int(P*ns*F8)
            if degree8 >= args.threshold:
                continue
            below_threshold += 1

            root_part = vector(QQ, p)*Groot.inverse()*Rroots
            degree_formula = (
                QQ(n*H6*n) + linear*n + QQ(O6*ns*F8)
                + correction + QQ(root_part*ns*F8)
            )
            assert degree_formula == degree8
            assert degree8 >= 0

            c = d13_coords(P)
            z = vector(ZZ, c[-4:])
            record = {
                "q6_word": [n1, n2, n3],
                "q6_height": str(n*H6*n),
                "q6_P_dot_O": int(P*ns*O6),
                "q6_component_correction": str(correction),
                "q8_degree": degree8,
                "q8_zero_intersection": int(P*ns*O8),
                "d13_mw": list(map(int, z)),
                "source_h3_ns_raw_presentation": list(map(int, P)),
            }

            if z[3] != 0:
                best_fourth.append(record)

            if z[1] == -1 and z[3] == 1:
                residual = q24 - z
                assert residual[1] == 0 and residual[3] == 0
                record["q24_correction_G1"] = int(residual[0])
                record["q24_correction_G3"] = int(residual[2])
                best_compatible.append(record)

            if z == q24:
                exact_q24.append(record)

best_fourth.sort(key=lambda r: (r["q8_degree"], sum(abs(v) for v in r["q6_word"])))
best_compatible.sort(key=lambda r: (r["q8_degree"], sum(abs(v) for v in r["q6_word"])))
exact_q24.sort(key=lambda r: r["q8_degree"])

print(
    f"Q6RAWSEARCH|tested_box={tested_box}|below_threshold={below_threshold}|"
    f"threshold={args.threshold}|fourth_candidates={len(best_fourth)}|"
    f"q24_compatible={len(best_compatible)}|exact_q24={len(exact_q24)}",
    flush=True,
)


def emit(prefix, records):
    for rank, record in enumerate(records[:args.top], 1):
        correction_text = ""
        if "q24_correction_G1" in record:
            correction_text = (
                f"q24_add={record['q24_correction_G1']}*G1+"
                f"{record['q24_correction_G3']}*G3|"
            )
        print(
            f"{prefix}|rank={rank}|"
            f"word={','.join(map(str,record['q6_word']))}|"
            f"q8_degree={record['q8_degree']}|"
            f"d13_mw={','.join(map(str,record['d13_mw']))}|"
            f"q6_height={record['q6_height']}|"
            f"q6_PdotO={record['q6_P_dot_O']}|"
            + correction_text
            + "status=CANDIDATE",
            flush=True,
        )

emit("Q6RAWSEARCH_FOURTH", best_fourth)
emit("Q6RAWSEARCH_COMPAT", best_compatible)
emit("Q6RAWSEARCH_EXACT", exact_q24)

status = (
    "PASS_FOUND_RAW_Q6_Q24_BRIDGE"
    if best_compatible
    else "PASS_NO_RAW_Q6_Q24_BRIDGE_BELOW_THRESHOLD"
)
print(
    "Q6RAWSEARCH_RESULT|status={}|best_compatible_degree={}|"
    "threshold={}|exhaustive_below_threshold=1".format(
        status,
        best_compatible[0]["q8_degree"] if best_compatible else "none",
        args.threshold,
    ),
    flush=True,
)

payload = {
    "schema": "elkies-k3.h92-q6-raw-mw-q24-bridge-search.v2",
    "status": status,
    "presentation": {
        "q6_fibre": "raw q6 isotropic class used by explicit child",
        "q6_zero": "source old zero transported backward through q6 Weyl record",
        "physical_roots": "consumed from q8-target-component-nef artifact",
        "presentation_consistency_checked": True,
    },
    "regressions": {
        "old_E7_7_q6_word": list(map(int, word_e77)),
        "old_E7_7_reconstructed_exactly": True,
        "old_E7_7_q8_degree": 2,
        "old_E7_7_D13_AJ": [1, 0, 1, 0],
        "G1": [1, 0, 0, 0],
        "G3": [0, 0, 1, 0],
    },
    "search": {
        "threshold": args.threshold,
        "tested_box": tested_box,
        "below_threshold": below_threshold,
        "completed_square_center": [str(v) for v in center],
        "global_coordinate_bounds": [
            [int(global_lower[i]), int(global_upper[i])] for i in range(3)
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
        "q24_anchored_D13_mw": [2, -1, -1, 1],
        "explicit_correction_basis": {
            "G1": [1, 0, 0, 0],
            "G3": [0, 0, 1, 0],
        },
        "compatibility_condition": "z[1]=-1 and z[3]=1",
    },
    "best_fourth_direction": best_fourth[:args.top],
    "best_q24_compatible": best_compatible[:args.top],
    "exact_q24": exact_q24[:args.top],
    "boundary": (
        "This exhaustively searches actual section classes in the raw q6 "
        "presentation below the stated q8-degree threshold. It does not "
        "perform Abel-Jacobi reduction of any degree>1 bridge on the q8 quartic."
    ),
}
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(f"OUTPUT|{OUTPUT}")
