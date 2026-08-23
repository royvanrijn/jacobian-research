#!/usr/bin/env sage -python
"""
Profile already-explicit old-H3 curves as sections/multisections of the
anchored q8/D13 fibration.

The anchored q8 zero is II*_E8_1 and IV*_E6_1 is the explicit height-3/4
generator G1=(1,0,0,0).  The previous old-q6 MW-image audit showed that

    <G1, A, B, C>  has index 3 in MW(D13),

with
    A=(5,-2,11,0),
    B=(3,-1,7,0),
    C=(22,-1,45,1).

Hence the unique missing saturation class may be represented by
    G3=(0,0,1,0),
because 3*G3 = 2*B - A - G1.

This script tests the rational curves we already have explicitly from the
original E7 fibre:
  * transported old zero;
  * transported old E7 simple components 1..7;
  * transported old E7 affine component.

For each curve it computes q8 degree, anchored D13 MW coordinates, height,
intersection with the anchored zero, and its class in the index-3 quotient.

If q8_degree=1, the curve itself is already an explicit D13 section.  A
nonzero quotient class then saturates the missing mod-3 direction without any
new section search.

Run:
  sage -python ~/Downloads/audit_h92_q8_explicit_old_curves.sage
"""

import argparse
import json
from pathlib import Path

from sage.all import (
    GF, QQ, ZZ, block_diagonal_matrix, identity_matrix, lcm, matrix, pari, vector
)

Q6_REFLECTIONS = (
    1, 2, 4, 3, 5, 4, 2, 6, 5, 4, 3, 1,
    7, 6, 5, 4, 2, 3, 4, 5, 6, 7,
)


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
        if (
            (candidate / "elkies-k3" / "scripts").is_dir()
            and (candidate / "artifacts" / "generated-results").is_dir()
        ):
            return candidate
    raise SystemExit("Could not locate jacobian-research; pass --repo PATH")


def load_gram(path):
    return matrix(ZZ, [
        [ZZ(v) for v in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ])


def reflect(row, gram, root):
    row = vector(ZZ, list(row))
    root = vector(ZZ, list(root))
    assert root * gram * root == -2
    return row + (row * gram * root) * root


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
        candidate = vector(
            ZZ,
            [(i + 1) ** 2 + shift * (i + 1) + 1 for i in range(gram.nrows())],
        )
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
    lll = matrix(ZZ, pari((scale * H).change_ring(ZZ)).qflllgram())
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


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo", type=Path)
parser.add_argument("--output", type=Path)
args = parser.parse_args()

ROOT = locate_repo(args.repo)
LOCAL = ROOT / "artifacts/local/elkies-k3"
FRAME = ROOT / "elkies-k3/data/fibrations/kumar_e7e8_mw2_frame_3.txt"
TARGET = LOCAL / "q8-target-component-nef.json"
BRANCH = LOCAL / "q8-d13-branch-anchor.json"
OLDMW = LOCAL / "q8-old-q6-mw-images.json"
OUTPUT = (
    args.output.resolve()
    if args.output
    else LOCAL / "q8-explicit-old-curves.json"
)

for path in (FRAME, TARGET, BRANCH, OLDMW):
    if not path.exists():
        raise SystemExit(f"Missing prerequisite: {path}")

frame = load_gram(FRAME)
target = json.loads(TARGET.read_text())
branch = json.loads(BRANCH.read_text())
oldmw = json.loads(OLDMW.read_text())

assert target["status"] == "PASS_EXACT_Q6_CHILD_Q8_PHYSICAL_ROOT_TARGET"
assert branch["status"] == "PASS_EXACT_D13_BRANCH_ANCHOR"
assert oldmw["status"] == "PASS_EXACT_Q8_OLD_Q6_MW_IMAGES"

U2 = matrix(ZZ, ((0,1),(1,0)))
ns = block_diagonal_matrix(U2, -frame)
F0 = vector(ZZ, [1,0] + [0]*17)
O0 = vector(ZZ, [-1,1] + [0]*17)

source_simple = tuple(
    vector(ZZ, [0,0] + [ZZ(i == node) for i in range(17)])
    for node in range(15)
)
E7_simple = matrix(ZZ, [list(source_simple[i]) for i in range(7)])
E7_cartan = -(E7_simple * ns * E7_simple.transpose())
assert E7_cartan.det() == 2

e7_highest = highest_root(E7_cartan)
E7_affine = F0 - e7_highest * E7_simple
assert E7_affine * ns * E7_affine == -2

# Equation-level q6 curves are obtained by inverse Weyl transport, exactly as
# derive_h92_q6_child_q8_ivstar_orientation.sage does for old E7_7.
def inverse_q6_weyl(curve):
    result = vector(ZZ, list(curve))
    for node1 in reversed(Q6_REFLECTIONS):
        result = reflect(result, ns, source_simple[node1 - 1])
    return result

curves = [("old_zero", inverse_q6_weyl(O0))]
for i in range(7):
    curves.append((f"old_E7_{i+1}", inverse_q6_weyl(source_simple[i])))
curves.append(("old_E7_affine", inverse_q6_weyl(E7_affine)))

for name, curve in curves:
    assert curve * ns * curve == -2

# Anchored q8 frame.
F8 = vector(ZZ, target["selected_q8"]["source_h3_ns_vector"])
E6 = matrix(
    ZZ, target["selected_q8"]["E6"]["simple_root_vectors_in_source_h3_ns"]
)
E8 = matrix(
    ZZ, target["selected_q8"]["E8"]["simple_root_vectors_in_source_h3_ns"]
)
O8 = vector(ZZ, E8.row(0))
G1_curve = vector(ZZ, E6.row(0))
assert O8 * ns * F8 == 1
assert G1_curve * ns * F8 == 1

child8, Bzero = child_frame_with_zero(ns, F8, O8)
A13, adapted, H = d13_root_adaptation(child8)
Badapt = block_diagonal_matrix(identity_matrix(ZZ,2), A13) * Bzero
Gadapt = block_diagonal_matrix(U2, -adapted)
assert Badapt * ns * Badapt.transpose() == Gadapt

def coordinates(curve):
    value = vector(QQ, curve) * Badapt.inverse()
    assert all(v in ZZ for v in value)
    return vector(ZZ, value)

G1_coords = coordinates(G1_curve)
G1 = vector(ZZ, G1_coords[-4:])
assert G1 == vector(ZZ, branch["iv_sections"]["lattice"]["E6_1_mw"])
assert G1 == vector(ZZ, (1,0,0,0))

# Reconfirm old index-3 lattice and derive a quotient functional mod 3.
A = vector(ZZ, oldmw["old_q6_generators"][0]["mw_coordinates"])
B = vector(ZZ, oldmw["old_q6_generators"][1]["mw_coordinates"])
C = vector(ZZ, oldmw["old_q6_generators"][2]["mw_coordinates"])
L = matrix(ZZ, [list(G1), list(A), list(B), list(C)])
assert abs(L.det()) == 3
assert 2*B - A - G1 == vector(ZZ, (0,0,3,0))

# Find a nonzero linear functional ell on (Z/3)^4 annihilating L.
Fp = GF(3)
Lm = matrix(Fp, L)
kernel = Lm.right_kernel()
assert kernel.dimension() == 1
ell = vector(Fp, kernel.basis()[0])
assert ell and all(vector(Fp,row) * ell == 0 for row in L.rows())

def residue(z):
    return int(Fp(vector(Fp, z) * ell))

records = []
section_candidates = []
for name, curve in curves:
    c = coordinates(curve)
    degree = int(c[1])
    z = vector(ZZ, c[-4:])
    height = QQ(z * H * z)
    Oint = int(curve * ns * O8)
    res = residue(z)

    # If this is a genuine q8 section, Shioda's local correction is determined.
    correction = None
    if degree == 1:
        correction = QQ(4 + 2*Oint) - height
        assert correction in (QQ(0), QQ(1), QQ(13)/4)
        section_candidates.append((name, z, height, Oint, correction, res))

    record = {
        "curve": name,
        "q8_degree": degree,
        "q8_zero_intersection": Oint,
        "mw_coordinates": list(map(int, z)),
        "height": str(height),
        "index3_residue": res,
        "section_local_correction": None if correction is None else str(correction),
        "source_h3_ns": list(map(int, curve)),
    }
    records.append(record)

    print(
        "Q8OLDCURVE|"
        f"curve={name}|q8_degree={degree}|O8={Oint}|"
        f"mw={','.join(map(str,z))}|height={height}|"
        f"mod3={res}|"
        f"section_correction={'NA' if correction is None else correction}",
        flush=True,
    )

nonzero_sections = [
    item for item in section_candidates if item[-1] != 0
]
print(
    "Q8OLDCURVE_SUMMARY|"
    f"q8_sections={len(section_candidates)}|"
    f"missing_coset_sections={len(nonzero_sections)}|"
    "names={}".format(
        ",".join(item[0] for item in nonzero_sections) or "none"
    ),
    flush=True,
)

# Exact lattice relations implied by the previous audit.
G3 = vector(ZZ, (0,0,1,0))
G2 = vector(ZZ, (0,1,0,0))
G4 = vector(ZZ, (0,0,0,1))
assert 3*G3 == 2*B - A - G1
assert G2 == 3*G1 + 7*G3 - B
assert G4 == C - 22*G1 + G2 - 45*G3

# The q24 target for the branch-point zero was certified previously.
q24 = vector(ZZ, (2,-1,-1,1))
assert q24 == C - 20*G1 - 46*G3
assert q24 * H * q24 == 52

payload = {
    "schema": "elkies-k3.h92-q8-explicit-old-curves.v1",
    "status": "PASS_EXACT_Q8_EXPLICIT_OLD_CURVE_PROFILE",
    "anchor": {
        "zero": "II*_E8_1",
        "G1": list(map(int, G1)),
        "height_gram": [[str(v) for v in row] for row in H.rows()],
    },
    "old_index3_lattice": {
        "A": list(map(int, A)),
        "B": list(map(int, B)),
        "C": list(map(int, C)),
        "determinant": int(L.det()),
        "quotient_functional_mod3": [int(v) for v in ell],
        "relations": {
            "3G3": "2B-A-G1",
            "G2": "3G1+7G3-B",
            "G4": "C-22G1+G2-45G3",
            "q24_branch_zero": "C-20G1-46G3",
        },
    },
    "curves": records,
    "q8_section_candidates": [
        {
            "curve": name,
            "mw_coordinates": list(map(int,z)),
            "height": str(height),
            "P_dot_O": Oint,
            "local_correction": str(correction),
            "index3_residue": res,
        }
        for name,z,height,Oint,correction,res in section_candidates
    ],
    "missing_coset_q8_sections": [
        {
            "curve": name,
            "mw_coordinates": list(map(int,z)),
            "height": str(height),
            "P_dot_O": Oint,
            "local_correction": str(correction),
            "index3_residue": res,
        }
        for name,z,height,Oint,correction,res in nonzero_sections
    ],
    "boundary": (
        "These are exact NS classes of already-explicit old curves. "
        "If a missing-coset candidate has q8 degree one, its rational "
        "Weierstrass coordinates can be obtained from the existing q6 "
        "section artifact and the anchored quartic map; this script itself "
        "does not perform that coordinate substitution."
    ),
}

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print("Q8OLDCURVE_RESULT|status=PASS_EXACT_Q8_EXPLICIT_OLD_CURVE_PROFILE")
print(f"OUTPUT|{OUTPUT}")
