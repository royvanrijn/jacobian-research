#!/usr/bin/env sage -python
"""
Profile the original H3 rational sections -P1 and P2 in the anchored q8/D13
fibration.

This is a cheap search for the missing fourth MW direction.  The original H3
surface already has two explicit section classes:
    -P1 = [5,1,-2,-3,-4,-6,-5,-4,-3,0,...,0,1,0]
     P2 = [22,1,0,...,0,1]
in the canonical H3 NS frame.

We compute their q8 multisection degrees and anchored D13 MW coordinates with
zero II*_E8_1.  We also profile source O and a few degree-zero combinations at
the MW-coordinate level.

Run:
  sage -python ~/Downloads/audit_h92_q8_original_h3_sections.sage
"""

import argparse
import json
from pathlib import Path

from sage.all import QQ, ZZ, block_diagonal_matrix, identity_matrix, lcm, matrix, pari, vector


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
    mate = zero + fibre
    orth = matrix(
        ZZ, [list(fibre * ns), list(mate * ns)]
    ).right_kernel_matrix()
    basis = matrix(
        ZZ, [list(fibre), list(mate)] + [list(row) for row in orth.rows()]
    )
    assert abs(basis.det()) == 1
    child = -(orth * ns * orth.transpose())
    U = matrix(ZZ, ((0,1),(1,0)))
    assert basis * ns * basis.transpose() == block_diagonal_matrix(U, -child)
    return child, basis


def roots_and_data(gram):
    result = pari(gram).qfminim(2)
    count = ZZ(result[0])
    half = [vector(ZZ, c) for c in matrix(ZZ, result[2]).columns()]
    roots = tuple(half + [-r for r in half])
    rb = matrix(ZZ, [list(r) for r in roots]).row_module().basis_matrix()
    rg = rb * gram * rb.transpose()
    return roots, rb, (rb.rank(), count, abs(ZZ(rg.det())))


def deterministic_simple_roots(gram):
    roots, unused, data = roots_and_data(gram)
    rank = data[0]
    regular = None
    for shift in range(1,1000):
        candidate = vector(ZZ, [
            (i+1)**2 + shift*(i+1) + 1 for i in range(gram.nrows())
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
    assert invariants == (13,312,4), invariants
    simple, cartan = deterministic_simple_roots(child)
    assert cartan.det() == 4

    smith, left, right = root_basis.smith_form()
    assert smith == left * root_basis * right
    completion = right.inverse()
    initial = simple.stack(completion[13:])
    assert abs(initial.det()) == 1

    adapted = initial * child * initial.transpose()
    root = adapted[:13,:13]
    coupling = adapted[:13,13:]
    tail = adapted[13:,13:]
    H = tail - coupling.transpose()*root.inverse()*coupling
    scale = ZZ(1)
    for v in H.list():
        scale = lcm(scale, ZZ(QQ(v).denominator()))
    lll = matrix(ZZ, pari((scale*H).change_ring(ZZ)).qflllgram())
    change = block_diagonal_matrix(identity_matrix(ZZ,13), lll.transpose())
    basis = change * initial
    adapted = basis * child * basis.transpose()

    root = adapted[:13,:13]
    coupling = adapted[:13,13:]
    tail = adapted[13:,13:]
    H = tail - coupling.transpose()*root.inverse()*coupling
    assert H.det() == 237
    return basis, adapted, H


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo", type=Path)
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
    else LOCAL / "q8-original-h3-sections.json"
)

for path in (FRAME, TARGET, BRANCH, G3FILE):
    if not path.exists():
        raise SystemExit(f"Missing prerequisite: {path}")

frame = load_gram(FRAME)
target = json.loads(TARGET.read_text())
branch = json.loads(BRANCH.read_text())
g3data = json.loads(G3FILE.read_text())
assert branch["status"] == "PASS_EXACT_D13_BRANCH_ANCHOR"
assert g3data["status"] == "PASS_EXACT_D13_G3_FROM_E77_BISECTION"

U2 = matrix(ZZ, ((0,1),(1,0)))
ns = block_diagonal_matrix(U2, -frame)
F0 = vector(ZZ, [1,0] + [0]*17)
O0 = vector(ZZ, [-1,1] + [0]*17)

minus_p1 = vector(ZZ, [5,1] + [
    -2,-3,-4,-6,-5,-4,-3,
    0,0,0,0,0,0,0,0,
    1,0
])
p2 = vector(ZZ, [22,1] + [0]*16 + [1])

for name, P in (("source_O",O0),("minus_P1",minus_p1),("P2",p2)):
    assert len(P) == 19, (name,len(P))
    assert P * ns * P == -2, (name,P*ns*P)
    assert P * ns * F0 == 1, (name,P*ns*F0)

F8 = vector(ZZ, target["selected_q8"]["source_h3_ns_vector"])
E6 = matrix(ZZ, target["selected_q8"]["E6"]["simple_root_vectors_in_source_h3_ns"])
E8 = matrix(ZZ, target["selected_q8"]["E8"]["simple_root_vectors_in_source_h3_ns"])
O8 = vector(ZZ, E8.row(0))
G1curve = vector(ZZ, E6.row(0))

child8, Bzero = child_frame_with_zero(ns,F8,O8)
A13, adapted, H = d13_root_adaptation(child8)
Badapt = block_diagonal_matrix(identity_matrix(ZZ,2),A13) * Bzero
Gadapt = block_diagonal_matrix(U2,-adapted)
assert Badapt * ns * Badapt.transpose() == Gadapt

def coords(curve):
    c = vector(QQ,curve) * Badapt.inverse()
    assert all(v in ZZ for v in c)
    return vector(ZZ,c)

G1 = vector(ZZ, coords(G1curve)[-4:])
assert G1 == vector(ZZ,(1,0,0,0))
G3 = vector(ZZ,(0,0,1,0))
assert G3 * H * G3 == QQ(11)/4

records=[]
for name,P in (("source_O",O0),("minus_P1",minus_p1),("P2",p2)):
    c=coords(P)
    d=ZZ(c[1])
    z=vector(ZZ,c[-4:])
    height=QQ(z*H*z)
    rec={
        "curve":name,
        "q8_degree":int(d),
        "q8_zero_intersection":int(P*ns*O8),
        "mw_coordinates":list(map(int,z)),
        "height":str(height),
        "source_h3_ns":list(map(int,P)),
    }
    records.append(rec)
    print(
        "Q8H3SECTION|"
        f"curve={name}|q8_degree={d}|O8={P*ns*O8}|"
        f"mw={','.join(map(str,z))}|height={height}",
        flush=True,
    )

# Search small integer combinations of the three KNOWN anchored vectors
# G1, G3, and each source-section AJ vector for a pure fourth-coordinate
# direction.  This is only lattice arithmetic; no section effectiveness claim.
candidates=[]
for rec in records:
    z=vector(ZZ,rec["mw_coordinates"])
    for a in range(-12,13):
        for b in range(-12,13):
            w=z-a*G1-b*G3
            if w[0]==0 and w[1]==0 and w[2]==0 and w[3]!=0:
                candidates.append((abs(w[3]),rec["curve"],a,b,tuple(w)))
candidates.sort()

if candidates:
    best=candidates[0]
    print(
        "Q8H3SECTION_FOURTH|"
        f"curve={best[1]}|subtract={best[2]}*G1+{best[3]}*G3|"
        f"residual={','.join(map(str,best[4]))}|status=FOUND",
        flush=True,
    )
else:
    print("Q8H3SECTION_FOURTH|status=NONE_SMALL",flush=True)

payload={
    "schema":"elkies-k3.h92-q8-original-h3-sections.v1",
    "status":"PASS_EXACT_Q8_ORIGINAL_H3_SECTION_PROFILE",
    "anchor":{"G1":[1,0,0,0],"G3":[0,0,1,0],
              "height_gram":[[str(v) for v in row] for row in H.rows()]},
    "sections":records,
    "pure_fourth_candidates":[
        {"abs_coefficient":int(absw),"curve":curve,
         "subtract_G1":int(a),"subtract_G3":int(b),
         "residual":list(map(int,w))}
        for absw,curve,a,b,w in candidates
    ],
    "boundary":"Exact NS profiles only; degree>1 curves require Abel-Jacobi reduction for equation coordinates."
}
OUTPUT.parent.mkdir(parents=True,exist_ok=True)
OUTPUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print("Q8H3SECTION_RESULT|status=PASS_EXACT_Q8_ORIGINAL_H3_SECTION_PROFILE")
print(f"OUTPUT|{OUTPUT}")
