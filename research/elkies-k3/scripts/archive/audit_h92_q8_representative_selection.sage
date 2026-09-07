#!/usr/bin/env sage -python
"""
Audit the H3 q8 representative selection.

Compare, in the SAME original H3 NS coordinates:

  A) dominant D13/MW4 orbit hit used by derive_h92_q8_generic_rr_ambient.sage;
  B) q8.nef_representative from classify_h3_q6_child_q8_orbits.sage;
  C) the repository's current degree-18 "source-nef" generic-ambient class;
  D) the later experimental degree-16 horizontal/fibre chamber reduction.

This is a coordinate identity audit only.  It does not use any local RR
compiler assumptions.

Run:
  sage -python ~/Downloads/audit_h92_q8_representative_selection.sage
"""

import argparse
import json
from pathlib import Path

from sage.all import ZZ, block_diagonal_matrix, matrix, vector


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
    for candidate in candidates:
        if (
            (candidate / "elkies-k3" / "scripts").is_dir()
            and (candidate / "artifacts" / "generated-results").is_dir()
        ):
            return candidate.resolve()
    raise SystemExit("Could not locate jacobian-research; pass --repo PATH")


def load_gram(path):
    return matrix(ZZ, [
        [ZZ(v) for v in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ])


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo", type=Path)
args = parser.parse_args()

ROOT = locate_repo(args.repo)
GEN = ROOT / "artifacts" / "generated-results"
FRAME = ROOT / "elkies-k3" / "data" / "fibrations" / "kumar_e7e8_mw2_frame_3.txt"
ORBITS = GEN / "elkies-k3-h3-q6-q8-orbits.json"
AMBIENT = GEN / "elkies-k3-h92-q8-generic-rr-ambient.json"
CHAMBER = GEN / "zz-h92-q8-complete-chamber-reduction.json"

for path in (FRAME, ORBITS, AMBIENT):
    if not path.exists():
        raise SystemExit(f"Missing prerequisite: {path}")

frame = load_gram(FRAME)
NS = block_diagonal_matrix(matrix(ZZ, ((0,1),(1,0))), -frame)
F = vector(ZZ, [1,0] + [0]*17)
simple = tuple(
    vector(ZZ, [0,0] + [ZZ(i == j) for i in range(17)])
    for j in range(15)
)
highest_e7 = (2,2,3,4,3,2,1)
highest_e8 = (2,3,4,6,5,4,3,2)
zero = vector(ZZ, [0]*19)
affine_e7 = F - sum((c*simple[i] for i,c in enumerate(highest_e7)), zero)
affine_e8 = F - sum((c*simple[7+i] for i,c in enumerate(highest_e8)), zero)

orbits = json.loads(ORBITS.read_text())
ambient = json.loads(AMBIENT.read_text())
assert orbits["status"] == "PASS_H3_Q6_CHILD_Q8_WEYL_CLASSIFICATION"
assert ambient["status"] == "PASS_EXACT_Q8_GENERIC_RR_AMBIENT"

dominant_record = next(
    hit for hit in orbits["q8"]["d13_mw4_hits"]
    if hit["mw_projection"] == [0,-2,0]
)
nef_record = orbits["q8"]["nef_representative"]
assert nef_record["mw_projection"] == [0,-2,0]

classes = {
    "dominant_raw": vector(ZZ, dominant_record["fiber_source_h3_ns"]),
    "classifier_nef": vector(ZZ, nef_record["fiber_source_h3_ns"]),
    "ambient_degree18": vector(ZZ, ambient["source_q8_lattice_class"]),
}
if CHAMBER.exists():
    chamber = json.loads(CHAMBER.read_text())
    classes["experimental_degree16"] = vector(ZZ, chamber["final_class"])


def invariants(v):
    return {
        "square": int(v*NS*v),
        "old_F_degree": int(v*NS*F),
        "simple_pairings": [int(v*NS*r) for r in simple],
        "affine_pairings": [
            int(v*NS*affine_e7),
            int(v*NS*affine_e8),
        ],
        "vector": [int(x) for x in v],
    }


for name,v in classes.items():
    inv = invariants(v)
    print(
        "Q8REP|"
        f"name={name}|square={inv['square']}|oldF={inv['old_F_degree']}|"
        f"simple={','.join(map(str,inv['simple_pairings']))}|"
        f"affine={','.join(map(str,inv['affine_pairings']))}|"
        f"vector={','.join(map(str,inv['vector']))}",
        flush=True,
    )

names = list(classes)
for i in range(len(names)):
    for j in range(i+1,len(names)):
        a,b = names[i], names[j]
        diff = classes[a]-classes[b]
        print(
            "Q8REP_DIFF|"
            f"a={a}|b={b}|equal={int(diff==0)}|"
            f"diff_square={int(diff*NS*diff)}|diff_F={int(diff*NS*F)}|"
            f"pairing={int(classes[a]*NS*classes[b])}|"
            f"diff={','.join(map(str,diff))}",
            flush=True,
        )

if "experimental_degree16" in classes:
    print(
        "Q8REP_KEYCHECK|"
        f"experimental_eq_classifier_nef={int(classes['experimental_degree16']==classes['classifier_nef'])}|"
        f"ambient18_eq_classifier_nef={int(classes['ambient_degree18']==classes['classifier_nef'])}|"
        f"dominant_eq_classifier_nef={int(classes['dominant_raw']==classes['classifier_nef'])}",
        flush=True,
    )

# Reproduce exactly the finite-root reduction performed by the generic ambient
# and see where each source representative lands.
def reduce_finite_roots(v):
    v = vector(ZZ, v)
    steps = []
    for _ in range(1000):
        neg = [
            (i, int(v*NS*r))
            for i,r in enumerate(simple)
            if v*NS*r < 0
        ]
        if not neg:
            return v, steps
        i,p = neg[0]
        v += p*simple[i]
        steps.append((i+1,p))
    raise RuntimeError("finite-root reduction did not terminate")

for name in ("dominant_raw","classifier_nef"):
    reduced, steps = reduce_finite_roots(classes[name])
    print(
        "Q8REP_FINITEREDUCE|"
        f"name={name}|steps={len(steps)}|oldF={int(reduced*NS*F)}|"
        f"equals_ambient18={int(reduced==classes['ambient_degree18'])}|"
        f"equals_classifier_nef={int(reduced==classes['classifier_nef'])}|"
        f"vector={','.join(map(str,reduced))}",
        flush=True,
    )
