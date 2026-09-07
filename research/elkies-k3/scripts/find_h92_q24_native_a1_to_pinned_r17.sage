#!/usr/bin/env sage -python
"""
Search successive q6 rootless children of the q24-native A1/MW16 frame until
one is integrally isometric to the pinned recovered R17 lattice.

A rootless determinant-948 child is NOT sufficient: several integral classes
may occur. Each rootless hit is therefore tested with PARI qfisom against
data/lattice/rank17_gram.txt. Nonmatching rootless classes are recorded and
the deterministic stream resumes after their orbit index.
"""

import json
import subprocess
from pathlib import Path

from sage.all import ZZ, block_diagonal_matrix, identity_matrix, matrix, pari, vector

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "elkies-k3" / "scripts"
LOCAL = ROOT / "artifacts" / "local" / "elkies-k3" / "q24-native-suffix"

SEARCH = SCRIPTS / "search_root_adapted_weyl_neighbors_targeted.sage"
SUFFIX = LOCAL / "q24-native-d12-to-a1.json"
PINNED_PATH = ROOT / "elkies-k3" / "data" / "lattice" / "rank17_gram.txt"

CACHE = LOCAL / "step10-rootless-q6-mw-vectors-cap10000.json"
OUT = LOCAL / "q24-native-a1-to-pinned-r17.json"
CANDIDATES = LOCAL / "q24-native-a1-rootless-candidates.json"

U2 = matrix(ZZ, ((0,1),(1,0)))


def load_gram(path):
    return matrix(ZZ, [
        [ZZ(x) for x in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ])


def rows(M):
    return [[int(x) for x in row] for row in M.rows()]


def run(cmd):
    print("+", " ".join(map(str, cmd)), flush=True)
    subprocess.run([str(x) for x in cmd], cwd=ROOT, check=True)


def exact_isometry_child_to_pinned(child, pinned):
    # PARI: qfisom(G,H) returns S with G = S^t H S.
    raw = pari(child).qfisom(pari(pinned))
    if str(raw) != "0":
        S = matrix(ZZ, raw)
        assert child == S.transpose() * pinned * S
        Q = S.inverse().transpose()
        if Q.change_ring(ZZ) == Q:
            Q = Q.change_ring(ZZ)
            assert Q * child * Q.transpose() == pinned
            assert abs(ZZ(Q.det())) == 1
            return Q, "qfisom(child,pinned)"

    # Symmetric fallback: pinned = S^t child S.
    raw = pari(pinned).qfisom(pari(child))
    if str(raw) != "0":
        S = matrix(ZZ, raw)
        assert pinned == S.transpose() * child * S
        Q = S.transpose()
        assert Q * child * Q.transpose() == pinned
        assert abs(ZZ(Q.det())) == 1
        return Q, "qfisom(pinned,child)"

    return None, None


for path in (SEARCH, SUFFIX, PINNED_PATH):
    if not path.exists():
        raise SystemExit(f"missing prerequisite: {path}")

suffix = json.loads(SUFFIX.read_text())
assert suffix["status"] == "PASS_Q24_NATIVE_D12_TO_A1"
a1_path = ROOT / suffix["final_frame"]
if not a1_path.exists():
    raise SystemExit(f"missing native A1 frame: {a1_path}")

A1 = load_gram(a1_path)
pinned = load_gram(PINNED_PATH)
assert A1.dimensions() == pinned.dimensions() == (17,17)
assert A1.det() == pinned.det() == 948

qf = pari(A1).qfminim(2)
assert (matrix(ZZ,qf[2]).rank(), ZZ(qf[0])) == (1,2)

print(
    "Q24NATIVE_R17_PARENT|"
    f"frame={a1_path.relative_to(ROOT)}|root_data=1,2,2|MW=16|status=PASS",
    flush=True,
)

skip = 0
chunk = 5000
max_stream = 100000
attempt = 0
tested_rootless = []
found = None

while skip < max_stream:
    attempt += 1
    search_out = LOCAL / f"step10-rootless-q6-search-skip{skip:06d}.json"
    frames = LOCAL / f"step10-rootless-q6-frames-skip{skip:06d}"

    cmd = [
        "sage","-python",SEARCH,
        "--frame",a1_path,
        "--root-rank","1",
        "--q","6",
        "--degree","2",
        "--adapt-mw-at-least","17",
        "--rank-growth-only",
        "--stop-after-first-growth",
        "--stream-first-growth",
        "--target-root-rank","0",
        "--target-root-count","0",
        "--target-root-determinant","1",
        "--mw-vector-cap","10000",
        "--mw-vectors-cache",CACHE,
        "--stream-skip",str(skip),
        "--stream-limit",str(chunk),
        "--output",search_out,
        "--frames-dir",frames,
        "--stream-progress-every","500",
    ]
    run(cmd)

    data = json.loads(search_out.read_text())
    assert data["status"] == "PASS_ROOT_ADAPTED_WEYL_NEIGHBORS"
    summaries = data.get("summaries", [])
    assert len(summaries) == 1
    summary = summaries[0]

    hits = [
        r for r in data.get("neighbors", [])
        if int(r["q"]) == 6
        and int(r["old_fiber_degree"]) == 2
        and tuple(r["child_root_data"]) == (0,0,1)
        and r["child_ade"] == "rootless"
        and int(r["child_mw_rank"]) == 17
    ]

    if not hits:
        tested = int(summary.get("stream_tested", 0))
        next_skip = int(summary.get("screened_orbits", skip + tested))
        print(
            "Q24NATIVE_R17_CHUNK|"
            f"skip={skip}|tested={tested}|next_skip={next_skip}|"
            "rootless_hits=0|status=CONTINUE",
            flush=True,
        )
        if next_skip <= skip:
            raise SystemExit("stream made no progress")
        skip = next_skip
        continue

    assert len(hits) == 1
    rec = hits[0]
    orbit = int(rec["orbit_index"])
    child = matrix(ZZ, rec["child_root_adapted_frame"])
    assert child.det() == 948
    assert ZZ(pari(child).qfminim(2)[0]) == 0

    Q, method = exact_isometry_child_to_pinned(child, pinned)
    is_pinned = Q is not None

    tested_rootless.append({
        "orbit_index": orbit,
        "mw_projection": rec["mw_projection"],
        "dominant_labels": rec["dominant_labels"],
        "witness": rec["witness"],
        "pinned_r17_isometric": bool(is_pinned),
        "search_artifact": str(search_out.relative_to(ROOT)),
    })
    CANDIDATES.write_text(json.dumps({
        "schema":"elkies-k3.h3-q24-native-a1-rootless-candidates.v1",
        "status":"SEARCHING" if not is_pinned else "PASS_PINNED_FOUND",
        "parent":str(a1_path.relative_to(ROOT)),
        "tested_rootless":tested_rootless,
    }, indent=2, sort_keys=True) + "\n")

    print(
        "Q24NATIVE_R17_CANDIDATE|"
        f"orbit={orbit}|rootless=1|det=948|"
        f"pinned_isometric={int(is_pinned)}|"
        f"status={'PASS_PINNED_R17' if is_pinned else 'SKIP_OTHER_ROOTLESS_CLASS'}",
        flush=True,
    )

    if is_pinned:
        found = (rec, child, Q, method, search_out)
        break

    # Continue immediately after this deterministic orbit.
    skip = orbit

if found is None:
    raise SystemExit(
        f"No pinned-R17 rootless child found in first {max_stream} deterministic "
        "q6 witnesses. Increase max_stream / MW-vector cap."
    )

rec, child, Q, method, search_out = found
orbit = int(rec["orbit_index"])

raw_child = matrix(ZZ, rec["child_frame"])
adapt = matrix(ZZ, rec["child_root_adapted_basis"])
neighbor = matrix(ZZ, rec["neighbor_basis"])
assert adapt * raw_child * adapt.transpose() == child

transition = block_diagonal_matrix(identity_matrix(ZZ,2), adapt) * neighbor
NS_A1 = block_diagonal_matrix(U2, -A1)
NS_child = block_diagonal_matrix(U2, -child)
assert transition * NS_A1 * transition.transpose() == NS_child
assert abs(ZZ(transition.det())) == 1

endpoint = block_diagonal_matrix(identity_matrix(ZZ,2), Q)
NS_pinned = block_diagonal_matrix(U2, -pinned)
assert endpoint * NS_child * endpoint.transpose() == NS_pinned

total = endpoint * transition
assert total * NS_A1 * total.transpose() == NS_pinned
assert abs(ZZ(total.det())) == 1

payload = {
    "schema":"elkies-k3.h3-q24-native-a1-to-pinned-r17.v2",
    "status":"PASS_Q24_NATIVE_A1_TO_PINNED_R17",
    "parent":{
        "frame":str(a1_path.relative_to(ROOT)),
        "root_data":[1,2,2],
        "ade":"A1",
        "mw_rank":16,
    },
    "neighbor":{
        "q":6,
        "factor_order":[3,2],
        "old_fibre_degree":2,
        "discovered_orbit_index":orbit,
        "mw_projection":rec["mw_projection"],
        "dominant_labels":rec["dominant_labels"],
        "witness":rec["witness"],
        "fiber":rec["fiber"],
        "child_root_data":[0,0,1],
        "child_ade":"rootless",
        "child_mw_rank":17,
        "search_is_exhaustive":False,
        "mw_vector_cap":10000,
        "search_artifact":str(search_out.relative_to(ROOT)),
    },
    "rootless_candidates_tested":tested_rootless,
    "rootless_child_frame":rows(child),
    "pinned_rank17_frame":str(PINNED_PATH.relative_to(ROOT)),
    "isometry_method":method,
    "rootless_child_to_pinned_r17_isometry":rows(Q),
    "native_a1_to_rootless_ns_transport":rows(transition),
    "native_a1_to_pinned_r17_ns_transport":rows(total),
    "proof_boundary":(
        "Exact positive lattice/NS certificate from the q24-native A1 frame "
        "to pinned recovered R17. The bounded deterministic search only "
        "discovers candidate witnesses; the chosen q6 neighbor, rootlessness, "
        "integral endpoint isometry, and determinant-one NS transport are "
        "verified exactly. This does not execute the characteristic-zero "
        "A1->R17 Weierstrass pencil."
    ),
}

OUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{OUT}",flush=True)
print(
    "Q24NATIVE_R17_RESULT|"
    f"orbit={orbit}|rootless_candidates_tested={len(tested_rootless)}|"
    "final=rootless/MW17|pinned_R17=1|"
    "status=PASS_Q24_NATIVE_A1_TO_PINNED_R17",
    flush=True,
)
