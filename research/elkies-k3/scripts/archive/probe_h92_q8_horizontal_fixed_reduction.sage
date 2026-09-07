#!/usr/bin/env sage -python
"""
Reduce the H92 q8 source isotropic class against horizontal (-2)-curve
obstructions of the old E7+E8 fibration.

This continues the bounded exact nef scan after the discovery of five
degree-two roots with D.C=-1.  First, every discovered root is reduced against
all 17 old fibre components (finite and affine) to identify duplicate
vertical representatives of the same horizontal bisection.  Then the q8 class
is Weyl-reflected across a deterministic negative horizontal root and the
bounded search is repeated.

The search at each round is exact for every (-2)-class C with
1 <= C.F <= max_degree that can satisfy current_D.C <= 0.

Run:
  sage -python ~/Downloads/probe_h92_q8_horizontal_fixed_reduction.sage

Optional:
  --repo /path/to/jacobian-research
  --max-degree 24
  --max-rounds 12
"""

import argparse
import json
from pathlib import Path

from sage.all import ZZ, block_diagonal_matrix, matrix, pari, vector


def locate_repo(explicit=None):
    candidates = []
    if explicit:
        candidates.append(Path(explicit).expanduser())
    cwd = Path.cwd().resolve()
    candidates += [cwd, *cwd.parents]
    home = Path.home()
    candidates += [
        home / "jacobian-research",
        home / "src" / "jacobian-research",
        home / "git" / "jacobian-research",
        home / "projects" / "jacobian-research",
        home / "Documents" / "jacobian-research",
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
            and (candidate / "elkies-k3" / "data" / "fibrations").is_dir()
            and (candidate / "artifacts" / "generated-results").is_dir()
        ):
            return candidate
    raise SystemExit("Could not locate jacobian-research; pass --repo PATH")


def load_gram(path):
    return matrix(ZZ, [
        [ZZ(value) for value in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ])


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo", type=Path)
parser.add_argument("--max-degree", type=int, default=24)
parser.add_argument("--max-rounds", type=int, default=12)
args = parser.parse_args()

ROOT = locate_repo(args.repo)
GEN = ROOT / "artifacts" / "generated-results"
FRAME = ROOT / "elkies-k3" / "data" / "fibrations" / "kumar_e7e8_mw2_frame_3.txt"
AMBIENT = GEN / "elkies-k3-h92-q8-generic-rr-ambient.json"
AUDIT = GEN / "zz-h92-q8-source-nef-multisection-audit.json"
OUTPUT = GEN / "zz-h92-q8-horizontal-fixed-reduction.json"

ambient = json.loads(AMBIENT.read_text())
audit = json.loads(AUDIT.read_text())
assert ambient["status"] == "PASS_EXACT_Q8_GENERIC_RR_AMBIENT"
assert audit["status"] == "FOUND_NEGATIVE_EFFECTIVE_ROOTS"

G = load_gram(FRAME)
NS = block_diagonal_matrix(matrix(ZZ, ((0, 1), (1, 0))), -G)
F = vector(ZZ, [1, 0] + [0] * 17)

simple = tuple(
    vector(ZZ, [0, 0] + [ZZ(index == node) for index in range(17)])
    for node in range(15)
)
highest_e7 = (2, 2, 3, 4, 3, 2, 1)
highest_e8 = (2, 3, 4, 6, 5, 4, 3, 2)
zero19 = vector(ZZ, [0] * 19)
affine_e7 = F - sum(
    (c * simple[i] for i, c in enumerate(highest_e7)), zero19
)
affine_e8 = F - sum(
    (c * simple[7+i] for i, c in enumerate(highest_e8)), zero19
)
fiber_components = simple + (affine_e7, affine_e8)
fiber_names = tuple(f"E7_{i}" for i in range(1,8)) + tuple(
    f"E8_{i}" for i in range(1,9)
) + ("E7_affine", "E8_affine")

for C in fiber_components:
    assert C * NS * C == -2
    assert C * NS * F == 0


def fiber_pairings(C):
    return tuple(int(C * NS * E) for E in fiber_components)


def reduce_root_to_fiber_chamber(C):
    C = vector(ZZ, C)
    steps = []
    # At positive old-fibre degree the affine Weyl alcove is finite.
    for _ in range(10000):
        pairings = fiber_pairings(C)
        negative = [(i, p) for i, p in enumerate(pairings) if p < 0]
        if not negative:
            assert C * NS * C == -2 and C * NS * F > 0
            return C, steps
        i, p = min(negative, key=lambda item: (item[0], item[1]))
        C = C + p * fiber_components[i]
        steps.append((fiber_names[i], int(p)))
        assert C * NS * C == -2
    raise RuntimeError("fiber chamber reduction did not terminate")


initial_hits = []
for n_text, entries in audit["negative_roots"].items():
    for entry in entries:
        C = vector(ZZ, entry["class"])
        reduced, steps = reduce_root_to_fiber_chamber(C)
        initial_hits.append({
            "input": list(map(int, C)),
            "input_D_pairing": int(vector(ZZ, ambient["source_q8_lattice_class"]) * NS * C),
            "reduced": list(map(int, reduced)),
            "steps": [[name, p] for name, p in steps],
            "reduced_fiber_pairings": list(fiber_pairings(reduced)),
        })

unique_initial = {}
D0 = vector(ZZ, ambient["source_q8_lattice_class"])
for rec in initial_hits:
    key = tuple(rec["reduced"])
    unique_initial.setdefault(key, {
        "class": list(key),
        "multiplicity": 0,
        "D_pairing": int(D0 * NS * vector(ZZ, key)),
        "old_fibre_degree": int(vector(ZZ, key) * NS * F),
        "fiber_pairings": list(fiber_pairings(vector(ZZ, key))),
    })
    unique_initial[key]["multiplicity"] += 1

print(
    "Q8HORIZONTAL_INITIAL|"
    f"raw_hits={len(initial_hits)}|unique_fiber_nef={len(unique_initial)}",
    flush=True,
)
for i, rec in enumerate(unique_initial.values()):
    print(
        "Q8HORIZONTAL_CLASS|"
        f"index={i}|multiplicity={rec['multiplicity']}|"
        f"degree={rec['old_fibre_degree']}|Dpair={rec['D_pairing']}|"
        f"fiber_pairings={','.join(map(str, rec['fiber_pairings']))}|"
        f"class={','.join(map(str, rec['class']))}",
        flush=True,
    )


def enumerate_nonpositive_roots(D, max_degree):
    """
    Exact enumeration of roots C with 1<=C.F<=max_degree and D.C<=0.

    For D=(d0,r,d), C=(c0,n,v):
      2*r*n*(D.C) = Q(r*v-n*d) - 2*r^2.
    So D.C<=0 implies Q(r*v-n*d)<=2*r^2.
    """
    r = ZZ(D * NS * F)
    assert r > 0 and D * NS * D == 0
    d = vector(ZZ, D[2:])
    result = []

    for n in range(1, int(max_degree) + 1):
        nd = n * d
        Gd = G * nd.column()

        A = matrix(ZZ, 18, 18)
        A[:17, :17] = r*r * G
        cross = -r * Gd
        for i in range(17):
            A[i, 17] = cross[i, 0]
            A[17, i] = cross[i, 0]
        A[17,17] = ZZ(nd * G * nd) + 1
        assert A.is_positive_definite()

        # Qshift <= 2r^2 corresponds to augmented norm <= 2r^2+1.
        raw = pari(A).qfminim(2*r*r + 1)
        reps = matrix(ZZ, raw[2]).transpose() if int(raw[0]) else matrix(ZZ, 0, 18)

        roots = {}
        for row in reps.rows():
            k = ZZ(row[17])
            if abs(k) != 1:
                continue
            rr = vector(ZZ, row)
            if k == -1:
                rr = -rr
            if rr[17] != 1:
                continue
            v = vector(ZZ, rr[:17])
            qv = ZZ(v * G * v)
            if (qv - 2) % (2*n):
                continue
            c0 = (qv - 2) // (2*n)
            C = vector(ZZ, [c0, n] + list(v))
            assert C * NS * C == -2 and C * NS * F == n
            pair = ZZ(D * NS * C)
            if pair > 0:
                continue
            reduced, steps = reduce_root_to_fiber_chamber(C)
            # Record the fibre-nef representative; it is the useful horizontal
            # curve class rather than a root plus vertical components.
            rpair = ZZ(D * NS * reduced)
            if rpair > 0:
                continue
            key = tuple(reduced)
            roots[key] = {
                "degree": int(reduced * NS * F),
                "D_pairing": int(rpair),
                "class": list(map(int, reduced)),
                "fiber_pairings": list(fiber_pairings(reduced)),
                "reduction_steps": [[name, p] for name, p in steps],
            }
        result.extend(roots.values())

    # dedupe across n (degree is invariant so normally unnecessary)
    dedup = {tuple(rec["class"]): rec for rec in result}
    return sorted(
        dedup.values(),
        key=lambda rec: (
            rec["D_pairing"],
            rec["degree"],
            tuple(rec["class"]),
        ),
    )


history = []
D = vector(ZZ, D0)
for round_index in range(int(args.max_rounds)):
    roots = enumerate_nonpositive_roots(D, args.max_degree)
    negative = [rec for rec in roots if rec["D_pairing"] < 0]
    zero = [rec for rec in roots if rec["D_pairing"] == 0]

    print(
        "Q8HORIZONTAL_ROUND|"
        f"round={round_index}|DF={D*NS*F}|D2={D*NS*D}|"
        f"negative={len(negative)}|zero={len(zero)}",
        flush=True,
    )

    record = {
        "round": round_index,
        "class_before": list(map(int, D)),
        "old_fibre_degree_before": int(D * NS * F),
        "negative_roots": negative,
        "zero_root_count": len(zero),
    }

    if not negative:
        record["action"] = "stop"
        history.append(record)
        break

    chosen = negative[0]
    C = vector(ZZ, chosen["class"])
    p = ZZ(D * NS * C)
    assert p < 0
    Dnext = D + p * C  # Weyl reflection in the (-2)-root C
    assert Dnext * NS * Dnext == 0

    record["action"] = "reflect"
    record["chosen_root"] = chosen
    record["reflection_coefficient"] = int(p)
    record["class_after"] = list(map(int, Dnext))
    record["old_fibre_degree_after"] = int(Dnext * NS * F)
    history.append(record)

    print(
        "Q8HORIZONTAL_REFLECT|"
        f"round={round_index}|pairing={p}|curve_degree={C*NS*F}|"
        f"DF_before={D*NS*F}|DF_after={Dnext*NS*F}|"
        f"class={','.join(map(str, C))}",
        flush=True,
    )
    D = Dnext
else:
    print("Q8HORIZONTAL_WARNING|max_rounds_reached=1", flush=True)

final_roots = enumerate_nonpositive_roots(D, args.max_degree)
final_negative = [rec for rec in final_roots if rec["D_pairing"] < 0]

status = (
    "NO_NEGATIVE_HORIZONTAL_ROOTS_THROUGH_BOUND"
    if not final_negative
    else "NEGATIVE_HORIZONTAL_ROOTS_REMAIN"
)

OUTPUT.write_text(json.dumps({
    "schema": "elkies-k3.h92-q8-horizontal-fixed-reduction.v1",
    "status": status,
    "max_old_fibre_degree_searched": int(args.max_degree),
    "initial_q8_class": list(map(int, D0)),
    "initial_old_fibre_degree": int(D0 * NS * F),
    "initial_negative_hits": initial_hits,
    "initial_unique_fiber_nef_roots": list(unique_initial.values()),
    "history": history,
    "final_class": list(map(int, D)),
    "final_square": int(D * NS * D),
    "final_old_fibre_degree": int(D * NS * F),
    "final_negative_roots": final_negative,
    "interpretation": (
        "The original q8 class was only reduced against old fibre roots. "
        "Negative positive-degree (-2)-curves are horizontal fixed-component/"
        "Weyl-chamber obstructions. The final class is the bounded horizontal "
        "chamber reduction obtained by exact root enumeration through the "
        "declared old-fibre degree."
    ),
}, indent=2, sort_keys=True) + "\n")

print(
    "Q8HORIZONTAL_SUMMARY|"
    f"initial_DF={D0*NS*F}|final_DF={D*NS*F}|"
    f"rounds={sum(rec['action']=='reflect' for rec in history)}|"
    f"remaining_negative={len(final_negative)}|status={status}",
    flush=True,
)
print(f"OUTPUT|{OUTPUT}", flush=True)
