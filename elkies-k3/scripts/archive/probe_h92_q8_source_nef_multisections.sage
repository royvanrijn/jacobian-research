#!/usr/bin/env sage -python
"""
Audit whether the H92 q8 source class is truly nef against old-fibration
(-2)-curves, not only the 15 displayed fibre roots.

For a (-2)-class C=(c0,n,v) in U + (-Frame), with n=C.F>0,

    c0 = (Q(v)-2)/(2n).

For the q8 isotropic class D=(d0,18,d), D^2=0 gives

    D.C = Q(18*v - n*d)/(36*n) - 18/n.

Hence D.C < 0 iff

    Q(18*v - n*d) < 648.

For each old-fibre degree n this is an exact closest-vector/coset enumeration.
We convert it to a positive-definite rank-18 qfminim calculation by adding a
homogenizing coordinate k and enumerating

    Q(18*v - k*n*d) + k^2

through 648, then retain k=+/-1 and the integrality condition for c0.

Any retained C has C^2=-2 and C.F=n>0. By RR, -C cannot be effective because
F is nef and (-C).F<0, so C is an effective (-2)-curve. If D.C<0 it is a
genuine fixed component obstruction to nefness.

Run:
  sage -python ~/Downloads/probe_h92_q8_source_nef_multisections.sage

Optional:
  --repo /path/to/jacobian-research
  --max-degree 24
"""

import argparse
import json
from pathlib import Path

from sage.all import (
    ZZ, block_diagonal_matrix, matrix, pari, vector
)


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
args = parser.parse_args()

ROOT = locate_repo(args.repo)
FRAME = ROOT / "elkies-k3" / "data" / "fibrations" / "kumar_e7e8_mw2_frame_3.txt"
AMBIENT = ROOT / "artifacts" / "generated-results" / "elkies-k3-h92-q8-generic-rr-ambient.json"
OUTPUT = ROOT / "artifacts" / "generated-results" / "zz-h92-q8-source-nef-multisection-audit.json"

payload = json.loads(AMBIENT.read_text())
assert payload["status"] == "PASS_EXACT_Q8_GENERIC_RR_AMBIENT"

G = load_gram(FRAME)
assert G.nrows() == 17 and G.ncols() == 17 and G.is_positive_definite()

NS = block_diagonal_matrix(matrix(ZZ, ((0, 1), (1, 0))), -G)
F = vector(ZZ, [1, 0] + [0] * 17)
O = vector(ZZ, [-1, 1] + [0] * 17)
D = vector(ZZ, payload["source_q8_lattice_class"])

assert len(D) == 19
assert D * NS * D == 0
assert D * NS * F == 18
assert D[1] == 18

d0 = ZZ(D[0])
d = vector(ZZ, D[2:])
Qd = ZZ(d * G * d)
assert 36 * d0 == Qd

# Known section used by the generic-frame construction.
twice_minuscule = (2, 3, 4, 6, 5, 4, 3)
minus_P1 = vector(
    ZZ,
    [5, 1]
    + [-value for value in twice_minuscule]
    + [0] * 8
    + [1, 0],
)
assert minus_P1 * NS * minus_P1 == -2
assert minus_P1 * NS * F == 1

print(
    "Q8NEF_KNOWN|"
    f"D2={D*NS*D}|DF={D*NS*F}|"
    f"DO={D*NS*O}|DminusP1={D*NS*minus_P1}",
    flush=True,
)

negative = {}
zero = {}
counts = {}

# Exact bound:
# for k=1, augmented norm = Q(18v-n*d)+1.
# Negative D.C requires Q(...) < 648, hence augmented norm <=648.
ENUM_BOUND = ZZ(648)

for n in range(1, int(args.max_degree) + 1):
    nd = n * d
    Gd = G * nd.column()
    top = 324 * G
    cross = -18 * Gd
    bottom = ZZ(nd * G * nd) + 1

    A = matrix(ZZ, 18, 18)
    A[:17, :17] = top
    for i in range(17):
        A[i, 17] = cross[i, 0]
        A[17, i] = cross[i, 0]
    A[17, 17] = bottom
    assert A.is_positive_definite()

    # PARI returns one representative from each +/- pair as columns.
    raw = pari(A).qfminim(ENUM_BOUND)
    reps = matrix(ZZ, raw[2]).transpose() if int(raw[0]) else matrix(ZZ, 0, 18)

    candidates = {}
    for row in reps.rows():
        k = ZZ(row[17])
        if abs(k) != 1:
            continue
        rr = vector(ZZ, row)
        if k == -1:
            rr = -rr
            k = ZZ(1)
        assert rr[17] == 1
        v = vector(ZZ, rr[:17])

        qv = ZZ(v * G * v)
        numerator = qv - 2
        if numerator % (2 * n):
            continue
        c0 = numerator // (2 * n)
        C = vector(ZZ, [c0, n] + list(v))
        assert C * NS * C == -2
        assert C * NS * F == n

        w = 18 * v - n * d
        qw = ZZ(w * G * w)
        direct = ZZ(D * NS * C)
        formula_num = qw - 648
        assert formula_num % (36 * n) == 0
        formula = formula_num // (36 * n)
        assert direct == formula

        key = tuple(C)
        candidates[key] = {
            "degree": n,
            "D_pairing": int(direct),
            "Q_shift": int(qw),
            "class": list(map(int, C)),
        }

    neg_n = sorted(
        (entry for entry in candidates.values() if entry["D_pairing"] < 0),
        key=lambda e: (e["D_pairing"], e["Q_shift"], e["class"]),
    )
    zero_n = sorted(
        (entry for entry in candidates.values() if entry["D_pairing"] == 0),
        key=lambda e: (e["Q_shift"], e["class"]),
    )
    counts[n] = {
        "coset_root_candidates": len(candidates),
        "negative": len(neg_n),
        "zero": len(zero_n),
    }
    if neg_n:
        negative[n] = neg_n
    if zero_n:
        zero[n] = zero_n

    best = min(
        (entry["D_pairing"] for entry in candidates.values()),
        default=None,
    )
    print(
        "Q8NEF_DEGREE|"
        f"n={n}|candidates={len(candidates)}|"
        f"negative={len(neg_n)}|zero={len(zero_n)}|"
        f"best={best if best is not None else 'NA'}",
        flush=True,
    )

    # A negative root is already a decisive counterexample; print it fully,
    # but continue through max-degree so we see the local pattern.
    for index, entry in enumerate(neg_n[:10]):
        print(
            "Q8NEF_NEGATIVE|"
            f"n={n}|index={index}|Dpair={entry['D_pairing']}|"
            f"Qshift={entry['Q_shift']}|class={','.join(map(str, entry['class']))}",
            flush=True,
        )

status = (
    "FOUND_NEGATIVE_EFFECTIVE_ROOTS"
    if negative
    else "NO_NEGATIVE_ROOTS_THROUGH_BOUND"
)
OUTPUT.write_text(json.dumps({
    "schema": "elkies-k3.h92-q8-source-nef-multisection-audit.v1",
    "status": status,
    "max_old_fibre_degree": int(args.max_degree),
    "source_q8_class": list(map(int, D)),
    "known_pairings": {
        "D2": int(D * NS * D),
        "D_F": int(D * NS * F),
        "D_O": int(D * NS * O),
        "D_minus_P1": int(D * NS * minus_P1),
    },
    "degree_counts": {str(k): v for k, v in counts.items()},
    "negative_roots": {str(k): v for k, v in negative.items()},
    "zero_roots": {str(k): v for k, v in zero.items()},
    "logic": (
        "Every retained class C has C^2=-2 and C.F=n>0. "
        "By K3 Riemann-Roch one of +/-C is effective; -C cannot be effective "
        "because the old fibre F is nef and (-C).F<0. Hence C is effective. "
        "Any D.C<0 is therefore a genuine nefness/fixed-component obstruction."
    ),
}, indent=2, sort_keys=True) + "\n")

print(
    "Q8NEF_SUMMARY|"
    f"max_degree={args.max_degree}|"
    f"negative_degrees={','.join(map(str, sorted(negative))) if negative else '-'}|"
    f"status={status}",
    flush=True,
)
print(f"OUTPUT|{OUTPUT}", flush=True)
