#!/usr/bin/env sage -python
"""
Test the missing E7 identity-component regularity on the corrected 7D survivor.

The corrected q8 cusp-resolution calculations cover the seven nonidentity E7
components.  The old zero section O meets the remaining identity/affine
component at the smooth point at infinity, outside those affine blow-up charts.

In the common representative used by the global compiler, the global -10F is
placed at E8.  Hence the E7 identity component has divisor coefficient 0.

At its generic point:
    t=1/u is the base uniformizer,
    x and y are nonzero functions on the normalized special cubic,
    m=(y-yP)/(x-xP) is a unit,
    h(u) ~ h4*u^4.

Therefore
    u^d/h(u)^16 * x^a*m^b ~ t^(64-d) * unit,
and regularity along the omitted identity component is exactly d<=64 for all
16 generic families.

This script restricts that condition to the already-computed 7D survivor
after generic E7 + translated E7.

Run:
  sage -python ~/Downloads/probe_h92_q8_corrected1278_e7_identity.sage
"""

import argparse
import json
from pathlib import Path

from sage.all import GF, matrix


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
            and (candidate / "artifacts" / "generated-results").is_dir()
        ):
            return candidate
    raise SystemExit("Could not locate jacobian-research; pass --repo PATH")


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo", type=Path)
parser.add_argument("--prime", type=int, default=43)
args = parser.parse_args()

ROOT = locate_repo(args.repo)
GEN = ROOT / "artifacts" / "generated-results"
P = int(args.prime)
finite = GF(P)

AMBIENT = GEN / "zz-h92-q8-corrected1278-ambient.json"
GLOBAL = GEN / f"zz-h92-q8-corrected1278-global-kernel-mod-{P}.json"
TRANSLATED = GEN / f"zz-h92-q8-corrected1278-two-translated-divisors-mod-{P}.json"
OUTPUT = GEN / f"zz-h92-q8-corrected1278-e7-identity-mod-{P}.json"

for path in (AMBIENT, GLOBAL, TRANSLATED):
    if not path.exists():
        raise SystemExit(f"Missing prerequisite: {path}")

ambient = json.loads(AMBIENT.read_text())
global_kernel = json.loads(GLOBAL.read_text())
translated = json.loads(TRANSLATED.read_text())

assert int(ambient["ambient_dimension"]) == 1278
assert int(global_kernel["dimensions"]["kernel"]) == 14
assert int(translated["combined"]["remaining_dimension"]) == 7

basis = ambient["ambient_basis"]

K = matrix(
    finite,
    [[finite(v) for v in row] for row in global_kernel["kernel_basis_rows"]],
)
assert K.nrows() == 14 and K.ncols() == 1278

trows = []
for record in translated["divisors"]:
    trows.extend([[finite(v) for v in row] for row in record["row_space_basis"]])
T = matrix(finite, trows).row_space().basis_matrix()
assert T.rank() == 7 and T.ncols() == 14

C7 = T.right_kernel().basis_matrix()
assert C7.nrows() == 7
S7 = C7 * K
assert S7.rank() == 7

excluded = [
    i for i, entry in enumerate(basis)
    if int(entry["actual_u_power"]) > 64
]
assert len(excluded) == 238, len(excluded)

M = S7[:, excluded].transpose()
rank = int(M.rank())
remaining = 7 - rank

# Degree-by-degree fingerprint.
records = []
running = []
running_rank = 0
for d in sorted({
    int(basis[i]["actual_u_power"]) for i in excluded
}, reverse=True):
    cols = [
        i for i in excluded
        if int(basis[i]["actual_u_power"]) == d
    ]
    Md = S7[:, cols].transpose()
    rd = int(Md.rank())
    before = running_rank
    running.extend(cols)
    running_rank = int(S7[:, running].transpose().rank())
    gain = running_rank - before
    records.append({
        "degree": d,
        "columns": len(cols),
        "rank": rd,
        "gain": gain,
        "cumulative": running_rank,
    })
    print(
        f"E7IDENTITY_DEGREE|d={d}|cols={len(cols)}|rank={rd}|"
        f"gain={gain}|cumulative={running_rank}",
        flush=True,
    )

C = M.right_kernel().basis_matrix()
S = C * S7
assert S.nrows() == remaining and S.rank() == remaining

print(
    "Q8E7IDENTITY|"
    f"prime={P}|before=7|excluded_columns={len(excluded)}|"
    f"rank={rank}|remaining={remaining}|cap=64",
    flush=True,
)

OUTPUT.write_text(json.dumps({
    "schema": "elkies-k3.h92-q8-corrected1278-e7-identity-modp.v1",
    "status": "EXPERIMENTAL_MODULAR_CORRECTED_Q8_E7_IDENTITY_COMPONENT",
    "prime": P,
    "input_dimension": 7,
    "identity_component": {
        "divisor_coefficient_in_common_representative": 0,
        "generic_orders": {"t": 1, "x": 0, "m": 0},
        "coefficient_order": "ord_t(u^d/h^16)=64-d",
        "degree_cap": 64,
    },
    "excluded_columns": len(excluded),
    "rank": rank,
    "remaining": remaining,
    "degree_records": records,
    "survivor_basis_rows": [[int(v) for v in row] for row in S.rows()],
    "boundary": (
        "This is a modular restriction of the omitted E7 identity-component "
        "valuation to the corrected 7D survivor. If it cuts to dimension 2, "
        "the next step is to rebuild the global ambient with this component "
        "included from the start and repeat at another good prime."
    ),
}, indent=2, sort_keys=True) + "\n")

print(f"OUTPUT|{OUTPUT}", flush=True)
