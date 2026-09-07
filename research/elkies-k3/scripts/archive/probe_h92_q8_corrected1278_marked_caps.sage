#!/usr/bin/env sage -python
"""
Test the corrected marked-E7 degree caps directly on the current 7D survivor.

Current corrected1278 ambient used helper-shifted caps
    m^b:    d<=78, b=0..8
    x*m^b:  d<=80, b=0..6.

But on the resolved marked chart E7_5, the exact q6 certificate says the raw
frame <1,m> is valid.  For the movable q8 divisor the local E7 contribution is

    8*(-P1) + 6*E7_5

(with the global -10F represented at E8).  Since near -P1
    m ~ unit/W,
    x ~ t^2 * unit,
the local q8 generator is t^-6*m^8.  Therefore a true coefficient
u^d/h^16 ~ t^(64-d) must satisfy

    m-family:    64-d >= -6       => d<=70
    x*m-family:  64-d+2 >= -6     => d<=72.

This script reconstructs the 7D space after the corrected translated E7
condition and imposes only these sharper marked-point caps.  It reports the
rank of the 128 high-degree coefficients and the remaining dimension.

Run:
  sage -python ~/Downloads/probe_h92_q8_corrected1278_marked_caps.sage

Optional:
  --repo /path/to/jacobian-research
  --prime 43
"""

import argparse
import json
from pathlib import Path

from sage.all import GF, ZZ, matrix


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
OUTPUT = GEN / f"zz-h92-q8-corrected1278-marked-caps-mod-{P}.json"

for path in (AMBIENT, GLOBAL, TRANSLATED):
    if not path.exists():
        raise SystemExit(f"Missing prerequisite: {path}")

ambient = json.loads(AMBIENT.read_text())
global_kernel = json.loads(GLOBAL.read_text())
translated = json.loads(TRANSLATED.read_text())

assert int(ambient["ambient_dimension"]) == 1278
assert int(global_kernel["prime"]) == P
assert int(global_kernel["dimensions"]["kernel"]) == 14
assert int(translated["prime"]) == P
assert int(translated["common"]["global_survivor_dimension"]) == 14
assert int(translated["combined"]["restricted_rank"]) == 7
assert int(translated["combined"]["remaining_dimension"]) == 7

basis = ambient["ambient_basis"]
assert len(basis) == 1278

K = matrix(
    finite,
    [[finite(v) for v in row] for row in global_kernel["kernel_basis_rows"]],
)
assert K.nrows() == 14 and K.ncols() == 1278 and K.rank() == 14

translated_rows = []
for record in translated["divisors"]:
    translated_rows.extend(
        [[finite(v) for v in row] for row in record["row_space_basis"]]
    )
T = matrix(finite, translated_rows)
T = T.row_space().basis_matrix()
assert T.ncols() == 14 and T.rank() == 7

# Coordinates in the 14D global-kernel basis which survive translated E7.
C7 = T.right_kernel().basis_matrix()
assert C7.nrows() == 7 and C7.ncols() == 14

# Seven surviving actual ambient coefficient vectors.
S7 = C7 * K
assert S7.nrows() == 7 and S7.ncols() == 1278 and S7.rank() == 7

high_cols = []
high_by_family = []
for col, entry in enumerate(basis):
    a = int(entry["x_power"])
    b = int(entry["m_power"])
    d = int(entry["actual_u_power"])
    cap = 70 if a == 0 else 72
    if d > cap:
        high_cols.append(col)
        high_by_family.append((col, a, b, d))

assert len(high_cols) == 128, len(high_cols)

# Requiring all excluded coefficients to vanish gives a map S7 -> F_p^128.
Mhigh = S7[:, high_cols].transpose()
high_rank = int(Mhigh.rank())
remaining = 7 - high_rank

# Separate ranks for diagnostic purposes.
m_cols = [col for col,a,b,d in high_by_family if a == 0]
x_cols = [col for col,a,b,d in high_by_family if a == 1]
Mm = S7[:, m_cols].transpose()
Mx = S7[:, x_cols].transpose()
rank_m = int(Mm.rank())
rank_x = int(Mx.rank())
rank_union = high_rank

# Degree-by-degree incremental ranks.
degree_records = []
running_cols = []
running_rank = 0
for d in sorted({d for _,_,_,d in high_by_family}, reverse=True):
    cols_d = [col for col,a,b,dd in high_by_family if dd == d]
    Md = S7[:, cols_d].transpose()
    rank_d = int(Md.rank())
    before = running_rank
    running_cols.extend(cols_d)
    running_rank = int(S7[:, running_cols].transpose().rank())
    degree_records.append({
        "degree": d,
        "columns": len(cols_d),
        "rank_at_degree": rank_d,
        "incremental_rank_high_to_low": running_rank - before,
        "cumulative_rank": running_rank,
    })
    print(
        "MARKEDCAP_DEGREE|"
        f"d={d}|cols={len(cols_d)}|rank={rank_d}|"
        f"gain={running_rank-before}|cumulative={running_rank}",
        flush=True,
    )

# Surviving 2D (or whatever remains) actual ambient basis.
Cmarked = Mhigh.right_kernel().basis_matrix()
assert Cmarked.nrows() == remaining and Cmarked.ncols() == 7
Smarked = Cmarked * S7
assert Smarked.rank() == remaining

print(
    "Q8MARKEDCAPS|"
    f"prime={P}|before=7|excluded_columns={len(high_cols)}|"
    f"m_excluded={len(m_cols)}|x_excluded={len(x_cols)}|"
    f"rank_m={rank_m}|rank_x={rank_x}|combined_rank={rank_union}|"
    f"remaining={remaining}",
    flush=True,
)

OUTPUT.write_text(json.dumps({
    "schema": "elkies-k3.h92-q8-corrected1278-marked-caps-modp.v1",
    "status": "EXPERIMENTAL_MODULAR_CORRECTED_Q8_MARKED_E7_CAPS",
    "prime": P,
    "input_dimension": 7,
    "marked_local_geometry": {
        "component": "E7_5",
        "horizontal_multiplicity_at_minus_P1": 8,
        "vertical_E7_5_coefficient": 6,
        "m_behavior": "m~unit/W",
        "x_behavior": "x~t^2*unit",
        "coefficient_behavior": "u^d/h^16~t^(64-d)*unit",
        "local_generator": "t^-6*m^8",
        "m_degree_cap": 70,
        "xm_degree_cap": 72,
    },
    "excluded_columns": len(high_cols),
    "ranks": {
        "m_family": rank_m,
        "x_family": rank_x,
        "combined": high_rank,
        "remaining": remaining,
    },
    "degree_records": degree_records,
    "survivor_basis_rows": [
        [int(v) for v in row] for row in Smarked.rows()
    ],
    "boundary": (
        "This tests the direct resolved marked-point divisor condition on the "
        "existing corrected 7D modular survivor. A remaining dimension 2 is "
        "strong evidence for the q8 pencil but should be repeated at another "
        "good prime and incorporated into a rebuilt sharp ambient."
    ),
}, indent=2, sort_keys=True) + "\n")

print(f"OUTPUT|{OUTPUT}", flush=True)
