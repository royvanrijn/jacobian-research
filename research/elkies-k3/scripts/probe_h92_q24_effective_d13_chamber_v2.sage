#!/usr/bin/env sage -python
"""
Pure lattice/chamber diagnostic for H3 D13 --q24--> D12.

This script intentionally does NOT depend on the current ordinary-blowup
resolution artifact.  It uses the exact q24 translation script only to obtain
the deterministic D13 frame (`adapted`), the q24 child-frame coordinates
(`cD`), and the exact vertical root correction (`vr`).

It then:
  * Weyl-reduces the FULL q24 positive-frame vector to a dominant D13 chamber;
  * constructs the Weyl matrix W explicitly, avoiding convention ambiguity;
  * transports the deterministic D13 simple roots back to the effective q24
    chamber by W^{-1};
  * expresses the exact vertical correction in that effective component basis;
  * forms the forced vanishing cycle;
  * unloads it to the minimal anti-nef cycle;
  * reports the Du Val colength Z^2/2.

For a D13 rational double point the crepant-resolution colength of a complete
ideal represented by an anti-nef exceptional cycle Z is Z*C*Z/2, where C is
the positive D13 Cartan matrix.  The missing RR local codimension is expected
to be 8.
"""

import argparse
import contextlib
import io
import json
import sys
from pathlib import Path

from sage.all import QQ, ZZ, identity_matrix, matrix, vector


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
            (candidate / "elkies-k3/scripts").is_dir()
            and (candidate / "artifacts/generated-results").is_dir()
        ):
            return candidate
    raise SystemExit("Could not locate jacobian-research; pass --repo PATH")


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo", type=Path)
parser.add_argument("--output", type=Path)
args = parser.parse_args()

ROOT = locate_repo(args.repo)
LOCAL = ROOT / "artifacts/local/elkies-k3"
CLOSE = ROOT / "elkies-k3/scripts/close_h92_q8_q24_by_q6_translation.sage"
if not CLOSE.exists():
    raise SystemExit(f"Missing {CLOSE}")

# Run the exact q24 translation in a private scope.
saved_argv = list(sys.argv)
scope = {"__name__": "__embedded_q24_chamber_v2__"}
captured = io.StringIO()
try:
    sys.argv = [str(CLOSE)]
    with contextlib.redirect_stdout(captured):
        exec(compile(CLOSE.read_text(), str(CLOSE), "exec"), scope)
finally:
    sys.argv = saved_argv

required = ("adapted", "cD", "vr", "vf")
missing = [name for name in required if name not in scope]
if missing:
    raise SystemExit("q24 translation scope lacks: " + ",".join(missing))

adapted = matrix(ZZ, scope["adapted"])
cD = vector(ZZ, scope["cD"])
vr = vector(ZZ, scope["vr"])
vf = ZZ(scope["vf"])

assert adapted.dimensions() == (17, 17)
assert len(cD) == 19
assert len(vr) == 13
C = adapted[:13, :13]
assert C.det() == 4
assert all(C[i, i] == 2 for i in range(13))

l0 = vector(ZZ, cD[2:])
assert len(l0) == 17
norm0 = ZZ(l0 * adapted * l0)


def reflect_coordinate_row(row, simple_index):
    row = vector(ZZ, list(row))
    pairing = ZZ(row * adapted.column(simple_index))
    out = vector(ZZ, row)
    out[simple_index] -= pairing
    assert ZZ(out * adapted * out) == ZZ(row * adapted * row)
    return out, pairing


# Weyl reduce the full q24 positive-frame vector.
l = vector(ZZ, l0)
steps = []
for unused in range(10000):
    pairings = vector(ZZ, l * adapted[:, :13])
    negative = [i for i, value in enumerate(pairings) if value < 0]
    if not negative:
        break
    i = negative[0]
    l, pairing = reflect_coordinate_row(l, i)
    assert pairing < 0
    steps.append(i)
else:
    raise ArithmeticError("Weyl reduction did not terminate")

dominant_pairings = vector(ZZ, l * adapted[:, :13])
assert all(value >= 0 for value in dominant_pairings)
assert ZZ(l * adapted * l) == norm0

# Construct the exact row-action matrix W by applying the SAME sequence to the
# standard basis.  This removes all forward/reverse convention ambiguity:
#
#     l0 * W = l_dominant.
#
W = identity_matrix(ZZ, 17)
for i in steps:
    W = matrix(
        ZZ,
        [list(reflect_coordinate_row(row, i)[0]) for row in W.rows()],
    )
assert abs(ZZ(W.det())) == 1
assert W * adapted * W.transpose() == adapted
assert vector(ZZ, l0 * W) == l

Winv_q = W.inverse()
assert all(value in ZZ for value in Winv_q.list())
Winv = matrix(ZZ, Winv_q)

# beta_j = e_j * W^{-1}; then
#   <l0,beta_j> = <l0 W,e_j> = dominant_pairing_j.
B_eff = Winv[:13, :]
assert B_eff.dimensions() == (13, 17)
assert B_eff * adapted * B_eff.transpose() == C
replayed = vector(ZZ, l0 * adapted * B_eff.transpose())
assert replayed == dominant_pairings

print(
    "Q24CHAMBER_V2|"
    f"weyl_steps={len(steps)}|"
    f"dominant_pairings={','.join(map(str, dominant_pairings))}|"
    f"positive_hits={sum(value > 0 for value in dominant_pairings)}|"
    f"sum_hits={sum(dominant_pairings)}|"
    f"Wdet={W.det()}|replay=1|status=PASS",
    flush=True,
)

# Exact root correction in the deterministic basis is [vr,0,0,0,0].
root_vector = vector(ZZ, list(vr) + [0] * 4)

# Express it in the effective root basis beta_j.
coeff_q = matrix(QQ, B_eff.transpose()).solve_right(vector(QQ, root_vector))
assert all(value in ZZ for value in coeff_q)
coeff_beta = vector(ZZ, coeff_q)
assert coeff_beta * B_eff == root_vector

# Effective geometric components are C_j = -beta_j because
# D.C_j = <l0,beta_j> >= 0.
component_coeff = vector(ZZ, [-value for value in coeff_beta])

# div(f)+D >= 0 forces v_C(f) >= -D_C.
raw_cycle = vector(ZZ, [max(ZZ(0), value) for value in coeff_beta])

def positive_pairings(z):
    return vector(ZZ, z * C)

def unload_to_antinef(z):
    """Minimal componentwise enlargement with z*C >= 0."""
    z = vector(ZZ, z)
    record = []
    for unused in range(100000):
        p = positive_pairings(z)
        bad = [i for i, value in enumerate(p) if value < 0]
        if not bad:
            return z, record
        i = bad[0]
        # Add enough copies of E_i to make the i-th pairing nonnegative.
        # Adding E_i changes that pairing by +2.
        amount = ZZ((-p[i] + 1) // 2)
        if amount <= 0:
            amount = ZZ(1)
        z[i] += amount
        record.append((i, int(amount)))
    raise ArithmeticError("anti-nef unloading did not terminate")

anti_cycle, unload_record = unload_to_antinef(raw_cycle)
anti_pairings = positive_pairings(anti_cycle)
assert all(value >= 0 for value in anti_pairings)

raw_square = ZZ(raw_cycle * C * raw_cycle)
anti_square = ZZ(anti_cycle * C * anti_cycle)
raw_colength = QQ(raw_square) / 2
anti_colength = QQ(anti_square) / 2

print(
    "Q24VERTICAL_V2|"
    f"fibre_twist_represented_elsewhere={vf}|"
    f"deterministic_roots={','.join(map(str, vr))}|"
    f"effective_beta_coeff={','.join(map(str, coeff_beta))}|"
    f"component_coeff={','.join(map(str, component_coeff))}|"
    f"raw_cycle={','.join(map(str, raw_cycle))}|"
    f"raw_pairings={','.join(map(str, positive_pairings(raw_cycle)))}|"
    f"raw_colength={raw_colength}|status=PASS",
    flush=True,
)

print(
    "Q24ANTINEF_V2|"
    f"unload_steps={len(unload_record)}|"
    f"cycle={','.join(map(str, anti_cycle))}|"
    f"pairings={','.join(map(str, anti_pairings))}|"
    f"square={anti_square}|colength={anti_colength}|"
    f"expected=8|"
    f"status={'HIT_COLENGTH8' if anti_colength == 8 else 'DIAGNOSTIC'}",
    flush=True,
)

payload = {
    "schema": "elkies-k3.h3-q24-effective-d13-chamber-v2.v1",
    "status": (
        "HIT_EXPECTED_LOCAL_COLENGTH_8"
        if anti_colength == 8
        else "PASS_EFFECTIVE_D13_CHAMBER_V2"
    ),
    "q24": {
        "weyl_step_count": len(steps),
        "weyl_steps_1_based": [i + 1 for i in steps],
        "W_determinant": int(W.det()),
        "dominant_simple_pairings": list(map(int, dominant_pairings)),
        "vertical_fibre_coefficient_represented_elsewhere": int(vf),
        "deterministic_root_coefficients": list(map(int, vr)),
        "effective_beta_coefficients": list(map(int, coeff_beta)),
        "effective_component_coefficients": list(map(int, component_coeff)),
        "raw_required_cycle": list(map(int, raw_cycle)),
        "raw_cycle_pairings": list(map(int, positive_pairings(raw_cycle))),
        "raw_colength": str(raw_colength),
        "anti_nef_cycle": list(map(int, anti_cycle)),
        "anti_nef_pairings": list(map(int, anti_pairings)),
        "anti_nef_unloading": [[i + 1, amount] for i, amount in unload_record],
        "anti_nef_square_positive": int(anti_square),
        "anti_nef_colength": str(anti_colength),
    },
    "boundary": (
        "This is the abstract effective-D13 chamber and complete-ideal cycle "
        "diagnostic. It does not identify those effective roots with specific "
        "ordinary-blowup chart components."
    ),
}

OUT = (
    args.output.resolve()
    if args.output
    else LOCAL / "q24-effective-d13-chamber-v2.json"
)
OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(f"OUTPUT|{OUT}", flush=True)
print(
    "Q24CHAMBER_V2_RESULT|"
    f"weyl_steps={len(steps)}|"
    f"positive_hits={sum(v > 0 for v in dominant_pairings)}|"
    f"raw_colength={raw_colength}|anti_colength={anti_colength}|"
    f"status={payload['status']}",
    flush=True,
)
