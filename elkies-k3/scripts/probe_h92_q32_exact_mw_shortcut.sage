#!/usr/bin/env sage -python
import json
from pathlib import Path
from sage.all import ZZ, matrix, vector

ROOT=Path(__file__).resolve().parents[2]
LOCAL=ROOT/"artifacts/local/elkies-k3"

OLD=LOCAL/"q8-explicit-old-curves.json"
if not OLD.exists():
    raise SystemExit(f"missing {OLD}")

d=json.loads(OLD.read_text())
assert d["status"]=="PASS_EXACT_Q8_EXPLICIT_OLD_CURVE_PROFILE"

G1=vector(ZZ,d["anchor"]["G1"])
idx=d["old_index3_lattice"]
A=vector(ZZ,idx["A"])
B=vector(ZZ,idx["B"])
C=vector(ZZ,idx["C"])
G3=vector(ZZ,(0,0,1,0))

# This is the exact relation certified by the artifact's own construction.
target=C-20*G1-46*G3

sections=[]
for rec in d["curves"]:
    if int(rec["q8_degree"])==1:
        sections.append((rec["curve"],vector(ZZ,rec["mw_coordinates"])))

# G1 and G3 are already explicit equation-level sections by separate exact
# certificates. Add them as available generators.
available=[("G1",G1),("G3",G3)]
for name,z in sections:
    if all(z != zz for _,zz in available):
        available.append((name,z))

print(
    "Q32MWSHORT_INPUT|"
    f"target={','.join(map(str,target))}|"
    f"explicit_degree1={len(sections)}|"
    f"available={';'.join(name+':'+','.join(map(str,z)) for name,z in available)}|"
    "status=PASS",
    flush=True,
)

M=matrix(ZZ,[list(z) for _,z in available]).transpose()
rank=M.rank()

# Solve M*c = target over QQ; check integrality.
sol=None
if rank==4:
    # right_kernel based exact solve when rectangular
    aug=M.augment(matrix(ZZ,4,1,list(target)))
    if aug.rank()==rank:
        qsol=M.solve_right(target)
        if all(x in ZZ for x in qsol):
            sol=vector(ZZ,qsol)

print(
    "Q32MWSHORT_SPAN|"
    f"rank={rank}|generators={len(available)}|"
    f"integral_target={int(sol is not None)}|"
    f"coeffs={'NA' if sol is None else ','.join(map(str,sol))}|"
    "status="+("PASS_DIRECT_MW_SHORTCUT" if sol is not None else "NO_DIRECT_SHORTCUT"),
    flush=True,
)

# Also inspect the artifact's explicitly named q8 section candidates.
for rec in d.get("q8_section_candidates",[]):
    print(
        "Q32MWSHORT_SECTION|"
        f"curve={rec['curve']}|mw={','.join(map(str,rec['mw_coordinates']))}|"
        f"height={rec['height']}|PdotO={rec['P_dot_O']}|"
        f"correction={rec['local_correction']}|status=PASS",
        flush=True,
    )
