#!/usr/bin/env python3
"""
status: HISTORICAL_DIAGNOSTIC
claim: old q24 R17/A11 progress wrapper ends in rejected zero-pole recovery.
superseded-by: fixed q24 D42 resolved-RR construction route.

Run the q24 D12 -> R17-directed A11 progress pipeline at one good prime.

Existing expensive prerequisites are reused when present.
"""
import subprocess
import sys
from pathlib import Path

ROOT = Path.home() / "Documents" / "jacobian-research"
S = ROOT / "elkies-k3" / "scripts"
L = ROOT / "artifacts" / "local" / "elkies-k3"
D = L / "q24-downstream-lift"
p = int(sys.argv[1]) if len(sys.argv) > 1 else 100003

def run(cmd):
    print("+", " ".join(map(str,cmd)), flush=True)
    subprocess.run([str(x) for x in cmd], cwd=ROOT, check=True)

def ensure(path, cmd):
    if path.exists():
        print(f"REUSE|{path.relative_to(ROOT)}", flush=True)
    else:
        run(cmd)

# Closeout lattice/marking certificate from equation D13 to pinned R17.
ensure(
    L / "q24-equation-d13-to-pinned-r17.json",
    ["sage","-python",S/"certify_h92_q24_equation_d13_to_pinned_r17.sage"],
)

# Modular q24 D12 + pointed spinor prerequisites.
ensure(
    L / f"q24-orbit85-d12-signature-mod-{p}.json",
    ["sage","-python",S/"extract_h92_q24_d12_modp_signature.sage","--prime",p],
)
ensure(
    L / f"q24-d12-spinor-anchor-mod-{p}.json",
    ["sage","-python",S/"anchor_h92_q24_d12_spinor_modp.sage","--prime",p],
)
ensure(
    D / f"d12-to-a11-equation-friendly-p{p}.json",
    ["sage","-python",S/"scan_h92_q24_d12_a11_equation_friendly.sage","--prime",p],
)
ensure(
    D / f"pointed-d12-a11-profile-p{p}.json",
    ["sage","-python",S/"profile_h92_q24_pointed_d12_a11.sage","--prime",p],
)
ensure(
    D / f"a11-easy-spinor-shift-p{p}.json",
    ["sage","-python",S/"analyze_h92_q24_a11_easy_spinor_shift.sage","--prime",p],
)

# Exact R17 marking gate.
run(["sage","-python",S/"identify_h92_q24_r17_a11_pointed_target.sage","--prime",p])

# Recover the marked section from the pointed zero-pole basis.
run([
    "sage","-python",S/"archive/recover_h92_q24_r17_a11_zero_pole_sections.sage",
    "--prime",p,
    "--target-json",D/f"r17-a11-pointed-target-p{p}.json",
])
