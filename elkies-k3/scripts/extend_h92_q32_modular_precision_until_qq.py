#!/usr/bin/env python3
from pathlib import Path
import json
import subprocess
import sys

repo = Path.home() / "Documents" / "jacobian-research"
scripts = repo / "elkies-k3" / "scripts"
local = repo / "artifacts" / "local" / "elkies-k3"

# Fresh good-prime candidates. We stop as soon as projective LLL passes.
primes = [
    100297,100313,100333,100343,100357,100361,100363,100379,
    100391,100393,100403,100411,100417,100447,100459,100469,
    100483,100493,100501,100511,100517,100519,100523,100537,
]

def json_status(path):
    try:
        return json.loads(path.read_text()).get("status")
    except Exception:
        return None

def run(cmd, label):
    print(f"BATCH|{label}|cmd={' '.join(map(str,cmd))}|status=START", flush=True)
    subprocess.run([str(x) for x in cmd], cwd=repo, check=True)
    print(f"BATCH|{label}|status=PASS", flush=True)

def success_now():
    out = local / "q32-pointed-d12-qq-projective-lll.json"
    return json_status(out) == "PASS_EXACT_Q32_POINTED_D12_QQ_PROJECTIVE_LLL_HELDOUT"

required_scripts = [
    "recover_h92_q24_degree46_direct_global_modp_reuse_sign.sage",
    "audit_h92_q24_i9star_exceptional_components_modp.sage",
    "probe_h92_q32_generic_divisorial_rr_modp.sage",
    "extract_h92_q32_modp_signature.sage",
    "probe_h92_q32_spinor_quartic_section_anchor.sage",
    "probe_h92_q32_pointed_spinor_weierstrass.sage",
    "reconstruct_h92_q32_pointed_d12_projective_lll.sage",
]
for name in required_scripts:
    if not (scripts / name).exists():
        raise SystemExit(f"missing required script: {scripts/name}")

if success_now():
    print("BATCH_RESULT|status=ALREADY_PASS_EXACT_Q32_POINTED_D12_QQ_PROJECTIVE_LLL_HELDOUT", flush=True)
    raise SystemExit(0)

completed = 0
failures = 0

for p in primes:
    print(f"\nBATCH_PRIME|prime={p}|status=START", flush=True)

    q24 = local / f"q24-degree46-direct-global-mod-{p}.json"
    res = local / f"q24-i9star-resolution-mod-{p}.json"
    div = local / f"q32-d12-generic-divval-mod-{p}.json"
    sig = local / f"q32-signature-mod-{p}.json"

    try:
        if json_status(q24) != "PASS_MODULAR_Q24_FROM_DIRECT_DEGREE46_BRIDGE":
            run([
                "sage", "-python",
                scripts / "recover_h92_q24_degree46_direct_global_modp_reuse_sign.sage",
                "--prime", p,
                "--global-sign", 1,
                "--output", q24,
            ], f"prime={p}|stage=q24-horizontal")
        else:
            print(f"BATCH|prime={p}|stage=q24-horizontal|status=SKIP_EXISTING", flush=True)

        if json_status(res) != "PASS_EXPLICIT_MODP_I9STAR_D13_EXCEPTIONAL_COMPONENTS":
            run([
                "sage", "-python",
                scripts / "audit_h92_q24_i9star_exceptional_components_modp.sage",
                "--prime", p,
                "--output", res,
            ], f"prime={p}|stage=i9star-resolution")
        else:
            print(f"BATCH|prime={p}|stage=i9star-resolution|status=SKIP_EXISTING", flush=True)

        # Skip bad-reduction primes where the D13 I9* resolution does not
        # have the expected split E10 center with two spinor tangent factors.
        try:
            rdat = json.loads(res.read_text())
            e10 = next(c for c in rdat["centers"] if c["label"] == "E10")
            if len(e10.get("tangent_factors", [])) != 2:
                print(
                    f"BATCH_PRIME|prime={p}|stage=q32-divval|"
                    f"E10_factors={len(e10.get('tangent_factors', []))}|"
                    "status=SKIP_BAD_REDUCTION",
                    flush=True,
                )
                continue
        except Exception as exc:
            print(
                f"BATCH_PRIME|prime={p}|stage=q32-divval|"
                f"reason={type(exc).__name__}:{exc}|status=SKIP_BAD_RESOLUTION_ARTIFACT",
                flush=True,
            )
            continue

        if json_status(div) != "PASS_Q32_GENERIC_DIVISORIAL_RR_H0_TWO":
            run([
                "sage", "-python",
                scripts / "probe_h92_q32_generic_divisorial_rr_modp.sage",
                "--prime", p,
                "--output", div,
            ], f"prime={p}|stage=q32-divval")
        else:
            print(f"BATCH|prime={p}|stage=q32-divval|status=SKIP_EXISTING", flush=True)

        if json_status(sig) != "PASS_Q32_MODP_SIGNATURE":
            run([
                "sage", "-python",
                scripts / "extract_h92_q32_modp_signature.sage",
                "--prime", p,
            ], f"prime={p}|stage=q32-signature")
        else:
            print(f"BATCH|prime={p}|stage=q32-signature|status=SKIP_EXISTING", flush=True)

        if json_status(sig) != "PASS_Q32_MODP_SIGNATURE":
            raise RuntimeError(f"signature did not pass for {p}")

        completed += 1

        # Rebuild global spinor anchors from all available signatures.
        run([
            "sage", "-python",
            scripts / "probe_h92_q32_spinor_quartic_section_anchor.sage",
        ], f"prime={p}|stage=spinor-anchor")

        run([
            "sage", "-python",
            scripts / "probe_h92_q32_pointed_spinor_weierstrass.sage",
        ], f"prime={p}|stage=pointed-anchor")

        # Try exact projective reconstruction after every added prime.
        run([
            "sage", "-python",
            scripts / "reconstruct_h92_q32_pointed_d12_projective_lll.sage",
        ], f"prime={p}|stage=projective-lll")

        if success_now():
            print(
                f"BATCH_RESULT|added_primes={completed}|last_prime={p}|"
                "status=PASS_EXACT_Q32_POINTED_D12_QQ_PROJECTIVE_LLL_HELDOUT",
                flush=True,
            )
            raise SystemExit(0)

        print(
            f"BATCH_PRIME|prime={p}|added={completed}|"
            "status=MORE_PRECISION_NEEDED",
            flush=True,
        )

    except subprocess.CalledProcessError as exc:
        failures += 1
        print(
            f"BATCH_PRIME|prime={p}|returncode={exc.returncode}|"
            "status=SKIP_FAILED_PRIME",
            flush=True,
        )
        if failures >= 4:
            raise SystemExit("too many failed primes; inspect the first failure before continuing")

print(
    f"BATCH_RESULT|added_primes={completed}|failed_primes={failures}|"
    "status=EXHAUSTED_CANDIDATES_NEED_STRONGER_PIVOT",
    flush=True,
)
