#!/usr/bin/env python3
from pathlib import Path
import json
import subprocess
import sys

repo = Path.home() / "Documents" / "jacobian-research"
scripts = repo / "elkies-k3" / "scripts"
local = repo / "artifacts" / "local" / "elkies-k3"

# Continue beyond the currently explored range. Existing artifacts are reused.
primes = [
    100363,100379,100391,100393,100403,100411,100417,100447,
    100459,100469,100483,100493,100501,100511,100517,100519,
    100523,100537,100547,100549,100559,100591,100609,100613,
    100621,100649,100669,100673,100693,100699,
]

def jstatus(path):
    try:
        return json.loads(path.read_text()).get("status")
    except Exception:
        return None

def run(cmd, label):
    print(f"FULLNORM_BATCH|{label}|status=START", flush=True)
    subprocess.run([str(x) for x in cmd], cwd=repo, check=True)
    print(f"FULLNORM_BATCH|{label}|status=PASS", flush=True)

success_artifact = local / "q32-d12-fully-normalized-qq.json"

def success():
    return jstatus(success_artifact) == "PASS_EXACT_Q32_D12_FULLY_NORMALIZED_QQ_HELDOUT"

# Reuse the globally certified sign from prior successful q24 artifacts.
sign = None
for path in sorted(local.glob("q24-degree46-direct-global-mod-*.json")):
    try:
        d = json.loads(path.read_text())
        if d.get("status") != "PASS_MODULAR_Q24_FROM_DIRECT_DEGREE46_BRIDGE":
            continue
        s = int(d["direct_global_map"]["global_scale_root_sign"])
        if s not in (-1,1):
            continue
        if sign is None:
            sign = s
        elif sign != s:
            raise SystemExit(f"inconsistent global signs: saw {sign} and {s}")
    except Exception:
        continue
if sign is None:
    raise SystemExit("no certified q24 global sign found")

q24script = scripts / "recover_h92_q24_degree46_direct_global_modp_reuse_sign.sage"
if not q24script.exists():
    raise SystemExit(f"missing {q24script}; run the fresh-prime orientation patch first")

required = [
    scripts / "audit_h92_q24_i9star_exceptional_components_modp.sage",
    scripts / "probe_h92_q32_generic_divisorial_rr_modp.sage",
    scripts / "extract_h92_q32_modp_signature.sage",
    scripts / "probe_h92_q32_spinor_quartic_section_anchor.sage",
    scripts / "probe_h92_q32_pointed_spinor_weierstrass.sage",
    scripts / "reconstruct_h92_q32_fully_normalized_qq.sage",
]
for p in required:
    if not p.exists():
        raise SystemExit(f"missing required script {p}")

if success():
    print("FULLNORM_BATCH_RESULT|status=ALREADY_PASS", flush=True)
    raise SystemExit(0)

added = 0
bad = 0

for prime in primes:
    print(f"\nFULLNORM_PRIME|prime={prime}|status=START", flush=True)

    q24 = local / f"q24-degree46-direct-global-mod-{prime}.json"
    res = local / f"q24-i9star-resolution-mod-{prime}.json"
    div = local / f"q32-d12-generic-divval-mod-{prime}.json"
    sig = local / f"q32-signature-mod-{prime}.json"

    try:
        if jstatus(q24) != "PASS_MODULAR_Q24_FROM_DIRECT_DEGREE46_BRIDGE":
            run([
                "sage","-python",q24script,
                "--prime",prime,
                "--global-sign",sign,
                "--output",q24,
            ], f"prime={prime}|stage=q24")
        else:
            print(f"FULLNORM_BATCH|prime={prime}|stage=q24|status=SKIP_EXISTING", flush=True)

        if jstatus(res) != "PASS_EXPLICIT_MODP_I9STAR_D13_EXCEPTIONAL_COMPONENTS":
            run([
                "sage","-python",
                scripts/"audit_h92_q24_i9star_exceptional_components_modp.sage",
                "--prime",prime,
                "--output",res,
            ], f"prime={prime}|stage=i9star")
        else:
            print(f"FULLNORM_BATCH|prime={prime}|stage=i9star|status=SKIP_EXISTING", flush=True)

        # Bad reduction is expected occasionally. E10 must split into two spinors.
        try:
            rd = json.loads(res.read_text())
            e10 = next(c for c in rd["centers"] if c["label"] == "E10")
            nf = len(e10.get("tangent_factors", []))
        except Exception as exc:
            print(
                f"FULLNORM_PRIME|prime={prime}|reason={type(exc).__name__}:{exc}|"
                "status=SKIP_BAD_RESOLUTION",
                flush=True,
            )
            bad += 1
            continue

        if nf != 2:
            print(
                f"FULLNORM_PRIME|prime={prime}|E10_factors={nf}|"
                "status=SKIP_BAD_REDUCTION",
                flush=True,
            )
            bad += 1
            continue

        if jstatus(div) != "PASS_Q32_GENERIC_DIVISORIAL_RR_H0_TWO":
            run([
                "sage","-python",
                scripts/"probe_h92_q32_generic_divisorial_rr_modp.sage",
                "--prime",prime,
                "--output",div,
            ], f"prime={prime}|stage=divval")
        else:
            print(f"FULLNORM_BATCH|prime={prime}|stage=divval|status=SKIP_EXISTING", flush=True)

        if jstatus(div) != "PASS_Q32_GENERIC_DIVISORIAL_RR_H0_TWO":
            print(f"FULLNORM_PRIME|prime={prime}|status=SKIP_BAD_DIVVAL", flush=True)
            bad += 1
            continue

        if jstatus(sig) != "PASS_Q32_MODP_SIGNATURE":
            run([
                "sage","-python",
                scripts/"extract_h92_q32_modp_signature.sage",
                "--prime",prime,
            ], f"prime={prime}|stage=signature")
        else:
            print(f"FULLNORM_BATCH|prime={prime}|stage=signature|status=SKIP_EXISTING", flush=True)

        if jstatus(sig) != "PASS_Q32_MODP_SIGNATURE":
            bad += 1
            continue

        added += 1

        # Rebuild global marking artifacts from every currently available good prime.
        run([
            "sage","-python",
            scripts/"probe_h92_q32_spinor_quartic_section_anchor.sage",
        ], f"prime={prime}|stage=spinor")

        run([
            "sage","-python",
            scripts/"probe_h92_q32_pointed_spinor_weierstrass.sage",
        ], f"prime={prime}|stage=pointed")

        # This is now the only reconstruction criterion we care about.
        run([
            "sage","-python",
            scripts/"reconstruct_h92_q32_fully_normalized_qq.sage",
        ], f"prime={prime}|stage=fullnorm-qq")

        if success():
            d = json.loads(success_artifact.read_text())
            print(
                f"FULLNORM_BATCH_RESULT|added={added}|bad={bad}|last_prime={prime}|"
                f"modulus_bits={d.get('crt_modulus_bits')}|"
                "status=PASS_EXACT_Q32_D12_FULLY_NORMALIZED_QQ_HELDOUT",
                flush=True,
            )
            raise SystemExit(0)

        print(
            f"FULLNORM_PRIME|prime={prime}|added={added}|bad={bad}|"
            "status=MORE_PRECISION_NEEDED",
            flush=True,
        )

    except subprocess.CalledProcessError as exc:
        bad += 1
        print(
            f"FULLNORM_PRIME|prime={prime}|returncode={exc.returncode}|"
            "status=SKIP_FAILED_PRIME",
            flush=True,
        )
        continue
    except Exception as exc:
        bad += 1
        print(
            f"FULLNORM_PRIME|prime={prime}|"
            f"reason={type(exc).__name__}:{exc}|"
            "status=SKIP_FAILED_PRIME",
            flush=True,
        )
        continue

print(
    f"FULLNORM_BATCH_RESULT|added={added}|bad={bad}|"
    "status=EXHAUSTED_GOOD_PRIMES",
    flush=True,
)
