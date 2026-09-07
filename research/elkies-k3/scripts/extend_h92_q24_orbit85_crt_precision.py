#!/usr/bin/env python3
"""
Extend q24/orbit85 modular precision for QQ RR-kernel reconstruction.

For fresh primes this runner builds exactly the prerequisites needed by the
fixed q24/orbit85 construction stack:

  1. q24 degree-46 direct bridge, reusing the globally certified sign;
  2. explicit I9* exceptional resolution;
  3. skip bad reductions where the split center is not E10 with two factors;
  4. q24/orbit85 resolved-RR modular stack and compact signature;
  5. scalar CRT and compact-row LLL reconstruction attempts.

It does not run Hensel, SRR, or any alternate-neighbor search.
"""

import argparse
import json
import subprocess
from pathlib import Path


DEFAULT_PRIMES = [
    100393, 100403, 100411, 100417, 100447, 100459,
    100469, 100483, 100493, 100501, 100511, 100517,
    100519, 100523, 100537, 100547, 100549, 100559,
    100591, 100609, 100613, 100621, 100649, 100669,
    100673, 100693, 100699,
]

SINGULAR_FACTOR_LIMIT = 2**29

Q24_STATUS = "PASS_MODULAR_Q24_FROM_DIRECT_DEGREE46_BRIDGE"
RES_STATUS = "PASS_EXPLICIT_MODP_I9STAR_D13_EXCEPTIONAL_COMPONENTS"
SIG_STATUSES = {
    "PASS_H3_Q24_ORBIT85_D12_MODP_SIGNATURE",
    "CANDIDATE_H3_Q24_ORBIT85_D12_MODP_SIGNATURE",
}
CRT_PASS = "PASS_Q24_ORBIT85_RR_KERNEL_CRT_HELDOUT_VALIDATED"
LLL_PASS = "PASS_Q24_ORBIT85_COMPACT_ROWS_HELDOUT"


def locate_repo(explicit=None):
    candidates = []
    if explicit:
        candidates.append(Path(explicit).expanduser())
    cwd = Path.cwd().resolve()
    candidates += [cwd, *cwd.parents]
    home = Path.home()
    candidates += [
        home / "Documents/jacobian-research",
        home / "jacobian-research",
        home / "src/jacobian-research",
        home / "git/jacobian-research",
        home / "projects/jacobian-research",
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
            and (candidate / "artifacts/local/elkies-k3").is_dir()
        ):
            return candidate
    raise SystemExit("Could not locate jacobian-research")


def is_prime(n):
    n = int(n)
    if n < 2:
        return False
    small = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
    for p in small:
        if n == p:
            return True
        if n % p == 0:
            return False
    d = n - 1
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1
    # Deterministic for n < 2^64; massive overkill for our <2^29 primes.
    for a in [2, 3, 5, 7, 11, 13, 17]:
        if a >= n:
            continue
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False
    return True


def generated_primes(start, count):
    primes = []
    candidate = int(start)
    if candidate >= SINGULAR_FACTOR_LIMIT:
        candidate = SINGULAR_FACTOR_LIMIT - 1
    if candidate % 2 == 0:
        candidate -= 1
    while candidate >= 3 and len(primes) < int(count):
        if is_prime(candidate):
            primes.append(candidate)
        candidate -= 2
    if len(primes) < int(count):
        raise SystemExit(
            f"only generated {len(primes)} primes below {start}, needed {count}"
        )
    return primes


def parse_primes(raw, generated_start, generated_count):
    if generated_count:
        return generated_primes(generated_start, generated_count)
    if not raw:
        return list(DEFAULT_PRIMES)
    primes = []
    for token in str(raw).replace(",", " ").split():
        primes.append(int(token))
    return primes


def json_status(path):
    try:
        return json.loads(path.read_text()).get("status")
    except Exception:
        return None


def run(cmd, cwd, label):
    print(f"Q24O85EXTEND|{label}|status=START", flush=True)
    subprocess.run([str(part) for part in cmd], cwd=str(cwd), check=True)
    print(f"Q24O85EXTEND|{label}|status=PASS", flush=True)


def certified_global_sign(local):
    sign = None
    for path in sorted(local.glob("q24-degree46-direct-global-mod-*.json")):
        try:
            data = json.loads(path.read_text())
            if data.get("status") != Q24_STATUS:
                continue
            value = int(data["direct_global_map"]["global_scale_root_sign"])
            if value not in (-1, 1):
                continue
            if sign is None:
                sign = value
            elif sign != value:
                raise SystemExit(f"inconsistent q24 global signs: {sign} and {value}")
        except Exception:
            continue
    if sign is None:
        raise SystemExit("no certified q24 global sign found")
    return sign


def canonical_e10_split(resolution_path):
    try:
        data = json.loads(resolution_path.read_text())
        if data.get("status") != RES_STATUS:
            return False, f"resolution_status={data.get('status')}"
        e10 = next(item for item in data["centers"] if item["label"] == "E10")
        factor_count = len(e10.get("tangent_factors", []))
        return factor_count == 2, f"E10_factors={factor_count}"
    except Exception as exc:
        return False, f"{type(exc).__name__}:{exc}"


def success(local):
    crt = local / "q24-orbit85-rr-kernel-crt-qq.json"
    lll = local / "q24-orbit85-compact-row-lll.json"
    return json_status(crt) == CRT_PASS or json_status(lll) == LLL_PASS


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo", type=Path)
parser.add_argument("--primes", help="Comma/space separated candidate primes")
parser.add_argument(
    "--generate-large-primes",
    type=int,
    default=0,
    help="Generate this many descending primes below --large-prime-start",
)
parser.add_argument(
    "--large-prime-start",
    type=int,
    default=SINGULAR_FACTOR_LIMIT - 1,
    help="Upper bound for generated primes; must stay below 2^29 for Singular factorization",
)
parser.add_argument("--max-added", type=int, default=0, help="0 means no limit")
args = parser.parse_args()

repo = locate_repo(args.repo)
scripts = repo / "elkies-k3/scripts"
local = repo / "artifacts/local/elkies-k3"
primes = parse_primes(
    args.primes,
    args.large_prime_start,
    args.generate_large_primes,
)
sign = certified_global_sign(local)

required = [
    scripts / "recover_h92_q24_degree46_direct_global_modp_reuse_sign.sage",
    scripts / "audit_h92_q24_i9star_exceptional_components_modp.sage",
    scripts / "run_h92_q24_orbit85_modp_stack.py",
    scripts / "reconstruct_h92_q24_orbit85_rr_kernel_crt.sage",
    scripts / "reconstruct_h92_q24_orbit85_compact_rows_lll.sage",
]
for path in required:
    if not path.exists():
        raise SystemExit(f"missing required script: {path}")

if success(local):
    print("Q24O85EXTEND_RESULT|status=ALREADY_PASS", flush=True)
    raise SystemExit(0)

added = 0
bad = 0
failed = 0

for prime in primes:
    if args.max_added and added >= args.max_added:
        break

    print(f"\nQ24O85EXTEND_PRIME|prime={prime}|status=START", flush=True)
    q24 = local / f"q24-degree46-direct-global-mod-{prime}.json"
    res = local / f"q24-i9star-resolution-mod-{prime}.json"
    sig = local / f"q24-orbit85-d12-signature-mod-{prime}.json"

    try:
        if json_status(q24) != Q24_STATUS:
            run([
                "sage", "-python",
                scripts / "recover_h92_q24_degree46_direct_global_modp_reuse_sign.sage",
                "--prime", prime,
                "--global-sign", sign,
                "--output", q24,
            ], repo, f"prime={prime}|stage=q24-direct")
        else:
            print(
                f"Q24O85EXTEND|prime={prime}|stage=q24-direct|status=SKIP_EXISTING",
                flush=True,
            )

        if json_status(res) != RES_STATUS:
            run([
                "sage", "-python",
                scripts / "audit_h92_q24_i9star_exceptional_components_modp.sage",
                "--prime", prime,
                "--output", res,
            ], repo, f"prime={prime}|stage=i9star-resolution")
        else:
            print(
                f"Q24O85EXTEND|prime={prime}|stage=i9star-resolution|status=SKIP_EXISTING",
                flush=True,
            )

        good, reason = canonical_e10_split(res)
        if not good:
            bad += 1
            print(
                f"Q24O85EXTEND_PRIME|prime={prime}|reason={reason}|"
                "status=SKIP_BAD_REDUCTION",
                flush=True,
            )
            continue

        if json_status(sig) not in SIG_STATUSES:
            run([
                "python3",
                scripts / "run_h92_q24_orbit85_modp_stack.py",
                str(prime),
            ], repo, f"prime={prime}|stage=q24-orbit85-stack")
        else:
            print(
                f"Q24O85EXTEND|prime={prime}|stage=q24-orbit85-stack|status=SKIP_EXISTING",
                flush=True,
            )

        if json_status(sig) not in SIG_STATUSES:
            raise RuntimeError(f"signature did not pass for {prime}")

        added += 1

        run([
            "sage", "-python",
            scripts / "reconstruct_h92_q24_orbit85_rr_kernel_crt.sage",
        ], repo, f"prime={prime}|stage=scalar-crt")
        run([
            "sage", "-python",
            scripts / "reconstruct_h92_q24_orbit85_compact_rows_lll.sage",
        ], repo, f"prime={prime}|stage=compact-row-lll")

        if success(local):
            print(
                f"Q24O85EXTEND_RESULT|added={added}|bad={bad}|failed={failed}|"
                "status=PASS_RECONSTRUCTION",
                flush=True,
            )
            raise SystemExit(0)

    except subprocess.CalledProcessError as exc:
        failed += 1
        print(
            f"Q24O85EXTEND_PRIME|prime={prime}|exit={exc.returncode}|"
            "status=FAILED_SUBPROCESS",
            flush=True,
        )
    except Exception as exc:
        failed += 1
        print(
            f"Q24O85EXTEND_PRIME|prime={prime}|reason={type(exc).__name__}:{exc}|"
            "status=FAILED",
            flush=True,
        )

print(
    f"Q24O85EXTEND_RESULT|added={added}|bad={bad}|failed={failed}|"
    "status=NEED_MORE_PRECISION_OR_STRATEGY",
    flush=True,
)
