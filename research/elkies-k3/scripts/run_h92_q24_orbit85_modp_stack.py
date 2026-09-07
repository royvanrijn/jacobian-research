#!/usr/bin/env python3
"""
Run the deterministic H3 q24/orbit85 -> D12 modular construction stack.

This is intentionally a thin runner, not a second implementation of the
geometry.  It follows the fixed construction path:

    q24 direct bridge + explicit I9* resolution
      -> D12 RR preflight
      -> affine D13 component graph
      -> effective D13 transport
      -> resolved component-valuation RR compiler
      -> 2x56 plane / quartic / D12 signature

There is no Hensel lift, no SRR retry, and no alternative neighbor search
here.  The purpose is to collect compatible modular signatures for the small
RR kernel so the characteristic-zero q24/orbit85 construction can be recovered
by CRT/rational reconstruction if direct QQ linear algebra is awkward.
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path


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
            and (candidate / "artifacts/local/elkies-k3").is_dir()
        ):
            return candidate
    raise SystemExit("Could not locate jacobian-research")


def parse_prime_tokens(tokens):
    primes = []
    for token in tokens:
        for part in str(token).replace(",", " ").split():
            if part:
                primes.append(int(part))
    if not primes:
        raise SystemExit("At least one prime is required")
    return primes


def load_status(path):
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text()).get("status")
    except Exception as exc:
        raise SystemExit(f"Could not read status from {path}: {exc}")


def require_prerequisites(local, prime):
    required = [
        local / f"q24-degree46-direct-global-mod-{prime}.json",
        local / f"q24-i9star-resolution-mod-{prime}.json",
    ]
    missing = [path for path in required if not path.exists()]
    if missing:
        lines = "\n".join(f"  - {path}" for path in missing)
        raise SystemExit(
            f"Prime {prime} is missing q24/orbit85 modular prerequisites:\n{lines}"
        )

    resolution_path = local / f"q24-i9star-resolution-mod-{prime}.json"
    resolution = json.loads(resolution_path.read_text())
    split_centers = (
        resolution.get("split_tangent_cone_centers")
        or resolution.get("split_center_labels")
        or []
    )
    canonical_split_centers = [
        "C" + str(label)[1:]
        if str(label).startswith("E") and str(label)[1:].isdigit()
        else str(label)
        for label in split_centers
    ]
    if canonical_split_centers != ["C10"]:
        raise SystemExit(
            f"Prime {prime} has split center {canonical_split_centers}; "
            "the current component-valuation compiler is canonicalized for C10. "
            "Use a prime with split center E10/C10, or first generalize the "
            "resolved-cluster relabeling."
        )


def stream_command(cmd, cwd):
    print("COMMAND|" + " ".join(str(part) for part in cmd), flush=True)
    proc = subprocess.Popen(
        cmd,
        cwd=str(cwd),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )
    assert proc.stdout is not None
    for line in proc.stdout:
        print(line, end="", flush=True)
    return proc.wait()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("primes", nargs="+", help="Prime(s), comma or space separated")
    parser.add_argument("--repo", type=Path)
    parser.add_argument("--sage", default="sage", help="Sage executable")
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-run steps even when their expected artifact already passes",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the planned steps without running Sage",
    )
    args = parser.parse_args()

    root = locate_repo(args.repo)
    scripts = root / "elkies-k3/scripts"
    local = root / "artifacts/local/elkies-k3"
    primes = parse_prime_tokens(args.primes)

    steps = [
        {
            "name": "rr-preflight",
            "script": scripts / "probe_h92_q24_d12_degree2_rr_modp.sage",
            "artifact": "q24-d12-rr-preflight-mod-{prime}.json",
            "statuses": {"PASS_H3_Q24_D12_MODP_RR_PREFLIGHT"},
        },
        {
            "name": "component-graph",
            "script": scripts / "certify_h92_q24_i9star_component_graph_modp.sage",
            "artifact": "q24-i9star-component-graph-mod-{prime}.json",
            "statuses": {"PASS_H3_Q24_AFFINE_D13_COMPONENT_GRAPH"},
        },
        {
            "name": "effective-d13-transport",
            "script": scripts / "certify_h92_q24_effective_d13_transport.sage",
            "artifact": "q24-effective-d13-transport-mod-{prime}.json",
            "statuses": {"PASS_EXACT_H3_Q24_EFFECTIVE_D13_TRANSPORT"},
        },
        {
            "name": "component-valuation-rr",
            "script": scripts / "probe_h92_q24_d12_component_valuation_rr_modp.sage",
            "artifact": "q24-d12-component-valuation-rr-mod-{prime}.json",
            "statuses": {
                "CANDIDATE_H3_Q24_EFFECTIVE_D13_D12_MODP",
                "PASS_H3_Q24_EFFECTIVE_D13_D12_MODP",
            },
        },
        {
            "name": "orbit85-d12-signature",
            "script": scripts / "extract_h92_q24_d12_modp_signature.sage",
            "artifact": "q24-orbit85-d12-signature-mod-{prime}.json",
            "statuses": {
                "PASS_H3_Q24_ORBIT85_D12_MODP_SIGNATURE",
                "CANDIDATE_H3_Q24_ORBIT85_D12_MODP_SIGNATURE",
            },
        },
    ]

    print(
        "Q24O85STACK_START|"
        f"repo={root}|primes={','.join(map(str, primes))}|"
        f"force={int(args.force)}|dry_run={int(args.dry_run)}",
        flush=True,
    )

    for prime in primes:
        require_prerequisites(local, prime)
        print(f"Q24O85STACK_PRIME|prime={prime}|status=START", flush=True)

        for step in steps:
            artifact = local / step["artifact"].format(prime=prime)
            status = load_status(artifact)
            if status in step["statuses"] and not args.force:
                print(
                    "Q24O85STACK_STEP|"
                    f"prime={prime}|step={step['name']}|artifact={artifact}|"
                    f"status={status}|action=SKIP_EXISTING",
                    flush=True,
                )
                continue

            cmd = [
                args.sage,
                "-python",
                str(step["script"]),
                "--repo",
                str(root),
                "--prime",
                str(prime),
            ]
            print(
                "Q24O85STACK_STEP|"
                f"prime={prime}|step={step['name']}|artifact={artifact}|action=RUN",
                flush=True,
            )
            if args.dry_run:
                print("COMMAND|" + " ".join(cmd), flush=True)
                continue

            ret = stream_command(cmd, root)
            if ret != 0:
                raise SystemExit(
                    f"Step {step['name']} failed for prime {prime} with exit {ret}"
                )

            status = load_status(artifact)
            if status not in step["statuses"]:
                allowed = ",".join(sorted(step["statuses"]))
                raise SystemExit(
                    f"Unexpected status for {artifact}: {status}; expected {allowed}"
                )
            print(
                "Q24O85STACK_STEP|"
                f"prime={prime}|step={step['name']}|artifact={artifact}|"
                f"status={status}|action=PASS",
                flush=True,
            )

        print(f"Q24O85STACK_PRIME|prime={prime}|status=PASS", flush=True)

    print(
        "Q24O85STACK_RESULT|"
        f"primes={','.join(map(str, primes))}|"
        "status=PASS_Q24_ORBIT85_MODP_STACK",
        flush=True,
    )


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.exit("Interrupted")
