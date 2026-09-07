#!/usr/bin/env python3
"""
Run/reuse all four modular embeddings for one split prime and combine the
normalized j(V)=N(V)/D(V) coefficients into the basis 1,s,j,s*j mod p.

This script deliberately does NOT reinstall q80_modular_j_probe.sage, so a
locally patched/working worker is preserved.

Usage (from any working directory):
  python3 elkies-k3/scripts/q80-orbit1222-char0/modular_packet.py --prime 79
"""

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
MODULAR_DATA = (
    REPO_ROOT
    / "artifacts"
    / "generated-results"
    / "q80-orbit1222-char0"
    / "modular"
)
EMBEDDINGS = ((1, 1), (1, -1), (-1, 1), (-1, -1))


def filename(repo, p, ss, sj):
    return repo / f"q80_modj_p{p}_s{ss:+d}_j{sj:+d}.json"


def valid_record(path, p, ss, sj):
    if not path.exists():
        return False
    try:
        r = json.loads(path.read_text())
    except Exception:
        return False
    return (
        r.get("version") == 1
        and int(r.get("prime", -1)) == p
        and int(r.get("sign_s", 0)) == ss
        and int(r.get("sign_j", 0)) == sj
        and len(r.get("kernel", ())) == 50
    )


def load(path):
    return json.loads(path.read_text())


def modinv(a, p):
    return pow(a % p, -1, p)


def combine(values, rs, rj, p):
    """
    values[(eps_s,eps_j)] = sigma(c) for
      c = a + b*s + c*j + d*s*j.
    """
    fpp = values[(1, 1)] % p
    fpm = values[(1, -1)] % p
    fmp = values[(-1, 1)] % p
    fmm = values[(-1, -1)] % p

    inv4 = modinv(4, p)
    a = (fpp + fpm + fmp + fmm) * inv4 % p
    b = (fpp + fpm - fmp - fmm) * modinv(4 * rs, p) % p
    c = (fpp - fpm + fmp - fmm) * modinv(4 * rj, p) % p
    d = (fpp - fpm - fmp + fmm) * modinv(4 * rs * rj, p) % p
    return (a, b, c, d)


def evaluate_component(component, ss, sj, rs, rj, p):
    a, b, c, d = component
    return (
        a
        + b * (ss * rs)
        + c * (sj * rj)
        + d * (ss * sj * rs * rj)
    ) % p


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--modular-data", type=Path, default=MODULAR_DATA)
    ap.add_argument("--prime", type=int, default=79)
    ap.add_argument("--jobs", type=int, default=4)
    args = ap.parse_args()

    modular_data = args.modular_data.expanduser().resolve()
    modular_data.mkdir(parents=True, exist_ok=True)
    worker = HERE / "modular_j_probe.sage"
    if not worker.exists():
        raise SystemExit(
            f"missing {worker}; run/install the modular j probe first"
        )

    # Make the known serialization fix persistent without otherwise modifying
    # the working modular engine.
    text = worker.read_text()
    old = "json.dumps(record, indent=2, sort_keys=True)"
    new = "json.dumps(record, indent=2, sort_keys=True, default=int)"
    if old in text:
        worker.write_text(text.replace(old, new, 1))
        print(f"patched JSON serializer: {worker}")

    sage = shutil.which("sage") or "/usr/local/bin/sage"

    missing = [
        (ss, sj) for ss, sj in EMBEDDINGS
        if not valid_record(filename(modular_data, args.prime, ss, sj), args.prime, ss, sj)
    ]
    cached = 4 - len(missing)
    print(
        f"Q80MODPACKET|p={args.prime}|cached_embeddings={cached}|"
        f"missing_embeddings={len(missing)}",
        flush=True,
    )

    live = []
    failures = []
    pending = list(missing)
    jobs = max(1, min(args.jobs, 4))

    while pending or live:
        while pending and len(live) < jobs:
            ss, sj = pending.pop(0)
            out = filename(modular_data, args.prime, ss, sj)
            cmd = [
                sage, str(worker),
                "--prime", str(args.prime),
                "--sign-s", str(ss),
                "--sign-j", str(sj),
                "--out", str(out),
            ]
            print("+", " ".join(cmd), flush=True)
            live.append(((ss, sj), subprocess.Popen(cmd, cwd=REPO_ROOT)))

        emb, proc = live.pop(0)
        rc = proc.wait()
        if rc:
            failures.append((emb, rc))

    if failures:
        for emb, rc in failures:
            print(f"FAILED embedding={emb} rc={rc}", file=sys.stderr)
        raise SystemExit(1)

    records = {
        emb: load(filename(modular_data, args.prime, *emb))
        for emb in EMBEDDINGS
    }

    p = args.prime
    rs = int(records[(1, 1)]["s_root"]) % p
    rj = int(records[(1, 1)]["j_root"]) % p
    if rs == 0 or rj == 0 or rs*rs % p != (-6) % p or rj*rj % p != (-3) % p:
        raise ArithmeticError(f"bad positive roots rs={rs}, rj={rj}")

    # Every sign file must be using the same underlying positive root.
    for ss, sj in EMBEDDINGS:
        r = records[(ss, sj)]
        if int(r["s_root"]) % p != ss*rs % p:
            raise ArithmeticError(f"inconsistent s root for {(ss,sj)}")
        if int(r["j_root"]) % p != sj*rj % p:
            raise ArithmeticError(f"inconsistent j root for {(ss,sj)}")

    components = []
    for i in range(50):
        vals = {
            emb: int(records[emb]["kernel"][i]) % p
            for emb in EMBEDDINGS
        }
        comp = combine(vals, rs, rj, p)

        # Exact round-trip check against all four observed embeddings.
        for ss, sj in EMBEDDINGS:
            observed = vals[(ss, sj)]
            recovered = evaluate_component(comp, ss, sj, rs, rj, p)
            if recovered != observed:
                raise ArithmeticError(
                    f"coefficient {i} roundtrip failed at {(ss,sj)}: "
                    f"{recovered} != {observed}"
                )
        components.append(comp)

    names = ("1", "s", "j", "sj")
    support_counts = {
        names[k]: sum(1 for c in components if c[k] % p != 0)
        for k in range(4)
    }
    support_indices = {
        names[k]: [i for i,c in enumerate(components) if c[k] % p != 0]
        for k in range(4)
    }

    packet = {
        "version": 1,
        "prime": p,
        "s_positive_root": rs,
        "j_positive_root": rj,
        "normalization": "kernel numerator constant coefficient = 1",
        "basis": ["1", "s=sqrt(-6)", "j=sqrt(-3)", "s*j"],
        "kernel_basis_components_mod_p": [list(c) for c in components],
        "numerator_basis_components_mod_p": [list(c) for c in components[:25]],
        "denominator_basis_components_mod_p": [list(c) for c in components[25:]],
        "support_counts": support_counts,
        "support_indices": support_indices,
        "source_files": {
            f"s{ss:+d}_j{sj:+d}": str(filename(modular_data, p, ss, sj))
            for ss, sj in EMBEDDINGS
        },
    }

    out = modular_data / f"q80_modj_p{p}_packet.json"
    out.write_text(json.dumps(packet, indent=2, sort_keys=True) + "\n")

    print(
        "Q80MODPACKET|"
        f"p={p}|rs={rs}|rj={rj}|"
        f"support_1={support_counts['1']}|"
        f"support_s={support_counts['s']}|"
        f"support_j={support_counts['j']}|"
        f"support_sj={support_counts['sj']}|"
        "roundtrip=PASS|status=PASS_GALOIS_PACKET",
        flush=True,
    )
    print(
        "Q80MODPACKET|"
        f"j_indices={','.join(map(str,support_indices['j'])) or '-'}|"
        f"sj_indices={','.join(map(str,support_indices['sj'])) or '-'}",
        flush=True,
    )
    print(f"Q80MODPACKET|out={out}", flush=True)


if __name__ == "__main__":
    main()
