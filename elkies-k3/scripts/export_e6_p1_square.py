#!/usr/bin/env python3
from pathlib import Path
import argparse, subprocess, random, sys

ap = argparse.ArgumentParser()
ap.add_argument("--p", type=int, default=101)
ap.add_argument("--seed", type=int, default=1)
ap.add_argument("--slices", type=int, default=3)
ap.add_argument("--out", required=True)
args = ap.parse_args()

out = Path(args.out)
out.parent.mkdir(parents=True, exist_ok=True)

raw = out.with_name(out.stem + "-raw.ms")

# First let the known-good Sage exporter construct the 15-equation system.
cmd = [
    "sage",
    "elkies-k3/scripts/export_e6_p1_sliced.sage",
    "--p", str(args.p),
    "--seed", str(args.seed),
    "--slices", str(args.slices),
    "--out", str(raw),
]

r = subprocess.run(cmd, text=True, capture_output=True)

sys.stdout.write(r.stdout)
sys.stderr.write(r.stderr)

if r.returncode:
    raise SystemExit(r.returncode)

# Read msolve file. Important: equations are one per physical line.
lines = raw.read_text().splitlines()

if len(lines) < 3:
    raise SystemExit(
        f"E6SQUARE|error=raw_file_too_short|lines={len(lines)}|file={raw}"
    )

varline = lines[0].strip()
p = int(lines[1].strip())

vars_ = [x.strip() for x in varline.split(",") if x.strip()]
n = len(vars_)

eqs = []
for line in lines[2:]:
    e = line.strip()
    if not e:
        continue
    if e.endswith(","):
        e = e[:-1].rstrip()
    if e:
        eqs.append(e)

m = len(eqs)

print(
    f"E6SQUARE|stage=input|vars={n}|source_eqs={m}|p={p}|raw={raw}",
    flush=True,
)

if m < n:
    raise SystemExit(
        f"E6SQUARE|error=underdetermined_input|vars={n}|eqs={m}"
    )

# Construct n generic linear combinations over GF(p).
random.seed(args.seed ^ 0xE6E6A17)

C = []
for i in range(n):
    row = [random.randrange(p) for _ in range(m)]

    if all(c == 0 for c in row):
        row[i % m] = 1

    C.append(row)

# Ensure every original equation participates at least once.
for j in range(m):
    if all(C[i][j] == 0 for i in range(n)):
        C[j % n][j] = 1

combined = []

for i, row in enumerate(C):
    parts = []

    for c, e in zip(row, eqs):
        c %= p

        if c == 0:
            continue

        if c == 1:
            parts.append(f"({e})")
        else:
            parts.append(f"{c}*({e})")

    if not parts:
        raise RuntimeError(f"empty combination row {i}")

    combined.append("+".join(parts))

with out.open("w") as h:
    h.write(varline + "\n")
    h.write(str(p) + "\n")

    for i, e in enumerate(combined):
        h.write(e)
        if i + 1 < len(combined):
            h.write(",\n")
        else:
            h.write("\n")

meta = out.with_suffix(".square.meta.txt")

with meta.open("w") as h:
    h.write(f"p={p}\n")
    h.write(f"seed={args.seed}\n")
    h.write(f"vars={n}\n")
    h.write(f"source_eqs={m}\n")
    h.write(f"square_eqs={len(combined)}\n")
    h.write(f"raw={raw}\n")
    h.write("matrix=\n")

    for row in C:
        h.write(" ".join(map(str, row)) + "\n")

print(
    f"E6SQUARE|stage=done|vars={n}|source_eqs={m}|square_eqs={len(combined)}|out={out}",
    flush=True,
)
