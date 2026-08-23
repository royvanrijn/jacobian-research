#!/usr/bin/env sage -python
"""Inspect rational-reconstruction progress of a q6-third Hensel artifact.

Example:
  sage -python ~/Downloads/inspect_h92_q6_third_hensel.sage \
    artifacts/local/elkies-k3/q6-third-hensel-p256.json
"""

import argparse
import json
from pathlib import Path
from sage.all import ZZ

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("input", type=Path)
args = parser.parse_args()

data = json.loads(args.input.read_text())
assert data["schema"] == "elkies-k3.h92-q6-third-hensel-lift.v1"
p = ZZ(data["prime"])
precision = int(data["precision"])
residues = [ZZ(v) for v in data["residues"]]
assert len(residues) == 138

blocks = (("Z", 0, 21), ("X", 21, 68), ("Y", 68, 138))

def digits(n):
    n = abs(ZZ(n))
    return 1 if not n else len(str(n))

levels = sorted(set(
    [k for k in (16, 32, 64, 128, 256, 512, 1024) if k <= precision]
    + [precision]
))

for k in levels:
    modulus = p**k
    recovered = []
    failures = []
    max_num_digits = 0
    max_den_digits = 0
    for i, residue in enumerate(residues):
        try:
            q = ZZ(residue % modulus).rational_reconstruction(modulus)
            recovered.append((i, q))
            max_num_digits = max(max_num_digits, digits(q.numerator()))
            max_den_digits = max(max_den_digits, digits(q.denominator()))
        except (ArithmeticError, ValueError):
            failures.append(i)

    parts = []
    for name, lo, hi in blocks:
        count = sum(lo <= i < hi for i, q in recovered)
        parts.append(f"{name}={count}/{hi-lo}")

    print(
        f"Q6THIRDRR|precision={k}|recovered={len(recovered)}/138|"
        + "|".join(parts)
        + f"|max_num_digits={max_num_digits}|max_den_digits={max_den_digits}|"
        f"failures={len(failures)}"
    )

# At full precision, report failed coefficient labels.
modulus = p**precision
failed_labels = []
for name, lo, hi in blocks:
    for i in range(lo, hi):
        try:
            ZZ(residues[i] % modulus).rational_reconstruction(modulus)
        except (ArithmeticError, ValueError):
            failed_labels.append(
                f"{name}[{i-lo}]"
            )
print(
    "Q6THIRDRR_FAILED|"
    + (",".join(failed_labels) if failed_labels else "none")
)
