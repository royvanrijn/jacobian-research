#!/usr/bin/env sage -python
"""Measure cross-precision stability of rational reconstruction candidates.

A rational reconstruction at one modulus is only a candidate.  A genuine
small rational coefficient should reconstruct to the SAME value once the
precision is large enough.  Random near-sqrt(M) reconstructions typically
change when the modulus is doubled.

Example:
  sage -python ~/Downloads/check_h92_q6_third_rr_stability.sage \
    artifacts/local/elkies-k3/q6-third-hensel-p1024.json
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

blocks = (("Z",0,21),("X",21,68),("Y",68,138))

def lab(i):
    for name,lo,hi in blocks:
        if lo <= i < hi:
            return f"{name}[{i-lo}]"
    raise AssertionError(i)

def rr(residue, k):
    M = p**k
    try:
        return ZZ(residue % M).rational_reconstruction(M)
    except (ArithmeticError, ValueError):
        return None

levels = [k for k in (64,128,256,512,1024,2048,4096) if k <= precision]
values = {k:[rr(r,k) for r in residues] for k in levels}

for k in levels:
    got = [i for i,q in enumerate(values[k]) if q is not None]
    parts=[]
    for name,lo,hi in blocks:
        parts.append(f"{name}={sum(lo<=i<hi for i in got)}/{hi-lo}")
    print(f"Q6RRSTABLE_LEVEL|precision={k}|recovered={len(got)}/138|"+"|".join(parts))

for lo,hi in zip(levels, levels[1:]):
    same=[]; changed=[]; new=[]; lost=[]
    for i,(a,b) in enumerate(zip(values[lo], values[hi])):
        if a is not None and b is not None:
            (same if a==b else changed).append(i)
        elif a is None and b is not None:
            new.append(i)
        elif a is not None and b is None:
            lost.append(i)
    print(
        f"Q6RRSTABLE_PAIR|{lo}->{hi}|same={len(same)}|changed={len(changed)}|"
        f"new={len(new)}|lost={len(lost)}"
    )
    if changed:
        print("Q6RRSTABLE_CHANGED|"+",".join(lab(i) for i in changed))

# Strong candidates: unchanged on the last available doubling.
if len(levels) >= 2:
    lo,hi = levels[-2],levels[-1]
    stable=[i for i,(a,b) in enumerate(zip(values[lo],values[hi])) if a is not None and a==b]
    parts=[]
    for name,blo,bhi in blocks:
        parts.append(f"{name}={sum(blo<=i<bhi for i in stable)}/{bhi-blo}")
    print(
        f"Q6RRSTABLE_RESULT|pair={lo}->{hi}|stable={len(stable)}/138|"+"|".join(parts)
    )
    unstable=[i for i in range(138) if i not in stable]
    print("Q6RRSTABLE_UNSTABLE|"+(",".join(lab(i) for i in unstable) if unstable else "none"))
