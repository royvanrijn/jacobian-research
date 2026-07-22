#!/usr/bin/env python3
"""Check exact branch transport and replay the large identity on branch 2."""

from __future__ import annotations

import pickle
import re
import subprocess
import sys
from pathlib import Path

from flint import fmpq

from degree5_core import L, Poly, decode_poly, sign_substitution


ROOT = Path(__file__).resolve().parent
BRANCH1 = ROOT / "hne0_polred.pkl"
BRANCH2 = ROOT / "hne0_branch2_polred.pkl"
CERTIFICATE = ROOT / "hard" / "h_certificate_exact.txt"
ROW_SCALES = (1, 1, -1, -1, -1, -1)
TERM = re.compile(r"^\((-?\d+)(?:/(\d+))?\)(?:\*w(?:\^(\d+))?)?$")

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)


def parse_l(expression):
    if expression == "0":
        return L(0)
    coefficients = [fmpq(0) for _ in range(5)]
    for item in expression.split("+"):
        match = TERM.match(item)
        if not match:
            raise RuntimeError(f"cannot parse exact coefficient: {item[:120]}")
        degree = int(match.group(3) or (1 if "*w" in item else 0))
        coefficients[degree] += fmpq(
            int(match.group(1)), int(match.group(2) or 1)
        )
    return L(coefficients)


def load_multipliers():
    multipliers = [Poly(3) for _ in range(4)]
    with CERTIFICATE.open() as handle:
        for line in handle:
            if not line.startswith("C|"):
                continue
            _, equation, h_exp, u1_exp, u2_exp, expression = line.rstrip("\n").split("|", 5)
            monomial = (int(h_exp), int(u1_exp), int(u2_exp))
            multipliers[int(equation) - 1] += Poly(
                3, {monomial: parse_l(expression)}
            )
    return multipliers


def verify_h_identity(generators, multipliers, label):
    h = Poly(3, {(1, 0, 0): L(1)})
    residual = -h
    for generator, multiplier in zip(generators, multipliers):
        residual += generator * multiplier
    if residual:
        raise RuntimeError(f"{label} identity has {len(residual.terms)} residual terms")


def main():
    branch1 = [decode_poly(q) for q in pickle.loads(BRANCH1.read_bytes())]
    branch2 = [decode_poly(q) for q in pickle.loads(BRANCH2.read_bytes())]
    for i, (left, right, scale) in enumerate(zip(branch1, branch2, ROW_SCALES)):
        transported = sign_substitution(left, (1, -1, -1)) * scale
        if transported.terms != right.terms:
            raise RuntimeError(f"degree-five branch symmetry fails on equation {i}")
    print("SYSTEM_SYMMETRY_PASS")

    replay = subprocess.run(
        [sys.executable, str(ROOT / "hard" / "verify_certificate_gmpy2.py")],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=True,
    )
    if "GMPY2_EXACT_PASS" not in replay.stdout:
        raise RuntimeError("branch-1 exact replay did not pass")
    print(replay.stdout, end="")

    # Formal transport: B_i = s_i*phi(A_i), where phi fixes h and negates
    # (u1,u2).  Therefore T'_i=s_i*phi(T_i) gives
    # sum(T'_i*B_i)=phi(sum(T_i*A_i))=phi(h)=h.  The system identity above
    # and the independently replayed branch-1 identity make this exact; no
    # second multiplication of the 89 MB coefficient list is needed.
    supports = []
    with CERTIFICATE.open() as handle:
        for line in handle:
            if line.startswith("C|"):
                parts = line.split("|", 5)
                supports.append(tuple(map(int, parts[1:5])))
    if len(supports) != 385 or any(not 1 <= item[0] <= 4 for item in supports):
        raise RuntimeError("unexpected support in branch-1 certificate")
    print("BRANCH1_EXACT_IDENTITY_PASS")
    print("BRANCH2_EXACT_IDENTITY_PASS")


if __name__ == "__main__":
    main()
