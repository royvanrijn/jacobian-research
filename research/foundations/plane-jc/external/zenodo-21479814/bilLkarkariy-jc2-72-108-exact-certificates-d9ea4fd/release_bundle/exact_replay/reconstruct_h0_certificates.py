#!/usr/bin/env python3
"""Serialize and independently check the two pre-division h=0 certificates."""

from __future__ import annotations

import hashlib
import pickle
import re
import subprocess
from pathlib import Path

import sympy as sp
from flint import fmpq_poly

import exact_core as ec


ROOT = Path(__file__).resolve().parent
BRANCH1_SYSTEM = ROOT / "case1_branch1_after_w.pkl"
BRANCH2_SYSTEM = ROOT / "case1_branch2_after_w.pkl"
LIFT_SCRIPT = ROOT / "case1_branch1_h0_lift.sing"
BRANCH1_CERTIFICATE = ROOT / "h0_branch1_exact_certificate.pkl"
BRANCH2_CERTIFICATE = ROOT / "h0_exact_certificate.pkl"
ROW_SCALES = (-1, 1, 1, -1, -1, -1, -1)
RECORD = re.compile(r"^C\|(\d+)\|(\d+)\|(\d+)\|(\d+)\|(\d+)\|\((.*)\)$")


def decode_k(serialized):
    coefficients = [ec.fmpq(0) for _ in range(35)]
    for degree, (numerator, denominator) in serialized.items():
        coefficients[int(degree)] = ec.fmpq(int(numerator), int(denominator))
    return ec.K(fmpq_poly(coefficients))


def encode_k(value):
    return {
        i: (int(value.p[i].p), int(value.p[i].q))
        for i in range(len(value.p))
        if value.p[i]
    }


def decode_system(path):
    raw = pickle.loads(path.read_bytes())
    return [
        {tuple(m): decode_k(c) for m, c in equation.items()}
        for equation in raw
    ]


def specialize_h0(system):
    return [
        {m[1:]: c for m, c in equation.items() if m[0] == 0}
        for equation in system
    ]


def poly_add(left, right):
    out = dict(left)
    for monomial, coefficient in right.items():
        value = out.get(monomial, ec.K(0)) + coefficient
        if value:
            out[monomial] = value
        else:
            out.pop(monomial, None)
    return out


def poly_scale(poly, scalar):
    return {m: c * scalar for m, c in poly.items() if c * scalar}


def poly_mul(left, right):
    out = {}
    for m, c in left.items():
        for n, d in right.items():
            monomial = tuple(a + b for a, b in zip(m, n))
            out[monomial] = out.get(monomial, ec.K(0)) + c * d
    return {m: c for m, c in out.items() if c}


def sign_substitution(poly, signs):
    out = {}
    for monomial, coefficient in poly.items():
        parity = sum(e for e, sign in zip(monomial, signs) if sign == -1)
        out[monomial] = coefficient * (-1 if parity % 2 else 1)
    return out


def verify_identity(generators, multipliers):
    residual = {(0, 0, 0): ec.K(-1)}
    for generator, multiplier in zip(generators, multipliers):
        residual = poly_add(residual, poly_mul(generator, multiplier))
    if residual:
        monomial = next(iter(residual))
        raise RuntimeError(f"h=0 identity has a nonzero residual at {monomial}")


def run_and_parse_lift():
    completed = subprocess.run(
        ["Singular", "-q", str(LIFT_SCRIPT)],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=True,
    )
    output = completed.stdout
    if "LIFT_BEGIN" not in output or "LIFT_END" not in output:
        raise RuntimeError("Singular did not emit a complete h=0 lift")
    multipliers = [dict() for _ in range(7)]
    for line in output.splitlines():
        match = RECORD.match(line)
        if not match:
            continue
        equation, h_exp, u1_exp, u2_exp, u3_exp = map(int, match.groups()[:5])
        if h_exp:
            raise RuntimeError("specialized h=0 multiplier unexpectedly contains h")
        expression = sp.sympify(
            match.group(6).replace("^", "**"), locals={"u": ec.U}
        )
        coefficient = ec.K(ec.qpoly(expression))
        multipliers[equation - 1][(u1_exp, u2_exp, u3_exp)] = coefficient
    if not any(multipliers) or len(multipliers) != 7:
        raise RuntimeError("no Singular lift records were parsed")
    return output, multipliers


def encode_poly(poly):
    return {m: encode_k(c) for m, c in poly.items()}


def payload(branch, source, singular_hash, generators, multipliers):
    return {
        "format": "jc2-h0-nullstellensatz-v1",
        "branch": branch,
        "field_minpoly_coefficients": [
            (int(ec.MOD[i].p), int(ec.MOD[i].q)) for i in range(len(ec.MOD))
        ],
        "source_system": source.name,
        "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "singular_lift_stdout_sha256": singular_hash,
        "generators": [encode_poly(q) for q in generators],
        "multipliers": [encode_poly(q) for q in multipliers],
    }


def main():
    branch1_full = decode_system(BRANCH1_SYSTEM)
    branch2_full = decode_system(BRANCH2_SYSTEM)
    for i, (left, right, scale) in enumerate(zip(branch1_full, branch2_full, ROW_SCALES)):
        transported = poly_scale(sign_substitution(left, (1, -1, -1, 1)), scale)
        if transported != right:
            raise RuntimeError(f"pre-division branch symmetry fails on equation {i}")
    print("PRE_DIVISION_SYSTEM_SYMMETRY_PASS")

    branch1 = specialize_h0(branch1_full)
    branch2 = specialize_h0(branch2_full)
    singular_output, multipliers1 = run_and_parse_lift()
    verify_identity(branch1, multipliers1)

    multipliers2 = [
        poly_scale(sign_substitution(multiplier, (-1, -1, 1)), scale)
        for multiplier, scale in zip(multipliers1, ROW_SCALES)
    ]
    verify_identity(branch2, multipliers2)

    singular_hash = hashlib.sha256(singular_output.encode()).hexdigest()
    BRANCH1_CERTIFICATE.write_bytes(
        pickle.dumps(
            payload("s=c", BRANCH1_SYSTEM, singular_hash, branch1, multipliers1),
            protocol=4,
        )
    )
    BRANCH2_CERTIFICATE.write_bytes(
        pickle.dumps(
            payload("s=-c", BRANCH2_SYSTEM, singular_hash, branch2, multipliers2),
            protocol=4,
        )
    )
    print("H0_BRANCH1_EXACT_IDENTITY_PASS")
    print("H0_BRANCH2_EXACT_IDENTITY_PASS")
    print(f"wrote {BRANCH1_CERTIFICATE} ({BRANCH1_CERTIFICATE.stat().st_size} bytes)")
    print(f"wrote {BRANCH2_CERTIFICATE} ({BRANCH2_CERTIFICATE.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
