#!/usr/bin/env python3
"""Replay the three serialized characteristic-zero certificates exactly."""

from __future__ import annotations

import hashlib
import pickle
from pathlib import Path

from flint import fmpq_poly

import case2_exact_generate as c2
import exact_core as ec


ROOT = Path(__file__).resolve().parent


def decode_k(serialized):
    coefficients = [ec.fmpq(0) for _ in range(35)]
    for degree, (numerator, denominator) in serialized.items():
        coefficients[int(degree)] = ec.fmpq(int(numerator), int(denominator))
    return ec.K(fmpq_poly(coefficients))


def decode_pp(serialized):
    return ec.PP({tuple(m): decode_k(c) for m, c in serialized.items()})


def encode_k(value):
    return {
        i: (int(value.p[i].p), int(value.p[i].q))
        for i in range(len(value.p))
        if value.p[i]
    }


def verify_field(payload):
    expected = [
        (int(ec.MOD[i].p), int(ec.MOD[i].q)) for i in range(len(ec.MOD))
    ]
    if payload["field_minpoly_coefficients"] != expected:
        raise RuntimeError("certificate uses an unexpected degree-35 number field")


def verify_case2():
    path = ROOT / "case2_exact_certificate.pkl"
    payload = pickle.loads(path.read_bytes())
    if payload.get("format") != "jc2-case2-nullstellensatz-v1":
        raise RuntimeError("unsupported Case-2 certificate format")
    verify_field(payload)
    regenerated = [c2.final[i] for i in payload["generator_indices"]]
    saved_generators = [decode_pp(q) for q in payload["generators"]]
    if any(a.d != b.d for a, b in zip(regenerated, saved_generators)):
        raise RuntimeError("Case-2 saved generators differ from exact regeneration")
    multipliers = [decode_pp(q) for q in payload["multipliers"]]
    residual = ec.PP.const(-1)
    for generator, multiplier in zip(saved_generators, multipliers):
        residual += generator * multiplier
    if residual:
        raise RuntimeError("Case-2 serialized identity has a nonzero residual")
    print("CASE2_SERIALIZED_EXACT_PASS")


def poly_add(left, right):
    out = dict(left)
    for monomial, coefficient in right.items():
        value = out.get(monomial, ec.K(0)) + coefficient
        if value:
            out[monomial] = value
        else:
            out.pop(monomial, None)
    return out


def poly_mul(left, right):
    out = {}
    for m, c in left.items():
        for n, d in right.items():
            monomial = tuple(a + b for a, b in zip(m, n))
            out[monomial] = out.get(monomial, ec.K(0)) + c * d
    return {m: c for m, c in out.items() if c}


def decode_poly(serialized):
    return {tuple(m): decode_k(c) for m, c in serialized.items()}


def specialize_machine_system(path):
    raw = pickle.loads(path.read_bytes())
    return [
        {
            tuple(m[1:]): decode_k(c)
            for m, c in equation.items()
            if m[0] == 0
        }
        for equation in raw
    ]


def verify_h0(filename, expected_branch):
    certificate_path = ROOT / filename
    payload = pickle.loads(certificate_path.read_bytes())
    if payload.get("format") != "jc2-h0-nullstellensatz-v1":
        raise RuntimeError(f"unsupported certificate format in {filename}")
    if payload.get("branch") != expected_branch:
        raise RuntimeError(f"unexpected branch label in {filename}")
    verify_field(payload)
    source = ROOT / payload["source_system"]
    if hashlib.sha256(source.read_bytes()).hexdigest() != payload["source_sha256"]:
        raise RuntimeError(f"source hash mismatch for {filename}")
    regenerated = specialize_machine_system(source)
    generators = [decode_poly(q) for q in payload["generators"]]
    if regenerated != generators:
        raise RuntimeError(f"h=0 generators differ from regeneration in {filename}")
    multipliers = [decode_poly(q) for q in payload["multipliers"]]
    residual = {(0, 0, 0): ec.K(-1)}
    for generator, multiplier in zip(generators, multipliers):
        residual = poly_add(residual, poly_mul(generator, multiplier))
    if residual:
        raise RuntimeError(f"h=0 serialized identity fails in {filename}")
    print(f"{expected_branch} H0_SERIALIZED_EXACT_PASS")


def main():
    verify_case2()
    verify_h0("h0_branch1_exact_certificate.pkl", "s=c")
    verify_h0("h0_exact_certificate.pkl", "s=-c")
    print("ALL_SERIALIZED_EXACT_CERTIFICATES_PASS")


if __name__ == "__main__":
    main()
