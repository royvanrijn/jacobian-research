#!/usr/bin/env sage
"""One bounded equation-derived arithmetic profile worker.

Run under ``sage -python``.  All arithmetic is exact except logarithms used only
as descriptive size summaries.  Class-group mode is explicitly provisional
(`proof=False`) unless a later independent certification is supplied.
"""
from __future__ import annotations

import argparse
import json
import math
import re
import resource
from pathlib import Path
import sys

from sage.all import EllipticCurve, NumberField, PolynomialRing, QQ, ZZ


def stable(obj):
    return json.dumps(obj, sort_keys=True, indent=2, allow_nan=False) + "\n"


def write(path, obj):
    Path(path).write_text(stable(obj))


def log2int(n):
    n = abs(ZZ(n))
    if n == 0:
        return None
    # Accurate enough for a display statistic while avoiding float overflow.
    bits = n.nbits()
    if bits <= 53:
        return math.log2(int(n))
    shift = bits - 53
    top = int(n >> shift)
    return math.log2(top) + shift


def primitive_poly(poly):
    den = ZZ(1)
    for c in poly.list():
        den = den.lcm(ZZ(QQ(c).denominator()))
    coeffs = [ZZ(QQ(c) * den) for c in poly.list()]
    g = ZZ(0)
    for c in coeffs:
        g = g.gcd(c)
    if g:
        coeffs = [c // g for c in coeffs]
    if coeffs and coeffs[-1] < 0:
        coeffs = [-c for c in coeffs]
    R = poly.parent()
    return R(coeffs)


def kodaira_kind(symbol):
    s = str(symbol).replace(" ", "")
    # Sage prints I5, I0*, II, III, IV, etc.  Plain In is multiplicative;
    # starred I_n and all named symbols are additive.
    if "*" not in s and re.fullmatch(r"I_?\d+", s):
        return "multiplicative"
    return "additive"


def construct(row):
    E = EllipticCurve(QQ, [ZZ(x) for x in row["ainvs"]])
    try:
        Emin = E.global_minimal_model()
    except Exception:
        Emin = E.minimal_model()
    return E, Emin


def base_profile(row):
    E, Emin = construct(row)
    b2, b4, b6, _b8 = Emin.b_invariants()
    R = PolynomialRing(QQ, "x")
    x = R.gen()
    cubic = primitive_poly(4*x**3 + b2*x**2 + 2*b4*x + b6)
    roots = cubic.roots(QQ, multiplicities=False)
    Delta = ZZ(Emin.discriminant())
    c4 = ZZ(Emin.c_invariants()[0])
    c6 = ZZ(Emin.c_invariants()[1])
    return {
        "schema": "elliptic-curves.wide-arithmetic-base.v1",
        "status": "PASS",
        "curve_key": row["curve_key"],
        "minimal_ainvs": [str(v) for v in Emin.ainvs()],
        "c4": str(c4),
        "c6": str(c6),
        "minimal_discriminant": str(Delta),
        "delta_sign": -1 if Delta < 0 else 1,
        "log2_abs_minimal_discriminant": log2int(Delta),
        "two_division_cubic": [str(c) for c in cubic.list()],
        "two_division_cubic_discriminant": str(ZZ(cubic.discriminant())),
        "log2_abs_polynomial_discriminant": log2int(cubic.discriminant()),
        "rational_2torsion_rank": 0 if len(roots) == 0 else (1 if len(roots) == 1 else 2),
        "cubic_irreducible": bool(cubic.is_irreducible()),
    }


def local_profile(row):
    E, Emin = construct(row)
    b2, b4, b6, _b8 = Emin.b_invariants()
    R = PolynomialRing(QQ, "x")
    x = R.gen()
    cubic = primitive_poly(4*x**3 + b2*x**2 + 2*b4*x + b6)
    Delta = ZZ(Emin.discriminant())
    roots = cubic.roots(QQ, multiplicities=False)
    payload = {
        "schema": "elliptic-curves.wide-arithmetic-local.v1",
        "curve_key": row["curve_key"],
        "status": "PASS",
        "rational_2torsion_rank": 0 if len(roots) == 0 else (1 if len(roots) == 1 else 2),
        "bk_applicable": len(roots) == 0 and cubic.is_irreducible(),
    }
    if not cubic.is_irreducible():
        payload.update({
            "status": "PASS_BK_NOT_APPLICABLE",
            "reason": "2-division cubic reducible over Q",
        })
        return payload

    K = NumberField(cubic.monic(), "a")
    OK = K.maximal_order()
    field_disc = ZZ(K.discriminant())
    field_disc_fac = list(ZZ(abs(field_disc)).factor())
    payload.update({
        "field_discriminant": str(field_disc),
        "log2_abs_field_discriminant": log2int(field_disc),
        "field_signature": list(map(int, K.signature())),
        "field_ramified_primes": [str(p) for p, _e in field_disc_fac],
        "field_ramified_prime_count": len(field_disc_fac),
    })

    delta_fac = list(ZZ(abs(Delta)).factor())
    locals_out = []
    phi_m = []
    phi_a = []
    additive_split_counts = []
    conductor = ZZ(1)
    root_factors = []
    for p, valuation in delta_fac:
        p = ZZ(p)
        ld = Emin.local_data(p)
        symbol = str(ld.kodaira_symbol())
        kind = kodaira_kind(symbol)
        vdelta = int(ld.discriminant_valuation())
        fcond = int(ld.conductor_valuation())
        conductor *= p**fcond
        try:
            w = int(Emin.root_number(p))
        except Exception:
            w = None
        if w is not None:
            root_factors.append(w)
        nprimes = None
        if kind == "multiplicative":
            if vdelta % 2 == 0:
                phi_m.append(int(p))
        else:
            phi_a.append(int(p))
            # Number of primes of K above p, as in Brumer-Kramer.
            try:
                nprimes = len(K.primes_above(p))
            except Exception:
                nprimes = len(OK.ideal(p).factor())
            additive_split_counts.append(int(nprimes))
        locals_out.append({
            "p": str(p), "v_delta": vdelta, "kodaira": symbol,
            "kind": kind, "conductor_exponent": fcond,
            "root_number": w, "cubic_prime_count": nprimes,
        })

    u = 1 if Delta < 0 else 2
    nterm = len(phi_m) + sum(n - 1 for n in additive_split_counts)
    root_number = None
    try:
        root_number = int(Emin.root_number())
    except Exception:
        if len(root_factors) == len(delta_fac):
            root_number = -1
            for w in root_factors:
                root_number *= w
    payload.update({
        "minimal_discriminant_factorization": [[str(p), int(e)] for p, e in delta_fac],
        "local_data": locals_out,
        "phi_m": list(map(str, phi_m)),
        "phi_m_count": len(phi_m),
        "phi_a": list(map(str, phi_a)),
        "phi_a_count": len(phi_a),
        "additive_split_counts": additive_split_counts,
        "bk_u": u,
        "bk_n": nterm,
        "bk_local_term": u + nterm,
        "conductor": str(conductor),
        "log2_conductor": log2int(conductor),
        "root_number": root_number,
    })
    return payload


def class_profile(row):
    E, Emin = construct(row)
    b2, b4, b6, _b8 = Emin.b_invariants()
    R = PolynomialRing(QQ, "x")
    x = R.gen()
    cubic = primitive_poly(4*x**3 + b2*x**2 + 2*b4*x + b6)
    if not cubic.is_irreducible():
        return {
            "schema": "elliptic-curves.wide-arithmetic-class-probe.v1",
            "curve_key": row["curve_key"], "status": "NOT_APPLICABLE_REDUCIBLE_2DIVISION",
        }
    K = NumberField(cubic.monic(), "a")
    # Explicitly conditional/provisional.  This is never promoted to an
    # unconditional class-group or rank claim by the controller.
    C = K.class_group(proof=False)
    invariants = [ZZ(v) for v in C.invariants()]
    g = sum(1 for v in invariants if v % 2 == 0)
    return {
        "schema": "elliptic-curves.wide-arithmetic-class-probe.v1",
        "curve_key": row["curve_key"],
        "status": "PASS_PROVISIONAL_GRH",
        "class_group_invariants": [str(v) for v in invariants],
        "class_2rank": int(g),
        "boundary": "proof=False / provisional class-group computation; no unconditional certification",
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--mode", choices=("base", "local", "class"), required=True)
    ap.add_argument("--memory-gb", type=float, default=0.0)
    args = ap.parse_args()
    if args.memory_gb > 0:
        limit = int(args.memory_gb * (1024 ** 3))
        try:
            resource.setrlimit(resource.RLIMIT_AS, (limit, limit))
        except Exception:
            pass
    row = json.loads(Path(args.input).read_text())
    try:
        if args.mode == "base":
            result = base_profile(row)
        elif args.mode == "local":
            result = local_profile(row)
        else:
            result = class_profile(row)
    except Exception as exc:
        result = {
            "schema": f"elliptic-curves.wide-arithmetic-{args.mode}-error.v1",
            "curve_key": row.get("curve_key"),
            "status": "UNKNOWN_ARITHMETIC_FAILURE",
            "error_type": type(exc).__name__,
            "error": str(exc),
        }
    write(args.output, result)
    return 0 if result.get("status", "").startswith(("PASS", "NOT_APPLICABLE")) else 2


if __name__ == "__main__":
    raise SystemExit(main())
