#!/usr/bin/env sage -python
"""Replay the 11 exact hidden generic sections on the finite-minimal family."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

from sage.all import FractionField, PolynomialRing, QQ

ROOTS = (-47, -43, -31, 30, 45, 46)


def load_builder():
    candidates = [
        Path("elliptic-curves/cas/newfamily/newfamily_rank11_minimal_common.py"),
        Path("elliptic-curves/cas/newfamily/archive/newfamily_rank11_minimal_common.py"),
        Path("/tmp/newfamily_rank11_minimal_common.py"),
    ]
    for path in candidates:
        if path.exists():
            sys.path.insert(0, str(path.parent.resolve()))
            from newfamily_rank11_minimal_common import build_finite_minimal_family
            return build_finite_minimal_family
    raise SystemExit("missing newfamily_rank11_minimal_common.py")


def q(value):
    if isinstance(value, tuple):
        return QQ(value[0]) / QQ(value[1])
    return QQ(value)


def polynomial(ring, coefficients):
    return ring([q(value) for value in coefficients])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--data",
        default="elliptic-curves/cas/newfamily/hidden_sections_data.py",
    )
    args = parser.parse_args()

    data_path = Path(args.data)
    if not data_path.exists():
        raise SystemExit(
            f"missing {data_path}; generate it with export_hidden_sections.py first"
        )
    sys.path.insert(0, str(data_path.parent.resolve()))
    module_name = data_path.stem
    data = __import__(module_name)

    build = load_builder()
    family = build(ROOTS)

    R = PolynomialRing(QQ, "T")
    K = FractionField(R)
    A = K(family["Amin"])
    B = K(family["Bmin"])

    if len(data.SECTIONS) != 11:
        raise AssertionError("expected 11 section records")

    print("NEWFAMILY_HIDDEN_SECTIONS_REPLAY_V1")
    print(f"A_degree={family['Amin'].degree()} B_degree={family['Bmin'].degree()}")

    for expected, record in enumerate(data.SECTIONS):
        if record["index"] != expected:
            raise AssertionError("section records are out of order")
        xn = polynomial(R, record["x_num"])
        xd = polynomial(R, record["x_den"])
        yn = polynomial(R, record["y_num"])
        yd = polynomial(R, record["y_den"])
        x = K(xn) / K(xd)
        y = K(yn) / K(yd)
        identity = y**2 - (x**3 + A * x + B)
        ok = identity == 0
        print(
            "U%d x=(%d,%d) y=(%d,%d) identity=%s" %
            (expected, xn.degree(), xd.degree(), yn.degree(), yd.degree(), ok)
        )
        if not ok:
            raise AssertionError(f"section U{expected} failed the curve identity")

    print("VERIFIED_GENERIC_SECTIONS=11/11")
    print("DONE")


if __name__ == "__main__":
    main()
