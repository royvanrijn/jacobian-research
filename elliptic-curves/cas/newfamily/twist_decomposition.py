#!/usr/bin/env sage -python
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from sage.all import *

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(HERE))

from newfamily_rank11_minimal_common import build_finite_minimal_family

DEFAULT_ROOTS = (-47, -43, -31, 30, 45, 46)

def parse_roots(text):
    roots = tuple(ZZ(x.strip()) for x in text.split(","))
    if len(roots) != 6:
        raise ValueError("need six roots")
    return roots

def even_descend(poly):
    RT = poly.parent()
    T = RT.gen()
    RU = PolynomialRing(QQ, "U")
    U = RU.gen()
    out = RU(0)
    for i,c in enumerate(poly.list()):
        if c != 0 and i % 2:
            raise ArithmeticError(f"polynomial is not even: nonzero T^{i}")
        if c != 0:
            out += QQ(c)*U**(i//2)
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--roots", default=",".join(map(str, DEFAULT_ROOTS)))
    ap.add_argument("--output")
    args = ap.parse_args()

    roots = parse_roots(args.roots)
    D = build_finite_minimal_family(roots)
    A, B, Delta = D["Amin"], D["Bmin"], D["Deltamin"]

    A0 = even_descend(A)
    B0 = even_descend(B)
    Delta0 = even_descend(Delta)

    if (A0.degree(), B0.degree(), Delta0.degree()) != (4,6,10):
        raise ArithmeticError(
            f"unexpected rational-surface degrees "
            f"{(A0.degree(), B0.degree(), Delta0.degree())}"
        )

    if gcd(Delta0, Delta0.derivative()).degree() != 0:
        raise ArithmeticError("base finite discriminant is not squarefree")
    if gcd(Delta0, A0).degree() != 0:
        raise ArithmeticError("base finite discriminant meets c4=0")

    U = A0.parent().gen()
    Atw = U**2 * A0
    Btw = U**3 * B0
    Dtw = -16*(4*Atw**3 + 27*Btw**2)

    result = {
        "roots": [int(x) for x in roots],
        "K3": {
            "A_T": str(A),
            "B_T": str(B),
            "Delta_T": str(Delta),
            "degrees": [int(A.degree()), int(B.degree()), int(Delta.degree())],
        },
        "quadratic_base_change": "U=T^2",
        "rational_surface_E0": {
            "equation": f"y^2=x^3+({A0})*x+({B0})",
            "A_U": str(A0),
            "B_U": str(B0),
            "Delta_U": str(Delta0),
            "degrees": [int(A0.degree()), int(B0.degree()), int(Delta0.degree())],
            "finite_fibers": "10 I1",
            "infinity_fiber": "I2",
            "geometric_Picard_rank": 10,
            "trivial_lattice_rank": 3,
            "geometric_MW_rank": 7,
        },
        "quadratic_twist_by_U": {
            "equation": f"y^2=x^3+({Atw})*x+({Btw})",
            "A_U": str(Atw),
            "B_U": str(Btw),
            "Delta_U": str(Dtw),
            "relation": (
                "rank MW(K3 over Qbar(T)) = 7 + "
                "rank MW(E0^(U) over Qbar(U))"
            ),
            "current_implication": (
                "if K3 MW rank is 11 or 12, twist MW rank is 4 or 5"
            ),
        },
    }

    print("BASE_CHANGE|U=T^2")
    print(f"E0|degA={A0.degree()}|degB={B0.degree()}|degDelta={Delta0.degree()}")
    print("E0|fibers=10I1+I2|rho=10|trivial=3|MW=7")
    print("TWIST|K3_MW=7+twist_MW|current_twist_rank=4_or_5")
    print()
    print("A0 =", A0)
    print("B0 =", B0)

    if args.output:
        path = Path(args.output)
    else:
        path = (
            REPO / "artifacts/local/elliptic-curves/newfamily"
            / "quadratic_base_change.json"
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(f"OUTPUT|{path}")

if __name__ == "__main__":
    main()
