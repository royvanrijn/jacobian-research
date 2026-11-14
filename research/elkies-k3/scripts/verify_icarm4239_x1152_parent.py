#!/usr/bin/env sage -python
"""Check the published ICARM4239 parent and its retrospective specialization.

The generic model export omits the exceptional parameter and points. Full
generic rank 17 is an external claim until a basis or Picard proof is imported.
"""
import argparse
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
from sage.all import EllipticCurve, PolynomialRing, QQ, ZZ

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT/"elliptic-curves/data/external/icarm_curve4239_20261113.json"


def check():
    R = PolynomialRing(QQ, "t")
    t = R.gen()
    L = 466548*t*t-828066*t+32416
    p = -478380*t*t-1076788*t+103008
    q = 82612*t*t-180205*t+81622
    D = -26666192028*t**4-91335649060*t**3+11566653808*t*t+239872553168*t-53374213248
    F = -47423161625*t**4+118525002574*t**3-104259187849*t*t+76235905664*t-24372268764
    B = p*F+q*D+p*q*(L-p-q)
    E = EllipticCurve(R.fraction_field(), [L, D+F, B, D*F, 0])
    # These are four exact sections, with NO independence claim.
    points = [E(0, 0), E(-D, 0), E(-F, 0), E(-p*q, -p*(F-p*q))]
    assert points[0]+points[1]+points[2] == E(0)
    data = json.loads(SOURCE.read_text())
    parameter = QQ(-2868768)/362209
    fiber = EllipticCurve(QQ, [a(parameter) for a in E.a_invariants()])
    target = EllipticCurve(QQ, [ZZ(a) for a in data["ainvs"]])
    iso = fiber.isomorphism_to(target)
    assert fiber.j_invariant() == target.j_invariant()
    # Inverse transport all witness points back to the published parent fiber.
    transported = [(~iso)(target(*map(QQ, P))) for P in data["points"]]
    assert len(transported) == 32
    A, C = R(-E.c4()/48), R(-E.c6()/864)
    spec = importlib.util.spec_from_file_location("record_mestre_gate", Path(__file__).with_name("pilot_hill1152_transfer.py"))
    pilot = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(pilot)
    gate = pilot.mestre_gate(A, C, {"mestre_prime_bound_exclusive": 2000})
    return {
        "schema": "icarm4239-parent-replay-v1", "status": "PASS",
        "source_url": "https://elliptic-rank.icarm.cloud/curve/4239.json",
        "source_sha256": sha256(SOURCE.read_bytes()).hexdigest(),
        "generic_model": {
            "a_invariants_low_to_high": [[str(v) for v in R(a).list()] for a in E.a_invariants()],
            "short_A_low_to_high": [str(v) for v in A.list()],
            "short_B_low_to_high": [str(v) for v in C.list()],
            "external_generic_rank_lower_bound": 17, "generic_rank_independently_certified": False,
            "exact_sections_verified": 4, "independence_asserted": False,
            "hill_parent_fibration_or_surface_identity": "UNKNOWN",
        },
        "retrospective_control": {
            "parameter": str(parameter), "fiber_isomorphism_urst": [str(v) for v in iso.tuple()],
            "public_points_transported": len(transported), "source_label": "KNOWN_RECORD_CONTROL",
        },
        "mestre_gate": gate,
        "boundary": "The record parent and Hill parent both have reported determinant1152; equality of this number is not a surface or fibration identification. No exceptional control selects prospective candidates.",
    }


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", required=True, type=Path)
    p.add_argument("--generic-export", type=Path)
    args = p.parse_args()
    result = check()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x") as f:
        json.dump(result, f, indent=2, sort_keys=True)
        f.write("\n")
    if args.generic_export:
        with args.generic_export.open("x") as f:
            json.dump(result["generic_model"], f, indent=2, sort_keys=True)
            f.write("\n")
    print(json.dumps({"status": "PASS", "public_points_transported": 32,
                      "mestre_gate": result["mestre_gate"]}, sort_keys=True))
