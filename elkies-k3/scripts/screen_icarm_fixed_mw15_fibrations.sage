#!/usr/bin/env sage-python
"""Screen exact fixed-corridor MW15 fibrations against pinned ICARM targets.

This is the first fail-closed layer after the complete A1/MW16 atlas.  It
currently consumes the exact fixed-corridor 2A1/MW15 Jacobian obtained from
the source-identified q12/orbit5867 rootless endpoint.  For each requested
target it first applies projective modular j-preimage tests at a declared
prime chain.  Only modular survivors are factored over QQ, and every rational
j-match is checked for QQ-isomorphism rather than accepted up to twist.

The certificate is intentionally bounded by its explicit model list.  It is
not an exhaustion of all A2/2A1 fibrations on the K3.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing, QQ, ZZ
from sage.env import SAGE_VERSION


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_MODEL = ROOT / "artifacts/local/elkies-k3/fixed-reverse-2a1-rr-qq.json"
DEFAULT_SNAPSHOT = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-norm12-icarm-database-sweep-v2.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-icarm-curve302-273-fixed-2a1-mw15-screen-v1.json"
)
DEFAULT_PRIMES = (
    1009, 1013, 1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061,
    1063, 1069, 1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123,
    1129, 1151, 1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213,
    1217, 1223, 1229, 1231,
)


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path.resolve())


def normalized_ainvs(values):
    if len(values) == 2:
        values = [0, 0, 0, *values]
    if len(values) != 5:
        raise ValueError("target curve has neither two nor five a-invariants")
    return tuple(Fraction(str(value)) for value in values)


def invariants(values):
    a1, a2, a3, a4, a6 = normalized_ainvs(values)
    b2 = a1**2 + 4 * a2
    b4 = 2 * a4 + a1 * a3
    b6 = a3**2 + 4 * a6
    b8 = a1**2 * a6 + 4 * a2 * a6 - a1 * a3 * a4 + a2 * a3**2 - a4**2
    c4 = b2**2 - 24 * b4
    c6 = -b2**3 + 36 * b2 * b4 - 216 * b6
    delta = -b2**2 * b8 - 8 * b4**3 - 27 * b6**2 + 9 * b2 * b4 * b6
    if not delta or c4**3 - c6**2 != 1728 * delta:
        raise ArithmeticError("target curve is singular or has inconsistent invariants")
    return c4, c6, delta


def target_record(snapshot, curve_id):
    rows = snapshot.get("snapshot", {}).get("curves")
    if not isinstance(rows, list):
        raise ValueError("target snapshot has no pinned curve list")
    matches = [row for row in rows if int(row.get("id", -1)) == curve_id]
    if len(matches) != 1:
        raise ValueError(f"target snapshot has {len(matches)} rows for curve {curve_id}")
    return matches[0]


def reduce_rational(value, field):
    value = Fraction(str(value))
    if value.denominator % field.characteristic() == 0:
        raise ZeroDivisionError("coefficient denominator vanishes modulo the selected prime")
    return field(value.numerator) / field(value.denominator)


def has_projective_root(poly, projective_degree, field):
    if not poly:
        return True, True, True
    ring = poly.parent()
    affine = poly.gcd(ring.gen() ** field.order() - ring.gen()).degree() > 0
    infinity = poly[projective_degree] == 0
    return bool(affine or infinity), bool(affine), bool(infinity)


def rational_roots(factorization):
    answer = []
    for factor, multiplicity in factorization:
        if factor.degree() == 1:
            answer.append({
                "parameter": str(-factor[0] / factor[1]),
                "multiplicity": int(multiplicity),
            })
    return answer


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--target-snapshot", type=Path, default=DEFAULT_SNAPSHOT)
    parser.add_argument("--curve-ids", default="302,273")
    parser.add_argument("--primes", default=",".join(map(str, DEFAULT_PRIMES)))
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    curve_ids = tuple(int(value) for value in args.curve_ids.split(",") if value.strip())
    primes = tuple(int(value) for value in args.primes.split(",") if value.strip())
    if not curve_ids or len(set(curve_ids)) != len(curve_ids):
        parser.error("--curve-ids must be a nonempty list without duplicates")
    if not primes or any(p <= 3 or not ZZ(p).is_prime() for p in primes):
        parser.error("--primes must contain primes greater than three")

    model_path = args.model.resolve()
    snapshot_path = args.target_snapshot.resolve()
    output_path = args.output.resolve()
    model = json.loads(model_path.read_text())
    if model.get("status") != "PASS_EXACT_QQ_FIXED_REVERSE_2A1_RR_JACOBIAN":
        raise ArithmeticError("input is not the exact fixed-corridor 2A1 Jacobian")
    if model.get("child", {}).get("fibre_configuration") not in (None, "2I2+20I1"):
        raise ArithmeticError("input child is not the declared 2A1 fibration")

    ring = PolynomialRing(QQ, "s")
    s = ring.gen()
    child = model["child"]
    child_a = ring(child["minimal_A_coefficients_low_to_high"])
    child_b = ring(child["minimal_B_coefficients_low_to_high"])
    if child_a.degree() != 8 or child_b.degree() != 12:
        raise ArithmeticError("fixed 2A1 model lost its K3 degree profile")
    projective_degree = 24

    snapshot = json.loads(snapshot_path.read_text())
    records = []
    for curve_id in curve_ids:
        pinned = target_record(snapshot, curve_id)
        c4_fraction, c6_fraction, delta_fraction = invariants(pinned["ainvs"])
        target_a = -27 * QQ(str(c4_fraction))
        target_b = -54 * QQ(str(c6_fraction))
        target_curve = EllipticCurve(QQ, [target_a, target_b])
        comparison = child_a**3 * target_b**2 - target_a**3 * child_b**2
        modular = []
        survives = True
        first_exclusion = None
        for prime in primes:
            field = GF(prime)
            try:
                if reduce_rational(delta_fraction, field) == 0:
                    raise ZeroDivisionError("target has bad reduction")
                a_mod = PolynomialRing(field, "s")(
                    [reduce_rational(value, field) for value in child_a]
                )
                b_mod = a_mod.parent()(
                    [reduce_rational(value, field) for value in child_b]
                )
                if a_mod.degree() != 8 or b_mod.degree() != 12:
                    raise ZeroDivisionError("model degree drops")
                ta = reduce_rational(target_a, field)
                tb = reduce_rational(target_b, field)
                if not ta or not tb:
                    raise ZeroDivisionError("target short coefficient vanishes")
                comparison_mod = a_mod**3 * tb**2 - ta**3 * b_mod**2
                keep, affine, infinity = has_projective_root(
                    comparison_mod, projective_degree, field
                )
                modular.append({
                    "prime": prime,
                    "usable": True,
                    "affine_root_exists": affine,
                    "infinity_root_exists": infinity,
                    "survives": keep,
                })
                if not keep:
                    survives = False
                    first_exclusion = prime
                    break
            except ZeroDivisionError as error:
                modular.append({"prime": prime, "usable": False, "reason": str(error)})

        exact = None
        if survives:
            factorization = comparison.factor()
            roots = rational_roots(factorization)
            infinity = comparison[projective_degree] == 0
            specializations = []
            for root_record in roots:
                parameter = QQ(root_record["parameter"])
                specialized = EllipticCurve(
                    QQ, [QQ(child_a(parameter)), QQ(child_b(parameter))]
                )
                specializations.append({
                    "parameter": str(parameter),
                    "isomorphic_to_target_over_Q": bool(
                        specialized.is_isomorphic(target_curve)
                    ),
                })
            if infinity:
                specialized = EllipticCurve(QQ, [child_a[8], child_b[12]])
                specializations.append({
                    "parameter": "infinity",
                    "isomorphic_to_target_over_Q": bool(
                        specialized.is_isomorphic(target_curve)
                    ),
                })
            exact = {
                "comparison_degree": int(comparison.degree()),
                "projective_degree": projective_degree,
                "factor_degrees_with_multiplicity": [
                    [int(factor.degree()), int(multiplicity)]
                    for factor, multiplicity in factorization
                ],
                "rational_affine_parameters": roots,
                "rational_parameter_at_infinity": bool(infinity),
                "specializations": specializations,
                "qq_isomorphic_count": sum(
                    item["isomorphic_to_target_over_Q"] for item in specializations
                ),
            }
        records.append({
            "curve_id": curve_id,
            "snapshot_rank_lower_bound": int(pinned["snapshot_rank_lower_bound"]),
            "modular_tests": modular,
            "first_exclusion_prime": first_exclusion,
            "survived_declared_prime_chain": survives,
            "exact": exact,
            "status": (
                "PASS_MODULAR_EXCLUDED_FIXED_2A1_FIBRATION"
                if not survives
                else (
                    "PASS_EXACT_QQ_ISOMORPHIC_FIXED_2A1_FIBRATION"
                    if exact["qq_isomorphic_count"]
                    else "PASS_EXACT_NO_QQ_ISOMORPHIC_FIXED_2A1_PARAMETER"
                )
            ),
        })
        print(
            f"ICARMMW15|curve={curve_id}|survives={int(survives)}"
            f"|exclude={first_exclusion}|qq_isomorphic="
            f"{0 if exact is None else exact['qq_isomorphic_count']}",
            flush=True,
        )

    payload = {
        "schema": "elkies-k3.icarm-fixed-mw15-screen.v1",
        "status": "PASS_EXACT_BOUNDED_FIXED_2A1_MW15_TARGET_SCREEN",
        "scope": {
            "fibration_count": 1,
            "root_type": "2A1",
            "generic_mw_rank": 15,
            "curve_ids": list(curve_ids),
            "prime_chain": list(primes),
            "complete_for_all_A2_or_2A1_fibrations": False,
        },
        "records": records,
        "proof_boundary": (
            "This is an exact projective modular screen, followed when necessary by "
            "exact QQ factorization and twist-sensitive isomorphism tests, for one "
            "source-identified fixed-corridor 2A1/MW15 fibration. It is not an atlas "
            "or exhaustion of all A2/2A1/MW15 fibrations."
        ),
        "inputs": {
            relative(path): digest(path)
            for path in (Path(__file__).resolve(), model_path, snapshot_path)
        },
        "software": {"sage_version": SAGE_VERSION},
        "reproducing_command": (
            "sage -python elkies-k3/scripts/screen_icarm_fixed_mw15_fibrations.sage "
            f"--curve-ids {','.join(map(str, curve_ids))} --output {relative(output_path)}"
        ),
    }
    serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not output_path.is_file() or output_path.read_text() != serialized:
            raise ArithmeticError("stored fixed MW15 screen differs from replay")
    else:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(serialized)
    print(f"ICARMMW15|status={payload['status']}|output={relative(output_path)}")


if __name__ == "__main__":
    main()
