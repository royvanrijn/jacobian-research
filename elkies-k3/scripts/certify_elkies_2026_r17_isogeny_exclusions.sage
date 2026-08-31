#!/usr/bin/env sage
"""Exclude small-prime isogeny images of published-R17 fibres.

For each pinned target curve and ell in {3,5,7,11}, form the reduction of

    Phi_ell(j_R17(t), j(E)) = 0

after clearing the power of the denominator of ``j_R17``.  A rational
characteristic-zero parameter would reduce to a point of P^1(F_p) at every
clean prime.  The first clean prime with no finite root and no root at
infinity is therefore an exact exclusion witness.

This certificate concerns only the published R17 fibration.  It neither
constructs nor tests the equation-open alternate Q80-derived endpoint.
"""

import argparse
import hashlib
import json
import platform
from pathlib import Path

from sage.all import (
    GF,
    PolynomialRing,
    QQ,
    ZZ,
    gcd,
    lcm,
    prime_range,
)
from sage.schemes.elliptic_curves.mod_poly import classical_modular_polynomial
from sage.version import version as sage_version


ROOT = Path(__file__).resolve().parents[2]
MODEL = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_model.json"
TARGETS = ROOT / "elliptic-curves/data/elkies_2026_r17_j_recognition_targets.json"
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves/"
    "elkies_2026_published_r17_isogeny_exclusions_v1.json"
)
TARGET_ORDER = ("rank29", "curve398", "curve399", "curve400", "curve273", "curve302")
ELLS = (3, 5, 7, 11)


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--model", type=Path, default=MODEL)
parser.add_argument("--targets", type=Path, default=TARGETS)
parser.add_argument("--output", type=Path, default=OUTPUT)
parser.add_argument("--prime-bound", type=int, default=5000)
parser.add_argument("--check", action="store_true")
args = parser.parse_args()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def generalized_invariants(ainvs):
    a1, a2, a3, a4, a6 = map(QQ, ainvs)
    b2 = a1**2 + 4 * a2
    b4 = a1 * a3 + 2 * a4
    b6 = a3**2 + 4 * a6
    b8 = a1**2 * a6 + 4 * a2 * a6 - a1 * a3 * a4 + a2 * a3**2 - a4**2
    c4 = b2**2 - 24 * b4
    c6 = -b2**3 + 36 * b2 * b4 - 216 * b6
    delta = -b2**2 * b8 - 8 * b4**3 - 27 * b6**2 + 9 * b2 * b4 * b6
    if delta == 0:
        raise ArithmeticError("singular target curve")
    if c4**3 - c6**2 != 1728 * delta:
        raise ArithmeticError("target invariant identity failed")
    return {"c4": c4, "c6": c6, "discriminant": delta, "j": c4**3 / delta}


def primitive_integer_pair(numerator, denominator):
    values = list(numerator) + list(denominator)
    common_denominator = ZZ(1)
    for value in values:
        common_denominator = lcm(common_denominator, QQ(value).denominator())
    integers = [ZZ(QQ(value) * common_denominator) for value in values]
    content = gcd(integers)
    integers = [value // content for value in integers]
    numerator_integers = integers[: len(numerator)]
    denominator_integers = integers[len(numerator) :]
    if denominator_integers[-1] < 0:
        numerator_integers = [-value for value in numerator_integers]
        denominator_integers = [-value for value in denominator_integers]
    return numerator_integers, denominator_integers


def coefficient_hash(polynomial):
    raw = json.dumps(list(map(str, polynomial.list())), separators=(",", ":")).encode()
    return hashlib.sha256(raw).hexdigest()


def factor_degrees(polynomial):
    return sorted(
        [int(factor.degree()) for factor, multiplicity in polynomial.factor() for _ in range(multiplicity)]
    )


model_payload = json.loads(args.model.read_text())
target_payload = json.loads(args.targets.read_text())
if model_payload.get("schema") != "elkies-k3.elkies-2026-published-r17-model.v1":
    raise ValueError("unexpected published-R17 model schema")

rational_ring = PolynomialRing(QQ, "t")
t = rational_ring.gen()
A = rational_ring(model_payload["A_coefficients_low_to_high"])
B = rational_ring(model_payload["B_coefficients_low_to_high"])
if A.degree() > 8 or B.degree() > 12:
    raise ArithmeticError("published model is outside elliptic-K3 degree bounds")
c4 = -48 * A
c6 = -864 * B
surface_delta = -16 * (4 * A**3 + 27 * B**2)
if c4**3 - c6**2 != 1728 * surface_delta:
    raise ArithmeticError("surface invariant identity failed")
common = (c4**3).gcd(surface_delta)
j_numerator_QQ = (c4**3) // common
j_denominator_QQ = surface_delta // common
j_numerator_ZZ, j_denominator_ZZ = primitive_integer_pair(
    j_numerator_QQ.list(), j_denominator_QQ.list()
)
integer_ring = PolynomialRing(ZZ, "t")
j_numerator = integer_ring(j_numerator_ZZ)
j_denominator = integer_ring(j_denominator_ZZ)
map_degree = max(j_numerator.degree(), j_denominator.degree())
if map_degree != 24 or j_numerator.gcd(j_denominator) != 1:
    raise ArithmeticError("expected a reduced degree-24 published j-map")

target_rows = target_payload.get("targets")
if not isinstance(target_rows, list):
    raise ValueError("target input has no target list")
target_by_label = {row["label"]: row for row in target_rows}
missing = [label for label in TARGET_ORDER if label not in target_by_label]
if missing:
    raise ValueError(f"missing targets: {missing}")

modular_polynomials = {ell: classical_modular_polynomial(ell) for ell in ELLS}
modular_hashes = {
    str(ell): hashlib.sha256(str(modular_polynomials[ell]).encode()).hexdigest()
    for ell in ELLS
}

records = []
for label in TARGET_ORDER:
    target = target_by_label[label]
    invariants = generalized_invariants(target["ainvs"])
    target_j = invariants["j"]
    ell_records = []
    for ell in ELLS:
        phi = modular_polynomials[ell]
        X, Y = phi.parent().gens()
        x_degree = ell + 1
        if phi.degree(X) != x_degree or phi.degree(Y) != x_degree:
            raise ArithmeticError(f"unexpected Phi_{ell} degree")
        expected_degree = map_degree * x_degree
        witness = None
        clean_tested = 0
        for prime in prime_range(13, args.prime_bound):
            if prime == ell:
                continue
            finite = GF(prime)
            try:
                target_j_mod = finite(target_j)
                target_delta_mod = finite(invariants["discriminant"])
                numerator_mod = PolynomialRing(finite, "t")(j_numerator)
                denominator_mod = numerator_mod.parent()(j_denominator)
                surface_delta_mod = numerator_mod.parent()(surface_delta)
            except (ZeroDivisionError, TypeError, ValueError):
                continue
            if target_delta_mod == 0:
                continue
            if (
                numerator_mod.degree() != j_numerator.degree()
                or denominator_mod.degree() != j_denominator.degree()
                or numerator_mod.gcd(denominator_mod) != 1
                or surface_delta_mod.degree() != surface_delta.degree()
                or surface_delta_mod.gcd(surface_delta_mod.derivative()) != 1
            ):
                continue

            phi_mod = phi.change_ring(finite)
            finite_X, finite_Y = phi_mod.parent().gens()
            phi_x = phi_mod(X=finite_X, Y=target_j_mod)
            coefficients = [
                finite(phi_x.monomial_coefficient(finite_X**index))
                for index in range(x_degree + 1)
            ]
            recognition = numerator_mod.parent().zero()
            for index, coefficient in enumerate(coefficients):
                recognition += (
                    coefficient
                    * numerator_mod**index
                    * denominator_mod ** (x_degree - index)
                )
            if not recognition:
                continue
            # Homogenization has declared degree ``expected_degree``.  A drop
            # in affine degree is precisely a projective root at infinity, so
            # it cannot be used as a no-P^1-root witness.
            infinity_root = recognition.degree() < expected_degree
            if recognition.degree() > expected_degree:
                raise ArithmeticError("recognition degree exceeded the projective bound")
            clean_tested += 1
            degrees = factor_degrees(recognition)
            finite_linear_factors = degrees.count(1)
            if finite_linear_factors == 0 and not infinity_root:
                witness = {
                    "prime": int(prime),
                    "projective_degree": int(expected_degree),
                    "affine_degree": int(recognition.degree()),
                    "factor_degrees_with_multiplicity": degrees,
                    "finite_linear_factor_count": 0,
                    "infinity_root": False,
                    "squarefree": bool(recognition.gcd(recognition.derivative()) == 1),
                    "coefficient_sha256_low_to_high": coefficient_hash(recognition),
                    "clean_primes_tested_through_witness": clean_tested,
                }
                break
        if witness is None:
            raise ArithmeticError(
                f"no no-projective-root witness for {label}, ell={ell} below {args.prime_bound}"
            )
        ell_records.append({"ell": ell, "witness": witness})
        print(
            f"R17ISOGENY|target={label}|ell={ell}|p={witness['prime']}|"
            f"degree={witness['affine_degree']}|linear=0|infinity=0|status=EXCLUDED",
            flush=True,
        )
    records.append(
        {
            "label": label,
            "ainvs": list(map(str, target["ainvs"])),
            "target_j": str(target_j),
            "target_discriminant": str(invariants["discriminant"]),
            "isogeny_degrees": ell_records,
            "all_excluded": True,
        }
    )

output = {
    "schema": "elkies-k3.elkies-2026-published-r17-isogeny-exclusions.v1",
    "status": "PASS_PUBLISHED_R17_SMALL_PRIME_ISOGENY_EXCLUSIONS",
    "claim": (
        "For every listed target and ell in {3,5,7,11}, the displayed clean "
        "prime gives no point of P^1(F_p) satisfying the reduced modular-"
        "polynomial relation. Hence no rational characteristic-zero published-"
        "R17 fibre parameter can satisfy Phi_ell(j_R17(t),j(E))=0."
    ),
    "claim_boundary": (
        "This excludes only cyclic isogenies of degrees 3,5,7,11 from rational "
        "published-R17 fibres. It does not test composite degrees, algebraic "
        "fibre parameters, or the equation-open alternate Q80-derived fibration."
    ),
    "published_j_map": {
        "coordinate": "t",
        "degree": int(map_degree),
        "numerator_degree": int(j_numerator.degree()),
        "denominator_degree": int(j_denominator.degree()),
        "numerator_sha256_low_to_high": coefficient_hash(j_numerator),
        "denominator_sha256_low_to_high": coefficient_hash(j_denominator),
    },
    "clean_prime_gate": [
        "p>11 and p!=ell",
        "target discriminant and target-j denominator nonzero modulo p",
        "published j numerator and denominator retain degree and are coprime",
        "published surface discriminant retains degree and is squarefree",
        "recognition polynomial is nonzero and has projective degree at most 24*(ell+1)",
    ],
    "method": (
        "Specialize Phi_ell(X,Y) at Y=j(E) modulo p, substitute X=N(t)/D(t), "
        "and clear D(t)^(ell+1). Factor the resulting polynomial and separately "
        "test the projective point at infinity."
    ),
    "modular_polynomial_string_sha256": modular_hashes,
    "targets": records,
    "inputs": {
        "model": {"path": str(args.model.relative_to(ROOT)), "sha256": sha256(args.model)},
        "targets": {"path": str(args.targets.relative_to(ROOT)), "sha256": sha256(args.targets)},
        "worker": {
            "path": str(Path(__file__).resolve().relative_to(ROOT)),
            "sha256": sha256(Path(__file__).resolve()),
        },
    },
    "software": {
        "sage": str(sage_version),
        "python": platform.python_version(),
    },
    "prime_search_bound_exclusive": args.prime_bound,
}

args.output.parent.mkdir(parents=True, exist_ok=True)
serialized = json.dumps(output, indent=2, sort_keys=True) + "\n"
if args.check:
    if not args.output.exists() or args.output.read_text() != serialized:
        raise SystemExit(f"stale or missing artifact: {args.output}")
else:
    args.output.write_text(serialized)
print(
    f"R17ISOGENY|targets={len(records)}|ells={len(ELLS)}|"
    f"output={args.output}|status=PASS",
    flush=True,
)
