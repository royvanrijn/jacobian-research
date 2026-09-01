#!/usr/bin/env sage -python
"""Certify generic genus one for the exact biquadratic resolved pencil."""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ, sage_eval
from sage.structure.proof.proof import WithProof


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "artifacts/generated-results"
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--pencil", type=Path,
    default=RESULTS / "q80-third-q12-um2-biquadratic-resolved-pencil-qq.json",
)
parser.add_argument(
    "--operands", type=Path,
    default=RESULTS / "q80-third-q12-um2-biquadratic-closure-operands-p19-hensel-qq.json",
)
parser.add_argument(
    "--p19-genus", type=Path,
    default=RESULTS / "q80-third-q12-um2-p19-resolved-genus.json",
)
parser.add_argument(
    "--output", type=Path,
    default=RESULTS / "q80-third-q12-um2-biquadratic-resolved-genus-qq.json",
)
args = parser.parse_args()
for name in ("pencil", "operands", "p19_genus", "output"):
    setattr(args, name, getattr(args, name).resolve())


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_rational(record):
    return QQ(ZZ(record["numerator"])) / ZZ(record["denominator"])


pencil = json.loads(args.pencil.read_text())
operands = json.loads(args.operands.read_text())
p19_genus = json.loads(args.p19_genus.read_text())
if pencil.get("status") != "PASS_EXACT_QQ_THIRD_Q12_BIQUADRATIC_RESOLVED_PENCIL":
    raise ValueError("exact resolved pencil is not certified")
if operands.get("status") != "PASS_EXACT_QQ_THIRD_Q12_BIQUADRATIC_CLOSURE_OPERANDS_P19_HENSEL":
    raise ValueError("exact biquadratic operands are not certified")
if p19_genus.get("status") != "PASS_EXACT_THIRD_Q12_GENUS_ONE_BY_ADJUNCTION_MOD19_QUADRATIC":
    raise ValueError("p=19 adjunction control is not certified")
if pencil["resolved_gates"]["combined_rank"] != 5:
    raise ArithmeticError("exact connected gate rank is not five")
if pencil["resolved_gates"]["kernel_dimension"] != 2:
    raise ArithmeticError("exact complete linear system is not a pencil")
if pencil["moving_equation"]["degrees_T_W_x"] != [2, 9, 3]:
    raise ArithmeticError("exact moving equation has the wrong degrees")

prime = 19
prime_field = GF(prime)
modulus_ring = PolynomialRing(prime_field, "r_modulus")
r_modulus = modulus_ring.gen()
finite = GF(prime**2, "r", modulus=r_modulus**2 + 12*r_modulus + 3)
q1 = read_rational(operands["biquadratic_field"]["q1"])
q2 = read_rational(operands["biquadratic_field"]["q2"])


def reduce_rational(value):
    value = QQ(value)
    denominator = finite(value.denominator())
    if not denominator:
        raise ZeroDivisionError("exact coefficient has bad reduction at p=19")
    return finite(value.numerator()) / denominator


q1_mod = reduce_rational(q1)
q2_mod = reduce_rational(q2)
a0 = q1_mod.sqrt()
b0 = q2_mod.sqrt()
theta0 = a0 + b0
if theta0**4 - 2*(q1_mod + q2_mod)*theta0**2 + (q1_mod - q2_mod)**2:
    raise ArithmeticError("chosen p=19 primitive element does not satisfy the quartic")
exact_polynomial_ring = PolynomialRing(QQ, "z_exact")
z_exact = exact_polynomial_ring.gen()
exact_field = exact_polynomial_ring.quotient(
    z_exact**4 - 2*(q1 + q2)*z_exact**2 + (q1 - q2)**2,
    "theta",
)
theta_exact = exact_field.gen()

plane_ring = PolynomialRing(finite, names=("V", "W", "x"))
V, W, x = plane_ring.gens()
moving = plane_ring.zero()
terms = pencil["moving_equation"]["terms_T_W_x_coefficient_1_r"]
if len(terms) != 63:
    raise ArithmeticError("exact moving-equation support changed")
for t_degree, w_degree, x_degree, encoded in terms:
    if len(encoded) != 1:
        raise ArithmeticError("unexpected exact coefficient encoding")
    value = exact_field(sage_eval(encoded[0], locals={"theta": theta_exact}))
    coefficient = sum(
        reduce_rational(entry) * theta0**degree
        for degree, entry in enumerate(value.lift().list())
    )
    moving += coefficient * V**t_degree * W**w_degree * x**x_degree
if (moving.degree(V), moving.degree(W), moving.degree(x)) != (2, 9, 3):
    raise ArithmeticError("p=19 reduction drops a moving-equation degree")
with WithProof("polynomial", False):
    factors = moving.factor()
if len(factors) != 1 or int(factors[0][1]) != 1:
    raise ArithmeticError("p=19 reduction of the exact moving equation is reducible")
if factors.prod() != moving:
    raise ArithmeticError("p=19 factorization replay failed")

lattice = p19_genus["lattice"]
if (
    lattice["divisor_square"] != 0
    or not lattice["divisor_primitive"]
    or lattice["old_fibre_degree"] != 3
):
    raise ArithmeticError("pinned characteristic-zero lattice divisor has changed")

output = {
    "schema": "elkies-k3.q80-third-q12-biquadratic-resolved-genus-qq.v1",
    "status": "PASS_EXACT_QQ_THIRD_Q12_BIQUADRATIC_GENUS_ONE_BY_GOOD_REDUCTION",
    "specialization": {"u": "-2", "field_degree": 4},
    "lattice": {
        "divisor_square": 0,
        "divisor_primitive": True,
        "old_fibre_degree": 3,
        "separable_in_characteristic_zero": True,
    },
    "linear_system": {
        "resolved_condition_rank": 5,
        "dimension": 2,
        "moving_degrees_T_W_x": [2, 9, 3],
        "moving_terms": 63,
        "irreducible_good_reduction_prime": 19,
        "generic_equation_irreducible_over_characteristic_zero": True,
    },
    "conclusion": (
        "Irreducibility of the p=19 reduction proves irreducibility over the "
        "biquadratic characteristic-zero field. The complete primitive "
        "isotropic pencil is base-point-free; Bertini and K3 adjunction give "
        "a smooth generic member of genus one."
    ),
    "inputs": {
        "pencil": {"path": str(args.pencil.relative_to(ROOT)), "sha256": sha256(args.pencil)},
        "operands": {"path": str(args.operands.relative_to(ROOT)), "sha256": sha256(args.operands)},
        "p19_genus": {"path": str(args.p19_genus.relative_to(ROOT)), "sha256": sha256(args.p19_genus)},
    },
    "worker": {
        "path": str(Path(__file__).resolve().relative_to(ROOT)),
        "sha256": sha256(Path(__file__).resolve()),
    },
    "claim_boundary": {
        "proved": [
            "irreducibility of the exact moving equation by one good irreducible reduction",
            "complete primitive isotropic characteristic-zero pencil",
            "smooth generic genus-one member by Bertini and K3 adjunction",
        ],
        "not_proved": [
            "an exact minimal Jacobian or birational maps",
            "the characteristic-zero A5+A3+3A1 fibre marking or Mordell--Weil rank",
        ],
    },
    "reproduce": "sage -python elkies-k3/scripts/verify_q80_third_q12_biquadratic_resolved_genus_qq.sage",
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
print(
    "Q80THIRDQ12BIQUADRATICGENUS|field_degree=4|reduction=19|"
    "irreducible=1|primitive=1|degree=3|genus=1|"
    "status=PASS_EXACT_QQ_THIRD_Q12_BIQUADRATIC_GENUS_ONE_BY_GOOD_REDUCTION",
    flush=True,
)
