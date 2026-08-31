#!/usr/bin/env sage
"""Certify the second q=4 Q80 child over the exact coefficient curve.

The equation is recorded compactly as the binary-quartic Jacobian over
``QQ(d,p,q,e)`` together with the already certified rational substitution
``(d,p,q,e)=(d(u),p(u),q(u),e(u))``.  We nevertheless perform the complete
substitution into ``QQ(u)[W]`` here: the discriminant factorization is checked
in that function field, while one rational value of ``u`` certifies all open
minimality and squarefreeness conditions.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import *
from sage.misc.persist import load


ROOT = Path(__file__).resolve().parents[2]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--parameter",
    type=Path,
    default=ROOT / "artifacts/generated-results/"
    "q80-cm24-slope-8-87-qq-PDQE-parameter.json",
)
parser.add_argument(
    "--output",
    type=Path,
    default=ROOT / "artifacts/generated-results/"
    "q80-unmarked-second-q4-qq.json",
)
parser.add_argument(
    "--sample-cache",
    type=Path,
    default=ROOT / "artifacts/generated-results/"
    "q80-unmarked-second-q4-repeated-root-samples-qq.json",
)
parser.add_argument("--witness-u", type=int, default=0)
args = parser.parse_args()

load(str(ROOT / "elkies-k3/scripts/derive_q80_second_q4_pencil.sage"))


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


payload = json.loads(args.parameter.read_text())
if payload.get("schema") != "q80-cm24-qq-PDQE-parameter-v1":
    raise ValueError("unexpected Q80 parameter schema")

parameter_ring = PolynomialRing(QQ, "u")
u = parameter_ring.gen()
parameter_field = parameter_ring.fraction_field()


def parameter_function(name):
    record = payload["original_functions"][name]
    numerator = parameter_ring(record["numerator"])
    denominator = parameter_ring(record["denominator"])
    if denominator == 0 or numerator.gcd(denominator) != 1:
        raise ArithmeticError(f"{name}(u) is not reduced")
    if [int(numerator.degree()), int(denominator.degree())] != record["degrees"]:
        raise ArithmeticError(f"{name}(u) has inconsistent recorded degrees")
    return parameter_field(numerator / denominator)


parameter_values = tuple(
    parameter_function(name) for name in ("d", "p", "q", "e")
)
specialization = parameters.hom(parameter_values, parameter_field)


def specialize_scalar(value):
    value = K(value)
    numerator = specialization(parameters(value.numerator()))
    denominator = specialization(parameters(value.denominator()))
    if denominator == 0:
        raise ZeroDivisionError("coefficient curve lies on a pencil denominator")
    return parameter_field(numerator / denominator)


function_W_ring = PolynomialRing(parameter_field, "W")
function_W = function_W_ring.gen()


def specialize_W_polynomial(polynomial):
    return function_W_ring(
        [specialize_scalar(coefficient) for coefficient in KW(polynomial).list()]
    )


quartic_coefficients = tuple(
    specialize_W_polynomial(second_curve[index]) for index in range(5)
)
q0, q1, q2, q3, q4 = quartic_coefficients
invariant_i = 12*q4*q0 - 3*q3*q1 + q2**2
invariant_j = (
    72*q4*q2*q0 + 9*q3*q2*q1 - 27*q4*q1**2
    - 27*q3**2*q0 - 2*q2**3
)
child_a = -27*invariant_i
child_b = -27*invariant_j
child_delta = 4*child_a**3 + 27*child_b**2
if child_a.degree() != 6 or child_b.degree() != 9:
    raise ArithmeticError("unexpected second-child Weierstrass degrees")

# General multivariate gcd/factor routines are disproportionately expensive
# for these large rational functions.  Recover the repeated root from exact
# rational specializations, reserve eight values as holdouts, and then check
# the resulting rational function by exact substitution below.
sample_ring = PolynomialRing(QQ, "W")


def specialize_at(value, point):
    value = parameter_field(value)
    numerator = parameter_ring(value.numerator())(point)
    denominator = parameter_ring(value.denominator())(point)
    if denominator == 0:
        raise ZeroDivisionError
    return QQ(numerator / denominator)


parameter_hash = sha256(args.parameter)
root_samples = []
if args.sample_cache.is_file():
    cache = json.loads(args.sample_cache.read_text())
    if (
        cache.get("schema") == "q80-unmarked-second-q4-root-samples-qq-v1"
        and cache.get("parameter_sha256") == parameter_hash
    ):
        root_samples = [
            (QQ(record["u"]), QQ(record["root"])) for record in cache["samples"]
        ]
if not root_samples:
    candidate_points = [QQ(0)]
    for integer in range(1, 100):
        candidate_points.extend((QQ(integer), QQ(-integer)))
    for point in candidate_points:
        try:
            specialized_a = sample_ring(
                [specialize_at(coefficient, point) for coefficient in child_a.list()]
            )
        except ZeroDivisionError:
            continue
        repeated = specialized_a.gcd(specialized_a.derivative()).monic()
        if repeated.degree() != 1:
            continue
        root_samples.append((point, -repeated[0]))
        if len(root_samples) == 64:
            break
    args.sample_cache.parent.mkdir(parents=True, exist_ok=True)
    args.sample_cache.write_text(json.dumps({
        "schema": "q80-unmarked-second-q4-root-samples-qq-v1",
        "parameter_sha256": parameter_hash,
        "samples": [
            {"u": str(point), "root": str(value)} for point, value in root_samples
        ],
    }, indent=2, sort_keys=True) + "\n")
if len(root_samples) != 64:
    raise ArithmeticError("could not collect 64 regular repeated-root samples")

reconstruction_samples = root_samples[:-8]
holdout_samples = root_samples[-8:]


def modular_degree_profile(prime):
    """Discover the Padé degree profile cheaply in one finite field."""
    finite = GF(prime)
    finite_ring = PolynomialRing(finite, "u")
    finite_u = finite_ring.gen()
    samples = [(finite(point), finite(value)) for point, value in root_samples]
    training = samples[:-8]
    interpolation = finite_ring.zero()
    modulus = finite_ring.one()
    for point, value in training:
        interpolation += ((value-interpolation(point))/modulus(point))*modulus
        modulus *= finite_u-point
    interpolation %= modulus
    r0, r1 = modulus, interpolation
    t0, t1 = finite_ring.zero(), finite_ring.one()
    candidates = []
    while r1:
        if t1:
            common = r1.gcd(t1)
            numerator = r1 // common
            denominator = t1 // common
            if all(
                denominator(point) != 0
                and numerator(point) == value*denominator(point)
                for point, value in samples
            ):
                candidates.append((
                    numerator.degree()+denominator.degree(),
                    numerator.degree(), denominator.degree(),
                ))
        quotient, r2 = r0.quo_rem(r1)
        t2 = t0-quotient*t1
        r0, r1 = r1, r2
        t0, t1 = t1, t2
    if not candidates:
        raise ArithmeticError(f"no modular Padé candidate at prime {prime}")
    return min(candidates)


profiles = tuple(
    modular_degree_profile(prime) for prime in (1000000007, 1000000009)
)
if profiles[0][1:] != profiles[1][1:]:
    raise ArithmeticError(f"modular Padé profiles disagree: {profiles}")
numerator_degree, denominator_degree = profiles[0][1:]
unknown_count = numerator_degree+1+denominator_degree
if unknown_count > len(reconstruction_samples):
    raise ArithmeticError("insufficient exact samples for recovered Padé profile")

# With a monic denominator the interpolation equations are linear and square.
# Solve once over QQ, then replay all unused reconstruction and holdout values.
rows = []
right_hand_side = []
for point, value in reconstruction_samples[:unknown_count]:
    rows.append(
        [point**degree for degree in range(numerator_degree+1)]
        + [-value*point**degree for degree in range(denominator_degree)]
    )
    right_hand_side.append(value*point**denominator_degree)
solution = Matrix(QQ, rows).solve_right(vector(QQ, right_hand_side))
solution_values = list(solution)
root_numerator = parameter_ring(solution_values[:numerator_degree+1])
root_denominator = parameter_ring(
    solution_values[numerator_degree+1:] + [QQ.one()]
)
if root_numerator.gcd(root_denominator).degree() != 0:
    raise ArithmeticError("exact repeated-root interpolant is not reduced")
if any(
    root_denominator(point) == 0
    or root_numerator(point) != value*root_denominator(point)
    for point, value in root_samples
):
    raise ArithmeticError("exact repeated-root interpolant misses a replay sample")
repeated_root = parameter_field(root_numerator/root_denominator)
linear_factor = function_W-repeated_root

# Exact function-field identities promote the interpolation from evidence to
# a certificate.  The later witness remains independent because it was not
# needed to choose among candidates: every one of the 64 samples was replayed.
if child_delta % linear_factor**7 or child_delta % linear_factor**8 == 0:
    raise ArithmeticError("generic repeated place does not have exponent seven")
generic_residual = child_delta // linear_factor**7
if generic_residual.degree() != 8:
    raise ArithmeticError("generic residual discriminant does not have degree eight")
factor_degrees = [(1, 7), (8, 1)]
finite_a_valuation = child_a.valuation(linear_factor)
finite_b_valuation = child_b.valuation(linear_factor)
if (finite_a_valuation, finite_b_valuation) != (2, 3):
    raise ArithmeticError("the exponent-seven place is not minimal I1*")

infinity_valuations = (
    int(8-child_a.degree()),
    int(12-child_b.degree()),
    int(24-child_delta.degree()),
)
if infinity_valuations != (2, 3, 9):
    raise ArithmeticError("the infinity place is not I3*")

witness = QQ(args.witness_u)
witness_map = parameter_ring.hom((witness,), QQ)


def at_witness(value):
    value = parameter_field(value)
    numerator = witness_map(parameter_ring(value.numerator()))
    denominator = witness_map(parameter_ring(value.denominator()))
    if denominator == 0:
        raise ZeroDivisionError("witness lies on an equation denominator")
    return QQ(numerator / denominator)


witness_ring = PolynomialRing(QQ, "W")
witness_a = witness_ring([at_witness(value) for value in child_a.list()])
witness_b = witness_ring([at_witness(value) for value in child_b.list()])
witness_delta = 4*witness_a**3 + 27*witness_b**2
witness_factors = tuple(witness_delta.factor())
if sorted((factor.degree(), int(exponent)) for factor, exponent in witness_factors) != [
    (1, 7), (8, 1)
]:
    raise ArithmeticError("witness discriminant has the wrong factorization")
witness_residual = next(factor for factor, exponent in witness_factors if exponent == 1)
if witness_residual.gcd(witness_residual.derivative()) != 1:
    raise ArithmeticError("witness residual discriminant is not squarefree")

# Store the substantially smaller ambient quartic formula.  Together with the
# hashed exact parameter artifact this is an explicit equation over QQ(u),
# while avoiding a second copy of its very large expanded coefficients.
ambient_coefficients = [str(second_curve[index]) for index in range(5)]
output = {
    "schema": "q80-unmarked-second-q4-qq-v1",
    "status": "PASS_EXACT_UNMARKED_SECOND_Q4",
    "inputs": {
        "parameter": {
            "path": str(args.parameter.relative_to(ROOT)),
            "sha256": sha256(args.parameter),
        },
        "pencil_derivation": {
            "path": "elkies-k3/scripts/derive_q80_second_q4_pencil.sage",
            "sha256": sha256(
                ROOT / "elkies-k3/scripts/derive_q80_second_q4_pencil.sage"
            ),
        },
    },
    "coordinates": {
        "first": "U=(x-T)/T^2",
        "local": "v=U-d+1",
        "second": "W=(X-3*v^3-x1*v-x0)/v^2",
        "x0": "-3*B1(d-1)/(2*A1(d-1))",
        "x1": "-A1'(d-1)/(6*x0)",
    },
    "binary_quartic": {
        "variable": "v",
        "base": "W",
        "coefficients_q0_through_q4_over_QQ_d_p_q_e": ambient_coefficients,
        "substitution": "(d,p,q,e)=(d(u),p(u),q(u),e(u))",
        "jacobian": {
            "I": "12*q4*q0-3*q3*q1+q2^2",
            "J": "72*q4*q2*q0+9*q3*q2*q1-27*q4*q1^2-27*q3^2*q0-2*q2^3",
            "A": "-27*I",
            "B": "-27*J",
        },
    },
    "checks": {
        "equation_ring": "QQ(u)[W]",
        "degrees_A_B_Delta": [
            int(child_a.degree()), int(child_b.degree()), int(child_delta.degree())
        ],
        "discriminant_factor_degrees_and_exponents": [
            list(value) for value in factor_degrees
        ],
        "finite_I1star_valuations_A_B_Delta": [
            int(finite_a_valuation), int(finite_b_valuation), 7
        ],
        "infinity_I3star_valuations_A_B_Delta": list(infinity_valuations),
        "witness_u": str(witness),
        "witness_residual_degree": int(witness_residual.degree()),
        "witness_residual_squarefree": True,
    },
    "generic_fibres": ["I3*", "I1*", "8 I1"],
    "generic_root_lattice": "D7+D5",
    "claim_boundary": {
        "proved": [
            "explicit second-q4 Jacobian over QQ(u)",
            "generic reducible fibre configuration I3*+I1*",
            "eight remaining generic nodal fibres",
        ],
        "not_proved": [
            "global marked Mordell-Weil sections",
            "the following q12 pencils over QQ(u)",
            "a rootless MW17 equation",
        ],
    },
}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")

print(
    "Q80UNMARKEDSECONDQ4|equation=QQ(u)[W]|"
    "Delta=(linear)^7*(squarefree_degree8)|finite=I1*|infinity=I3*|"
    "ADE=D7+D5|status=PASS_EXACT_UNMARKED_SECOND_Q4",
    flush=True,
)
