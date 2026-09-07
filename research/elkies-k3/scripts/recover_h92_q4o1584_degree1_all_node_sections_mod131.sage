#!/usr/bin/env sage -python
"""Recover q4/orbit1584 degree-one all-node sections modulo 131.

Write a section with P.O=1 as x=X/Z^2, y=Y/Z^3, where Z=t-z,
deg(X)<=5, and deg(Y)<=7.  Requiring the section to meet the singular
point of every finite reducible fibre fixes the four values of X and makes
Y divisible by their support polynomial L.  For each z the X interpolation
space has only two free parameters; exhaust those parameters and test

    (X^3+A*X*Z^4+B*Z^6)/L^2

for a polynomial square.  This is a finite-field seed search only.  It uses
neither elimination nor a Groebner basis.
"""

import hashlib
import json
import time
from itertools import product
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
MODEL = LOCAL / "q4o1584-compact-weierstrass-qq.json"
OUTPUT = LOCAL / "q4o1584-degree1-all-node-sections-mod131.json"
PRIME = 131

started = time.monotonic()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def log(stage, **fields):
    tail = "|".join(f"{key}={value}" for key, value in fields.items())
    print(
        f"Q4O1584D1NODES|stage={stage}|elapsed={time.monotonic()-started:.3f}"
        + (f"|{tail}" if tail else ""),
        flush=True,
    )


compact = json.loads(MODEL.read_text())
assert compact["status"] == "PASS_EXACT_QQ_Q4O1584_COMPACT_WEIERSTRASS_NORMALIZATION"

k = GF(PRIME)
R = PolynomialRing(k, "t")
t = R.gen()


def reduce_qq(value):
    value = QQ(value)
    return k(value.numerator()) / k(value.denominator())


A = R([reduce_qq(value) for value in compact["compact_model"]["A_coefficients_low_to_high"]])
B = R([reduce_qq(value) for value in compact["compact_model"]["B_coefficients_low_to_high"]])
records = compact["compact_model"]["finite_reducible_fibres"]
assert [record["kodaira"] for record in records] == ["I2", "I2", "I2", "I4"]
supports = [reduce_qq(record["support"]) for record in records]
assert len(set(supports)) == 4

RX = PolynomialRing(k, "x")
x = RX.gen()
nodes = []
for support in supports:
    cubic = x**3 + A(support) * x + B(support)
    repeated = cubic.gcd(cubic.derivative())
    assert repeated.degree() == 1
    node = -repeated[0] / repeated[1]
    assert A(support) == -3 * node**2 and B(support) == 2 * node**3
    nodes.append(node)

L = R.prod(t - support for support in supports)
L2 = L**2
evaluation = matrix(k, [[support**degree for degree in range(6)] for support in supports])
kernel = evaluation.right_kernel().basis()
assert len(kernel) == 2
# A fixed right inverse turns the z-dependent four target values into one
# degree-five interpolation polynomial without solving a system per z.
right_inverse = evaluation.solve_right(matrix.identity(k, 4))
assert evaluation * right_inverse == matrix.identity(k, 4)
log("LOAD", supports=",".join(map(str, supports)), nodes=",".join(map(str, nodes)))


def square_roots(polynomial):
    """Return both polynomial square roots, using leading coefficients."""
    polynomial = R(polynomial)
    if not polynomial:
        return (R.zero(),)
    degree = polynomial.degree()
    if degree % 2:
        return ()
    root_degree = degree // 2
    leading = polynomial[degree]
    if not leading.is_square():
        return ()
    answers = []
    for top in leading.sqrt(all=True):
        coefficients = [k.zero()] * (root_degree + 1)
        coefficients[root_degree] = top
        for target_degree in range(2 * root_degree - 1, root_degree - 1, -1):
            unknown_degree = target_degree - root_degree
            known = sum(
                (
                    coefficients[left] * coefficients[target_degree-left]
                    for left in range(unknown_degree + 1, root_degree)
                    if 0 <= target_degree - left <= root_degree
                ),
                k.zero(),
            )
            coefficients[unknown_degree] = (polynomial[target_degree] - known) / (2 * top)
        candidate = R(coefficients)
        if candidate**2 == polynomial:
            answers.append(candidate)
    return tuple(answers)


answers = []
tests = 0
divisibility_failures = 0
for z in k:
    # If Z vanishes at a reducible support, these projective coordinates meet
    # the section at infinity there rather than the prescribed finite node.
    if z in supports:
        log("SKIP_POLE_AT_NODE", z=int(z))
        continue
    Z = t - z
    targets = vector(k, [node * (support - z)**2 for support, node in zip(supports, nodes)])
    particular = right_inverse * targets
    for free in product(k, repeat=2):
        coefficients = particular + free[0] * kernel[0] + free[1] * kernel[1]
        X = R(list(coefficients))
        rhs = X**3 + A * X * Z**4 + B * Z**6
        quotient, remainder = rhs.quo_rem(L2)
        tests += 1
        if remainder:
            divisibility_failures += 1
            continue
        for Q in square_roots(quotient):
            Y = L * Q
            assert Y**2 == rhs
            assert X.degree() <= 5 and Y.degree() <= 7
            assert all(
                X(support) == node * Z(support)**2 and Y(support) == 0
                for support, node in zip(supports, nodes)
            )
            answers.append({
                "z": int(z),
                "free_parameters": [int(value) for value in free],
                "Z_coefficients_low_to_high": [int(value) for value in Z.list()],
                "X_coefficients_low_to_high": [int(value) for value in X.list()],
                "Y_coefficients_low_to_high": [int(value) for value in Y.list()],
                "Q_coefficients_low_to_high": [int(value) for value in Q.list()],
            })
    log("Z", z=int(z), tests=tests, sections=len(answers))

# The two signs of Y are distinct sections.  Guard against an accidental
# duplicate caused by a degenerate zero square root.
unique = {
    (tuple(record["Z_coefficients_low_to_high"]),
     tuple(record["X_coefficients_low_to_high"]),
     tuple(record["Y_coefficients_low_to_high"])): record
    for record in answers
}
answers = [unique[key] for key in sorted(unique)]

payload = {
    "schema": "elkies-k3.q4o1584-degree1-all-node-sections-mod131.v1",
    "status": "PASS_MOD131_Q4O1584_DEGREE1_ALL_NODE_SECTION_SEARCH",
    "prime": PRIME,
    "ansatz": {
        "coordinates": "x=X/Z^2, y=Y/Z^3",
        "Z": "t-z",
        "maximum_degrees_X_Y": [5, 7],
        "finite_node_supports": [int(value) for value in supports],
        "finite_node_x_coordinates": [int(value) for value in nodes],
        "support_polynomial_L_coefficients_low_to_high": [int(value) for value in L.list()],
    },
    "search": {
        "parameters": ["z", "two coordinates in the X-interpolation kernel"],
        "tests": tests,
        "expected_tests": (PRIME - len(supports)) * PRIME**2,
        "known_L_squared_divisibility_failures": divisibility_failures,
        "sections_found_including_sign": len(answers),
        "sections": answers,
        "excluded_denominator_poles_at_reducible_supports": [int(value) for value in supports],
        "exhaustive_for_displayed_ansatz": tests == (PRIME - len(supports)) * PRIME**2,
    },
    "method": {
        "large_Groebner_required": False,
        "polynomial_factorization_required": False,
        "runtime_seconds": time.monotonic() - started,
    },
    "proof_boundary": (
        "This exhaustive finite-field search supplies modular branches for exact lifting; "
        "it does not by itself prove a characteristic-zero section. Every returned branch "
        "is checked by literal substitution modulo 131 and meets all four prescribed nodes."
    ),
    "next_required": (
        "Lift a regular resolved branch to QQ, identify its exact lattice class, and transport "
        "it through the q4/orbit164 quartic map."
    ),
    "inputs": {
        "paths": [str(MODEL.relative_to(ROOT))],
        "sha256": {str(MODEL.relative_to(ROOT)): sha256(MODEL)},
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "Q4O1584D1NODES|tests={}|sections={}|divisibility_failures={}|status={}|output={}".format(
        tests, len(answers), divisibility_failures, payload["status"], OUTPUT,
    ),
    flush=True,
)
