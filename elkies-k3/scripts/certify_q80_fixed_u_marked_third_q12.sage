#!/usr/bin/env sage
"""Search fixed rational Q80 fibres for the marked third-q12 horizontal.

This worker deliberately specializes before attempting a ``QQ(u)`` lift.  It
has three fail-closed layers:

1. specialize the exact rational Q80 coefficient curve and construct the
   source, first-q4, and second-q4 equations over ``QQ``;
2. test the forced first Q80 marking exactly and, independently, enumerate
   polynomial sections of the ``D7+D5`` child at several good primes;
3. combine the modular polynomial shell with decoded direct-section charts,
   and independently export denominator-one or denominator-two schemes;
   retain shell sections only when they have
   the equation-side fingerprint of the third-q12 horizontal:

       P.O=2, height 8, identity at both I*_n fibres.

The script does not call the CM24 third-q12 seed and does not accept a CM24
horizontal.  A modular candidate is reconnaissance, not an exact section.
The terminal success status is therefore intentionally reserved for a later
exact reconstruction/connected-quotient/child compilation stage.  In
particular, this worker never reports ``A5+A3+3A1/MW6`` merely from its lattice
label.

The declared parameter list starts with the projective point ``u=infinity``
and is then ordered by affine projective height.  The first affine entry
``u=0`` is unusually small on the reconstructed coefficient curve; subsequent
entries are the smallest reduced rationals, with negative values first at a
fixed height.  Use ``--parameter-limit`` to enlarge the finite experiment,
and ``--prime`` repeatedly to change the modular ensemble.
"""

import argparse
import hashlib
import itertools
import json
import shutil
import subprocess
import time
from pathlib import Path

from sage.all import (
    EllipticCurve,
    GF,
    Integer,
    PolynomialRing,
    QuadraticForm,
    QQ,
    ZZ,
    lcm,
    matrix,
    prod,
    vector,
)
from sage.misc.persist import load


ROOT = Path(__file__).resolve().parents[2]
PARAMETER = (
    ROOT
    / "artifacts/generated-results/q80-cm24-slope-8-87-qq-PDQE-parameter.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/q80-fixed-u-marked-third-q12-search.json"
)
GLOBAL_FIRST_MARKING_COVER = (
    ROOT
    / "artifacts/generated-results/q80-slope-8-87-first-marked-cover-qq.json"
)
SIMPLIFIED_FIRST_MARKING_COVER = (
    ROOT
    / "artifacts/generated-results/q80-first-marked-cover-simplified-qq.json"
)
RANK_FIVE_LATTICE_TARGET = (
    ROOT / "artifacts/generated-results/q80-d7d5-mw5-height-lattice.json"
)

# The projective point at infinity and the first forty affine reduced
# rationals in the declared deterministic height ordering.  This is data, not
# a dynamically changing search boundary.
PREDECLARED_U = (None,) + tuple(
    map(
        QQ,
        (
            "0", "-1", "1", "-2", "2", "-1/2", "1/2", "-3", "3",
            "-3/2", "3/2", "-1/3", "1/3", "-2/3", "2/3", "-4", "4",
            "-4/3", "4/3", "-1/4", "1/4", "-3/4", "3/4", "-5", "5",
            "-5/2", "5/2", "-5/3", "5/3", "-5/4", "5/4", "-1/5",
            "1/5", "-2/5", "2/5", "-3/5", "3/5", "-4/5", "4/5", "-6",
        ),
    )
)

TARGET_ROOTS = "A5+A3+3A1"
TARGET_MW_RANK = 6
TARGET_PO = 2
TARGET_HEIGHT = QQ(8)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def coefficient_strings(poly):
    return [str(value) for value in poly.list()]


def finite_coefficients(poly):
    return [int(value) for value in poly.list()]


def rational_record(value):
    value = QQ(value)
    return {
        "value": str(value),
        "numerator": str(value.numerator()),
        "denominator": str(value.denominator()),
    }


def parameter_value(value):
    if value is None:
        return None
    if str(value).lower() in ("infinity", "inf", "oo"):
        return None
    return QQ(value)


def parameter_label(value):
    return "infinity" if value is None else str(QQ(value))


def valuation(poly, factor):
    if not poly:
        return None
    answer = 0
    while poly % factor == 0:
        poly //= factor
        answer += 1
    return answer


def square_root_polynomial(poly):
    """Return the monic-free exact square root, or ``None``."""
    if not poly:
        return poly.parent().zero()
    unit = poly.factor().unit()
    if not unit.is_square():
        return None
    answer = poly.parent()(unit.sqrt())
    for factor, exponent in poly.factor():
        if int(exponent) % 2:
            return None
        answer *= factor ** (int(exponent) // 2)
    if answer**2 != poly:
        return None
    return answer


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--parameter", type=Path, default=PARAMETER)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument(
    "--parameter-limit",
    type=int,
    default=1,
    help="number of entries from the fixed PREDECLARED_U list to process",
)
parser.add_argument(
    "--u",
    action="append",
    default=[],
    help="explicit rational parameter; repeat to override the predeclared prefix",
)
parser.add_argument(
    "--prime",
    action="append",
    type=int,
    default=[],
    help="good-prime candidate; repeat (default: 23)",
)
parser.add_argument(
    "--word-length",
    type=int,
    default=3,
    help="maximum signed polynomial-shell word length",
)
parser.add_argument(
    "--max-polynomial-pairs",
    type=int,
    default=12,
    help="fail closed if the unsigned polynomial shell is larger than this",
)
parser.add_argument(
    "--max-height-shell",
    type=int,
    default=10000,
    help="fail closed if the unsigned height-eight lattice shell exceeds this cap",
)
parser.add_argument(
    "--marking-only",
    action="store_true",
    help="run exact equations and the forced first-marking gate, but no modular shell",
)
parser.add_argument(
    "--prime-audit-only",
    action="store_true",
    help="check reduction signatures without enumerating polynomial sections",
)
parser.add_argument(
    "--height-shell",
    action="store_true",
    help=(
        "recover the exact modular Shioda lattice from fourth multiples and "
        "enumerate its complete height-eight shell"
    ),
)
parser.add_argument(
    "--po1-slice-certificate",
    type=Path,
    default=None,
    help=(
        "decoded output of scan_q80_po1_msolve_slices.sage; verified modular "
        "sections are adjoined before height-lattice recovery"
    ),
)
parser.add_argument(
    "--direct-msolve-dir",
    type=Path,
    default=None,
    help=(
        "export the direct rational-section schemes here; "
        "one file is written for each good (u,prime,relative-sign)"
    ),
)
parser.add_argument(
    "--direct-chart",
    choices=("auxiliary", "recursive"),
    default="auxiliary",
    help=(
        "direct section scheme: sparse auxiliary-M coefficients (default) or "
        "dense top-down recursive elimination"
    ),
)
parser.add_argument(
    "--direct-pole-order",
    type=int,
    choices=(0, 1, 2),
    default=2,
    help="section pole order for the exported direct scheme (default: 2)",
)
parser.add_argument(
    "--direct-pole-location",
    choices=("finite", "infinity"),
    default="finite",
    help="location of the unique pole when --direct-pole-order=1",
)
parser.add_argument(
    "--run-msolve",
    action="store_true",
    help="run msolve on every exported direct section scheme",
)
parser.add_argument("--msolve-threads", type=int, default=4)
parser.add_argument("--msolve-timeout", type=int, default=300)
args = parser.parse_args()

if args.parameter_limit < 1 or args.parameter_limit > len(PREDECLARED_U):
    raise ValueError("parameter-limit lies outside the predeclared finite list")
if args.word_length < 1 or args.word_length > 20:
    raise ValueError("word-length must lie between one and twenty")
if args.run_msolve and args.direct_msolve_dir is None:
    raise ValueError("--run-msolve requires --direct-msolve-dir")
if args.msolve_threads < 1 or args.msolve_timeout < 1:
    raise ValueError("msolve threads and timeout must be positive")
if args.max_height_shell < 1:
    raise ValueError("max-height-shell must be positive")
if args.direct_pole_order != 1 and args.direct_pole_location != "finite":
    raise ValueError("--direct-pole-location applies only to P.O=1")
if args.run_msolve and shutil.which("msolve") is None:
    raise RuntimeError("--run-msolve requested but msolve is unavailable")

selected_u = (
    tuple(parameter_value(value) for value in args.u)
    if args.u
    else PREDECLARED_U[: args.parameter_limit]
)
primes = tuple(ZZ(value) for value in (args.prime or (23,)))
if any(not prime.is_prime() or prime in (2, 3) for prime in primes):
    raise ValueError("modular primes must be prime and different from 2,3")

started = time.monotonic()
parameter_payload = json.loads(args.parameter.read_text(encoding="utf-8"))
if parameter_payload.get("schema") != "q80-cm24-qq-PDQE-parameter-v1":
    raise ValueError("unexpected exact Q80 coefficient-curve schema")

# Load the universal q4/q4 construction once.  This deliberately uses no
# marked CM24 third-q12 file.
load(str(ROOT / "elkies-k3/scripts/derive_q80_second_q4_pencil.sage"))

parameter_u_ring = PolynomialRing(QQ, "u_parameter")
u_parameter = parameter_u_ring.gen()
parameter_u_field = parameter_u_ring.fraction_field()


def read_parameter_function(name):
    record = parameter_payload["original_functions"][name]
    numerator = parameter_u_ring(record["numerator"].replace("u", "u_parameter"))
    denominator = parameter_u_ring(record["denominator"].replace("u", "u_parameter"))
    if not denominator or numerator.gcd(denominator) != 1:
        raise ArithmeticError(f"stored {name}(u) is not reduced")
    if [int(numerator.degree()), int(denominator.degree())] != record["degrees"]:
        raise ArithmeticError(f"stored {name}(u) has inconsistent degrees")
    return parameter_u_field(numerator / denominator)


parameter_functions = {
    name: read_parameter_function(name) for name in ("d", "p", "q", "e")
}


def coefficient_values(point):
    values = {}
    for name, function in parameter_functions.items():
        if point is None:
            numerator = function.numerator().leading_coefficient()
            denominator = function.denominator().leading_coefficient()
        else:
            numerator = function.numerator()(point)
            denominator = function.denominator()(point)
        if not denominator:
            raise ZeroDivisionError(
                f"u={parameter_label(point)} is a pole of {name}(u)"
            )
        values[name] = QQ(numerator / denominator)
    return values


def specialize_scalar(value, homomorphism):
    value = K(value)
    denominator = homomorphism(parameters(value.denominator()))
    if not denominator:
        raise ZeroDivisionError("specialization hits a universal-q4 denominator")
    return QQ(homomorphism(parameters(value.numerator())) / denominator)


def specialize_polynomial(poly, source_ring, target_ring, homomorphism):
    return target_ring(
        [specialize_scalar(value, homomorphism) for value in source_ring(poly).list()]
    )


def source_equation(values):
    map_to_qq = parameters.hom(
        tuple(values[name] for name in ("d", "p", "q", "e")), QQ
    )
    source_ring = PolynomialRing(QQ, "T")
    source_A = specialize_polynomial(A, KT, source_ring, map_to_qq)
    source_B = specialize_polynomial(B, KT, source_ring, map_to_qq)
    return source_ring, source_A, source_B, map_to_qq


def first_child_equation(map_to_qq):
    ring = PolynomialRing(QQ, "U")
    first_A = specialize_polynomial(jacobian_a, KU, ring, map_to_qq)
    first_B = specialize_polynomial(jacobian_b, KU, ring, map_to_qq)
    delta = 4 * first_A**3 + 27 * first_B**2
    factors = sorted((factor.degree(), int(exponent)) for factor, exponent in delta.factor())
    if factors != [(1, 5), (8, 1)]:
        # On the rank-19 coefficient curve the ambient I4 collides once with
        # the residual nonic, giving finite I5.  Any further collision is bad.
        raise ArithmeticError(("bad first-q4 specialization", factors))
    infinity = (8 - first_A.degree(), 12 - first_B.degree(), 24 - delta.degree())
    if infinity != (2, 3, 11):
        raise ArithmeticError(("bad first-q4 infinity", infinity))
    return ring, first_A, first_B, delta


def second_child_equation(map_to_qq):
    ring = PolynomialRing(QQ, "W")
    local_ring = PolynomialRing(ring, "v")
    curve = local_ring(
        [
            ring(
                [specialize_scalar(value, map_to_qq) for value in KW(coefficient).list()]
            )
            for coefficient in second_curve.list()
        ]
    )
    q0, q1, q2, q3, q4 = [curve[index] for index in range(5)]
    invariant_i = 12 * q4 * q0 - 3 * q3 * q1 + q2**2
    invariant_j = (
        72 * q4 * q2 * q0
        + 9 * q3 * q2 * q1
        - 27 * q4 * q1**2
        - 27 * q3**2 * q0
        - 2 * q2**3
    )
    second_A = -27 * invariant_i
    second_B = -27 * invariant_j
    delta = 4 * second_A**3 + 27 * second_B**2
    factors = tuple(delta.factor())
    signature = sorted((factor.degree(), int(exponent)) for factor, exponent in factors)
    if signature != [(1, 7), (8, 1)]:
        raise ArithmeticError(("bad second-q4 specialization", signature))
    star_factor = next(factor.monic() for factor, exponent in factors if int(exponent) == 7)
    if (
        valuation(second_A, star_factor),
        valuation(second_B, star_factor),
        valuation(delta, star_factor),
    ) != (2, 3, 7):
        raise ArithmeticError("finite exponent-seven place is not minimal I1*")
    residual = next(factor for factor, exponent in factors if int(exponent) == 1)
    if not residual.is_squarefree():
        raise ArithmeticError("second-q4 residual discriminant is not squarefree")
    infinity = (8 - second_A.degree(), 12 - second_B.degree(), 24 - delta.degree())
    if infinity != (2, 3, 9):
        raise ArithmeticError(("bad second-q4 infinity", infinity))
    return ring, second_A, second_B, delta, star_factor


def forced_first_marking(source_ring, source_A, source_B, values):
    """Test the exact forced P1 on the source Q80 model."""
    T_source = source_ring.gen()
    forced_x = T_source + (values["d"] - 1) * T_source**2
    radicand = forced_x**3 + source_A * forced_x + source_B
    forced_y = square_root_polynomial(radicand)
    record = {
        "x_coefficients_low_to_high": coefficient_strings(forced_x),
        "radicand_factor_degrees_and_exponents": [
            [int(factor.degree()), int(exponent)] for factor, exponent in radicand.factor()
        ],
        "rational_over_QQ": forced_y is not None,
    }
    if forced_y is not None:
        if forced_y**2 != forced_x**3 + source_A * forced_x + source_B:
            raise ArithmeticError("forced P1 square-root replay failed")
        record["y_coefficients_low_to_high"] = coefficient_strings(forced_y)
    else:
        # The T^4 coefficient is the first marking-cover squareclass and is a
        # compact exact obstruction to this orientation over QQ.
        record["leading_squareclass_T4"] = rational_record(radicand[4])
    return record


def reduce_qq_polynomial(poly, finite, target_ring):
    result = target_ring.zero()
    variable = target_ring.gen()
    for degree, value in enumerate(poly.list()):
        value = QQ(value)
        denominator = finite(ZZ(value.denominator()))
        if not denominator:
            raise ZeroDivisionError("prime divides an equation denominator")
        result += finite(ZZ(value.numerator())) / denominator * variable**degree
    return result


def canonical_polynomial_point(x_value, y_value):
    """Choose one representative of the pair {P,-P}."""
    positive = tuple(map(int, y_value.list()))
    negative = tuple(map(int, (-y_value).list()))
    return (x_value, y_value) if positive <= negative else (x_value, -y_value)


def singular_x_coordinate(cubic):
    """Return the unique singular x-coordinate of a reduced cubic fibre."""
    common = cubic.gcd(cubic.derivative())
    roots = common.roots(multiplicities=False)
    if len(roots) != 1:
        return None
    return roots[0]


def polynomial_shell(
    second_A, second_B, prime, maximum_pairs, enumerate_sections=True
):
    finite = GF(prime)
    ring = PolynomialRing(finite, "W")
    W_finite = ring.gen()
    finite_A = reduce_qq_polynomial(second_A, finite, ring)
    finite_B = reduce_qq_polynomial(second_B, finite, ring)
    delta = 4 * finite_A**3 + 27 * finite_B**2
    factorization = tuple(delta.factor())
    signature = sorted((factor.degree(), int(exponent)) for factor, exponent in factorization)
    star_factors = [
        factor.monic()
        for factor, exponent in factorization
        if int(exponent) == 7 and factor.degree() == 1
    ]
    residual_factors = [
        factor.monic()
        for factor, exponent in factorization
        if int(exponent) == 1
    ]
    if (
        len(star_factors) != 1
        or any(int(exponent) not in (1, 7) for factor, exponent in factorization)
        or sum(factor.degree() for factor in residual_factors) != 8
    ):
        raise ValueError(f"bad reduction signature {signature}")
    star_factor = star_factors[0]
    if (
        valuation(finite_A, star_factor),
        valuation(finite_B, star_factor),
        valuation(delta, star_factor),
    ) != (2, 3, 7):
        raise ValueError("finite I1* valuations changed after reduction")
    residual = prod(residual_factors, ring.one())
    if not residual.is_squarefree():
        raise ValueError("residual discriminant is not squarefree after reduction")

    if not enumerate_sections:
        return {
            "field": finite,
            "ring": ring,
            "A": finite_A,
            "B": finite_B,
            "delta": delta,
            "star_factor": star_factor,
            "curve": None,
            "shell": tuple(),
            "points": tuple(),
            "signature": signature,
        }

    signed = []
    # A P.O=0 section on this K3 chart has polynomial degrees x<=4,y<=6.
    for entries in itertools.product(finite, repeat=5):
        x_value = ring(entries)
        right_side = x_value**3 + finite_A * x_value + finite_B
        if not right_side.is_square():
            continue
        y_value = right_side.sqrt()
        if y_value.degree() > 6:
            continue
        signed.append(canonical_polynomial_point(x_value, y_value))
    shell = {}
    for x_value, y_value in signed:
        key = (tuple(map(int, x_value.list())), tuple(map(int, y_value.list())))
        shell[key] = (x_value, y_value)
    shell = tuple(shell[key] for key in sorted(shell))
    if len(shell) > maximum_pairs:
        raise ArithmeticError(
            f"unsigned polynomial shell {len(shell)} exceeds declared cap {maximum_pairs}"
        )

    function_field = ring.fraction_field()
    curve = EllipticCurve(
        function_field,
        [0, 0, 0, function_field(finite_A), function_field(finite_B)],
    )
    points = tuple(
        curve(function_field(x_value), function_field(y_value))
        for x_value, y_value in shell
    )
    return {
        "field": finite,
        "ring": ring,
        "A": finite_A,
        "B": finite_B,
        "delta": delta,
        "star_factor": star_factor,
        "curve": curve,
        "shell": shell,
        "points": points,
        "signature": signature,
    }


def direct_pole_two_system(modular, relative_sign, chart):
    """Build an unsliced P.O=2 section scheme over one good finite field.

    Write ``h=W^2+h1*W+h0``, ``x=N/h^2`` and ``y=M/h^3`` with
    ``deg(N)=8`` and ``deg(M)=12``.  The leading coefficients are normalized
    to ``l^2`` and ``relative_sign*l^3``.  The top twelve coefficients of

        M^2 = N^3 + A*N*h^4 + B*h^6

    either remain as sparse auxiliary variables or are solved recursively for
    ``M_11,...,M_0``.  This is the direct modular problem for the target
    horizontal, rather than a word search in the polynomial-section shell.
    """
    finite = modular["field"]
    base_names = (
        "h0",
        "h1",
        "l",
        "n0",
        "n1",
        "n2",
        "n3",
        "n4",
        "n5",
        "n6",
        "n7",
    )
    names = base_names + (
        tuple(f"m{degree}" for degree in range(12)) if chart == "auxiliary" else tuple()
    ) + ("sat",)
    scheme_ring = PolynomialRing(finite, names=names, order="degrevlex")
    variables = scheme_ring.gens_dict()
    fraction = scheme_ring.fraction_field()
    polynomial_ring = PolynomialRing(fraction, "W_direct")
    W_direct = polynomial_ring.gen()

    def lift_polynomial(poly):
        return polynomial_ring([fraction(value) for value in poly.list()])

    direct_A = lift_polynomial(modular["A"])
    direct_B = lift_polynomial(modular["B"])
    h0, h1, leading = map(
        fraction, (variables["h0"], variables["h1"], variables["l"])
    )
    h = W_direct**2 + h1 * W_direct + h0
    numerator = sum(
        fraction(variables[f"n{degree}"]) * W_direct**degree
        for degree in range(8)
    ) + leading**2 * W_direct**8
    square = numerator**3 + direct_A * numerator * h**4 + direct_B * h**6

    y_coefficients = [fraction.zero() for _ in range(13)]
    y_coefficients[12] = fraction(relative_sign) * leading**3
    if chart == "auxiliary":
        for degree in range(12):
            y_coefficients[degree] = fraction(variables[f"m{degree}"])
    else:
        for degree in range(23, 11, -1):
            index = degree - 12
            partial = sum(
                y_coefficients[j] * W_direct**j for j in range(13)
            )
            known_coefficient = (partial**2)[degree]
            y_coefficients[index] = (
                square[degree] - known_coefficient
            ) / (2 * y_coefficients[12])
    y_numerator = sum(
        y_coefficients[index] * W_direct**index for index in range(13)
    )
    identity = y_numerator**2 - square
    if chart == "recursive":
        if any(identity[index] for index in range(12, 25)):
            raise ArithmeticError("top-down denominator-two recursion did not close")
        equation_degrees = range(12)
    else:
        equation_degrees = range(25)
    equations = [
        scheme_ring(identity[index].numerator()) for index in equation_degrees
    ]

    star_factor = modular["star_factor"]
    star_root = finite(-star_factor[0] / star_factor[1])
    node_ring = PolynomialRing(finite, "x_node_direct")
    x_node_direct = node_ring.gen()
    nodal_cubic = (
        x_node_direct**3
        + modular["A"](star_root) * x_node_direct
        + modular["B"](star_root)
    )
    singular_x = singular_x_coordinate(nodal_cubic)
    if singular_x is None:
        raise ArithmeticError("finite star fibre has no unique rational singular x")

    # These open conditions enforce exact degree at infinity and the identity
    # component at the finite I1* fibre.  Lower-P.O cancellation solutions are
    # intentionally left in the exported closure and must be filtered after
    # solving; including the two large resultants here makes the system much
    # less useful computationally.
    open_product = (
        leading
        * h(star_root)
        * (numerator(star_root) - fraction(singular_x) * h(star_root) ** 2)
    )
    saturation = fraction(variables["sat"]) * open_product - 1
    equations.append(scheme_ring(saturation.numerator()))
    equations = tuple(equation for equation in equations if equation)
    return {
        "names": names,
        "ring": scheme_ring,
        "equations": equations,
        "star_root": int(star_root),
        "singular_x": int(singular_x),
        "relative_sign": int(relative_sign),
        "chart": chart,
        "pole_order": 2,
        "pole_location": "finite",
        "star_chart": "regular_nonnode_at_finite_I1star",
        "coefficient_term_counts": [len(equation.dict()) for equation in equations],
    }


def direct_pole_one_system(modular, relative_sign, chart):
    """Build the direct P.O=1 section scheme over one good finite field."""
    finite = modular["field"]
    base_names = (
        "z",
        "l",
        "n0",
        "n1",
        "n2",
        "n3",
        "n4",
        "n5",
    )
    names = base_names + (
        tuple(f"m{degree}" for degree in range(9)) if chart == "auxiliary" else tuple()
    ) + ("sat",)
    scheme_ring = PolynomialRing(finite, names=names, order="degrevlex")
    variables = scheme_ring.gens_dict()
    fraction = scheme_ring.fraction_field()
    polynomial_ring = PolynomialRing(fraction, "W_direct")
    W_direct = polynomial_ring.gen()

    def lift_polynomial(poly):
        return polynomial_ring([fraction(value) for value in poly.list()])

    direct_A = lift_polynomial(modular["A"])
    direct_B = lift_polynomial(modular["B"])
    pole = fraction(variables["z"])
    leading = fraction(variables["l"])
    h = W_direct - pole
    numerator = sum(
        fraction(variables[f"n{degree}"]) * W_direct**degree
        for degree in range(6)
    ) + leading**2 * W_direct**6
    square = numerator**3 + direct_A * numerator * h**4 + direct_B * h**6

    y_coefficients = [fraction.zero() for _ in range(10)]
    y_coefficients[9] = fraction(relative_sign) * leading**3
    if chart == "auxiliary":
        for degree in range(9):
            y_coefficients[degree] = fraction(variables[f"m{degree}"])
    else:
        for degree in range(17, 8, -1):
            index = degree - 9
            partial = sum(
                y_coefficients[j] * W_direct**j for j in range(10)
            )
            known_coefficient = (partial**2)[degree]
            y_coefficients[index] = (
                square[degree] - known_coefficient
            ) / (2 * y_coefficients[9])
    y_numerator = sum(
        y_coefficients[index] * W_direct**index for index in range(10)
    )
    identity = y_numerator**2 - square
    if chart == "recursive":
        if any(identity[index] for index in range(9, 19)):
            raise ArithmeticError("top-down denominator-one recursion did not close")
        equation_degrees = range(9)
    else:
        equation_degrees = range(19)
    equations = [
        scheme_ring(identity[index].numerator()) for index in equation_degrees
    ]

    # N(z) != 0 prevents cancellation of the unique pole.  No star-fibre
    # condition is imposed here: these sections are marking generators, and
    # a pole over a reducible fibre is a legitimate separate local chart.
    open_product = leading * numerator(pole)
    saturation = fraction(variables["sat"]) * open_product - 1
    equations.append(scheme_ring(saturation.numerator()))
    equations = tuple(equation for equation in equations if equation)
    return {
        "names": names,
        "ring": scheme_ring,
        "equations": equations,
        "star_root": None,
        "singular_x": None,
        "relative_sign": int(relative_sign),
        "chart": chart,
        "pole_order": 1,
        "pole_location": "finite",
        "star_chart": None,
        "coefficient_term_counts": [len(equation.dict()) for equation in equations],
    }


def direct_pole_one_infinity_system(modular, relative_sign, chart):
    """Build the P.O=1 chart whose sole intersection with O is at infinity."""
    finite = modular["field"]
    base_names = ("l", "n0", "n1", "n2", "n3", "n4", "n5")
    names = base_names + (
        tuple(f"m{degree}" for degree in range(9)) if chart == "auxiliary" else tuple()
    ) + ("sat",)
    scheme_ring = PolynomialRing(finite, names=names, order="degrevlex")
    variables = scheme_ring.gens_dict()
    fraction = scheme_ring.fraction_field()
    polynomial_ring = PolynomialRing(fraction, "W_direct")
    W_direct = polynomial_ring.gen()

    def lift_polynomial(poly):
        return polynomial_ring([fraction(value) for value in poly.list()])

    direct_A = lift_polynomial(modular["A"])
    direct_B = lift_polynomial(modular["B"])
    leading = fraction(variables["l"])
    numerator = sum(
        fraction(variables[f"n{degree}"]) * W_direct**degree
        for degree in range(6)
    ) + leading**2 * W_direct**6
    square = numerator**3 + direct_A * numerator + direct_B

    y_coefficients = [fraction.zero() for _ in range(10)]
    y_coefficients[9] = fraction(relative_sign) * leading**3
    if chart == "auxiliary":
        for degree in range(9):
            y_coefficients[degree] = fraction(variables[f"m{degree}"])
    else:
        for degree in range(17, 8, -1):
            index = degree - 9
            partial = sum(
                y_coefficients[j] * W_direct**j for j in range(10)
            )
            known_coefficient = (partial**2)[degree]
            y_coefficients[index] = (
                square[degree] - known_coefficient
            ) / (2 * y_coefficients[9])
    y_numerator = sum(
        y_coefficients[index] * W_direct**index for index in range(10)
    )
    identity = y_numerator**2 - square
    if chart == "recursive":
        if any(identity[index] for index in range(9, 19)):
            raise ArithmeticError("top-down infinity-pole recursion did not close")
        equation_degrees = range(9)
    else:
        equation_degrees = range(19)
    equations = [
        scheme_ring(identity[index].numerator()) for index in equation_degrees
    ]
    saturation = fraction(variables["sat"]) * leading - 1
    equations.append(scheme_ring(saturation.numerator()))
    equations = tuple(equation for equation in equations if equation)
    return {
        "names": names,
        "ring": scheme_ring,
        "equations": equations,
        "star_root": None,
        "singular_x": None,
        "relative_sign": int(relative_sign),
        "chart": chart,
        "pole_order": 1,
        "pole_location": "infinity",
        "star_chart": None,
        "coefficient_term_counts": [len(equation.dict()) for equation in equations],
    }


def direct_polynomial_system(modular, relative_sign, chart):
    """Build an affine P.O=0 section scheme over the finite field."""
    finite = modular["field"]
    if chart == "recursive":
        names = ("a", "n0", "n1", "n2", "n3", "sat")
        scheme_ring = PolynomialRing(finite, names=names, order="degrevlex")
        variables = scheme_ring.gens_dict()
        fraction = scheme_ring.fraction_field()
        polynomial_ring = PolynomialRing(fraction, "W_direct")
        W_direct = polynomial_ring.gen()
        direct_A = polynomial_ring([fraction(value) for value in modular["A"].list()])
        direct_B = polynomial_ring([fraction(value) for value in modular["B"].list()])
        leading = fraction(variables["a"])
        x_value = sum(
            fraction(variables[f"n{degree}"]) * W_direct**degree
            for degree in range(4)
        ) + leading**2 * W_direct**4
        square = x_value**3 + direct_A * x_value + direct_B
        y_coefficients = [fraction.zero() for _ in range(7)]
        y_coefficients[6] = fraction(relative_sign) * leading**3
        for degree in range(11, 5, -1):
            index = degree - 6
            partial = sum(
                y_coefficients[j] * W_direct**j for j in range(7)
            )
            y_coefficients[index] = (
                square[degree] - (partial**2)[degree]
            ) / (2 * y_coefficients[6])
        y_value = sum(
            y_coefficients[index] * W_direct**index for index in range(7)
        )
        identity = y_value**2 - square
        if any(identity[index] for index in range(6, 13)):
            raise ArithmeticError("top-down polynomial-section recursion did not close")
        equations = tuple(
            scheme_ring(identity[index].numerator())
            for index in range(6)
            if identity[index]
        )
        equations += (variables["sat"] * variables["a"] - 1,)
        return {
            "names": names,
            "ring": scheme_ring,
            "equations": equations,
            "star_root": None,
            "singular_x": None,
            "relative_sign": int(relative_sign),
            "chart": "identity_at_infinity_recursive",
            "pole_order": 0,
            "pole_location": None,
            "star_chart": None,
            "open_condition": "a != 0 (enforced by sat*a-1)",
            "coefficient_term_counts": [len(equation.dict()) for equation in equations],
        }

    names = tuple(f"n{degree}" for degree in range(5)) + tuple(
        f"m{degree}" for degree in range(7)
    )
    scheme_ring = PolynomialRing(finite, names=names, order="degrevlex")
    variables = scheme_ring.gens_dict()
    polynomial_ring = PolynomialRing(scheme_ring, "W_direct")
    W_direct = polynomial_ring.gen()
    direct_A = polynomial_ring(modular["A"])
    direct_B = polynomial_ring(modular["B"])
    x_value = sum(
        variables[f"n{degree}"] * W_direct**degree for degree in range(5)
    )
    y_value = sum(
        variables[f"m{degree}"] * W_direct**degree for degree in range(7)
    )
    identity = y_value**2 - x_value**3 - direct_A * x_value - direct_B
    equations = tuple(identity[index] for index in range(13) if identity[index])
    return {
        "names": names,
        "ring": scheme_ring,
        "equations": equations,
        "star_root": None,
        "singular_x": None,
        "relative_sign": None,
        "chart": "complete_polynomial",
        "pole_order": 0,
        "pole_location": None,
        "star_chart": None,
        "open_condition": None,
        "coefficient_term_counts": [len(equation.dict()) for equation in equations],
    }


def safe_parameter_tag(point):
    if point is None:
        return "uinf"
    point = QQ(point)
    if point == 0:
        return "u0"
    prefix = "um" if point < 0 else "u"
    point = abs(point)
    return f"{prefix}{point.numerator()}d{point.denominator()}"


def export_direct_msolve(system, path, prime):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write(",".join(system["names"]) + "\n")
        handle.write(str(prime) + "\n")
        for index, equation in enumerate(system["equations"]):
            handle.write(str(equation).replace("**", "^"))
            handle.write(",\n" if index + 1 < len(system["equations"]) else "\n")


def direct_msolve_probe(modular, point, prime, relative_sign):
    system_started = time.monotonic()
    if args.direct_pole_order == 0:
        system = direct_polynomial_system(modular, relative_sign, args.direct_chart)
    elif args.direct_pole_order == 1 and args.direct_pole_location == "infinity":
        system = direct_pole_one_infinity_system(
            modular, relative_sign, args.direct_chart
        )
    elif args.direct_pole_order == 1:
        system = direct_pole_one_system(modular, relative_sign, args.direct_chart)
    else:
        system = direct_pole_two_system(modular, relative_sign, args.direct_chart)
    tag = safe_parameter_tag(point)
    input_path = (
        args.direct_msolve_dir
        / (
            f"q80-third-q12-{tag}-p{prime}-po{args.direct_pole_order}-"
            f"{system.get('pole_location') or 'polynomial'}-"
            f"{args.direct_chart}"
            f"-sign{relative_sign:+d}.ms"
        )
    )
    export_direct_msolve(system, input_path, prime)
    record = {
        "relative_sign": (
            int(relative_sign) if system["relative_sign"] is not None else None
        ),
        "chart": system["chart"],
        "pole_order": int(system["pole_order"]),
        "pole_location": system.get("pole_location"),
        "star_chart": system["star_chart"],
        "open_condition": system.get("open_condition"),
        "input_path": str(input_path),
        "input_sha256": sha256(input_path),
        "variables": len(system["names"]),
        "equations": len(system["equations"]),
        "coefficient_term_counts": system["coefficient_term_counts"],
        "finite_I1star_root": system["star_root"],
        "finite_I1star_singular_x": system["singular_x"],
        "status": f"EXPORTED_DIRECT_PO{system['pole_order']}_SCHEME",
    }
    if args.run_msolve:
        executable = shutil.which("msolve")
        if executable is None:
            raise RuntimeError("--run-msolve requested but msolve is unavailable")
        solution_path = input_path.with_suffix(".solve")
        log_path = input_path.with_suffix(".log")
        command = [
            executable,
            "-t",
            str(args.msolve_threads),
            "-v",
            "2",
            "-f",
            str(input_path),
            "-o",
            str(solution_path),
        ]
        solve_started = time.monotonic()
        try:
            completed = subprocess.run(
                command,
                text=True,
                capture_output=True,
                timeout=args.msolve_timeout,
                check=False,
            )
            log_path.write_text(
                completed.stdout + completed.stderr, encoding="utf-8"
            )
            record.update(
                {
                    "status": (
                        "MSOLVE_FINISHED"
                        if completed.returncode == 0
                        else "MSOLVE_FAILED"
                    ),
                    "msolve_returncode": int(completed.returncode),
                    "msolve_runtime_seconds": time.monotonic() - solve_started,
                    "solution_path": str(solution_path),
                    "solution_exists": solution_path.exists(),
                    "solution_size_bytes": (
                        solution_path.stat().st_size if solution_path.exists() else 0
                    ),
                    "log_path": str(log_path),
                }
            )
        except subprocess.TimeoutExpired as error:
            stdout = error.stdout or ""
            stderr = error.stderr or ""
            if isinstance(stdout, bytes):
                stdout = stdout.decode(errors="replace")
            if isinstance(stderr, bytes):
                stderr = stderr.decode(errors="replace")
            log_path.write_text(stdout + stderr, encoding="utf-8")
            record.update(
                {
                    "status": "MSOLVE_TIMEOUT",
                    "msolve_runtime_seconds": time.monotonic() - solve_started,
                    "log_path": str(log_path),
                }
            )
    record["runtime_seconds"] = time.monotonic() - system_started
    return record


def rational_function_key(value):
    value = value.parent()(value)
    numerator = value.numerator()
    denominator = value.denominator()
    return (
        tuple(map(int, numerator.list())),
        tuple(map(int, denominator.list())),
    )


def point_key(point):
    return (rational_function_key(point[0]), rational_function_key(point[1]))


def po_from_x(point):
    if point.is_zero():
        return None
    denominator = point[0].denominator()
    square_root = square_root_polynomial(denominator)
    if square_root is None:
        return None
    finite_poles = ZZ(square_root.degree())
    numerator_degree = ZZ(point[0].numerator().degree())
    twice_intersection = max(2 * finite_poles, numerator_degree - 4)
    if twice_intersection < 0 or twice_intersection % 2:
        return None
    return int(twice_intersection // 2)


def target_profile(point, modular):
    """Exact equation-side profile of the generic height-eight target."""
    x_value = point[0]
    y_value = point[1]
    if po_from_x(point) != TARGET_PO:
        return None
    if y_value.denominator().degree() != 3 * TARGET_PO:
        return None
    # Identity component at infinity: after xbar=s^4*x,ybar=s^6*y the
    # section reduces to a nonzero smooth point of ybar^2=xbar^3.
    if (
        x_value.numerator().degree() - x_value.denominator().degree() != 4
        or y_value.numerator().degree() - y_value.denominator().degree() != 6
    ):
        return None

    star_factor = modular["star_factor"]
    if star_factor.degree() != 1:
        return None
    root = -star_factor[0] / star_factor[1]
    if not x_value.denominator()(root):
        return None
    node_ring = PolynomialRing(modular["field"], "x")
    node_x = node_ring.gen()
    cubic = node_x**3 + modular["A"](root) * node_x + modular["B"](root)
    singular_x = singular_x_coordinate(cubic)
    if singular_x is None:
        return None
    if x_value(root) == singular_x:
        return None

    # With identity components at both reducible fibres, Shioda's formula is
    # height=4+2(P.O)=8.  This is the exact lattice profile being matched.
    height = QQ(4 + 2 * TARGET_PO)
    if height != TARGET_HEIGHT:
        raise ArithmeticError("internal target-height profile mismatch")
    return {
        "P_dot_O": TARGET_PO,
        "height": str(height),
        "finite_I1star_identity": True,
        "infinity_I3star_identity": True,
        "x_numerator_degree": int(x_value.numerator().degree()),
        "x_denominator_degree": int(x_value.denominator().degree()),
        "y_numerator_degree": int(y_value.numerator().degree()),
        "y_denominator_degree": int(y_value.denominator().degree()),
    }


def intersection_fingerprint(point, shell_points):
    values = []
    for shell_point in shell_points:
        difference = point - shell_point
        values.append(po_from_x(difference))
    return sorted(-1 if value is None else int(value) for value in values)


def star_pair_canonical_height(point):
    """Exact Shioda height on the D7+D5 K3 using fourth multiples.

    Both component groups have exponent four.  Hence ``4*point`` is on the
    identity component at both reducible fibres and its local correction is
    zero.  Shioda's K3 formula then reads

        16*h(point) = 4 + 2*(4*point).O.

    The eighth multiple supplies an independent exact replay.
    """
    if point.is_zero():
        return QQ(0)
    fourth = 4 * point
    if fourth.is_zero():
        return QQ(0)
    fourth_po = po_from_x(fourth)
    if fourth_po is None:
        raise ArithmeticError("fourth multiple has a non-section pole divisor")
    height = QQ(4 + 2 * fourth_po) / 16
    eighth = 2 * fourth
    if eighth.is_zero():
        if height:
            raise ArithmeticError("inconsistent torsion height from fourth multiple")
        return height
    eighth_po = po_from_x(eighth)
    if eighth_po is None:
        raise ArithmeticError("eighth multiple has a non-section pole divisor")
    if QQ(4 + 2 * eighth_po) / 64 != height:
        raise ArithmeticError("fourth/eighth multiple Shioda heights disagree")
    return height


def point_linear_combination(curve, coefficients, points):
    answer = curve(0)
    for coefficient, point in zip(coefficients, points):
        if coefficient:
            answer += ZZ(coefficient) * point
    return answer


def torsion_closure(curve, generators, maximum=64):
    zero = curve(0)
    closure = {point_key(zero): zero}
    for generator in generators:
        if generator.is_zero():
            continue
        while True:
            additions = {}
            for point in closure.values():
                candidate = point + generator
                key = point_key(candidate)
                if key not in closure:
                    additions[key] = candidate
            if not additions:
                break
            closure.update(additions)
            if len(closure) > maximum:
                raise ArithmeticError("height-zero residual subgroup exceeds torsion cap")
    return tuple(closure[key] for key in sorted(closure, key=str))


def modular_height_eight_shell(modular, maximum_vectors):
    """Enumerate the complete height-eight shell in the generated MW lattice."""
    points = modular["points"]
    count = len(points)
    if not points:
        raise ArithmeticError("height-shell search needs polynomial generators")

    generator_heights = tuple(star_pair_canonical_height(point) for point in points)
    height_gram = matrix(QQ, count)
    for left in range(count):
        height_gram[left, left] = generator_heights[left]
        for right in range(left):
            sum_height = star_pair_canonical_height(points[left] + points[right])
            pairing = (
                sum_height - generator_heights[left] - generator_heights[right]
            ) / 2
            height_gram[left, right] = height_gram[right, left] = pairing
    rank = int(height_gram.rank())
    if rank < 1:
        raise ArithmeticError("polynomial shell has zero canonical-height rank")

    independent = None
    for indices in itertools.combinations(range(count), rank):
        principal = height_gram.matrix_from_rows_and_columns(indices, indices)
        if principal.det():
            independent = tuple(indices)
            break
    if independent is None:
        raise ArithmeticError("could not find a positive principal height minor")
    principal = height_gram.matrix_from_rows_and_columns(independent, independent)
    coordinates = matrix(
        QQ,
        [
            vector(QQ, [height_gram[row, column] for column in independent])
            * principal.inverse()
            for row in range(count)
        ],
    )
    if coordinates * principal * coordinates.transpose() != height_gram:
        raise ArithmeticError("height-lattice coordinate reconstruction failed")

    coordinate_denominator = lcm(
        value.denominator() for value in coordinates.list()
    )
    integral_coordinates = (
        coordinate_denominator * coordinates
    ).change_ring(ZZ)
    coordinate_hnf, generator_transform = integral_coordinates.hermite_form(
        transformation=True, include_zero_rows=False
    )
    if generator_transform * integral_coordinates != coordinate_hnf:
        raise ArithmeticError("height-lattice HNF transformation replay failed")
    lattice_basis_coordinates = coordinate_hnf / coordinate_denominator
    lattice_height_gram = (
        lattice_basis_coordinates
        * principal
        * lattice_basis_coordinates.transpose()
    )
    if lattice_height_gram.rank() != rank:
        raise ArithmeticError("generated height lattice is not positive definite")

    lattice_points = tuple(
        point_linear_combination(modular["curve"], row, points)
        for row in generator_transform.rows()
    )
    inverse_lattice_coordinates = lattice_basis_coordinates.inverse()
    torsion_residuals = []
    generator_lattice_coordinates = []
    for index, point in enumerate(points):
        lattice_coordinates = vector(QQ, coordinates[index]) * inverse_lattice_coordinates
        if not all(value in ZZ for value in lattice_coordinates):
            raise ArithmeticError("generator is nonintegral in recovered lattice basis")
        lattice_coordinates = vector(ZZ, lattice_coordinates)
        generator_lattice_coordinates.append(tuple(map(int, lattice_coordinates)))
        lifted = point_linear_combination(
            modular["curve"], lattice_coordinates, lattice_points
        )
        residual = point - lifted
        if star_pair_canonical_height(residual):
            raise ArithmeticError("height-lattice generator residual is not torsion")
        torsion_residuals.append(residual)
    torsion = torsion_closure(modular["curve"], torsion_residuals)

    height_denominator = lcm(
        value.denominator() for value in lattice_height_gram.list()
    )
    scaled_height = (height_denominator * lattice_height_gram).change_ring(ZZ)
    quadratic_form = QuadraticForm(ZZ, 2 * scaled_height)
    target_value = ZZ(height_denominator * TARGET_HEIGHT)
    shells = quadratic_form.short_vector_list_up_to_length(
        int(target_value + 1), up_to_sign_flag=True
    )
    half_vectors = tuple(vector(ZZ, row) for row in shells[int(target_value)])
    if len(half_vectors) > maximum_vectors:
        raise ArithmeticError(
            f"height-eight half-shell {len(half_vectors)} exceeds cap {maximum_vectors}"
        )

    candidate_records = []
    evaluated = {}
    for half_vector in half_vectors:
        for lattice_vector in (half_vector, -half_vector):
            base_point = point_linear_combination(
                modular["curve"], lattice_vector, lattice_points
            )
            original_coefficients = vector(ZZ, lattice_vector) * generator_transform
            for torsion_index, torsion_point in enumerate(torsion):
                candidate = base_point + torsion_point
                if candidate.is_zero():
                    continue
                key = point_key(candidate)
                if key in evaluated:
                    continue
                evaluated[key] = True
                if star_pair_canonical_height(candidate) != TARGET_HEIGHT:
                    raise ArithmeticError("enumerated lattice vector has wrong height")
                profile = target_profile(candidate, modular)
                if profile is None:
                    continue
                x_value, y_value = candidate[0], candidate[1]
                candidate_records.append(
                    {
                        "lattice_vector": list(map(int, lattice_vector)),
                        "polynomial_generator_coefficients": list(
                            map(int, original_coefficients)
                        ),
                        "torsion_index": int(torsion_index),
                        "profile": profile,
                        "intersection_fingerprint": intersection_fingerprint(
                            candidate, points
                        ),
                        "x": {
                            "numerator_coefficients_low_to_high": finite_coefficients(
                                x_value.numerator()
                            ),
                            "denominator_coefficients_low_to_high": finite_coefficients(
                                x_value.denominator()
                            ),
                        },
                        "y": {
                            "numerator_coefficients_low_to_high": finite_coefficients(
                                y_value.numerator()
                            ),
                            "denominator_coefficients_low_to_high": finite_coefficients(
                                y_value.denominator()
                            ),
                        },
                    }
                )

    return {
        "generator_heights": list(map(str, generator_heights)),
        "generator_height_gram": [
            list(map(str, row)) for row in height_gram.rows()
        ],
        "rank": rank,
        "independent_generator_indices_one_based": [index + 1 for index in independent],
        "generator_coordinates_in_lattice_basis": generator_lattice_coordinates,
        "lattice_height_gram": [
            list(map(str, row)) for row in lattice_height_gram.rows()
        ],
        "lattice_height_determinant": str(lattice_height_gram.det()),
        "torsion_subgroup_order_in_generated_shell": len(torsion),
        "height_eight_half_shell_count": len(half_vectors),
        "height_eight_points_evaluated": len(evaluated),
        "target_profile_candidates": candidate_records,
        "target_profile_candidate_count": len(candidate_records),
    }


def certified_po1_points(modular, prime):
    """Import decoded finite-pole sections and replay them on this surface."""
    if args.po1_slice_certificate is None:
        return tuple(), []
    certificate = json.loads(args.po1_slice_certificate.read_text())
    if certificate.get("schema") != "elkies-k3.q80-fixed-u-po1-msolve-slices.v1":
        raise ValueError("unexpected P.O=1 slice-certificate schema")
    if int(certificate["prime"]) != int(prime):
        raise ValueError("P.O=1 slice-certificate prime mismatch")
    imported = []
    records = []
    seen_up_to_sign = set()
    function_field = modular["curve"].base_ring()
    polynomial_ring = function_field.ring()
    for candidate in certificate.get("candidates", []):
        decoded = candidate.get("decoded_degree_one_solution")
        if not decoded or not decoded.get("scheme_equations_replayed"):
            continue

        def rational_function(record):
            numerator = polynomial_ring(
                record["numerator_coefficients_low_to_high"]
            )
            denominator = polynomial_ring(
                record["denominator_coefficients_low_to_high"]
            )
            return function_field(numerator) / function_field(denominator)

        point = modular["curve"](
            rational_function(decoded["x"]), rational_function(decoded["y"])
        )
        if po_from_x(point) != 1:
            raise ArithmeticError("imported denominator-one section has wrong P.O")
        key = min(point_key(point), point_key(-point), key=str)
        if key in seen_up_to_sign:
            continue
        seen_up_to_sign.add(key)
        height = star_pair_canonical_height(point)
        imported.append(point)
        records.append(
            {
                "pole_root": int(candidate["pole_root"]),
                "leading_coefficient": int(candidate["leading_coefficient"]),
                "canonical_height": str(height),
                "P_dot_O": 1,
                "x": decoded["x"],
                "y": decoded["y"],
                "literal_curve_substitution": True,
            }
        )
    return tuple(imported), records


def modular_horizontal_candidates(modular, maximum_word_length):
    points = modular["points"]
    signed_generators = tuple(
        (sign * point, sign * (index + 1))
        for index, point in enumerate(points)
        for sign in (1, -1)
    )
    zero = modular["curve"](0)
    current = {point_key(zero): (zero, tuple())}
    seen = dict(current)
    accepted = {}
    new_points_by_word_length = []
    for word_length in range(1, maximum_word_length + 1):
        following = {}
        for point, word in current.values():
            for generator, label in signed_generators:
                candidate = point + generator
                if candidate.is_zero():
                    continue
                key = point_key(candidate)
                candidate_word = word + (label,)
                if key not in seen and key not in following:
                    following[key] = (candidate, candidate_word)
                    profile = target_profile(candidate, modular)
                    if profile is not None:
                        accepted[key] = {
                            "point": candidate,
                            "word": candidate_word,
                            "profile": profile,
                        }
        seen.update(following)
        current = following
        new_points_by_word_length.append(len(following))
        print(
            "Q80FIXEDUTHIRDQ12|stage=modular_word_shell|"
            f"length={word_length}|new_points={len(following)}|"
            f"seen={len(seen)}|target_profiles={len(accepted)}",
            flush=True,
        )
        if not current:
            break

    records = []
    for key in sorted(accepted, key=str):
        candidate = accepted[key]
        point = candidate["point"]
        x_value, y_value = point[0], point[1]
        records.append(
            {
                "word": list(map(int, candidate["word"])),
                "profile": candidate["profile"],
                "intersection_fingerprint": intersection_fingerprint(point, points),
                "x": {
                    "numerator_coefficients_low_to_high": finite_coefficients(
                        x_value.numerator()
                    ),
                    "denominator_coefficients_low_to_high": finite_coefficients(
                        x_value.denominator()
                    ),
                },
                "y": {
                    "numerator_coefficients_low_to_high": finite_coefficients(
                        y_value.numerator()
                    ),
                    "denominator_coefficients_low_to_high": finite_coefficients(
                        y_value.denominator()
                    ),
                },
            }
        )
    return {
        "candidates": records,
        "new_points_by_word_length": new_points_by_word_length,
        "distinct_points_seen_including_zero": len(seen),
        "shell_subgroup_saturated": bool(
            new_points_by_word_length and new_points_by_word_length[-1] == 0
        ),
    }


parameter_records = []
for point in selected_u:
    parameter_started = time.monotonic()
    record = {"u": parameter_label(point), "status": "STARTED"}
    try:
        values = coefficient_values(point)
        source_ring, source_A, source_B, map_to_qq = source_equation(values)
        first_ring, first_A, first_B, first_delta = first_child_equation(map_to_qq)
        second_ring, second_A, second_B, second_delta, second_star = (
            second_child_equation(map_to_qq)
        )
        marking = forced_first_marking(source_ring, source_A, source_B, values)
        record.update(
            {
                "coefficient_values": {
                    name: rational_record(value) for name, value in values.items()
                },
                "exact_equations": {
                    "source": {
                        "base": "T",
                        "A_coefficients_low_to_high": coefficient_strings(source_A),
                        "B_coefficients_low_to_high": coefficient_strings(source_B),
                        "fibres": "I1*+I4+IV*+5I1",
                        "roots": "D5+A3+E6",
                    },
                    "first_q4": {
                        "coordinate": "U=(x-T)/T^2",
                        "A_coefficients_low_to_high": coefficient_strings(first_A),
                        "B_coefficients_low_to_high": coefficient_strings(first_B),
                        "fibres": "I5*+I5+8I1",
                        "roots": "D9+A4",
                    },
                    "second_q4": {
                        "coordinate": "W=(X-3*v^3-x1*v-x0)/v^2; v=U-d+1",
                        "A_coefficients_low_to_high": coefficient_strings(second_A),
                        "B_coefficients_low_to_high": coefficient_strings(second_B),
                        "finite_I1star_factor_coefficients_low_to_high": coefficient_strings(
                            second_star
                        ),
                        "fibres": "I3*+I1*+8I1",
                        "roots": "D7+D5",
                    },
                    "selected_pencil_coordinates_retained": True,
                    "parent_child_jacobian_maps_retained": False,
                },
                "forced_first_marking": marking,
                "modular": [],
            }
        )

        if not args.marking_only:
            for prime in primes:
                modular_started = time.monotonic()
                modular_record = {"prime": int(prime)}
                try:
                    modular = polynomial_shell(
                        second_A,
                        second_B,
                        prime,
                        args.max_polynomial_pairs,
                        enumerate_sections=(
                            not args.prime_audit_only or args.height_shell
                        ),
                    )
                    shell_search = (
                        {
                            "candidates": [],
                            "new_points_by_word_length": [],
                            "distinct_points_seen_including_zero": 0,
                            "shell_subgroup_saturated": False,
                        }
                        if args.prime_audit_only
                        else modular_horizontal_candidates(modular, args.word_length)
                    )
                    candidates = shell_search["candidates"]
                    po1_points, po1_records = certified_po1_points(modular, prime)
                    height_modular = dict(modular)
                    height_modular["points"] = modular["points"] + po1_points
                    height_shell = (
                        modular_height_eight_shell(
                            height_modular, args.max_height_shell
                        )
                        if args.height_shell
                        else None
                    )
                    height_candidates = (
                        height_shell["target_profile_candidates"]
                        if height_shell is not None
                        else []
                    )
                    direct_probes = []
                    if args.direct_msolve_dir is not None:
                        # Changing the sign of y exchanges the two possible
                        # leading signs, so +1 is the only required chart.
                        direct_probes.append(
                            direct_msolve_probe(modular, point, prime, 1)
                        )
                    modular_record.update(
                        {
                            "status": (
                                "PASS_GOOD_REDUCTION_AUDIT"
                                if args.prime_audit_only
                                and args.direct_msolve_dir is None
                                else "PASS_GOOD_REDUCTION_MODULAR_SHELL"
                            ),
                            "reduction_signature": [
                                list(map(int, value)) for value in modular["signature"]
                            ],
                            "unsigned_polynomial_section_pairs": len(
                                modular["shell"]
                            ),
                            "polynomial_shell": [
                                {
                                    "x_coefficients_low_to_high": finite_coefficients(x),
                                    "y_coefficients_low_to_high": finite_coefficients(y),
                                }
                                for x, y in modular["shell"]
                            ],
                            "target_profile_candidates": candidates,
                            "target_profile_candidate_count": (
                                len(candidates) + len(height_candidates)
                            ),
                            "shell_search": {
                                key: value
                                for key, value in shell_search.items()
                                if key != "candidates"
                            },
                            "height_eight_shell": height_shell,
                            "imported_finite_pole_P_dot_O_1_sections": po1_records,
                            "direct_section_schemes": direct_probes,
                            "direct_P_dot_O_2_schemes": (
                                direct_probes if args.direct_pole_order == 2 else []
                            ),
                        }
                    )
                except (ArithmeticError, ValueError, ZeroDivisionError) as error:
                    modular_record.update(
                        {
                            "status": "SKIP_BAD_OR_OVERSIZED_REDUCTION",
                            "reason": f"{type(error).__name__}: {error}",
                        }
                    )
                modular_record["runtime_seconds"] = time.monotonic() - modular_started
                record["modular"].append(modular_record)

        good_modular = [
            item
            for item in record["modular"]
            if item["status"]
            in (
                "PASS_GOOD_REDUCTION_AUDIT",
                "PASS_GOOD_REDUCTION_MODULAR_SHELL",
            )
        ]
        candidate_modular = [
            item for item in good_modular if item["target_profile_candidate_count"]
        ]
        if marking["rational_over_QQ"]:
            record["status"] = "PASS_EXACT_RATIONAL_FIRST_MARKING_CANDIDATE"
        elif candidate_modular:
            record["status"] = "PASS_MODULAR_THIRD_Q12_HORIZONTAL_CANDIDATES"
        else:
            record["status"] = "NO_RATIONAL_FIRST_MARKING_OR_MODULAR_HORIZONTAL_YET"
    except (ArithmeticError, ValueError, ZeroDivisionError) as error:
        record.update(
            {
                "status": "SKIP_BAD_RATIONAL_SPECIALIZATION",
                "reason": f"{type(error).__name__}: {error}",
            }
        )
    record["runtime_seconds"] = time.monotonic() - parameter_started
    parameter_records.append(record)
    print(
        "Q80FIXEDUTHIRDQ12|u={}|status={}|seconds={:.3f}".format(
            parameter_label(point), record["status"], record["runtime_seconds"]
        ),
        flush=True,
    )

exact_marking_hits = [
    record
    for record in parameter_records
    if record["status"] == "PASS_EXACT_RATIONAL_FIRST_MARKING_CANDIDATE"
]
modular_hits = [
    record
    for record in parameter_records
    if record["status"] == "PASS_MODULAR_THIRD_Q12_HORIZONTAL_CANDIDATES"
]
direct_probes = [
    probe
    for record in parameter_records
    for modular in record.get("modular", [])
    for probe in modular.get(
        "direct_section_schemes",
        modular.get("direct_P_dot_O_2_schemes", []),
    )
]

# This version stops before exact CRT/LLL reconstruction.  Keep the desired
# equation-level terminal status unreachable until the horizontal, complete
# connected D7+D5 quotient, child maps, and child fibre frame all pass.
status = (
    "PASS_FIXED_U_MODULAR_HORIZONTAL_RECONNAISSANCE"
    if modular_hits
    else (
        f"DIRECT_PO{args.direct_pole_order}_SCHEME_EXPORTED_NO_SECTION_YET"
        if direct_probes
        else "NO_FIXED_U_MARKING_IN_SELECTED_LIST"
    )
)
output = {
    "schema": "elkies-k3.q80-fixed-u-marked-third-q12-search.v1",
    "status": status,
    "target": {
        "edge": "D7+D5/MW5 --q12--> A5+A3+3A1/MW6",
        "parent_roots": "D7+D5",
        "parent_mw_rank": 5,
        "horizontal_P_dot_O": TARGET_PO,
        "horizontal_height": str(TARGET_HEIGHT),
        "horizontal_component_profile": "identity at finite I1* and infinity I3*",
        "child_roots": TARGET_ROOTS,
        "child_mw_rank": TARGET_MW_RANK,
    },
    "inputs": {
        "worker": {
            "path": "elkies-k3/scripts/certify_q80_fixed_u_marked_third_q12.sage",
            "sha256": sha256(Path(__file__).resolve()),
        },
        "parameter": {
            "path": str(args.parameter.relative_to(ROOT)),
            "sha256": sha256(args.parameter),
        },
        "universal_first_q4": {
            "path": "elkies-k3/scripts/derive_q80_first_q4_pencil.sage",
            "sha256": sha256(
                ROOT / "elkies-k3/scripts/derive_q80_first_q4_pencil.sage"
            ),
        },
        "universal_second_q4": {
            "path": "elkies-k3/scripts/derive_q80_second_q4_pencil.sage",
            "sha256": sha256(
                ROOT / "elkies-k3/scripts/derive_q80_second_q4_pencil.sage"
            ),
        },
        "global_first_marking_cover": (
            {
                "path": str(GLOBAL_FIRST_MARKING_COVER.relative_to(ROOT)),
                "sha256": sha256(GLOBAL_FIRST_MARKING_COVER),
                "use": (
                    "boundary only: the forced source section is defined on a "
                    "squarefree degree-six genus-two cover, not over QQ(u)"
                ),
            }
            if GLOBAL_FIRST_MARKING_COVER.exists()
            else None
        ),
        "simplified_first_marking_cover": (
            {
                "path": str(SIMPLIFIED_FIRST_MARKING_COVER.relative_to(ROOT)),
                "sha256": sha256(SIMPLIFIED_FIRST_MARKING_COVER),
                "use": (
                    "fallback boundary: exact simplified genus-two model and "
                    "completed bounded rational-point search"
                ),
            }
            if SIMPLIFIED_FIRST_MARKING_COVER.exists()
            else None
        ),
        "rank_five_lattice_target": (
            {
                "path": str(RANK_FIVE_LATTICE_TARGET.relative_to(ROOT)),
                "sha256": sha256(RANK_FIVE_LATTICE_TARGET),
                "use": (
                    "acceptance target only: saturated MW5 Gram determinant "
                    "237/4 and pinned height-eight coordinates (-1,1,-1,1,0)"
                ),
            }
            if RANK_FIVE_LATTICE_TARGET.exists()
            else None
        ),
    },
    "declared_search": {
        "predeclared_u": [parameter_label(point) for point in PREDECLARED_U],
        "selected_u": [parameter_label(point) for point in selected_u],
        "primes": list(map(int, primes)),
        "maximum_signed_word_length": int(args.word_length),
        "maximum_unsigned_polynomial_pairs": int(args.max_polynomial_pairs),
        "maximum_height_eight_half_shell": int(args.max_height_shell),
        "marking_only": bool(args.marking_only),
        "prime_audit_only": bool(args.prime_audit_only),
        "height_shell": bool(args.height_shell),
        "po1_slice_certificate": (
            {
                "path": str(args.po1_slice_certificate),
                "sha256": sha256(args.po1_slice_certificate),
            }
            if args.po1_slice_certificate is not None
            else None
        ),
        "direct_msolve_dir": (
            str(args.direct_msolve_dir) if args.direct_msolve_dir is not None else None
        ),
        "direct_chart": args.direct_chart,
        "direct_pole_order": int(args.direct_pole_order),
        "direct_pole_location": args.direct_pole_location,
        "run_msolve": bool(args.run_msolve),
        "msolve_threads": int(args.msolve_threads),
        "msolve_timeout": int(args.msolve_timeout),
    },
    "parameters": parameter_records,
    "summary": {
        "processed": len(parameter_records),
        "exact_rational_first_marking_hits": len(exact_marking_hits),
        "modular_horizontal_parameter_hits": len(modular_hits),
        "direct_section_schemes_exported": len(direct_probes),
        "direct_section_pole_order": int(args.direct_pole_order),
        "direct_msolve_finished": sum(
            probe["status"] == "MSOLVE_FINISHED" for probe in direct_probes
        ),
        "direct_msolve_timeouts": sum(
            probe["status"] == "MSOLVE_TIMEOUT" for probe in direct_probes
        ),
        "exact_marked_third_q12_equations": 0,
    },
    "acceptance_gate": {
        "required_terminal_status": "PASS_EXACT_FIXED_U_MARKED_Q12_A5_A3_3A1_MW6",
        "currently_passed": False,
        "remaining": [
            "recover two non-polynomial P.O>=2 directions completing the modular rank-five Shioda lattice, then match its height-eight identity/identity shell to the pinned rank-five coordinates",
            "align one horizontal across several good primes by component/intersection fingerprints",
            "projective CRT/LLL reconstruction and literal QQ(W) section substitution",
            "complete connected D7+D5 resolved quotient and exact h0=2 calculation",
            "quartic/Jacobian compilation with retained parent-child rational maps",
            "minimal fibre proof A5+A3+3A1 and Picard-rank-19 upper bound",
        ],
    },
    "fallback": {
        "triggered_for_selected_search": not exact_marking_hits and not modular_hits,
        "route": (
            "Reconstruct the marked horizontal over the persisted genus-two first-"
            "marking cover (or the smaller finite cover cut out by the horizontal), "
            "then interpolate and verify at held-out primes. Do not substitute the "
            "CM24 third-q12 boundary seed."
        ),
        "reason": (
            "The specialize-first experiment has not produced a rational forced "
            "source marking or a modular target-profile horizontal. At u=-2,p=19 "
            "the polynomial subgroup has rank three and empty height-eight shell; "
            "the complete P.O=1 finite/infinity charts add no missing direction. "
            "This is a routing decision toward P.O>=2, not a nonexistence claim."
        ),
    },
    "claim_boundary": (
        "Exact QQ specialization and both q4 equations are proved for every retained "
        "parameter. Modular shells and group-law words are finite-field reconnaissance. "
        "No third-q12 equation, generic family marking, Picard-rank upper bound, or "
        "A5+A3+3A1/MW6 child is claimed by this artifact."
    ),
    "runtime_seconds": time.monotonic() - started,
}

args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
print(
    "Q80FIXEDUTHIRDQ12|processed={}|exact_marking_hits={}|modular_hits={}|"
    "output={}|status={}".format(
        len(parameter_records),
        len(exact_marking_hits),
        len(modular_hits),
        args.output,
        status,
    ),
    flush=True,
)
