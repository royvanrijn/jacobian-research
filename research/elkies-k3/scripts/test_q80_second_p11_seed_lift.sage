#!/usr/bin/env sage
"""Test lifting of selected fully profile/pair-gated q=80 seeds.

Unlike the disproved standard-P1 branch, all four ambient parameters are free.
P3 is encoded directly with its visible simple pole.  The test is exact over
ZZ and reports Jacobian rank and compatibility through 11^2; it does not claim
a characteristic-zero model unless subsequent lifts and reconstruction pass.
"""

from sage.all import *
import argparse
import random
from itertools import product as itertools_product


parser = argparse.ArgumentParser()
parser.add_argument("--seed", choices=("correct7", "second11"), default="second11")
parser.add_argument("--jacobian-only", action="store_true")
parser.add_argument("--first-obstruction", action="store_true")
parser.add_argument("--greedy", action="store_true")
parser.add_argument("--greedy-after", type=int, default=None)
parser.add_argument("--beam-after", type=int, default=None)
parser.add_argument("--beam-size", type=int, default=20)
parser.add_argument("--random-seed", type=int, default=1)
parser.add_argument("--exponent", type=int, default=7)
parser.add_argument("--stop-after", type=int, default=None)
args = parser.parse_args()
random.seed(int(args.random_seed))
prime = 7 if args.seed == "correct7" else 11
field = GF(prime)
use_pair_gates = args.seed == "correct7"
names = (
    ("d", "p", "q", "r", "b1", "b2", "b3", "b4", "e", "rho", "w")
    + ("a", "b", "u0", "u1", "u2")
    + ("k2", "l2", "x21", "x22", "x23", "y21", "y22", "y23", "y24", "y25")
    + ("z",)
    + tuple(f"x3{index}" for index in range(7))
    + tuple(f"y3{index}" for index in range(10))
    + (
        tuple(f"c12_{index}" for index in range(2))
        + tuple(f"c13_{index}" for index in range(4))
        + tuple(f"c23_{index}" for index in range(6))
        if use_pair_gates else ()
    )
)
ring = PolynomialRing(ZZ, names=names)
v = ring.gens_dict()
polynomials = PolynomialRing(ring, "T")
T = polynomials.gen()

A = T**2*(-3+v["p"]*T+v["q"]*T**2+v["r"]*T**3)
B = T**3*(
    2+v["b1"]*T+v["b2"]*T**2+v["b3"]*T**3+v["b4"]*T**4+v["e"]*T**5
)
discriminant = 4*A**3+27*B**2
equations = [A(1)+3*v["d"]**2, B(1)-2*v["d"]**3]
equations += [discriminant.derivative(order)(1) for order in range(1, 4)]
equations += [
    v["w"]**3+A(v["rho"])*v["w"]+B(v["rho"]),
    3*v["w"]**2+A(v["rho"]),
    discriminant.derivative()(v["rho"]),
]

X1 = T*(v["a"]+v["b"]*T+(v["d"]-v["a"]-v["b"])*T**2)
Y1 = T**2*(T-1)*(v["u0"]+v["u1"]*T+v["u2"]*T**2)
equations += (Y1**2-X1**3-A*X1-B).list()
equations += [X1(v["rho"])-v["w"], Y1(v["rho"])]

X2 = v["k2"]**2+v["x21"]*T+v["x22"]*T**2+v["x23"]*T**3+v["l2"]**2*T**4
Y2 = (
    v["k2"]**3+v["y21"]*T+v["y22"]*T**2+v["y23"]*T**3
    +v["y24"]*T**4+v["y25"]*T**5+v["l2"]**3*T**6
)
equations += (Y2**2-X2**3-A*X2-B).list()

Z3 = T-v["z"]
X3 = sum(v[f"x3{index}"]*T**index for index in range(7))
Y3 = sum(v[f"y3{index}"]*T**index for index in range(10))
equations += (Y3**2-X3**3-A*X3*Z3**4-B*Z3**6).list()
equations += [
    X3(1)-v["d"]*Z3(1)**2,
    Y3(1),
    X3(v["rho"])-v["w"]*Z3(v["rho"])**2,
    Y3(v["rho"]),
]


def chord_data(XP, YP, ZP, XQ, YQ, ZQ):
    D = XP*ZQ**2-XQ*ZP**2
    S = YP*ZQ**3+YQ*ZP**3
    H = D*ZP*ZQ
    N = S**2-D**2*(XP*ZQ**2+XQ*ZP**2)
    return H, N


if use_pair_gates:
    one = polynomials.one()
    for label, degree, P, Q in (
        ("12", 2, (X1, Y1, one), (X2, Y2, one)),
        ("13", 4, (X1, Y1, one), (X3, Y3, Z3)),
        ("23", 6, (X2, Y2, one), (X3, Y3, Z3)),
    ):
        cancellation = T**degree+sum(
            v[f"c{label}_{index}"]*T**index for index in range(degree)
        )
        H, N = chord_data(*P, *Q)
        remainder_H = H.quo_rem(cancellation)[1]
        remainder_N = N.quo_rem(cancellation**2)[1]
        equations += remainder_H.list()+remainder_N.list()
equations = tuple(ring(equation) for equation in equations if equation)

core_seed_values = {
    "correct7": (
        (3, 4, 3, 4, 0, 2, 6, 0, 2, 4, 1)
        + (1, 2, 2, 3, 0)
        + (6, 2, 2, 2, 6, 4, 4, 5, 4, 4)
        + (6,)
        + (1, 4, 3, 0, 6, 1, 4)
        + (1, 6, 2, 5, 3, 4, 1, 2, 3, 1)
    ),
    "second11": (
        (3, 6, 2, 1, 5, 7, 10, 7, 1, 4, 6)
        + (9, 5, 7, 1, 0)
        + (10, 7, 2, 10, 2, 8, 7, 10, 7, 10)
        + (9,)
        + (4, 7, 6, 10, 3, 7, 1)
        + (3, 1, 2, 4, 6, 1, 5, 5, 5, 1)
    ),
}[args.seed]


def finite_polynomial(coefficients):
    finite_ring = PolynomialRing(field, "t")
    t = finite_ring.gen()
    return sum(field(value)*t**index for index, value in enumerate(coefficients))


seed_values = core_seed_values
if use_pair_gates:
    finite_ring = PolynomialRing(field, "t")
    t = finite_ring.gen()
    P1 = (
        finite_polynomial((0, 1, 2, 0)),
        finite_polynomial((0, 0, 5, 6, 3, 0)),
        finite_ring.one(),
    )
    P2 = (
        finite_polynomial((1, 2, 2, 6, 4)),
        finite_polynomial((6, 4, 4, 5, 4, 4, 1)),
        finite_ring.one(),
    )
    P3 = (
        finite_polynomial((1, 4, 3, 0, 6, 1, 4)),
        finite_polynomial((1, 6, 2, 5, 3, 4, 1, 2, 3, 1)),
        t-6,
    )

    def finite_cancellation(P, Q, expected_degree):
        XP, YP, ZP = P
        XQ, YQ, ZQ = Q
        D = XP*ZQ**2-XQ*ZP**2
        S = YP*ZQ**3+YQ*ZP**3
        H = D*ZP*ZQ
        N = S**2-D**2*(XP*ZQ**2+XQ*ZP**2)
        first = gcd(H, N)
        cancellation = gcd(H, N//first).monic()
        assert cancellation.degree() == expected_degree
        return tuple(ZZ(cancellation[index]) for index in range(expected_degree))

    seed_values += (
        finite_cancellation(P1, P2, 2)
        + finite_cancellation(P1, P3, 4)
        + finite_cancellation(P2, P3, 6)
    )
seed = vector(ZZ, seed_values)
assert len(seed) == len(names)
substitution = dict(zip(ring.gens(), seed))
residuals = vector(ZZ, [equation(*seed) for equation in equations])
assert all(value % prime == 0 for value in residuals)

substitution_field = dict(zip(ring.gens(), map(field, seed)))
jacobian = matrix(
    field,
    [[equation.derivative(variable).subs(substitution_field) for variable in ring.gens()]
     for equation in equations],
)
rhs = vector(field, [field(-(value//prime)) for value in residuals])
rank = jacobian.rank()
augmented_rank = jacobian.augment(rhs).rank()
liftable = rank == augmented_rank
print(
    f"Q80SEEDLIFT|seed={args.seed}|prime={prime}|variables={len(names)}|equations={len(equations)}|"
    f"terms={sum(len(equation.dict()) for equation in equations)}|"
    f"max_terms={max(len(equation.dict()) for equation in equations)}|"
    f"rank={rank}|kernel={len(names)-rank}|liftable_to_121={ZZ(liftable)}",
    flush=True,
)
if liftable:
    kernel = jacobian.right_kernel().basis()
    for index, tangent in enumerate(kernel):
        support = tuple(
            f"{name}:{ZZ(value)}"
            for name, value in zip(names, tangent)
            if value
        )
        print(
            f"Q80SEEDLIFT|seed={args.seed}|stage=kernel|index={index}|"
            f"support={support}",
            flush=True,
        )
    if args.jacobian_only:
        args.exponent = 1

    integer_polynomials = PolynomialRing(ZZ, "t")
    t_integer = integer_polynomials.gen()
    base_equation_count = 52

    def values_at(point):
        base_values = [equation(*point) for equation in equations[:base_equation_count]]
        if not use_pair_gates:
            return vector(ZZ, base_values)
        point_by_name = dict(zip(names, map(ZZ, point)))
        d_value = point_by_name["d"]
        X1_value = t_integer*(
            point_by_name["a"]+point_by_name["b"]*t_integer
            +(d_value-point_by_name["a"]-point_by_name["b"])*t_integer**2
        )
        Y1_value = t_integer**2*(t_integer-1)*(
            point_by_name["u0"]+point_by_name["u1"]*t_integer
            +point_by_name["u2"]*t_integer**2
        )
        X2_value = (
            point_by_name["k2"]**2+point_by_name["x21"]*t_integer
            +point_by_name["x22"]*t_integer**2
            +point_by_name["x23"]*t_integer**3
            +point_by_name["l2"]**2*t_integer**4
        )
        Y2_value = (
            point_by_name["k2"]**3+point_by_name["y21"]*t_integer
            +point_by_name["y22"]*t_integer**2
            +point_by_name["y23"]*t_integer**3
            +point_by_name["y24"]*t_integer**4
            +point_by_name["y25"]*t_integer**5
            +point_by_name["l2"]**3*t_integer**6
        )
        Z3_value = t_integer-point_by_name["z"]
        X3_value = sum(
            point_by_name[f"x3{index}"]*t_integer**index for index in range(7)
        )
        Y3_value = sum(
            point_by_name[f"y3{index}"]*t_integer**index for index in range(10)
        )
        one_value = integer_polynomials.one()
        pair_values = []
        for label, degree, P, Q in (
            ("12", 2, (X1_value, Y1_value, one_value), (X2_value, Y2_value, one_value)),
            ("13", 4, (X1_value, Y1_value, one_value), (X3_value, Y3_value, Z3_value)),
            ("23", 6, (X2_value, Y2_value, one_value), (X3_value, Y3_value, Z3_value)),
        ):
            cancellation = t_integer**degree+sum(
                point_by_name[f"c{label}_{index}"]*t_integer**index
                for index in range(degree)
            )
            H_value, N_value = chord_data(*P, *Q)
            remainder_H = H_value.quo_rem(cancellation)[1]
            remainder_N = N_value.quo_rem(cancellation**2)[1]
            pair_values += [remainder_H[index] for index in range(degree)]
            pair_values += [remainder_N[index] for index in range(2*degree)]
        answer = vector(ZZ, base_values+pair_values)
        assert len(answer) == len(equations)
        return answer

    assert values_at(seed) == residuals

    def compatible(rhs_value):
        return not any(jacobian.left_kernel_matrix()*rhs_value)

    if args.first_obstruction:
        # The obstruction to extending a first correction
        # particular+s*kernel[0]+t*kernel[1] from p^2 to p^3 is quadratic in
        # (s,t).  Recover its exact vector-valued polynomial over GF(p), then
        # row-reduce the coefficient vectors to print only independent scalar
        # obstruction equations.
        left_kernel = jacobian.left_kernel_matrix()
        particular = jacobian.solve_right(rhs)
        samples = []
        obstruction_values = []
        for s_value in field:
            for t_value in field:
                delta = vector(field, particular+s_value*kernel[0]+t_value*kernel[1])
                candidate = seed+prime*vector(ZZ, map(ZZ, delta))
                candidate_values = values_at(candidate)
                assert all(value % prime**2 == 0 for value in candidate_values)
                next_rhs = vector(
                    field, [field(-(value//prime**2)) for value in candidate_values]
                )
                samples.append((s_value, t_value))
                obstruction_values.append(left_kernel*next_rhs)
        evaluation = matrix(
            field,
            [
                [1, s_value, t_value, s_value**2, s_value*t_value, t_value**2]
                for s_value, t_value in samples
            ],
        )
        coefficient_rows = []
        for obstruction_index in range(left_kernel.nrows()):
            values = vector(
                field,
                [value[obstruction_index] for value in obstruction_values],
            )
            coefficients = evaluation.solve_right(values)
            assert evaluation*coefficients == values
            coefficient_rows.append(coefficients)
        independent = matrix(field, coefficient_rows).row_space().basis_matrix()
        zero_locus = []
        for sample, value in zip(samples, obstruction_values):
            if not any(value):
                zero_locus.append(tuple(map(ZZ, sample)))
        monomial_names = ("1", "s", "t", "s^2", "s*t", "t^2")
        for index, row in enumerate(independent.rows(), 1):
            terms = tuple(
                f"{ZZ(coefficient)}*{name}"
                for coefficient, name in zip(row, monomial_names)
                if coefficient
            )
            print(
                f"Q80SEEDLIFT|seed={args.seed}|stage=first_obstruction|"
                f"equation={index}|polynomial={'+'.join(terms) if terms else '0'}",
                flush=True,
            )
        print(
            f"Q80SEEDLIFT|seed={args.seed}|stage=first_obstruction|"
            f"rank={independent.nrows()}|zero_locus={tuple(zero_locus)}",
            flush=True,
        )
        args.exponent = 1

    states = [seed]
    target_exponent = args.exponent
    for exponent in range(1, target_exponent):
        modulus = prime**exponent
        next_states = []
        greedy_now = args.greedy or (
            args.greedy_after is not None and exponent >= args.greedy_after
        )
        beam_now = args.beam_after is not None and exponent >= args.beam_after
        state_order = list(states)
        if beam_now:
            random.shuffle(state_order)
        diversify_beam = beam_now and len(state_order) > args.beam_size
        for state in state_order:
            values = values_at(state)
            assert all(value % modulus == 0 for value in values)
            rhs_state = vector(
                field, [field(-(value//modulus)) for value in values]
            )
            if not compatible(rhs_state):
                continue
            particular = jacobian.solve_right(rhs_state)

            def candidate_for(coefficients):
                delta = vector(field, particular)
                for coefficient, tangent in zip(coefficients, kernel):
                    delta += coefficient*tangent
                return state+modulus*vector(
                    ZZ, [ZZ(value) for value in delta]
                )

            coefficient_tuples = tuple(
                itertools_product(field, repeat=len(kernel))
            )
            next_modulus = modulus*prime
            needs_lookahead = exponent+1 < target_exponent
            if needs_lookahead and exponent >= 2:
                zero_coefficients = tuple(field.zero() for _ in kernel)
                zero_candidate = candidate_for(zero_coefficients)
                zero_values = values_at(zero_candidate)
                assert all(value % next_modulus == 0 for value in zero_values)
                zero_rhs = vector(
                    field, [field(-(value//next_modulus)) for value in zero_values]
                )
                slopes = []
                for index in range(len(kernel)):
                    unit = list(zero_coefficients)
                    unit[index] = field.one()
                    unit_values = values_at(candidate_for(tuple(unit)))
                    assert all(value % next_modulus == 0 for value in unit_values)
                    unit_rhs = vector(
                        field,
                        [field(-(value//next_modulus)) for value in unit_values],
                    )
                    slopes.append(unit_rhs-zero_rhs)
                for coefficients in coefficient_tuples:
                    lookahead_rhs = vector(field, zero_rhs)
                    for coefficient, slope in zip(coefficients, slopes):
                        lookahead_rhs += coefficient*slope
                    if not compatible(lookahead_rhs):
                        continue
                    next_states.append(candidate_for(coefficients))
                    if greedy_now or diversify_beam:
                        break
            else:
                for coefficients in coefficient_tuples:
                    candidate = candidate_for(coefficients)
                    candidate_values = values_at(candidate)
                    assert all(value % next_modulus == 0 for value in candidate_values)
                    if needs_lookahead:
                        lookahead_rhs = vector(
                            field,
                            [field(-(value//next_modulus)) for value in candidate_values],
                        )
                        if not compatible(lookahead_rhs):
                            continue
                    next_states.append(candidate)
                    if greedy_now or diversify_beam:
                        break
            if greedy_now and next_states:
                break
            if beam_now and len(next_states) >= args.beam_size:
                break
        states = next_states
        print(
            f"Q80SEEDLIFT|seed={args.seed}|stage=lift|exponent={exponent+1}|"
            f"survivors={len(states)}",
            flush=True,
        )
        if not states:
            break
        if args.stop_after is not None and exponent+1 >= args.stop_after:
            break
        if len(states) > 5000:
            print(f"Q80SEEDLIFT|seed={args.seed}|stage=lift|status=POSITIVE_DIMENSIONAL_BEAM_CAP")
            break
    print(
        f"Q80SEEDLIFT|seed={args.seed}|SUMMARY|status="
        + (
            "JACOBIAN_ONLY"
            if args.jacobian_only
            else ("LIFTS_THROUGH_P_ADIC_DIGITS" if states else "HIGHER_OBSTRUCTION")
        ),
        flush=True,
    )
else:
    print(f"Q80SEEDLIFT|seed={args.seed}|SUMMARY|status=FIRST_LIFT_OBSTRUCTED", flush=True)
