#!/usr/bin/env python3
"""Prove the exact stable W_3(F_2) lift degree is nineteen by exact SAT.

The system retains every coefficient of an arbitrary degree-18 first Witt
correction that can affect the z^0 or z^1 determinant digits.  At the second
Witt digit, degree-18 corrections reach plane degree at most 27 in z^0 and
26 in z^1.  Requiring all terms above those bounds to vanish gives an exact
necessary Boolean system.  Z3 proves that system unsatisfiable.

With ``--degree 19 --seek-lift``, the script instead solves the complete
determinant equations in a z-linear first-correction ansatz and directly
replays the returned degree-19 polynomial map modulo eight.

The fixed-first second correction is also audited as an affine F_2 system.
Its independent z0, z1, and z2 blocks can be minimized separately with exact
pseudo-Boolean bisection.  Incidence-component decomposition makes the
fixed-first minimum exact by reducing it to small independent decoding jobs.

With ``--require-w4-degree 19``, the same nonlinear master is coupled to an
exact affine W_4 completion subproblem.  Inconsistent subproblems yield
replayable dual cuts; singleton zero rows are batched as codomain-hole
constraints.  A bounded or timed-out run remains an experiment, not a proof
of UNSAT.  ``--w4-quotient-compiler`` factors selected coefficients through
shared projected two-by-two Jacobian minors; structural seeds can impose the
complete monomial complement or one output z-layer.

The mathematical derivation and scope boundary are in
verified/HUQ_KURUVILLA_PLANE_W2_OBSTRUCTION.md.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from functools import reduce
import hashlib
import itertools
import json
from pathlib import Path

import z3

from affine_support_decoder import (
    AffineComponent,
    affine_inconsistency_certificate,
    affine_component_statistics,
    affine_components,
    minimize_affine_components,
    rank_affine_rows,
    solve_affine_rows,
)


Monomial = tuple[int, int]
BooleanPolynomial = dict[Monomial, list[z3.BoolRef]]
Monomial3 = tuple[int, int, int]
Polynomial3 = dict[Monomial3, int]
SymbolicPolynomial3 = dict[
    Monomial3, list[tuple[int, z3.BoolRef | None]]
]
BitVectorPolynomial3 = dict[Monomial3, z3.BitVecRef]

parser = argparse.ArgumentParser()
parser.add_argument(
    "--degree",
    type=int,
    default=18,
    help="first-correction degree to test (18 is the proved theorem)",
)
parser.add_argument(
    "--show-model",
    action="store_true",
    help="print the true coefficient supports when the system is satisfiable",
)
parser.add_argument(
    "--seek-lift",
    action="store_true",
    help=(
        "solve the complete determinant equations in the ansatz where the "
        "first correction is at most linear in z"
    ),
)
parser.add_argument(
    "--write-certificate",
    type=Path,
    help="write the returned degree-19 support model as a generated JSON certificate",
)
parser.add_argument(
    "--replay-certificate",
    type=Path,
    help="replay the pinned degree-19 JSON support certificate without solving",
)
parser.add_argument(
    "--minimize-second",
    type=Path,
    metavar="CERTIFICATE",
    help=(
        "fix the first Witt digit from CERTIFICATE and prove the minimum "
        "support of the second correction by cardinality bisection"
    ),
)
parser.add_argument(
    "--second-support-bound",
    type=int,
    help="with --minimize-second, test this support bound before bisection",
)
parser.add_argument(
    "--second-support-layer",
    choices=("all", "z0", "z1", "z2"),
    default="all",
    help="with --minimize-second, minimize only the selected independent z-layer",
)
parser.add_argument(
    "--audit-second-linear",
    action="store_true",
    help=(
        "with --minimize-second, compute the exact rank and nullity of the "
        "fixed-first affine F_2 completion system"
    ),
)
parser.add_argument(
    "--component-minimize-second",
    action="store_true",
    help=(
        "with --minimize-second, minimize every affine incidence component "
        "in all layers or in the selected z-layer"
    ),
)
parser.add_argument(
    "--component-determinant",
    choices=("one", "any"),
    default="one",
    help=(
        "with --component-minimize-second, require determinant one or allow "
        "any odd constant modulo eight"
    ),
)
parser.add_argument(
    "--first-support-bound",
    type=int,
    help=(
        "with --degree 19 --seek-lift, bound the six first-correction "
        "supports in the full nonlinear system"
    ),
)
parser.add_argument(
    "--decode-found-first",
    action="store_true",
    help=(
        "after a successful --first-support-bound search, minimize the "
        "second correction componentwise for that returned first digit"
    ),
)
parser.add_argument(
    "--require-w4-class-zero",
    action="store_true",
    help=(
        "add every odd-odd-odd coefficient equation for vanishing of the "
        "next determinant-digit class in H^3_dR"
    ),
)
parser.add_argument(
    "--audit-w4-class-certificate",
    type=Path,
    metavar="CERTIFICATE",
    help=(
        "compile the next-class equations and replay their values on a "
        "degree-19 W_3 support certificate without solving"
    ),
)
parser.add_argument(
    "--w4-fix-first-certificate",
    type=Path,
    metavar="CERTIFICATE",
    help=(
        "with --require-w4-class-zero, fix only the first Witt digit from "
        "CERTIFICATE before solving"
    ),
)
parser.add_argument(
    "--require-w4-degree",
    type=int,
    choices=(19,),
    help=(
        "add exact necessary equations for a W_4 correction of this degree; "
        "currently the degree-19 boundary system is supported"
    ),
)
parser.add_argument(
    "--w4-cut-limit",
    type=int,
    default=16,
    help="maximum certified completion cuts in the degree-19 joint search",
)
parser.add_argument(
    "--w4-seed-cuts",
    type=int,
    default=8,
    help="number of singleton codomain-hole cuts to compile initially",
)
parser.add_argument(
    "--w4-seed-source",
    choices=("pinned", "structural"),
    default="pinned",
    help=(
        "seed singleton holes activated by the pinned witness or from the "
        "entire structural target complement"
    ),
)
parser.add_argument(
    "--w4-structural-z-layer",
    type=int,
    choices=(0, 1, 2, 3),
    help="with structural seeds, restrict codomain holes to this output z-layer",
)
parser.add_argument(
    "--write-w4-dimacs",
    type=Path,
    help="write the seeded degree-19 joint Boolean system as DIMACS and exit",
)
parser.add_argument(
    "--w4-quotient-compiler",
    action="store_true",
    help=(
        "compile selected W_4 quotient coefficients through shared factored "
        "two-by-two minors instead of cubic term expansion"
    ),
)
parser.add_argument(
    "--solve-equations-first",
    action="store_true",
    help="run Z3 solve-eqs before bit-blasting the Boolean/BitVec master",
)
parser.add_argument("--timeout-ms", type=int, default=60_000)
parser.add_argument(
    "--random-seed",
    type=int,
    default=0,
    help="SAT branching seed for reproducible Pareto-frontier sampling",
)
args = parser.parse_args()
assert args.degree >= 1
assert args.write_certificate is None or args.seek_lift
assert args.replay_certificate is None or not args.seek_lift
assert args.minimize_second is None or (args.seek_lift and args.degree == 19)
assert args.second_support_bound is None or args.minimize_second is not None
assert not args.audit_second_linear or args.minimize_second is not None
assert args.second_support_layer == "all" or args.minimize_second is not None
assert args.second_support_layer == "all" or not args.audit_second_linear
assert not args.component_minimize_second or args.minimize_second is not None
assert not args.component_minimize_second or not args.audit_second_linear
assert args.component_determinant == "one" or args.component_minimize_second
assert args.first_support_bound is None or args.seek_lift
assert args.first_support_bound is None or args.minimize_second is None
assert not args.decode_found_first or args.first_support_bound is not None
assert not args.require_w4_class_zero or (args.seek_lift and args.degree == 19)
assert args.audit_w4_class_certificate is None or (
    args.seek_lift and args.degree == 19
)
assert args.w4_fix_first_certificate is None or args.require_w4_class_zero
assert args.w4_fix_first_certificate is None or args.minimize_second is None
assert args.audit_w4_class_certificate is None or not args.require_w4_class_zero
assert not args.require_w4_class_zero or args.minimize_second is None
assert not args.require_w4_class_zero or not args.decode_found_first
assert args.require_w4_degree is None or (args.seek_lift and args.degree == 19)
assert args.require_w4_degree is None or args.minimize_second is None
assert args.require_w4_degree is None or not args.decode_found_first
assert args.w4_cut_limit >= 1
assert 1 <= args.w4_seed_cuts <= args.w4_cut_limit
assert args.write_w4_dimacs is None or args.require_w4_degree == 19
assert not args.w4_quotient_compiler or args.require_w4_degree == 19
assert args.w4_structural_z_layer is None or args.w4_seed_source == "structural"
if args.require_w4_degree is not None:
    # Every W_4 completion must first have zero top de Rham obstruction.
    args.require_w4_class_zero = True

AFFINE_DECODER_SHA256 = (
    "cd829ec30aca757c052fadf936c5cfb55abbf6622e09917e7bae841f263e834f"
)
assert hashlib.sha256(
    Path(__file__).with_name("affine_support_decoder.py").read_bytes()
).hexdigest() == AFFINE_DECODER_SHA256
PREFERRED_W4_TOP_CLASS_SHA256 = (
    "e219eeab3de8badeaf76c9cb393a1b1f0d8a791ae794ef66e97fa4b70c77b9fb"
)

FROZEN_CERTIFICATE_SHA256 = {
    "huq_kuruvilla_w3_degree19_witness.json": (
        "6f84a06494efbaa653e3b6471e7023d1ac2fada8e0513b6f84a12cf1be0c5246"
    ),
    "huq_kuruvilla_w3_degree19_witness_second_reduced.json": (
        "00892eac3da54c4ffb59da1d90746b44e02c61e72727854d9ae07b0c126c32e0"
    ),
    "huq_kuruvilla_w3_degree19_witness_second_280.json": (
        "4eb9a9db42beb7f7b17d55aa17d4e103f97bfcb70d2b2ad082da35492484119e"
    ),
    "huq_kuruvilla_w3_degree19_witness_second_180.json": (
        "4a0807522c4bbd7f75b540391b6058e59bbb1e69f4de6665a51a1ca0299a65af"
    ),
    "huq_kuruvilla_w3_degree19_witness_second_177.json": (
        "bd3d960e2134eb111d2fda977115a431ea53f838c8d8a03ee7fc80a47e7708fb"
    ),
    "huq_kuruvilla_w3_degree19_witness_second_172.json": (
        "17540744b69fc43c1dd9b27cb81cdd524fe6b47ff9a34c8d4919b30fe1a5109c"
    ),
    "huq_kuruvilla_w3_degree19_witness_second_165.json": (
        "ccdb0b58ddecc0ac3a2eec73c3d9d1a97ecce8d5a49072f83803de90f2944150"
    ),
    "huq_kuruvilla_w3_degree19_witness_first_280.json": (
        "e3d08551772afcbb90217302f91ab022ebb60b72b740ae104fec6bf525ad7f67"
    ),
    "huq_kuruvilla_w3_degree19_witness_first_280_second_160.json": (
        "a79984550854ce01d903783156baa6f7d4720f56ec2824bd5823cd13088a5d7f"
    ),
    "huq_kuruvilla_w3_degree19_w4_class_zero.json": (
        "9ad15068593af7cca87169c25eed2ff53068cc466d183ce320c2fc7d0e2c1aaa"
    ),
}

K_SUPPORT = {
    (1, 1),
    (5, 2),
    (7, 1),
    (9, 3),
    (10, 1),
    (11, 2),
    (12, 3),
    (13, 4),
}
K_INTEGER = {
    (1, 1): 1,
    (3, 0): 2,
    (5, 2): 3,
    (6, 0): -2,
    (7, 1): -1,
    (9, 0): 2,
    (9, 3): -1,
    (10, 1): -1,
    (11, 2): 3,
    (12, 3): -1,
    (13, 4): 1,
}
J_INTEGER = {(0, 0): 1}
for _monomial, _coefficient in K_INTEGER.items():
    J_INTEGER[_monomial] = J_INTEGER.get(_monomial, 0) + 2 * _coefficient

PX = ((0, 0),)
PY = ((2, 0),)
QX = ((4, 0), (6, 2))
QY = ((0, 0), (6, 0), (8, 2))

PX_INTEGER = {(0, 0): 1, (1, 1): 2, (3, 0): 4, (5, 2): 6}
PY_INTEGER = {(2, 0): 1, (6, 1): 2}
QX_INTEGER = {(4, 0): 5, (5, 1): 6, (6, 2): 7, (7, 3): 8}
QY_INTEGER = {(0, 0): 1, (6, 0): 1, (7, 1): 2, (8, 2): 3}


def monomials(maximum_degree: int) -> list[Monomial]:
    return [
        (i, total - i)
        for total in range(maximum_degree + 1)
        for i in range(total + 1)
    ]


def xor(terms: list[z3.BoolRef]) -> z3.BoolRef:
    if not terms:
        return z3.BoolVal(False)
    return reduce(z3.Xor, terms)


def affine_boolean_support(
    expression: z3.BoolRef,
    variable_index: dict[int, int],
) -> tuple[bool, set[int]]:
    """Expand a simplified affine Boolean expression as constant XOR vars."""

    if z3.is_true(expression):
        return True, set()
    if z3.is_false(expression):
        return False, set()
    index = variable_index.get(expression.get_id())
    if index is not None:
        return False, {index}

    kind = expression.decl().kind()
    if kind == z3.Z3_OP_NOT:
        constant, support = affine_boolean_support(
            expression.arg(0), variable_index
        )
        return not constant, support
    if kind in (z3.Z3_OP_XOR, z3.Z3_OP_EQ):
        constant = kind == z3.Z3_OP_EQ
        support: set[int] = set()
        for child in expression.children():
            child_constant, child_support = affine_boolean_support(
                child, variable_index
            )
            constant ^= child_constant
            support.symmetric_difference_update(child_support)
        return constant, support
    raise AssertionError(f"non-affine residual expression: {expression}")


def audit_affine_boolean_system(
    assertions: list[z3.BoolRef],
    variables_to_audit: list[z3.BoolRef],
    substitutions: list[tuple[z3.BoolRef, z3.BoolRef]],
    variable_layers: list[str],
) -> tuple[
    int,
    int,
    int,
    dict[str, tuple[int, int, int]],
    dict[str, list[tuple[int, int, int, int]]],
    list[AffineComponent],
]:
    """Return equation count, GF(2) rank, and affine solution dimension."""

    assert len(variables_to_audit) == len(variable_layers)
    variable_index = {
        variable.get_id(): index
        for index, variable in enumerate(variables_to_audit)
    }
    combined = z3.simplify(
        z3.substitute(z3.And(*assertions), *substitutions)
    )
    reduced_assertions = (
        list(combined.children())
        if combined.decl().kind() == z3.Z3_OP_AND
        else [combined]
    )
    rows: list[tuple[int, int]] = []
    rows_by_layer: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for reduced in reduced_assertions:
        constant, support = affine_boolean_support(reduced, variable_index)
        if not support:
            assert constant, reduced
            continue
        left = sum(1 << index for index in support)
        row = (left, int(not constant))
        rows.append(row)
        layers = {variable_layers[index] for index in support}
        assert len(layers) == 1, f"mixed second-correction layers: {layers}"
        rows_by_layer[layers.pop()].append(row)

    rank = rank_affine_rows(rows)
    layer_statistics: dict[str, tuple[int, int, int]] = {}
    for layer in sorted(set(variable_layers)):
        variable_count = variable_layers.count(layer)
        layer_rows = rows_by_layer[layer]
        layer_rank = rank_affine_rows(layer_rows)
        layer_statistics[layer] = (
            len(layer_rows),
            layer_rank,
            variable_count - layer_rank,
        )
    assert sum(statistics[1] for statistics in layer_statistics.values()) == rank
    components = affine_components(rows, variable_layers)
    return (
        len(rows),
        rank,
        len(variables_to_audit) - rank,
        layer_statistics,
        affine_component_statistics(components),
        components,
    )


def variables(prefix: str, basis: list[Monomial]) -> dict[Monomial, z3.BoolRef]:
    return {
        monomial: z3.Bool(f"{prefix}_{monomial[0]}_{monomial[1]}")
        for monomial in basis
    }


def second_binary_digit(
    terms: list[tuple[int, z3.BoolRef | None]],
) -> z3.BoolRef:
    """Return whether an integer linear form is 2 modulo 4."""

    total = z3.BitVecVal(0, 2)
    for coefficient, variable in terms:
        value = z3.BitVecVal(coefficient % 4, 2)
        total += value if variable is None else z3.If(
            variable, value, z3.BitVecVal(0, 2)
        )
    return total == z3.BitVecVal(2, 2)


def add_term(
    polynomial: Polynomial3,
    monomial: Monomial3,
    coefficient: int,
    modulus: int = 8,
) -> None:
    coefficient = (polynomial.get(monomial, 0) + coefficient) % modulus
    if coefficient:
        polynomial[monomial] = coefficient
    else:
        polynomial.pop(monomial, None)


def derivative(
    polynomial: Polynomial3,
    axis: int,
    modulus: int = 8,
) -> Polynomial3:
    output: Polynomial3 = {}
    for monomial, coefficient in polynomial.items():
        exponent = monomial[axis]
        if not exponent:
            continue
        result = list(monomial)
        result[axis] -= 1
        add_term(output, tuple(result), coefficient * exponent, modulus)
    return output


def product(
    left: Polynomial3,
    right: Polynomial3,
    modulus: int = 8,
) -> Polynomial3:
    output: Polynomial3 = {}
    for (i, j, k), left_coefficient in left.items():
        for (u, v, w), right_coefficient in right.items():
            add_term(
                output,
                (i + u, j + v, k + w),
                left_coefficient * right_coefficient,
                modulus,
            )
    return output


def product3(
    first: Polynomial3,
    second: Polynomial3,
    third: Polynomial3,
    modulus: int = 8,
) -> Polynomial3:
    return product(product(first, second, modulus), third, modulus)


def add_symbolic_term(
    polynomial: SymbolicPolynomial3,
    monomial: Monomial3,
    coefficient: int,
    variable: z3.BoolRef | None,
) -> None:
    if coefficient:
        polynomial.setdefault(monomial, []).append((coefficient, variable))


def symbolic_derivative(
    polynomial: SymbolicPolynomial3,
    axis: int,
) -> SymbolicPolynomial3:
    output: SymbolicPolynomial3 = {}
    for monomial, terms in polynomial.items():
        exponent = monomial[axis]
        if not exponent:
            continue
        result = list(monomial)
        result[axis] -= 1
        for coefficient, variable in terms:
            add_symbolic_term(
                output,
                tuple(result),
                coefficient * exponent,
                variable,
            )
    return output


def build_symbolic_witt_jacobian(
    first: tuple[
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
    ],
    second: tuple[
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
    ],
) -> list[list[list[SymbolicPolynomial3]]]:
    p, q, c, a, b, r = first
    u0, v0, w1, u1, v1, u2, v2, w3 = second
    coordinate_digits: list[list[SymbolicPolynomial3]] = [
        [{}, {}, {}],
        [{}, {}, {}],
        [{}, {}, {}],
    ]
    for coordinate, base in enumerate(
        (
            {
                (1, 0, 0): 1,
                (2, 1, 0): 1,
                (4, 0, 0): 1,
                (6, 2, 0): 1,
            },
            {
                (0, 1, 0): 1,
                (5, 0, 0): 1,
                (6, 1, 0): 1,
                (7, 2, 0): 1,
                (8, 3, 0): 1,
            },
            {(0, 0, 1): 1},
        )
    ):
        for monomial, coefficient in base.items():
            add_symbolic_term(
                coordinate_digits[0][coordinate],
                monomial,
                coefficient,
                None,
            )

    for polynomial, coordinate, z_degree in (
        (p, 0, 0),
        (q, 1, 0),
        (c, 2, 0),
        (a, 0, 1),
        (b, 1, 1),
        (r, 2, 1),
    ):
        for (i, j), variable in polynomial.items():
            add_symbolic_term(
                coordinate_digits[1][coordinate],
                (i, j, z_degree),
                1,
                variable,
            )
    for polynomial, coordinate, z_degree in (
        (u0, 0, 0),
        (v0, 1, 0),
        (w1, 2, 1),
        (u1, 0, 1),
        (v1, 1, 1),
        (u2, 0, 2),
        (v2, 1, 2),
        (w3, 2, 3),
    ):
        for (i, j), variable in polynomial.items():
            add_symbolic_term(
                coordinate_digits[2][coordinate],
                (i, j, z_degree),
                1,
                variable,
            )

    return [
        [
            [symbolic_derivative(coordinate_digits[digit][row], axis)
             for axis in range(3)]
            for row in range(3)
        ]
        for digit in range(3)
    ]


def compile_w4_top_class_coefficients(
    first: tuple[
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
    ],
    second: tuple[
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
    ],
) -> tuple[dict[Monomial3, z3.BitVecRef], int]:
    """Compile odd-odd-odd coefficients of the determinant modulo 16.

    Write the Jacobian matrix as M0+2*M1+4*M2.  On an odd-z monomial the
    base determinant is zero.  Terms of 2-adic weights one and two determine
    the top de Rham class of the next digit.  Weight-three terms are sums of
    wedges of exact one-forms (the M1-M2 cross terms and det(M1)), so their
    top de Rham projection vanishes identically and need not be expanded.
    """

    jacobian = build_symbolic_witt_jacobian(first, second)
    raw: dict[
        Monomial3,
        dict[tuple[int, ...], tuple[int, tuple[z3.BoolRef, ...]]],
    ] = defaultdict(dict)

    def record(
        monomial: Monomial3,
        coefficient: int,
        variables_in_term: tuple[z3.BoolRef, ...],
    ) -> None:
        coefficient %= 16
        if not coefficient or not all(exponent % 2 for exponent in monomial):
            return
        by_id = {variable.get_id(): variable for variable in variables_in_term}
        key = tuple(sorted(by_id))
        variables_key = tuple(by_id[index] for index in key)
        old_coefficient = raw[monomial].get(key, (0, variables_key))[0]
        new_coefficient = (old_coefficient + coefficient) % 16
        if new_coefficient:
            raw[monomial][key] = (new_coefficient, variables_key)
        else:
            raw[monomial].pop(key, None)

    signed_permutations = (
        (1, (0, 1, 2)),
        (1, (1, 2, 0)),
        (1, (2, 0, 1)),
        (-1, (2, 1, 0)),
        (-1, (1, 0, 2)),
        (-1, (0, 2, 1)),
    )
    digit_choices = [
        choice
        for choice in itertools.product(range(3), repeat=3)
        if sum(choice) in (1, 2)
    ]
    for sign, columns in signed_permutations:
        for digits in digit_choices:
            digit_sum = sum(digits)
            entries = [
                jacobian[digits[row]][row][columns[row]]
                for row in range(3)
            ]
            if any(not entry for entry in entries):
                continue
            entries.sort(key=len)
            for monomial0, terms0 in entries[0].items():
                for monomial1, terms1 in entries[1].items():
                    partial = tuple(
                        monomial0[index] + monomial1[index]
                        for index in range(3)
                    )
                    for monomial2, terms2 in entries[2].items():
                        target = tuple(
                            partial[index] + monomial2[index]
                            for index in range(3)
                        )
                        if not all(exponent % 2 for exponent in target):
                            continue
                        for coefficient0, variable0 in terms0:
                            for coefficient1, variable1 in terms1:
                                for coefficient2, variable2 in terms2:
                                    variables_in_term = tuple(
                                        variable
                                        for variable in (
                                            variable0,
                                            variable1,
                                            variable2,
                                        )
                                        if variable is not None
                                    )
                                    record(
                                        target,
                                        sign
                                        * (1 << digit_sum)
                                        * coefficient0
                                        * coefficient1
                                        * coefficient2,
                                        variables_in_term,
                                    )

    expressions: dict[Monomial3, z3.BitVecRef] = {}
    term_count = 0
    for monomial, terms in raw.items():
        expression = z3.BitVecVal(0, 4)
        for coefficient, variables_in_term in terms.values():
            condition = (
                z3.BoolVal(True)
                if not variables_in_term
                else z3.And(*variables_in_term)
            )
            expression += z3.If(
                condition,
                z3.BitVecVal(coefficient, 4),
                z3.BitVecVal(0, 4),
            )
            term_count += 1
        expression = z3.simplify(expression)
        if not (z3.is_bv_value(expression) and expression.as_long() == 0):
            expressions[monomial] = expression
    return dict(sorted(expressions.items())), term_count


def compile_w4_target_coefficients(
    first: tuple[
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
    ],
    second: tuple[
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
    ],
    targets: set[Monomial3],
) -> tuple[dict[Monomial3, z3.BitVecRef], int]:
    """Compile selected determinant coefficients modulo sixteen exactly."""

    jacobian = build_symbolic_witt_jacobian(first, second)
    raw: dict[
        Monomial3,
        dict[tuple[int, ...], tuple[int, tuple[z3.BoolRef, ...]]],
    ] = {target: {} for target in targets}

    def record(
        target: Monomial3,
        coefficient: int,
        variables_in_term: tuple[z3.BoolRef, ...],
    ) -> None:
        coefficient %= 16
        if not coefficient:
            return
        by_id = {variable.get_id(): variable for variable in variables_in_term}
        key = tuple(sorted(by_id))
        variables_key = tuple(by_id[index] for index in key)
        old_coefficient = raw[target].get(key, (0, variables_key))[0]
        new_coefficient = (old_coefficient + coefficient) % 16
        if new_coefficient:
            raw[target][key] = (new_coefficient, variables_key)
        else:
            raw[target].pop(key, None)

    signed_permutations = (
        (1, (0, 1, 2)),
        (1, (1, 2, 0)),
        (1, (2, 0, 1)),
        (-1, (2, 1, 0)),
        (-1, (1, 0, 2)),
        (-1, (0, 2, 1)),
    )
    digit_choices = [
        choice
        for choice in itertools.product(range(3), repeat=3)
        if sum(choice) <= 3
    ]
    for sign, columns in signed_permutations:
        for digits in digit_choices:
            digit_sum = sum(digits)
            entries = [
                jacobian[digits[row]][row][columns[row]]
                for row in range(3)
            ]
            if any(not entry for entry in entries):
                continue
            entries.sort(key=len)
            for monomial0, terms0 in entries[0].items():
                for monomial1, terms1 in entries[1].items():
                    partial = tuple(
                        monomial0[index] + monomial1[index]
                        for index in range(3)
                    )
                    matches: list[
                        tuple[Monomial3, list[tuple[int, z3.BoolRef | None]]]
                    ] = []
                    if len(targets) <= len(entries[2]):
                        for target in targets:
                            needed = tuple(
                                target[index] - partial[index]
                                for index in range(3)
                            )
                            if min(needed) < 0:
                                continue
                            terms2 = entries[2].get(needed)
                            if terms2 is not None:
                                matches.append((target, terms2))
                    else:
                        for monomial2, terms2 in entries[2].items():
                            target = tuple(
                                partial[index] + monomial2[index]
                                for index in range(3)
                            )
                            if target in targets:
                                matches.append((target, terms2))
                    for target, terms2 in matches:
                        for coefficient0, variable0 in terms0:
                            for coefficient1, variable1 in terms1:
                                for coefficient2, variable2 in terms2:
                                    variables_in_term = tuple(
                                        variable
                                        for variable in (
                                            variable0,
                                            variable1,
                                            variable2,
                                        )
                                        if variable is not None
                                    )
                                    record(
                                        target,
                                        sign
                                        * (1 << digit_sum)
                                        * coefficient0
                                        * coefficient1
                                        * coefficient2,
                                        variables_in_term,
                                    )

    expressions: dict[Monomial3, z3.BitVecRef] = {}
    term_count = 0
    for target in sorted(targets):
        expression = z3.BitVecVal(0, 4)
        for coefficient, variables_in_term in raw[target].values():
            condition = (
                z3.BoolVal(True)
                if not variables_in_term
                else z3.And(*variables_in_term)
            )
            expression += z3.If(
                condition,
                z3.BitVecVal(coefficient, 4),
                z3.BitVecVal(0, 4),
            )
            term_count += 1
        expressions[target] = z3.simplify(expression)
    return expressions, term_count


def build_weighted_symbolic_witt_jacobian(
    first: tuple[
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
    ],
    second: tuple[
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
    ],
) -> list[list[BitVectorPolynomial3]]:
    """Combine M0+2*M1+4*M2 into a modulo-sixteen symbolic Jacobian."""

    digit_jacobian = build_symbolic_witt_jacobian(first, second)
    weighted: list[list[BitVectorPolynomial3]] = [
        [{}, {}, {}],
        [{}, {}, {}],
        [{}, {}, {}],
    ]
    for row in range(3):
        for axis in range(3):
            monomials_in_entry = set().union(
                *(set(digit_jacobian[digit][row][axis]) for digit in range(3))
            )
            for monomial in monomials_in_entry:
                expression = z3.BitVecVal(0, 4)
                for digit in range(3):
                    for coefficient, variable in digit_jacobian[digit][row][axis].get(
                        monomial,
                        [],
                    ):
                        value = z3.BitVecVal(
                            (1 << digit) * coefficient,
                            4,
                        )
                        expression += (
                            value
                            if variable is None
                            else z3.If(
                                variable,
                                value,
                                z3.BitVecVal(0, 4),
                            )
                        )
                expression = z3.simplify(expression)
                if not (
                    z3.is_bv_value(expression) and expression.as_long() == 0
                ):
                    weighted[row][axis][monomial] = expression
    return weighted


def compile_w4_factored_target_coefficients(
    first: tuple[
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
    ],
    second: tuple[
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
    ],
    targets: set[Monomial3],
) -> tuple[dict[Monomial3, z3.BitVecRef], dict[str, int]]:
    """Compile target coefficients through shared projected two-by-two minors.

    For each determinant permutation, choose the entry to contract last.
    Only pair coefficients that can reach a requested quotient target are
    materialized.  Their BitVec expressions are shared by every final target,
    avoiding the legacy expansion into separate cubic Boolean conjunctions.
    """

    jacobian = build_weighted_symbolic_witt_jacobian(first, second)
    expressions = {
        target: z3.BitVecVal(0, 4)
        for target in targets
    }
    statistics = {
        "minor_coefficients": 0,
        "minor_products": 0,
        "final_products": 0,
    }
    for sign, columns in (
        (1, (0, 1, 2)),
        (1, (1, 2, 0)),
        (1, (2, 0, 1)),
        (-1, (2, 1, 0)),
        (-1, (1, 0, 2)),
        (-1, (0, 2, 1)),
    ):
        entries = [jacobian[row][columns[row]] for row in range(3)]
        if any(not entry for entry in entries):
            continue

        factorings: list[
            tuple[
                int,
                int,
                int,
                set[Monomial3],
            ]
        ] = []
        for last_index in range(3):
            pair_indices = [index for index in range(3) if index != last_index]
            last_entry = entries[last_index]
            needed_pairs: set[Monomial3] = set()
            final_lookup_count = 0
            for target in targets:
                for last_monomial in last_entry:
                    needed = tuple(
                        target[index] - last_monomial[index]
                        for index in range(3)
                    )
                    if min(needed) >= 0:
                        needed_pairs.add(needed)
                        final_lookup_count += 1
            pair_lookup_count = (
                len(needed_pairs)
                * min(
                    len(entries[pair_indices[0]]),
                    len(entries[pair_indices[1]]),
                )
            )
            factorings.append(
                (
                    pair_lookup_count + final_lookup_count,
                    pair_indices[0],
                    pair_indices[1],
                    needed_pairs,
                )
            )
        _, left_index, right_index, needed_pairs = min(
            factorings,
            key=lambda factoring: factoring[0],
        )
        last_index = ({0, 1, 2} - {left_index, right_index}).pop()
        left_entry = entries[left_index]
        right_entry = entries[right_index]
        last_entry = entries[last_index]
        if len(left_entry) > len(right_entry):
            left_entry, right_entry = right_entry, left_entry

        pair_coefficients: BitVectorPolynomial3 = {}
        for pair_target in needed_pairs:
            products: list[z3.BitVecRef] = []
            for left_monomial, left_expression in left_entry.items():
                right_monomial = tuple(
                    pair_target[index] - left_monomial[index]
                    for index in range(3)
                )
                if min(right_monomial) < 0:
                    continue
                right_expression = right_entry.get(right_monomial)
                if right_expression is not None:
                    products.append(left_expression * right_expression)
            if products:
                pair_coefficients[pair_target] = sum(
                    products,
                    z3.BitVecVal(0, 4),
                )
                statistics["minor_coefficients"] += 1
                statistics["minor_products"] += len(products)

        for target in targets:
            products = []
            for last_monomial, last_expression in last_entry.items():
                pair_target = tuple(
                    target[index] - last_monomial[index]
                    for index in range(3)
                )
                if min(pair_target) < 0:
                    continue
                pair_expression = pair_coefficients.get(pair_target)
                if pair_expression is not None:
                    products.append(pair_expression * last_expression)
            if products:
                contribution = sum(products, z3.BitVecVal(0, 4))
                expressions[target] += (
                    contribution if sign == 1 else -contribution
                )
                statistics["final_products"] += len(products)

    return {
        target: z3.simplify(expression)
        for target, expression in expressions.items()
    }, statistics


def structural_w4_target_universe(
    first: tuple[
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
    ],
    second: tuple[
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
        dict[Monomial, z3.BoolRef],
    ],
) -> set[Monomial3]:
    """Return every structurally possible determinant target modulo 16."""

    jacobian = build_symbolic_witt_jacobian(first, second)

    def minkowski(
        left: set[Monomial3],
        right: set[Monomial3],
    ) -> set[Monomial3]:
        return {
            tuple(left_monomial[index] + right_monomial[index] for index in range(3))
            for left_monomial in left
            for right_monomial in right
        }

    universe: set[Monomial3] = set()
    for columns in (
        (0, 1, 2),
        (1, 2, 0),
        (2, 0, 1),
        (2, 1, 0),
        (1, 0, 2),
        (0, 2, 1),
    ):
        for digits in itertools.product(range(3), repeat=3):
            if sum(digits) > 3:
                continue
            entries = [
                set(jacobian[digits[row]][row][columns[row]])
                for row in range(3)
            ]
            if any(not entry for entry in entries):
                continue
            universe.update(minkowski(minkowski(entries[0], entries[1]), entries[2]))
    return universe


SUPPORT_LAYOUT = {
    "p": (2, 0, 0),
    "q": (2, 1, 0),
    "c": (2, 2, 0),
    "a": (2, 0, 1),
    "b": (2, 1, 1),
    "r": (2, 2, 1),
    "u0": (4, 0, 0),
    "v0": (4, 1, 0),
    "w1": (4, 2, 1),
    "u1": (4, 0, 1),
    "v1": (4, 1, 1),
    "u2": (4, 0, 2),
    "v2": (4, 1, 2),
    "w3": (4, 2, 3),
}


def supports_from_model(
    model: z3.ModelRef,
    named_variables: list[tuple[str, dict[Monomial, z3.BoolRef]]],
) -> dict[str, list[Monomial]]:
    return {
        name: [
            monomial
            for monomial, variable in polynomial.items()
            if z3.is_true(model.eval(variable, model_completion=True))
        ]
        for name, polynomial in named_variables
    }


def verify_direct_supports(
    supports: dict[str, list[Monomial]],
    maximum_degree: int,
) -> int:
    """Reconstruct support data and recompute its determinant modulo eight."""

    assert set(supports) == set(SUPPORT_LAYOUT)

    coordinates: list[Polynomial3] = [
        {
            (1, 0, 0): 1,
            (2, 1, 0): 1,
            (4, 0, 0): 1,
            (6, 2, 0): 1,
        },
        {
            (0, 1, 0): 1,
            (5, 0, 0): 1,
            (6, 1, 0): 1,
            (7, 2, 0): 1,
            (8, 3, 0): 1,
        },
        {(0, 0, 1): 1},
    ]

    for name, support in supports.items():
        scale, coordinate, z_degree = SUPPORT_LAYOUT[name]
        assert len(support) == len(set(support))
        for i, j in support:
            assert i >= 0 and j >= 0
            add_term(coordinates[coordinate], (i, j, z_degree), scale)

    assert all(
        max(map(sum, coordinate), default=0) <= maximum_degree
        for coordinate in coordinates
    )
    reductions = [
        {monomial for monomial, coefficient in coordinate.items() if coefficient % 2}
        for coordinate in coordinates
    ]
    assert reductions == [
        {(1, 0, 0), (2, 1, 0), (4, 0, 0), (6, 2, 0)},
        {(0, 1, 0), (5, 0, 0), (6, 1, 0), (7, 2, 0), (8, 3, 0)},
        {(0, 0, 1)},
    ]

    jacobian = [
        [derivative(coordinate, axis) for axis in range(3)]
        for coordinate in coordinates
    ]
    determinant: Polynomial3 = {}
    for sign, indices in (
        (1, (0, 1, 2)),
        (1, (1, 2, 0)),
        (1, (2, 0, 1)),
        (-1, (2, 1, 0)),
        (-1, (1, 0, 2)),
        (-1, (0, 2, 1)),
    ):
        term = product3(
            jacobian[0][indices[0]],
            jacobian[1][indices[1]],
            jacobian[2][indices[2]],
        )
        for monomial, coefficient in term.items():
            add_term(determinant, monomial, sign * coefficient)

    assert set(determinant) == {(0, 0, 0)}, determinant
    constant = determinant[(0, 0, 0)]
    assert constant % 2 == 1
    return constant


def coordinates_from_supports(
    supports: dict[str, list[Monomial]],
    modulus: int,
) -> list[Polynomial3]:
    """Reconstruct the three stabilized coordinates at a chosen modulus."""

    coordinates: list[Polynomial3] = [
        {
            (1, 0, 0): 1,
            (2, 1, 0): 1,
            (4, 0, 0): 1,
            (6, 2, 0): 1,
        },
        {
            (0, 1, 0): 1,
            (5, 0, 0): 1,
            (6, 1, 0): 1,
            (7, 2, 0): 1,
            (8, 3, 0): 1,
        },
        {(0, 0, 1): 1},
    ]
    assert set(supports) == set(SUPPORT_LAYOUT)
    for name, support in supports.items():
        scale, coordinate, z_degree = SUPPORT_LAYOUT[name]
        for i, j in support:
            add_term(
                coordinates[coordinate],
                (i, j, z_degree),
                scale,
                modulus,
            )
    return coordinates


def determinant_polynomial(
    coordinates: list[Polynomial3],
    modulus: int,
) -> Polynomial3:
    jacobian = [
        [derivative(coordinate, axis, modulus) for axis in range(3)]
        for coordinate in coordinates
    ]
    determinant: Polynomial3 = {}
    for sign, indices in (
        (1, (0, 1, 2)),
        (1, (1, 2, 0)),
        (1, (2, 0, 1)),
        (-1, (2, 1, 0)),
        (-1, (1, 0, 2)),
        (-1, (0, 2, 1)),
    ):
        term = product3(
            jacobian[0][indices[0]],
            jacobian[1][indices[1]],
            jacobian[2][indices[2]],
            modulus,
        )
        for monomial, coefficient in term.items():
            add_term(determinant, monomial, sign * coefficient, modulus)
    return determinant


def determinant_error_digit(
    supports: dict[str, list[Monomial]],
    constant_mod_8: int,
) -> Polynomial3:
    """Return (det-c)/8 modulo two, asserting det=c modulo eight."""

    determinant = determinant_polynomial(
        coordinates_from_supports(supports, 16),
        16,
    )
    error: Polynomial3 = {}
    for monomial in set(determinant) | {(0, 0, 0)}:
        residue = (
            determinant.get(monomial, 0)
            - constant_mod_8 * int(monomial == (0, 0, 0))
        ) % 16
        assert residue in (0, 8), (monomial, residue)
        if residue == 8:
            error[monomial] = 1
    return error


def full_w4_completion_rows(
    error: Polynomial3,
    degree: int,
) -> tuple[
    list[tuple[int, int]],
    list[tuple[str, Monomial3]],
    list[Monomial3],
]:
    """Encode D_F(R,S,T)=error for corrections of total degree <= degree."""

    contributions: dict[Monomial3, int] = {}
    variable_keys: list[tuple[str, Monomial3]] = []

    def add_variable(
        name: str,
        source: Monomial3,
        targets: set[Monomial3],
    ) -> None:
        if not targets:
            return
        index = len(variable_keys)
        variable_keys.append((name, source))
        for target in targets:
            contributions[target] = contributions.get(target, 0) ^ (1 << index)

    for z_degree in sorted({monomial[2] for monomial in error}):
        plane_bound = degree - z_degree
        if plane_bound >= 0:
            for i, j in monomials(plane_bound):
                r_targets: set[Monomial3] = set()
                if i % 2:
                    r_targets.update(
                        (i - 1 + u, j + v, z_degree) for u, v in QY
                    )
                if j % 2:
                    r_targets.update(
                        (i + u, j - 1 + v, z_degree) for u, v in QX
                    )
                add_variable("R", (i, j, z_degree), r_targets)

                s_targets: set[Monomial3] = set()
                if j % 2:
                    s_targets.update(
                        (i + u, j - 1 + v, z_degree) for u, v in PX
                    )
                if i % 2:
                    s_targets.update(
                        (i - 1 + u, j + v, z_degree) for u, v in PY
                    )
                add_variable("S", (i, j, z_degree), s_targets)

        t_plane_bound = degree - z_degree - 1
        if z_degree % 2 == 0 and t_plane_bound >= 0:
            for i, j in monomials(t_plane_bound):
                add_variable(
                    "T",
                    (i, j, z_degree + 1),
                    {(i, j, z_degree)},
                )

    targets = sorted(set(contributions) | set(error))
    nonconstant_targets = [
        target for target in targets if target != (0, 0, 0)
    ]
    rows = [
        (contributions.get(target, 0), int(target in error))
        for target in nonconstant_targets
    ]
    return rows, variable_keys, nonconstant_targets


def load_certificate(path: Path) -> tuple[dict, dict[str, list[Monomial]]]:
    raw = path.read_bytes()
    if path.name in FROZEN_CERTIFICATE_SHA256:
        assert (
            hashlib.sha256(raw).hexdigest()
            == FROZEN_CERTIFICATE_SHA256[path.name]
        )
    data = json.loads(raw)
    assert data["schema_version"] == 1
    assert data["claim"] == "HKM2W3"
    assert data["degree"] == 19
    assert data["modulus"] == 8
    assert data["system"] == {"constraints": 4513, "variables": 2685}
    supports = {
        name: [tuple(monomial) for monomial in support]
        for name, support in data["supports"].items()
    }
    counts = {name: len(support) for name, support in supports.items()}
    assert counts == data["support_counts"]
    assert sum(counts.values()) == data["total_support"]
    return data, supports


def replay_certificate(path: Path) -> None:
    data, supports = load_certificate(path)
    if path.name == "huq_kuruvilla_w3_degree19_witness_second_180.json":
        assert set(supports["u2"]) == {(5, 2), (6, 3)}
        assert set(supports["v2"]) == {
            (4, 5),
            (7, 2),
            (7, 8),
            (8, 9),
        }
        assert set(supports["w3"]) == {(0, 0), (0, 4), (8, 4)}
    if path.name == "huq_kuruvilla_w3_degree19_witness_second_172.json":
        assert data["total_support"] == 471
        assert data["jacobian_constant_mod_8"] == 1
        assert data["support_counts"] == {
            "a": 17,
            "b": 38,
            "c": 61,
            "p": 51,
            "q": 80,
            "r": 52,
            "u0": 50,
            "u1": 30,
            "u2": 0,
            "v0": 24,
            "v1": 18,
            "v2": 3,
            "w1": 44,
            "w3": 3,
        }
        assert supports["u2"] == []
        assert set(supports["v2"]) == {(3, 2), (4, 5), (8, 3)}
        assert set(supports["w3"]) == {(0, 0), (0, 4), (8, 4)}
    if path.name == "huq_kuruvilla_w3_degree19_witness_second_165.json":
        assert data["total_support"] == 464
        assert data["jacobian_constant_mod_8"] == 1
        assert data["optimization"]["proved_minimum"] is True
        assert data["optimization"]["second_support"] == 165
        assert data["optimization"]["second_support_lower_bound"] == 165
        assert data["support_counts"] == {
            "a": 17,
            "b": 38,
            "c": 61,
            "p": 51,
            "q": 80,
            "r": 52,
            "u0": 46,
            "u1": 30,
            "u2": 0,
            "v0": 25,
            "v1": 18,
            "v2": 3,
            "w1": 40,
            "w3": 3,
        }
        assert supports["u2"] == []
        assert set(supports["v2"]) == {(3, 2), (4, 5), (8, 3)}
        assert set(supports["w3"]) == {(0, 0), (0, 4), (8, 4)}
    if path.name == (
        "huq_kuruvilla_w3_degree19_witness_first_280_second_160.json"
    ):
        assert data["total_support"] == 440
        assert data["jacobian_constant_mod_8"] == 1
        assert data["optimization"]["proved_minimum"] is True
        assert data["optimization"]["second_support"] == 160
        assert data["optimization"]["second_support_lower_bound"] == 160
        assert data["support_counts"] == {
            "a": 27,
            "b": 61,
            "c": 52,
            "p": 41,
            "q": 54,
            "r": 45,
            "u0": 41,
            "u1": 28,
            "u2": 2,
            "v0": 12,
            "v1": 26,
            "v2": 5,
            "w1": 39,
            "w3": 7,
        }
        assert set(supports["u2"]) == {(1, 0), (2, 5)}
        assert set(supports["v2"]) == {
            (3, 0),
            (7, 0),
            (3, 6),
            (11, 0),
            (3, 14),
        }
        assert set(supports["w3"]) == {
            (0, 2),
            (2, 4),
            (6, 0),
            (2, 8),
            (2, 10),
            (0, 14),
            (2, 12),
        }
    if path.name == "huq_kuruvilla_w3_degree19_w4_class_zero.json":
        assert data["total_support"] == 818
        assert data["jacobian_constant_mod_8"] == 5
        assert data["w4_top_class_zero"] == {
            "active_constraints": 4754,
            "collected_symbolic_terms": 38760,
            "equations": 241,
            "proved": True,
        }
        assert data["support_counts"] == {
            "a": 22,
            "b": 53,
            "c": 51,
            "p": 49,
            "q": 97,
            "r": 55,
            "u0": 71,
            "u1": 64,
            "u2": 47,
            "v0": 71,
            "v1": 53,
            "v2": 40,
            "w1": 89,
            "w3": 56,
        }
    constant = verify_direct_supports(supports, data["degree"])
    assert constant == data["jacobian_constant_mod_8"]
    assert constant * constant % 8 == 1
    print(
        "PASS: frozen degree-19 witness has "
        f"{data['total_support']} correction coefficients"
    )
    print(f"PASS: direct certificate replay gives Jacobian {constant} modulo 8")
    print("PASS: target scaling gives determinant one at unchanged degree")
    if path.name == "huq_kuruvilla_w3_degree19_witness_second_180.json":
        print("PASS: compact u2, v2, w3 tail support identities")
    if path.name == "huq_kuruvilla_w3_degree19_witness_second_172.json":
        print("PASS: six-term z2 tail support identities")
    if path.name == "huq_kuruvilla_w3_degree19_witness_second_165.json":
        print("PASS: component-minimized support identities")
    if path.name == (
        "huq_kuruvilla_w3_degree19_witness_first_280_second_160.json"
    ):
        print("PASS: preferred sparse-first component-minimized identities")
    if path.name == "huq_kuruvilla_w3_degree19_w4_class_zero.json":
        print("PASS: pinned W_4-class-zero degree-19 support identities")


def add_plane_variation(
    output: BooleanPolynomial,
    variable: z3.BoolRef,
    i: int,
    j: int,
    coordinate: str,
) -> None:
    """Add D_F(variable*x^i*y^j,0) or D_F(0,...)."""

    if coordinate == "A":
        if i % 2:
            for u, v in QY:
                output[(i - 1 + u, j + v)].append(variable)
        if j % 2:
            for u, v in QX:
                output[(i + u, j - 1 + v)].append(variable)
    else:
        if j % 2:
            output[(i, j - 1)].append(variable)
        if i % 2:
            output[(i + 1, j)].append(variable)


def add_integer_derivative_product(
    output: dict[Monomial, list[tuple[int, z3.BoolRef]]],
    polynomial: dict[Monomial, z3.BoolRef],
    axis: int,
    multiplier: dict[Monomial, int],
    sign: int,
) -> None:
    for (i, j), variable in polynomial.items():
        derivative_coefficient = i if axis == 0 else j
        if not derivative_coefficient:
            continue
        exponent = (i - 1, j) if axis == 0 else (i, j - 1)
        for (u, v), coefficient in multiplier.items():
            output[(exponent[0] + u, exponent[1] + v)].append(
                (sign * derivative_coefficient * coefficient, variable)
            )


def add_integer_product(
    output: dict[Monomial, list[tuple[int, z3.BoolRef]]],
    polynomial: dict[Monomial, z3.BoolRef],
    multiplier: dict[Monomial, int],
) -> None:
    for (i, j), variable in polynomial.items():
        for (u, v), coefficient in multiplier.items():
            output[(i + u, j + v)].append((coefficient, variable))


if args.replay_certificate is not None:
    replay_certificate(args.replay_certificate)
    raise SystemExit


degree_main = monomials(args.degree)
degree_odd = monomials(args.degree - 1)

# Write the first correction as
#   A=p+z*a+..., B=q+z*b+..., C=c+z*r+...
# Higher z-layers cannot affect the z^0 and z^1 constraints used below.
p = variables("p", degree_main)
q = variables("q", degree_main)
r = variables("r", degree_odd)
a = variables("a", degree_odd)
b = variables("b", degree_odd)
c = variables("c", degree_main)

solver_tactics = [z3.Tactic("simplify")]
if args.solve_equations_first:
    solver_tactics.append(z3.Tactic("solve-eqs"))
solver_tactics.extend((z3.Tactic("bit-blast"), z3.Tactic("sat")))
solver = z3.Then(*solver_tactics).solver()
solver.set(timeout=args.timeout_ms, random_seed=args.random_seed)

# First Witt digit.  The z^0 equation is K+D_F(p,q)+r=constant, while the
# odd z-layer has D_F(a,b)=0 because d(z^2)/dz vanishes in characteristic two.
variation_zero: BooleanPolynomial = defaultdict(list)
variation_one: BooleanPolynomial = defaultdict(list)
for (i, j), variable in p.items():
    add_plane_variation(variation_zero, variable, i, j, "A")
for (i, j), variable in q.items():
    add_plane_variation(variation_zero, variable, i, j, "B")
for (i, j), variable in a.items():
    add_plane_variation(variation_one, variable, i, j, "A")
for (i, j), variable in b.items():
    add_plane_variation(variation_one, variable, i, j, "B")

for monomial in set(variation_zero) | set(r) | K_SUPPORT:
    if monomial == (0, 0):
        continue
    terms = list(variation_zero[monomial])
    if monomial in r:
        terms.append(r[monomial])
    if monomial in K_SUPPORT:
        terms.append(z3.BoolVal(True))
    solver.add(xor(terms) == z3.BoolVal(False))
for terms in variation_one.values():
    solver.add(xor(terms) == z3.BoolVal(False))

# The first-Witt equations themselves are consistent (the canonical
# correction is one witness); the contradiction must come from the next
# determinant digit.
assert solver.check() == z3.sat

# Integral z^0 first variation.  On the first-digit solution space its
# nonconstant coefficients are even; their second binary digit is the
# Bockstein contribution to the W_3 error.
linear_zero: dict[Monomial, list[tuple[int, z3.BoolRef | None]]] = defaultdict(list)
for monomial, coefficient in K_INTEGER.items():
    linear_zero[monomial].append((coefficient, None))
add_integer_derivative_product(linear_zero, p, 0, QY_INTEGER, 1)
add_integer_derivative_product(linear_zero, q, 1, PX_INTEGER, 1)
add_integer_derivative_product(linear_zero, p, 1, QX_INTEGER, -1)
add_integer_derivative_product(linear_zero, q, 0, PY_INTEGER, -1)
add_integer_product(linear_zero, r, J_INTEGER)

# Quadratic z^0 determinant digit.
quadratic_zero: BooleanPolynomial = defaultdict(list)
for (i, j), p_variable in p.items():
    for (u, v), q_variable in q.items():
        if i % 2 and v % 2:
            quadratic_zero[(i - 1 + u, j + v - 1)].append(
                z3.And(p_variable, q_variable)
            )
        if j % 2 and u % 2:
            quadratic_zero[(i + u - 1, j - 1 + v)].append(
                z3.And(p_variable, q_variable)
            )

# The part D_F(p,q)*r.
for (i, j), terms in variation_zero.items():
    coefficient = xor(terms)
    for (u, v), r_variable in r.items():
        quadratic_zero[(i + u, j + v)].append(z3.And(coefficient, r_variable))

# The cross term from z*a,z*b and the z-independent third correction c.
for (ci, cj), c_variable in c.items():
    if ci % 2:
        for (i, j), a_variable in a.items():
            for u, v in QY:
                quadratic_zero[(i + ci - 1 + u, j + cj + v)].append(
                    z3.And(a_variable, c_variable)
                )
        for (i, j), b_variable in b.items():
            for u, v in PY:
                quadratic_zero[(i + ci - 1 + u, j + cj + v)].append(
                    z3.And(b_variable, c_variable)
                )
    if cj % 2:
        for (i, j), a_variable in a.items():
            for u, v in QX:
                quadratic_zero[(i + ci + u, j + cj - 1 + v)].append(
                    z3.And(a_variable, c_variable)
                )
        for (i, j), b_variable in b.items():
            for u, v in PX:
                quadratic_zero[(i + ci + u, j + cj - 1 + v)].append(
                    z3.And(b_variable, c_variable)
                )

# A degree-d second correction reaches plane degree at most d+9 in z^0.  Every
# higher coefficient of the existing second digit must therefore vanish.
for monomial in set(linear_zero) | set(quadratic_zero):
    if sum(monomial) <= args.degree + 9:
        continue
    bockstein = second_binary_digit(linear_zero[monomial])
    solver.add(xor([bockstein] + quadratic_zero[monomial]) == z3.BoolVal(False))

# Integral and quadratic z^1 digits.  Corrections in this layer have plane
# coefficient degree at most d-1 and hence reach plane degree at most d+8.
linear_one: dict[Monomial, list[tuple[int, z3.BoolRef]]] = defaultdict(list)
add_integer_derivative_product(linear_one, a, 0, QY_INTEGER, 1)
add_integer_derivative_product(linear_one, b, 1, PX_INTEGER, 1)
add_integer_derivative_product(linear_one, a, 1, QX_INTEGER, -1)
add_integer_derivative_product(linear_one, b, 0, PY_INTEGER, -1)

quadratic_one: BooleanPolynomial = defaultdict(list)
for (i, j), p_variable in p.items():
    for (u, v), b_variable in b.items():
        if i % 2 and v % 2:
            quadratic_one[(i - 1 + u, j + v - 1)].append(
                z3.And(p_variable, b_variable)
            )
        if j % 2 and u % 2:
            quadratic_one[(i + u - 1, j - 1 + v)].append(
                z3.And(p_variable, b_variable)
            )
for (i, j), a_variable in a.items():
    for (u, v), q_variable in q.items():
        if i % 2 and v % 2:
            quadratic_one[(i - 1 + u, j + v - 1)].append(
                z3.And(a_variable, q_variable)
            )
        if j % 2 and u % 2:
            quadratic_one[(i + u - 1, j - 1 + v)].append(
                z3.And(a_variable, q_variable)
            )

for (ri, rj), r_variable in r.items():
    if ri % 2:
        for (i, j), a_variable in a.items():
            for u, v in QY:
                quadratic_one[(i + ri - 1 + u, j + rj + v)].append(
                    z3.And(a_variable, r_variable)
                )
        for (i, j), b_variable in b.items():
            for u, v in PY:
                quadratic_one[(i + ri - 1 + u, j + rj + v)].append(
                    z3.And(b_variable, r_variable)
                )
    if rj % 2:
        for (i, j), a_variable in a.items():
            for u, v in QX:
                quadratic_one[(i + ri + u, j + rj - 1 + v)].append(
                    z3.And(a_variable, r_variable)
                )
        for (i, j), b_variable in b.items():
            for u, v in PX:
                quadratic_one[(i + ri + u, j + rj - 1 + v)].append(
                    z3.And(b_variable, r_variable)
                )

for monomial in set(linear_one) | set(quadratic_one):
    if sum(monomial) <= args.degree + 8:
        continue
    bockstein = second_binary_digit(linear_one[monomial])
    solver.add(xor([bockstein] + quadratic_one[monomial]) == z3.BoolVal(False))

second_variables: list[tuple[str, dict[Monomial, z3.BoolRef]]] = []
determinant_one_constraint: z3.BoolRef | None = None
if args.seek_lift:
    assert args.degree >= 3

    # Complete the restricted ansatz with
    #   U=u0+z*u1+z^2*u2, V=v0+z*v1+z^2*v2,
    #   W=z*w1+z^3*w3.
    # A satisfying model is therefore an actual lift, not merely a survivor
    # of the high-degree necessary equations.
    u0 = variables("u0", monomials(args.degree))
    v0 = variables("v0", monomials(args.degree))
    w1 = variables("w1", monomials(args.degree - 1))
    u1 = variables("u1", monomials(args.degree - 1))
    v1 = variables("v1", monomials(args.degree - 1))
    u2 = variables("u2", monomials(args.degree - 2))
    v2 = variables("v2", monomials(args.degree - 2))
    w3 = variables("w3", monomials(args.degree - 3))
    second_variables = [
        ("u0", u0),
        ("v0", v0),
        ("w1", w1),
        ("u1", u1),
        ("v1", v1),
        ("u2", u2),
        ("v2", v2),
        ("w3", w3),
    ]
    second_zero: BooleanPolynomial = defaultdict(list)
    second_one: BooleanPolynomial = defaultdict(list)
    second_two: BooleanPolynomial = defaultdict(list)
    for polynomial, coordinate, output in (
        (u0, "A", second_zero),
        (v0, "B", second_zero),
        (u1, "A", second_one),
        (v1, "B", second_one),
        (u2, "A", second_two),
        (v2, "B", second_two),
    ):
        for (i, j), variable in polynomial.items():
            add_plane_variation(output, variable, i, j, coordinate)

    for monomial in (
        set(linear_zero) | set(quadratic_zero) | set(second_zero) | set(w1)
    ):
        terms = [second_binary_digit(linear_zero[monomial])]
        terms.extend(quadratic_zero[monomial])
        terms.extend(second_zero[monomial])
        if monomial in w1:
            terms.append(w1[monomial])
        equation = xor(terms) == z3.BoolVal(False)
        if monomial == (0, 0):
            determinant_one_constraint = equation
        else:
            solver.add(equation)

    for monomial in (
        set(linear_one) | set(quadratic_one) | set(second_one)
    ):
        terms = [second_binary_digit(linear_one[monomial])]
        terms.extend(quadratic_one[monomial])
        terms.extend(second_one[monomial])
        solver.add(xor(terms) == z3.BoolVal(False))

    quadratic_two: BooleanPolynomial = defaultdict(list)
    for (i, j), a_variable in a.items():
        for (u, v), b_variable in b.items():
            if i % 2 and v % 2:
                quadratic_two[(i - 1 + u, j + v - 1)].append(
                    z3.And(a_variable, b_variable)
                )
            if j % 2 and u % 2:
                quadratic_two[(i + u - 1, j - 1 + v)].append(
                    z3.And(a_variable, b_variable)
                )

    for monomial in set(quadratic_two) | set(second_two) | set(w3):
        terms = list(quadratic_two[monomial])
        terms.extend(second_two[monomial])
        if monomial in w3:
            terms.append(w3[monomial])
        solver.add(xor(terms) == z3.BoolVal(False))

w3_constraint_count = len(solver.assertions())
w4_top_class_coefficients: dict[Monomial3, z3.BitVecRef] = {}
if args.require_w4_class_zero or args.audit_w4_class_certificate is not None:
    w4_top_class_coefficients, w4_symbolic_term_count = (
        compile_w4_top_class_coefficients(
            (p, q, c, a, b, r),
            (u0, v0, w1, u1, v1, u2, v2, w3),
        )
    )
    print(
        "PASS: compiled next W_4 top class into "
        f"{len(w4_top_class_coefficients)} odd-odd-odd coefficient "
        f"expressions with {w4_symbolic_term_count} collected terms"
    )

if args.audit_w4_class_certificate is not None:
    certificate_data, certificate_supports = load_certificate(
        args.audit_w4_class_certificate
    )
    named_variables = [
        ("p", p),
        ("q", q),
        ("c", c),
        ("a", a),
        ("b", b),
        ("r", r),
        *second_variables,
    ]
    substitutions = [
        (variable, z3.BoolVal(monomial in set(certificate_supports[name])))
        for name, polynomial in named_variables
        for monomial, variable in polynomial.items()
    ]
    nonzero_class: list[Monomial3] = []
    for monomial, expression in w4_top_class_coefficients.items():
        value = z3.simplify(z3.substitute(expression, *substitutions))
        assert z3.is_bv_value(value), (monomial, value)
        assert value.as_long() in (0, 8), (monomial, value)
        if value.as_long() == 8:
            nonzero_class.append(monomial)
    encoding = json.dumps(nonzero_class, separators=(",", ":")).encode()
    digest = hashlib.sha256(encoding).hexdigest()
    if args.audit_w4_class_certificate.name == (
        "huq_kuruvilla_w3_degree19_witness_first_280_second_160.json"
    ):
        assert certificate_data["jacobian_constant_mod_8"] == 1
        assert len(nonzero_class) == 48
        assert (1, 1, 1) in nonzero_class
        assert digest == PREFERRED_W4_TOP_CLASS_SHA256
    print(
        "PASS: compiled next-class replay has "
        f"{len(nonzero_class)} nonzero monomials and SHA256 {digest}"
    )
    raise SystemExit

if args.require_w4_class_zero:
    for expression in w4_top_class_coefficients.values():
        solver.add(expression == z3.BitVecVal(0, 4))

w4_degree_boundary_coefficients: dict[Monomial3, z3.BitVecRef] = {}
w4_degree_boundary_compilation = ""
if args.require_w4_degree == 19:
    # Seed the joint search with several singleton dual obstructions.  Each
    # is a determinant target outside the degree-19 image of D_F.  The
    # targets are discovered independently from the pinned class-zero model,
    # then compiled symbolically for the whole master solution space.
    pinned_path = (
        Path(__file__).resolve().parents[1]
        / "artifacts"
        / "generated-results"
        / "huq_kuruvilla_w3_degree19_w4_class_zero.json"
    )
    pinned_data, pinned_supports = load_certificate(pinned_path)
    pinned_error = determinant_error_digit(
        pinned_supports,
        pinned_data["jacobian_constant_mod_8"],
    )
    pinned_rows, _, pinned_targets = full_w4_completion_rows(
        pinned_error,
        args.require_w4_degree,
    )
    pinned_singleton_targets = [
        target
        for (left, right), target in zip(pinned_rows, pinned_targets)
        if left == 0 and right == 1
    ]
    assert pinned_singleton_targets[0] == (0, 19, 1)

    structural_universe = structural_w4_target_universe(
        (p, q, c, a, b, r),
        (u0, v0, w1, u1, v1, u2, v2, w3),
    )
    structural_rows, structural_keys, structural_targets = (
        full_w4_completion_rows(
            {target: 1 for target in structural_universe},
            args.require_w4_degree,
        )
    )
    structural_holes = sorted(
        target
        for (left, _), target in zip(structural_rows, structural_targets)
        if target in structural_universe and left == 0
    )
    structural_universe_encoding = json.dumps(
        sorted(structural_universe),
        separators=(",", ":"),
    ).encode()
    structural_hole_encoding = json.dumps(
        structural_holes,
        separators=(",", ":"),
    ).encode()
    assert len(structural_universe) == 5396
    assert len(structural_holes) == 4340
    assert hashlib.sha256(structural_universe_encoding).hexdigest() == (
        "14d9435ce8a8fce698194c8a2f86f72bfd962bbfbbfb1bf17d7600478451d454"
    )
    assert hashlib.sha256(structural_hole_encoding).hexdigest() == (
        "de442207ad627a8202168496c37fcd2b9af7bb8cf03cbeb96bf90a662097ab99"
    )
    expected_layer_statistics = {
        0: (500, 275, 1540, 1223),
        1: (270, 177, 1485, 1228),
        2: (405, 229, 1431, 1165),
        3: (216, 145, 940, 724),
    }
    for z_degree, expected in expected_layer_statistics.items():
        variable_count = sum(
            1
            for name, source in structural_keys
            if (source[2] - 1 if name == "T" else source[2]) == z_degree
        )
        layer_rows = [
            (left, 0)
            for (left, _), target in zip(structural_rows, structural_targets)
            if target[2] == z_degree
        ]
        layer_rank = rank_affine_rows(layer_rows)
        universe_count = sum(
            target[2] == z_degree for target in structural_universe
        )
        hole_count = sum(target[2] == z_degree for target in structural_holes)
        actual_statistics = (
            variable_count,
            layer_rank,
            universe_count,
            hole_count,
        )
        assert actual_statistics == expected, (
            z_degree,
            actual_statistics,
            expected,
        )
    print(
        "PASS: structural degree-19 W_4 target universe has 5396 monomials; "
        "4340 are singleton codomain holes"
    )
    print(
        "AUDIT: structural-hole SHA256 "
        "de442207ad627a8202168496c37fcd2b9af7bb8cf03cbeb96bf90a662097ab99"
    )
    if args.w4_seed_source == "pinned":
        boundary_pool = pinned_singleton_targets
    else:
        boundary_pool = [
            target
            for target in structural_holes
            if args.w4_structural_z_layer is None
            or target[2] == args.w4_structural_z_layer
        ]
    assert args.w4_seed_cuts <= len(boundary_pool), (
        args.w4_seed_cuts,
        len(boundary_pool),
    )
    boundary_targets = set(boundary_pool[:args.w4_seed_cuts])
    if args.w4_quotient_compiler:
        (
            w4_degree_boundary_coefficients,
            quotient_statistics,
        ) = compile_w4_factored_target_coefficients(
            (p, q, c, a, b, r),
            (u0, v0, w1, u1, v1, u2, v2, w3),
            boundary_targets,
        )
        w4_degree_boundary_compilation = (
            f"{quotient_statistics['minor_coefficients']} shared minor "
            f"coefficients, {quotient_statistics['minor_products']} minor "
            f"products, and {quotient_statistics['final_products']} final "
            "products"
        )
    else:
        (
            w4_degree_boundary_coefficients,
            w4_degree_boundary_term_count,
        ) = compile_w4_target_coefficients(
            (p, q, c, a, b, r),
            (u0, v0, w1, u1, v1, u2, v2, w3),
            boundary_targets,
        )
        w4_degree_boundary_compilation = (
            f"{w4_degree_boundary_term_count} collected terms"
        )
    named_variables = [
        ("p", p),
        ("q", q),
        ("c", c),
        ("a", a),
        ("b", b),
        ("r", r),
        *second_variables,
    ]
    pinned_substitutions = [
        (variable, z3.BoolVal(monomial in set(pinned_supports[name])))
        for name, polynomial in named_variables
        for monomial, variable in polynomial.items()
    ]
    replay_targets = set(sorted(boundary_targets)[:min(8, len(boundary_targets))])
    for target, boundary_expression in w4_degree_boundary_coefficients.items():
        if target in replay_targets:
            pinned_value = z3.simplify(
                z3.substitute(boundary_expression, *pinned_substitutions)
            )
            assert z3.is_bv_value(pinned_value), (target, pinned_value)
            expected_pinned_value = 8 * int(target in pinned_error)
            assert pinned_value.as_long() == expected_pinned_value, (
                target,
                pinned_value,
                expected_pinned_value,
            )
        solver.add(
            z3.Extract(3, 3, boundary_expression) == z3.BitVecVal(0, 1)
        )
    if args.w4_quotient_compiler:
        preferred_path = (
            Path(__file__).resolve().parents[1]
            / "artifacts"
            / "generated-results"
            / "huq_kuruvilla_w3_degree19_witness_first_280_second_160.json"
        )
        _, preferred_supports = load_certificate(preferred_path)
        preferred_substitutions = [
            (
                variable,
                z3.BoolVal(monomial in set(preferred_supports[name])),
            )
            for name, polynomial in named_variables
            for monomial, variable in polynomial.items()
        ]
        preferred_determinant = determinant_polynomial(
            coordinates_from_supports(preferred_supports, 16),
            16,
        )
        for target in sorted(boundary_targets)[:8]:
            replayed = z3.simplify(
                z3.substitute(
                    w4_degree_boundary_coefficients[target],
                    *preferred_substitutions,
                )
            )
            assert z3.is_bv_value(replayed), (target, replayed)
            assert replayed.as_long() == preferred_determinant.get(target, 0), (
                target,
                replayed,
                preferred_determinant.get(target, 0),
            )
        print(
            "PASS: factored quotient coefficients agree with independent "
            "direct expansion on two pinned W_3 representatives"
        )
    print(
        f"PASS: compiled {len(boundary_targets)} unreachable W_4 boundary "
        f"coefficients from {w4_degree_boundary_compilation}"
    )
    print(
        f"PASS: {len(replay_targets)} seeded boundary coefficients agree "
        "with direct expansion on the pinned top-class-zero witness"
    )
    sorted_boundary_targets = sorted(boundary_targets)
    if len(sorted_boundary_targets) <= 64:
        print(f"AUDIT: seeded singleton targets {sorted_boundary_targets}")
    else:
        boundary_encoding = json.dumps(
            sorted_boundary_targets,
            separators=(",", ":"),
        ).encode()
        print(
            "AUDIT: seeded singleton-target SHA256 "
            f"{hashlib.sha256(boundary_encoding).hexdigest()}"
        )

if args.w4_fix_first_certificate is not None:
    _, fixed_w4_supports = load_certificate(args.w4_fix_first_certificate)
    for name, polynomial in (
        ("p", p),
        ("q", q),
        ("c", c),
        ("a", a),
        ("b", b),
        ("r", r),
    ):
        support = set(fixed_w4_supports[name])
        for monomial, variable in polynomial.items():
            solver.add(variable == z3.BoolVal(monomial in support))

variable_count = sum(map(len, (p, q, r, a, b, c))) + sum(
    len(polynomial) for _, polynomial in second_variables
)
constraint_count = len(solver.assertions())
selected_model: z3.ModelRef | None = None
optimization: dict | None = None
w4_completion_solution: set[int] | None = None
w4_completion_variable_keys: list[tuple[str, Monomial3]] = []
w4_completion_error: Polynomial3 = {}
w4_cut_limit_reached = False
w4_cut_records: list[list[Monomial3]] = (
    [[target] for target in sorted(boundary_targets)]
    if args.require_w4_degree == 19
    else []
)

if args.first_support_bound is not None:
    assert args.degree == 19
    first_support_variables = [
        variable
        for polynomial in (p, q, c, a, b, r)
        for variable in polynomial.values()
    ]
    solver.add(
        z3.PbLe(
            [(variable, 1) for variable in first_support_variables],
            args.first_support_bound,
        )
    )

if args.write_w4_dimacs is not None:
    target = args.write_w4_dimacs.resolve()
    assert not target.exists(), f"refusing to overwrite {target}"
    cnf_goal = z3.Goal()
    cnf_goal.add(*solver.assertions())
    cnf_subgoals = z3.Then(
        z3.Tactic("simplify"),
        z3.Tactic("bit-blast"),
        z3.Tactic("tseitin-cnf"),
    )(cnf_goal)
    assert len(cnf_subgoals) == 1
    dimacs = cnf_subgoals[0].dimacs()
    target.write_text(dimacs)
    header = next(line for line in dimacs.splitlines() if line.startswith("p "))
    print(f"PASS: wrote joint W_3/W_4 Boolean system to {target}")
    print(f"PASS: DIMACS {header}")
    raise SystemExit

if args.minimize_second is not None:
    certificate_data, fixed_supports = load_certificate(args.minimize_second)
    first_variables = [
        ("p", p),
        ("q", q),
        ("c", c),
        ("a", a),
        ("b", b),
        ("r", r),
    ]
    first_substitutions: list[tuple[z3.BoolRef, z3.BoolRef]] = []
    for name, polynomial in first_variables:
        support = set(fixed_supports[name])
        for monomial, variable in polynomial.items():
            value = z3.BoolVal(monomial in support)
            first_substitutions.append((variable, value))
            solver.add(variable == value)

    second_boolean_variables = [
        variable
        for _, polynomial in second_variables
        for variable in polynomial.values()
    ]
    layer_for_name = {
        "u0": "z0",
        "v0": "z0",
        "w1": "z0",
        "u1": "z1",
        "v1": "z1",
        "u2": "z2",
        "v2": "z2",
        "w3": "z2",
    }
    second_variable_layers = [
        layer_for_name[name]
        for name, polynomial in second_variables
        for _ in polynomial.values()
    ]
    layer_names = {
        "all": {name for name, _ in second_variables},
        "z0": {"u0", "v0", "w1"},
        "z1": {"u1", "v1"},
        "z2": {"u2", "v2", "w3"},
    }[args.second_support_layer]
    minimization_variables = [
        variable
        for name, polynomial in second_variables
        if name in layer_names
        for variable in polynomial.values()
    ]
    known_second_support = sum(
        certificate_data["support_counts"][name]
        for name in layer_names
    )
    affine_components_data: list[AffineComponent] = []
    if args.audit_second_linear or args.component_minimize_second:
        audit_assertions = list(solver.assertions())
        if (
            args.component_minimize_second
            and args.component_determinant == "one"
        ):
            assert determinant_one_constraint is not None
            audit_assertions.append(determinant_one_constraint)
        (
            equation_count,
            rank,
            nullity,
            layer_statistics,
            component_statistics,
            affine_components_data,
        ) = audit_affine_boolean_system(
            audit_assertions,
            second_boolean_variables,
            first_substitutions,
            second_variable_layers,
        )
        print(
            "PASS: fixed-first second correction is an affine F_2 system "
            f"with {len(second_boolean_variables)} variables, "
            f"{equation_count} nontrivial equations, rank {rank}, and "
            f"nullity {nullity}"
        )
        for layer, statistics in layer_statistics.items():
            layer_equations, layer_rank, layer_nullity = statistics
            print(
                f"PASS: {layer} block has {layer_equations} equations, "
                f"rank {layer_rank}, and nullity {layer_nullity}"
            )
            nontrivial_components = [
                component
                for component in component_statistics[layer]
                if component[1]
            ]
            largest = max(nontrivial_components, default=(0, 0, 0, 0))
            print(
                f"PASS: {layer} incidence graph has "
                f"{len(nontrivial_components)} nontrivial components; "
                "largest (variables,equations,rank,nullity) is "
                f"{largest}"
            )
            if layer == "z0" and args.audit_second_linear:
                print(f"AUDIT: z0 component statistics {nontrivial_components}")
    if args.audit_second_linear:
        replay_certificate(args.minimize_second)
        raise SystemExit

    if args.component_minimize_second:
        second_variable_keys = [
            (name, monomial)
            for name, polynomial in second_variables
            for monomial in polynomial
        ]
        known_true_indices = {
            index
            for index, (name, monomial) in enumerate(second_variable_keys)
            if monomial in set(fixed_supports[name])
        }
        target_layers = (
            {"z0", "z1", "z2"}
            if args.second_support_layer == "all"
            else {args.second_support_layer}
        )
        selected_indices = set(known_true_indices)
        component_records: dict[str, list[tuple[int, int, int, int]]] = {}
        total_lower = 0
        total_upper = 0
        all_exact = True
        for layer in sorted(target_layers):
            (
                layer_selected,
                layer_lower,
                layer_upper,
                layer_exact,
                records,
            ) = minimize_affine_components(
                affine_components_data,
                layer,
                known_true_indices,
                args.timeout_ms,
            )
            selected_indices.difference_update(
                index
                for index, variable_layer in enumerate(second_variable_layers)
                if variable_layer == layer
            )
            selected_indices.update(layer_selected)
            component_records[layer] = records
            total_lower += layer_lower
            total_upper += layer_upper
            all_exact &= layer_exact and layer_lower == layer_upper
            status = "exact" if layer_exact and layer_lower == layer_upper else "bounded"
            print(
                f"PASS: componentwise {layer} support is {status} at "
                f"{layer_lower}..{layer_upper}"
            )

        supports = {
            name: [tuple(monomial) for monomial in support]
            for name, support in fixed_supports.items()
        }
        for name, _ in second_variables:
            if layer_for_name[name] in target_layers:
                supports[name] = []
        for index in sorted(selected_indices):
            name, monomial = second_variable_keys[index]
            if layer_for_name[name] in target_layers:
                supports[name].append(monomial)

        support_counts = {
            name: len(support) for name, support in supports.items()
        }
        total_support = sum(support_counts.values())
        constant = verify_direct_supports(supports, args.degree)
        if args.component_determinant == "one":
            assert constant == 1
        second_support = sum(
            support_counts[name] for name, _ in second_variables
        )
        proved_global_minimum = target_layers == {"z0", "z1", "z2"} and all_exact
        optimization = {
            "fixed_first_certificate_sha256": hashlib.sha256(
                args.minimize_second.read_bytes()
            ).hexdigest(),
            "second_support": second_support,
            "second_support_lower_bound": (
                total_lower if proved_global_minimum else 0
            ),
            "proved_minimum": proved_global_minimum,
            "determinant": args.component_determinant,
            "componentwise_layer_bounds": {
                layer: {
                    "lower": sum(record[2] for record in records),
                    "upper": sum(record[3] for record in records),
                    "components": len(records),
                }
                for layer, records in component_records.items()
            },
        }
        print(
            f"PASS: componentwise witness uses {total_support} correction "
            f"coefficients with second support {second_support}"
        )
        print(
            "PASS: direct polynomial replay gives determinant "
            f"{constant} modulo 8"
        )

        if args.write_certificate is not None:
            repository_root = Path(__file__).resolve().parents[1]
            generated_directory = (
                repository_root / "artifacts" / "generated-results"
            ).resolve()
            target = args.write_certificate.resolve()
            assert target.parent == generated_directory
            assert not target.exists(), f"refusing to overwrite {target}"
            certificate = {
                "schema_version": 1,
                "claim": "HKM2W3",
                "degree": args.degree,
                "modulus": 8,
                "solver": f"Z3 {z3.get_version_string()} componentwise SAT",
                "system": {
                    "variables": variable_count,
                    "constraints": constraint_count,
                },
                "ansatz": certificate_data["ansatz"],
                "support_counts": support_counts,
                "total_support": total_support,
                "jacobian_constant_mod_8": constant,
                "supports": {
                    name: [list(monomial) for monomial in support]
                    for name, support in supports.items()
                },
                "optimization": optimization,
            }
            target.write_text(
                json.dumps(certificate, indent=2, sort_keys=True) + "\n"
            )
            print(f"PASS: wrote frozen support certificate {target}")
        raise SystemExit

    solver.push()
    solver.add(
        z3.PbLe(
            [(variable, 1) for variable in minimization_variables],
            known_second_support,
        )
    )
    result = solver.check()
    assert result == z3.sat, f"frozen witness failed to seed minimization: {result}"
    selected_model = solver.model()
    upper = sum(
        z3.is_true(selected_model.eval(variable, model_completion=True))
        for variable in minimization_variables
    )
    solver.pop()
    lower = 0
    proved_minimum = True

    if args.second_support_bound is not None:
        bound = args.second_support_bound
        assert 0 <= bound < upper
        solver.push()
        solver.add(
                z3.PbLe(
                    [(variable, 1) for variable in minimization_variables],
                    bound,
            )
        )
        trial = solver.check()
        if trial == z3.sat:
            selected_model = solver.model()
            upper = sum(
                z3.is_true(
                    selected_model.eval(variable, model_completion=True)
                )
                for variable in minimization_variables
            )
            print(f"TARGET second support <= {bound}: SAT ({upper} used)")
        elif trial == z3.unsat:
            lower = bound + 1
            print(f"TARGET second support <= {bound}: UNSAT")
        else:
            proved_minimum = False
            print(f"TARGET second support <= {bound}: UNKNOWN")
        solver.pop()

    while lower < upper and proved_minimum:
        bound = (lower + upper) // 2
        solver.push()
        solver.add(
            z3.PbLe(
                [(variable, 1) for variable in minimization_variables],
                bound,
            )
        )
        trial = solver.check()
        if trial == z3.sat:
            selected_model = solver.model()
            upper = sum(
                z3.is_true(
                    selected_model.eval(variable, model_completion=True)
                )
                for variable in minimization_variables
            )
            print(f"MINIMIZE second support <= {bound}: SAT ({upper} used)")
        elif trial == z3.unsat:
            lower = bound + 1
            print(f"MINIMIZE second support <= {bound}: UNSAT")
        else:
            proved_minimum = False
            print(f"MINIMIZE second support <= {bound}: UNKNOWN")
            solver.pop()
            break
        solver.pop()

    result = z3.sat
    optimization = {
        "fixed_first_certificate_sha256": hashlib.sha256(
            args.minimize_second.read_bytes()
        ).hexdigest(),
        "second_support": upper,
        "second_support_lower_bound": lower,
        "proved_minimum": proved_minimum and lower == upper,
    }
    status = "exact minimum" if optimization["proved_minimum"] else "best found"
    print(
        f"PASS: second-correction {args.second_support_layer} support "
        f"{status} is {upper}"
    )
else:
    if args.require_w4_degree == 19:
        named_joint_variables = [
            ("p", p),
            ("q", q),
            ("c", c),
            ("a", a),
            ("b", b),
            ("r", r),
            *second_variables,
        ]
        seen_cuts = {tuple(record) for record in w4_cut_records}
        while True:
            result = solver.check()
            if result != z3.sat:
                break
            selected_model = solver.model()
            candidate_supports = supports_from_model(
                selected_model,
                named_joint_variables,
            )
            candidate_constant = verify_direct_supports(
                candidate_supports,
                args.degree,
            )
            candidate_error = determinant_error_digit(
                candidate_supports,
                candidate_constant,
            )
            completion_rows, completion_keys, completion_targets = (
                full_w4_completion_rows(
                    candidate_error,
                    args.require_w4_degree,
                )
            )
            completion_solution = solve_affine_rows(completion_rows)
            print(
                "SEARCH: degree-19 W_4 completion subproblem is "
                f"{'SAT' if completion_solution is not None else 'UNSAT'} "
                f"({len(completion_keys)} variables, "
                f"{len(completion_rows)} equations, error support "
                f"{len(candidate_error)})"
            )
            if completion_solution is not None:
                w4_completion_solution = completion_solution
                w4_completion_variable_keys = completion_keys
                w4_completion_error = candidate_error
                break

            remaining_cut_slots = args.w4_cut_limit - len(w4_cut_records)
            singleton_targets = [
                target
                for (left, right), target in zip(
                    completion_rows,
                    completion_targets,
                )
                if left == 0
                and right == 1
                and (target,) not in seen_cuts
            ]
            if singleton_targets and remaining_cut_slots:
                selected_targets = sorted(singleton_targets)[:min(
                    8,
                    remaining_cut_slots,
                )]
                if args.w4_quotient_compiler:
                    singleton_coefficients, singleton_statistics = (
                        compile_w4_factored_target_coefficients(
                            (p, q, c, a, b, r),
                            (u0, v0, w1, u1, v1, u2, v2, w3),
                            set(selected_targets),
                        )
                    )
                    singleton_compilation = (
                        f"{singleton_statistics['minor_coefficients']} shared "
                        "minor coefficients and "
                        f"{singleton_statistics['minor_products'] + singleton_statistics['final_products']} "
                        "factored products"
                    )
                else:
                    singleton_coefficients, singleton_term_count = (
                        compile_w4_target_coefficients(
                            (p, q, c, a, b, r),
                            (u0, v0, w1, u1, v1, u2, v2, w3),
                            set(selected_targets),
                        )
                    )
                    singleton_compilation = (
                        f"{singleton_term_count} collected determinant terms"
                    )
                singleton_substitutions = [
                    (
                        variable,
                        z3.BoolVal(
                            monomial in set(candidate_supports[name])
                        ),
                    )
                    for name, polynomial in named_joint_variables
                    for monomial, variable in polynomial.items()
                ]
                for target in selected_targets:
                    bit = (
                        z3.Extract(3, 3, singleton_coefficients[target])
                        == z3.BitVecVal(1, 1)
                    )
                    violated = z3.simplify(
                        z3.substitute(bit, *singleton_substitutions)
                    )
                    assert z3.is_true(violated), (target, violated)
                    solver.add(bit == z3.BoolVal(False))
                    seen_cuts.add((target,))
                    w4_cut_records.append([target])
                print(
                    f"CUTS {len(w4_cut_records) - len(selected_targets) + 1}"
                    f"..{len(w4_cut_records)}: added "
                    f"{len(selected_targets)} singleton codomain-hole "
                    f"obstructions with {singleton_compilation}"
                )
                print(f"AUDIT: singleton targets {selected_targets}")
                continue

            dual_rows = affine_inconsistency_certificate(completion_rows)
            assert dual_rows is not None
            dual_left = 0
            dual_right = 0
            for row_index in dual_rows:
                left, right = completion_rows[row_index]
                dual_left ^= left
                dual_right ^= right
            assert dual_left == 0 and dual_right == 1
            cut_targets = sorted(completion_targets[index] for index in dual_rows)
            cut_key = tuple(cut_targets)
            assert cut_key not in seen_cuts
            if len(w4_cut_records) >= args.w4_cut_limit:
                w4_cut_limit_reached = True
                print(
                    "TARGET: certified cut limit reached with an unexcluded "
                    "master candidate"
                )
                break

            if args.w4_quotient_compiler:
                cut_coefficients, cut_statistics = (
                    compile_w4_factored_target_coefficients(
                        (p, q, c, a, b, r),
                        (u0, v0, w1, u1, v1, u2, v2, w3),
                        set(cut_targets),
                    )
                )
                cut_compilation = (
                    f"{cut_statistics['minor_coefficients']} shared minor "
                    "coefficients and "
                    f"{cut_statistics['minor_products'] + cut_statistics['final_products']} "
                    "factored products"
                )
            else:
                cut_coefficients, cut_term_count = compile_w4_target_coefficients(
                    (p, q, c, a, b, r),
                    (u0, v0, w1, u1, v1, u2, v2, w3),
                    set(cut_targets),
                )
                cut_compilation = f"{cut_term_count} collected determinant terms"
            cut_bits = [
                z3.Extract(3, 3, cut_coefficients[target])
                == z3.BitVecVal(1, 1)
                for target in cut_targets
            ]
            violated = z3.simplify(
                z3.substitute(xor(cut_bits), *[
                    (
                        variable,
                        z3.BoolVal(
                            monomial in set(candidate_supports[name])
                        ),
                    )
                    for name, polynomial in named_joint_variables
                    for monomial, variable in polynomial.items()
                ])
            )
            assert z3.is_true(violated), (cut_targets, violated)
            solver.add(xor(cut_bits) == z3.BoolVal(False))
            seen_cuts.add(cut_key)
            w4_cut_records.append(cut_targets)
            print(
                f"CUT {len(w4_cut_records)}: added a {len(cut_targets)}-"
                f"coefficient dual obstruction with {cut_compilation}"
            )
            print(f"AUDIT: dual targets {cut_targets}")
    else:
        result = solver.check()
        if result == z3.sat:
            selected_model = solver.model()

if args.degree == 18 and not args.seek_lift:
    assert variable_count == 1083
    assert w3_constraint_count == 1639
    assert result == z3.unsat, f"expected unsat, received {result}"
    print(
        "PASS: exact degree-18 W_3 necessary system is UNSAT "
        f"({variable_count} variables, {constraint_count} constraints)"
    )
    print(
        "PASS: every stable W_3(F_2) lift has maximum coordinate degree "
        "at least 19"
    )
elif args.degree == 19 and args.seek_lift:
    assert variable_count == 2685
    assert w3_constraint_count == 4513
    if args.require_w4_degree is not None:
        print(
            "TARGET joint degree-19 W_3/W_4 search: "
            f"{result} ({len(w4_cut_records)} certified cuts; "
            f"{len(solver.assertions())} final constraints)"
        )
        if result == z3.unsat:
            print(
                "PASS: no representative in the degree-19 W_3 existence "
                "ansatz has a degree-19 W_4 correction"
            )
        elif w4_completion_solution is not None:
            print(
                "PASS: a joint degree-19 W_3/W_4 witness was found"
            )
        elif w4_cut_limit_reached:
            print(
                "EXPERIMENT: cut limit reached before SAT/UNSAT was decided"
            )
        else:
            print("EXPERIMENT: solver returned unknown before a decision")
    elif args.require_w4_class_zero:
        scope = (
            "fixed-first"
            if args.w4_fix_first_certificate is not None
            else "unrestricted degree-19 ansatz"
        )
        print(
            f"TARGET W_4 top-class-zero {scope}: {result} "
            f"({constraint_count} active constraints)"
        )
    elif args.first_support_bound is None:
        assert result == z3.sat, f"expected sat, received {result}"
        print(
            "PASS: exact degree-19 W_3 lift system is SAT "
            f"({variable_count} variables, {constraint_count} constraints)"
        )
    else:
        print(
            "EXPERIMENT: degree-19 first support <= "
            f"{args.first_support_bound}: {result}"
        )
else:
    print(
        f"EXPERIMENT: degree {args.degree}: {result} "
        f"({variable_count} variables, {constraint_count} constraints)"
    )

if args.seek_lift and result == z3.sat:
    assert selected_model is not None
    model = selected_model
    named_variables = [
        ("p", p),
        ("q", q),
        ("c", c),
        ("a", a),
        ("b", b),
        ("r", r),
        *second_variables,
    ]
    supports = supports_from_model(model, named_variables)
    if args.decode_found_first:
        first_substitutions = [
            (variable, z3.BoolVal(monomial in set(supports[name])))
            for name, polynomial in (
                ("p", p),
                ("q", q),
                ("c", c),
                ("a", a),
                ("b", b),
                ("r", r),
            )
            for monomial, variable in polynomial.items()
        ]
        second_boolean_variables = [
            variable
            for _, polynomial in second_variables
            for variable in polynomial.values()
        ]
        layer_for_name = {
            "u0": "z0",
            "v0": "z0",
            "w1": "z0",
            "u1": "z1",
            "v1": "z1",
            "u2": "z2",
            "v2": "z2",
            "w3": "z2",
        }
        second_variable_layers = [
            layer_for_name[name]
            for name, polynomial in second_variables
            for _ in polynomial.values()
        ]
        audit_assertions = list(solver.assertions())
        assert determinant_one_constraint is not None
        audit_assertions.append(determinant_one_constraint)
        (*_, affine_components_data) = audit_affine_boolean_system(
            audit_assertions,
            second_boolean_variables,
            first_substitutions,
            second_variable_layers,
        )
        second_variable_keys = [
            (name, monomial)
            for name, polynomial in second_variables
            for monomial in polynomial
        ]
        known_true_indices = {
            index
            for index, (name, monomial) in enumerate(second_variable_keys)
            if monomial in set(supports[name])
        }
        selected_indices: set[int] = set()
        component_bounds: dict[str, dict[str, int]] = {}
        second_minimum = 0
        for layer in ("z0", "z1", "z2"):
            (
                layer_selected,
                layer_lower,
                layer_upper,
                layer_exact,
                records,
            ) = minimize_affine_components(
                affine_components_data,
                layer,
                known_true_indices,
                args.timeout_ms,
            )
            assert layer_exact and layer_lower == layer_upper
            selected_indices.update(layer_selected)
            second_minimum += layer_upper
            component_bounds[layer] = {
                "lower": layer_lower,
                "upper": layer_upper,
                "components": len(records),
            }
            print(f"PASS: decoded returned {layer} minimum is {layer_upper}")
        for name, _ in second_variables:
            supports[name] = []
        for index in sorted(selected_indices):
            name, monomial = second_variable_keys[index]
            supports[name].append(monomial)
        optimization = {
            "first_support_bound": args.first_support_bound,
            "first_support_random_seed": args.random_seed,
            "second_support": second_minimum,
            "second_support_lower_bound": second_minimum,
            "proved_minimum": True,
            "determinant": "one",
            "componentwise_layer_bounds": component_bounds,
        }
    if args.minimize_second is not None and args.second_support_layer != "all":
        for name, _ in second_variables:
            if name not in layer_names:
                supports[name] = fixed_supports[name]
    support_counts = {
        name: len(support) for name, support in supports.items()
    }
    if optimization is not None and args.second_support_layer != "all":
        total_second_support = sum(
            support_counts[name] for name, _ in second_variables
        )
        prior_lower = certificate_data.get("optimization", {}).get(
            "second_support_lower_bound", 0
        )
        optimization = {
            "fixed_first_certificate_sha256": hashlib.sha256(
                args.minimize_second.read_bytes()
            ).hexdigest(),
            "second_support": total_second_support,
            "second_support_lower_bound": prior_lower,
            "proved_minimum": False,
            "proved_layer_minimum": {
                args.second_support_layer: upper,
            },
        }
    total_support = sum(support_counts.values())
    first_support = sum(
        support_counts[name] for name in ("p", "q", "c", "a", "b", "r")
    )
    second_support = total_support - first_support
    constant = verify_direct_supports(
        supports,
        args.degree,
    )
    print(f"PASS: witness uses {total_support} correction coefficients")
    print(
        f"PASS: support splits as first {first_support} plus second "
        f"{second_support}"
    )
    print(
        "PASS: direct polynomial replay gives constant Jacobian "
        f"{constant} modulo 8"
    )
    assert constant * constant % 8 == 1
    print(
        "PASS: scaling one target coordinate by that odd constant gives "
        "determinant one"
    )
    if w4_completion_solution is not None:
        w4_coordinates = coordinates_from_supports(supports, 16)
        coordinate_for_name = {"R": 0, "S": 1, "T": 2}
        w4_support_counts = {"R": 0, "S": 0, "T": 0}
        for index in sorted(w4_completion_solution):
            name, monomial = w4_completion_variable_keys[index]
            add_term(
                w4_coordinates[coordinate_for_name[name]],
                monomial,
                8,
                16,
            )
            w4_support_counts[name] += 1
        assert max(
            max(map(sum, coordinate), default=0)
            for coordinate in w4_coordinates
        ) <= args.require_w4_degree
        w4_determinant = determinant_polynomial(w4_coordinates, 16)
        assert set(w4_determinant) == {(0, 0, 0)}, w4_determinant
        w4_constant = w4_determinant[(0, 0, 0)]
        assert w4_constant % 8 == constant
        print(
            "PASS: direct joint replay gives constant Jacobian "
            f"{w4_constant} modulo 16"
        )
        print(
            "PASS: degree-19 W_4 correction support is "
            f"R:{w4_support_counts['R']}, S:{w4_support_counts['S']}, "
            f"T:{w4_support_counts['T']}"
        )

    if args.write_certificate is not None:
        repository_root = Path(__file__).resolve().parents[1]
        generated_directory = (
            repository_root / "artifacts" / "generated-results"
        ).resolve()
        target = args.write_certificate.resolve()
        assert target.parent == generated_directory
        assert not target.exists(), f"refusing to overwrite {target}"
        certificate = {
            "schema_version": 1,
            "claim": "HKM2W3",
            "degree": args.degree,
            "modulus": 8,
            "solver": f"Z3 {z3.get_version_string()} SAT",
            "system": {
                "variables": variable_count,
                "constraints": w3_constraint_count,
            },
            "ansatz": {
                "first": "A=p+z*a, B=q+z*b, C=c+z*r",
                "second": (
                    "U=u0+z*u1+z^2*u2, V=v0+z*v1+z^2*v2, "
                    "W=z*w1+z^3*w3"
                ),
            },
            "support_counts": support_counts,
            "total_support": total_support,
            "jacobian_constant_mod_8": constant,
            "supports": {
                name: [list(monomial) for monomial in support]
                for name, support in supports.items()
            },
        }
        if optimization is not None:
            certificate["optimization"] = optimization
        if args.require_w4_class_zero:
            certificate["w4_top_class_zero"] = {
                "equations": len(w4_top_class_coefficients),
                "collected_symbolic_terms": w4_symbolic_term_count,
                "active_constraints": constraint_count,
                "proved": True,
            }
        target.write_text(json.dumps(certificate, indent=2, sort_keys=True) + "\n")
        print(f"PASS: wrote frozen support certificate {target}")

if args.show_model and result == z3.sat:
    if not args.seek_lift:
        model = solver.model()
        named_variables = [
            ("p", p),
            ("q", q),
            ("r", r),
            ("a", a),
            ("b", b),
            ("c", c),
            *second_variables,
        ]
        supports = supports_from_model(model, named_variables)
    for name, support in supports.items():
        print(f"MODEL {name} {support}")
