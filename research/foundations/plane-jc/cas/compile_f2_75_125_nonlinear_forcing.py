#!/usr/bin/env python3
"""Compile the nonlinear F2 carrier forcing as sparse arithmetic circuits.

This is the coupled continuation after the exact carrier specialization.  It
substitutes the ten fixed-endpoint pivot circuits into the genuine source-band
polynomials, forms every zero Laurent row on descents 12..35 and 37, and
projects the resulting nonlinear forcing into the *full* Laurent-row
cokernel over ``QQ(rho)``.  It then appends the seven target, six layer-zero
Hermite, five descent-eight incidence, and one quadratic defect equations.

The calculation deliberately uses a hash-consed arithmetic DAG.  Expanding
the 927-variable degree-eight presentation into an ordinary multivariate
coefficient list would recreate the coefficient explosion that the
kernel/cokernel implementation is meant to avoid.

The rational squarefree carrier is not on this component.  A separate audit
routes it into the later first-defect spacings 9..90 and records the exact
pre-target/target-row census without pretending that those branches have
already been compiled.

This produces an exact component presentation, not a Groebner computation,
unit ideal, F2 exclusion, or Keller map.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from math import comb, factorial
from pathlib import Path
from typing import Iterable

import sympy as sp
from sympy.polys.matrices import DomainMatrix

from classify_f2_75_125_layers import (
    band_factor_data,
    make_band,
    nonlinear_first_defect_audit,
)
from reduce_f2_75_125_endpoint_system import normalized_layer_support_audit
from verify_f2_75_125_carrier_specializations import (
    RHO_FIELD,
    RHO_POLYNOMIAL,
    RHO_SYMBOL,
    evaluation_row,
    field_element,
    field_matrix,
    full_q_operator,
    matrix_digest,
    polynomial_coefficient_matrix,
    reduced_left_cokernel,
)


ROOT = Path(__file__).resolve().parents[2]
ARTIFACT = (
    ROOT
    / "artifacts/generated-results/jc2_f2_75_125_nonlinear_forcing.json"
)


def _rational_key(value: object) -> tuple[int, int]:
    return int(value.numerator), int(value.denominator)


def field_key(value: object) -> tuple[tuple[int, int], tuple[int, int]]:
    coefficients = list(value.to_list())
    coefficients = [RHO_FIELD.dom.zero] * (2 - len(coefficients)) + coefficients
    rho_coefficient, constant = coefficients
    return _rational_key(constant), _rational_key(rho_coefficient)


class CircuitDAG:
    """A small canonical arithmetic DAG over ``QQ(rho)``."""

    def __init__(self) -> None:
        self.nodes: list[tuple[object, ...]] = []
        self.degrees: list[int] = []
        self.cache: dict[tuple[object, ...], int] = {}
        self.field_value_cache: dict[
            tuple[tuple[int, int], tuple[int, int]], object
        ] = {}
        self.field_key_cache: dict[
            object, tuple[tuple[int, int], tuple[int, int]]
        ] = {}
        self.zero = self.constant(RHO_FIELD.zero)
        self.one = self.constant(RHO_FIELD.one)

    def _intern(self, key: tuple[object, ...], degree: int) -> int:
        existing = self.cache.get(key)
        if existing is not None:
            return existing
        index = len(self.nodes)
        self.cache[key] = index
        self.nodes.append(key)
        self.degrees.append(degree)
        return index

    def constant(self, value: object) -> int:
        return self._intern(("const", self._field_key(value)), 0)

    def _field_key(
        self, value: object
    ) -> tuple[tuple[int, int], tuple[int, int]]:
        known = self.field_key_cache.get(value)
        if known is not None:
            return known
        key = field_key(value)
        self.field_key_cache[value] = key
        self.field_value_cache.setdefault(key, value)
        return key

    def variable(self, name: str) -> int:
        return self._intern(("var", name), 1)

    def scale(self, coefficient: object, node: int) -> int:
        if not coefficient or node == self.zero:
            return self.zero
        if coefficient == RHO_FIELD.one:
            return node
        kind = self.nodes[node][0]
        if kind == "const":
            node_value = self._field_from_key(self.nodes[node][1])
            return self.constant(coefficient * node_value)
        if kind == "scale":
            old_coefficient = self._field_from_key(self.nodes[node][1])
            return self.scale(coefficient * old_coefficient, self.nodes[node][2])
        return self._intern(
            ("scale", self._field_key(coefficient), node), self.degrees[node]
        )

    def linear_combination(
        self, terms: Iterable[tuple[object, int]]
    ) -> int:
        return self.add(
            *(self.scale(coefficient, node) for coefficient, node in terms)
        )

    def add(self, *nodes: int) -> int:
        active = [node for node in nodes if node != self.zero]
        if not active:
            return self.zero

        def add_pair(left: int, right: int) -> int:
            if left == self.zero:
                return right
            if right == self.zero:
                return left
            if left == right:
                return self.scale(RHO_FIELD.convert(2), left)
            if (
                self.nodes[right][0] == "scale"
                and self._field_from_key(self.nodes[right][1])
                == RHO_FIELD.convert(-1)
                and self.nodes[right][2] == left
            ) or (
                self.nodes[left][0] == "scale"
                and self._field_from_key(self.nodes[left][1])
                == RHO_FIELD.convert(-1)
                and self.nodes[left][2] == right
            ):
                return self.zero
            left, right = sorted((left, right))
            return self._intern(
                ("add", left, right),
                max(self.degrees[left], self.degrees[right]),
            )

        while len(active) > 1:
            active = [
                add_pair(active[index], active[index + 1])
                if index + 1 < len(active)
                else active[index]
                for index in range(0, len(active), 2)
            ]
        return active[0]

    def multiply(self, *nodes: int) -> int:
        active: list[int] = []
        coefficient = RHO_FIELD.one
        for node in nodes:
            if node == self.zero:
                return self.zero
            kind = self.nodes[node][0]
            if kind == "const":
                coefficient *= self._field_from_key(self.nodes[node][1])
            else:
                active.append(node)
        if not coefficient:
            return self.zero
        if not active:
            return self.constant(coefficient)

        def multiply_pair(left: int, right: int) -> int:
            left, right = sorted((left, right))
            return self._intern(
                ("mul", left, right),
                self.degrees[left] + self.degrees[right],
            )

        while len(active) > 1:
            active = [
                multiply_pair(active[index], active[index + 1])
                if index + 1 < len(active)
                else active[index]
                for index in range(0, len(active), 2)
            ]
        base = active[0]
        return self.scale(coefficient, base)

    def power(self, node: int, exponent: int) -> int:
        if exponent < 0:
            raise ValueError("circuit powers must be nonnegative")
        if exponent == 0:
            return self.one
        return self.multiply(*([node] * exponent))

    def _field_from_key(
        self, key: tuple[tuple[int, int], tuple[int, int]]
    ) -> object:
        known = self.field_value_cache.get(key)
        if known is not None:
            return known
        (constant_p, constant_q), (rho_p, rho_q) = key
        value = (
            RHO_FIELD.convert(sp.Rational(constant_p, constant_q))
            + RHO_FIELD.convert(sp.Rational(rho_p, rho_q)) * RHO_FIELD.unit
        )
        self.field_value_cache[key] = value
        return value

    def expression_node(
        self, expression: sp.Expr, symbol_nodes: dict[sp.Symbol, int]
    ) -> int:
        expression = sp.sympify(expression)
        if expression.is_Integer or expression.is_Rational:
            return self.constant(field_element(expression))
        if expression.is_Symbol:
            return symbol_nodes[expression]
        if expression.is_Add:
            return self.add(
                *(self.expression_node(term, symbol_nodes) for term in expression.args)
            )
        if expression.is_Mul:
            return self.multiply(
                *(self.expression_node(term, symbol_nodes) for term in expression.args)
            )
        if expression.is_Pow and expression.exp.is_Integer:
            return self.power(
                self.expression_node(expression.base, symbol_nodes),
                int(expression.exp),
            )
        raise TypeError(f"unsupported endpoint expression: {expression!r}")

    def structural_hash(self, node: int, memo: dict[int, str]) -> str:
        known = memo.get(node)
        if known is not None:
            return known
        record = self.nodes[node]
        kind = record[0]
        if kind in ("const", "var"):
            payload = record
        elif kind == "scale":
            payload = (
                "scale",
                record[1],
                self.structural_hash(record[2], memo),
            )
        elif kind in ("add", "mul"):
            payload = (
                kind,
                self.structural_hash(record[1], memo),
                self.structural_hash(record[2], memo),
            )
        else:
            raise AssertionError("unknown circuit node")
        digest = hashlib.sha256(
            json.dumps(payload, separators=(",", ":")).encode()
        ).hexdigest()
        memo[node] = digest
        return digest

    def equation_digest(self, equations: list[int]) -> str:
        memo: dict[int, str] = {}
        hashes = [self.structural_hash(node, memo) for node in equations]
        return hashlib.sha256("\n".join(hashes).encode()).hexdigest()

    def summary(self) -> dict[str, object]:
        operations = Counter(str(node[0]) for node in self.nodes)
        variables = [node[1] for node in self.nodes if node[0] == "var"]
        return {
            "node_count": len(self.nodes),
            "nodes_by_operation": dict(sorted(operations.items())),
            "variable_count": len(variables),
            "variable_name_digest_sha256": hashlib.sha256(
                "\n".join(sorted(variables)).encode()
            ).hexdigest(),
            "maximum_total_degree": max(self.degrees, default=0),
        }


Polynomial = list[int]


def poly_trim(dag: CircuitDAG, polynomial: Polynomial) -> Polynomial:
    result = list(polynomial)
    while len(result) > 1 and result[-1] == dag.zero:
        result.pop()
    return result


def poly_add(dag: CircuitDAG, left: Polynomial, right: Polynomial) -> Polynomial:
    length = max(len(left), len(right))
    return poly_trim(
        dag,
        [
            dag.add(
                left[index] if index < len(left) else dag.zero,
                right[index] if index < len(right) else dag.zero,
            )
            for index in range(length)
        ],
    )


def poly_scale(dag: CircuitDAG, scalar: object, polynomial: Polynomial) -> Polynomial:
    return poly_trim(dag, [dag.scale(scalar, node) for node in polynomial])


def poly_multiply(
    dag: CircuitDAG, left: Polynomial, right: Polynomial
) -> Polynomial:
    result = [dag.zero] * (len(left) + len(right) - 1)
    for left_degree, left_node in enumerate(left):
        if left_node == dag.zero:
            continue
        for right_degree, right_node in enumerate(right):
            if right_node == dag.zero:
                continue
            degree = left_degree + right_degree
            result[degree] = dag.add(
                result[degree], dag.multiply(left_node, right_node)
            )
    return poly_trim(dag, result)


def poly_derivative(dag: CircuitDAG, polynomial: Polynomial) -> Polynomial:
    if len(polynomial) == 1:
        return [dag.zero]
    return poly_trim(
        dag,
        [
            dag.scale(RHO_FIELD.convert(degree), polynomial[degree])
            for degree in range(1, len(polynomial))
        ],
    )


def poly_shift(dag: CircuitDAG, polynomial: Polynomial, shift: int) -> Polynomial:
    if shift >= 0:
        return [dag.zero] * shift + polynomial
    removed = -shift
    if any(node != dag.zero for node in polynomial[:removed]):
        raise AssertionError("a negative Kummer shift did not cancel")
    return poly_trim(dag, polynomial[removed:])


def poly_pad(
    dag: CircuitDAG, polynomial: Polynomial, length: int
) -> Polynomial:
    if len(polynomial) > length and any(
        node != dag.zero for node in polynomial[length:]
    ):
        raise AssertionError("a circuit polynomial exceeded its ambient degree")
    return polynomial[:length] + [dag.zero] * max(0, length - len(polynomial))


def poly_divmod_fixed(
    dag: CircuitDAG,
    dividend: Polynomial,
    monic_divisor: Polynomial,
) -> tuple[Polynomial, Polynomial]:
    """Divide a circuit polynomial by a fixed monic field polynomial."""

    divisor = poly_trim(dag, monic_divisor)
    if divisor[-1] != dag.one:
        raise AssertionError("the fixed circuit divisor is not monic")
    divisor_degree = len(divisor) - 1
    if divisor_degree <= 0:
        raise ValueError("the fixed divisor must have positive degree")
    remainder = list(dividend)
    quotient = [dag.zero] * max(1, len(remainder) - divisor_degree)
    while len(remainder) - 1 >= divisor_degree:
        shift = len(remainder) - 1 - divisor_degree
        leading = remainder[-1]
        quotient[shift] = dag.add(quotient[shift], leading)
        for degree, divisor_node in enumerate(divisor):
            target = shift + degree
            remainder[target] = dag.add(
                remainder[target],
                dag.scale(RHO_FIELD.convert(-1), dag.multiply(leading, divisor_node)),
            )
        remainder = poly_trim(dag, remainder)
    return poly_trim(dag, quotient), poly_pad(dag, remainder, divisor_degree)


def fixed_polynomial(dag: CircuitDAG, expression: sp.Expr, variable: sp.Symbol) -> Polynomial:
    polynomial = sp.Poly(sp.cancel(expression), variable)
    return [
        dag.constant(field_element(polynomial.coeff_monomial(variable**degree)))
        for degree in range(polynomial.degree() + 1)
    ]


def delta_basis_polynomial(
    dag: CircuitDAG, coefficients: list[int]
) -> Polynomial:
    result = [dag.zero]
    for order, coefficient in enumerate(coefficients):
        if coefficient == dag.zero:
            continue
        basis = [
            dag.constant(field_element(comb(order, degree) * (-1) ** (order - degree)))
            for degree in range(order + 1)
        ]
        result = poly_add(
            dag,
            result,
            [dag.multiply(coefficient, node) for node in basis],
        )
    return poly_trim(dag, result)


def evaluate_polynomial(
    dag: CircuitDAG, polynomial: Polynomial, point: object, order: int = 0
) -> int:
    terms: list[tuple[object, int]] = []
    for degree in range(order, len(polynomial)):
        coefficient = factorial(degree) // factorial(degree - order)
        scalar = RHO_FIELD.convert(coefficient) * point ** (degree - order)
        terms.append((scalar, polynomial[degree]))
    return dag.linear_combination(terms)


def bracket_pair(
    dag: CircuitDAG,
    p_layer: int,
    p_u: int,
    p_factor: Polynomial,
    q_layer: int,
    q_u: int,
    q_factor: Polynomial,
) -> Polynomial:
    q_term = poly_add(
        dag,
        poly_scale(dag, RHO_FIELD.convert(q_u), q_factor),
        poly_shift(
            dag,
            poly_scale(dag, RHO_FIELD.convert(5), poly_derivative(dag, q_factor)),
            1,
        ),
    )
    p_term = poly_add(
        dag,
        poly_scale(dag, RHO_FIELD.convert(p_u), p_factor),
        poly_shift(
            dag,
            poly_scale(dag, RHO_FIELD.convert(5), poly_derivative(dag, p_factor)),
            1,
        ),
    )
    result = poly_add(
        dag,
        poly_scale(
            dag,
            RHO_FIELD.convert(p_layer),
            poly_multiply(dag, p_factor, q_term),
        ),
        poly_scale(
            dag,
            RHO_FIELD.convert(-q_layer),
            poly_multiply(dag, p_term, q_factor),
        ),
    )
    shift = (p_u + q_u - 1) // 5
    return poly_shift(dag, result, shift)


def endpoint_sympy_circuits(
    p_bands: dict[int, object], q_bands: dict[int, object]
) -> dict[str, object]:
    """Rebuild the ten exact endpoint solutions for circuit substitution."""

    delta = sp.Symbol("delta")
    local_w = 1 + delta
    jet_requirements: dict[str, dict[int, int]] = {"P": {}, "Q": {}}

    def require(side: str, layer: int, maximum_jet: int) -> None:
        jet_requirements[side][layer] = max(
            jet_requirements[side].get(layer, -1), maximum_jet
        )

    target_pairs: list[tuple[int, int]] = []
    for p_layer in range(-21, 16):
        q_layer = 4 - p_layer
        _, p_vanishing, _ = band_factor_data(p_bands[p_layer])
        _, q_vanishing, _ = band_factor_data(q_bands[q_layer])
        if p_vanishing + q_vanishing > 5:
            continue
        maximum_jet = 5 - p_vanishing - q_vanishing
        require("P", p_layer, maximum_jet)
        require("Q", q_layer, maximum_jet)
        target_pairs.append((p_layer, q_layer))

    layer_zero_pairs: list[tuple[int, int]] = []
    for p_layer in range(-25, 16):
        if p_layer == 0:
            continue
        q_layer = -p_layer
        _, p_vanishing, _ = band_factor_data(p_bands[p_layer])
        _, q_vanishing, _ = band_factor_data(q_bands[q_layer])
        if p_vanishing + q_vanishing > 5:
            continue
        maximum_jet = 5 - p_vanishing - q_vanishing
        require("P", p_layer, maximum_jet)
        require("Q", q_layer, maximum_jet)
        layer_zero_pairs.append((p_layer, q_layer))

    if len(target_pairs) != 22 or len(layer_zero_pairs) != 25:
        raise AssertionError("the endpoint dependency cone changed")

    coefficient_lists: dict[tuple[str, int], tuple[sp.Symbol, ...]] = {}
    local_polynomials: dict[tuple[str, int], sp.Expr] = {}
    for side in ("P", "Q"):
        for layer, maximum_jet in sorted(jet_requirements[side].items()):
            layer_name = f"m{-layer}" if layer < 0 else str(layer)
            coefficients = sp.symbols(
                f"ep_{side.lower()}{layer_name}_0:{maximum_jet + 1}"
            )
            coefficient_lists[side, layer] = coefficients
            local_polynomials[side, layer] = sum(
                coefficient * delta**order
                for order, coefficient in enumerate(coefficients)
            )

    fixed_values = {
        coefficient_lists["P", 3][0]: sp.Rational(1, 5),
        coefficient_lists["Q", 1][0]: sp.Integer(-1),
        coefficient_lists["Q", 13][0]: -sp.Rational(3, 5**5),
    }

    def target_pair_expression(p_layer: int, q_layer: int) -> sp.Expr:
        p_u, p_vanishing, _ = band_factor_data(p_bands[p_layer])
        q_u, q_vanishing, _ = band_factor_data(q_bands[q_layer])
        p_factor = delta**p_vanishing * local_polynomials["P", p_layer]
        q_factor = delta**q_vanishing * local_polynomials["Q", q_layer]
        shift = (p_u + q_u - 1) // 5
        return sp.expand(
            local_w**shift
            * (
                p_layer
                * p_factor
                * (q_u * q_factor + 5 * local_w * sp.diff(q_factor, delta))
                - q_layer
                * (p_u * p_factor + 5 * local_w * sp.diff(p_factor, delta))
                * q_factor
            )
        )

    def layer_zero_expression(p_layer: int, q_layer: int) -> sp.Expr:
        p_u, p_vanishing, _ = band_factor_data(p_bands[p_layer])
        q_u, q_vanishing, _ = band_factor_data(q_bands[q_layer])
        shift = (p_u + q_u) // 5
        return sp.expand(
            p_layer
            * local_w**shift
            * delta ** (p_vanishing + q_vanishing)
            * local_polynomials["P", p_layer]
            * local_polynomials["Q", q_layer]
        )

    full_target = sp.expand(
        sum(target_pair_expression(*pair) for pair in target_pairs).subs(
            fixed_values
        )
    )
    full_integral = sp.expand(
        sum(layer_zero_expression(*pair) for pair in layer_zero_pairs).subs(
            fixed_values
        )
    )
    h_at_zero = sp.Symbol("H_at_zero")
    rows = [full_target.coeff(delta, order) for order in range(1, 5)] + [
        full_integral.coeff(delta, 0) - h_at_zero,
        *[full_integral.coeff(delta, order) for order in range(1, 6)],
    ]
    pivots = (
        *coefficient_lists["P", 3][1:5],
        *coefficient_lists["P", -1][:6],
    )
    pivot_zero = {pivot: 0 for pivot in pivots}
    matrix = sp.Matrix([[sp.diff(row, pivot) for pivot in pivots] for row in rows])
    forcing = sp.Matrix([sp.expand(row.subs(pivot_zero)) for row in rows])
    if sp.factor(matrix.det()) != 75000:
        raise AssertionError("the endpoint circuit determinant changed")
    solutions = [sp.expand(value) for value in (-matrix.inv() * forcing)]
    if sum(len(sp.Add.make_args(value)) for value in solutions) != 1489:
        raise AssertionError("the endpoint circuit term census changed")
    return {
        "coefficient_lists": coefficient_lists,
        "fixed_values": fixed_values,
        "pivots": pivots,
        "solutions": solutions,
        "H_at_zero": h_at_zero,
    }


def build_band_circuits(dag: CircuitDAG) -> dict[str, object]:
    """Construct every source band and substitute the ten endpoint circuits."""

    w = sp.Symbol("w")
    p_bands = {
        layer: make_band("P", 75, 3, layer) for layer in range(-75, 16)
    }
    q_bands = {
        layer: make_band("Q", 125, 5, layer) for layer in range(-125, 26)
    }
    endpoint = endpoint_sympy_circuits(p_bands, q_bands)

    delta_coefficients: dict[tuple[str, int], list[int]] = {}
    fixed_endpoint_values = {
        ("P", 3, 0): field_element(sp.Rational(1, 5)),
        ("Q", 1, 0): field_element(-1),
        ("Q", 13, 0): field_element(-sp.Rational(3, 5**5)),
    }
    pivot_orders = {
        ("P", 3): set(range(1, 5)),
        ("P", -1): set(range(6)),
    }
    for side, bands, interval in (
        ("P", p_bands, range(-25, 16)),
        ("Q", q_bands, range(-15, 26)),
    ):
        for layer in interval:
            _, _, degree = band_factor_data(bands[layer])
            coefficients: list[int] = []
            for order in range(degree + 1):
                fixed = fixed_endpoint_values.get((side, layer, order))
                if fixed is not None:
                    coefficients.append(dag.constant(fixed))
                elif order in pivot_orders.get((side, layer), set()):
                    coefficients.append(dag.zero)
                else:
                    coefficients.append(dag.variable(f"{side}_{layer}_d{order}"))
            delta_coefficients[side, layer] = coefficients

    R = (w - RHO_SYMBOL) ** 2 / (25 * (1 - RHO_SYMBOL) ** 2)
    top_factors = {
        ("P", 15): fixed_polynomial(dag, (w - 1) ** 6 * R**3, w),
        ("Q", 25): fixed_polynomial(
            dag, -sp.Rational(9, 5) * (w - 1) ** 10 * R**5, w
        ),
    }

    def factor_from_coefficients(side: str, layer: int) -> Polynomial:
        top = top_factors.get((side, layer))
        if top is not None:
            return top
        bands = p_bands if side == "P" else q_bands
        _, vanishing, _ = band_factor_data(bands[layer])
        K = delta_basis_polynomial(dag, delta_coefficients[side, layer])
        vanishing_factor = fixed_polynomial(dag, (w - 1) ** vanishing, w)
        return poly_multiply(dag, vanishing_factor, K)

    baseline_factors = {
        (side, layer): factor_from_coefficients(side, layer)
        for side, interval in (("P", range(-25, 16)), ("Q", range(-15, 26)))
        for layer in interval
    }

    # Descents 1..11 have already been eliminated by the exact upper
    # triangular block.  On the common-power locus their full tangent
    # solution is Q_follow=-3*C0^2*P_follow.  The resonant descents 5 and 10
    # retain the one-dimensional centralizers C0^4 and C0^3.
    C0 = fixed_polynomial(dag, (w - 1) ** 2 * R, w)
    C0_squared = poly_multiply(dag, C0, C0)

    def substitute_upper_tangents(
        factor_table: dict[tuple[str, int], Polynomial]
    ) -> None:
        for descent in range(1, 12):
            p_layer = 15 - descent
            q_layer = 25 - descent
            q_factor = poly_scale(
                dag,
                RHO_FIELD.convert(-3),
                poly_multiply(
                    dag, C0_squared, factor_table["P", p_layer]
                ),
            )
            if descent in (5, 10):
                exponent = 5 - descent // 5
                centralizer = dag.variable(f"upper_C0_power_{exponent}")
                centralizer_factor = fixed_polynomial(
                    dag, ((w - 1) ** 2 * R) ** exponent, w
                )
                q_factor = poly_add(
                    dag,
                    q_factor,
                    [
                        dag.multiply(centralizer, coefficient)
                        for coefficient in centralizer_factor
                    ],
                )
            factor_table["Q", q_layer] = q_factor

    substitute_upper_tangents(baseline_factors)

    # H(0) is independent of the ten follower pivots because every follower
    # is divisible by w^2.  Compute it on the pivot-zero baseline.
    H0 = dag.zero
    for p_layer in range(-25, 16):
        if p_layer == 0:
            continue
        q_layer = -p_layer
        p_u, _, _ = band_factor_data(p_bands[p_layer])
        q_u, _, _ = band_factor_data(q_bands[q_layer])
        product = poly_multiply(
            dag,
            baseline_factors["P", p_layer],
            baseline_factors["Q", q_layer],
        )
        product = poly_shift(dag, product, (p_u + q_u) // 5)
        H0 = dag.add(
            H0,
            dag.scale(
                RHO_FIELD.convert(p_layer),
                evaluate_polynomial(dag, product, RHO_FIELD.zero),
            ),
        )

    symbol_nodes: dict[sp.Symbol, int] = {
        endpoint["H_at_zero"]: H0
    }
    for (side, layer), symbols in endpoint["coefficient_lists"].items():
        coefficients = delta_coefficients[side, layer]
        for order, symbol in enumerate(symbols):
            fixed = endpoint["fixed_values"].get(symbol)
            symbol_nodes[symbol] = (
                dag.constant(field_element(fixed))
                if fixed is not None
                else coefficients[order]
            )

    solution_nodes = [
        dag.expression_node(expression, symbol_nodes)
        for expression in endpoint["solutions"]
    ]
    K_polynomials = {
        (side, layer): delta_basis_polynomial(
            dag, delta_coefficients[side, layer]
        )
        for side, interval in (("P", range(-25, 16)), ("Q", range(-15, 26)))
        for layer in interval
    }
    for (side, layer, orders), nodes in (
        (("P", 3, range(1, 5)), solution_nodes[:4]),
        (("P", -1, range(6)), solution_nodes[4:]),
    ):
        for order, solution in zip(orders, nodes):
            lift_expression = (
                (w - 1) ** order
                + (order - 7) * (-1) ** order * (w - 1) ** 6
                + (order - 6) * (-1) ** order * (w - 1) ** 7
            )
            lift = fixed_polynomial(dag, sp.expand(lift_expression), w)
            K_polynomials[side, layer] = poly_add(
                dag,
                K_polynomials[side, layer],
                [dag.multiply(solution, coefficient) for coefficient in lift],
            )

    factors = dict(baseline_factors)
    for side, layer in (("P", 3), ("P", -1)):
        _, vanishing, _ = band_factor_data(p_bands[layer])
        factors[side, layer] = poly_multiply(
            dag,
            fixed_polynomial(dag, (w - 1) ** vanishing, w),
            K_polynomials[side, layer],
        )

    full_H0 = dag.zero
    for p_layer in range(-25, 16):
        if p_layer == 0:
            continue
        q_layer = -p_layer
        p_u, _, _ = band_factor_data(p_bands[p_layer])
        q_u, _, _ = band_factor_data(q_bands[q_layer])
        product = poly_multiply(
            dag, factors["P", p_layer], factors["Q", q_layer]
        )
        product = poly_shift(dag, product, (p_u + q_u) // 5)
        full_H0 = dag.add(
            full_H0,
            dag.scale(
                RHO_FIELD.convert(p_layer),
                evaluate_polynomial(dag, product, RHO_FIELD.zero),
            ),
        )
    if full_H0 != H0:
        raise AssertionError("the endpoint follower circuits changed H(0)")

    return {
        "p_bands": p_bands,
        "q_bands": q_bands,
        "delta_coefficients": delta_coefficients,
        "K_polynomials": K_polynomials,
        "factors": factors,
        "endpoint_solution_nodes": solution_nodes,
        "H0": H0,
        "C0": C0,
    }


def project_rows(
    dag: CircuitDAG, rows: DomainMatrix, vector: Polynomial
) -> list[int]:
    if rows.shape[1] != len(vector):
        raise AssertionError("a pinned functional has the wrong ambient width")
    return [
        dag.linear_combination(zip(row, vector))
        for row in rows.to_list()
    ]


def coupled_quotient_operator(
    descent: int,
    maximum_full_degree: int,
) -> tuple[DomainMatrix, sp.Expr, int, int]:
    """Return the endpoint-domain quotient operator and forced divisor.

    A negative Kummer shift is handled by multiplying the entire Laurent row
    by ``w`` before division.  At the two fixed endpoint coordinates the
    movable Q factor has one extra ``w-1``.  These conventions give exactly
    the 53 pinned quotient coordinates, with the complementary 294
    divisibility coordinates retained separately.
    """

    w = sp.Symbol("w")
    p_layer = 15
    q_layer = 25 - descent
    p_band = make_band("P", 75, 3, p_layer)
    q_band = make_band("Q", 125, 5, q_layer)
    p_u, _, _ = band_factor_data(p_band)
    q_u, q_vanishing, q_degree = band_factor_data(q_band)
    shift = (p_u + q_u - 1) // 5
    endpoint_fixed = q_layer in (1, 13)
    orders = list(range(q_degree + 1))
    if endpoint_fixed:
        orders.remove(0)

    divisor = sp.expand(
        w ** max(shift, 0)
        * (w - 1) ** (q_vanishing + 5 + int(endpoint_fixed))
        * (w - RHO_SYMBOL) ** 5
    )
    quotient_maximum_degree = (
        maximum_full_degree + max(-shift, 0) - sp.degree(divisor, w)
    )
    # Build the factor-reduced columns from the three-term normalized
    # operator.  This avoids repeatedly dividing high-degree expressions in
    # QQ(rho)[w].  Changing from w^k to the endpoint basis (w-1)^j is an
    # exact triangular basis change.  If K(1) is fixed, the remaining image
    # is divided by its additional factor w-1.
    monomial_columns: list[sp.Expr] = []
    for degree in range(q_degree + 1):
        A = RHO_SYMBOL * (5 * degree + q_u)
        B = (
            (RHO_SYMBOL + 1) * (2 * q_layer - 5 * degree - q_u)
            - 5 * q_vanishing * RHO_SYMBOL
        )
        C = 5 * degree + q_u + 5 * q_vanishing - 4 * q_layer
        monomial_columns.append(
            A * w**degree + B * w ** (degree + 1) + C * w ** (degree + 2)
        )

    columns: list[sp.Expr] = []
    for order in orders:
        column = sp.expand(
            sum(
                comb(order, degree)
                * (-1) ** (order - degree)
                * monomial_columns[degree]
                for degree in range(order + 1)
            )
        )
        if endpoint_fixed:
            column, remainder = sp.div(column, w - 1, w)
            if remainder != 0:
                raise AssertionError("the fixed endpoint quotient changed")
        columns.append(sp.expand(column))

    operator = polynomial_coefficient_matrix(
        columns, w, int(quotient_maximum_degree)
    )
    return operator, divisor, shift, int(quotient_maximum_degree)


def coupled_forcing_audit(
    dag: CircuitDAG, bands: dict[str, object]
) -> dict[str, object]:
    """Compile descents 12..35 and 37 into 294+53 exact equations."""

    w = sp.Symbol("w")
    p_bands = bands["p_bands"]
    q_bands = bands["q_bands"]
    factors = bands["factors"]
    support = normalized_layer_support_audit(p_bands, q_bands)
    support_by_layer = {int(row["layer"]): row for row in support["rows"]}

    all_divisibility: list[int] = []
    all_quotient: list[int] = []
    records: list[dict[str, object]] = []
    for descent in [*range(12, 36), 37]:
        layer = 40 - descent
        q_new_layer = 25 - descent
        row = support_by_layer[layer]
        minimum_degree, maximum_degree = map(int, row["w_degree_interval"])

        forcing = [dag.zero]
        for p_layer in range(max(-22, layer - 25), 16):
            q_layer = layer - p_layer
            if not -12 <= q_layer <= 25:
                continue
            if p_layer == 15 and q_layer == q_new_layer:
                continue
            p_u, _, _ = band_factor_data(p_bands[p_layer])
            q_u, _, _ = band_factor_data(q_bands[q_layer])
            forcing = poly_add(
                dag,
                forcing,
                bracket_pair(
                    dag,
                    p_layer,
                    p_u,
                    factors["P", p_layer],
                    q_layer,
                    q_u,
                    factors["Q", q_layer],
                ),
            )

        # Q_13(1) and Q_1(1) are constants, not movable new-Q coordinates.
        # Their top-pair contribution belongs to the forcing.
        fixed_new_q = {
            13: -sp.Rational(3, 5**5),
            1: sp.Integer(-1),
        }.get(q_new_layer)
        if fixed_new_q is not None:
            q_u, q_vanishing, _ = band_factor_data(q_bands[q_new_layer])
            fixed_factor = fixed_polynomial(
                dag, fixed_new_q * (w - 1) ** q_vanishing, w
            )
            p_u, _, _ = band_factor_data(p_bands[15])
            forcing = poly_add(
                dag,
                forcing,
                bracket_pair(
                    dag,
                    15,
                    p_u,
                    factors["P", 15],
                    q_new_layer,
                    q_u,
                    fixed_factor,
                ),
            )

        forcing = poly_pad(dag, forcing, maximum_degree + 1)
        if any(node != dag.zero for node in forcing[:minimum_degree]):
            raise AssertionError("a coupled forcing escaped its support interval")

        operator, divisor_expression, shift, quotient_maximum = (
            coupled_quotient_operator(descent, maximum_degree)
        )
        cokernel = reduced_left_cokernel(operator)
        transformed = poly_shift(dag, forcing, max(-shift, 0))
        divisor = fixed_polynomial(dag, divisor_expression, w)
        quotient, remainder = poly_divmod_fixed(dag, transformed, divisor)
        quotient = poly_pad(dag, quotient, quotient_maximum + 1)

        # For shift=-1 the inserted leading w makes remainder coefficient zero
        # tautological; it is not a Laurent compatibility coordinate.
        remainder_start = max(-shift, 0)
        divisibility_coordinates = remainder[remainder_start:]
        quotient_coordinates = project_rows(dag, cokernel, quotient)
        full_cokernel_count = (
            maximum_degree - minimum_degree + 1 - operator.rank()
        )
        if (
            len(divisibility_coordinates) + len(quotient_coordinates)
            != full_cokernel_count
        ):
            raise AssertionError("the 294+53 Laurent cokernel split changed")

        all_divisibility.extend(divisibility_coordinates)
        all_quotient.extend(quotient_coordinates)
        equations = divisibility_coordinates + quotient_coordinates
        records.append(
            {
                "descent": descent,
                "Laurent_layer": layer,
                "new_Q_band": q_new_layer,
                "forcing_support": [minimum_degree, maximum_degree],
                "Kummer_shift": shift,
                "forced_divisor": str(sp.factor(divisor_expression)),
                "divisibility_coordinate_count": len(divisibility_coordinates),
                "quotient_operator_shape": list(operator.shape),
                "quotient_operator_rank": operator.rank(),
                "quotient_cokernel_coordinate_count": len(quotient_coordinates),
                "full_Laurent_cokernel_coordinate_count": full_cokernel_count,
                "quotient_operator_digest_sha256": matrix_digest(operator),
                "quotient_cokernel_basis_digest_sha256": matrix_digest(cokernel),
                "equation_digest_sha256": dag.equation_digest(equations),
                "identically_zero_coordinates_after_substitution": sum(
                    node == dag.zero for node in equations
                ),
                "maximum_circuit_total_degree": max(
                    (dag.degrees[node] for node in equations), default=0
                ),
            }
        )

    if (len(all_divisibility), len(all_quotient)) != (294, 53):
        raise AssertionError("the coupled 294+53 coordinate census changed")
    return {
        "coefficient_algebra": (
            "QQ(rho)[y]/(27*y^2-9*y+1), the quartic compositum presented "
            "relatively over QQ(rho)"
        ),
        "projection_order": [[12, 35], [37, 37]],
        "upper_tangent_substitution": (
            "Q_(25-delta)=-3*C0^2*P_(15-delta), delta=1..11, "
            "with C0^4 and C0^3 resonant parameters"
        ),
        "rows": records,
        "divisibility_coordinate_count": len(all_divisibility),
        "pinned_quotient_cokernel_coordinate_count": len(all_quotient),
        "full_Laurent_cokernel_coordinate_count": (
            len(all_divisibility) + len(all_quotient)
        ),
        "divisibility_equation_digest_sha256": dag.equation_digest(
            all_divisibility
        ),
        "quotient_equation_digest_sha256": dag.equation_digest(all_quotient),
        "all_equation_nodes": all_divisibility + all_quotient,
    }


def target_and_layer_zero_audit(
    dag: CircuitDAG, bands: dict[str, object]
) -> dict[str, object]:
    """Append the seven target and six residual Hermite functionals."""

    w = sp.Symbol("w")
    p_bands = bands["p_bands"]
    q_bands = bands["q_bands"]
    factors = bands["factors"]

    target_polynomial = [dag.zero]
    for p_layer in range(-21, 16):
        q_layer = 4 - p_layer
        p_u, _, _ = band_factor_data(p_bands[p_layer])
        q_u, _, _ = band_factor_data(q_bands[q_layer])
        target_polynomial = poly_add(
            dag,
            target_polynomial,
            bracket_pair(
                dag,
                p_layer,
                p_u,
                factors["P", p_layer],
                q_layer,
                q_u,
                factors["Q", q_layer],
            ),
        )
    target_polynomial[0] = dag.add(
        target_polynomial[0], dag.scale(RHO_FIELD.convert(-1), dag.one)
    )
    target_polynomial = poly_pad(dag, target_polynomial, 34)

    target_jets = [
        evaluate_polynomial(
            dag, target_polynomial, RHO_FIELD.unit, derivative_order
        )
        for derivative_order in range(5)
    ]
    target_divisor = fixed_polynomial(
        dag, w**2 * (w - 1) ** 5 * (w - RHO_SYMBOL) ** 5, w
    )
    target_quotient, _ = poly_divmod_fixed(
        dag, target_polynomial, target_divisor
    )
    target_quotient = poly_pad(dag, target_quotient, 22)

    E = (w - 1) * (w - RHO_SYMBOL)
    target_columns = [
        sp.expand(
            5 * w * E * sp.diff(w**degree, w)
            + (11 * E + 22 * w * sp.diff(E, w)) * w**degree
        )
        for degree in range(20)
    ]
    target_operator = polynomial_coefficient_matrix(target_columns, w, 21)
    target_cokernel = reduced_left_cokernel(target_operator)
    if target_cokernel.shape != (2, 22):
        raise AssertionError("the two target quotient residues changed")
    target_residues = project_rows(dag, target_cokernel, target_quotient)
    target_equations = target_jets + target_residues

    H = [dag.zero]
    for p_layer in range(-25, 16):
        if p_layer == 0:
            continue
        q_layer = -p_layer
        p_u, _, _ = band_factor_data(p_bands[p_layer])
        q_u, _, _ = band_factor_data(q_bands[q_layer])
        product = poly_multiply(
            dag, factors["P", p_layer], factors["Q", q_layer]
        )
        product = poly_shift(dag, product, (p_u + q_u) // 5)
        H = poly_add(
            dag,
            H,
            poly_scale(dag, RHO_FIELD.convert(p_layer), product),
        )
    H = poly_pad(dag, H, 34)
    H0 = bands["H0"]
    hermite_equations = [
        dag.add(
            evaluate_polynomial(dag, H, RHO_FIELD.unit),
            dag.scale(RHO_FIELD.convert(-1), H0),
        )
    ] + [
        evaluate_polynomial(dag, H, RHO_FIELD.unit, derivative_order)
        for derivative_order in range(1, 6)
    ]
    equations = target_equations + hermite_equations
    if len(equations) != 13:
        raise AssertionError("the final 7+6 functional count changed")
    return {
        "target": {
            "movable_rho_jet_orders": list(range(5)),
            "quotient_residue_count": len(target_residues),
            "equation_count": len(target_equations),
            "quotient_operator_digest_sha256": matrix_digest(target_operator),
            "quotient_cokernel_basis_digest_sha256": matrix_digest(
                target_cokernel
            ),
            "equation_digest_sha256": dag.equation_digest(target_equations),
        },
        "layer_zero": {
            "Hermite_profile": "H(rho)-H(0), H^(1..5)(rho)",
            "equation_count": len(hermite_equations),
            "equation_digest_sha256": dag.equation_digest(hermite_equations),
        },
        "combined_equation_count": len(equations),
        "all_equation_nodes": equations,
    }


def descent_eight_incidence_audit(
    dag: CircuitDAG, bands: dict[str, object]
) -> dict[str, object]:
    """Compile the five incidence rows and the relative quadratic relation."""

    K = bands["K_polynomials"]
    p_seven = K["P", 7]
    q_one = K["Q", 1]
    p_minus_one = K["P", -1]
    rho = RHO_FIELD.unit
    incidence = [
        evaluate_polynomial(dag, p_seven, rho, order) for order in range(3)
    ]
    incidence.append(evaluate_polynomial(dag, q_one, rho))
    a = dag.scale(
        RHO_FIELD.convert(sp.Rational(1, 6)),
        evaluate_polynomial(dag, p_seven, rho, 3),
    )
    y = dag.variable("descent8_y")
    scale = rho * (rho - RHO_FIELD.one) ** 6 * RHO_FIELD.convert(25**3)
    incidence.append(
        dag.add(
            evaluate_polynomial(dag, p_minus_one, rho),
            dag.scale(
                -scale,
                dag.multiply(y, dag.power(a, 2)),
            ),
        )
    )
    defect = dag.add(
        dag.scale(RHO_FIELD.convert(27), dag.power(y, 2)),
        dag.scale(RHO_FIELD.convert(-9), y),
        dag.one,
    )
    equations = incidence + [defect]
    return {
        "relative_base": "QQ(rho), rho^2-3*rho+1=0",
        "quartic_compositum_presentation": (
            "QQ(rho)[y]/(27*y^2-9*y+1)"
        ),
        "incidence_rows": [
            "K_P7(rho)",
            "K_P7'(rho)",
            "K_P7''(rho)",
            "K_Q1(rho)",
            (
                "K_P-1(rho)-25^3*rho*(rho-1)^6*y*"
                "(K_P7'''(rho)/6)^2"
            ),
        ],
        "incidence_equation_count": len(incidence),
        "defect_relation": "27*y^2-9*y+1",
        "equation_digest_sha256": dag.equation_digest(equations),
        "all_equation_nodes": equations,
    }


def squarefree_first_defect_ledger() -> dict[str, object]:
    """Route the squarefree carrier through every later defect position."""

    p_bands = {
        layer: make_band("P", 75, 3, layer) for layer in range(-75, 16)
    }
    q_bands = {
        layer: make_band("Q", 125, 5, layer) for layer in range(-125, 26)
    }
    audit = nonlinear_first_defect_audit(p_bands, q_bands)
    ledger = audit["first_defect_position_ledger"]
    later_rows = [
        row for row in ledger["exact_rows"] if int(row["spacing"]) >= 9
    ]
    if len(later_rows) != 82:
        raise AssertionError("the squarefree later-defect ledger changed")
    regime_counts = Counter(
        tuple(row["primitive_zero_multiples_before_target"])
        for row in later_rows
    )
    expected = {(2, 3): 3, (2,): 6, (): 73}
    if dict(regime_counts) != expected:
        raise AssertionError("the squarefree later-defect regimes changed")
    encoded = json.dumps(later_rows, sort_keys=True, separators=(",", ":"))
    return {
        "carrier": "R=(w^2-3*w+3)/25",
        "carrier_discriminant": "-3/625",
        "carrier_root_profile": "two simple roots; no movable double prime",
        "descent_eight_route": (
            "incompatible with the movable-double-root incidence prime"
        ),
        "later_spacing_interval": [9, 90],
        "later_spacing_count": len(later_rows),
        "regime_counts": {
            "spacings_9_to_11_with_multiples_2_3": regime_counts[(2, 3)],
            "spacings_12_to_17_with_multiple_2": regime_counts[(2,)],
            "spacings_18_to_90_with_no_pre_target_multiple": regime_counts[()],
        },
        "exact_rows": later_rows,
        "ledger_digest_sha256": hashlib.sha256(encoded.encode()).hexdigest(),
        "status": (
            "routed but not eliminated; each spacing still requires its own "
            "target-and-tail Fitting compiler"
        ),
    }


def reachable_circuit_summary(
    dag: CircuitDAG, roots: list[int]
) -> dict[str, object]:
    reachable: set[int] = set()
    stack = list(roots)
    while stack:
        node = stack.pop()
        if node in reachable:
            continue
        reachable.add(node)
        record = dag.nodes[node]
        if record[0] == "scale":
            stack.append(record[2])
        elif record[0] in ("add", "mul"):
            stack.extend((record[1], record[2]))
    operations = Counter(str(dag.nodes[node][0]) for node in reachable)
    variables = sorted(
        str(dag.nodes[node][1])
        for node in reachable
        if dag.nodes[node][0] == "var"
    )
    return {
        "reachable_node_count": len(reachable),
        "reachable_nodes_by_operation": dict(sorted(operations.items())),
        "reachable_variable_count": len(variables),
        "reachable_variable_name_digest_sha256": hashlib.sha256(
            "\n".join(variables).encode()
        ).hexdigest(),
        "maximum_total_degree": max(
            (dag.degrees[node] for node in roots), default=0
        ),
    }


def build_payload() -> dict[str, object]:
    dag = CircuitDAG()
    bands = build_band_circuits(dag)
    coupled = coupled_forcing_audit(dag, bands)
    coupled_nodes = coupled.pop("all_equation_nodes")
    final_functionals = target_and_layer_zero_audit(dag, bands)
    final_nodes = final_functionals.pop("all_equation_nodes")
    incidence = descent_eight_incidence_audit(dag, bands)
    incidence_nodes = incidence.pop("all_equation_nodes")
    all_equations = coupled_nodes + final_nodes + incidence_nodes
    if len(all_equations) != 366:
        raise AssertionError("the complete nonlinear equation ledger changed")

    return {
        "schema": "plane-jc.f2-75-125-nonlinear-forcing.v1",
        "status": "exact-arithmetic-circuit-presentation;ideal-test-open",
        "double_carrier": {
            "R": "(w-rho)^2/(25*(1-rho)^2)",
            "rho_polynomial": str(RHO_POLYNOMIAL.as_expr()),
            "quartic_primitive_polynomial": (
                "729*theta^4-4860*theta^3+10341*theta^2-"
                "7470*theta+1756"
            ),
        },
        "endpoint_circuits": {
            "pivot_count": len(bands["endpoint_solution_nodes"]),
            "upper_tangent_descents": [1, 11],
            "resonant_upper_centralizers": ["C0^4", "C0^3"],
            "H_at_zero_is_zero_after_full_band_specialization": (
                bands["H0"] == dag.zero
            ),
        },
        "coupled_Laurent_forcing": coupled,
        "final_thirteen_functionals": final_functionals,
        "descent_eight_incidence": incidence,
        "equation_ledger": {
            "coupled_full_Laurent": len(coupled_nodes),
            "target_and_layer_zero": len(final_nodes),
            "incidence_and_relative_field_relation": len(incidence_nodes),
            "total": len(all_equations),
            "all_equation_digest_sha256": dag.equation_digest(all_equations),
        },
        "arithmetic_circuit": reachable_circuit_summary(dag, all_equations),
        "squarefree_carrier_later_first_defect": (
            squarefree_first_defect_ledger()
        ),
        "decision": {
            "unit_ideal_obtained": False,
            "counterexample_obtained": False,
            "what_is_now_exact": (
                "the ten endpoint circuits, upper tangent substitutions, "
                "294 divisibility coordinates, 53 pinned quotient cokernel "
                "coordinates, final 7+6 functionals, and five incidence rows"
            ),
            "next_exact_operation": (
                "base-change the circuit ideal to each quartic embedding or "
                "a finite field of good reduction and run modular ideal/"
                "Jacobian tests before any characteristic-zero Groebner lift"
            ),
        },
        "reproduction_command": (
            ".venv/bin/python plane-jc/cas/"
            "compile_f2_75_125_nonlinear_forcing.py"
        ),
        "software": {"sympy": sp.__version__},
    }


def artifact_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=ARTIFACT)
    parser.add_argument("--refresh", action="store_true")
    args = parser.parse_args()

    payload = build_payload()
    artifact = args.artifact.resolve()
    if args.refresh:
        artifact.parent.mkdir(parents=True, exist_ok=True)
        artifact.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        try:
            display = artifact.relative_to(ROOT)
        except ValueError:
            display = artifact
        print(f"WROTE {display}")
    else:
        expected = json.loads(artifact.read_text())
        current_claim = {key: value for key, value in payload.items() if key != "software"}
        pinned_claim = {key: value for key, value in expected.items() if key != "software"}
        if current_claim != pinned_claim:
            raise AssertionError(
                "the pinned nonlinear-forcing artifact is stale; inspect "
                "the change before using --refresh"
            )

    print("F2_ENDPOINT_CIRCUIT_SUBSTITUTION_PASS")
    print("F2_UPPER_TANGENT_PROPAGATION_PASS")
    print("F2_COUPLED_DIVISIBILITY_294_PASS")
    print("F2_COUPLED_PINNED_COKERNEL_53_PASS")
    print("F2_FINAL_TARGET_HERMITE_7_PLUS_6_PASS")
    print("F2_DESCENT_EIGHT_INCIDENCE_PASS")
    print("F2_SQUAREFREE_LATER_DEFECT_9_TO_90_ROUTED")
    print("F2_NONLINEAR_IDEAL_TEST_OPEN")
    print(f"F2_NONLINEAR_FORCING_ARTIFACT_SHA256={artifact_sha256(artifact)}")


if __name__ == "__main__":
    main()
