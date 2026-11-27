"""Exact degree-six Keller chart checks for diagonal one-parameter actions."""

from collections import defaultdict
from itertools import permutations, product
from math import gcd

from sympy import Poly, symbols

from verify_sextic_graded_112 import bracket, verify_unit_ideal


DEGREE = 6
LINEAR = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
BALANCED = (1, -1, -1)
KNOWN_EXAMPLE_WEIGHT = (1, -1, -2)
FIXED_LINE_WEIGHT = (1, 0, -1)
NONLINEAR = tuple(
    exponent for exponent in product(range(DEGREE + 1), repeat=3)
    if 2 <= sum(exponent) <= DEGREE
)
RESONANCES = tuple(
    tuple(
        (exponent, tuple(exponent[j] - LINEAR[i][j] for j in range(3)))
        for exponent in NONLINEAR
    )
    for i in range(3)
)
DETERMINANT_SIGNS = tuple(
    (order, (-1)**sum(order[i] > order[j] for i in range(3) for j in range(i + 1, 3)))
    for order in permutations(range(3))
)


def cross(left, right):
    a, b, c = left
    x, y, z = right
    return b*z - c*y, c*x - a*z, a*y - b*x


def candidate_weights():
    other_resonances = {difference for terms in RESONANCES for _, difference in terms}
    positive_resonances = {
        difference for exponent, difference in RESONANCES[0]
        if exponent[0] >= 2 and exponent[1] + exponent[2] >= 1
    }
    weights = set()
    for first in positive_resonances:
        for second in other_resonances:
            a, b, c = cross(first, second)
            if a < 0:
                a, b, c = -a, -b, -c
            if a > 0 >= b >= c and c < 0:
                common_factor = gcd(a, b, c)
                weights.add((a // common_factor, b // common_factor, c // common_factor))
    return weights


def supports_for(weight):
    supports = tuple(
        tuple(exponent for exponent, difference in terms
              if sum(w*d for w, d in zip(weight, difference)) == 0)
        for terms in RESONANCES
    )
    if weight[1] != 0:
        return supports
    # The restriction and first derivative on the fixed line are the identity.
    return tuple(
        tuple(exponent for exponent in terms if exponent[0] + exponent[2] > int(i != 1))
        for i, terms in enumerate(supports)
    )


def determinant_coefficients(components):
    gradients = [
        [[(exponent, parameter, exponent[j]) for exponent, parameter in terms if exponent[j]]
         for j in range(3)]
        for terms in components
    ]
    residual = defaultdict(lambda: defaultdict(int))
    for order, sign in DETERMINANT_SIGNS:
        for selected in product(*(gradients[i][j] for i, j in enumerate(order))):
            # Every column occurs once, so every exponent is lowered by one.
            exponent = tuple(sum(term[0][j] for term in selected) - 1 for j in range(3))
            if exponent == (0, 0, 0):
                continue
            scalar = sign
            for _, _, derivative in selected:
                scalar *= derivative
            monomial = "*".join(parameter for _, parameter, _ in selected if parameter != "1")
            residual[exponent][monomial] += scalar
    equations = []
    for exponent in sorted(residual):
        terms = residual[exponent]
        summands = [
            f"{coefficient:+d}{'*' + monomial if monomial else ''}"
            for monomial, coefficient in sorted(terms.items()) if coefficient
        ]
        if summands:
            equations.append("".join(summands).lstrip("+"))
    return equations


def check_chart(supports, chart):
    parameters = []
    components = []
    for i, terms in enumerate(supports):
        component = [(LINEAR[i], "1")]
        for j, exponent in enumerate(terms):
            if i == 0 and j < chart:
                continue
            if i == 0 and j == chart:
                component.append((exponent, "1"))
            else:
                parameter = f"a{i}b{j}"
                component.append((exponent, parameter))
                parameters.append(parameter)
        components.append(component)
    equations = determinant_coefficients(components)
    if not equations:
        raise RuntimeError("The normalized chart has an identically zero determinant residual")
    verify_unit_ideal(equations, parameters or ["unused"])


def check_balanced():
    u, v, c = symbols("u v c")
    first = symbols("a1:8")
    second = symbols("b1:8")
    monomials = (u**2, u*v, v**2, u**3, u**2*v, u*v**2, v**3)
    B = u + sum(coeff * term for coeff, term in zip(first, monomials))
    C = v + sum(coeff * term for coeff, term in zip(second, monomials))
    choices = (1 + u + v**2 + c*u**2, 1 + u + u*v, 1 + u + c*u**2, 1 + u*v, 1 + u**2)
    for multiplier in choices:
        residual = (
            multiplier * bracket(B, C, u, v) + B * bracket(multiplier, C, u, v)
            + C * bracket(B, multiplier, u, v) - 1
        )
        parameters = (*first, *second, c) if c in multiplier.free_symbols else (*first, *second)
        verify_unit_ideal(Poly(residual, u, v).coeffs(), parameters)
    return len(choices)


def main():
    control = [[(exponent, "1")] for exponent in LINEAR]
    if determinant_coefficients(control):
        raise RuntimeError("The identity map gave a nonzero determinant residual")
    control[0].append(((2, 0, 1), "1"))
    if determinant_coefficients(control) != ["2"]:
        raise RuntimeError("The determinant residual for x -> x + x^2*z is incorrect")

    counts = {False: 0, True: 0}
    for weight in sorted(candidate_weights() - {BALANCED, KNOWN_EXAMPLE_WEIGHT}):
        supports = supports_for(weight)
        for chart in range(len(supports[0])):
            chart_supports = supports
            if weight == FIXED_LINE_WEIGHT and chart == 0:
                # Axis translation and axis exchange leave this chart only
                # when both coefficients of x*z in the multipliers are constant.
                chart_supports = (
                    tuple(term for term in supports[0] if term[2] != 1 or term[1] == 0),
                    supports[1],
                    tuple(term for term in supports[2] if term[0] != 1 or term[1] == 0),
                )
            try:
                check_chart(chart_supports, chart)
            except Exception as error:
                raise RuntimeError(f"Failed weight {weight}, chart {chart}") from error
            counts[weight[1] == 0] += 1
    balanced_count = check_balanced()
    if (counts[False], counts[True], balanced_count) != (349, 20, 5):
        raise RuntimeError("The exact chart census changed")
    print(f"{counts[False]} nonzero-weight charts, {counts[True]} fixed-line charts, "
          f"and {balanced_count} balanced charts: unit ideals over Q.")
    print("The four charts for (1,-1,-2) are covered by verify_sextic_graded_112.py.")


if __name__ == "__main__":
    main()
