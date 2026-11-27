#!/usr/bin/env python3
"""Replay exact unit ideals for the ordinary sextic (1,-1,-2) Keller charts."""

from subprocess import run

from sympy import Poly, diff, expand, symbols


def bracket(left, right, u, v):
    return diff(left, u) * diff(right, v) - diff(left, v) * diff(right, u)


def keller_polynomial(multiplier, first_mate, second_mate, u, v):
    return (
        multiplier * bracket(first_mate, second_mate, u, v)
        + 2 * second_mate * bracket(first_mate, multiplier, u, v)
        + first_mate * bracket(multiplier, second_mate, u, v)
    )


def verify_unit_ideal(coefficients, variables):
    variable_list = ",".join(map(str, variables))
    equation_list = ",\n".join(str(coefficient).replace("**", "^") for coefficient in coefficients)
    program = f"""
ring R=0,({variable_list}),dp;
ideal I={equation_list};
option(redSB);
ideal G=slimgb(I);
if (size(G)!=1 || G[1]!=1) {{ ERROR("The chart did not give a unit ideal"); }}
print("Exact unit ideal verified");
quit;
"""
    result = run(
        ["Singular", "-q"], input=program, text=True, capture_output=True, timeout=60, check=True
    )
    if result.stdout.strip() != "Exact unit ideal verified" or result.stderr:
        raise RuntimeError(f"Invalid Singular result:\n{result.stdout}{result.stderr}")


def main():
    u, v = symbols("u v")
    multipliers = symbols("l1 l2 l3 l4")
    first_mate_coefficients = symbols("b2 b3 b4 b5 b6 b7")
    second_mate_coefficients = symbols("a1 a2 a3 a4 a5 a6 a7")
    parameters = multipliers + first_mate_coefficients + second_mate_coefficients

    l1, l2, l3, l4 = multipliers
    b2, b3, b4, b5, b6, b7 = first_mate_coefficients
    a1, a2, a3, a4, a5, a6, a7 = second_mate_coefficients
    multiplier = 1 + l1*u + l2*u**2 + l3*v + l4*u*v
    first_mate = u + b2*u**2 + b3*u**3 + (b4 + b5*u + b6*u**2)*v + b7*v**2
    second_mate = v + a1*u**2 + a2*u**3 + a3*u**4 + (a4*u + a5*u**2)*v + (a6 + a7*u)*v**2
    residual = keller_polynomial(multiplier, first_mate, second_mate, u, v) - 1

    if expand(keller_polynomial(1, u, v + u**2, u, v)) != 1:
        raise RuntimeError("The triangular control failed")
    seed = v*(1 + u)**2 + u**2*(4 + 3*u)
    known_multiplier = (2 - 3*u - v) / 2
    known_first_mate = u + 3*seed
    known_second_mate = (1 + u)*seed
    if expand(keller_polynomial(known_multiplier, known_first_mate, known_second_mate, u, v)) != 1:
        raise RuntimeError("The degree-seven control failed")
    truncated = known_second_mate - u**3*v
    if expand(keller_polynomial(known_multiplier, known_first_mate, truncated, u, v)) == 1:
        raise RuntimeError("The determinant failed to detect truncation")

    earlier = {}
    for chart in (l4, l2, l3, l1):
        substitutions = {**earlier, chart: 1}
        variables = [parameter for parameter in parameters if parameter not in substitutions]
        coefficients = Poly(residual.subs(substitutions), u, v).coeffs()
        verify_unit_ideal(coefficients, variables)
        prefix = ", ".join(f"{parameter} = 0" for parameter in earlier)
        print(f"{prefix + ', ' if prefix else ''}{chart} = 1: exact unit ideal", flush=True)
        earlier[chart] = 0


if __name__ == "__main__":
    main()
