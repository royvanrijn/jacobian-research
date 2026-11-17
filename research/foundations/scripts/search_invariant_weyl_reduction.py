#!/usr/bin/env python3
"""Bounded search for an invariant position-type Weyl pair in the A_3 lift.

Let Phi_F be the inverse-Jacobian Weyl endomorphism attached to the
foundational Keller map F.  Suppose an invariant copy of A_1 has generators

    h(X) and L = sum_i v_i(X) D_i + c(X),    [L,h] = 1.

Normal ordering in (h,L) and the differential-order filtration imply

    k<h,L> intersect k[X] = k[h].

Therefore Phi_F(k<h,L>) contained in k<h,L> forces the necessary polynomial
semiconjugacy h(F)=g(h).  In characteristic zero this implies

    dh wedge d(h(F)) = 0.

The script forms every coefficient equation in that wedge for a completely
general h of total degree at most ``--max-degree``.  The equations are
homogeneous quadrics in the nonconstant coefficients of h.  A Groebner basis
over F_32003 proves that their projective zero locus is empty.  Since the
equations are integral and homogeneous, any characteristic-zero projective
solution would have a nonzero reduction over the algebraic closure of
F_32003.  Thus the modular calculation is a characteristic-zero no-go
certificate for this bounded ansatz.

This does not search Weyl pairs for which both generators have positive
differential order, nor higher-degree h.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess

import sympy as sp


PRIME = 32003


def monomials_through_degree(variables: tuple[sp.Symbol, ...], degree: int):
    """Return all nonconstant monomials in three variables through degree."""

    x, y, z = variables
    result = []
    for total in range(1, degree + 1):
        for x_degree in range(total + 1):
            for y_degree in range(total - x_degree + 1):
                z_degree = total - x_degree - y_degree
                result.append(x**x_degree * y**y_degree * z**z_degree)
    return tuple(result)


def foundational_map(x: sp.Symbol, y: sp.Symbol, z: sp.Symbol):
    u = 1 + x * y
    return sp.Matrix(
        [
            u**3 * z + y**2 * u * (4 + 3 * x * y),
            y + 3 * x * u**2 * z + 3 * x * y**2 * (4 + 3 * x * y),
            2 * x - 3 * x**2 * y - x**3 * z,
        ]
    )


def coefficient_equations(max_degree: int):
    x, y, z = sp.symbols("x y z")
    variables = (x, y, z)
    monomials = monomials_through_degree(variables, max_degree)
    coefficients = sp.symbols(f"a0:{len(monomials)}")
    h = sum(
        coefficient * monomial
        for coefficient, monomial in zip(coefficients, monomials, strict=True)
    )
    F = foundational_map(x, y, z)
    h_after_F = sp.expand(
        h.subs(dict(zip(variables, F, strict=True)), simultaneous=True)
    )
    dh = sp.Matrix([sp.diff(h, variable) for variable in variables])
    dh_after_F = sp.Matrix(
        [sp.diff(h_after_F, variable) for variable in variables]
    )

    equations = []
    for expression in dh.cross(dh_after_F):
        equations.extend(sp.Poly(expression, *variables).coeffs())

    # SymPy can return the same coefficient from different wedge components.
    unique = []
    seen = set()
    for equation in equations:
        equation = sp.expand(equation)
        key = sp.srepr(equation)
        if equation != 0 and key not in seen:
            seen.add(key)
            unique.append(equation)
    return variables, coefficients, monomials, F, tuple(unique)


def singular_polynomial(expression: sp.Expr) -> str:
    return str(expression).replace("**", "^")


def modular_projective_certificate(
    coefficients: tuple[sp.Symbol, ...], equations: tuple[sp.Expr, ...]
):
    singular = shutil.which("Singular")
    if singular is None:
        raise RuntimeError(
            "Singular is required for the modular Groebner certificate"
        )

    names = ",".join(map(str, coefficients))
    generators = ",".join(singular_polynomial(equation) for equation in equations)
    program = [
        f"ring r={PRIME},({names}),dp;",
        f"ideal I={generators};",
        "ideal G=std(I);",
        'if (dim(G)==0) { print("PASS_DIM_ZERO"); } else { print("FAIL_DIM"); };',
        'print("VDIM");',
        "vdim(G);",
        'print("GB_SIZE");',
        "size(G);",
    ]
    # Nilpotence of every homogeneous coordinate proves that the radical is
    # the irrelevant maximal ideal.  Fourth powers cover the checked bounds.
    for coefficient in coefficients:
        program.append(
            f'if (reduce({coefficient}^4,G)==0) '
            f'{{ print("PASS_NILPOTENT_{coefficient}"); }} '
            f'else {{ print("FAIL_NILPOTENT_{coefficient}"); }};'
        )

    result = subprocess.run(
        [singular, "-q"],
        input="\n".join(program) + "\n",
        text=True,
        capture_output=True,
        check=True,
        timeout=240,
    )
    output = result.stdout
    assert "PASS_DIM_ZERO" in output and "FAIL_DIM" not in output
    assert "FAIL_NILPOTENT" not in output
    for coefficient in coefficients:
        assert f"PASS_NILPOTENT_{coefficient}" in output
    return output


def check_canonical_image_pairs(variables, F):
    """The three most immediate pairs <F_i,delta_i> are not invariant."""

    points = ((0, 0, 0), (0, 0, 1), (0, 0, 0))
    expected_witnesses = ((0, 2, 0), (0, 0, -6), (0, -4, 0))
    witnesses = []
    substitution_F = dict(zip(variables, F, strict=True))
    for coordinate, point, expected in zip(
        F, points, expected_witnesses, strict=True
    ):
        iterate = sp.expand(coordinate.subs(substitution_F, simultaneous=True))
        first_gradient = sp.Matrix(
            [sp.diff(coordinate, variable) for variable in variables]
        )
        iterate_gradient = sp.Matrix(
            [sp.diff(iterate, variable) for variable in variables]
        )
        witness = tuple(
            value.subs(dict(zip(variables, point, strict=True)))
            for value in first_gradient.cross(iterate_gradient)
        )
        assert witness == expected
        witnesses.append((point, witness))
    return tuple(witnesses)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-degree", type=int, default=4, choices=range(1, 5))
    args = parser.parse_args()

    variables, coefficients, monomials, F, equations = coefficient_equations(
        args.max_degree
    )
    output = modular_projective_certificate(coefficients, equations)
    witnesses = check_canonical_image_pairs(variables, F)

    lines = [line for line in output.splitlines() if line.strip()]
    vdim = lines[lines.index("VDIM") + 1]
    basis_size = lines[lines.index("GB_SIZE") + 1]
    print(
        f"PASS: {len(equations)} wedge equations in {len(coefficients)} "
        f"coefficients have empty projective zero locus over F_{PRIME}"
    )
    print(
        f"PASS: quotient vector-space dimension={vdim}; "
        f"Groebner-basis size={basis_size}; every coefficient is nilpotent"
    )
    print(
        "PASS: no nonconstant h of total degree <= "
        f"{args.max_degree} satisfies dh wedge d(h o F)=0 in characteristic zero"
    )
    for index, (point, witness) in enumerate(witnesses, start=1):
        print(
            f"PASS: F_{index} o F is not in k[F_{index}] "
            f"(witness {witness} at {point})"
        )
    print(
        "SCOPE: excludes invariant position/first-order Weyl pairs only; "
        "positive-order pairs and higher-degree positions remain open"
    )


if __name__ == "__main__":
    main()
