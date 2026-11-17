#!/usr/bin/env python3
"""Search for a one-element depth certificate for the degree-42 residual.

For the residual ideal ``I`` and support ideal

    k = (w0, w1*w2, A*B*w2),

an equality ``I:f = I`` for one ``f`` in ``k`` proves that ``f`` is a
non-zero-divisor on ``R/I``.  This is equivalent to positive k-grade and
therefore implies ``I:k^infinity = I``.
"""

from __future__ import annotations

import argparse
import os
import shutil
import signal
import subprocess

import sympy as sp

from verify_degree42_transported_27_normal_jets import (
    serialize,
    transformed_problem,
)


def support_factors(
    e1: sp.Symbol,
    e2: sp.Symbol,
    translation: sp.Symbol,
) -> tuple[sp.Expr, sp.Expr]:
    """Return the two exceptional Kuranishi factors A and B."""

    factor_a = (
        -e1**2 * e2**2
        + e2**3 * translation
        + 4 * e1**3
        - 6 * e1 * e2 * translation
    )
    factor_b = (
        e1**4 * e2**4
        - 2 * e1**2 * e2**5 * translation
        + e2**6 * translation**2
        - 8 * e1**5 * e2**2
        + 20 * e1**3 * e2**3 * translation
        - 12 * e1 * e2**4 * translation**2
        + 16 * e1**6
        - 48 * e1**4 * e2 * translation
        + 27 * e1**2 * e2**2 * translation**2
        + 9 * e2**3 * translation**3
        + 36 * e1**3 * translation**2
        - 54 * e1 * e2 * translation**3
        + 27 * translation**4
    )
    return factor_a, factor_b


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--coefficients",
        default="1,1,1",
        help="comma-separated integer coefficients alpha,beta,gamma",
    )
    parser.add_argument("--prime", type=int, default=0)
    parser.add_argument("--timeout", type=int, default=900)
    parser.add_argument(
        "--algorithm",
        choices=("std", "slimgb"),
        default="slimgb",
    )
    parser.add_argument(
        "--order",
        default=None,
        help=(
            "Singular monomial order; defaults to the normal/base block "
            "order, or (dp(5),dp(6)) with --full-presentation"
        ),
    )
    parser.add_argument(
        "--full-presentation",
        action="store_true",
        help="skip the three exact unit-pivot eliminations (diagnostic only)",
    )
    parser.add_argument(
        "--method",
        choices=("pivots", "height", "colon"),
        default="height",
        help=(
            "pivots certifies only the exact two-normal-variable reduction; "
            "height tests a codimension-two perfect presentation; colon "
            "performs the direct ideal quotient"
        ),
    )
    args = parser.parse_args()

    coefficients = tuple(
        int(value.strip()) for value in args.coefficients.split(",")
    )
    assert len(coefficients) == 3
    assert any(coefficients)
    assert not (
        args.full_presentation and args.method in {"pivots", "height"}
    )

    normals, bases, residuals, _defect = transformed_problem()
    e1, e2, translation, w0, w1, w2 = bases
    factor_a, factor_b = support_factors(e1, e2, translation)
    alpha, beta, gamma = coefficients
    regular_candidate = (
        alpha * w0
        + beta * w1 * w2
        + gamma * factor_a * factor_b * w2
    )

    singular = shutil.which("Singular")
    assert singular is not None, "Singular is required for the certificate"
    variables = normals + bases
    characteristic = args.prime if args.prime else 0
    basis = args.algorithm
    if args.full_presentation:
        monomial_order = args.order or "(dp(5),dp(6))"
        program = f"""
ring q={characteristic},({",".join(map(str, variables))}),{monomial_order};
ideal I={",".join(serialize(item) for item in residuals)};
poly f={serialize(regular_candidate)};
print("DEPTH_STAGE_PRESENTATION");
ideal G={basis}(I);
print("DEPTH_STAGE_BASIS");
ideal Q={basis}(quotient(G,ideal(f)));
print("DEPTH_STAGE_COLON");
int pivots=1;
int generated=1;
int expected_dimensions=1;
int forward=(reduce(Q,G)==0);
int reverse=(reduce(G,Q)==0);
"""
    else:
        monomial_order = args.order or "(dp(2),dp(6))"
        normal_images = ",".join(
            [str(normals[0]), str(normals[1]), "0", "0", "0"]
        )
        mapped_residuals = ",".join(
            f"M[1,{index}]" for index in range(1, len(residuals) + 1)
        )
        program = f"""
ring source={characteristic},({",".join(map(str, variables))}),(dp(5),dp(6));
ideal I={",".join(serialize(item) for item in residuals)};
matrix residual_matrix[1][{len(residuals)}]=I;
poly p3=subst(I[5],{normals[2]},0);
int pivot3=(I[5]==-{normals[2]}+p3);
poly pivot4_poly=subst(I[11],{normals[2]},p3);
poly p4=subst(pivot4_poly,{normals[3]},0);
int pivot4=(pivot4_poly==-{normals[3]}+p4);
poly pivot5_poly=subst(
  subst(I[17],{normals[2]},p3),{normals[3]},p4);
poly p5=subst(pivot5_poly,{normals[4]},0);
int pivot5=(pivot5_poly==-{normals[4]}+p5);
matrix M3=subst(residual_matrix,{normals[2]},p3);
matrix M4=subst(M3,{normals[3]},p4);
matrix M5=subst(M4,{normals[4]},p5);
poly pivot_product=pivot3*pivot4*pivot5;

ring q={characteristic},({normals[0]},{normals[1]},
  {",".join(map(str, bases))}),{monomial_order};
map reduce_normals=source,{normal_images},{",".join(map(str, bases))};
matrix M=reduce_normals(M5);
ideal I={mapped_residuals};
ideal terminal=M[1,18],M[1,19];
poly f={serialize(regular_candidate)};
print("DEPTH_STAGE_PRESENTATION");
int pivots=(imap(source,pivot_product)==1);
"""
        if args.method == "pivots":
            program += """
ideal G=ideal(0);
ideal Q=G;
int generated=1;
int expected_dimensions=1;
int forward=1;
int reverse=1;
"""
        elif args.method == "height":
            program += f"""
ideal G={basis}(terminal);
int terminal_generated=(reduce(I,G)==0);
if (!terminal_generated)
{{
  G={basis}(I);
}}
print("DEPTH_STAGE_BASIS");
int generated=(reduce(I,G)==0);
resolution L=mres(G,0);
ideal H={basis}(I+ideal(f));
print("DEPTH_STAGE_HEIGHT");
int expected_dimensions=(
  dim(G)==6 and dim(H)==5 and size(L)==3
);
ideal Q=G;
int forward=1;
int reverse=1;
"""
        else:
            program += f"""
ideal G={basis}(I);
print("DEPTH_STAGE_BASIS");
ideal Q={basis}(quotient(G,ideal(f)));
print("DEPTH_STAGE_COLON");
int generated=1;
int expected_dimensions=1;
int forward=(reduce(Q,G)==0);
int reverse=(reduce(G,Q)==0);
"""
    program += """
print("DEGREE42_DEPTH_CERTIFICATE");
print(pivots);
print(generated);
print(expected_dimensions);
print(forward);
print(reverse);
print(size(G));
print(size(Q));
"""
    process = subprocess.Popen(
        [singular, "-q"],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        start_new_session=True,
    )
    try:
        stdout, stderr = process.communicate(program, timeout=args.timeout)
    except subprocess.TimeoutExpired as error:
        os.killpg(process.pid, signal.SIGTERM)
        try:
            process.wait(timeout=5)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGKILL)
            process.wait()
        partial = error.stdout or ""
        if isinstance(partial, bytes):
            partial = partial.decode(errors="replace")
        raise TimeoutError(
            f"Singular did not finish within {args.timeout}s"
            + (f"; partial output:\n{partial}" if partial else "")
        ) from error
    if process.returncode:
        raise RuntimeError(stdout + stderr)
    compact = stdout.split()
    if "DEGREE42_DEPTH_CERTIFICATE" not in compact:
        raise RuntimeError(stdout + stderr)
    marker = compact.index("DEGREE42_DEPTH_CERTIFICATE")
    (
        pivots,
        generated,
        expected_dimensions,
        forward,
        reverse,
    ) = compact[marker + 1 : marker + 6]
    assert (
        pivots
        == generated
        == expected_dimensions
        == forward
        == reverse
        == "1"
    ), stdout + stderr
    field = f"GF({characteristic})" if characteristic else "QQ"
    if args.method == "pivots":
        print(
            "PASS: residuals 5, 11, and 17 are global unit pivots; "
            "the degree-42 residual algebra has an exact presentation "
            "with only two normal variables"
        )
    else:
        reason = (
            "after three unit pivots the reduced ideal is codimension-two "
            "perfect, and adjoining f raises its height by one"
            if args.method == "height"
            else "I:f=I"
        )
        print(
            "PASS: f="
            f"{alpha}*w0+{beta}*w1*w2+{gamma}*A*B*w2 is regular on R/I "
            f"over {field}; {reason} (basis sizes "
            f"{compact[marker + 6]}/{compact[marker + 7]})"
        )


if __name__ == "__main__":
    main()
