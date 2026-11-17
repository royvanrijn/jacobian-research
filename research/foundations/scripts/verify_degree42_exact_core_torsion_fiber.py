#!/usr/bin/env python3
"""Verify exact, untruncated support torsion on a degree-42 base fiber.

This checker is deliberately weaker than a generic saturation computation.
It specializes ``(e1,e2,t)`` to rational values, performs the three global
unit-pivot eliminations, and works in the full polynomial ring

    QQ[u,v,w0,w1,w2]

without imposing a normal-jet cutoff.  On the default fiber it verifies

    h = v*w0^2*(e2*u - (t + e1*e2)*v)

is nonzero modulo the residual ideal, while

    w0*h, w0^2*h, w0*w2*h, w2^2*h

vanish.  Hence ``h`` is killed by ``(w0,w2)^2`` on that exact fiber.
The calculation does not by itself imply generic torsion, since saturation
need not commute with specialization.
"""

from __future__ import annotations

import argparse
import os
import shutil
import signal
import subprocess

import sympy as sp

from verify_degree42_depth_certificate import support_factors
from verify_degree42_transported_27_normal_jets import (
    serialize,
    transformed_problem,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--base-values",
        default="1,2,3",
        help="comma-separated rational values of e1,e2,t",
    )
    parser.add_argument("--timeout", type=int, default=600)
    args = parser.parse_args()

    values = tuple(
        sp.Rational(value.strip()) for value in args.base_values.split(",")
    )
    assert len(values) == 3

    normals, bases, residuals, _defect = transformed_problem()
    variables = normals + bases
    e1, e2, translation, w0, w1, w2 = bases
    substitution = dict(zip((e1, e2, translation), values, strict=True))
    factor_a, factor_b = support_factors(e1, e2, translation)
    specialized_a = sp.factor(factor_a.subs(substitution))
    specialized_b = sp.factor(factor_b.subs(substitution))
    assert specialized_a != 0 and specialized_b != 0, (
        "the selected fiber lies on A*B=0, where k does not reduce to "
        "(w0,w2)"
    )

    program = f"""
ring source=0,({",".join(map(str, variables))}),(dp(5),dp(6));
ideal I={",".join(serialize(item) for item in residuals)};
poly p3=subst(I[5],{normals[2]},0);
ideal I3=subst(I,{normals[2]},p3);
poly p4=subst(I3[11],{normals[3]},0);
ideal I4=subst(I3,{normals[3]},p4);
poly p5=subst(I4[17],{normals[4]},0);
ideal I5=subst(I4,{normals[4]},p5);

ring q=0,(u,v,{w0},{w1},{w2},H),dp;
map phi=source,u,v,0,0,0,{values[0]},{values[1]},{values[2]},
  {w0},{w1},{w2};
ideal Core=phi(I5);
ideal HomogeneousCore=homog(Core,H);
ideal HomogeneousBasis=groebner(HomogeneousCore);
ideal G=std(subst(HomogeneousBasis,H,1));
poly h=v*{w0}^2*(
  {values[1]}*u-({values[2]}+{values[0]}*{values[1]})*v
);

int h_nonzero=(reduce(h,G)!=0);
int w0_kills=(reduce({w0}*h,G)==0);
int w2_does_not_kill=(reduce({w2}*h,G)!=0);
int square_kills=(
  reduce({w0}^2*h,G)==0
  and reduce({w0}*{w2}*h,G)==0
  and reduce({w2}^2*h,G)==0
);

print("DEGREE42_EXACT_CORE_TORSION_FIBER");
print(size(HomogeneousBasis));
print(h_nonzero);
print(w0_kills);
print(w2_does_not_kill);
print(square_kills);
print(reduce(h,G));
"""

    singular = shutil.which("Singular")
    assert singular is not None, "Singular is required"
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
        raise TimeoutError(
            f"Singular did not finish within {args.timeout}s"
        ) from error
    if process.returncode or "? error occurred" in stdout:
        raise RuntimeError(stdout + stderr)

    marker = "DEGREE42_EXACT_CORE_TORSION_FIBER"
    assert marker in stdout, stdout + stderr
    result = stdout[stdout.index(marker) :].splitlines()
    assert result[2] == result[3] == result[5] == "1", stdout + stderr
    print(
        "PASS: exact untruncated fiber "
        f"(e1,e2,t)=({values[0]},{values[1]},{values[2]}) has "
        "(w0,w2)^2-torsion represented by "
        "v*w0^2*(e2*u-(t+e1*e2)*v) "
        f"(basis size {result[1]})"
    )


if __name__ == "__main__":
    main()
