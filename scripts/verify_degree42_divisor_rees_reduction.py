#!/usr/bin/env python3
"""Certify the two-normal-variable Rees reduction on ``w0=0``."""

from __future__ import annotations

import shutil
import subprocess

import sympy as sp

from verify_degree42_transported_27_normal_jets import (
    serialize,
    transformed_problem,
)


def main() -> None:
    normals, bases, residuals, _defect = transformed_problem()
    variables = normals + bases
    e1, e2, translation, w0, w1, _w2 = bases
    zero_normals = ",".join(f"{variable},0" for variable in normals)
    singular = shutil.which("Singular")
    assert singular is not None, "Singular is required for the certificate"

    program = f"""
ring q=0,({",".join(map(str, variables))}),(dp(5),dp(6));
ideal I={",".join(serialize(item) for item in residuals)};
matrix J=jacob(I);
matrix C[nrows(J)][5];
int i;
int j;
for (i=1;i<=nrows(J);i++)
{{
  for (j=1;j<=5;j++)
  {{
    C[i,j]=subst(J[i,j],{zero_normals},{w0},0);
  }}
}}
int tangent_two=1;
for (i=1;i<=nrows(C);i++)
{{
  if (C[i,1]!=0 or C[i,2]!=0) {{ tangent_two=0; }}
}}
int triangular=(
  C[5,3]==-1 and
  C[11,3]==0 and C[11,4]==-1 and
  C[17,3]==0 and C[17,4]==0 and C[17,5]==-1
);
int terminal_zero=1;
for (i=18;i<=19;i++)
{{
  for (j=1;j<=5;j++)
  {{
    if (C[i,j]!=0) {{ terminal_zero=0; }}
  }}
}}
poly q11=subst(diff(diff(I[18],{normals[0]}),{normals[0]}),
  {zero_normals},{w0},0)/2;
poly q12=subst(diff(diff(I[18],{normals[0]}),{normals[1]}),
  {zero_normals},{w0},0);
poly q22=subst(diff(diff(I[18],{normals[1]}),{normals[1]}),
  {zero_normals},{w0},0)/2;
poly r11=subst(diff(diff(I[19],{normals[0]}),{normals[0]}),
  {zero_normals},{w0},0)/2;
poly r12=subst(diff(diff(I[19],{normals[0]}),{normals[1]}),
  {zero_normals},{w0},0);
poly r22=subst(diff(diff(I[19],{normals[1]}),{normals[1]}),
  {zero_normals},{w0},0)/2;
int quadrics=(
  q11==-9/2*{e2}*{w1} and
  q12==(9*{e2}^2+3/2*{e1})*{w1} and
  q22==(-9/2*{e2}^3+15/4*{e1}*{e2}+3/4*{translation})*{w1} and
  r11==-3/4*{w1} and
  r12==3/2*{e2}*{w1} and
  r22==(-3/4*{e2}^2+3/4*{e1})*{w1}
);
print("DEGREE42_DIVISOR_REES");
print(tangent_two);
print(triangular);
print(terminal_zero);
print(quadrics);
"""
    process = subprocess.run(
        [singular, "-q"],
        input=program,
        text=True,
        capture_output=True,
        timeout=300,
        check=True,
    )
    compact = process.stdout.split()
    marker = compact.index("DEGREE42_DIVISOR_REES")
    assert compact[marker + 1 : marker + 5] == ["1", "1", "1", "1"], (
        process.stdout + process.stderr
    )

    u, v = sp.symbols("u v")
    q = (
        -sp.Rational(9, 2) * e2 * u**2
        + (9 * e2**2 + sp.Rational(3, 2) * e1) * u * v
        + (
            -sp.Rational(9, 2) * e2**3
            + sp.Rational(15, 4) * e1 * e2
            + sp.Rational(3, 4) * translation
        )
        * v**2
    )
    r = (
        -sp.Rational(3, 4) * u**2
        + sp.Rational(3, 2) * e2 * u * v
        + (-sp.Rational(3, 4) * e2**2 + sp.Rational(3, 4) * e1)
        * v**2
    )
    discriminant = (translation + e1 * e2) ** 2 - 4 * e1**3
    assert sp.factor(sp.resultant(q.subs(v, 1), r.subs(v, 1), u)) == (
        sp.Rational(81, 256) * discriminant
    )
    print(
        "PASS: on w0=0 three unit pivots leave two normal variables; "
        "their quadratic Rees resultant is (81/256)*w1^4*Delta"
    )


if __name__ == "__main__":
    main()
