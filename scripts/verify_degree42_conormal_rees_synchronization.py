#!/usr/bin/env python3
"""Certify the conceptual conormal/Rees form of degree-42 synchronization.

On the transported ``{2,7}`` power chart, let ``M`` be the ideal of the
five refinement-normal coordinates and let ``I`` be the pulled-back Hessian
residual ideal.  This checker proves that the residual conormal map

    I -> M/M^2

is surjective precisely off ``w0=0``.  More exactly, its maximal-minor ideal
is ``(w0^2)``; on ``w0=0`` its rank is exactly three.  Thus over ``D(w0)``
the equality ``M=I+M^2`` and complete Nakayama give ``M^ = I^``.  Since the
synchronization defect lies in ``M``, it vanishes to every Rees order there.
"""

from __future__ import annotations

import shutil
import subprocess

from verify_degree42_transported_27_normal_jets import (
    serialize,
    transformed_problem,
)


def main() -> None:
    normals, bases, residuals, defect = transformed_problem()
    variables = normals + bases
    w0 = bases[3]
    zero_substitutions = ",".join(f"{variable},0" for variable in normals)
    singular = shutil.which("Singular")
    assert singular is not None, "Singular is required for the certificate"

    program = f"""
ring q=0,({",".join(map(str, variables))}),(dp(5),dp(6));
ideal I={",".join(serialize(item) for item in residuals)};
poly delta={serialize(defect)};
ideal M={",".join(map(str, normals))};
ideal GM=std(M);
int onzero=1;
int i;
for (i=1;i<=size(I);i++)
{{
  if (reduce(I[i],GM)!=0) {{ onzero=0; }}
}}
if (reduce(delta,GM)!=0) {{ onzero=0; }}

matrix J=jacob(I);
matrix C[nrows(J)][5];
int j;
for (i=1;i<=nrows(J);i++)
{{
  for (j=1;j<=5;j++)
  {{
    C[i,j]=subst(J[i,j],{zero_substitutions});
  }}
}}

ideal G5=std(minor(C,5));
ideal D4=minor(C,4);
D4=D4,{w0};
ideal G4=std(D4);
ideal D3=minor(C,3);
D3=D3,{w0};
ideal G3=std(D3);

print("DEGREE42_CONORMAL_REES");
print(onzero);
print(size(G5)==1 and G5[1]=={w0}^2);
print(size(G4)==1 and G4[1]=={w0});
print(size(G3)==1 and G3[1]==1);
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
    marker = compact.index("DEGREE42_CONORMAL_REES")
    assert compact[marker + 1 : marker + 5] == ["1", "1", "1", "1"], (
        process.stdout + process.stderr
    )
    print(
        "PASS: the degree-42 residual conormal map has maximal-minor "
        "ideal (w0^2), has rank exactly three on w0=0, and gives "
        "all-order formal synchronization on D(w0)"
    )


if __name__ == "__main__":
    main()
