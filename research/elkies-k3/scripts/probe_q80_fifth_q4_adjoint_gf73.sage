#!/usr/bin/env sage
"""Bounded adjoint-ideal probe for the fifth q=4 compact curve.

The raw line through the selected section leaves a degree-13 residual because
the stored bidegree-(14,3) model is highly singular.  This probe asks
Singular only for the conductor/adjoint ideal of its generic projective plane
closure over GF(73)(t).  It has a strict timeout and does not normalize the
curve or construct a Jacobian.
"""

import argparse
import hashlib
import json
import os
import re
import subprocess
from pathlib import Path

from sage.all import FunctionField, GF, PolynomialRing


ROOT = Path(__file__).resolve().parents[2]
ARTIFACT = (
    ROOT / "artifacts/generated-results/"
    "q80-fourth-q12-cm24-moving-cubic-gf73.json"
)
artifact_bytes = ARTIFACT.read_bytes()
assert hashlib.sha256(artifact_bytes).hexdigest() == (
    "c6560b3db2d1232866e9996fc727924090aa46293c2482885cf9f9dbf4c21c89"
)
data = json.loads(artifact_bytes)

parser = argparse.ArgumentParser()
parser.add_argument("--timeout", type=int, default=30)
parser.add_argument("--specialize", type=int)
arguments = parser.parse_args()

finite = GF(73, impl="modn")
if arguments.specialize is None:
    parameter_field = FunctionField(finite, "t")
    t = parameter_field.gen()
    coefficient_ring = parameter_field
    singular_ring = "ring R=(73,t),(v,x,z),dp;"
else:
    t = finite(arguments.specialize)
    coefficient_ring = finite
    singular_ring = "ring R=73,(v,x,z),dp;"
projective_ring = PolynomialRing(coefficient_ring, names=("v", "x", "z"))
v, x, z = projective_ring.gens()
degree = 14
projective = projective_ring(sum(
    coefficient_ring(coefficient)*t**t_degree
    *v**v_degree*x**x_degree*z**(degree-v_degree-x_degree)
    for t_degree, v_degree, x_degree, coefficient
    in data["moving_terms_T_v_x_coefficient"]
))
assert projective.is_homogeneous()
assert projective.total_degree() == degree

sage_local = Path(os.environ.get("SAGE_LOCAL", ""))
singular_path = sage_local/"bin/Singular"
if not singular_path.is_file():
    singular_path = Path("/opt/homebrew/bin/Singular")

program = f'''
LIB "paraplanecurves.lib";
{singular_ring}
poly f={projective};
ideal AI=adjointIdeal(f,4);
"Q80FIFTHADJOINT|generators="+string(size(AI));
int i;
for(i=1;i<=size(AI);i++)
{{
  "Q80FIFTHADJOINTGEN|"+string(i)+"|"+string(deg(AI[i]))+"|"+string(AI[i]);
}}
'''
try:
    completed = subprocess.run(
        [str(singular_path), "-q"],
        input=program,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=arguments.timeout,
        check=False,
    )
except subprocess.TimeoutExpired as error:
    partial = error.stdout or ""
    if isinstance(partial, bytes):
        partial = partial.decode(errors="replace")
    print(partial, end="", flush=True)
    print(
        "Q80FIFTHADJOINT|"
        f"timeout_seconds={arguments.timeout}|status=BOUNDED_TIMEOUT",
        flush=True,
    )
    raise SystemExit(2)

print(completed.stdout, end="", flush=True)
if completed.returncode:
    raise SystemExit(completed.returncode)
lines = tuple(
    line for line in completed.stdout.splitlines()
    if line.startswith("Q80FIFTHADJOINTGEN|")
)
if not lines:
    raise RuntimeError("Singular returned no adjoint generators")
degrees = tuple(int(line.split("|", 3)[2]) for line in lines)
print(
    "Q80FIFTHADJOINT|"
    f"specialize={arguments.specialize}|projective_degree={degree}|"
    f"generator_degrees={degrees}|"
    "normalization_used=0|jacobian_used=0|status=PASS_BOUNDED_ADJOINT",
    flush=True,
)
