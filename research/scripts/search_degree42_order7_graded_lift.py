#!/usr/bin/env python3
"""Extract the graded correction in the degree-42 order-seven lift.

At the specialization (e1,e2,t)=(1,2,3), let c6 be the first support
torsion witness modulo (u,v)^6.  If a saturated order-seven class c7 maps
to c6, then

    c7 - c6 = q + h,  q in I,  h in (u,v)^6.

Modulo (u,v)^7 only the normal-constant coefficients of the seven standard
generators of (u,v)^6 matter.  This script extracts that pure degree-six
correction and verifies c6+h directly at order seven.
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
    parser.add_argument("--prime", type=int, default=101)
    parser.add_argument("--timeout", type=int, default=900)
    args = parser.parse_args()

    normals, bases, residuals, _defect = transformed_problem()
    variables = normals + bases
    e1, e2, translation, w0, w1, w2 = bases
    specialization = {e1: sp.Integer(1), e2: sp.Integer(2), translation: sp.Integer(3)}
    factor_a, factor_b = support_factors(e1, e2, translation)
    factor_a = factor_a.subs(specialization)
    factor_b = factor_b.subs(specialization)
    characteristic = args.prime if args.prime else 0

    program = f"""
LIB "elim.lib";
ring source={characteristic},({",".join(map(str, variables))}),(dp(5),dp(6));
ideal I={",".join(serialize(item) for item in residuals)};
poly p3=subst(I[5],{normals[2]},0);
ideal I3=subst(I,{normals[2]},p3);
poly p4=subst(I3[11],{normals[3]},0);
ideal I4=subst(I3,{normals[3]},p4);
poly p5=subst(I4[17],{normals[4]},0);
ideal I5=subst(I4,{normals[4]},p5);

ring q={characteristic},(u,v,{w0},{w1},{w2}),(dp(2),dp(3));
map phi=source,u,v,0,0,0,1,2,3,{w0},{w1},{w2};
ideal Core=phi(I5);
ideal M=u,v;
ideal J6=Core,M^6;
ideal G6=slimgb(J6);
ideal M6=M^6;
matrix T6;
ideal L6=liftstd(J6,T6);
poly A={serialize(factor_a)};
poly B={serialize(factor_b)};
ideal K={w0},{w1}*{w2},A*B*{w2};
list Sat6=sat(G6,K);
ideal S6=slimgb(Sat6[1]);

int i;
poly c6=0;
for (i=1;i<=size(S6);i++)
{{
  if (c6==0)
  {{
    c6=reduce(S6[i],G6);
  }}
}}

ideal J7=Core,M^7;
ideal G7=slimgb(J7);
list Sat7=sat(G7,K);
ideal S7=slimgb(Sat7[1]);
ideal LiftIdeal=slimgb(S7+G6);
matrix L=lift(S7+G6,ideal(c6));
poly c7=0;
for (i=1;i<=size(S7);i++)
{{
  c7=c7+S7[i]*L[i,1];
}}
c7=reduce(c7,G7);

matrix DL=lift(L6,ideal(c7-c6));
matrix D=T6*DL;
poly h=0;
int coreSize=size(Core);
int coreSlots=ncols(matrix(Core));
poly reconstruction=0;
for (i=1;i<=nrows(D);i++)
{{
  reconstruction=reconstruction+J6[i]*D[i,1];
}}
for (i=1;i<=7;i++)
{{
  poly coefficient=D[coreSlots+i,1];
  coefficient=subst(subst(coefficient,u,0),v,0);
  h=h+J6[coreSlots+i]*coefficient;
}}
int reconstructs=(reconstruction==c7-c6);
int rawCorrection=(reduce(c7-c6-h,G7)==0);
poly gradedCorrection=h;
h=reduce(h,G7);
poly candidate=reduce(c6+h,G7);

print("DEGREE42_ORDER7_GRADED_LIFT");
print(size(Core));
print(coreSlots);
print(size(J6));
print(size(L6));
print(nrows(T6));
print(ncols(T6));
print(nrows(DL));
print(ncols(DL));
print(nrows(D));
print(reconstructs);
print(rawCorrection);
print(reduce(c6,G6)!=0);
print(reduce(c7,G7)!=0);
print(reduce(c7-c6,G6)==0);
print(reduce(c7-candidate,G7)==0);
print(reduce({w0}*candidate,G7)==0);
print(reduce({w2}^2*candidate,G7)==0);
print(candidate);
print(gradedCorrection);
print(h);
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
    except subprocess.TimeoutExpired:
        os.killpg(process.pid, signal.SIGTERM)
        try:
            process.wait(timeout=5)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGKILL)
            process.wait()
        raise
    if process.returncode or "? error occurred" in stdout:
        raise RuntimeError(stdout + stderr)
    marker = "DEGREE42_ORDER7_GRADED_LIFT"
    assert marker in stdout, stdout + stderr
    print(stdout[stdout.index(marker) :].strip())
    if stderr.strip():
        print(stderr.strip())


if __name__ == "__main__":
    main()
