#!/usr/bin/env python3
"""Verify an order-seven degree-42 lift using the single w0-colon.

The modular order-seven support-torsion class is already annihilated by w0.
Thus a lift of the exact order-six class can be found without computing the
full support saturation: if Jn=I+(u,v)^n, it is enough to test

    c6 in (J7 : w0) + J6.

When this membership holds, the component in J7:w0 gives an explicit
order-seven representative congruent to c6 modulo J6.
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
from verify_degree42_order7_known_witness import C6


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime", type=int, default=0)
    parser.add_argument("--timeout", type=int, default=1200)
    args = parser.parse_args()

    normals, bases, residuals, _defect = transformed_problem()
    variables = normals + bases
    e1, e2, translation, w0, w1, w2 = bases
    characteristic = args.prime if args.prime else 0
    specialization = (sp.Integer(1), sp.Integer(2), sp.Integer(3))

    program = f"""
LIB "elim.lib";
LIB "primdec.lib";
ring source={characteristic},({",".join(map(str, variables))}),(dp(5),dp(6));
ideal I={",".join(serialize(item) for item in residuals)};
poly p3=subst(I[5],{normals[2]},0);
ideal I3=subst(I,{normals[2]},p3);
poly p4=subst(I3[11],{normals[3]},0);
ideal I4=subst(I3,{normals[3]},p4);
poly p5=subst(I4[17],{normals[4]},0);
ideal I5=subst(I4,{normals[4]},p5);

ring q={characteristic},(u,v,{w0},{w1},{w2}),(dp(2),dp(3));
map phi=source,u,v,0,0,0,{specialization[0]},{specialization[1]},
  {specialization[2]},{w0},{w1},{w2};
ideal Core=phi(I5);
ideal M=u,v;
ideal J6=Core,M^6;
ideal G6=slimgb(J6);
poly c6={C6};
int i;

ideal J7=Core,M^7;
ideal G7=slimgb(J7);
ideal H7=std(G7);
int directBefore=(reduce({w0}*c6,H7)==0);
ideal Q7=slimgb(quotient(G7,ideal({w0})));
ideal Transition=slimgb(Q7+G6);
int lifts=(reduce(c6,Transition)==0);

poly c7=0;
if (lifts)
{{
  matrix L=lift(Q7+G6,ideal(c6));
  for (i=1;i<=size(Q7);i++)
  {{
    c7=c7+Q7[i]*L[i,1];
  }}
  c7=reduce(c7,G7);
}}

print("DEGREE42_ORDER7_COLON_LIFT");
print(size(G6));
print(size(G7));
print(size(Q7));
print(directBefore);
print(lifts);
print(reduce(c6,G6)!=0);
print(reduce(c7,G7)!=0);
print(reduce(c7-c6,G6)==0);
print(reduce({w0}*c7,G7)==0);
print(c6);
print(c7);
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
    marker = "DEGREE42_ORDER7_COLON_LIFT"
    assert marker in stdout, stdout + stderr
    print(stdout[stdout.index(marker) :].strip())
    if stderr.strip():
        print(stderr.strip())


if __name__ == "__main__":
    main()
