#!/usr/bin/env python3
"""Search the degree-42 residual for support torsion on finite normal jets.

This is a discovery driver, not an all-order certificate.  It performs the
three exact unit-pivot eliminations, adds the requested power of the
two-variable normal ideal, and computes

    (I + m^n) : k^infinity,

where ``k=(w0,w1*w2,A*B*w2)``.  If the saturation grows, the driver extracts
the first generator-by-generator normal form; ideal-level ``reduce`` is
deliberately avoided because Singular can discard this quotient witness in
bulk reductions.
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
    parser.add_argument("--normal-power", type=int, default=6)
    parser.add_argument("--timeout", type=int, default=600)
    parser.add_argument(
        "--base-values",
        help="optional comma-separated specialization of e1,e2,t",
    )
    parser.add_argument(
        "--transition",
        action="store_true",
        help="test whether the first class lifts to the next normal power",
    )
    parser.add_argument(
        "--print-lift",
        action="store_true",
        help="with --transition, print one explicit next-order lift",
    )
    parser.add_argument(
        "--analyze-witness",
        action="store_true",
        help="also compute the cyclic annihilator and its radical",
    )
    args = parser.parse_args()
    assert args.normal_power >= 2
    assert not args.print_lift or args.transition

    normals, bases, residuals, _defect = transformed_problem()
    variables = normals + bases
    e1, e2, translation, w0, w1, w2 = bases
    factor_a, factor_b = support_factors(e1, e2, translation)
    normal_images = "u,v,0,0,0"
    characteristic = args.prime if args.prime else 0
    if args.base_values:
        base_values = tuple(
            sp.Rational(value.strip())
            for value in args.base_values.split(",")
        )
        assert len(base_values) == 3
        base_substitution = dict(
            zip((e1, e2, translation), base_values, strict=True)
        )
        factor_a = factor_a.subs(base_substitution)
        factor_b = factor_b.subs(base_substitution)
        coefficient_images = ",".join(map(str, base_values))
        target_variables = f"u,v,{w0},{w1},{w2}"
        target_order = "(dp(2),dp(3))"
        map_images = (
            f"{normal_images},{coefficient_images},{w0},{w1},{w2}"
        )
    else:
        target_variables = f"u,v,{','.join(map(str, bases))}"
        target_order = "(dp(2),dp(6))"
        map_images = (
            f"{normal_images},{','.join(map(str, bases))}"
        )

    analysis = ""
    if args.analyze_witness:
        analysis = """
ideal Annih=slimgb(quotient(G,ideal(witness)));
ideal Rad=std(radical(Annih));
ideal BaseRad=std(eliminate(Rad,u*v));
print(size(Annih));
print(size(Rad));
print(BaseRad);
"""
    transition = ""
    if args.transition:
        lift_output = ""
        if args.print_lift:
            lift_output = """
matrix LiftRep=lift(Snext+G,ideal(witness));
poly liftedWitness=0;
for (i=1;i<=size(Snext);i++)
{
  liftedWitness=liftedWitness+Snext[i]*LiftRep[i,1];
}
liftedWitness=reduce(liftedWitness,Gnext);
print(liftedWitness);
"""
        transition = f"""
ideal Jnext=phi(I5),M^{args.normal_power + 1};
ideal Gnext=slimgb(Jnext);
list saturationNext=sat(Gnext,K);
ideal Snext=slimgb(saturationNext[1]);
ideal Lift=slimgb(Snext+G);
int lifts=(reduce(witness,Lift)==0);
print(size(Gnext));
print(size(Snext));
print(lifts);
{lift_output}
"""

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

ring q={characteristic},({target_variables}),{target_order};
map phi=source,{map_images};
ideal M=u,v;
ideal J=phi(I5),M^{args.normal_power};
ideal G=slimgb(J);
poly A={serialize(factor_a)};
poly B={serialize(factor_b)};
ideal K={w0},{w1}*{w2},A*B*{w2};
list saturation=sat(G,K);
ideal S=slimgb(saturation[1]);

int i;
int first=0;
poly witness=0;
poly remainder=0;
for (i=1;i<=size(S);i++)
{{
  remainder=reduce(S[i],G);
  if (remainder!=0 and first==0)
  {{
    first=i;
    witness=remainder;
  }}
}}

print("DEGREE42_SUPPORT_SATURATION_JET");
print(size(G));
print(size(S));
print(first);
print(witness);
{analysis}
{transition}
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
    if process.returncode:
        raise RuntimeError(stdout + stderr)
    if "? error occurred" in stdout:
        raise RuntimeError(stdout + stderr)
    marker = "DEGREE42_SUPPORT_SATURATION_JET"
    assert marker in stdout, stdout + stderr
    print(stdout[stdout.index(marker) :].strip())
    if stderr.strip():
        print(stderr.strip())


if __name__ == "__main__":
    main()
