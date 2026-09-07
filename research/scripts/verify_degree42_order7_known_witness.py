#!/usr/bin/env python3
"""Directly test the known exact degree-42 witness at normal order seven."""

from __future__ import annotations

import argparse
import os
import shutil
import signal
import subprocess

from verify_degree42_transported_27_normal_jets import (
    serialize,
    transformed_problem,
)


C6 = (
    "6344815050*v^4*sync42p_w2"
    "+285348825*u^3*sync42p_w2^2"
    "-1644952050*u^3*sync42p_w1"
    "-151067025*u^2*v*sync42p_w2^2"
    "+2960913690*u^2*v*sync42p_w1"
    "-4209975*u*v^2*sync42p_w1*sync42p_w2"
    "+50355675*u*v^2*sync42p_w2^2"
    "+89503785*u*v^2*sync42p_w0"
    "-7754773950*u*v^2*sync42p_w1"
    "+67195575*v^3*sync42p_w1*sync42p_w2"
    "-839261250*v^3*sync42p_w2^2"
    "-212157225*v^3*sync42p_w0"
    "+657980820*v^3*sync42p_w1"
    "+4019868*u^2*sync42p_w1^2"
    "-8544420*u^2*sync42p_w0*sync42p_w2"
    "-80569080*u^2*sync42p_w1*sync42p_w2"
    "+563983560*u^2*sync42p_w0"
    "+1321164*u*v*sync42p_w0*sync42p_w1"
    "-29042064*u*v*sync42p_w1^2"
    "+10954980*u*v*sync42p_w0*sync42p_w2"
    "+684837180*u*v*sync42p_w1*sync42p_w2"
    "-281991780*u*v*sync42p_w0"
    "-666198*v^2*sync42p_w0*sync42p_w1"
    "+5578308*v^2*sync42p_w1^2"
    "-18594360*v^2*sync42p_w0*sync42p_w2"
    "-60426810*v^2*sync42p_w1*sync42p_w2"
    "+3806889030*v^2*sync42p_w0"
    "-660572*u*sync42p_w0^2"
    "+8641728*u*sync42p_w0*sync42p_w1"
    "-241707240*u*sync42p_w0*sync42p_w2"
    "-72784*v*sync42p_w0^2"
    "+5359824*v*sync42p_w0*sync42p_w1"
    "-107425440*v*sync42p_w0*sync42p_w2"
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime", type=int, default=0)
    parser.add_argument("--timeout", type=int, default=1200)
    args = parser.parse_args()

    normals, bases, residuals, _defect = transformed_problem()
    variables = normals + bases
    _e1, _e2, _translation, w0, w1, w2 = bases
    characteristic = args.prime if args.prime else 0

    program = f"""
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
ideal J7=Core,M^7;
ideal G7=slimgb(J7);
ideal H7=std(G7);
poly c6={C6};

print("DEGREE42_ORDER7_KNOWN_WITNESS");
print(size(G7));
print(size(H7));
print(reduce(c6,G7)!=0);
print(reduce({w0}*c6,H7)==0);
print(reduce({w0}*c6,G7)==0);
print(reduce({w2}^2*c6,G7)==0);
print(reduce({w2}*c6,G7)==0);
print(reduce(c6,G7));
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
    marker = "DEGREE42_ORDER7_KNOWN_WITNESS"
    assert marker in stdout, stdout + stderr
    print(stdout[stdout.index(marker) :].strip())
    if stderr.strip():
        print(stderr.strip())


if __name__ == "__main__":
    main()
