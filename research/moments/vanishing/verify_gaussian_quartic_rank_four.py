"""Exact rank-four source, kernel and principal-part controls.

The companion proof supplies independence and all-order completeness.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import platform
import sympy as sp
from gaussian_quartic_rank_four import (
    Z,W,Y,U,T,CONNECTION,row,laurent_period,reduce_marker,reduce_time_amplitude,
    connection_derivative,is_in_full_image,truncated_series,original_projection,
    original_projection_allows,
)


def zero(expr):
    assert sp.cancel(sp.together(expr))==0,expr


def direct_coefficient(a,b,j,m):
    if (m+j)%2 or not 0<=m+b-a<=2*m:
        return sp.Integer(0)
    return (sp.factorial(m+b)/sp.factorial(m)*sp.binomial(2*m,m+b-a)
            *sp.factorial2(m+j-1))


def checks():
    P=W*(1+Z)**2*Y
    assert reduce_marker(1)==row([0,1,0,0])
    assert reduce_marker(P)==row([0,0,1,0])
    assert reduce_marker(1-Z)==row([0,0,0,1])
    assert reduce_time_amplitude(1-2*T*(1+Z)*Y)==row([1,0,0,0])
    assert reduce_marker(2*(1+Z)*Y)==row([-1/T,1/T,0,0])

    markers=[(0,0,0),(1,0,0),(2,0,0),(3,0,0),(0,1,0),
             (0,0,1),(0,0,2),(0,0,3),(1,1,0),(1,0,2),
             (1,1,1),(0,2,0)]
    for a,b,j in markers:
        coordinates=reduce_marker(Z**a*W**b*Y**j)
        actual=truncated_series(coordinates,8)
        for m in range(9):
            zero(actual.coeff(T,m)-direct_coefficient(a,b,j,m))
        assert is_in_full_image(coordinates)
        assert original_projection_allows(coordinates)

    kernel=(1+Z)**2*(Z*W-Y*Y)
    assert reduce_marker(kernel)==row([0,0,0,0])
    assert reduce_marker(kernel+1)==row([0,1,0,0])
    for m in range(9):
        actual=0
        for (a,b,j),coefficient in sp.Poly(sp.expand(kernel),Z,W,Y).terms():
            actual+=coefficient*direct_coefficient(a,b,j,m)
        zero(actual)

    # Exact recurrence and frame determinant, not a rank inference from jets.
    for n in range(-5,9):
        assert (2*laurent_period(n+3)-3*laurent_period(n+2)+laurent_period(n+1)
                -n*T*T*laurent_period(n-1)).applyfunc(sp.cancel)==row([0,0,0,0])
    frame=sp.Matrix.vstack(*(laurent_period(n) for n in [-1,0,1,2]))
    zero(frame.det()-2*T**3)
    h2=reduce_marker(Z*Z)
    assert h2==row([1/(2*T*T),1-1/(2*T*T),8*T,2])
    rows=[h2]
    for _ in range(3):
        rows.append(connection_derivative(rows[-1]))
    cyclic=sp.factor(sp.Matrix.vstack(*rows).det())
    zero(cyclic-3*(60*T*T+1)/(64*T**9))

    # The two decisive rational-ODE coefficient polynomials in the written proof.
    N=sp.symbols("N")
    zero(N*(N-1)-3*N+sp.Rational(15,4)-(N-sp.Rational(3,2))*(N-sp.Rational(5,2)))
    zero(N/sp.Integer(16)-sp.Rational(3,16)-(N-3)/16)

    # M membership uses only the certified finite principal part of the input.
    assert is_in_full_image([1,0,0,0])
    assert not original_projection_allows([1,0,0,0])
    assert not is_in_full_image([1/T,0,0,0])
    assert is_in_full_image([-1/T,1/T,0,0])
    assert is_in_full_image([-1/T,0,0,1/T])
    assert not is_in_full_image([1/(1-T),0,0,0])
    for n in range(1,13):
        projection=original_projection(laurent_period(n))
        assert sp.degree(projection,T)<=((n-1)//2)
    for n in range(-8,0):
        projection=original_projection(laurent_period(n))
        if projection:
            numerator,denominator=sp.fraction(projection)
            assert sp.degree(numerator,T)-sp.degree(denominator,T)<0
    return {
        "status":"passed",
        "actual_marker_controls":{"monomials":len(markers),"spatial_degree_max":3,"moment_orders":[0,8]},
        "exact_kernel":{"zero":"(1+Z)^2(ZW-Y^2)","rejected_corruption":"kernel+1"},
        "recurrence_indices":[-5,8],"frame_determinant":"2t^3",
        "cyclic_Z_squared_determinant":str(cyclic),
        "full_image_membership":"Laurent coordinates plus finite certified principal part",
        "literal_V_boundary":"1 is in M but its constant-line projection excludes it from V",
        "scope":"Exact controls; the written rational-ODE and source arguments prove rank four, global completeness and the original projection image.",
    }


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path)
    args=parser.parse_args()
    result=checks()
    here=Path(__file__).resolve()
    utility=here.with_name("gaussian_quartic_rank_four.py")
    proof=here.with_name("GAUSSIAN_QUARTIC_BASELINE_RANK_FOUR.md")
    receipt={"status":"passed","checked_on":"2026-11-17",
             "environment":{"python":platform.python_version(),"sympy":sp.__version__},
             "checker_sha256":hashlib.sha256(here.read_bytes()).hexdigest(),
             "utility_sha256":hashlib.sha256(utility.read_bytes()).hexdigest(),
             "proof_sha256":hashlib.sha256(proof.read_bytes()).hexdigest(),
             "checks":result,"independent_replay":False,"formal_verification":False,"external_review":False}
    if args.output:
        args.output.write_text(json.dumps(receipt,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(receipt,indent=2))

if __name__=="__main__":
    main()