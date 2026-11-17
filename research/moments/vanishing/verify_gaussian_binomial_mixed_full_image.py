"""Bounded exact controls for all-m faithful rank and global-image proof."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import platform
import sympy as sp
from gaussian_binomial_mixed_full_image import (
    Z,W,Y,T,Q,row,source_period,connection_derivative,reduce_marker,
    primitive_marker,unit_operator,unit_primitive_frame,truncated_series,is_in_full_image,
)


def zero(value):
    assert sp.cancel(value)==0,value


def zero_row(value,m):
    assert value.applyfunc(sp.cancel)==row([0]*(2*m),m)


def direct_coefficient(a,b,j,m,n):
    if (n+j)%2 or not 0<=n+b-a<=m*n:
        return sp.Integer(0)
    return (sp.factorial(n+b)/sp.factorial(n)*sp.binomial(m*n,n+b-a)
            *sp.factorial2(n+j-1))


def controls():
    markers=[(0,0,0),(1,0,0),(2,0,0),(0,1,0),(0,0,1),
             (0,0,2),(1,1,0),(1,0,2),(1,1,1),(0,2,0)]
    for m in range(1,5):
        for a,b,j in markers:
            coordinates=reduce_marker(Z**a*W**b*Y**j,m)
            actual=truncated_series(coordinates,m,6)
            for n in range(7):
                zero(actual.coeff(T,n)-direct_coefficient(a,b,j,m,n))
            assert is_in_full_image(coordinates,m)
        for s in range(-m-2,m+2):
            zero_row(m*source_period(s+2*m,m)-(2*m-1)*source_period(s+2*m-1,m)
                     +(m-1)*source_period(s+2*m-2,m)
                     -(s+m-1)*T*T*source_period(s,m),m)
        one=m*source_period(1,m)-(m-1)*source_period(0,m)
        zero_row(connection_derivative(one,m),m)
        zero_row(source_period(0,m)-T*reduce_marker(m*(1+Z)**(m-1)*Y,m)-one,m)
        assert is_in_full_image(one,m)
        assert not is_in_full_image(one/T,m)
        assert is_in_full_image((source_period(0,m)-one)/T,m)
        assert not is_in_full_image(one/(1-T),m)

    determinants={}
    for m in range(1,4):
        determinant=sp.factor(unit_primitive_frame(m).det())
        assert determinant!=0
        numerator,denominator=sp.fraction(determinant)
        assert len(sp.Poly(numerator,T).terms())==len(sp.Poly(denominator,T).terms())==1
        determinants[str(m)]=str(determinant)
    for m in range(2,6):
        c,n,d=unit_operator(m)
        for q in range(1,6):
            a=sp.factorial(2*m*q)/(2**q*sp.factorial(q)*sp.factorial(2*(m-1)*q))
            previous=sp.factorial(2*m*(q-1))/(2**(q-1)*sp.factorial(q-1)*sp.factorial(2*(m-1)*(q-1)))
            zero(a/previous-c*n.subs(Q,q)/d.subs(Q,q))
        assert sp.gcd(n,d)==1
        assert n.subs(Q,sp.Rational(1,2))!=0 and d.subs(Q,sp.Rational(1,2))!=0
        marker=primitive_marker(m)
        for order in range(9):
            coefficient=sp.expand(marker*(1+Z)**(m*order)).coeff(Z,order)
            zero(coefficient-sp.binomial(m*order,order)/sp.Integer(order+1))
        k=reduce_marker(marker,m)
        zero_row(T*connection_derivative(k,m)+k-source_period(0,m),m)
    # Exceptional m=1 inhomogeneous and homogeneous equations, exactly in coordinates.
    m=1
    f=source_period(0,m)
    one=source_period(1,m)
    zero_row(T**3*connection_derivative(f,m)+(T*T-1)*f+one,m)
    fp=connection_derivative(f,m)
    zero_row(T**3*connection_derivative(fp,m)+(4*T*T-1)*fp+2*T*f,m)
    return {
        "status":"passed","direct_original_controls":{"m_values":[1,2,3,4],"monomials":len(markers),"moment_orders":[0,6]},
        "frame_determinants":determinants,"unit_ratio_and_primitive_m_values":[2,3,4,5],
        "checks":["all-integer recurrence instances","exact connection and Jacobian unit","finite principal-part membership",
                  "unit factorial quotient and disjoint half-integer orbit","closed polynomial primitive identity",
                  "exceptional m1 inhomogeneous and homogeneous equations"],
        "scope":"Bounded exact controls. The written extreme-pole and differential-lattice proofs establish all-m actual ranks, global lattices and analytic boundaries.",
    }


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path)
    args=parser.parse_args()
    result=controls()
    here=Path(__file__).resolve()
    files=[here,here.with_name("gaussian_binomial_mixed_full_image.py"),
           here.with_name("GAUSSIAN_BINOMIAL_MIXED_FULL_IMAGE.md")]
    receipt={"status":"passed","checked_on":"2026-11-17",
             "environment":{"python":platform.python_version(),"sympy":sp.__version__},
             "sha256":{path.name:hashlib.sha256(path.read_bytes()).hexdigest() for path in files},
             "controls":result,"independent_replay":False,"formal_verification":False,"external_review":False}
    if args.output:
        args.output.write_text(json.dumps(receipt,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(receipt,indent=2))


if __name__=="__main__":
    main()
