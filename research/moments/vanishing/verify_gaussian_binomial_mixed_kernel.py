"""Bounded exact controls for the all-m original Stein-kernel proof."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import platform
import sympy as sp
from gaussian_binomial_mixed_kernel import (
    Z,W,Y,Q,delta_one,delta_two,normal_form,coefficient_classes,
    canonical_monomials,is_in_kernel,
)


def zero(value):
    assert sp.cancel(value)==0,value


def gaussian_moment(n):
    return sp.Integer(0) if n%2 else sp.factorial2(n-1)


def actual_coefficient(amplitude,m,n):
    answer=0
    for (a,b,j),coefficient in sp.Poly(sp.expand(amplitude),Z,W,Y).terms():
        if 0<=n+b-a<=m*n:
            answer+=(coefficient*sp.factorial(n+b)/sp.factorial(n)
                     *sp.binomial(m*n,n+b-a)*gaussian_moment(n+j))
    return answer


def controls():
    count=0
    for m in range(1,6):
        for a in range(4):
            for b in range(4-a):
                for j in range(4-a-b):
                    amplitude=Z**a*W**b*Y**j
                    normal,first,second=normal_form(amplitude,m)
                    zero(amplitude-normal-delta_one(first)-delta_two(second,m))
                    assert sp.Poly(normal,Z,W,Y).total_degree()<=a+b+j
                    ratios=coefficient_classes(amplitude,m)
                    for n in range(9):
                        q=n//2
                        base=(sp.factorial2(2*q-1)*sp.binomial(2*m*q,2*q) if n%2==0
                              else sp.factorial2(2*q+1)*sp.binomial(m*(2*q+1),2*q+1))
                        zero(actual_coefficient(amplitude,m,n)-base*ratios[n%2].subs(Q,q))
                    count+=1
        for degree in range(5):
            canonical=canonical_monomials(degree,m)
            expected=2*degree+1 if m==1 else (1 if degree==0 else 4*degree)
            assert len(canonical)==expected
            for parity in range(2):
                ratios=[coefficient_classes(value,m)[parity] for value in canonical]
                ratios=[value for value in ratios if value]
                if ratios:
                    matrix=sp.Matrix([[value.subs(Q,q) for value in ratios]
                                      for q in range(len(ratios))])
                    assert matrix.det()!=0
        for marker in [1,Z,W,Y,Z*W,Y**2,Z+W+Y]:
            assert is_in_kernel(delta_one(marker),m)
            assert is_in_kernel(delta_two(marker,m),m)
        assert not is_in_kernel(1,m)
        for n in range(9):
            # Polynomial-time Jacobian identity Phi[1-t m(1+Z)^(m-1)Y]=1.
            coefficient=actual_coefficient(1,m,n)
            if n:
                coefficient-=actual_coefficient(m*(1+Z)**(m-1)*Y,m,n-1)
            zero(coefficient-(1 if n==0 else 0))
    # Nonzero literal-kernel example shared with the faithful rank-four model.
    kernel=(1+Z)**2*(Z*W-Y**2)
    assert is_in_kernel(kernel,2)
    assert not is_in_kernel(kernel+1,2)
    from gaussian_quartic_rank_four import reduce_marker
    for marker in [1,Z,W,Y,Z*W,Y**2,Z**3,W**2*Y,kernel]:
        normal,_,_=normal_form(marker,2)
        assert reduce_marker(marker)==reduce_marker(normal)
    return {
        "status":"passed","m_values":[1,2,3,4,5],
        "original_monomials":count,"spatial_degree_max":3,"moment_orders":[0,8],
        "filtered_dimension_controls":{"degree_range":[0,4],"m_ge_2":"4d for d>=1; 1 for d=0","m_1":"2d+1"},
        "checks":["exact operator reconstruction","direct Wick versus rational coefficient classes",
                  "bounded canonical-basis independence","two operator images in exact kernel",
                  "time-dependent Jacobian unit identity","faithful quartic specialization"],
        "scope":"Bounded exact controls only. The written degree descent and rational-pole proof establish all m, all markers and the complete literal kernel/image.",
    }


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path)
    args=parser.parse_args()
    result=controls()
    here=Path(__file__).resolve()
    files=[here,here.with_name("gaussian_binomial_mixed_kernel.py"),
           here.with_name("GAUSSIAN_BINOMIAL_MIXED_LITERAL_KERNEL.md"),
           here.with_name("gaussian_quartic_rank_four.py")]
    receipt={"status":"passed","checked_on":"2026-11-17",
             "environment":{"python":platform.python_version(),"sympy":sp.__version__},
             "sha256":{path.name:hashlib.sha256(path.read_bytes()).hexdigest() for path in files},
             "controls":result,"independent_replay":False,"formal_verification":False,"external_review":False}
    if args.output:
        args.output.write_text(json.dumps(receipt,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(receipt,indent=2))


if __name__=="__main__":
    main()
