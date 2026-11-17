#!/usr/bin/env python3
"""Finite exact controls for QSO3A; universal coverage uses the written proof."""
from __future__ import annotations

import argparse
from math import comb, factorial
from pathlib import Path
import hashlib
import json
import tempfile

import sympy as sp

HERE = Path(__file__).resolve().parent
PDF = HERE/'sources/beals-sattinger-solv-int-9801010v1.pdf'
PDF_SHA = '05e7df65dac43b343791d5a6ee7a02f0cd5ed662738a2c1530a105f2a21ac738'
x = sp.symbols('x')
u, w = sp.Function('u')(x), sp.Function('w')(x)


def clean(op):
    return {i: sp.expand(c) for i,c in op.items() if sp.expand(c) != 0}


def mul(left, right, minimum=0):
    out = {}
    for i, a in left.items():
        for j, b in right.items():
            for k in range(max(-1, i+j-minimum)+1):
                if i >= 0 and k > i:
                    break
                degree = i+j-k
                out[degree] = out.get(degree, 0)+sp.binomial(i,k)*a*sp.diff(b,x,k)
    return clean(out)


def subtract(a,b):
    return clean({i:a.get(i,0)-b.get(i,0) for i in set(a)|set(b)})


def comm(a,b):
    return subtract(mul(a,b),mul(b,a))


def power(op, n, minimum=0):
    out = {0:sp.Integer(1)}
    for j in range(1,n+1):
        out = mul(out,op,minimum-(n-j))
    return out


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    assert sha(PDF)==PDF_SHA

    P={3:sp.Integer(1),1:u,0:w}
    root={1:sp.Integer(1)}
    for k in range(1,5):
        trial=mul(mul(root,root,1-k),root,2-k)
        target={1:u,0:w}.get(2-k,0)
        root[-k]=sp.expand((target-trial.get(2-k,0))/3)
    assert root[-1]==u/3
    assert sp.expand(root[-2]-(w-sp.diff(u,x))/3)==0
    assert {i:c for i,c in power(root,3).items() if i>=0}==P
    flow_controls=[]
    for j in [1,2,4,5]:
        Bj=power(root,j)
        actual=comm(P,Bj)
        assert set(actual)<={0,1},(j,actual)
        flow_controls.append({'flow':j,'commutator_orders':sorted(actual)})

    # Universal positive coefficient formula for an arbitrary negative part.
    r1,r2,r3=[sp.Function('r'+str(i))(x) for i in [1,2,3]]
    negative={-1:r1,-2:r2,-3:r3}
    actual=subtract(mul(negative,P),mul(P,negative))
    expected={1:-3*sp.diff(r1,x),0:-3*(sp.diff(r1,x,2)+sp.diff(r2,x))}
    assert subtract(actual,expected)=={}
    f=sp.Function('f')(x)
    for j in range(9):
        assert sp.expand(comm(P,{j:f})[j+2]-3*sp.diff(f,x))==0
    corrupted=power(root,5)
    corrupted[3]+=f
    assert not set(comm(P,corrupted))<={0,1}

    # Positive Airy case and removal of a commuting cubic-order polynomial.
    Airy={3:sp.Integer(1),0:x}
    mate={1:-sp.Integer(1)}
    assert comm(Airy,mate)=={0:sp.Integer(1)}
    mixed=clean({i:mate.get(i,0)+power(Airy,2).get(i,0)
                 for i in set(mate)|set(power(Airy,2))})
    assert max(mixed)==6 and comm(Airy,mixed)=={0:sp.Integer(1)}
    assert subtract(mixed,power(Airy,2))==mate
    assert -3*(x/3)/(-1)==x

    # Any finite scalar Q gives a spectral matrix polynomial in lambda.
    lam=sp.symbols('lambda')
    companion=sp.Matrix([[0,1,0],[0,0,1],[lam-w,-u,0]])
    rows=[sp.Matrix([[1,0,0]])]
    for _ in range(9):
        rows.append(rows[-1].diff(x)+rows[-1]*companion)
    qs=[sp.Function('q'+str(j))(x) for j in range(8)]
    row=sum((qs[j]*rows[j] for j in range(8)),sp.zeros(1,3))
    second=row.diff(x)+row*companion
    third=second.diff(x)+second*companion
    spectral=sp.Matrix.vstack(row,second,third)
    spectral_degrees=[sp.Poly(sp.expand(c),lam).degree() for c in spectral]
    assert max(spectral_degrees)==3
    assert all(d>=0 for d in spectral_degrees)

    # Correct left-ratio bounds force a unipotent lower triangular first row.
    z=sp.symbols('z',nonzero=True)
    a,b,c,d=sp.symbols('a b c d')
    G=sp.Matrix([[1,0,0],[a,1,0],[b*z+c,d,1]])
    diag=sp.diag(1,z,z**2)
    scaled=diag.inv()*G*diag-sp.eye(3)
    assert all(sp.limit(t,z,sp.oo)==0 for t in scaled)
    assert list(G[0,:])==[1,0,0]
    wrong=sp.eye(3);wrong[0,1]=1
    assert sp.limit((diag.inv()*wrong*diag)[0,1],z,sp.oo)==sp.oo

    # Direct residue expansion against the terminating hypergeometric/Jacobi form.
    U,W,Z=sp.symbols('U W Z')
    residue_controls=[]
    for m in [1,2,4,5,7,8,10,11,13,14,16,17,19,20,25,26,31,32]:
        N=m+4;s=sp.Rational(m+3,3)
        H=sp.Integer(0)
        for aa in range(N//2+1):
            for bb in range(N//3+1):
                if 2*aa+3*bb==N:
                    H+=sp.binomial(s,aa+bb)*comb(aa+bb,aa)*U**aa*W**bb
        choices=[(aa,bb,(N-2*aa-3*bb)//6) for aa in [0,1] for bb in [0,1]
                 if N>=2*aa+3*bb and (N-2*aa-3*bb)%6==0]
        assert len(choices)==1
        aa,bb,l=choices[0]
        nu=sp.Rational(2+2*aa,3);eps=sp.Rational(1-2*bb,2)
        C0=sp.binomial(s,aa+bb+2*l)*comb(aa+bb+2*l,aa)
        poly=sum((sp.rf(-l,k)*sp.rf(-l+eps,k)/(sp.rf(nu,k)*factorial(k))*Z**k
                  for k in range(l+1)),sp.Integer(0))
        candidate=C0*U**aa*W**(bb+2*l)*poly.subs(Z,-4*U**3/(27*W**2))
        assert sp.cancel(H-candidate)==0,m
        jacobi=((1-Z)**l*factorial(l)/sp.rf(nu,l)
                *sp.jacobi(l,nu-1,-eps,1-2*Z/(Z-1)))
        assert sp.cancel(poly-jacobi)==0,m
        assert nu-1>-1 and -eps>-1
        if m<=14:
            result=sp.Poly(sp.resultant(sp.diff(H,U),sp.diff(H,W),U),W)
            assert len(result.terms())==1 and result.as_expr()!=0,m
        residue_controls.append({'m':m,'a':aa,'b':bb,'jacobi_degree':l,
                                 'alpha':str(nu-1),'beta':str(-eps),
                                 'residue_matches':True,'Pfaff_matches':True})

    names=['ALL_CUBIC_ALGEBRAIC_STRING_OBSTRUCTION.md',
           'check_all_cubic_algebraic_string.py','all_cubic_source_audit.json']
    receipt={
        'claim':'QSO3A','status':'PASS','flow_controls':flow_controls,
        'generic_negative_part_time_formula':True,
        'constant_leading_coefficient_orders_checked':list(range(9)),
        'corrupted_flow_rejected':True,
        'positive_Airy_and_order6_polynomial_subtraction':True,
        'spectral_matrix_order7_polynomial_degrees':spectral_degrees,
        'left_ratio_triangular_bound_control':True,
        'terminating_residue_Jacobi_controls':residue_controls,
        'source_pdf_sha256':sha(PDF),
        'hashes':{name:sha(HERE/name) for name in names},
        'coverage_boundary':(
            'Finite exact fractional-power flows1,2,4,5, universal negative-part '
            'positive coefficients, a positive Airy pair, order7 polynomial '
            'spectral reduction, left-ratio bounds, terminating residue/Pfaff '
            'identities in declared m, and selected monomial gradient resultants. '
            'All-order Stokes injectivity application, algebraic continuation, '
            'Jacobi root simplicity and infinity exclusion use the written '
            'proof and primary-source audit, not this finite set.'),
        'independent_replay':False,'formal_verification':False,'external_review':False,
    }
    output=args.output or Path(tempfile.mkdtemp(prefix='qso3a-check-'))/'all_cubic_controls.json'
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'claim':'QSO3A','status':'PASS','receipt':str(output)},indent=2))


if __name__=='__main__':
    main()
