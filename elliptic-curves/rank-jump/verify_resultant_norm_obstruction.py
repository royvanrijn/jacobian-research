#!/usr/bin/env python3
"""Integer Sylvester/Bareiss replay of the constant norm obstruction."""
import argparse
import sys
from fractions import Fraction as Q
from math import gcd,lcm,isqrt
from pathlib import Path
import retrospective as r
import second_generic_presentation as source
import resultant_norm_obstruction as resultants
import verify_second_generic_presentation as rational
import verify_j_preimage_sign_trace as sign

OUTPUT=r.OUT/'rank_jump_resultant_norm_obstruction_verification_v1.json'


def primitive(coeff):
    a=list(map(Q,coeff));d=lcm(*(x.denominator for x in a));ints=[int(d*x) for x in a];c=gcd(*ints)
    return [x//c for x in ints],Q(c,d)


def determinant(matrix):
    A=[row[:] for row in matrix];n=len(A);den=1;sgn=1
    for k in range(n-1):
        p=next((j for j in range(k,n) if A[j][k]),None)
        if p is None:return 0
        if p!=k:A[k],A[p]=A[p],A[k];sgn=-sgn
        pivot=A[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                num=A[i][j]*pivot-A[i][k]*A[k][j];v,rem=divmod(num,den);assert rem==0;A[i][j]=v
        for i in range(k+1,n):A[i][k]=0
        den=pivot
    return sgn*A[-1][-1]


def sylvester(a,b):
    a,ca=primitive(a);b,cb=primitive(b);m,n=len(a)-1,len(b)-1;a=a[::-1];b=b[::-1];rows=[]
    for i in range(n):rows.append([0]*i+a+[0]*(n-1-i))
    for i in range(m):rows.append([0]*i+b+[0]*(m-1-i))
    return ca**n*cb**m*determinant(rows)


def compute():
    # Historical coordinate scalings have certified integers beyond Python's
    # default 4300-digit conversion limit; this changes no arithmetic bound.
    sys.set_int_max_str_digits(100000)
    inputs=r.read(source.INPUT);data=r.read(resultants.OUTPUT);old=r.read(source.OUTPUT);v=r.read(sign.OUTPUT)
    for obj in (inputs,data,old,v,r.read(rational.OUTPUT)):
        for p,sha in obj['bindings'].items():assert r.digest((r.ROOT/p).read_bytes())==sha,p
    families={x['family']:x for x in inputs['families']};rows=[];byfamily={}
    for x in data['families']:
        f=families[x['family']];R=sylvester(f['A'],f['B']);assert str(R)==x['resultant'] and R>0
        a,b=int(x['numerator_floor_sqrt']),int(x['denominator_floor_sqrt'])
        assert a*a<=R.numerator<(a+1)**2 and b*b<=R.denominator<(b+1)**2
        assert a*a!=R.numerator or b*b!=R.denominator;assert x['rational_square'] is False
        rows.append({'family':x['family'],'sylvester_dimension':20,'resultant_nonsquare':True});byfamily[x['family']]=R
    # Norm squareclasses are invariant under the independently proved base changes.
    duplicate=r.read(rational.OUTPUT)['global_base_equivalences']
    for d in duplicate:assert rational.square_root(byfamily[d['left_family']]/byfamily[d['right_family']]) is not None
    cells=[];sign_rows={(x['token'],x['family']):x for x in v['rows']}
    for row in old['rows']:
        assert row['finite_degree']==24 and row['infinity_multiplicity']==0
        L=Q(row['j_equation_ascending'][-1]);R=byfamily[row['family']]
        predicted_norm=L**4/R
        for pre in row['preimages']:
            assert pre['rationally_isomorphic'] and pre['multiplicity']==1
            u2=Q(pre['scaling_square']);assert rational.square_root(u2) is not None
            predicted_norm/=u2
        # Check the consequent nonsquare also as an exact rational integer-square test.
        assert rational.square_root(predicted_norm) is None
        s=sign_rows[row['token'],row['family']];assert s['rational_irreducibility']=='PASS'
        assert s['preimage_degree']==24-len(row['preimages'])
        cells.append({'token':row['token'],'family':row['family'],'irreducible_residual_degree':s['preimage_degree'],
            'residual_scaling_norm_nonsquare':True,'local_sign_witness_agrees':True})
    return {'schema':'rank-jump.resultant-norm-obstruction-verification.v1','status':'PASS','families':rows,'cells':cells,
        'bindings':source.bindings([Path(__file__),resultants.OUTPUT,source.INPUT,source.OUTPUT,rational.OUTPUT,sign.OUTPUT]),
        'boundary':'Exact integer Sylvester determinants and square-root bounds independently verify all seven nonsquare resultants; the proved norm identity reproduces all112 sign obstructions. The identity is an algebraic theorem proved in the note, not a separate numerical norm calculation in112 number fields.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);args=p.parse_args();result=compute()
    if args.mode=='build':r.write_new(OUTPUT,result)
    else:assert result==r.read(OUTPUT)
    print('PASS seven integer Sylvester determinants; 112 norm obstructions')
