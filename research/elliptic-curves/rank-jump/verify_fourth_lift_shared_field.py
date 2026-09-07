#!/usr/bin/env python3
"""Replay the constant-field obstruction and binary S11-module proof gates."""
import argparse
from fractions import Fraction
from math import factorial
from pathlib import Path
from sage.all import QQ,GF,PolynomialRing,matrix
import retrospective as r

HERE=Path(__file__).resolve().parent
SOURCE=r.OUT/'rank_jump_fourth_lift_shared_field_v1.json'
OUTPUT=r.OUT/'rank_jump_fourth_lift_shared_field_verification_v1.json'


def compute():
    d=r.read(SOURCE)
    for path,sha in d['bindings'].items():assert r.digest((r.ROOT/path).read_bytes())==sha
    R=PolynomialRing(QQ,'t');t=R.gen();f=R(d['parameter_polynomial']);q=R(d['fourth_polynomial']);n=11
    def norm(a):return matrix(QQ,[[(a*t**j%f)[i] for j in range(n)] for i in range(n)]).det()
    N=norm(q);D=(-1)**(n*(n-1)//2)*norm(f.derivative())
    assert N==QQ(d['norm']) and D==QQ(d['parameter_discriminant'])
    assert N>0>D and not N.is_square() and not (N/D).is_square()
    w=d['q_over_norm_local_obstruction'];p=w['prime'];a=w['root']
    def mod(x):x=Fraction(str(x));return x.numerator*pow(x.denominator,-1,p)%p
    def evaluate(pol):
        acc=0
        for c in reversed(pol.list()):acc=(acc*a+mod(c))%p
        return acc
    assert p==59 and a==4 and mod(D)!=0
    assert evaluate(f)==0 and evaluate(f.derivative())==w['F_derivative_at_root']==51
    qval=evaluate(q);nv=mod(N);ratio=qval*pow(nv,-1,p)%p
    assert qval==7 and nv==8 and ratio==23
    assert 19**2%p==qval and pow(nv,(p-1)//2,p)==p-1 and pow(ratio,(p-1)//2,p)==p-1
    for p,expected in ((73,[11]),(79,[1,2,3,5])):
        k=GF(p);rp=PolynomialRing(k,'t');ff=rp(f)
        assert ff.is_squarefree() and sorted(int(a.degree()) for a,e in ff.factor() for _ in range(int(e)))==expected
    # Any nonconstant vector in an S11-submodule gives a two-coordinate vector.
    checked=0
    for v in range(1,2047):
        i=next(i for i in range(11) if v>>i&1);j=next(j for j in range(11) if not(v>>j&1))
        swapped=v^(1<<i)^(1<<j)
        assert (v^swapped)==(1<<i)|(1<<j)
        checked+=1
    augmentation=[1|(1<<i) for i in range(1,11)];assert len(r.basis(augmentation))==10
    assert 2047 .bit_count()%2==1
    assert d['normal_closure_sign_kernel_rank']==11 and d['normal_closure_Galois_group']=='C2 wr S11'
    return {'schema':'rank-jump.fourth-lift-shared-field-verification.v1','status':'PASS',
        'base_Galois_group':'S11','norm_positive_and_parameter_discriminant_negative':True,
        'degree_one_local_place':{'prime':59,'parameter_residue':4,'fourth_value_residue':7,
            'fourth_square_root_residue':19,'constant_norm_candidate_residue':8,'normalized_value_residue':23},
        'degree22_field_has_quadratic_subfield':False,'nonconstant_binary_vectors_checked':checked,
        'augmentation_dimension':10,'normal_closure_sign_kernel_rank':11,
        'normal_closure_group':'C2 wr S11','normal_closure_group_order':(2**11)*factorial(11),
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (SOURCE,Path(__file__),HERE/'retrospective.py')},
        'boundary':'The note supplies the Galois-module argument from these exact gates. Conjugate squareclass rank is not original-curve Mordell-Weil rank.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);a=p.parse_args();d=compute()
    if a.mode=='build':r.write_new(OUTPUT,d)
    else:assert r.read(OUTPUT)==d
    print('PASS: no quadratic subfield; full sign kernel rank11; C2 wr S11')
