#!/usr/bin/env python3
"""Independent exact section, height-input and branch-obstruction replay.

Uses rational polynomial arithmetic and Rabin irreducibility tests over finite
fields. Sage's factorization declarations and proposed height matrix are not
accepted as proof inputs. The local height and valuation arguments are given
in the canonical note; this is an arithmetic replay, not formal verification.
"""
import argparse
from fractions import Fraction as Q
from hashlib import sha256
import importlib.util
import json
from math import isqrt
from pathlib import Path
import time

ROOT=Path(__file__).resolve().parents[2]
DEFAULT=ROOT/'artifacts/generated-results/elkies-k3-r17-mestre-shared-twist-v1'
SOURCE=ROOT/'artifacts/generated-results/elkies-k3-r17-correlated-genus-one-v1/input.json'
SPEC=importlib.util.spec_from_file_location('mestre_polynomial_replay',ROOT/'elkies-k3/scripts/verify_r17_correlated_genus_one_pencils.py')
V=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(V)
require=V.require


def trim(a):
    a=list(a)
    while a and not a[-1]:a.pop()
    return a


def ff_sub(a,b,p):
    out=list(a)+[0]*max(0,len(b)-len(a))
    for i,c in enumerate(b):out[i]=(out[i]-c)%p
    return trim(out)


def ff_mul(a,b,p):
    if not a or not b:return []
    out=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]=(out[i+j]+x*y)%p
    return trim(out)


def ff_rem(a,b,p):
    a=trim(a);require(b,'zero modular divisor');inv=pow(b[-1],-1,p)
    while len(a)>=len(b):
        shift=len(a)-len(b);c=a[-1]*inv%p
        for j,x in enumerate(b):a[shift+j]=(a[shift+j]-c*x)%p
        a=trim(a)
    return a


def ff_gcd(a,b,p):
    while b:a,b=b,ff_rem(a,b,p)
    return [x*pow(a[-1],-1,p)%p for x in a] if a else []


def ff_power(a,e,f,p):
    out=[1]
    while e:
        if e&1:out=ff_rem(ff_mul(out,a,p),f,p)
        a=ff_rem(ff_mul(a,a,p),f,p);e//=2
    return out


def derivative(f):return V.trim([i*c for i,c in enumerate(f)][1:])


def reduction(f,p):
    require(p>=2 and all(p%d for d in range(2,isqrt(p)+1)),'nonprime finite certificate')
    require(f and all(c.denominator%p for c in f),'nonintegral finite polynomial')
    out=[c.numerator*pow(c.denominator,-1,p)%p for c in f]
    require(out[-1],'degree drop in finite certificate')
    return out


def irreducible(f,p):
    n=len(f)-1
    if n<1:return False
    divisors=[d for d in range(2,n+1) if n%d==0 and all(d%e for e in range(2,isqrt(d)+1))]
    tests={n//d for d in divisors};x=[0,1];z=x
    for i in range(1,n+1):
        z=ff_power(z,p,f,p)
        difference=ff_sub(z,x,p)
        if i in tests and len(ff_gcd(f,difference,p))>1:return False
    return not ff_rem(ff_sub(z,x,p),f,p)


def evaluate(f,x,p):
    out=0
    for c in reversed(f):out=(out*x+c)%p
    return out


def same_fraction(a,b):return V.mul(a[0],b[1])==V.mul(b[0],a[1])


def verify(path):
    packet=json.loads((path/'input.json').read_text());result=json.loads((path/'result.json').read_text());gate=json.loads((path/'branch-gate.json').read_text())
    require(packet['schema']=='r17-mestre-shared-twist-input-v1','input schema')
    source=json.loads(SOURCE.read_text())
    require(packet['source_sha256']==sha256(SOURCE.read_bytes()).hexdigest(),'generic source binding')
    require(packet['A']==source['A'] and packet['B']==source['B'] and packet['u']=='2','generic-only equation and fixed auxiliary selection')
    require(result['input_sha256']==gate['input_sha256']==sha256((path/'input.json').read_bytes()).hexdigest(),'input binding')
    V.source_projection(source)
    A,B=V.poly(packet['A']),V.poly(packet['B'])
    core=V.add(V.scale(V.power(B,2),9261),V.scale(V.power(A,3),400))
    D=V.scale(V.mul(V.mul(A,B),core),-5)
    Delta=V.add(V.scale(V.power(A,3),4),V.scale(V.power(B,2),27))
    require(D==V.poly(result['D']) and [len(A)-1,len(B)-1,len(core)-1,len(D)-1,len(Delta)-1]==[8,12,24,44,24],'branch polynomial identity and degrees')
    p=gate['smooth_disjoint_branch_prime'];d,e=reduction(D,p),reduction(Delta,p)
    for f in (d,e):
        df=trim([i*c%p for i,c in enumerate(f)][1:])
        require(len(ff_gcd(f,df,p))==1,'singular branch or original discriminant')
    require(len(ff_gcd(d,e,p))==1,'branch meets an original singular fibre')
    expected=[((V.scale(B,-21),V.scale(A,5)),([Q(1)],V.scale(V.power(A,2),25))),
              ((V.scale(B,-21),V.scale(A,20)),([Q(1)],V.scale(V.power(A,2),200)))]
    require(len(result['points'])==2,'point count')
    for row,(x,y) in zip(result['points'],expected):
        actual_x=(V.poly(row['x_numerator']),V.poly(row['x_denominator']))
        actual_y=(V.poly(row['y_radical_numerator']),V.poly(row['y_radical_denominator']))
        require(same_fraction(actual_x,x) and same_fraction(actual_y,y),'displayed section formula')
        xn,xd=x;yn,yd=y
        left=V.mul(V.mul(D,V.power(yn,2)),V.power(xd,3))
        right=V.mul(V.power(yd,2),V.add(V.add(V.power(xn,3),V.mul(A,V.mul(xn,V.power(xd,2)))),V.mul(B,V.power(xd,3))))
        require(left==right,'new section Weierstrass identity')
    # x(P+Q)=(-7/15*B^2-20/81*A^3)/(A*B), derived independently
    # from the group law, gives the pole divisor needed for the height.
    plus_n=V.add(V.scale(V.power(B,2),Q(-7,15)),V.scale(V.power(A,3),Q(-20,81)))
    plus_d=V.mul(A,B)
    (xn,xd),(yn,yd)=expected[0]
    slope_n=V.scale(V.mul(V.mul(D,V.power(yn,2)),V.power(xd,2)),49)
    slope_d=V.scale(V.mul(V.power(yd,2),V.power(xn,2)),36)
    sum_n=V.sub(V.mul(slope_n,xd),V.scale(V.mul(xn,slope_d),Q(5,4)))
    sum_d=V.mul(slope_d,xd)
    require(same_fraction((sum_n,sum_d),(plus_n,plus_d)),'addition identity for the height proof')
    require(len(ff_gcd(reduction(plus_n,p),reduction(plus_d,p),p))==1,'cancellation in the sum pole divisor')
    require(len(plus_n)-len(plus_d)<=4 and len(B)-len(A)<=4,'pole at infinity')
    pole_single=len(A)-1;pole_sum=len(A)+len(B)-2
    heights=[8+2*pole_single]*2;cross=Q((8+2*pole_sum)-sum(heights),2)
    require(heights==[24,24] and cross==0,'height Gram')
    fields=[]
    require([r['polynomial'] for r in gate['root_field_certificates']]==['A','B'],'root field coverage')
    for row,f in zip(gate['root_field_certificates'],[A,B]):
        ell=row['irreducible_reduction_prime'];require(irreducible(reduction(f,ell),ell),'root field irreducibility')
        p0=row['degree_one_prime']['prime'];root=row['degree_one_prime']['root'];g=reduction(f,p0)
        dg=trim([i*c%p0 for i,c in enumerate(g)][1:])
        require(evaluate(g,root,p0)==0 and evaluate(dg,root,p0)!=0,'simple degree-one local place')
        require(p0%4==3 and (row['polynomial']!='A' or p0%3==2),'cyclotomic exclusion congruence')
        fields.append({'polynomial':row['polynomial'],'degree':len(f)-1,'irreducibility_prime':ell,'degree_one_prime':p0})
    return {'status':'PASS_TWO_CORRELATED_SECTIONS_AND_MESTRE_GENUS_GATE',
            'quadratic_cover_genus':21,'height_matrix':[[24,0],[0,24]],'new_independent_directions':2,
            'generic_rank_lower_bound_with_inherited_parent':19,'parent_rank_is_inherited_dependency':True,
            'unrestricted_rational_auxiliary_genus_lower_bound':9,'forced_branch_degree':20,
            'root_field_certificates':fields,'goal_complete':False,
            'scope':'The displayed two sections are independent on one genus21 quadratic cover. Every rational auxiliary function in this Mestre identity forces genus at least9 on this fixed parent. The infinite-rational-base requirement remains unmet; other identities and parents are not excluded.'}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--input-dir',type=Path,default=DEFAULT);p.add_argument('--output',type=Path)
    args=p.parse_args();started=time.monotonic();result=verify(args.input_dir);result['elapsed_seconds']=round(time.monotonic()-started,3)
    if args.output:
        with args.output.open('x') as f:json.dump(result,f,indent=2,sort_keys=True);f.write('\n')
    print(json.dumps(result,sort_keys=True))


if __name__=='__main__':main()
