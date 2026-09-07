#!/usr/bin/env python3
"""Bounded original-fibre secant squareclasses; no exceptional-point input."""
import argparse
from itertools import combinations
from pathlib import Path
import subprocess
import sys
import retrospective as r
from retrospective_secant_pencils import rational_square

PROTOCOL=Path(__file__).with_name('ORIGINAL_SECANT_CLASSES_PROTOCOL.json')
INPUT=r.OUT/'rank_jump_original_secant_class_inputs_v1.json'
OUTPUT=r.OUT/'rank_jump_original_secant_classes_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-original-secant-classes-v1'


def export():
    import bad_prime_support as bad
    rows=[]
    for i,source in enumerate(bad.cases()):
        model,points=r.short(source['model'],source['generic_points'])
        rows.append({'case_index':i,'short_model':model,'generic_points':points})
    r.write_new(INPUT,{'schema':'rank-jump.original-secant-class-inputs.v1',
        'source_sha256':r.digest(r.INPUT.read_bytes()),'rows':rows})


def bindings():
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in
        (Path(__file__),PROTOCOL,INPUT,Path(__file__).with_name('retrospective_secant_pencils.py'))}


def squareclass_signature(value,primes):
    """Real sign plus complete odd-p squareclass (valuation parity, unit bit)."""
    value=r.F(value);assert value
    bits=[int(value<0)]
    for p in primes:
        n,d=value.numerator,value.denominator;v=0
        while n%p==0:n//=p;v+=1
        while d%p==0:d//=p;v-=1
        unit=n*pow(d,-1,p)%p
        bits.extend([v%2,int(pow(unit,(p-1)//2,p)==p-1)])
    return r.pack(bits)


def intercept(A,B,P,Q,sign):
    a,p=map(r.F,P);b,q=map(r.F,Q);q*=sign
    assert p*p==a**3+A*a+B and q*q==b**3+A*b+B
    if a==b:return {'status':'DEGENERATE','reason':'equal abscissae'}
    slope=(q-p)/(b-a)
    if not slope:return {'status':'DEGENERATE','reason':'horizontal secant'}
    x0=a-p/slope;C=x0**3+A*x0+B
    if not C:return {'status':'RATIONAL_TWO_TORSION','x0':str(x0),'C':'0'}
    return {'status':'PASS','x0':str(x0),'C':str(C),'rationally_soluble':rational_square(C)}


def compute(index):
    old=next(x for x in r.read(INPUT)['rows'] if x['case_index']==index)
    A,B=map(r.F,old['short_model'][3:]);points=old['generic_points']
    primes=r.read(PROTOCOL)['limits']['fingerprint_primes']
    rows=[];groups=[];buckets={};ratio_tests=0
    for i,j in combinations(range(len(points)),2):
        for sign in [1,-1]:
            row={'pair':[i,j],'relative_sign':sign,**intercept(A,B,points[i],points[j],sign)}
            if row['status']=='PASS':
                C=r.F(row['C']);sig=squareclass_signature(C,primes);row['local_signature']=sig
                found=None
                for group in buckets.get(sig,[]):
                    ratio_tests+=1
                    if rational_square(C/r.F(group['representative_C'])):found=group;break
                if found is None:
                    found={'class_index':len(groups),'representative_C':str(C),'local_signature':sig,'members':[]}
                    groups.append(found);buckets.setdefault(sig,[]).append(found)
                row['class_index']=found['class_index'];found['members'].append(len(rows))
            rows.append(row)
    return {'bindings':bindings(),'case_index':index,'status':'PASS','rows':rows,'classes':groups,
        'summary':{'secants':len(rows),'nondegenerate':sum(x['status']=='PASS' for x in rows),
            'rationally_soluble':sum(x.get('rationally_soluble',False) for x in rows),
            'rational_two_torsion_hits':sum(x['status']=='RATIONAL_TWO_TORSION' for x in rows),
            'distinct_quadratic_classes':len(groups),'repeated_classes':sum(len(x['members'])>1 for x in groups),
            'largest_class_multiplicity':max([len(x['members']) for x in groups],default=0),
            'exact_ratio_tests_after_local_filter':ratio_tests}}


def capture(check=False):
    WORK.mkdir(parents=True,exist_ok=True);rows=[]
    for i in r.read(PROTOCOL)['cases']:
        if check:
            assert compute(i)==next(x for x in r.read(OUTPUT)['rows'] if x['case_index']==i)
            print('PASS original secant replay',i,flush=True);continue
        path=WORK/f'case-{i}.json'
        if not path.exists():
            with (WORK/f'case-{i}.log').open('x') as log:
                try:
                    proc=subprocess.run([sys.executable,str(Path(__file__).resolve()),'worker','--index',str(i),
                        '--destination',str(path)],cwd=r.ROOT,stdout=log,stderr=log,timeout=30)
                    failure=None if proc.returncode==0 else 'worker failure'
                except subprocess.TimeoutExpired:failure='30-second timeout'
                if failure and not path.exists():r.write_new(path,{'bindings':bindings(),'case_index':i,'status':'UNKNOWN','reason':failure})
        row=r.read(path);assert row['bindings']==bindings();rows.append(row)
        print('checkpoint',i,row.get('summary',row['status']),flush=True)
    if not check:
        assert sum(x.get('summary',{}).get('secants',0) for x in rows)<=r.read(PROTOCOL)['limits']['secants_total']
        r.write_new(OUTPUT,{'schema':'rank-jump.original-secant-classes.v1','bindings':bindings(),'rows':rows})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['export','capture','worker','check'])
    p.add_argument('--index',type=int);p.add_argument('--destination',type=Path);args=p.parse_args()
    if args.mode=='export':export()
    elif args.mode=='worker':r.write_new(args.destination,compute(args.index))
    else:capture(args.mode=='check')
