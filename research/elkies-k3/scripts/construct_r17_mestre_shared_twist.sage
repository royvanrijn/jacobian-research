#!/usr/bin/env sage-python
"""Exact prospective Mestre two-point identity and covering-genus gate.

The auxiliary value2 is the first positive integer after the degenerate1.
No exceptional point or specialized rank is read. Factorizations are small
generic polynomial calculations; the worker has a30-second CPU cap.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import resource
import time

from sage.all import PolynomialRing, QQ, GF, prime_range

ROOT=Path(__file__).resolve().parents[2]
SOURCE=ROOT/'artifacts/generated-results/elkies-k3-r17-correlated-genus-one-v1/input.json'
DEFAULT=ROOT/'artifacts/generated-results/elkies-k3-r17-mestre-shared-twist-v1'


def coeffs(p):return [str(x) for x in p.list()]


def write_new(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('x') as f:json.dump(data,f,indent=2,sort_keys=True);f.write('\n')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=DEFAULT);args=p.parse_args()
    resource.setrlimit(resource.RLIMIT_CPU,(30,35));resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3))
    started=time.process_time();source=json.loads(SOURCE.read_text())
    packet={'schema':'r17-mestre-shared-twist-input-v1','source_sha256':sha256(SOURCE.read_bytes()).hexdigest(),
            'A':source['A'],'B':source['B'],'u':'2','limits':{'cpu_seconds':30,'address_space_gib':4},
            'selection':'The first nondegenerate positive integer auxiliary parameter, using only the published generic equation. The branch obstruction is then tested for arbitrary rational auxiliary functions.'}
    write_new(args.output/'input.json',packet)
    R=PolynomialRing(QQ,'t');t=R.gen();K=R.fraction_field()
    A,B=R(packet['A']),R(packet['B']);u=QQ(packet['u']);v=u*u;k=v*v+v+1
    core=B*B*k**3+A**3*v*v*(v+1)**2
    D=-A*B*(v+1)*core
    x1=-K(B)*k/(A*(v+1));x2=x1/v
    z1=K(1)/(A*A*(v+1)**2);z2=z1/u**3
    assert D*z1*z1==x1**3+A*x1+B
    assert D*z2*z2==x2**3+A*x2+B
    assert x1!=x2 and z1!=z2
    delta=4*A**3+27*B*B
    smooth=D.gcd(D.derivative()).degree()==0
    branch_good=D.gcd(delta).degree()==0
    rows={}
    for name,f in [('A',A),('B',B),('core',core),('D',D)]:
        rows[name]={'degree':int(f.degree()),'squarefree':bool(f.gcd(f.derivative()).degree()==0),
                    'factors':[{'coefficients':coeffs(g),'multiplicity':int(e)} for g,e in f.factor()]}
    assert smooth and branch_good and rows['A']['degree']==8 and rows['B']['degree']==12
    # A simple degree-one reduction at p=2 mod3 rules out zeta3 in
    # an irreducible characteristic-zero root field. Retain the full
    # irreducible-factor record; the independent check must certify it.
    witnesses=[]
    for factor in rows['A']['factors']:
        f=R(factor['coefficients'])
        for ell in prime_range(5,200):
            if ell%3!=2 or f.denominator()%ell==0:continue
            F=PolynomialRing(GF(ell),'t')(f)
            if F.degree()!=f.degree() or F.gcd(F.derivative()).degree()!=0:continue
            roots=F.roots(multiplicities=False)
            if roots:
                witnesses.append({'factor':coeffs(f),'prime':int(ell),'root':int(roots[0])});break
        else:witnesses.append({'factor':coeffs(f),'status':'NO_WITNESS_IN_FIXED_PRIME_LIST'})
    point=lambda x,z:{'x_numerator':coeffs(x.numerator()),'x_denominator':coeffs(x.denominator()),
                       'y_radical_numerator':coeffs(z.numerator()),'y_radical_denominator':coeffs(z.denominator())}
    result={'schema':'r17-mestre-shared-twist-result-v1','input_sha256':sha256((args.output/'input.json').read_bytes()).hexdigest(),
            'D':coeffs(D),'points':[point(x1,z1),point(x2,z2)],'factorizations':rows,
            'branch_disjoint_from_surface_discriminant':bool(branch_good),
            'cover_genus':int((D.degree()-2)//2),'A_root_field_witnesses':witnesses,
            'height_certificate_inputs':{'zero_intersections':[8,8],'mutual_intersection_at_A':8,'mutual_intersection_at_B':12,
                    'height_matrix':[[24,0],[0,24]],'requires_written_local_intersection_proof':True},
            'cpu_seconds':time.process_time()-started,
            'goal_status':'NOT_COMPLETE: the displayed cover has genus21 and cannot supply infinitely many rational base points.'}
    write_new(args.output/'result.json',result)
    print(json.dumps({'cover_genus':result['cover_genus'],'factor_degrees':{n:[len(g['coefficients'])-1 for g in z['factors']] for n,z in rows.items()},
                      'A_root_field_witnesses':witnesses,'cpu_seconds':result['cpu_seconds']}),flush=True)


if __name__=='__main__':main()
