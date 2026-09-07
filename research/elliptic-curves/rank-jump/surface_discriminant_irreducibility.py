#!/usr/bin/env python3
"""Bounded modular completion for the original surface discriminant."""
import argparse
from math import gcd
from pathlib import Path
import subprocess
import retrospective as r
import parameter_cover_capacity as source
import verify_parameter_cover_capacity as previous

PROTOCOL=Path(__file__).with_name('SURFACE_DISCRIMINANT_IRREDUCIBILITY_PROTOCOL.json')
OUTPUT=r.OUT/'rank_jump_surface_discriminant_irreducibility_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-surface-discriminant-irreducibility-v1'


def compute():
    from sage.all import QQ,GF,PolynomialRing
    assert r.read(previous.OUTPUT)['independent_irreducibility_status']=='UNKNOWN'
    coeff=list(map(int,r.read(source.OUTPUT)['discriminant_coefficients_ascending']))
    content=gcd(*coeff);coeff=[c//content for c in coeff]
    R=PolynomialRing(QQ,'u');D=R(coeff);spec=r.read(PROTOCOL)['limits']
    candidates=set(range(1,24));rows=[]
    for p in r.primes(spec['largest_prime']):
        if p<spec['smallest_prime']:continue
        f=D.change_ring(GF(p))
        if f.degree()!=24 or f.gcd(f.derivative())!=1:
            rows.append({'prime':p,'status':'SKIP_BAD_REDUCTION'});continue
        fac=f.factor();degrees=[int(g.degree()) for g,e in fac]
        assert all(e==1 for g,e in fac)
        candidates &= previous.possible_degrees(degrees)
        rows.append({'prime':p,'factor_degrees':degrees,
            'monic_factor_coefficients_ascending':[[int(c) for c in g.monic()] for g,e in fac],
            'remaining_proper_factor_degrees':sorted(candidates)})
        if not candidates:break
    return {'schema':'rank-jump.surface-discriminant-irreducibility.v1',
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (Path(__file__),PROTOCOL,source.OUTPUT,previous.OUTPUT)},
        'integer_content':str(content),'primitive_coefficients_ascending':list(map(str,coeff)),
        'status':'PASS_IRREDUCIBLE' if not candidates else 'UNKNOWN',
        'remaining_proper_factor_degrees':sorted(candidates),'modular_certificate':rows}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['capture','worker','check']);a=p.parse_args()
    if a.mode=='worker':r.write_new(WORK/'checkpoint.json',compute())
    elif a.mode=='check':assert r.read(OUTPUT)==compute();print('PASS surface discriminant modular replay')
    else:
        WORK.mkdir(parents=True,exist_ok=True);path=WORK/'checkpoint.json'
        if not path.exists():
            with (WORK/'worker.log').open('x') as log:
                try:
                    proc=subprocess.run(['sage','-python',str(Path(__file__).resolve()),'worker'],
                        cwd=r.ROOT,stdout=log,stderr=log,timeout=30)
                    reason=None if proc.returncode==0 else 'worker failure'
                except subprocess.TimeoutExpired:reason='30-second timeout'
                if reason and not path.exists():r.write_new(path,{'status':'UNKNOWN','reason':reason})
        data=r.read(path);r.write_new(OUTPUT,data);print(data['status'],data.get('modular_certificate',[]))
