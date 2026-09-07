#!/usr/bin/env python3
"""Four fixed rigorous first-derivative intervals; no elliptic point search."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import retrospective as r

HERE=Path(__file__).resolve().parent
PROTOCOL=HERE/'ZERO_CT_ANALYTIC_GATE_PROTOCOL.json'
PRIOR=r.OUT/'rank_jump_full_small_governing_block_v1.json'
OUTPUT=r.OUT/'rank_jump_zero_ct_analytic_gate_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-zero-ct-analytic-gate-v1'


def bindings():
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (Path(__file__),PROTOCOL,PRIOR,Path(r.__file__))}


def worker(index):
    from sage.all import QQ,ZZ,EllipticCurve,RealBallField,pari
    from sage.version import version
    spec=r.read(PROTOCOL);d=spec['twists'][index]
    curve=[0,-11*d,0,-14*d*d,-d**3];E=EllipticCurve(QQ,curve);minimal=E.minimal_model()
    conductor=ZZ(E.conductor());square=conductor.is_square();assert square
    root=int(E.root_number());assert root==-1
    factors=list(conductor.factor());local=[]
    for p,e in factors:
        data=minimal.local_data(p)
        local.append({'prime':int(p),'conductor_exponent':int(e),
                      'kodaira':str(data.kodaira_symbol()),'tamagawa_number':int(data.tamagawa_number()),
                      'local_root_number':int(pari.ellrootno(pari(E),p))})
    assert -__import__('math').prod(x['local_root_number'] for x in local)==root
    K=int(spec['limits']['terms_per_sqrt_conductor']*conductor.sqrt());assert K<=spec['limits']['maximum_terms']
    coefficients=list(map(int,E.anlist(K,python_ints=True)))
    assert len(coefficients)==K+1 and coefficients[0]==0 and coefficients[1]==1
    coefficient_bytes=json.dumps(coefficients,separators=(',',':')).encode()
    R=RealBallField(spec['limits']['precision_bits']);c=2*R.pi()/R(conductor).sqrt();zero=R(0)
    total=R(0)
    for n in range(1,K+1):
        if coefficients[n]:total+=2*R(coefficients[n])/n*zero.gamma_inc(c*n)
    tail=4*(-c*(K+1)).exp()/(c*(K+1)*(1-(-c).exp()))
    enclosure=total.add_error(tail.upper())
    enc=lambda x:{'lower':str(x.lower().exact_rational()),'upper':str(x.upper().exact_rational()),'display':str(x)}
    interval=enc(enclosure);nonzero=QQ(interval['lower'])>0 or QQ(interval['upper'])<0
    return {'schema':'rank-jump.zero-ct-analytic-case.v1','bindings':bindings(),'twist':d,'role':spec['roles'][index],
            'model':curve,'minimal_model':list(map(int,minimal.ainvs())),'conductor':str(conductor),'local_data':local,
            'root_number':root,'terms':K,'coefficient_sha256':r.digest(coefficient_bytes),
            'first_coefficients':coefficients[:25],'partial_sum':enc(total),'tail_bound':enc(tail),'L_derivative_interval':interval,
            'nonzero_derivative_certified':nonzero,'analytic_rank':1 if nonzero else 'UNKNOWN',
            'algebraic_rank':1 if nonzero else 'UNKNOWN','Sha_finite':True if nonzero else 'UNKNOWN',
            'rank_theorem':'Modularity and Gross-Zagier/Kolyvagin when root number -1 and the derivative interval excludes zero.',
            'software':{'sage':version,'pari':str(pari.version())},
            'boundary':'No rational point is constructed. Containing zero is inconclusive. The coefficient and local arithmetic use exact CAS computations; the sum and tail include directed interval error.'}


def capture():
    WORK.mkdir(parents=True,exist_ok=True);rows=[]
    for i,d in enumerate(r.read(PROTOCOL)['twists']):
        path=WORK/f'twist-{d}.json'
        if not path.exists():
            with (WORK/f'twist-{d}.log').open('x') as log:
                try:
                    process=subprocess.run([sys.executable,str(Path(__file__)),'worker','--index',str(i)],
                                           stdout=log,stderr=log,timeout=r.read(PROTOCOL)['limits']['seconds_per_worker'])
                    error=None if process.returncode==0 else 'worker failed'
                except subprocess.TimeoutExpired:error='bounded worker timed out'
            if error:raise RuntimeError(f'{error}; inspect {WORK}/twist-{d}.log')
        row=r.read(path);assert row['bindings']==bindings();rows.append(row)
        print('checkpoint',d,row['L_derivative_interval']['display'],row['analytic_rank'],flush=True)
    r.write_new(OUTPUT,{'schema':'rank-jump.zero-ct-analytic-gate.v1','status':'PASS','bindings':bindings(),'rows':rows,
                        'boundary':'Rank-one conclusions only where a rigorous derivative interval excludes zero. No point search, new twist population, or numerical higher-rank claim.'})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['capture','worker']);p.add_argument('--index',type=int);args=p.parse_args()
    if args.mode=='worker':
        d=r.read(PROTOCOL)['twists'][args.index];r.write_new(WORK/f'twist-{d}.json',worker(args.index))
    else:capture()
