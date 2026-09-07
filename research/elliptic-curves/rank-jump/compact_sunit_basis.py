#!/usr/bin/env python3
"""One bounded compact-unit attempt with a deliberately smaller factor base."""
import argparse
from pathlib import Path
import subprocess
import sys
import time
import retrospective as r

PROTOCOL=Path(__file__).with_name('COMPACT_SUNIT_BASIS_PROTOCOL.json')
PRIOR=r.OUT/'rank_jump_matched103b2_class_boundary_v1.json'
INPUT=r.OUT/'rank_jump_compact_sunit_basis_inputs_v1.json'
OUTPUT=r.OUT/'rank_jump_compact_sunit_basis_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-compact-sunit-basis-v1'


def binding(paths):return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in paths}


def export():
    cases=[]
    for row in r.read(PRIOR)['rows']:
        red=row['reduction'];assert red['status']=='PASS' and row['boundary']['status']=='PASS'
        cases.append({'token':row['token'],'cubic':red['reduced_cubic_ascending'],
            'basis':red['transported_maximal_order_basis'],'field_discriminant':red['field_discriminant'],
            'S_finite':row['boundary']['S_finite']})
    assert [c['token'] for c in cases]==r.read(PROTOCOL)['cases']
    r.write_new(INPUT,{'schema':'rank-jump.compact-sunit-basis-inputs.v1','cases':cases,
        'bindings':binding([PROTOCOL,PRIOR,Path(__file__)])})


def worker(token):
    from sage.all import QQ,PolynomialRing,pari
    pari.allocatemem(64000000,268435456,silent=True);pari.setrand(20260907)
    R=PolynomialRing(QQ,'z');tech=pari('[0.03,0.03,-1]')
    if token=='control':
        nf=pari.nfinit(pari(R([-11,0,0,1])));ps=[2,3,11]
    else:
        c=next(c for c in r.read(INPUT)['cases'] if c['token']==token)
        nf=pari.nfinit([pari(R(c['cubic'])),[pari(R(b)) for b in c['basis']]])
        assert str(nf.disc())==c['field_discriminant'];ps=c['S_finite']
    primes=[P for p in ps for P in pari.idealprimedec(nf,p)]
    r.write_new(WORK/f'{token}-setup.json',{'status':'PASS','field_discriminant':str(nf.disc()),
        'S_prime_ideals':len(primes),'signature':list(map(int,nf.nf_get_sign())),
        'pari_version':str(pari.version()),'tech':str(tech)})
    print('SETUP_COMPLETE',token,flush=True)
    pari.default('debug',1)
    try:
        bnf=pari.bnfinit(nf,1,tech)
        r.write_new(WORK/f'{token}-bnf.json',{'status':'RETURNED_UNCERTIFIED',
            'cyclic_factors':list(map(str,bnf.bnf_get_cyc())),'factor_base_size':len(bnf[4]),
            'compact_bnf':str(bnf)})
        print('BNF_RETURNED',token,flush=True)
        units=pari.bnfunits(bnf,primes)
        r.write_new(WORK/f'{token}-units.json',{'status':'RETURNED_UNCERTIFIED','compact_units':str(units),
            'generator_count':len(units[0]),'prime_count':len(primes)})
        print('UNITS_RETURNED',token,flush=True)
        certified=int(pari.bnfcertify(bnf))==1
        result={'status':'CERTIFIED' if certified else 'UNKNOWN','certified':certified,'generator_count':len(units[0])}
    except Exception as exc:
        result={'status':'UNKNOWN','exception_type':type(exc).__name__,'reason':str(exc)}
    r.write_new(WORK/f'{token}-terminal.json',result)


def capture():
    WORK.mkdir(parents=True,exist_ok=True);rows=[]
    for token in ['control']+r.read(PROTOCOL)['cases']:
        terminal=WORK/f'{token}-terminal.json';log=WORK/f'{token}.log'
        if not terminal.exists():
            start=time.monotonic()
            with log.open('x') as stream:
                try:
                    proc=subprocess.run([sys.executable,str(Path(__file__).resolve()),'worker','--token',token],stdout=stream,stderr=stream,timeout=30)
                    error=None if proc.returncode==0 else 'worker failure'
                except subprocess.TimeoutExpired:error='30-second timeout'
            if error:r.write_new(terminal,{'status':'UNKNOWN','reason':error})
            r.write_new(WORK/f'{token}-execution.json',{'elapsed_seconds':time.monotonic()-start,'process_error':error})
        result=r.read(terminal)
        stages={stage:r.read(WORK/f'{token}-{stage}.json') for stage in ('setup','bnf','units') if (WORK/f'{token}-{stage}.json').exists()}
        rows.append({'token':token,'terminal':result,'stages':stages,'execution':r.read(WORK/f'{token}-execution.json'),
            'log':log.read_text(),'log_sha256':r.digest(log.read_bytes())})
        print(token,result['status'],result.get('reason',''),flush=True)
        if token=='control' and result['status']!='CERTIFIED':break
    r.write_new(OUTPUT,{'schema':'rank-jump.compact-sunit-basis.v1','rows':rows,
        'bindings':binding([Path(__file__),PROTOCOL,INPUT]),'boundary':r.read(PROTOCOL)['boundary']})


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('mode',choices=['export','capture','worker']);parser.add_argument('--token');args=parser.parse_args()
    if args.mode=='worker':worker(args.token)
    else:globals()[args.mode]()
