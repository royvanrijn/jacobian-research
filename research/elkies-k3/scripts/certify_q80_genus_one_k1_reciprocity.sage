#!/usr/bin/env sage
"""Source-bound finite gates for the written Q80 genus-one k1 exclusion."""
from sage.all import GF, PolynomialRing, QQ
from collections import Counter
from hashlib import sha256
from pathlib import Path
import argparse
import itertools
import json
import resource
import subprocess
import tempfile
import time

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'artifacts/generated-results/elkies-k3-q80-genus-one-k1-reciprocity-v1'
SOURCE='artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json'
DEGREE='artifacts/generated-results/elkies-k3-q80-degree-two-reciprocity-v1'
BOUNDARY='artifacts/generated-results/elkies-k3-q80-genus-one-two-disagreement-v1'
CPP='elkies-k3/scripts/certify_q80_cubic_norm_census.cpp'
REPLAY='elkies-k3/scripts/replay_q80_cubic_norm_census.cpp'
CHECKER='elkies-k3/scripts/verify_q80_genus_one_k1_reciprocity.py'
TEST='tests/test_q80_genus_one_k1_reciprocity.py'
P=131
def read(p):return json.loads(p.read_text())
def digest(p):return sha256(p.read_bytes()).hexdigest()
def write(name,data):
    with (OUT/name).open('x') as f:json.dump(data,f,indent=2,sort_keys=True);f.write('\n')
def bindings():
    return [SOURCE,CPP,REPLAY,CHECKER,TEST,str(Path(__file__).relative_to(ROOT)),
            *[DEGREE+'/'+n for n in ['input.json','result.json','rational-fibres.json','quadratic-fibres.json','independent-replay.json']],
            *[BOUNDARY+'/'+n for n in ['input.json','result.json','character-gate.json','independent-replay.json']]]
def freeze():
    packet={'schema':1,'prime':P,'cubic_modulus_low_to_high':[3,1,0,1],
            'cubic_orbits':(P**3-P)//3,'factorization_trial_limit':64,'driver_cpu_seconds':30,
            'census_cpu_seconds':90,'driver_memory_bytes':2*1024**3,'census_memory_bytes':512*1024**2,
            'bindings':{p:digest(ROOT/p) for p in bindings()},
            'preserved_preflight':{str(p.relative_to(ROOT)):digest(p) for p in (OUT/'preflight').rglob('*') if p.is_file()},
            'scope':'Complete rational/quadratic/cubic norm gate; polynomial k1 integrality and the unbounded exclusion are written arguments.'}
    write('input.json',packet);print(json.dumps({'status':'FROZEN','input_sha256':digest(OUT/'input.json')}),flush=True)
def run():
    resource.setrlimit(resource.RLIMIT_CPU,(30,95));resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    start=time.process_time();wall=time.monotonic();packet=read(OUT/'input.json')
    assert set(packet['bindings'])==set(bindings())
    for p,h in {**packet['bindings'],**packet['preserved_preflight']}.items():assert digest(ROOT/p)==h
    source=read(ROOT/SOURCE);F=GF(P);R=PolynomialRing(F,'t');t=R.gen()
    poly=lambda v:R([F(QQ(a)) for a in v])
    A=poly(source['weierstrass_model']['A_coefficients_low_to_high']);B=poly(source['weierstrass_model']['B_coefficients_low_to_high'])
    delta=4*A**3+27*B**2;factors=delta.factor()
    assert [(f.degree(),e) for f,e in factors]==[(2,1),(3,1),(19,1)]
    assert all(delta(b) for b in F) and 4*A[8]**3+27*B[12]**2
    assert R([3,1,0,1]).is_irreducible()
    write('fibre-hypotheses.json',{'discriminant_factors':[[[int(v) for v in f.list()],int(e)] for f,e in factors],
          'all_projective_rational_fibres_smooth':True,'cubic_modulus_irreducible':True})
    rat=read(ROOT/DEGREE/'rational-fibres.json')['rows'];quad=read(ROOT/DEGREE/'quadratic-fibres.json')['rows']
    agreement={}
    for r in rat:
        for e,c in zip(r['rational_roots'],r['rational_codes']):
            assert c and c not in agreement;agreement[c]=(r['t'],e)
    assert len(agreement)==110
    split=[r for r in rat if len(r['rational_roots'])==3];assert len(split)==16
    counts=Counter();candidates=[]
    for rows in itertools.combinations(split,3):
        for choices in itertools.product(*[list(itertools.permutations(zip(r['rational_roots'],r['rational_codes']))) for r in rows]):
            counts['rational_triple_choices']+=1;codes=[0,0,0]
            for cc in choices:
                for j,(_,c) in enumerate(cc):codes[j]^=c
            if codes[2]:continue
            counts['rational_triples_unused_zero']+=1;assert codes[0]==codes[1]
            if codes[0] in agreement and agreement[codes[0]][0] not in [r['t'] for r in rows]:
                candidates.append({'type':'111','agreement':agreement[codes[0]],'d_bases':[r['t'] for r in rows]})
    for r in split:
        for s in quad:
            if len(s['roots'])!=3:continue
            for rr in itertools.permutations(zip(r['rational_roots'],r['rational_codes'])):
                for ss in set(itertools.permutations(zip(s['roots'],s['norm_codes']))):
                    counts['rational_quadratic_choices']+=1
                    codes=[rr[j][1]^ss[j][1] for j in range(3)]
                    if codes[2]:continue
                    counts['rational_quadratic_unused_zero']+=1;assert codes[0]==codes[1]
                    if codes[0] in agreement and agreement[codes[0]][0]!=r['t']:
                        candidates.append({'type':'12','agreement':agreement[codes[0]],'d_bases':[r['t'],s['t']]})
    assert counts['rational_triple_choices']==120960 and counts['rational_quadratic_choices']==838944
    assert not candidates
    write('composite-norms.json',{'counts':dict(counts),'candidate_rows':candidates,
          'agreement_points':110,'full_split_rational_sites':[r['t'] for r in split]})
    lines=[]
    def emit(f):
        v=[int(c) for c in f.list()];lines.append(str(len(v))+' '+' '.join(map(str,v)))
    emit(A);emit(B);lines.append('17')
    for section in source['sections']['records']:
        x=section['X'];num=poly(x['numerator_coefficients_low_to_high']);den=poly(x['denominator_coefficients_low_to_high'])
        g=num.gcd(den);emit(num//g);emit(den//g)
    lines.append('110')
    for c,(b,e) in agreement.items():lines.append(str(c)+' '+str(-2 if b is None else b))
    text_input='\n'.join(lines)+'\n'
    with (OUT/'cubic-input.txt').open('x') as f:f.write(text_input)
    with tempfile.TemporaryDirectory(prefix='q80-k1-cubic-') as tmp:
        exe=Path(tmp)/'census'
        subprocess.run(['g++','-O3','-std=c++17',str(ROOT/CPP),'-o',str(exe)],check=True,timeout=30)
        run=subprocess.run([str(exe)],input=text_input,text=True,capture_output=True,check=True,timeout=95)
    cubic=json.loads(run.stdout)
    assert cubic['status']=='COMPLETE' and cubic['orbits']==749320
    assert (cubic['smooth_split'],cubic['nodal'])==(125595,1)
    assert len(cubic['zero_norm_rows'])==8 and not cubic['candidate_rows']
    # Check every retained special row again using Sage's extension field.
    K=GF(P**3,name='theta',modulus=R([3,1,0,1]));theta=K.gen()
    decode=lambda n:K(n%P)+K(n//P%P)*theta+K(n//(P*P))*theta**2
    local=[]
    for row in cubic['zero_norm_rows']:
        tau=decode(row['t']);assert tau**P!=tau
        roots=[decode(e) for e in row['roots']]
        assert all(e**3+A(tau)*e+B(tau)==0 for e in roots)
        assert sum(roots)==0 and sum(roots[i]*roots[j] for i in range(3) for j in range(i+1,3))==A(tau)
        assert -roots[0]*roots[1]*roots[2]==B(tau)
        codes=[]
        for e in roots:
            code=0
            for j,section in enumerate(source['sections']['records']):
                x=section['X'];num=poly(x['numerator_coefficients_low_to_high']);den=poly(x['denominator_coefficients_low_to_high'])
                g=num.gcd(den);num//=g;den//=g
                if not den(tau):assert num(tau);continue
                value=num(tau)/den(tau)-e
                if not value:value=3*e*e+A(tau)
                assert value
                code |= int(not value.is_square())<<j
            codes.append(code)
        assert codes==row['norm_codes']
        assert all(c==0 or c not in agreement for c in codes)
        local.append({'t':row['t'],'norm_codes':codes})
    write('cubic-census.json',cubic)
    write('cubic-row-checks.json',{'sage_checked_rows':local})
    result={'status':'PASS','input_sha256':digest(OUT/'input.json'),
            'records':{p:digest(OUT/p) for p in ['fibre-hypotheses.json','composite-norms.json','cubic-input.txt','cubic-census.json','cubic-row-checks.json']},
            'finite_candidate_count':0,'written_integrality_argument_required':True,
            'excluded_allocations':['k1,j0','k1,j1','k1,j2','k1,j3'],
            'genus1_allocations_remaining':5,'positive_mw17_target_complete':False,
            'formal_verification':False,'old_norm8_replay_upgraded':False}
    write('result.json',result)
    write('execution.json',{'status':'PASS','input_sha256':digest(OUT/'input.json'),'result_sha256':digest(OUT/'result.json'),
          'sage_cpu_seconds':time.process_time()-start,'census_cpu_seconds':cubic['cpu_seconds'],'wall_seconds':time.monotonic()-wall,
          'driver_cpu_limit':30,'census_cpu_limit':90,'driver_memory_bytes':2*1024**3,'census_memory_bytes':512*1024**2})
    print(json.dumps(result,sort_keys=True),flush=True)
if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('action',choices=['freeze','run']);args=parser.parse_args()
    freeze() if args.action=='freeze' else run()
