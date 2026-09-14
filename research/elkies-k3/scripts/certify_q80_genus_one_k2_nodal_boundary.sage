#!/usr/bin/env sage
"""Frozen arithmetic gates for Q80 genus1 k2; the nodal denominator boundary stays open."""
from sage.all import GF, Integers, PolynomialRing, QQ, matrix, vector
from pathlib import Path
from hashlib import sha256
from collections import defaultdict
import argparse
import itertools
import json
import resource
import subprocess
import tempfile
import time

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'artifacts/generated-results/elkies-k3-q80-genus-one-two-disagreement-v1'
SOURCE='artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json'
DEGREE='artifacts/generated-results/elkies-k3-q80-degree-two-reciprocity-v1'
NORM4='artifacts/generated-results/elkies-k3-q80-genus-zero-integral-reduction-v1'
CHECKER='elkies-k3/scripts/verify_q80_genus_one_k2_nodal_boundary.py'
CPP='elkies-k3/scripts/certify_q80_nodal_quadratic_census.cpp'
OLD_CHECKER='elkies-k3/scripts/verify_q80_genus_zero_closure.py'
def read(path):return json.loads(path.read_text())
def digest(path):return sha256(path.read_bytes()).hexdigest()
def write(name,data):
    with (OUT/name).open('x') as f:json.dump(data,f,indent=2,sort_keys=True);f.write('\n')
def files():
    return [SOURCE,CHECKER,CPP,OLD_CHECKER,str(Path(__file__).relative_to(ROOT)),
            *[DEGREE+'/'+n for n in ['input.json','result.json','rational-fibres.json','quadratic-fibres.json','independent-replay.json']],
            *[NORM4+'/'+n for n in ['input.json','result.json','norm4-sections.json','complete-replay-input.json','complete-independent-replay.json']]]
def freeze():
    preview=[str(p.relative_to(ROOT)) for p in OUT.rglob('*') if p.is_file()]
    packet={'schema':1,'prime':131,'cpu_seconds':40,'memory_bytes':4*1024**3,
            'nodal_census_count':131**3,'etale_coefficient_count':131,'split_contact_pairs':21,
            'bindings':{s:digest(ROOT/s) for s in files()},
            'preserved_preflight':{s:digest(ROOT/s) for s in preview},
            'scope':'Complete integral genus1 k2 exclusion and necessary nodal denominator reduction; no remaining-case emptiness or positive MW17 cover.'}
    write('input.json',packet)
    print(json.dumps({'status':'FROZEN','input_sha256':digest(OUT/'input.json')}),flush=True)
def run():
    packet=read(OUT/'input.json')
    assert set(packet['bindings'])==set(files())
    for s,h in packet['bindings'].items():assert digest(ROOT/s)==h
    resource.setrlimit(resource.RLIMIT_CPU,(40,45));resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3))
    start=time.process_time();wall=time.monotonic();p=131
    source=read(ROOT/SOURCE);F=GF(p);R=PolynomialRing(F,'t');t=R.gen()
    A=R([F(QQ(v)) for v in source['weierstrass_model']['A_coefficients_low_to_high']])
    B=R([F(QQ(v)) for v in source['weierstrass_model']['B_coefficients_low_to_high']])
    discriminant=4*A**3+27*B**2;factors=discriminant.factor()
    assert [(int(q.degree()),int(e)) for q,e in factors]==[(2,1),(3,1),(19,1)]
    q=t*t+62*t+88;assert factors[0][0]==q
    node=F(38)/31*t+F(122)-F(38)/31*100
    simple=F(55)/31*t+F(18)-F(55)/31*100
    assert (node**3+A*node+B)%q==0 and (3*node**2+A)%q==0
    assert (simple**3+A*simple+B)%q==0 and (3*simple**2+A)%q!=0
    rational=read(ROOT/DEGREE/'rational-fibres.json')['rows']
    quadratic=read(ROOT/DEGREE/'quadratic-fibres.json')['rows']
    root_points=[(r['t'],e,c) for r in rational for e,c in zip(r['rational_roots'],r['rational_codes'])]
    assert len(root_points)==len({r[2] for r in root_points})==110
    assert all(r[2] for r in root_points)
    zero_rows=[r for r in quadratic if 0 in r['norm_codes']]
    assert len(zero_rows)==1 and zero_rows[0]['t']==4161
    pairs13412=[[list(a[:2]),list(b[:2])] for a,b in itertools.combinations_with_replacement(root_points,2) if a[2]^b[2]==13412]
    norm13412=[{'t':r['t'],'root':e,'smooth':r['smooth']} for r in rational+quadratic for e,c in zip(r['roots'],r['norm_codes']) if c==13412]
    assert pairs13412==[[[35,114],[75,11]]]
    assert norm13412==[{'t':4161,'root':5100,'smooth':False}]*2
    split_sites=[r['t'] for r in rational if len(r['rational_roots'])==3]
    character={'rational_root_codes':110,'distinct_rational_root_codes':110,'quadratic_zero_rows':zero_rows,
               'rational_pairs_norm13412':pairs13412,'quadratic_norm13412':norm13412,'split_sites':split_sites,
               'discriminant_factors':[[[int(v) for v in h.list()],int(e)] for h,e in factors],
               'node_root':[int(v) for v in node.list()],'simple_root':[int(v) for v in simple.list()]}
    write('character-gate.json',character)
    g=(t-35)*(t-75);D=q*g
    v35=(F(114)-node(35))/q(35);v75=(F(11)-node(75))/q(75)
    x0=node+q*(v35+(v75-v35)/F(40)*(t-35));assert x0.degree()<4
    etale=[]
    for c in F:
        X=x0+c*D;S,rem=(X**3+A*X+B).quo_rem(D);assert not rem and S
        etale.append({'c':int(c),'quotient':[int(v) for v in S.list()],
                      'square_up_to_scalar':bool(S.monic().is_square())})
    assert not any(r['square_up_to_scalar'] for r in etale)
    write('etale-coefficients.json',{'D':[int(v) for v in D.list()],'x0':[int(v) for v in x0.list()],'rows':etale})
    points=read(ROOT/NORM4/'norm4-sections.json')['records'];incidence=defaultdict(list)
    for point in points:
        X,Y=R(point['x']),R(point['y'])
        for b in list(F)+[None]:
            xx,yy=(X[4],Y[6]) if b is None else (X(b),Y(b))
            if not yy:incidence[None if b is None else int(b)].append((point['index'],int(xx)))
    by_pair=defaultdict(dict)
    for b,rows in incidence.items():
        for k,(i,e) in enumerate(rows):
            for j,f in rows[k+1:]:by_pair[(i,j)][b]=(e==f)
    mixed=[{'pair':list(ij),'agree':[b for b,v in sites.items() if v], 'disagree':[b for b,v in sites.items() if not v]}
           for ij,sites in by_pair.items() if any(sites.values()) and not all(sites.values())]
    assert len(mixed)==21
    Z=Integers(p*p);R2=PolynomialRing(Z,'t')
    def lift(r):return R2([int(v) for v in r.list()])
    def poly2(key):return R2([Z(QQ(v).numerator())/Z(QQ(v).denominator()) for v in source['weierstrass_model'][key+'_coefficients_low_to_high']])
    A2,B2=poly2('A'),poly2('B');lift_rows=[]
    for pair in mixed:
        for b in pair['disagree']:
            for c in pair['agree']:
                assert b in split_sites and b!=c
                H=(t-b)*(R(1) if c is None else t-c);D=H*H;fixed=2 if c is None else 4
                J=matrix(F,26,24);errors=[]
                for n,index in enumerate(pair['pair']):
                    X,Y=R(points[index]['x']),R(points[index]['y']);S,rem=Y.quo_rem(H)
                    assert not rem and S.degree()<=4
                    local=[t**k*(3*X*X+A) for k in range(5)]+[-2*D*S*t**k for k in range(5)]
                    dcols=[-S*S*t**k for k in range(5) if k!=fixed]
                    for j,f in enumerate(local):
                        for k in range(13):J[13*n+k,10*n+j]=f[k]
                    for j,f in enumerate(dcols):
                        for k in range(13):J[13*n+k,20+j]=f[k]
                    error=lift(X)**3+A2*lift(X)+B2-lift(D)*lift(S)**2
                    assert all(int(error[k])%p==0 for k in range(13))
                    errors += [int(F(-int(error[k])//p)) for k in range(13)]
                rank=int(J.rank());aug=int(J.augment(vector(F,errors).column()).rank());assert (rank,aug)==(24,25)
                lift_rows.append({'pair':pair['pair'],'disagree':b,'agree':c,'fixed_D_coefficient':fixed,
                                  'error_target':errors,'jacobian_rank':rank,'augmented_rank':aug})
    assert len(lift_rows)==21;write('split-lifts.json',{'rows':lift_rows})
    cpp_input='\n'.join(' '.join(str(int(v)) for v in row) for row in [A.list(),B.list(),simple.list()])+'\n'
    with tempfile.TemporaryDirectory(prefix='q80-nodal-census-') as temp:
        executable=Path(temp)/'census'
        subprocess.run(['g++','-O3','-std=c++17',str(ROOT/CPP),'-o',str(executable)],check=True,timeout=30)
        completed=subprocess.run([str(executable)],input=cpp_input,text=True,capture_output=True,check=True,timeout=30)
    census=json.loads(completed.stdout);assert census['status']=='PASS' and census['tried']==p**3 and census['accepted']==2
    for row in census['records']:
        X,S=R(row['x']),R(row['s']);assert X**3+A*X+B==F(row['scalar'])*q*S*S
        assert X%q==simple
        row['rational_s_roots']=[int(b) for b in F if not S(b)]
        assert not any(b in split_sites for b in row['rational_s_roots'])
    write('nodal-census.json',census)
    result={'status':'PASS','input_sha256':digest(OUT/'input.json'),
            'records':{s:digest(OUT/s) for s in ['character-gate.json','etale-coefficients.json','split-lifts.json','nodal-census.json']},
            'integral_k2_solutions':'EXCLUDED over Q_131 by the written classification and finite gates',
            'necessary_D_reduction':[int(v) for v in (q*q).list()],
            'both_abscissas_require_negative_gauss_valuation':True,
            'remaining_nodal_denominator_boundary':'UNKNOWN',
            'positive_mw17_target_complete':False,'old_norm8_replay_upgraded':False}
    write('result.json',result)
    write('execution.json',{'status':'PASS','input_sha256':digest(OUT/'input.json'),
          'result_sha256':digest(OUT/'result.json'),'sage_cpu_seconds':time.process_time()-start,
          'census_cpu_seconds':census['cpu_seconds'],'wall_seconds':time.monotonic()-wall,
          'producer_cpu_limit':40,'census_cpu_limit':30,'memory_limit_bytes':4*1024**3})
    print(json.dumps(result,sort_keys=True),flush=True)
if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('action',choices=['freeze','run']);args=parser.parse_args()
    freeze() if args.action=='freeze' else run()
