#!/usr/bin/env sage
"""Certify the complete norm4/repeated-contact input to Q80 genus0 closure."""
from sage.all import GF, QQ, ZZ, PolynomialRing, matrix, vector
from collections import Counter, defaultdict
from functools import lru_cache
from hashlib import sha256
from pathlib import Path
import argparse
import json
import resource
import time

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'artifacts/generated-results/elkies-k3-q80-genus-zero-integral-reduction-v1'
SOURCE='artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json'
CENSUS='artifacts/generated-results/elkies-k3-q80-alternate-rootless-bisection-orbits.json'
ALT='artifacts/generated-results/q80-alternate-fifth-q6-rootless-transport.json'
FROB='artifacts/generated-results/elliptic-curves/rank_jump_residual_quartic_at_131_v1.json'
FROB_REPLAY='artifacts/generated-results/elliptic-curves/rank_jump_residual_quartic_at_131_verification_v1.json'
CHECKER='elkies-k3/scripts/verify_q80_genus_zero_collision_closure.py'
PRIOR='artifacts/generated-results/elkies-k3-q80-degree-two-reciprocity-v1'
HELPERS=['elkies-k3/scripts/verify_q80_degree_two_reciprocity.py',
 'elkies-k3/scripts/verify_q80_single_branch_reciprocity.py',
 'elkies-k3/scripts/verify_q80_complete_genus_one_pairs.py',
 'elkies-k3/scripts/verify_r17_mestre_shared_twist.py',
 'elkies-k3/scripts/verify_r17_correlated_genus_one_pencils.py']
def read(path):return json.loads(path.read_text())
def digest(path):return sha256(path.read_bytes()).hexdigest()
def write(name,data):
    with (OUT/name).open('x') as f:json.dump(data,f,indent=2,sort_keys=True);f.write('\n')
def freeze():
    files=[SOURCE,CENSUS,ALT,FROB,FROB_REPLAY,CHECKER,str(Path(__file__).relative_to(ROOT))]+HELPERS
    files += [PRIOR+'/'+name for name in ['input.json','result.json','rational-fibres.json','independent-replay.json']]
    preflight=read(OUT/'preflight-norm4.json')
    packet={'schema':1,'prime':131,'unoriented_norm4_count':1313,
        'words':[r['word'] for r in preflight['rows']],
        'collision_sites':read(ROOT/PRIOR/'result.json')['rational_full_splitting_values'],
        'bindings':{s:digest(ROOT/s) for s in files},'cpu_seconds':40,'memory_bytes':4*1024**3,
        'preflight_sha256':digest(OUT/'preflight-norm4.json'),
        'scope':'Complete finite norm4 and repeated-contact proof input; written integrality, Hensel and constant-twist arguments close genus0 k0. No positive MW17 endpoint.'}
    write('input.json',packet);print(json.dumps({'status':'FROZEN','input_sha256':digest(OUT/'input.json')}),flush=True)
def run():
    packet=read(OUT/'input.json')
    for s,h in packet['bindings'].items():assert digest(ROOT/s)==h
    resource.setrlimit(resource.RLIMIT_CPU,(40,45));resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3))
    start=time.process_time();source=read(ROOT/SOURCE)
    F=GF(131);R=PolynomialRing(F,'t');t=R.gen();K=R.fraction_field()
    def poly(cs):return R([F(QQ(v)) for v in cs])
    A=poly(source['weierstrass_model']['A_coefficients_low_to_high']);B=poly(source['weierstrass_model']['B_coefficients_low_to_high'])
    G=matrix(ZZ,source['sections']['height_gram']);assert G.det()==948
    C=matrix(ZZ,source['sections']['coordinate_matrix_in_compiled_frame'])
    H=matrix(ZZ,source['frame_certificate']['frame_gram'])
    U=matrix(ZZ,source['frame_certificate']['integral_isometry_to_alternate_Q80'])
    HA=matrix(ZZ,read(ROOT/ALT)['rootless_frame'])
    assert abs(C.det())==abs(U.det())==1 and C*H*C.transpose()==G and U*HA*U.transpose()==H
    census=read(ROOT/CENSUS);assert census['streaming_enumeration']['signed_shell_counts']['4']==2626
    assert census['input']['rootless_frame_sha256']==read(ROOT/ALT)['rootless_frame_sha256']
    basis=[]
    for row in source['sections']['records']:
        pair=[]
        for key in ['X','Y']:
            v=row[key];pair.append(K(poly(v['numerator_coefficients_low_to_high']))/poly(v['denominator_coefficients_low_to_high']))
        assert pair[1]**2==pair[0]**3+A*pair[0]+B
        basis.append(tuple(pair))
    def plus(P,Q):
        if P is None:return Q
        if Q is None:return P
        x,y=P;u,v=Q
        if x==u:
            if y==-v:return None
            s=(3*x*x+A)/(2*y)
        else:s=(v-y)/(u-x)
        z=s*s-x-u;return z,s*(x-z)-y
    def height(P):
        if P is None:return 0
        x=P[0];return max(4+x.denominator().degree(),x.numerator().degree())
    reduced_gram=[]
    for i,P in enumerate(basis):
        row=[]
        for j,Q in enumerate(basis):
            value=height(P) if i==j else (height(plus(P,Q))-height(P)-height(Q))//2
            assert value==G[i,j];row.append(int(value))
        reduced_gram.append(row)
    codes=[c for r in read(ROOT/PRIOR/'rational-fibres.json')['rows'] for c in r['rational_codes']]
    assert matrix(GF(2),[[c>>i&1 for i in range(17)] for c in codes]).rank()==17
    old,replay=read(ROOT/FROB),read(ROOT/FROB_REPLAY)
    for data in [old,replay]:
        assert data['status']=='PASS'
        for s,h in data['bindings'].items():assert digest(ROOT/s)==h
    assert old['Frobenius_traces']==[1884,319520]
    T1,T2=1884-17*131,319520-17*131**2
    assert (T1,T2)==(-343,27783) and T1+131==-212
    assert (T1+131)**2==T2+131**2
    @lru_cache(None)
    def point(word):
        if not any(word):return None
        v=vector(ZZ,word);gv=G*v
        i=max((2*(1 if c>0 else -1)*gv[i]-G[i,i],i) for i,c in enumerate(word) if c)[1]
        sign=1 if word[i]>0 else -1;rest=list(word);rest[i]-=sign
        x,y=basis[i];return plus(point(tuple(rest)),(x,sign*y))
    words=packet['words'];assert len(words)==len({tuple(w) for w in words})==1313
    rows=[];repeated=[];incidences=0;determinants=[]
    for index,word in enumerate(words):
        w=vector(ZZ,word);assert w*G*w==4
        P=point(tuple(word));assert P is not None
        x,y=P;assert x.denominator()==y.denominator()==1
        x,y=R(x),R(y);a=3*x*x+A
        assert x.degree()<=4 and y.degree()<=6 and y*y==x*x*x+A*x+B
        assert a.gcd(y)==1 and (y[6]!=0 or a[8]!=0)
        roots=[]
        for b in packet['collision_sites']:
            if y(b):continue
            order=1;d=y.derivative()
            while not d(b):order+=1;d=d.derivative();assert d
            u=t-b;s=R(y/u)
            columns=[a*t**j for j in range(5)]+[-2*u*u*s*t**j for j in range(6)]+[-s*s*t**j for j in range(2)]
            det=int(matrix(F,13,13,lambda i,j:columns[j][i]).det())
            assert bool(det)==(order==1)
            rr={'b':b,'root':int(x(b)),'order':order,'jacobian_determinant':det}
            roots.append(rr);incidences+=1
            if order>1:repeated.append({'index':index,**rr})
            else:determinants.append(det)
        rows.append({'index':index,'word':word,'x':[int(v) for v in x.list()],'y':[int(v) for v in y.list()],'collision_roots':roots})
    assert len({tuple(r['x']) for r in rows})==1313
    fibres=defaultdict(set)
    for r in repeated:fibres[r['b']].add(r['root'])
    assert all(len(v)==1 for v in fibres.values())
    assert incidences==443 and len(repeated)==4 and len(determinants)==439
    write('norm4-sections.json',{'prime':131,'records':rows})
    write('finite-lattice.json',{'height_gram':reduced_gram,'determinant':948,'kummer_character_rank':17,
        'inherited_norm4_signed_count':2626,'residual_H2_eigenvalue_polynomial_coefficients':[2248091,44933,343,1],
        'constant_nonsquare_twist_rank_upper_bound':1})
    result={'status':'PASS','input_sha256':digest(OUT/'input.json'),
        'records':{s:digest(OUT/s) for s in ['norm4-sections.json','finite-lattice.json']},
        'unoriented_norm4_sections':1313,'collision_branch_incidences':443,
        'simple_root_hensel_exclusions':439,'repeated_root_records':repeated,
        'repeated_roots_by_fibre':{str(b):sorted(v) for b,v in sorted(fibres.items())},
        'repeated_root_pairs_selecting_distinct_torsion':0,
        'genus0_k0_all_five_allocations':'EXCLUDED over Q_131 by the written proof',
        'genus1_allocations':'UNKNOWN','positive_mw17_target_complete':False,
        'old_norm8_singular_replay_upgraded':False}
    write('result.json',result)
    write('execution.json',{'status':'PASS','input_sha256':digest(OUT/'input.json'),'result_sha256':digest(OUT/'result.json'),
        'cpu_seconds':time.process_time()-start,'cpu_limit_seconds':40,'memory_limit_bytes':4*1024**3})
    print(json.dumps(result,sort_keys=True),flush=True)
if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('action',choices=['freeze','run']);args=parser.parse_args()
    freeze() if args.action=='freeze' else run()
