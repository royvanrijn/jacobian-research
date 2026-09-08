#!/usr/bin/env python3
"""Four independent rank26/27 proofs with global minimal models and novelty checks."""
import argparse,json
from math import gcd
from pathlib import Path
import certify_compact_r17_candidates as cert
import memory_rank_certificate as memory
from mod2_reduction_independence import _is_prime
from research_runtime.store import checkpoint
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';D=ROOT/'artifacts/local/elliptic-curves/fresh-r17-paired-pari-v1';OUT=ART/'paired_high_rank_minimal_proofs_v2.json';OLD=ART/'compact_r17_wide_results_v1.json';DB=ROOT/'artifacts/local/elliptic-curves/next12-current-catalogue-v1/database.json';IDS=('074d9-007','074d9-017','07ca9-025','11952-010')

def sources():
    names=['certify_paired_high_rank_minimal_v2.py','memory_rank_certificate.py','certify_compact_r17_candidates.py','elliptic_candidate_record.py','mod2_reduction_independence.py','research_runtime/finite_reduction.py','research_runtime/memory_store.py','research_runtime/store.py']
    return {str(p.relative_to(ROOT)):cert.hashed(p) for p in (ROOT/'elliptic-curves/cas'/n for n in names)}

def valuation(a,p):
    a=abs(int(a));v=0
    if not a:raise ArithmeticError('zero invariant outside protocol')
    while a%p==0:a//=p;v+=1
    return v

def minimality(model):
    if any(x.denominator!=1 for x in model):raise ArithmeticError('integral equation required')
    inv=cert.weierstrass_invariants(model);g=gcd(abs(int(inv['c4'])),abs(int(inv['c6'])))
    if not 1<=g<=10**8:raise ArithmeticError('gcd outside cheap exact factor gate')
    remaining=g;factors=[];p=2
    while p*p<=remaining:
        if remaining%p==0:
            v=0
            while remaining%p==0:remaining//=p;v+=1
            factors.append([p,v])
        p+=1
    if remaining>1:factors.append([remaining,1])
    rows=[]
    for p,e in factors:
        if not _is_prime(p):raise ArithmeticError('gcd prime failed')
        v=[valuation(inv[k],p) for k in ('c4','c6','discriminant')]
        if not (v[0]<4 or v[1]<6 or v[2]<12):raise ArithmeticError('nonminimality not excluded at a gcd prime; needs separate Tate proof')
        rows.append({'prime':p,'valuations_c4_c6_discriminant':v})
    return {'invariant_gcd':g,'gcd_factorization':factors,'local_exclusion_rows':rows,'argument':'A nonminimal integral model at p requires v_p(c4)>=4, v_p(c6)>=6 and v_p(discriminant)>=12. Every prime dividing gcd(c4,c6) fails a necessary condition; other primes are excluded by the gcd.'}

def integral(short):
    inv=cert.weierstrass_invariants(short)
    for a1 in (0,1):
        for a2 in (-1,0,1):
            for a3 in (0,1):
                b2=cert.F(a1*a1+4*a2);b4=(b2*b2-inv['c4'])/24;b6=(-b2**3+36*b2*b4-inv['c6'])/216;model=tuple(map(cert.F,(a1,a2,a3,(b4-a1*a3)/2,(b6-a3*a3)/4)))
                if all(x.denominator==1 for x in model):return model
    raise ArithmeticError('small normalized integral model not found')

def check_row(r,catalogue,previous):
    model=tuple(map(cert.F,r['minimal_curve']));short=tuple(map(cert.F,r['discovery_curve']));inv=cert.weierstrass_invariants(model);points=[tuple(map(cert.F,P)) for P in r['points']]
    if r['minimality']!=minimality(model) or short!=(0,0,0,-inv['c4']/48,-inv['c6']/864):raise ArithmeticError('minimal model or short transport differs')
    a1,a2,a3,a4,a6=model;back=[(x+inv['b2']/12,y+(a1*x+a3)/2) for x,y in points]
    if [list(map(str,P)) for P in back]!=r['discovery_points'] or any(y*y+a1*x*y+a3*y!=x*x*x+a2*x*x+a4*x+a6 for x,y in points):raise ArithmeticError('exact point transport differs')
    p=r['rank_certificate'];actual=memory.checked_rank(short,back,[s['prime'] for s in p['signatures']],p['no_rational_2_torsion_prime'])
    if json.dumps(p,sort_keys=True)!=json.dumps(actual,sort_keys=True) or r['rank_lower_bound']!=len(points):raise ArithmeticError('finite rank proof differs')
    if r['icarm_matches']!=[q['id'] for q in catalogue if cert.isomorphic(model,q['ainvs'])] or r['previous_matches']!=[q['address'] for q in previous if cert.isomorphic(model,q['curve'])]:raise ArithmeticError('catalogue/prior comparison differs')

def build():
    if OUT.exists():raise FileExistsError('preserve minimal high-rank proofs')
    ledger=cert.read(D/'ledger.json')
    if ledger['status']!='COMPLETE_FIXED_BATCH_ATTEMPTS':raise ArithmeticError('post-batch comparison only')
    old=cert.read(OLD);previous=old['previous_equations']+[{'address':OLD.name+':'+r['family']+':'+r['parameter'],'curve':r['curve']} for r in old['curves']];projection=[{'id':r['id'],'ainvs':r['ainvs']} for r in cert.read(DB)['curves']];rows=[];inputs={}
    for identifier in IDS:
        path=D/identifier/'result.json';d=cert.read(path);entry=next(e for e in ledger['rows'] if e['id']==identifier)
        if entry['result_sha256']!=cert.hashed(path) or d['status']!='COMPLETE_DECLARED_POINT_ATTEMPT':raise ArithmeticError('terminal result binding differs')
        state=d['final_state']['state'];short=tuple(map(cert.F,d['curve']));points=[tuple(map(cert.F,P)) for P in state['reductions']['points']];proof=memory.checked_rank(short,points,state['reductions']['primes'],state['no_two_torsion_prime']);model=integral(short);inv=cert.weierstrass_invariants(model);transported=[]
        for x,y in points:
            X=x-inv['b2']/12;transported.append((X,y-(model[0]*X+model[2])/2))
        r={'id':identifier,'family':d['family'],'parameter':d['parameter'],'rank_lower_bound':len(points),'minimal_curve':list(map(str,model)),'points':[list(map(str,P)) for P in transported],'discovery_curve':d['curve'],'discovery_points':state['reductions']['points'],'rank_certificate':proof,'minimality':minimality(model),'icarm_matches':[q['id'] for q in projection if cert.isomorphic(model,q['ainvs'])],'previous_matches':[q['address'] for q in previous if cert.isomorphic(model,q['curve'])]};check_row(r,projection,previous);rows.append(r);inputs[str(path.relative_to(ROOT))]=cert.hashed(path);print('PROVED MINIMAL PAIRED CURVE',identifier,'rank >=',len(points),'gcd',r['minimality']['invariant_gcd'],flush=True)
    if any(cert.isomorphic(r['minimal_curve'],s['minimal_curve']) for i,r in enumerate(rows) for s in rows[:i]):raise ArithmeticError('four curves are not mutually distinct')
    checkpoint(OUT,{'schema':'elliptic-curves.paired-high-rank-minimal-proofs.v2','sources':sources(),'inputs':{**inputs,**{str(p.relative_to(ROOT)):cert.hashed(p) for p in (D/'ledger.json',OLD,DB)}},'curves':rows,'catalogue':{'url':'https://elliptic-rank.icarm.cloud/database.json','raw_sha256':cert.hashed(DB),'equations':projection},'previous_equations':previous,'claim_boundary':'Four mutually Q-nonisomorphic globally minimal elliptic curves, two with27 and two with26 exactly independent rational points. Three are unmatched in the pinned586-row catalogue and299 earlier equations: two with27 points and one with26. The fourth is catalogue390, an independently rediscovered known curve, not a novelty claim. These stand-alone lower bounds do not depend on a completed full search-history replay or imply exact rank, exact conductor or universal novelty.'})

def check():
    d=cert.read(OUT)
    if d['sources']!=sources() or tuple(r['id'] for r in d['curves'])!=IDS:raise ArithmeticError('fixed minimal-proof roster differs')
    for r in d['curves']:check_row(r,d['catalogue']['equations'],d['previous_equations'])
    if any(cert.isomorphic(r['minimal_curve'],s['minimal_curve']) for i,r in enumerate(d['curves']) for s in d['curves'][:i]):raise ArithmeticError('mutual distinctness differs')
    print('REPLAYED FOUR MINIMAL RANK26/27 PROOFS',flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();check() if a.check else build()
