#!/usr/bin/env python3
"""Exact integral models and point transports for unmatched higher-height R17 follow-up curves>=25."""
import argparse
from pathlib import Path
import certify_compact_r17_candidates as cert
import certify_discarded_rank26_minimal as arithmetic
from research_runtime.store import checkpoint
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';INPUT=ART/'scaled13_24_r17_results_v1.json';OUT=ART/'scaled13_24_rank25_minimal_proof_v1.json'
def sources():
    return {**arithmetic.sources(),str(Path(__file__).resolve().relative_to(ROOT)):cert.hashed(Path(__file__).resolve())}
def expected():
    d=cert.read(INPUT);selected=[r for r in d['curves'] if r['rank_lower_bound']>=25 and not r['icarm_matches'] and not r['previous_matches']];rows=[]
    if not selected:raise ArithmeticError('no unmatched>=25 cohort')
    for q in selected:
        short=tuple(map(cert.F,q['curve']));model=arithmetic.integral(short);inv=cert.weierstrass_invariants(model);points=[]
        for raw in q['points']:
            x,y=map(cert.F,raw);X=x-inv['b2']/12;points.append([str(X),str(y-(model[0]*X+model[2])/2)])
        r={'id':q['id'],'family':q['family'],'parameter':q['parameter'],'rank_lower_bound':q['rank_lower_bound'],'minimal_curve':list(map(str,model)),'points':points,'discovery_curve':q['curve'],'discovery_points':q['points'],'rank_certificate':q['rank_certificate'],'minimality':arithmetic.minimality(model),'icarm_matches':q['icarm_matches'],'previous_matches':q['previous_matches']}
        arithmetic.check_row(r,d['catalogue']['equations'],d['previous_equations']);rows.append(r)
    if any(cert.isomorphic(r['minimal_curve'],s['minimal_curve']) for i,r in enumerate(rows) for s in rows[:i]):raise ArithmeticError('selected high-rank curves not distinct')
    return {'schema':'elliptic-curves.scaled13_24-rank25-minimal.v1','status':'PASS','sources':sources(),'input_sha256':cert.hashed(INPUT),'curves':rows,'catalogue':d['catalogue'],'previous_equations':d['previous_equations'],'claim_boundary':'Exact independent-point lower bounds, proved global minimal integral models and exact point transports for the unmatched>=25 portion of the fixed24 higher-height cohort. Absence from593 pinned catalogue equations and472 earlier measured equations is not universal novelty. No exact rank or conductor is asserted.'}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();r=expected()
    if a.check:
        if cert.read(OUT)!=r:raise ArithmeticError('minimal model proof differs')
    else:
        if OUT.exists():raise FileExistsError('preserve high-rank minimal proof')
        checkpoint(OUT,r)
    print('EXACT NEW higher-height R17 HIGH-RANK MINIMAL PROOFS',[(q['id'],q['rank_lower_bound'],q['minimality']['invariant_gcd']) for q in r['curves']],flush=True)
