#!/usr/bin/env python3
"""One bounded exact incidence audit of the84 inventory additions."""
import argparse
from pathlib import Path
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';ART=ROOT/'artifacts/generated-results/elliptic-curves'
D=ROOT/'artifacts/local/elliptic-curves/inventory185-added-incidence-v1'
OUT=ART/'inventory185_added_cross_family_j_incidence_v1.json'
REPLAY=ART/'inventory185_added_cross_family_j_incidence_replay_v1.json'
SAGE='/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python'

def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve fixed84 incidence protocol')
    old=cert.read(ART/'inventory101_incidence_v1.json');index=cert.read(ART/'new_high_rank_curve_index_v13.json')
    rows=[r['id'] for r in index['curves'] if int(r['id'].split('-')[-1])>101]
    if old['status']!='PASS' or old['targets_checked']!=101 or len(index['curves'])!=185 or len(rows)!=84:raise ArithmeticError('exact prior101 and fixed84 additions required')
    paths=[Path(__file__).resolve(),CAS/'audit_inventory185_added_incidence.sage',CAS/'replay_inventory185_added_incidence.py',ART/'inventory101_incidence_v1.json',ART/'new_high_rank_curve_index_v13.json',ART/'compact_six_r17_atlas_v1.json',ART/'compact_five_mw16_atlas_v1.json',ROOT/'elkies-k3/data/fibrations/elkies_2026_published_r17_model.json']
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.inventory185-added-incidence.v1','sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},'rows':rows,'presentations':12,'pairs':1008,'audit_wall_seconds':1200,'replay_wall_seconds':600,'rss_bytes':2147483648,'maximum_workers':1,'modular_image_prime_bound':251,'residual_root_prime_bound':997,'gate':'The inventory contains84 new exact distinct curves beyond the completed101-curve incidence audit. Their original construction family is known, but alternate rational presentations are not yet classified. Compute the complete rational j-preimage problem in the same twelve recorded presentations to test this missing geometry; any genuinely additional presentation needs a separate rational-isomorphism and subgroup transport proof before it can motivate points.','scope':'Exactly84 targets and1008 family pairs, with finite projective-image exclusions first and exact factorization only for surviving pairs. Preserve every coefficient, modular exclusion, rational root and unresolved residual factor. Any factorization timeout or missing residual-root certificate remains UNKNOWN; no retry or larger prime/parameter search. Equal j alone proves neither rational isomorphism nor additional generic sections, higher rank or a new curve.'})

def launch():
    p=cert.read(D/'protocol.json')
    if any(cert.hashed(ROOT/n)!=h for n,h in p['sources'].items()):raise ArithmeticError('frozen incidence inputs changed')
    if (D/'ledger.json').exists():raise FileExistsError('preserve incidence attempt')
    ledger={'status':'RUNNING','rows':[]};checkpoint(D/'ledger.json',ledger)
    for name,executable,script,args,seconds in [('audit',SAGE,'audit_inventory185_added_incidence.sage',['--output',str(OUT)],p['audit_wall_seconds']),('replay','/usr/bin/python3','replay_inventory185_added_incidence.py',['--input',str(OUT),'--output',str(REPLAY)],p['replay_wall_seconds'])]:
        s=run([executable,str(CAS/script),*args],limits=Limits(seconds,p['rss_bytes']),log_path=D/(name+'.log'),checkpoint_path=D/(name+'.supervisor.json'),cwd=ROOT)
        ok=s['outcome']=='completed' and s['returncode']==0;ledger['rows'].append({'name':name,'status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':s});checkpoint(D/'ledger.json',ledger);print('ADDED84 INCIDENCE',name,ledger['rows'][-1]['status'],flush=True)
        if not ok:
            ledger['status']='FAILED_OR_CENSORED';checkpoint(D/'ledger.json',ledger);raise ArithmeticError('bounded incidence stopped; retain unresolved evidence')
    ledger['status']='PASS';checkpoint(D/'ledger.json',ledger)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','launch']);a=p.parse_args();globals()[a.stage]()
