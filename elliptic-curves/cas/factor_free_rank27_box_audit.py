#!/usr/bin/env python3
"""Freeze/replay factor-free maps on all seven retained catalogue-unmatched27 seeds."""
import argparse,sys
from pathlib import Path
import certify_compact_r17_candidates as cert
import pari_pointed_backend as backend
from memory_rank_certificate import checked_rank
from research_runtime.store import checkpoint,digest
from research_runtime.supervisor import run,Limits
from research_runtime.projective_box_change import classify
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';ART=ROOT/'artifacts/generated-results/elliptic-curves';LOCAL=ROOT/'artifacts/local/elliptic-curves'
D=LOCAL/'factor-free-rank27-box-audit-v1';OUT=ART/'factor_free_rank27_box_audit_v1.json'

def sources():
    names=['factor_free_rank27_box_audit.py','prepare_factor_free_rank27_boxes.sage','factor_free_pari_mapping.sage','research_runtime/projective_box_change.py','memory_rank_certificate.py']
    return {**backend.sources(),**{str((CAS/n).relative_to(ROOT)):cert.hashed(CAS/n) for n in names}}

def protocol():
    p=cert.read(D/'protocol.json')
    assert p['sources']==sources() and all(cert.hashed(ROOT/n)==h for n,h in p['inputs'].items())
    return p

def freeze():
    if (D/'protocol.json').exists():raise FileExistsError('preserve map audit')
    index=ART/'new_high_rank_curve_index_v22.json';inventory=cert.read(index)
    ids={'new-20260906-'+n for n in ('40','41','48','71','72','90','186')}
    candidates={r['id']:r for r in inventory['curves'] if r['rank_lower_bound']==27}
    assert set(candidates)==ids
    paths=[index,ART/'factor_free_known28_control_v1.json'];rows=[]
    for ident in sorted(ids):
        group='full11952-specialized-followup-v1' if ident.endswith('-186') else 'new27-specialized-parity-six-v1'
        folder=LOCAL/group/ident;seed=cert.read(folder/'seed.json');maps=cert.read(folder/'maps.json');data=cert.read(folder/'result.json')
        assert len(seed['points'])==27 and seed['curve']==candidates[ident]['curve']
        assert maps['status']=='COMPLETE_DECLARED_MAPS' and len(maps['rows'])==len(data['charts'])==49
        assert data['maps_sha256']==cert.hashed(folder/'maps.json') and data['rank_lower_bound']==27
        assert all(c['search']['status']=='bounded_search_complete' and c['search']['height_bound']==125000 for c in data['charts'])
        proof=seed['rank_certificate'];actual=checked_rank(tuple(map(cert.F,seed['curve'])),[tuple(map(cert.F,P)) for P in seed['points']],[s['prime'] for s in proof['signatures']],proof['no_rational_2_torsion_prime']);assert digest(actual)==digest(proof)
        paths.extend([folder/'seed.json',folder/'maps.json',folder/'result.json'])
        rows.append({'id':ident,'family':seed['family'],'parameter':seed['parameter'],'folder':str(folder.relative_to(ROOT))})
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.factor-free-rank27-box-audit.v1','sources':sources(),'inputs':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},'rows':rows,'height':125000,'wall_seconds_per_curve':120,'rss_bytes':1610612736,'maximum_workers':1,'scope':'Exactly the existing49 own27 centres on all seven V22 rank27 equations, with no new parameter, mask, height fit or point search. Construct exact factor-free maps and compare finite coordinate boxes. A coordinate witness need not be a quartic square or new elliptic point. No automatic point campaign.'})

def launch():
    p=protocol();ledger=D/'ledger.json'
    if ledger.exists():raise FileExistsError('preserve geometry execution')
    data={'status':'RUNNING','rows':[]};checkpoint(ledger,data)
    for i,row in enumerate(p['rows']):
        folder=D/row['id'];s=run(['/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python',str(CAS/'prepare_factor_free_rank27_boxes.sage'),str(i)],limits=Limits(p['wall_seconds_per_curve'],p['rss_bytes']),log_path=folder/'maps.log',checkpoint_path=folder/'maps.supervisor.json',cwd=ROOT)
        ok=s['outcome']=='completed' and s['returncode']==0;data['rows'].append({'id':row['id'],'supervision':s});checkpoint(ledger,data)
        print(row['id'],s['outcome'],s['returncode'],flush=True)
        if not ok:data['status']='FAILED_OR_CENSORED';checkpoint(ledger,data);raise ArithmeticError('geometry failed')
    data['status']='PASS';checkpoint(ledger,data)

def expected():
    p=protocol();paths=[D/'protocol.json',D/'ledger.json'];ledger=cert.read(paths[-1]);assert ledger['status']=='PASS';rows=[]
    for metadata in p['rows']:
        path=D/metadata['id']/'maps.json';paths.append(path);data=cert.read(path);old=cert.read(ROOT/metadata['folder']/'maps.json')
        assert data['status']=='PASS' and data['protocol_hash']==digest(p) and len(data['rows'])==49
        findings=[]
        for original,new in zip(old['rows'],data['rows']):
            assert original['centre']==new['mapping']['centre']
            actual=classify(original['matrix'],new['mapping']['matrix'],p['height']);assert actual==new['box_comparison'];findings.append(actual['status'])
        rows.append({**metadata,'counts':{name:findings.count(name) for name in sorted(set(findings))}})
    return {'schema':'elliptic-curves.factor-free-rank27-box-summary.v1','status':'PASS','sources':{str(x.relative_to(ROOT)):cert.hashed(x) for x in paths},'rows':rows,'charts':343,'point_boxes_attempted':0,'boundary':p['scope']}

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('mode',choices=['freeze','launch','report','check']);mode=a.parse_args().mode
    if mode in ('freeze','launch'):globals()[mode]()
    else:
        r=expected()
        if mode=='check':assert r==cert.read(OUT)
        else:
            if OUT.exists():raise FileExistsError('preserve result')
            checkpoint(OUT,r)
        print('PASS',[(x['id'],x['counts']) for x in r['rows']])
