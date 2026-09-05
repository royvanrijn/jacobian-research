#!/usr/bin/env python3
"""Post-batch exact point certificates and catalogue comparison for the MW16 pilot."""
import argparse
import json
from pathlib import Path
import certify_compact_r17_candidates as cert
import compact_mw16_specialization as spec
ROOT=Path(__file__).resolve().parents[2]
PREVIOUS=('compact_r17_new_curves_v1.json','compact_r17_wide_new_curves_v1.json','compact_r17_top64_interim_curves_v1.json',
          'compact_r17_largest_gain_curve_v1.json','compact_atlas_new_curves_v1.json')

def sources():
    paths=(Path(__file__).resolve(),Path(cert.__file__).resolve(),Path(spec.__file__).resolve(),spec.ATLAS,
           ROOT/'elliptic-curves/cas/mod2_reduction_independence.py',ROOT/'elliptic-curves/cas/elliptic_candidate_record.py')
    return {str(p.relative_to(ROOT)):cert.hashed(p) for p in paths}

def verify(row,families,projection,previous):
    model=tuple(map(cert.F,row['curve']));points=tuple(tuple(map(cert.F,p)) for p in row['points']);proof=row['rank_certificate']
    actual=cert.checked_rank(model,points,[s['prime'] for s in proof['signatures']],proof['no_rational_2_torsion_prime'])
    if json.dumps(actual,sort_keys=True)!=json.dumps(proof,sort_keys=True):raise ArithmeticError('rank proof changed')
    original,generic=spec.specialize(families[row['family']],row['parameter']);u=cert.F(row['family_to_curve_scale_u'])
    if not u or model!=(cert.F(0),cert.F(0),cert.F(0),original[3]/u**4,original[4]/u**6):raise ArithmeticError('model transport changed')
    if tuple((x/u**2,y/u**3) for x,y in generic)!=tuple(tuple(map(cert.F,p)) for p in row['generic_points']):raise ArithmeticError('generic point transport changed')
    matches=[r['id'] for r in projection if cert.isomorphic(model,r['ainvs'])]
    own=[r['address'] for r in previous if cert.isomorphic(model,r['curve'])]
    if matches!=row['icarm_matches'] or own!=row['previous_matches']:raise ArithmeticError('post-batch novelty mismatch')
    if row['rank_lower_bound']!=len(points):raise ArithmeticError('rank lower bound changed')

def build(directory,output):
    if output.exists():raise FileExistsError('preserve certificate')
    ledger=cert.read(directory/'point-ledger.json')
    if ledger['status']!='COMPLETE_FIXED_BATCH_ATTEMPTS' or len(ledger['rows'])!=20 or any(r['status']=='PENDING' for r in ledger['rows']):raise ArithmeticError('catalogue boundary: fixed batch not finished')
    # This is the first catalogue access in this pilot's pipeline.
    database=cert.read(cert.DATABASE);projection=[{'id':r['id'],'ainvs':r['ainvs']} for r in database['curves']]
    previous=[];previous_hashes={}
    for name in PREVIOUS:
        p=ROOT/'artifacts/generated-results/elliptic-curves'/name;previous_hashes[str(p.relative_to(ROOT))]=cert.hashed(p)
        previous.extend({'address':name+':'+r.get('family','published-R17')+':'+r['parameter'],'curve':r['curve']} for r in cert.read(p)['curves'])
    families={f['fibration_id']:f for f in cert.read(spec.ATLAS)['families']};rows=[];missing=[]
    for entry in ledger['rows']:
        path=ROOT/entry['result_path']
        if not path.exists():missing.append(entry);continue
        if cert.hashed(path)!=entry['result_sha256']:raise ArithmeticError('terminal measurement changed')
        d=cert.read(path);model=tuple(map(cert.F,d['curve']));state=d['final_state']['state'];points=tuple(tuple(map(cert.F,p)) for p in state['reductions']['points'])
        proof=cert.checked_rank(model,points,state['reductions']['primes'],state['no_two_torsion_prime'])
        row={'family':d['family'],'parameter':d['parameter'],'curve':d['curve'],'points':[list(map(str,p)) for p in points],
            'generic_points':d['generic_points'],'family_to_curve_scale_u':d['family_to_curve_scale_u'],'rank_certificate':proof,'rank_lower_bound':len(points),
            'icarm_matches':[r['id'] for r in projection if cert.isomorphic(model,r['ainvs'])],
            'previous_matches':[r['address'] for r in previous if cert.isomorphic(model,r['curve'])],
            'completed_charts':len(d['charts']),'search_status':d['status'],'supervision_status':entry['status'],
            'discovery_witness':{'path':str(path.relative_to(ROOT)),'sha256':cert.hashed(path)}}
        verify(row,families,projection,previous);rows.append(row)
        print('CERTIFIED PROSPECTIVE MW16',row['family'],row['parameter'],'rank >=',len(points),'catalogue matches',row['icarm_matches'],flush=True)
    duplicates=[]
    for i,r in enumerate(rows):
        for j,s in enumerate(rows[:i]):
            if cert.isomorphic(r['curve'],s['curve']):duplicates.append([j,i])
    cert.write(output,{'schema':'elliptic-curves.prospective-mw16-results.v1','sources':sources(),'curves':rows,'missing_measurements':missing,
        'within_batch_isomorphic_pairs':duplicates,'point_ledger_sha256':cert.hashed(directory/'point-ledger.json'),
        'selection_protocol_sha256':cert.hashed(directory/'protocol.json'),'point_protocol_sha256':cert.hashed(directory/'point-protocol.json'),
        'catalogue':{'date':'2026-09-05','url':'https://elliptic-rank.icarm.cloud/database.json','raw_sha256':cert.hashed(cert.DATABASE),'equations':projection,
                     'acknowledgement':'ICARM, supported by NSF Grant DMS2425401'},
        'previous_equations':previous,'previous_sources':previous_hashes,
        'claim_boundary':'Independent points give exact lower bounds only. Novelty means no Q-isomorphism in the pinned584-equation catalogue and21 earlier discovered equations, not universal novelty. Censored point attempts are not upper bounds. All catalogue reads occurred after the fixed20-address batch.'})

def check(path):
    d=cert.read(path)
    if d['sources']!=sources():raise ArithmeticError('certificate source changed')
    families={f['fibration_id']:f for f in cert.read(spec.ATLAS)['families']}
    for row in d['curves']:
        verify(row,families,d['catalogue']['equations'],d['previous_equations'])
        print('REPLAYED PROSPECTIVE MW16 CERTIFICATE',row['family'],row['parameter'],'rank >=',row['rank_lower_bound'],flush=True)
    pairs=[[j,i] for i,r in enumerate(d['curves']) for j,s in enumerate(d['curves'][:i]) if cert.isomorphic(r['curve'],s['curve'])]
    if pairs!=d['within_batch_isomorphic_pairs']:raise ArithmeticError('batch deduplication changed')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);g=p.add_mutually_exclusive_group(required=True);g.add_argument('--directory',type=Path);g.add_argument('--check',type=Path);p.add_argument('--output',type=Path);a=p.parse_args()
    check(a.check) if a.check else build(a.directory.resolve(),a.output)
