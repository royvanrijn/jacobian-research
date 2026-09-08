#!/usr/bin/env python3
"""Post-terminal exact point certificates and Q-isomorphism comparisons."""
import argparse,json
from pathlib import Path
import certify_compact_r17_candidates as cert
import compact_atlas_specialization as spec
from memory_rank_certificate import checked_rank
import next24_r17_pari_batch as batch
ROOT=batch.ROOT;ART=batch.ART;D=batch.D;OUT=ART/'next24_r17_results_v1.json';DATABASE=ROOT/'artifacts/local/elliptic-curves/next12-current-catalogue-v1/database.json';OLD=ART/'fresh_r17_paired_results_v1.json'

def verify(row, families, projection, previous):
    model = tuple(map(cert.F, row['curve']))
    points = [tuple(map(cert.F, p)) for p in row['points']]
    proof = row['rank_certificate']
    actual = checked_rank(model, points, [s['prime'] for s in proof['signatures']], proof['no_rational_2_torsion_prime'])
    if json.dumps(actual, sort_keys=True) != json.dumps(proof, sort_keys=True) or row['rank_lower_bound'] != len(points):
        raise ArithmeticError('independent point proof differs')
    original, generic = spec.specialize(families[row['family']], row['parameter'])
    u = cert.F(row['family_to_curve_scale_u'])
    if not u or model != (cert.F(0), cert.F(0), cert.F(0), original[3]/u**4, original[4]/u**6):
        raise ArithmeticError('family model transport differs')
    if tuple((x/u**2, y/u**3) for x, y in generic) != tuple(tuple(map(cert.F, p)) for p in row['generic_points']):
        raise ArithmeticError('generic point transport differs')
    matches = [r['id'] for r in projection if cert.isomorphic(model, r['ainvs'])]
    old = [r['address'] for r in previous if cert.isomorphic(model, r['curve'])]
    if matches != row['icarm_matches'] or old != row['previous_matches']:
        raise ArithmeticError('post-batch Q-isomorphism comparison differs')

def sources():
    return {str(p.relative_to(ROOT)):cert.hashed(p) for p in (Path(__file__).resolve(),Path(cert.__file__).resolve(),Path(spec.__file__).resolve(),spec.ATLAS,ROOT/'elliptic-curves/cas/memory_rank_certificate.py',ROOT/'elliptic-curves/cas/mod2_reduction_independence.py',ROOT/'elliptic-curves/cas/elliptic_candidate_record.py')}

def build():
    if OUT.exists():raise FileExistsError('preserve next24 certificate')
    p=batch.protocol();ledger=cert.read(D/'ledger.json');verified=cert.read(D/'verification-ledger.json')
    if ledger['status']!='COMPLETE_FIXED_BATCH_ATTEMPTS' or verified['status']!='PASS' or [r['id'] for r in verified['rows']]!=[r['id'] for r in p['rows']]:raise ArithmeticError('terminal attempts/replays required before catalogue access')
    projection=[{'id':r['id'],'ainvs':r['ainvs']} for r in cert.read(DATABASE)['curves']];old=cert.read(OLD);previous=old['previous_equations']+[{'address':OLD.name+':'+r['family']+':'+r['parameter'],'curve':r['curve']} for r in old['curves']];families={f['family']:f for f in cert.read(spec.ATLAS)['families']};rows=[]
    for entry in ledger['rows']:
        folder=D/entry['id'];path=folder/'result.json';data=cert.read(path);v=cert.read(folder/'verification.json');cloud_path=ROOT/v['cloud_certificate'];cloud=cert.read(cloud_path)
        if v['status']!='PASS' or v['input_sha256']!=cert.hashed(path) or entry['result_sha256']!=cert.hashed(path) or cloud['input_sha256']!=cert.hashed(path) or cloud['curve']!=data['curve']:raise ArithmeticError('terminal cloud binding differs')
        model=tuple(map(cert.F,data['curve']));row={'id':entry['id'],'arms':entry['arms'],'retained_index':entry['retained_index'],'family':data['family'],'parameter':data['parameter'],'curve':data['curve'],'points':cloud['independent_points'],'generic_points':data['generic_points'],'family_to_curve_scale_u':data['family_to_curve_scale_u'],'rank_certificate':cloud['rank_certificate'],'rank_lower_bound':cloud['rank_lower_bound'],'icarm_matches':[r['id'] for r in projection if cert.isomorphic(model,r['ainvs'])],'previous_matches':[r['address'] for r in previous if cert.isomorphic(model,r['curve'])],'attempted_charts':len(data['charts']),'declared_charts':len(data['centres']),'completed_boxes':sum(c['search']['status']=='bounded_search_complete' for c in data['charts']),'search_status':data['status'],'supervision_status':entry['status'],'discovery_witness':{'path':str(path.relative_to(ROOT)),'sha256':cert.hashed(path)},'complete_cloud_certificate':{'path':str(cloud_path.relative_to(ROOT)),'sha256':cert.hashed(cloud_path)}};verify(row,families,projection,previous);rows.append(row);print('CERTIFIED NEXT24 R17',row['id'],row['parameter'],'rank >=',row['rank_lower_bound'],'catalogue',row['icarm_matches'],'previous',len(row['previous_matches']),flush=True)
    pairs=[[j,i] for i,r in enumerate(rows) for j,s in enumerate(rows[:i]) if cert.isomorphic(r['curve'],s['curve'])]
    cert.write(OUT,{'schema':'elliptic-curves.next24-r17-results.v1','sources':sources(),'curves':rows,'within_batch_isomorphic_pairs':pairs,'point_protocol_sha256':cert.hashed(D/'protocol.json'),'point_ledger_sha256':cert.hashed(D/'ledger.json'),'verification_ledger_sha256':cert.hashed(D/'verification-ledger.json'),'previous_equations':previous,'previous_source':{'path':str(OLD.relative_to(ROOT)),'sha256':cert.hashed(OLD)},'catalogue':{'url':'https://elliptic-rank.icarm.cloud/database.json','raw_sha256':cert.hashed(DATABASE),'curve_count':len(projection),'equations':projection,'acknowledgement':'ICARM, supported by NSF Grant DMS2425401'},'claim_boundary':'Exact point lower bounds for the frozen next24 candidate cohort. Q-isomorphism and catalogue comparison occurs only after all attempts/replays. No exact rank, global minimality, exact conductor or universal novelty is inferred. Validation-prime enrichment and discovery yield are finite experimental outcomes.'})

def check():
    d=cert.read(OUT)
    if d['sources']!=sources():raise ArithmeticError('next24 certificate source differs')
    families={f['family']:f for f in cert.read(spec.ATLAS)['families']}
    for r in d['curves']:verify(r,families,d['catalogue']['equations'],d['previous_equations'])
    pairs=[[j,i] for i,r in enumerate(d['curves']) for j,s in enumerate(d['curves'][:i]) if cert.isomorphic(r['curve'],s['curve'])]
    if pairs!=d['within_batch_isomorphic_pairs']:raise ArithmeticError('within-batch Q-isomorphism comparison differs')
    print('REPLAYED NEXT24 R17 CERTIFICATES',len(d['curves']),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();check() if a.check else build()
