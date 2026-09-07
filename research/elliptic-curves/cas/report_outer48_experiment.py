#!/usr/bin/env python3
"""Bind deeper compact point exposure to its selection and exact certificates."""
import argparse,json
from collections import Counter
from pathlib import Path
import certify_compact_r17_candidates as cert
import outer48_r17_pari_batch as batch
import extend_outer131072_r17 as selection
from research_runtime.store import checkpoint
ROOT=batch.ROOT;ART=batch.ART;OUT=ART/'outer48_experiment_v1.json'

def expected():
    p=batch.protocol();paths={Path(__file__).resolve(),batch.D/'protocol.json',selection.D/'protocol.json'}
    def read(path):paths.add(path);return cert.read(path)
    def terminal(path):
        s=read(path)
        if s['outcome']!='completed' or s['returncode']!=0:raise ArithmeticError('required terminal stage failed or censored')
    selected=read(selection.D/'result.json')
    if selected['status']!='COMPLETE_FROZEN_TRACE_EXTENSION' or selected['selection']!=selection.selection(selected['rows']):raise ArithmeticError('fixed outer48 selection differs')
    terminal(selection.D/'replay.supervisor.json')
    maps=read(batch.D/'maps-ledger.json');point=read(batch.D/'ledger.json');v=read(batch.D/'verification-ledger.json');odd=read(batch.D/'odd-cloud-audit/ledger.json');proof=read(ART/'outer48_r17_results_v1.json')
    for label in ('certify','proof-replay','geometry'):terminal(batch.D/'post-batch'/(label+'.supervisor.json'))
    roster=[r['id'] for r in p['rows']]
    if len(roster)!=48 or len(set(roster))!=48 or maps['status']!='PASS' or point['status']!='COMPLETE_FIXED_BATCH_ATTEMPTS' or v['status']!='PASS' or odd['status']!='PASS':raise ArithmeticError('all48 required stages must be terminal')
    for d in (maps,point,v,odd):
        if [r['id'] for r in d['rows']]!=roster:raise ArithmeticError('fixed48 order differs')
    if [r['id'] for r in proof['curves']]!=roster or len(proof['previous_equations'])!=741 or proof['catalogue']['curve_count']!=593:raise ArithmeticError('point proof roster or pinned comparison population differs')
    rows=[]
    for r,m,w,c,o,q in zip(p['rows'],maps['rows'],point['rows'],v['rows'],odd['rows'],proof['curves']):
        folder=batch.D/r['id'];raw=read(folder/'result.json');cloud=read(ROOT/c['cloud_certificate']);ell=read(ROOT/o['output']);declared=len(p['generic_masks'][r['family']]);attempted=len(raw['charts']);completed=sum(a['search']['status']=='bounded_search_complete' for a in raw['charts'])
        if cert.hashed(folder/'maps.json')!=m['maps_sha256'] or cert.hashed(folder/'result.json')!=w['result_sha256'] or w['result_sha256']!=c['input_sha256'] or cert.hashed(ROOT/c['cloud_certificate'])!=c['cloud_sha256'] or cert.hashed(ROOT/o['output'])!=o['output_sha256']:raise ArithmeticError('immutable map/point/cloud binding differs')
        if q['rank_lower_bound']!=cloud['rank_lower_bound'] or o['mod2_lower_bound']!=cloud['rank_lower_bound'] or {str(a['modulus']):a['finite_column_rank'] for a in ell['audits']}!=o['odd_lower_bounds']:raise ArithmeticError('point lower bounds differ')
        if any(a['search']['height_bound']!=125000 or a['search']['timeout_seconds']!=10 for a in raw['charts']):raise ArithmeticError('fixed point budget differs')
        if raw['status']=='COMPLETE_DECLARED_POINT_ATTEMPT':
            if attempted!=declared:raise ArithmeticError('short complete point attempt')
        elif raw['status']=='TARGET_REACHED_PENDING_REPLAY':
            if raw['rank_lower_bound']<28:raise ArithmeticError('unsupported target stop')
        else:raise ArithmeticError('nonterminal point status')
        rows.append({'id':r['id'],'family':r['family'],'parameter':r['parameter'],'rank_lower_bound':q['rank_lower_bound'],'odd_modulus_lower_bounds':o['odd_lower_bounds'],'declared_charts':declared,'attempted_charts':attempted,'completed_boxes':completed,'retained_points':len(cloud['points']),'catalogue_matches':q['icarm_matches'],'previous_matches':q['previous_matches']})
    if sum(r['declared_charts'] for r in rows)!=2160:raise ArithmeticError('fixed2160 chart cap differs')
    return {'schema':'elliptic-curves.outer48-experiment.v1','status':'PASS','sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in sorted(paths)},'rows':rows,'saved_short_score_rows':6144,'outer_extended_score_rows':len(selected['rows']),'declared_point_boxes':2160,'attempted_point_boxes':sum(r['attempted_charts'] for r in rows),'completed_point_boxes':sum(r['completed_boxes'] for r in rows),'rank_lower_bound_counts':dict(sorted(Counter(r['rank_lower_bound'] for r in rows).items())),'catalogue_snapshot_curves':593,'previous_address_equations':741,'within_batch_isomorphic_pairs':proof['within_batch_isomorphic_pairs'],'catalogue_unmatched_candidates_at_least22':sum(r['rank_lower_bound']>=22 and not r['catalogue_matches'] and not r['previous_matches'] for r in rows),'claim_boundary':'A fixed point-exposure cohort of48 equations in the outer32768through131072 parameter envelope, selected by unchanged S1 with validation excluded from all ties. Previous and within-roster equation exclusions do not assert random sampling. Exact point lower bounds and pinned rational-isomorphism comparisons require separate inventory promotion. These results do not establish whole-curve rank, upper bounds, saturation, score optimality, rank density, point absence or universal novelty.'}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();d=expected()
    if a.check:
        if cert.read(OUT)!=json.loads(json.dumps(d)):raise ArithmeticError('outer48 aggregate differs')
    else:
        if OUT.exists():raise FileExistsError('preserve outer48 experiment aggregate')
        checkpoint(OUT,d)
    print('OUTER48 EXACT EXPERIMENT',d['rank_lower_bound_counts'],d['completed_point_boxes'],'completed boxes;',d['catalogue_unmatched_candidates_at_least22'],'unmatched >=22 candidates',flush=True)
