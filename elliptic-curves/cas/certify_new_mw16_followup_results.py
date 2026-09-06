#!/usr/bin/env python3
"""Exact post-terminal follow-up point proofs, preserving the original cohort."""
import argparse
from pathlib import Path
import certify_compact_r17_candidates as cert
import certify_extended20_mw16_results as cohort
import followup_new_mw16_rank26 as follow
from research_runtime.store import checkpoint
ROOT=follow.ROOT;ART=follow.ART;INPUT=ART/'extended20_mw16_union_results_v1.json';OUT=ART/'new_mw16_followup_results_v1.json'
def sources():
    return {**cohort.sources(),str(Path(__file__).resolve().relative_to(ROOT)):cert.hashed(Path(__file__).resolve()),str(Path(follow.__file__).resolve().relative_to(ROOT)):cert.hashed(Path(follow.__file__).resolve())}
def expected():
    old=cert.read(INPUT);ledger=cert.read(follow.LOCAL/'mw16-new26-followups-v1/ledger.json')
    if ledger['status']!='PASS' or [r['id'] for r in ledger['rows']]!=list(follow.IDS):raise ArithmeticError('three terminal exact histories required')
    rows=[];inputs={str(INPUT.relative_to(ROOT)):cert.hashed(INPUT)};families={f['fibration_id']:f for f in cert.read(cohort.spec.ATLAS)['families']}
    for identifier in follow.IDS:
        follow.configure(identifier);base=next(r for r in old['curves'] if r['id']==identifier);result=follow.D/'result.json';data=cert.read(result);cloudpath=ART/('mw16_new26_'+identifier.replace('-','_')+'_all_retained_mod2_v1.json');cloud=cert.read(cloudpath);ip=follow.D/'all-retained-point-cloud-only.json';v=follow.D/'cloud-verification-ledger.json'
        if cert.read(v)['status']!='PASS' or data['curve']!=base['curve'] or cloud['curve']!=data['curve'] or cloud['input_sha256']!=cert.hashed(ip):raise ArithmeticError('union-cloud binding differs')
        row={**base,'points':cloud['independent_points'],'rank_certificate':cloud['rank_certificate'],'rank_lower_bound':cloud['rank_lower_bound'],'initial_cohort_rank_lower_bound':base['rank_lower_bound'],'discovery_witness':{'path':str(result.relative_to(ROOT)),'sha256':cert.hashed(result)},'complete_cloud_certificate':{'path':str(cloudpath.relative_to(ROOT)),'sha256':cert.hashed(cloudpath)},'point_union_input':{'path':str(ip.relative_to(ROOT)),'sha256':cert.hashed(ip)},'attempted_adaptive_charts':len(data['charts']),'completed_adaptive_boxes':sum(c['search']['status']=='bounded_search_complete' for c in data['charts']),'adaptive_search_status':data['status']};cohort.verify(row,families,old['catalogue']['equations'],old['previous_equations']);rows.append(row)
        inputs.update({str(p.relative_to(ROOT)):cert.hashed(p) for p in [result,cloudpath,ip,v,follow.D/'protocol.json',follow.D/'maps.json']})
    return {'schema':'elliptic-curves.new-mw16-followup-results.v1','status':'PASS','sources':sources(),'inputs':inputs,'curves':rows,'catalogue':old['catalogue'],'previous_equations':old['previous_equations'],'claim_boundary':'Three already catalogue-unmatched curves, with exact point lower bounds from their initial and separately frozen301-centre follow-ups. New bounds do not create new curve IDs. Initial selector/exposure outcomes remain in the original cohort certificate. No exact rank, rank upper bound, saturation, conductor or universal novelty.'}
def main(check):
    result=expected()
    if check:
        if cert.read(OUT)!=result:raise ArithmeticError('follow-up result certificate differs')
    else:
        if OUT.exists():raise FileExistsError('preserve follow-up certificate')
        checkpoint(OUT,result)
    print('EXACT NEW MW16 FOLLOWUP RESULTS',[(r['id'],r['rank_lower_bound']) for r in result['curves']],flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();main(a.check)
