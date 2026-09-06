#!/usr/bin/env python3
"""Exact retrospective population-incidence audit of saved public R17 controls."""
import argparse
from pathlib import Path
from collections import Counter
import certify_compact_r17_candidates as cert
import compact_atlas_specialization as spec
from research_runtime.store import checkpoint
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';INPUT=ROOT/'artifacts/generated-results/elkies-k3-r17-norm12-icarm-calibration-dataset-v1.json';DB=ROOT/'artifacts/local/elliptic-curves/retention24-current-catalogue-v1/database.json';OUT=ART/'public_compact_parameter_heights_v1.json';BOUNDS=(1024,4096,16384,32768,65536,131072,524288)
def expected():
    dataset=cert.read(INPUT);atlas={r['family']:r for r in cert.read(spec.ATLAS)['families']};catalogue={r['id']:r for r in cert.read(DB)['curves']};rows=[]
    if len(dataset['rows'])!=69 or len(catalogue)!=593:raise ArithmeticError('frozen retrospective cohorts differ')
    for r in dataset['rows']:
        family=r['family'].split('-')[-1];public=catalogue[r['curve_id']];entry={'id':r['curve_id'],'native_family':r['family'],'native_parameter':r['parameter'],'snapshot_reported_rank_lower_bound':r['snapshot_rank_lower_bound'],'current_reported_rank_lower_bound':public['rank_lower_bound'],'public_curve':public['ainvs']}
        if family not in atlas:
            entry['status']='NO_PINNED_COMPACT_COORDINATE_MAP';rows.append(entry);continue
        f=atlas[family];a,b,c,d=map(cert.F,f['base_matrix_a_b_c_d']);u=cert.F(r['parameter'])
        if a*d==b*c:raise ArithmeticError('singular base matrix')
        den=a-c*u
        if den==0:raise ArithmeticError('projective infinity requires separate declared model')
        t=(d*u-b)/den
        if c*t+d==0 or (a*t+b)/(c*t+d)!=u:raise ArithmeticError('exact inverse base map failed')
        n,q=t.numerator,t.denominator;model=(cert.F(0),cert.F(0),cert.F(0),spec.polynomial(f['A_coefficients_low_to_high'],t)*q**8,spec.polynomial(f['B_coefficients_low_to_high'],t)*q**12)
        if not cert.isomorphic(model,public['ainvs']):raise ArithmeticError('compact coordinate does not recover the public Q-isomorphism class')
        h=max(abs(n),q);entry.update(status='PASS_EXACT_COMPACT_Q_ISOMORPHISM',family=family,compact_parameter=str(t),height=h,compact_specialization=list(map(str,model)),inside_boxes={str(B):h<=B for B in BOUNDS});rows.append(entry)
    matched=[r for r in rows if r['status']=='PASS_EXACT_COMPACT_Q_ISOMORPHISM'];coverage={}
    for lower in (24,25,26,27,28,29):
        cohort=[r for r in matched if r['current_reported_rank_lower_bound']>=lower];coverage[str(lower)]={'matched_public_curves':len(cohort),'inside_box_counts':{str(B):sum(r['height']<=B for r in cohort) for B in BOUNDS},'minimum_height':min((r['height'] for r in cohort),default=None)}
    paths=[Path(__file__).resolve(),Path(spec.__file__).resolve(),Path(cert.__file__).resolve(),INPUT,DB,spec.ATLAS]
    return {'schema':'elliptic-curves.public-compact-parameter-heights.v1','status':'PASS','sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},'rows':rows,'status_counts':dict(Counter(r['status'] for r in rows)),'reported_rank_coverage':coverage,'claim_boundary':'Retrospective exact coordinate inversion and Q-isomorphism checks on69 saved successful public R17 fibres. Reported ranks are pinned catalogue metadata, not newly replayed point proofs. Missing frame maps remain unknown. These public successes are a selected sample without a historical failed-trial denominator; height coverage is not a success probability, density theorem or fair family-quality comparison. No public parameter or point enters a prospective candidate list, no new scan, and no claim that larger height guarantees new near-record curves.'}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();d=expected()
    if a.check:
        if cert.read(OUT)!=d:raise ArithmeticError('retrospective height proof differs')
    else:
        if OUT.exists():raise FileExistsError('preserve height audit')
        checkpoint(OUT,d)
    print('EXACT PUBLIC COMPACT HEIGHT AUDIT',d['status_counts'],d['reported_rank_coverage'],flush=True)
