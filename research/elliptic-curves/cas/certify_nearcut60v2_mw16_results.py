#!/usr/bin/env python3
"""Export independently rechecked mod2 point bases for the completed near-finalist cohort."""
import argparse,json
from pathlib import Path
import certify_compact_r17_candidates as cert
import compact_mw16_specialization as spec
import certify_corrected60_mw16_results as existing
import nearcut60v2_mw16_pari_batch as batch
ROOT=batch.ROOT;ART=batch.ART;OUT=ART/'nearcut60v2_mw16_results_v1.json'
NOVELTY=ART/'nearcut60v2_mw16_novelty_v1.json'
REPORT=ART/'nearcut60v2_mw16_experiment_v1.json'


def sources():
    paths=[Path(__file__).resolve(),Path(existing.__file__),Path(cert.__file__),Path(spec.__file__),spec.ATLAS,
           ROOT/'elliptic-curves/cas/memory_rank_certificate.py',NOVELTY,REPORT]
    return {str(p.relative_to(ROOT)):cert.hashed(p) for p in paths}


def expected():
    gate=cert.read(batch.extension.D/'novelty-controller-v2/ledger.json')
    if gate['status']!='PASS':raise ArithmeticError('terminal novelty replay required')
    p=batch.protocol();novelty=cert.read(NOVELTY);report=cert.read(REPORT)
    if novelty['status']!='PASS' or report['status']!='COMPLETE_FIXED_RETAINED_TRIAL' or len(novelty['rows'])!=60 or novelty['within_cohort_isomorphic_pairs']:raise ArithmeticError('complete distinct certified cohort required')
    for data in (novelty,report):
        if any(cert.hashed(ROOT/n)!=h for n,h in data['sources'].items()):raise ArithmeticError('terminal source changed')
    families={f['fibration_id']:f for f in cert.read(spec.ATLAS)['families']};rows=[];bindings={}
    for chosen,n,r in zip(p['rows'],novelty['rows'],report['rows']):
        if chosen['id']!=n['id'] or n['id']!=r['id'] or n['rank_lower_bound'] is None:raise ArithmeticError('certified roster differs')
        bound=n['point_proof_binding'];path=ROOT/bound['cloud_path'];oddpath=ROOT/bound['odd_path']
        if cert.hashed(path)!=bound['cloud_sha256'] or cert.hashed(oddpath)!=bound['odd_sha256']:raise ArithmeticError('point proof changed')
        cloud=cert.read(path);odd=cert.read(oddpath);proofpath=batch.D/n['id']/'proof-input.json';data=cert.read(proofpath)
        if cloud['input_sha256']!=cert.hashed(proofpath) or cloud['curve']!=data['curve'] or not cert.isomorphic(n['curve'],cloud['curve']):raise ArithmeticError('rank proof equation provenance differs')
        if cloud['rank_lower_bound']!=r['rank_lower_bound'] or any(a['finite_column_rank']>cloud['rank_lower_bound'] for a in odd['audits']):raise ArithmeticError('stronger odd bound requires separate export; do not weaken claim')
        row={'id':n['id'],'family':n['family'],'parameter':n['parameter'],'band':n['band'],'arms':[n['arm']],
             'retained_index':chosen['retained_rank'],'curve':cloud['curve'],'points':cloud['independent_points'],
             'generic_points':data['generic_points'],'family_to_curve_scale_u':data['family_to_curve_scale_u'],
             'rank_certificate':cloud['rank_certificate'],'rank_lower_bound':cloud['rank_lower_bound'],
             'icarm_matches':n['catalogue_matches'],'previous_matches':n['previous_matches'],
             'attempted_charts':r['attempted_boxes'],'completed_boxes':r['completed_boxes'],'search_status':r['worker_status'],
             'complete_cloud_certificate':{'path':str(path.relative_to(ROOT)),'sha256':cert.hashed(path)}}
        existing.verify(row,families,novelty['catalogue']['equations'],novelty['previous_equations']);rows.append(row)
        for q in (path,oddpath,proofpath):bindings[str(q.relative_to(ROOT))]=cert.hashed(q)
    return {'schema':'elliptic-curves.nearcut60v2-results.v1','status':'PASS','sources':sources(),'point_bindings':bindings,
        'curves':rows,'within_batch_isomorphic_pairs':novelty['within_cohort_isomorphic_pairs'],
        'catalogue':novelty['catalogue'],'previous_equations':novelty['previous_equations'],
        'claim_boundary':'All60 retained near-finalist curves have rechecked independent mod2 point bases, exact generic transports and post-terminal Q-isomorphism comparisons with593 pinned catalogue and1405 prior equations. All odd bounds are covered by the exported mod2 bounds. This export follows the completed fixed experiment; it changes no allocation, score, budget or policy criterion. No exact rank or universal novelty.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();d=expected()
    if a.check:
        if cert.read(OUT)!=d:raise ArithmeticError('point export replay differs')
    else:
        if OUT.exists():raise FileExistsError('preserve near-finalist point export')
        cert.write(OUT,d)
    print('NEARCUT60V2 EXACT POINT EXPORT PASS',len(d['curves']),flush=True)
