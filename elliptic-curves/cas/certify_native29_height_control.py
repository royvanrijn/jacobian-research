#!/usr/bin/env python3
"""Exact finite proofs and honest coverage summary for the matched height controls."""
import argparse,json
from pathlib import Path
import certify_compact_r17_candidates as cert
from memory_rank_certificate import checked_rank
from research_runtime.store import checkpoint
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';LOCAL=ROOT/'artifacts/local/elliptic-curves';OUT=ART/'native29_height_control_v1.json'
FOLDERS={100000:LOCAL/'native11952-height-pair-v1/100000',125000:LOCAL/'native11952-height125-control-v1/125000',1000000:LOCAL/'native11952-height-pair-v1/1000000'}
def expected():
    replay=ART/'native11952_translated_visibility_replay_v1.json';audit=ART/'native11952_translated_visibility_v1.json';old=cert.read(audit);rp=cert.read(replay)
    if rp['status']!='PASS' or rp['sources'][str(audit.relative_to(ROOT))]!=cert.hashed(audit) or rp['exact_group_words_checked']!=196:raise ArithmeticError('exact196 group replay missing')
    inputs={str(p.relative_to(ROOT)):cert.hashed(p) for p in (Path(__file__).resolve(),replay,audit,ROOT/'elliptic-curves/cas/memory_rank_certificate.py')};rows=[];maps=[]
    for height,folder in FOLDERS.items():
        path=folder/'result.json';vp=folder/'verification.json';d,v=cert.read(path),cert.read(vp);cp=ART/f'native11952_height{height}_mod2_v1.json';c=cert.read(cp);proof=c['rank_certificate'];model=tuple(map(cert.F,c['curve']));points=[tuple(map(cert.F,p)) for p in c['independent_points']]
        if v['status']!='PASS' or c['input_sha256']!=cert.hashed(path) or d['height']!=height or len(d['charts'])!=49 or c['rank_lower_bound']!=len(points):raise ArithmeticError('terminal cloud binding differs')
        actual=checked_rank(model,points,[s['prime'] for s in proof['signatures']],proof['no_rational_2_torsion_prime'])
        if json.loads(json.dumps(actual))!=proof:raise ArithmeticError('control finite proof differs')
        for name in ('worker','replay'):
            sp=folder/(name+'.supervisor.json');s=cert.read(sp)
            if s['outcome']!='completed' or s['returncode']!=0 or cert.hashed(Path(s['log']))!=s['log_sha256']:raise ArithmeticError('history/supervision replay differs')
            inputs[str(sp.relative_to(ROOT))]=cert.hashed(sp)
        if d['generic_points']!=old['basis'][:17] or d['initial_dimension']!=17 or d['curve']!=old['curve']:raise ArithmeticError('generic-only starting input differs')
        counts={status:sum(r['search']['status']==status for r in d['charts']) for status in sorted({r['search']['status'] for r in d['charts']})};maps.append((d['maps_sha256'],d['centres']));rows.append({'height':height,'point_count':len(c['points']),'rank_lower_bound_from_this_attempt':c['rank_lower_bound'],'coverage_counts':counts,'worker_wall_seconds':cert.read(folder/'worker.supervisor.json')['wall_seconds'],'certificate':str(cp.relative_to(ROOT)),'certificate_sha256':cert.hashed(cp)})
        inputs.update({str(p.relative_to(ROOT)):cert.hashed(p) for p in (path,vp,cp)})
    if any(m!=maps[0] for m in maps):raise ArithmeticError('height arms changed maps/centres')
    if [r['rank_lower_bound_from_this_attempt'] for r in rows[:2]]!=[27,28] or any(r['coverage_counts']!={'bounded_search_complete':49} for r in rows[:2]):raise ArithmeticError('completed27-to28 control gate differs')
    return {'schema':'elliptic-curves.native29-height-control.v1','status':'PASS','sources':inputs,'exact_retrospective_translations':196,'best_missing_direction_heights':[r['visibility']['minimum_affine_height'] for r in rp['best_per_direction']],'identical_frozen_maps_and_centres':True,'generic_starting_rank':17,'rows':rows,'claim_boundary':'A known rank29 curve, not a new curve: exact retrospective group words justify the height experiment. Identical generic17-only maps certify27 at100000 and28 at125000, both49 boxes complete. The million-height arm has49 censored attempts and its retained generic17 certificate is not a rank loss, upper bound or search exclusion. These data justify a bounded prospective height test only; no universal recovery guarantee or exact rank.'}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();d=expected()
    if a.check:
        if cert.read(OUT)!=d:raise ArithmeticError('height-control certificate differs')
    else:
        if OUT.exists():raise FileExistsError('preserve height-control summary')
        checkpoint(OUT,d)
    print('REPLAYED NATIVE29 HEIGHT CONTROL',[(r['height'],r['rank_lower_bound_from_this_attempt'],r['coverage_counts']) for r in d['rows']],flush=True)
