#!/usr/bin/env python3
"""Bind known29 recovery and the fixed six-curve prospective parity-policy test."""
import argparse,json
from collections import Counter
from pathlib import Path
import certify_compact_r17_candidates as cert
import native28_specialized_parity_control as control
import new27_specialized_parity_six as batch
from research_runtime.store import checkpoint
ROOT=control.ROOT;ART=control.ART;OUT=ART/'specialized_parity_experiment_v1.json'

def expected():
    paths={Path(__file__).resolve()}
    def read(p):paths.add(p);return cert.read(p)
    def terminal(p):
        s=read(p)
        if s['outcome']!='completed' or s['returncode']!=0:raise ArithmeticError('required exact replay failed or censored')
    p=read(control.D/'protocol.json');m=read(control.D/'maps.json');r=read(control.D/'result.json');coverage=read(ART/'native28_specialized_parity_adaptive_coverage_v1.json')
    for label in ('replay','cloud-audit','geometry'):terminal(control.D/(label+'.supervisor.json'))
    if len(m['sample'])!=2048 or len(m['centres'])!=49 or r['rank_lower_bound']!=29 or len(r['charts'])!=18 or r['status']!='TARGET_REACHED_PENDING_INDEPENDENT_REPLAY' or any(c['search']['status']!='bounded_search_complete' for c in r['charts']):raise ArithmeticError('known29 control recovery differs')
    if coverage['mod2_lower_bound']!=29 or coverage['odd_modulus_lower_bounds']!={'3':29,'5':29} or coverage['retained_point_count']!=73:raise ArithmeticError('complete control cloud differs')
    for name,h in coverage['point_certificates'].items():
        read(ROOT/name)
        if cert.hashed(ROOT/name)!=h:raise ArithmeticError('control cloud certificate changed')
    bp=read(batch.BATCH/'protocol.json');v=read(batch.BATCH/'verification-ledger.json');terminal(batch.BATCH/'geometry.supervisor.json')
    if v['status']!='PASS' or len(v['rows'])!=6 or [r['id'] for r in bp['rows']]!=[r['id'] for r in v['rows']]:raise ArithmeticError('complete fixed six-curve proof roster required')
    index=read(batch.INDEX);expected_ids=sorted(r['id'] for r in index['curves'] if r['rank_lower_bound']==27)
    if [r['id'] for r in bp['rows']]!=expected_ids:raise ArithmeticError('all and only inventory27 curves required')
    rows=[]
    for row,q in zip(bp['rows'],v['rows']):
        d=batch.BATCH/row['id'];raw=read(d/'result.json');maps=read(d/'maps.json');seed=read(d/'seed.json');ell=read(ROOT/q['modl_certificate']);cloud=read(ROOT/q['mod2_certificate'])
        for label in ('maps','worker','replay'):terminal(d/(label+'.supervisor.json'))
        if q['input_sha256']!=cert.hashed(d/'result.json') or q['mod2_sha256']!=cert.hashed(ROOT/q['mod2_certificate']) or q['modl_sha256']!=cert.hashed(ROOT/q['modl_certificate']):raise ArithmeticError('prospective certificate input changed')
        n=len(raw['charts']);bound=cloud['rank_lower_bound'];completed=sum(c['search']['status']=='bounded_search_complete' for c in raw['charts'])
        if len(maps['sample'])!=2048 or len(maps['rows'])!=49 or not 1<=n<=49 or q['attempted_charts']!=n or q['completed_boxes']!=completed or q['rank_lower_bound']!=bound or raw['initial_dimension']!=27 or len(seed['points'])!=27:raise ArithmeticError('fixed sampled exposure differs')
        if any(c['search']['height_bound']!=125000 or c['search']['timeout_seconds']!=10 for c in raw['charts']):raise ArithmeticError('uniform prospective budget differs')
        if raw['status']=='COMPLETE_DECLARED_ADAPTIVE_ATTEMPT':
            if n!=49 or raw['rank_lower_bound']>=28:raise ArithmeticError('complete no-target status differs')
        elif raw['status']=='TARGET_REACHED_PENDING_INDEPENDENT_REPLAY':
            if bound<28 or any(c['rank_lower_bound']>=28 for c in raw['charts'][:-1]):raise ArithmeticError('first target stop differs')
        else:raise ArithmeticError('nonterminal prospective result')
        if {str(a['modulus']):a['finite_column_rank'] for a in ell['audits']}!=q['odd_modulus_lower_bounds']:raise ArithmeticError('odd-modulus bound differs')
        rows.append({**row,'initial_rank_lower_bound':27,'rank_lower_bound':bound,'attempted_charts':n,'completed_boxes':completed,'retained_points':len(cloud['points']),'odd_modulus_lower_bounds':q['odd_modulus_lower_bounds'],'status':raw['status']})
    return {'schema':'elliptic-curves.specialized-parity-experiment.v1','status':'PASS','sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in sorted(paths)},'control':{'initial_rank_lower_bound':28,'final_rank_lower_bound':29,'sampled_parities':2048,'selected_charts':49,'completed_boxes':18,'unsearched_after_target':31,'height':125000,'all_retained_points':73,'moduli':[2,3,5],'boundary':'A known29 curve selected using prior control history. Its public unrecovered point and oracle words were excluded from this new centre selection and execution. This is recovery, not new-curve discovery or independent population validation.'},'prospective_curves':rows,'prospective_sampled_parities':6*2048,'prospective_selected_charts':294,'prospective_attempted_charts':sum(r['attempted_charts'] for r in rows),'prospective_completed_boxes':sum(r['completed_boxes'] for r in rows),'rank_lower_bound_counts':dict(sorted(Counter(r['rank_lower_bound'] for r in rows).items())),'claim_boundary':'A finite sampled specialized-parity policy using each curve own certified point subgroup. All metric changes, representative parities and rational chart/point identities replay exactly; numerical canonical heights and CVP schedule choices without certifying global coset or covering optimality. The known control success does not guarantee prospective recovery. New bounds require their exact independent point certificates and separate inventory promotion. No rank upper bound, exact rank, saturation, point-absence or universal-novelty conclusion.'}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();d=expected()
    if a.check:
        if cert.read(OUT)!=json.loads(json.dumps(d)):raise ArithmeticError('specialized-parity report differs')
    else:
        if OUT.exists():raise FileExistsError('preserve specialized-parity report')
        checkpoint(OUT,d)
    print('SPECIALIZED PARITY REPORT PASS',d['rank_lower_bound_counts'],d['prospective_completed_boxes'],'completed boxes',flush=True)
