#!/usr/bin/env python3
"""Fail-closed certified discovery yield, completed exposure and matched policies."""
import argparse
from collections import Counter
from pathlib import Path
import strata60_mw16_pari_batch as batch
import verify_strata60_mw16_points as proof
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
ROOT=batch.ROOT; D=batch.D; ART=batch.ART
OUT=ART/'retained_mw16_score_strata_experiment_v1.json'
ARMS=('top','moderate','lower')


def summarize(rows):
    gains=[r['certified_gain'] for r in rows]
    complete=[r['completed_boxes'] for r in rows]
    gain=sum(x for x in gains if x is not None)
    seconds=sum(r['discovery_worker_seconds'] for r in rows)
    inclusive=seconds+sum(r['verification_seconds'] for r in rows)
    all_gains=all(x is not None for x in gains)
    all_exposure=all(x is not None for x in complete)
    boxes=sum(x for x in complete if x is not None)
    return {'allocated_curves':len(rows),'certified_gain_sum_known':gain,
        'unresolved_gain_curves':sum(x is None for x in gains),
        'certified_gain_distribution':dict(sorted(Counter(str(x) for x in gains if x is not None).items())),
        'completed_boxes_known':boxes,'unresolved_exposure_curves':sum(x is None for x in complete),
        'allocated_boxes':43*len(rows),'discovery_worker_seconds':seconds,
        'including_verification_seconds':inclusive,
        'gain_per_discovery_second':gain/seconds if all_gains and seconds>0 else None,
        'gain_per_inclusive_second':gain/inclusive if all_gains and inclusive>0 else None,
        'gain_per_completed_box':gain/boxes if all_gains and all_exposure and boxes>0 else None,
        'completion_fraction':boxes/(43*len(rows)) if all_exposure and rows else None,
        'censored_or_failed_workers':sum(r['worker_status']!='COMPLETE_DECLARED_POINT_ATTEMPT' for r in rows),
        'known_gain_sum_boundary':'Sum of certified directions in retained evidence. Unresolved curves are explicitly counted; zero certified gain does not establish absence of more points.'}


def comparison(rows):
    top=summarize([r for r in rows if r['arm']=='top']); pooled=summarize(rows)
    keys=('gain_per_discovery_second','gain_per_completed_box','completion_fraction')
    evaluable=all(a[k] is not None for a in (top,pooled) for k in keys)
    leave=[]
    for family in sorted({r['family'] for r in rows}):
        subset=[r for r in rows if r['family']!=family]
        t=summarize([r for r in subset if r['arm']=='top'])['gain_per_discovery_second']
        q=summarize(subset)['gain_per_discovery_second']
        leave.append({'omitted_family':family,'top_gain_per_second':t,'pooled_gain_per_second':q,
                      'pooled_advantage':t is not None and q is not None and q>t})
    diversify=evaluable and pooled[keys[0]]>top[keys[0]] and pooled[keys[1]]>top[keys[1]] and pooled[keys[2]]>=top[keys[2]] and bool(leave) and all(x['pooled_advantage'] for x in leave)
    top_favoured=evaluable and top[keys[0]]>pooled[keys[0]] and top[keys[1]]>=pooled[keys[1]] and top[keys[2]]>=pooled[keys[2]]
    return {'top_only':top,'equal_three_arm_portfolio':pooled,'leave_one_family_out':leave,
            'recommendation':'DIVERSIFY_WITHIN_RETAINED_POOL' if diversify else 'TOP_FAVOURED_IN_THIS_FINITE_SAMPLE' if top_favoured else 'INCONCLUSIVE',
            'diversification_criterion_met':bool(diversify),'all_endpoints_resolved':evaluable,
            'scope':'Overlapping policy estimates from the same matched curves, not independent samples or a rank-density theorem. No automatic follow-on sweep.'}


def expected():
    p=batch.protocol(); paths={Path(__file__).resolve(),D/'protocol.json',batch.extension.OUT}
    control=batch.extension.D/'point-controller/protocol.json';fixed=cert.read(control);paths.add(control)
    for name,h in fixed['sources'].items():
        source=ROOT/name;paths.add(source)
        if cert.hashed(source)!=h:raise ArithmeticError('frozen comparison source changed')
    def read(path):
        paths.add(path);return cert.read(path)
    def supervision(s):
        if s['outcome']=='running':raise ArithmeticError('terminal supervision required')
        path=Path(s['log']);paths.add(path)
        if cert.hashed(path)!=s['log_sha256']:raise ArithmeticError('execution log changed')
        if s['wall_seconds']<0:raise ArithmeticError('negative execution time')
        return s['wall_seconds']
    maps=read(D/'map-ledger.json')
    if maps['status']=='RUNNING':raise ArithmeticError('map stage not terminal')
    baselines=read(D/'baseline-ledger.json') if (D/'baseline-ledger.json').exists() else {'rows':[]}
    ledger=read(D/'ledger.json') if (D/'ledger.json').exists() else {'rows':[]}
    verifications=read(D/'verification-ledger.json') if (D/'verification-ledger.json').exists() else {'rows':[]}
    if any(x.get('status')=='RUNNING' for x in (baselines,ledger,verifications)):
        raise ArithmeticError('all allocated stages must be terminal')
    if ledger['rows'] and (ledger.get('status')!='COMPLETE_FIXED_BATCH_ATTEMPTS' or verifications.get('status')!='COMPLETE_ALL_ALLOCATED_VERIFICATIONS'):
        raise ArithmeticError('all terminal searches and proof attempts required')
    bymap={r['id']:r for r in maps['rows']};bybase={r['id']:r for r in baselines['rows']}
    bysearch={r['id']:r for r in ledger['rows']};byproof={r['id']:r for r in verifications['rows']}
    rows=[]
    for selected in p['rows']:
        row={k:selected[k] for k in ('id','arm','block','family','band','sign','parameter','retained_rank','j_height','parameter_height')}
        folder,cloud,odd=proof.paths(row);mapped=bymap[row['id']]
        map_seconds=supervision(mapped['supervision']);baseline_seconds=0
        if row['id'] in bybase:baseline_seconds=supervision(bybase[row['id']]['supervision'])
        searched=bysearch.get(row['id']);verified=byproof.get(row['id'])
        seconds=supervision(searched['point_supervision']) if searched else 0
        verify_seconds=baseline_seconds
        exposure_ok=False; gain=None; lower=None; charts=None; attempted=None; cpu=None
        verification_status='NOT_RUN'
        if verified:
            if read(folder/'verification.json')!=verified:raise ArithmeticError('verification ledger changed')
            verification_status=verified['status']
            for name,h in verified['inputs'].items():
                path=folder/name;paths.add(path)
                if cert.hashed(path)!=h:raise ArithmeticError('verified input changed')
            for step in verified['steps']:verify_seconds+=supervision(step['supervision'])
            exposure_ok=any(s['name']=='geometry' and s['status']=='PASS' for s in verified['steps'])
            if exposure_ok:
                exposure=read(folder/'exposure.json');charts=[r['search'] for r in exposure['charts'] if r['status']=='RETURNED'];attempted=len(exposure['charts'])
                cpu=sum(r['search_cpu_ms'] for r in charts if r['search_cpu_ms'] is not None)
            if verified['status']=='PASS':
                if len(verified['steps'])!=7 or any(s['status']!='PASS' for s in verified['steps']):raise ArithmeticError('all7 exact proof stages required')
                for path,h in ((folder/'proof-input.json',verified['proof_input_sha256']),(cloud,verified['cloud_sha256']),(odd,verified['odd_sha256'])):
                    paths.add(path)
                    if cert.hashed(path)!=h:raise ArithmeticError('point proof artifact changed')
                c=read(cloud);o=read(odd)
                if c['input_sha256']!=verified['proof_input_sha256'] or o['input_sha256']!=verified['cloud_sha256']:raise ArithmeticError('point proof chain differs')
                lower=max(c['rank_lower_bound'],*(a['finite_column_rank'] for a in o['audits']))
                if lower!=verified['rank_lower_bound'] or lower-16!=verified['certified_gain']:raise ArithmeticError('certified gain differs')
                gain=lower-16
        coefficients=None
        if mapped['status']=='PASS':
            path=folder/'maps.json';m=read(path)
            if cert.hashed(path)!=mapped['maps_sha256']:raise ArithmeticError('map artifact changed')
            coefficients=[{'numerator_bits':abs(cert.F(v).numerator).bit_length(),'denominator_bits':cert.F(v).denominator.bit_length()} for v in m['curve']]
        row.update(verification_status=verification_status,rank_lower_bound=lower,certified_gain=gain,
            allocated_boxes=43,attempted_boxes=attempted,returned_charts=len(charts) if charts is not None else None,
            completed_boxes=sum(r['status']=='bounded_search_complete' for r in charts) if charts is not None else None,
            chart_timeouts=sum(r['status']=='bounded_search_timeout' for r in charts) if charts is not None else None,
            backend_failures=sum(r['status']=='backend_failure' for r in charts) if charts is not None else None,
            in_flight_at_termination=attempted-len(charts) if charts is not None else None,
            unattempted_boxes=43-attempted if attempted is not None else None,
            worker_status=searched['status'] if searched else 'NOT_RUN_PRESEARCH_GATE',
            map_seconds=map_seconds,point_worker_seconds=seconds,discovery_worker_seconds=map_seconds+seconds,
            verification_seconds=verify_seconds,known_search_cpu_ms=cpu,
            charts_without_cpu_time=attempted-sum(r['search_cpu_ms'] is not None for r in charts) if charts is not None else None,
            exposure_verified=exposure_ok,search_model_coefficient_sizes=coefficients)
        rows.append(row)
    if len(rows)!=60 or Counter(r['arm'] for r in rows)!=Counter({a:20 for a in ARMS}):raise ArithmeticError('complete fixed allocation differs')
    return {'schema':'elliptic-curves.retained-score-strata-experiment.v1',
        'status':'COMPLETE_FIXED_COMPARISON' if all(r['certified_gain'] is not None for r in rows) else 'COMPLETE_ACCOUNTING_WITH_UNRESOLVED_OUTCOMES',
        'sources':{str(path.relative_to(ROOT)):cert.hashed(path) for path in sorted(paths)},'rows':rows,
        'arms':{a:summarize([r for r in rows if r['arm']==a]) for a in ARMS},
        'policy_comparison':comparison(rows),
        'families':{f:{a:summarize([r for r in rows if r['family']==f and r['arm']==a]) for a in ARMS} for f in sorted({r['family'] for r in rows})},
        'matched_triplets':{b:[r['id'] for r in rows if r['block']==b] for b in sorted({r['block'] for r in rows})},
        'shared_sunk_costs':{'retained_population_addresses':3059269468,'retained_candidates':1310720,'charged_to_arms':False,
                             'scope':'Previously completed trace-table construction and corrected population selection. Prospective matching/report orchestration costs are shared, not assigned to an arm.'},
        'validation':'Not computed for this comparison;65537..131071 remain separate and never affect allocation or the policy criterion.',
        'claim_boundary':'Certified lower-bound gains over the independent generic16 subgroup in this fixed matched sample. Unresolved gains and exposure remain null, never zero. Worker wall time includes search, startup, admission and interruption cleanup; map time is charged to discovery and baseline/terminal proofs to verification. GP CPU time is partial when a chart lacks its timing marker. Completed boxes trust pinned PARI execution. No absence, exact-rank, universal novelty, sampling-density or causal optimality theorem; no larger sweep authorized.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();data=expected()
    if a.check:
        if cert.read(OUT)!=data:raise ArithmeticError('comparison accounting replay differs')
    else:
        if OUT.exists():raise FileExistsError('preserve comparison outcome')
        checkpoint(OUT,data)
    print('STRATIFIED COMPARISON',data['status'],data['policy_comparison']['recommendation'],flush=True)
