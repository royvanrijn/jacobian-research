#!/usr/bin/env python3
"""Controlled MW18 anchor calibration, then a gated balanced candidate trial.

Freeze inputs and centre lists before discovery. Run each cell under the shared
supervisor, checkpoint every chart, and replay exact hits and incremental finite
independence. Neither a miss nor an incomplete cell supplies a rank upper bound.
"""
import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sys
from time import monotonic

from pointed_quartic_search import PointedQuarticSearch, ROOT, sources
from research_runtime.search_state import raw_state, reduction_cache
from research_runtime.mw_state import MWState
from research_runtime.store import checkpoint, digest
from research_runtime.supervisor import capture, Limits

INPUT = ROOT/'artifacts/generated-results/elkies-k3-r17-extreme-anchored-mw18-specializations-h1000-v1.json'
POLICIES = ('nearest_first', 'deepest', 'diverse_deep')
CHART_COUNT = 40


def check_protocol(document):
    row = dict(document); expected = row.pop('protocol_hash')
    if digest(row) != expected: raise ArithmeticError('frozen protocol changed')


def freeze(geometry_path, centres_dir, output):
    if output.exists(): raise FileExistsError('frozen protocols are immutable')
    geometry = json.loads(geometry_path.read_text())
    if geometry['status']!='PASS' or len(geometry['covers'])!=9: raise ArithmeticError('geometry certificate is incomplete')
    all_rows = json.loads(INPUT.read_text())['candidates']
    anchors = [r for r in all_rows if r['nagao']['is_certified_anchor']]
    if len(all_rows) != 178 or len(anchors) != 5: raise ArithmeticError('panel changed')
    geometries = {r['cover_label']: r for r in geometry['covers']}
    centre_records = {}
    cases = {}
    for row in anchors:
        label = row['cover_label']; cp = centres_dir/(label+'.json')
        c = json.loads(cp.read_text())
        if c['geometry_hash'] != digest(geometries[label]): raise ArithmeticError('centre geometry mismatch')
        if c['selection']['count']!=CHART_COUNT or any(len(v)!=CHART_COUNT for v in c['original_basis_centres'].values()):
            raise ArithmeticError('centre exposure differs from the planned experiment')
        centre_records[label] = c
        points = row['specialized_points']['generic_R17']+[row['specialized_points']['cover_section']]
        cases[row['candidate_id']] = {'curve': row['raw_short_model'], 'points': points,
            'cover_label': label, 'anchor_id': row['anchor_id'], 'centres': c['original_basis_centres']}
    # Freeze a prospective roster independently of recovery outcomes: three
    # height strata per cover, deterministic hash selection inside each stratum.
    covers = json.loads((ROOT/'artifacts/generated-results/elkies-k3-r17-extreme-anchored-mw18-covers-v1.json').read_text())
    anchor_t = {c['label']: Fraction(f['native_parameter']) for ch in covers['charts'] for f in ch['fibres'] for c in f['covers']}
    h=covers['historical_rank28_anchor']; anchor_t[h['label']]=Fraction(h['native_parameter'])
    prospective=[]
    for label in sorted(geometries):
        pool = sorted((r for r in all_rows if r['cover_label']==label and Fraction(r['base_t'])!=anchor_t[label]),
                      key=lambda r:(r['nagao']['projective_height'],r['candidate_id']))
        if len(pool)<3: raise ArithmeticError('balanced cover stratum is empty')
        for k in range(3):
            stratum=pool[k*len(pool)//3:(k+1)*len(pool)//3]
            selected=min(stratum,key=lambda r:digest({'seed':'mw18-deep-centres-v1','id':r['candidate_id']}))
            prospective.append({'candidate_id':selected['candidate_id'],'cover_label':label,'height_stratum':k,
                                'parameter_height':selected['nagao']['projective_height']})
    record = {'schema':'elliptic-curves.mw18-centre-experiment.v1', 'cases':cases,
        'policies':list(POLICIES), 'geometry':geometry, 'centre_records':centre_records,
        'source_hashes': {**sources(), **{str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in
            (Path(__file__).resolve(), ROOT/'elliptic-curves/cas/prepare_mw18_deep_centres.sage',
             ROOT/'elliptic-curves/cas/research_runtime/deep_centres.py', INPUT)}},
        'limits':{'height':100000, 'seconds_per_chart':20, 'chart_count':CHART_COUNT, 'cell_wall_seconds':1800,
                  'rss_bytes':2147483648, 'prime_bound':1000},
        'coordinate_policy':'metric:16', 'enumeration_backend':'gmp-pointed-sieve',
        'ranking_rule':['total_certified_gain_descending','minimum_anchor_gain_descending','policy_name'],
        'success_gate':{'minimum_total_gain':35,'minimum_each_anchor_gain':5,'all_boxes_complete':True},
        'demonstrated_remaining_directions_per_anchor':10,
        'prospective_roster':prospective,
        'prospective_selection':'Three equal-count parameter-height strata per cover, smallest seeded candidate hash in each; exclude both points above every known anchor t, including r=infinity.',
        'blinding':'Search inputs contain only the equation, eighteen sections, their generic geometry, and fixed centres. No public complement is read by search or selection.',
        'claim_boundary':'Five anchor presentations on four distinct curves. Historical calibration only; ten demonstrated remaining directions per presentation, never fourteen.'}
    record['protocol_hash']=digest(record); checkpoint(output,record)
    print('FROZEN|'+record['protocol_hash'],flush=True)


def cell(protocol, case_id, policy, output, verify=False):
    check_protocol(protocol); case=protocol['cases'][case_id]; limits=protocol['limits']
    centres=case['centres'][policy]
    if len(centres)!=limits['chart_count']: raise ArithmeticError('chart exposure changed')
    cache=reduction_cache(); previous=json.loads(output.read_text()) if output.exists() else None
    if verify and previous is None: raise FileNotFoundError('missing cell witnesses')
    identity={'protocol_hash':protocol['protocol_hash'],'case':case_id,'policy':policy}
    if previous:
        if any(previous[k]!=v for k,v in identity.items()): raise ArithmeticError('cell checkpoint changed')
        cache.store.import_snapshot(previous['arithmetic_facts'])
        initial=MWState.from_record(previous['initial_state'],cache=cache)
        expected=raw_state(case['curve'],case['points'],cache=cache,prime_bound=limits['prime_bound'])
        if initial!=expected: raise ArithmeticError('initial state does not match case')
    else: initial=raw_state(case['curve'],case['points'],cache=cache,prime_bound=limits['prime_bound'])
    if initial.rank!=18: raise ArithmeticError('eighteen input sections not certified')
    from mod2_reduction_independence import _is_prime
    primes=[p for p in range(3,limits['prime_bound']+1) if _is_prime(p)]
    state=initial
    result={**identity,'schema':'elliptic-curves.mw18-centre-cell.v1','initial_state':initial.record(),'charts':[],
            'status':'RUNNING','source_hashes':sources()}
    retained_charts=previous['charts'] if previous else []
    for i,centre in enumerate(centres):
        if verify and i>=len(retained_charts): raise ArithmeticError('incomplete retained cell')
        start=monotonic()
        search=PointedQuarticSearch(state=initial,centre=centre,coordinate_policy=protocol['coordinate_policy'])
        if i<len(retained_charts):
            old=retained_charts[i]
            if old['centre']!=centre: raise ArithmeticError('chart centre changed')
            outcome=search.verify_record(old['search'])
        else:
            outcome=search.search(limits['height'],limits['seconds_per_chart'],checkpoint_dir=output.parent/'charts')
        record=outcome.record
        if record['height_bound']!=limits['height'] or record['timeout_seconds']!=limits['seconds_per_chart']:
            raise ArithmeticError('chart budget differs from protocol')
        before=state
        for point in outcome.curve_points: state=state.adjoin(point,cache=cache,extra_primes=primes)
        row={'centre':centre,'search':record,'state_before':before.key,'state_after':state.key,
             'certified_gain':state.rank-before.rank}
        if i<len(retained_charts) and row!=retained_charts[i]: raise ArithmeticError('admission replay differs')
        result['charts'].append(row)
        result.update({'final_state':state.record(),'arithmetic_facts':cache.store.snapshot(),
                       'certified_gain':state.rank-18})
        if not verify: checkpoint(output,result)
        print(f"MW18_CELL|{case_id}|{policy}|chart={i+1}/{len(centres)}|gain={state.rank-18}|{record['status']}|wall={monotonic()-start:.2f}",flush=True)
    result['status']='COMPLETE' if all(r['search']['status']=='bounded_search_complete' for r in result['charts']) else 'INCOMPLETE_BOXES'
    if verify:
        if result['final_state']!=previous['final_state'] or result['status']!=previous['status']:
            raise ArithmeticError('cell final replay differs')
    else: checkpoint(output,result)
    return result


def summarize(protocol,directory):
    rows=[]
    for case in protocol['cases']:
        for policy in protocol['policies']:
            p=directory/'cells'/f'{case}--{policy}.json'
            if not p.exists():continue
            d=json.loads(p.read_text())
            if d['protocol_hash']!=protocol['protocol_hash']: raise ArithmeticError('foreign cell in summary')
            rows.append({'case':case,'anchor_id':protocol['cases'][case]['anchor_id'],'policy':policy,
                'status':d['status'],'chart_count':len(d['charts']),'certified_gain':d['certified_gain'],
                'complete_boxes':sum(r['search']['status']=='bounded_search_complete' for r in d['charts']),
                'result_path':str(p.relative_to(directory)),'result_sha256':sha256(p.read_bytes()).hexdigest()})
    complete=len(rows)==len(protocol['cases'])*len(protocol['policies']) and all(r['status']=='COMPLETE' for r in rows)
    ranking=[]
    if complete:
        for policy in protocol['policies']:
            pr=[r for r in rows if r['policy']==policy]
            groups={g: [r['certified_gain'] for r in pr if r['anchor_id']==g] for g in sorted({r['anchor_id'] for r in pr})}
            ranking.append({'policy':policy,'total_gain':sum(r['certified_gain'] for r in pr),
                'minimum_anchor_gain':min(r['certified_gain'] for r in pr),'gains_by_distinct_curve':groups})
        ranking.sort(key=lambda r:(-r['total_gain'],-r['minimum_anchor_gain'],r['policy']))
    gate=protocol['success_gate']
    passed=bool(complete and ranking[0]['total_gain']>=gate['minimum_total_gain'] and ranking[0]['minimum_anchor_gain']>=gate['minimum_each_anchor_gain'])
    result={'schema':'elliptic-curves.mw18-centre-comparison.v1','protocol_hash':protocol['protocol_hash'],
        'status':'COMPLETE' if complete else 'INCOMPLETE','cells':rows,'ranking':ranking,
        'success_gate_passed':passed,'frozen_winning_policy':ranking[0]['policy'] if passed else None,
        'prospective_trial_authorized_by_gate':passed,
        'claim_boundary':'Bounded anchor recoveries, not exact ranks or evidence of fourteen existing directions beyond MW18.'}
    if protocol.get('phase')=='balanced_prospective_trial':
        result={'schema':'elliptic-curves.mw18-balanced-trial.v1','protocol_hash':protocol['protocol_hash'],
            'status':'COMPLETE' if complete else 'INCOMPLETE','cells':rows,
            'frozen_policy':protocol['policies'][0],
            'total_certified_gain':sum(r['certified_gain'] for r in rows),
            'claim_boundary':'Frozen bounded prospective experiment; additional ranks on zero-gain fibres remain UNKNOWN.'}
    checkpoint(directory/'summary.json',result); return result


def run(protocol_path,directory,verify=False):
    protocol=json.loads(protocol_path.read_text());check_protocol(protocol)
    for case_id in protocol['cases']:
        for policy in protocol['policies']:
            output=directory/'cells'/f'{case_id}--{policy}.json'
            if not verify and output.exists() and json.loads(output.read_text())['status']=='COMPLETE':continue
            command=[sys.executable,str(Path(__file__).resolve()),'--protocol',str(protocol_path),'--cell',case_id,
                     '--policy',policy,'--output',str(output)]
            if verify:command.append('--verify')
            r=capture(command,limits=Limits(protocol['limits']['cell_wall_seconds'],protocol['limits']['rss_bytes']),
                      log_path=directory/('replay-logs' if verify else 'logs')/f'{case_id}--{policy}.log',check=False)
            print(f"MW18_EXPERIMENT|{case_id}|{policy}|returncode={r.returncode}",flush=True)
            summarize(protocol,directory)
            if r.returncode: raise RuntimeError('cell failed; retain its partial checkpoint and log')
    result=summarize(protocol,directory)
    if verify:
        checkpoint(directory/'replay.json',{'status':'PASS_EXACT_CHARTS_AND_INCREMENTAL_INDEPENDENCE',
            'protocol_hash':protocol['protocol_hash'],'cell_count':len(result['cells']),
            'chart_count':sum(r['chart_count'] for r in result['cells']),
            'summary_sha256':sha256((directory/'summary.json').read_bytes()).hexdigest(),
            'enumeration_repeated':False})
    return result


def freeze_prospective(protocol_path, directory, centres_dir, output):
    """Only a completed and independently replayed calibration opens this gate."""
    if output.exists(): raise FileExistsError('prospective protocols are immutable')
    parent=json.loads(protocol_path.read_text());check_protocol(parent)
    summary=summarize(parent,directory)
    replay=json.loads((directory/'replay.json').read_text())
    if (not summary['success_gate_passed'] or replay['protocol_hash']!=parent['protocol_hash']
        or replay['summary_sha256']!=sha256((directory/'summary.json').read_bytes()).hexdigest()):
        raise ArithmeticError('the replayed anchor success gate has not passed')
    winning=summary['frozen_winning_policy']; cases={}
    candidates={r['candidate_id']:r for r in json.loads(INPUT.read_text())['candidates']}
    geometries={r['cover_label']:r for r in parent['geometry']['covers']}
    for selected in parent['prospective_roster']:
        row=candidates[selected['candidate_id']]; label=row['cover_label']
        centres=json.loads((centres_dir/(label+'.json')).read_text())
        if centres['geometry_hash']!=digest(geometries[label]): raise ArithmeticError('prospective geometry changed')
        cases[row['candidate_id']]={'curve':row['raw_short_model'],
            'points':row['specialized_points']['generic_R17']+[row['specialized_points']['cover_section']],
            'cover_label':label,'anchor_id':row['anchor_id'],
            'centres':{winning:centres['original_basis_centres'][winning]}}
    # The same runner and limits apply; this is an outcome-free population
    # experiment. A prospective zero gain can never change the frozen winner.
    record={k:v for k,v in parent.items() if k not in ('protocol_hash','cases','policies','centre_records')}
    record.update({'cases':cases,'policies':[winning],'phase':'balanced_prospective_trial',
        'parent_protocol_hash':parent['protocol_hash'],'anchor_summary':summary,
        'claim_boundary':'Twenty-seven frozen non-anchor candidates, three height strata per cover. Bounded misses leave their additional rank UNKNOWN.'})
    record['protocol_hash']=digest(record);checkpoint(output,record)
    print('PROSPECTIVE_FROZEN|'+record['protocol_hash'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--freeze',type=Path);p.add_argument('--geometry',type=Path);p.add_argument('--centres',type=Path)
    p.add_argument('--protocol',type=Path);p.add_argument('--directory',type=Path)
    p.add_argument('--cell');p.add_argument('--policy',choices=POLICIES);p.add_argument('--output',type=Path)
    p.add_argument('--verify',action='store_true');p.add_argument('--summarize',action='store_true')
    p.add_argument('--freeze-prospective',type=Path)
    a=p.parse_args()
    if a.freeze:freeze(a.geometry,a.centres,a.freeze)
    elif a.freeze_prospective:freeze_prospective(a.protocol,a.directory,a.centres,a.freeze_prospective)
    elif a.cell:cell(json.loads(a.protocol.read_text()),a.cell,a.policy,a.output,a.verify)
    elif a.summarize:print(json.dumps(summarize(json.loads(a.protocol.read_text()),a.directory),indent=2))
    else:run(a.protocol,a.directory,a.verify)
