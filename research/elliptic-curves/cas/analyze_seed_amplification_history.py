#!/usr/bin/env python3
"""Retained-history analysis only. No Sage arithmetic, search or subprocesses.

Use the existing Sage Python (NumPy/SciPy/SymPy installed). `extract` binds
retained inputs; `analyze` uses the exported rows, never live search state.
"""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
import csv
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import gzip
import json
import math
from pathlib import Path
import sqlite3
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT/'artifacts/generated-results/elliptic-curves'
LOCAL = ROOT/'artifacts/local/elliptic-curves'
OUT = ART/'seed_amplification_history_v1'
INPUTS = {}
GAPS = []
EVENTS = []
LATTICES = {}
PRIME_FEATURES = []
LEGACY = '''compact192_r17 outer48_r17 joint64_r17 annulus64_r17
full11952_64_r17 full11952_late64_r17 retention24_r17 next24_r17
higher24_r17 productfirst24_r17 scaled13_24_r17 splitnode24_r17 skew8_r17
product22_r17 discarded12_r17 fresh60_mw16 outer60_mw16 higher60_mw16
broad60_mw16 corrected60_mw16 strata60_mw16 nearcut60v2_mw16
extended20_mw16 prospective_mw16 prospective_mw16_wide
prospective_mw16_next12 compact_r17_wide fresh_r17_paired'''.split()

def sha(b): return hashlib.sha256(b).hexdigest()
def rel(p):
    try: return str(Path(p).resolve().relative_to(ROOT))
    except ValueError: return str(p)
def read(p, expected=None):
    p = Path(p)
    b = p.read_bytes(); h = sha(b)
    if expected is not None and h != expected:
        raise ValueError(f'Binding changed: {p}: {h} != {expected}')
    if rel(p) in INPUTS: assert INPUTS[rel(p)] == h, f'Input changed during read: {p}'
    INPUTS[rel(p)] = h
    return json.loads(b)
@lru_cache(maxsize=64)
def small(p): return read(p)
def write(p, d):
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(d, indent=2, sort_keys=True, allow_nan=False)+'\n')
def source(p, expected=None):
    p = Path(p)
    if not p.is_file():
        GAPS.append({'path':rel(p),'status':'MISSING_RETAINED_ARTIFACT_NOT_REBUILT'})
        return None
    return read(p, expected)
def numeric_bits(x):
    x = F(x)
    return max(abs(x.numerator).bit_length(),x.denominator.bit_length())
def model_features(model):
    a1,a2,a3,a4,a6 = map(F,model)
    b2=a1*a1+4*a2; b4=2*a4+a1*a3; b6=a3*a3+4*a6
    b8=a1*a1*a6+4*a2*a6-a1*a3*a4+a2*a3*a3-a4*a4
    delta=-b2*b2*b8-8*b4**3-27*b6*b6+9*b2*b4*b6
    c4=b2*b2-24*b4
    assert delta
    # j groups are deliberately conservative leakage groups, NOT Q-isomorphism.
    # Twists in the same j group are withheld together, never deduplicated as curves.
    return {'displayed_model_bits':max(map(numeric_bits,model)),
            'displayed_discriminant_bits':numeric_bits(delta),
            'j_leakage_group':sha(str(c4**3/delta).encode())}

def lattice(family, gram, path, surface):
    import sympy as sp
    key = sha(json.dumps(gram,sort_keys=True).encode())
    mat = sp.Matrix([[sp.Rational(str(x)) for x in row] for row in gram])
    LATTICES[family] = {'family':family,'surface':surface,'source':rel(path),
        'gram_sha256':key,'gram':gram,'rank':len(gram),'determinant':str(mat.det()),
        'diagonal_min':str(min(mat.diagonal())),'diagonal_max':str(max(mat.diagonal())),
        'scope':'Recorded generic Mordell-Weil Gram; not the integral Neron-Severi frame.'}

def base_row(cohort, cid, family, parameter, model, generic, final):
    t = F(parameter)
    d = {'cohort':cohort,'id':f'{cohort}:{cid}','source_id':cid,'family':family,
         'surface':'X1092' if family.startswith('x1092') or family=='det1092' else 'X948',
         'parameter':str(t),'numerator':t.numerator,'denominator':t.denominator,
         'parameter_height':max(abs(t.numerator),t.denominator),
         'generic_rank':generic,'final_rank':final,'initial_rank':None,
         'final_gain':None if final is None or generic is None else final-generic,
         'rank_evidence':'RETAINED_CERTIFICATE_NOT_NEW_ARITHMETIC_REPLAY',
         'cpu_seconds':None,'elapsed_seconds':None,'timing_scope':'UNKNOWN',
         'known_equation':False,'previous_equation':False,'predictive_eligible':True,
         'all_frozen_features':{},'model':model,'calls':None,
         'trajectory_complete':False,'censored':None,'endpoint_status':'UNKNOWN',
         'adaptive_followup':'UNKNOWN','source':''}
    if model: d.update(model_features(model))
    if family in LATTICES:
        d['generic_gram_sha256']=LATTICES[family]['gram_sha256']
        d['generic_gram_source']=LATTICES[family]['source']
    return d

def attach_events(row, events, complete, times=None):
    """Same-call clouds cannot certify a later amplification event."""
    generic, final = row['generic_rank'], row['final_rank']
    if generic is None or final is None:
        row.update(seed=None,further_after_seed=None,deep_amplifier=None,
                   classification='UNKNOWN',post_seed_exposure=None)
        return
    row['seed']=int(final>generic)
    row['extra_direction_after_first']=int(final>=generic+2)
    row['endpoint_bucket']=('no seed' if final==generic else
        'seed only (endpoint +1/+2)' if final<=generic+2 else
        'endpoint >=24' if final>=24 else 'endpoint +3 or more')
    bycall=defaultdict(int)
    for e in events:
        assert e['after']>e['before']
        assert 0<=e['call']<=row['calls'], (row['id'],e,row['calls'])
        bycall[e['call']]+=e['after']-e['before']
    if complete: assert sum(bycall.values())==final-generic, row['id']
    row['trajectory_complete']=bool(complete)
    cursor=generic
    deep=False
    for call,gain in sorted(bycall.items()):
        event={'id':row['id'],'call':call,'rank_from_earlier_calls':cursor,
               'directions':gain,'rank_after_call_attribution':cursor+gain,
               'deep_event':cursor>=23,'attribution':'ORIGINATING_CALL; may be certified by later reconciliation'}
        if times and call in times: event.update(times[call])
        EVENTS.append(event);deep |= cursor>=23;cursor+=gain
    calls=sorted(bycall)
    row['observed_acquiring_calls']=len(calls)
    row['first_gain_call']=calls[0] if calls else None
    row['last_gain_call']=calls[-1] if calls else None
    row['post_seed_exposure']=(row['calls']-calls[0]) if calls else None
    row['calls_since_last_gain']=(row['calls']-calls[-1]) if calls else row['calls']
    if not complete:
        row.update(further_after_seed=None,deep_amplifier=None,classification='UNKNOWN_TRAJECTORY')
    else:
        row['further_after_seed']=int(len(calls)>1) if calls and row['post_seed_exposure']>0 else None
        row['deep_amplifier']=int(deep)
        row['classification']=('deep amplifier' if deep else 'amplifier' if len(calls)>1
                               else 'seed only' if calls else 'no seed')

def broad():
    folder=LOCAL/'broad-rank-v1'; rt=folder/'runtime/research'
    protocol=read(folder/'plan.json'); manifest=read(folder/'manifest.json')
    assert INPUTS[rel(folder/'plan.json')]==manifest['plan_sha256']
    review=read(folder/'COMPLETION_REVIEW.json'); queue=read(folder/'queue.json',review['queue_sha256'])
    for p in protocol['parents']:
        data=read(rt/p['path'],p['sha256'])
        lattice(p['id'],data['generic_height_gram'],rt/p['path'],
                'X948' if p['backend']=='native' else 'X1092')
        pool=rt/'broad-pools'/p['id']; selection=read(pool/'selection.json')
        receipt=read(pool/'result.json')
        assert INPUTS[rel(pool/'selection.json')]==receipt['selection_sha256']
        tables=read(pool/'tables.json',receipt['tables_sha256'])
        # The selected fields are already copied verbatim to the frozen queue.
        selected={z['id']:z for z in selection['rows']}
        for q in queue['rows']:
            if q['parent_id']==p['id']:
                assert all(q[k]==v for k,v in selected[q['id']].items())
                t=F(q['parameter']);primes=[];traces=[];smooth=[];units=[]
                for tab in tables['tables']:
                    prime=tab['prime'];den=t.denominator%prime
                    at=(t.numerator*pow(den,-1,prime))%prime if den else prime
                    primes.append(prime);traces.append(tab['traces'][at]);smooth.append(tab['smooth'][at]);units.append(tab['score_units'][at])
                assert sum(units)==q['score_units']
                PRIME_FEATURES.append({'id':'broad-rank-v1:'+q['id'],'primes':primes,'traces':traces,
                    'smooth':smooth,'score_units':units,'source':rel(pool/'tables.json'),
                    'source_sha256':receipt['tables_sha256'],'scope':'Lookup of frozen residues, no new finite-field point count'})
    rows=[]
    for i,q in enumerate(queue['rows']):
        case=rt/'broad-cases'/q['id']
        state=source(case/'state.json')
        row=base_row('broad-rank-v1',q['id'],q['parent_id'],q['parameter'],q['model'],17,
                     state.get('rank') if state else None)
        row.update(all_frozen_features={k:v for k,v in q.items() if k not in ['model']},
                   arm=q['arm'],detector=q['backend'],source=rel(case/'state.json'),
                   adaptive_followup='RANK_DEPENDENT_CONTINUATION_AND_HASH_RESCUE')
        if not state:
            row['endpoint_status']=q['intake_status'];row['predictive_eligible']=False
            attach_events(row,[],False);rows.append(row);continue
        row.update(calls=state['calls'],endpoint_status=state['status'],
                   censored=state['status'].startswith(('CENSORED','UNKNOWN')),
                   continuation_batches=state['round'],rescue_used=state['rescue_used'])
        events=[];times={};offset=0;cpu=0.;wall=0.;postseedphase=False
        for batch in sorted(case.glob('batch-*')):
            seal=read(batch/'seal.json'); assert seal['completed']
            result=read(batch/'result.json',seal['files']['result.json'])
            proofname='packet-verified.json' if q['backend']=='native' else 'verified.json'
            proof=read(batch/proofname,seal['files'][proofname]);assert proof['status']=='PASS_TWO_FINITE_IMPLEMENTATIONS'
            # Bind the packets, without rerunning elliptic arithmetic.
            packet=read(batch/'packet.json',result['packet_sha256'])
            assert proof['packet_sha256']==result['packet_sha256']==seal['files']['packet.json']
            assert len(packet['points'])==packet['rank_lower_bound']==result['rank_lower_bound']
            assert list(map(F,packet['curve']))==list(map(F,q['model']))
            if offset==0:
                row['initial_rank']=result['rank_lower_bound'];row['initial_calls']=result['calls']
                row['initial_seed_calls']=result.get('seed_calls')
            for e in result.get('gain_timeline',[]):
                call=offset+e['call'];events.append({**e,'call':call})
                times[call]={'batch_cpu_lower':cpu,'batch_cpu_upper':cpu+seal['cpu_seconds'],
                             'batch_elapsed_lower':wall,'batch_elapsed_upper':wall+seal['wall_seconds'],
                             'timing_scope':'Enclosing batch sums; includes preparation and verification. Not exact discovery CPU.'}
                postseedphase |= e.get('phase') not in ('seed','seed-cloud')
            offset+=result['calls'];cpu+=seal['cpu_seconds'];wall+=seal['wall_seconds']
        assert offset==state['calls'] and events==state['history']
        row.update(cpu_seconds=cpu,elapsed_seconds=wall,
                   timing_scope='SUM_OF_BATCH_PROCESS_TREE_CPU_AND_ELAPSED; excludes queue waiting and shared scoring',
                   gained_after_seed_phase=bool(postseedphase))
        attach_events(row,events,True,times);rows.append(row)
        if i%400==0:print(f'broad rows {i}/{len(queue["rows"])}',flush=True)
    assert Counter(str(x['final_rank']) for x in rows)==Counter(review['final_lower_bound_histogram'])
    return rows

def r17_panel():
    path=ART/'r17_60_panel_results_v1.json';d=read(path)
    roster={x['id']:x for x in read(ART/'r17_60_panel_roster_v1.json')['rows']}
    rows=[]
    for r in d['results']:
        s=roster[r['id']];p=r['packet']
        row=base_row('r17-60-panel-v1',r['id'],r['family'],r['parameter'],s['model'],17,r['rank_lower_bound'])
        assert r['status']=='PASS_INDEPENDENT_R17_60_CASE'
        row.update(calls=r['total_calls'],source=rel(path),all_frozen_features=s,
            censored=bool(r['point_timeouts'] or r['map_timeouts']),
            endpoint_status=r['stop_reason'],initial_rank=17+sum(e['after']-e['before'] for e in r['gain_timeline'] if e['phase'].startswith('seed')),
            adaptive_followup='ALL_SEEDED_FIBRES_AT_MOST_100_COMPLEMENT_CALLS',
            gained_after_seed_phase=any(e['phase']=='complement' for e in r['gain_timeline']))
        events=[{**e,'call':e['call']+(r['seed_calls'] if e['phase']=='complement' else 0)} for e in r['gain_timeline']]
        attach_events(row,events,True);rows.append(row)
    return rows

def protocol_features(folder, family, parameter, cid):
    p=folder/'protocol.json'
    if not p.is_file():
        # Early workers put results below a family subdirectory.
        if (folder.parent/'protocol.json').is_file():
            proto=small(folder.parent/'protocol.json')
            pop=folder/'population.json'
            records=small(pop).get('retained_candidates',[]) if pop.is_file() else proto.get('fixed_addresses',[])
            for r in records:
                if str(r.get('parameter'))==str(parameter):return r,proto
        return {},None
    d=small(p)
    records=d.get('rows',d.get('candidates',[]))
    if isinstance(records,list):
        for r in records:
            if isinstance(r,dict) and str(r.get('parameter'))==str(parameter) and r.get('family',family)==family:
                return r,d
    return {},d

def legacy():
    rows=[]
    # MW16 Gram is separate from its rank-17 integral frame.
    p=ART/'compact_five_mw16_atlas_v1.json'
    for f in read(p)['families']:
        lattice(f['fibration_id'],f['generic_height_gram'],p,'X948')
    for cohort in LEGACY:
        path=ART/(cohort+'_results_v1.json');d=source(path)
        if d is None: continue
        for i,r in enumerate(d['curves']):
            n=len(r['generic_points']);cid=r.get('id',str(i));family=r['family']
            row=base_row(cohort,cid,family,r['parameter'],r['curve'],n,r['rank_lower_bound'])
            row.update(source=rel(path),known_equation=bool(r.get('icarm_matches')),
                previous_equation=bool(r.get('previous_matches')),initial_rank=r['rank_lower_bound'],
                endpoint_status=r.get('search_status','UNKNOWN'),
                calls=r.get('attempted_charts',r.get('completed_charts')),
                adaptive_followup='FIXED_GENERIC_CHART_BANK; later targeted follow-ups not pooled here',
                all_frozen_features={},censored=True)
            witness=r.get('discovery_witness')
            if witness is None and cohort in ['strata60_mw16','nearcut60v2_mw16']:
                directory='strata60-mw16-pari-v1' if cohort=='strata60_mw16' else 'nearcut60-mw16-pari-v2'
                ledger=small(LOCAL/directory/'ledger.json')
                lr=next(x for x in ledger['rows'] if x['id']==cid)
                witness={'path':rel(LOCAL/directory/cid/'result.json'),'sha256':lr['result_sha256']}
            if witness is None:
                GAPS.append({'id':row['id'],'status':'NO_DISCOVERY_WITNESS'})
                attach_events(row,[],False);rows.append(row);continue
            wp=ROOT/witness['path'];w=source(wp,witness['sha256'])
            features,proto=protocol_features(wp.parent.parent,family,r['parameter'],cid)
            row['all_frozen_features']=features
            if proto:
                row['search_height']=proto.get('height')
                row['seconds_per_chart']=proto.get('seconds_per_chart')
                row['selection_protocol']=rel(next(p for p in
                    [wp.parent.parent/'protocol.json',wp.parent.parent.parent/'protocol.json'] if p.is_file()))
            if w is None:
                attach_events(row,[],False);rows.append(row);continue
            assert list(map(F,w['curve']))==list(map(F,row['model'])),row['id']
            charts=w.get('charts',[])
            if not isinstance(charts,list):charts=[]
            row['calls']=len(charts)
            prev=n;events=[];times={};wall=0.;cpu=0.;bounded=True;cpu_known=True
            for j,c in enumerate(charts,1):
                search=c.get('search',{})
                if isinstance(search,dict):
                    wall+=float(search.get('wall_seconds',0) or 0)
                    ms=search.get('search_cpu_ms')
                    if ms is None:cpu_known=False
                    else:cpu+=float(ms)/1000
                    bounded &= search.get('status')=='bounded_search_complete'
                rank=c.get('rank_lower_bound')
                if rank is not None and rank>prev:
                    events.append({'before':prev,'after':rank,'call':j});prev=rank
                    times[j]={'point_search_elapsed_sum':wall,
                        'point_search_cpu_sum':cpu if cpu_known else None,
                        'timing_scope':'Point-box sums only; excludes preparation/admission/verification'}
            row['point_search_elapsed_seconds']=wall
            row['point_search_cpu_seconds']=cpu if cpu_known else None
            row['censored']=not bounded or (r.get('declared_charts') is not None and len(charts)<r['declared_charts'])
            row['trajectory_rank']=prev
            row['trajectory_gap']=r['rank_lower_bound']-prev
            row['rank_evidence']='RETAINED_FULL_COHORT_CERTIFICATE; discovery hash checked; no arithmetic rerun'
            if proto and isinstance(proto.get('rows'),list):
                row['declared_cohort_size']=len(proto['rows'])
            attach_events(row,events,prev==r['rank_lower_bound'],times);rows.append(row)
        print(f'legacy {cohort}: {len(d["curves"])} rows',flush=True)
    return rows

def parent_foundry():
    folder=LOCAL/'parent-foundry-v4';path=folder/'state.json'
    if not path.is_file():return []
    state=read(path);roots=[folder/'runtime/research']
    ancestor=state.get('continuation_from')
    while ancestor:
        ancestor=Path(ancestor);prior=read(ancestor/'state.json');roots.append(ancestor/'runtime/research')
        ancestor=prior.get('continuation_from')
    def find(name):return next((rt/name for rt in roots if (rt/name).is_file()),None)
    parents={}
    for name,p in state['parents'].items():
        pp=find(p['path'])
        if pp:
            data=read(pp,p['sha256']);parents[name]=data
            if name not in LATTICES and 'generic_height_gram' in data:
                lattice(name,data['generic_height_gram'],pp,p['source_surface'].split(';')[0])
    jobs=defaultdict(list)
    for name,j in sorted(state['jobs'].items()):
        if j.get('fibre'):jobs[j['fibre']].append(j)
    rows=[]
    for cid,f in state['fibres'].items():
        packets=find(f['packet']) if f.get('packet') else None
        packet=read(packets,f['packet_sha256']) if packets else None
        model=packet['curve'] if packet else None
        row=base_row('parent-foundry-v1-v4',cid,f['family'],f['parameter'],model,
            f['generic_rank'] if f.get('rank') is not None else None,f.get('rank'))
        row.update(calls=f['calls'],source=rel(path),censored=True,
            endpoint_status='STOPPED_ADAPTIVE_PARENT_FOUNDRY',predictive_eligible=False,
            adaptive_followup='OUTCOME_DEPENDENT_PANEL_AND_EXPLOIT; administrative stop',
            all_frozen_features={'panel_index':f['panel_index']},
            elapsed_seconds=sum(j.get('finished_at',0)-j.get('started_at',0) for j in jobs[cid]),
            timing_scope='Summed retained job elapsed, not CPU; no timing for unretained attempts')
        events=[];offset=0;missing=False
        for job in jobs[cid]:
            rp=find(job['path']+'/result.json')
            if rp is None:missing=True;continue
            r=read(rp)
            if r.get('status')!='PASS_CERTIFIED_PARENT_EVALUATION':continue
            if row['initial_rank'] is None:row['initial_rank']=r['rank_lower_bound']
            events.extend({**e,'call':offset+e['call']} for e in r.get('gain_timeline',[]))
            offset+=r['calls']
        complete=not missing and offset==f['calls'] and f.get('rank') is not None
        attach_events(row,events,complete);rows.append(row)
    return rows

def record_scale():
    path=ART/'det1092_record_scale_points_v1.json';d=read(path)
    assert d['status']=='PASS' and d['certified_gains']==0 and d['certified_boxes']==2352
    sp=ART/'det1092_record_scale_intake_v1/selection-result.json.gz'
    raw=sp.read_bytes();INPUTS[rel(sp)]=sha(raw)
    selection=json.loads(gzip.decompress(raw));selected={r['id']:r for r in selection['selected']}
    rows=[]
    for r in d['rows']:
        s=selected[r['id']]
        assert r['status']=='PASS' and r['certified_rank_lower_bound']==17 and r['certified_boxes']==49
        row=base_row('det1092-record-scale-v1',r['id'],'det1092',r['parameter'],s['model'],17,17)
        row.update(calls=49,initial_rank=17,source=rel(path),censored=False,
            endpoint_status='COMPLETED_FIXED_EXPOSURE',all_frozen_features=s,
            elapsed_seconds=r['completed_stage_seconds'],timing_scope='Summed finished stage elapsed; not CPU',
            adaptive_followup='FOLLOWUP_CONDITIONAL_ON_GAIN; no seed observed')
        attach_events(row,[],True);rows.append(row)
    return rows

def foundry():
    path=LOCAL/'high-rank-foundry-v3/ledger.sqlite'
    if not path.is_file():return []
    before=sha(path.read_bytes())
    con=sqlite3.connect(f'file:{path}?mode=ro',uri=True)
    records=[json.loads(r[0]) for r in con.execute('select record from curves order by id')]
    jobs=defaultdict(list)
    for r, in con.execute('select record from jobs order by id'):
        j=json.loads(r);jobs[j['cid']].append(j)
    con.close();assert sha(path.read_bytes())==before
    INPUTS[rel(path)]=before
    rows=[]
    for r in records:
        rr=r.get('rank');row=base_row('high-rank-foundry-v1-v3',r['id'],r['family'],r['parameter'],r['model'],17 if rr is not None else None,rr)
        row.update(calls=r['total_calls'],all_frozen_features={k:r.get(k) for k in
            ['combined_good','combined_selection_units','parameter_height','selection_ordinal','lane','numerator','denominator']},
            source=rel(path),endpoint_status=r['state'],
            censored=r['state'] not in ['EXHAUSTED'],
            adaptive_followup='STOPPED_ADAPTIVE_FOUNDRY; outcome-dependent exposure',
            known_equation=bool(r.get('novelty',{}).get('matches')),
            elapsed_seconds=sum(j.get('charged_seconds',0) for j in jobs[r['id']]),
            timing_scope='Charged worker elapsed; CPU unavailable, not worker-hours')
        done=[j for j in jobs[r['id']] if j.get('rank_after') is not None]
        row['initial_rank']=done[0]['rank_after'] if done else None
        row['job_count']=len(jobs[r['id']])
        row['predictive_eligible']=False # administrative stop and adaptive policy, descriptive only
        attach_events(row,r.get('history',[]),rr is not None);rows.append(row)
    return rows

def csvfile(path, rows):
    keys=sorted({k for r in rows for k in r})
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader()
        for r in rows:w.writerow({k:json.dumps(v,sort_keys=True) if isinstance(v,(dict,list)) else v for k,v in r.items()})

def extract():
    start=time.monotonic();read(OUT/'protocol.json')
    rows=broad()+r17_panel()+legacy()+foundry()+parent_foundry()+record_scale()
    assert len({r['id'] for r in rows})==len(rows)
    # Hash j clusters are conservative fold exclusions; do not merge twists.
    counts=Counter(r.get('j_leakage_group') for r in rows)
    for r in rows:
        r['same_j_rows']=counts[r.get('j_leakage_group')]
        # A blind rediscovery is still part of the selected population. These
        # flags limit novelty, not outcome-based inclusion in predictive tables.
    OUT.mkdir(exist_ok=True,parents=True)
    (OUT/'fibres.jsonl').write_text(''.join(json.dumps(r,sort_keys=True,allow_nan=False)+'\n' for r in rows))
    csvfile(OUT/'fibres.csv',rows);csvfile(OUT/'gain_events.csv',EVENTS)
    (OUT/'prime_features.jsonl').write_text(''.join(json.dumps(r,sort_keys=True)+'\n' for r in PRIME_FEATURES))
    write(OUT/'lattices.json',LATTICES)
    write(OUT/'extraction.json',{'rows':len(rows),'events':len(EVENTS),'cohorts':dict(Counter(r['cohort'] for r in rows)),
        'inputs':INPUTS,'gaps':GAPS,'elapsed_seconds':time.monotonic()-start,
        'coverage_boundary':'Full denominators for listed cohorts. Early success-only exports and non-dispatched proposals do not supply complete prospective denominators; see coverage inventory. No points or missing artifacts are reconstructed.',
        'verification':'Checks joins, row accounting, retained hashes and batch packet/replay bindings; does not rerun elliptic arithmetic.',
        'outputs':{p.name:sha(p.read_bytes()) for p in [OUT/'fibres.jsonl',OUT/'fibres.csv',OUT/'gain_events.csv',OUT/'lattices.json',OUT/'prime_features.jsonl']}})
    print(json.dumps({'rows':len(rows),'events':len(EVENTS),'gaps':len(GAPS)}))
    enrich_features()

def enrich_features():
    """Join only selection files matching the frozen point-protocol hash."""
    manifest=read(OUT/'extraction.json')
    rows=[json.loads(l) for l in (OUT/'fibres.jsonl').read_text().splitlines()]
    wanted={}
    for row in rows:
        name=row.get('selection_protocol')
        if not name:continue
        p=ROOT/name
        if not p.is_file() and (p.parent.parent/'protocol.json').is_file():
            p=p.parent.parent/'protocol.json';row['selection_protocol']=rel(p)
        if p.is_file():
            proto=small(p)
            if proto.get('selection_sha256'):wanted.setdefault(proto['selection_sha256'],[]).append(row)
    candidates=sorted(set(LOCAL.glob('*select*/result.json'))|set(ART.glob('*selection*.json')))
    matched=[]
    for path in candidates:
        if path.stat().st_size>=30_000_000:continue
        b=path.read_bytes();h=sha(b)
        if h not in wanted:continue
        INPUTS[rel(path)]=h;d=json.loads(b)
        available=d.get('rows',d.get('selected',d.get('candidates',[])))
        if isinstance(available,dict):available=list(available.values())
        index={(r.get('family'),str(F(r['parameter']))):r for r in available
               if isinstance(r,dict) and r.get('parameter') is not None}
        count=0
        for row in wanted[h]:
            r=index.get((row['family'],row['parameter']))
            if r is None:continue
            prior=row['all_frozen_features']
            # Some older point protocols rename a later-cutoff score to the
            # same field used by an earlier selection stage. Preserve stages
            # separately; never replace the declared point-roster score.
            row['all_frozen_features']={**prior,'source_selection_record':r}
            row['full_selection_source']=rel(path);row['full_selection_sha256']=h;count+=1
        matched.append({'path':rel(path),'sha256':h,'rows_joined':count})
    (OUT/'fibres.jsonl').write_text(''.join(json.dumps(r,sort_keys=True,allow_nan=False)+'\n' for r in rows))
    csvfile(OUT/'fibres.csv',rows)
    for name in ['fibres.jsonl','fibres.csv']:manifest['outputs'][name]=sha((OUT/name).read_bytes())
    manifest['inputs'].update({k:v for k,v in INPUTS.items() if not k.startswith(rel(OUT)+'/')})
    manifest['selection_enrichment']=matched
    manifest['feature_boundary']='All exposed frozen row fields retained, plus hash-matched full selections and broad prime vectors. Further unmaterialized scoring pools stay available at their recorded sources; no scores are recomputed.'
    write(OUT/'extraction.json',manifest)
    print(json.dumps({'selection_rows_enriched':sum(r['rows_joined'] for r in matched)}))

FEATURES=['log2_height','log2_denominator','negative_parameter','model_bits',
          'discriminant_bits','score_97','score_257','score_997','local_bad_count']
OLDER_FEATURES=['log2_height','log2_denominator','negative_parameter','model_bits',
                'discriminant_bits','selection_score']

def feature(row, name):
    frozen=row['all_frozen_features']
    if name=='log2_height':return math.log2(row['parameter_height'])
    if name=='log2_denominator':return math.log2(row['denominator'])
    if name=='negative_parameter':return int(row['numerator']<0)
    if name=='model_bits':return frozen.get('model_bits',row.get('displayed_model_bits'))
    if name=='discriminant_bits':return frozen.get('discriminant_bits',row.get('displayed_discriminant_bits'))
    if name.startswith('score_'):
        v=frozen.get('scores_by_cutoff',{}).get(name[6:]);return v/1e12 if v is not None else None
    if name=='selection_score':
        v=frozen.get('combined_selection_units',frozen.get('selection_units',frozen.get('combined_late_units',frozen.get('score_units'))))
        return v/1e12 if v is not None else None
    if name=='local_bad_count':
        v=frozen.get('local_singular_primes');return len(v) if v is not None else None
    raise KeyError(name)

def train_transform(train, test, names):
    import numpy as np
    raw=np.array([[feature(r,n) if feature(r,n) is not None else np.nan for n in names] for r in train])
    ext=np.array([[feature(r,n) if feature(r,n) is not None else np.nan for n in names] for r in test])
    med=np.array([np.nanmedian(raw[:,j]) if np.isfinite(raw[:,j]).any() else 0 for j in range(len(names))])
    raw=np.where(np.isfinite(raw),raw,med);ext=np.where(np.isfinite(ext),ext,med)
    mean=raw.mean(axis=0);scale=raw.std(axis=0);scale[scale<1e-12]=1
    return (raw-mean)/scale,(ext-mean)/scale,mean,scale

def logistic(x,y,z):
    import numpy as np
    from scipy.optimize import minimize
    from scipy.special import expit
    xx=np.column_stack((np.ones(len(x)),x));zz=np.column_stack((np.ones(len(z)),z))
    prev=(sum(y)+.5)/(len(y)+1)
    if min(y)==max(y):return np.full(len(z),prev),{'constant':prev}
    def fun(w):
        eta=xx@w
        return float(np.sum(np.logaddexp(0,eta)-y*eta)+.5*np.dot(w[1:],w[1:]))
    def jac(w):
        gradient=xx.T@(expit(xx@w)-y);gradient[1:]+=w[1:];return gradient
    init=np.zeros(xx.shape[1]);init[0]=math.log(prev/(1-prev))
    fit=minimize(fun,init,jac=jac,method='L-BFGS-B',options={'maxiter':1000,'ftol':1e-12,'gtol':1e-7})
    assert fit.success, fit.message
    return expit(zz@fit.x),{'standardized_coefficients':fit.x.tolist(),'gradient_max':float(max(abs(jac(fit.x))))}

def tree(x,y,z,names,mean,scale):
    import numpy as np
    minleaf=max(10,math.ceil(.05*len(x)))
    def fit(indices,depth):
        yes=int(y[indices].sum());n=len(indices);leaf={'n':n,'yes':yes,'p':(yes+.5)/(n+1)}
        if depth==2 or n<2*minleaf or yes in (0,n):return leaf
        best=None
        for j in range(x.shape[1]):
            for threshold in sorted(set(float(a) for a in np.quantile(x[indices,j],np.arange(.1,1,.1)))):
                left=indices[x[indices,j]<=threshold];right=indices[x[indices,j]>threshold]
                if min(len(left),len(right))<minleaf:continue
                impurity=sum(float(y[a].sum())*(len(a)-float(y[a].sum()))/len(a) for a in [left,right])
                candidate=(impurity,j,threshold)
                if best is None or candidate<best[0]:best=(candidate,left,right)
        if best is None or best[0][0]>=yes*(n-yes)/n-1e-12:return leaf
        (_,j,threshold),left,right=best
        return {**leaf,'feature_index':j,'feature':names[j],'threshold_standardized':threshold,
            'threshold_original':threshold*scale[j]+mean[j],
            'left':fit(left,depth+1),'right':fit(right,depth+1)}
    model=fit(np.arange(len(x)),0)
    def prediction(row,node):
        if 'left' not in node:return node['p']
        return prediction(row,node['left'] if row[node['feature_index']]<=node['threshold_standardized'] else node['right'])
    return np.array([prediction(r,model) for r in z]),model

def metrics(y,p):
    import numpy as np
    from scipy.stats import rankdata
    y=np.asarray(y);p=np.clip(np.asarray(p),1e-12,1-1e-12)
    auc=None
    if 0<sum(y)<len(y):auc=float((rankdata(p)[y==1].sum()-sum(y)*(sum(y)+1)/2)/(sum(y)*(len(y)-sum(y))))
    return {'n':len(y),'yes':int(sum(y)),'brier':float(np.mean((p-y)**2)),
        'log_loss':float(-np.mean(y*np.log(p)+(1-y)*np.log1p(-p))),'auc':auc}

def contingency(a,b,c,d):
    from scipy.stats import fisher_exact
    rr=(a/(a+b))/(c/(c+d)) if a+b and c+d and c else None
    return {'table':[[a,b],[c,d]],'risk_ratio':rr,
            'fisher_one_sided':float(fisher_exact([[a,b],[c,d]],alternative='greater').pvalue)}

def descriptive(rows):
    selected=[r for r in rows if r.get('final_rank') is not None]
    seeds=[r for r in selected if r['seed']]
    followed=[r for r in seeds if r.get('further_after_seed') is not None]
    deep=[r for r in selected if r.get('deep_amplifier') is not None]
    return {'rows':len(rows),'certified_endpoints':len(selected),
        'seed':[sum(r['seed'] for r in selected),len(selected)],
        'later_acquiring_call_given_seed_and_exposure':[sum(r['further_after_seed'] for r in followed),len(followed)],
        'any_second_direction_given_seed':[sum(r['extra_direction_after_first'] for r in seeds),len(seeds)],
        'gain_after_seed_phase':sum(bool(r.get('gained_after_seed_phase')) for r in seeds),
        'deep':[sum(r['deep_amplifier'] for r in deep),len(deep)],
        'trajectories_unknown':sum(not r['trajectory_complete'] for r in rows),
        'censored_or_stopped':sum(r['censored'] is not False for r in rows),
        'classes':dict(Counter(r['classification'] for r in rows)),
        'ranks':dict(sorted(Counter(str(r['final_rank']) for r in selected).items())),
        'known_equations':sum(r['known_equation'] for r in rows),
        'previous_equations':sum(r['previous_equation'] for r in rows),
        'cpu_seconds':sum(r['cpu_seconds'] or 0 for r in rows) if any(r['cpu_seconds'] is not None for r in rows) else None}

def validate(rows, label, names, groupkey='family', external_train=None):
    import numpy as np
    folds=[];predictions=[]
    for target in ['seed','further_after_seed','deep_amplifier']:
        data=[r for r in rows if r.get(target) is not None and r['predictive_eligible'] and r['censored'] is False]
        training_data=data if external_train is None else [r for r in external_train
            if r.get(target) is not None and r['predictive_eligible'] and r['censored'] is False]
        groups=sorted({r[groupkey] for r in data})
        if len(groups)<2:continue
        for holdout in groups:
            test=[r for r in data if r[groupkey]==holdout]
            forbidden={r['j_leakage_group'] for r in test}
            train=[r for r in training_data if r[groupkey]!=holdout and r['j_leakage_group'] not in forbidden]
            if not train or not test:continue
            x,z,mean,scale=train_transform(train,test,names)
            y=np.array([r[target] for r in train],dtype=float);ty=[r[target] for r in test]
            prevalence=(sum(y)+.5)/(len(y)+1)
            base=metrics(ty,np.full(len(test),prevalence))
            for modelname in ['logistic','tree_depth2']:
                probs,model=logistic(x,y,z) if modelname=='logistic' else tree(x,y,z,names,mean,scale)
                nsel=math.ceil(len(test)/4)
                order=sorted(range(len(test)),key=lambda i:(-probs[i],test[i]['id']))
                chosen=set(order[:nsel]);a=sum(ty[i] for i in chosen);c=sum(ty)-a
                table=contingency(a,nsel-a,c,len(test)-nsel-c)
                item={'cohort':label,'target':target,'group_key':groupkey,'holdout':holdout,
                    'train_rows':len(train),'train_positives':int(sum(y)),'test':metrics(ty,probs),
                    'baseline':base,'model':modelname,'features':names,'fit':model,
                    'top_quarter':table,'selected':nsel,'training_only_mean':mean.tolist(),
                    'training_only_scale':scale.tolist()}
                folds.append(item)
                for i,r in enumerate(test):predictions.append({'id':r['id'],'cohort':label,'holdout':holdout,
                    'target':target,'model':modelname,'y':ty[i],'p':float(probs[i]),
                    'baseline_p':float(prevalence),'selected':i in chosen})
    # One correction over every reported model/fold test in this analysis set.
    order=sorted(range(len(folds)),key=lambda i:folds[i]['top_quarter']['fisher_one_sided'])
    previous=0.
    for j,i in enumerate(order):
        previous=max(previous,min(1.,(len(folds)-j)*folds[i]['top_quarter']['fisher_one_sided']))
        folds[i]['top_quarter']['holm_adjusted']=previous
    pooled=[]
    for target in ['seed','further_after_seed','deep_amplifier']:
        for modelname in ['logistic','tree_depth2']:
            p=[r for r in predictions if r['target']==target and r['model']==modelname]
            if not p:continue
            chosen=[r for r in p if r['selected']];rest=[r for r in p if not r['selected']]
            pooled.append({'cohort':label,'target':target,'model':modelname,
                'heldout':metrics([r['y'] for r in p],[r['p'] for r in p]),
                'baseline':metrics([r['y'] for r in p],[r['baseline_p'] for r in p]),
                'top_quarter_counts':[[sum(r['y'] for r in chosen),len(chosen)], [sum(r['y'] for r in rest),len(rest)]],
                'boundary':'Concatenated held-out predictions; folds share training data, no pooled significance test.'})
    return folds,predictions,pooled

def analyze():
    import numpy as np
    import scipy
    start=time.monotonic();extraction=read(OUT/'extraction.json')
    for name,h in extraction['outputs'].items():assert sha((OUT/name).read_bytes())==h
    rows=[json.loads(line) for line in (OUT/'fibres.jsonl').read_text().splitlines()]
    groups=[]
    for cohort in sorted({r['cohort'] for r in rows}):
        cohortrows=[r for r in rows if r['cohort']==cohort]
        groups.append({'cohort':cohort,'family':'ALL',**descriptive(cohortrows)})
        for family in sorted({r['family'] for r in cohortrows}):
            sub=[r for r in cohortrows if r['family']==family]
            groups.append({'cohort':cohort,'family':family,**descriptive(sub)})
            if cohort=='broad-rank-v1':
                for arm in ['ranked','control']:
                    groups.append({'cohort':cohort,'family':family,'arm':arm,
                        **descriptive([r for r in sub if r['arm']==arm])})
    studies=[('broad-native',[r for r in rows if r['cohort']=='broad-rank-v1' and r['surface']=='X948'],FEATURES,'family'),
        ('broad-all-parent-detector',[r for r in rows if r['cohort']=='broad-rank-v1'],FEATURES,'family'),
        ('broad-surface-diagnostic',[r for r in rows if r['cohort']=='broad-rank-v1'],FEATURES,'surface'),
        ('r17-60-panel',[r for r in rows if r['cohort']=='r17-60-panel-v1'],OLDER_FEATURES,'family')]
    for arm in ['ranked','control']:
        studies.append(('broad-native-'+arm,[r for r in rows if r['cohort']=='broad-rank-v1'
            and r['surface']=='X948' and r['arm']==arm],FEATURES,'family'))
    for cohort in LEGACY:
        r=[x for x in rows if x['cohort']==cohort]
        if len(r)>=20 and len({x['family'] for x in r})>=3:
            studies.append((cohort,r,OLDER_FEATURES,'family'))
    folds=[];predictions=[];pooled=[]
    for label,data,names,groupkey in studies:
        f,p,s=validate(data,label,names,groupkey);folds+=f;predictions+=p;pooled+=s
        print(f'validated {label}: {len(f)} model/fold records',flush=True)
    # One forward cohort-transfer check: earlier compact192 training, later
    # r17-60 test. The test fibration is absent from BOTH training populations.
    # Identical models and scoring cutoff; no rule/threshold tuning on this test.
    f,p,s=validate([r for r in rows if r['cohort']=='r17-60-panel-v1'],
        'compact192-to-r17-60',OLDER_FEATURES,external_train=[r for r in rows if r['cohort']=='compact192_r17'])
    folds+=f;predictions+=p;pooled+=s
    order=sorted(range(len(folds)),key=lambda i:folds[i]['top_quarter']['fisher_one_sided'])
    previous=0.
    for j,i in enumerate(order):
        previous=max(previous,min(1.,(len(folds)-j)*folds[i]['top_quarter']['fisher_one_sided']))
        folds[i]['top_quarter']['holm_all_reported_tests']=previous
    csvfile(OUT/'heldout_predictions.csv',predictions)
    csvfile(OUT/'cohort_summary.csv',groups)
    write(OUT/'heldout_folds.json',folds)
    write(OUT/'summary.json',{'schema':'seed-amplification-history.analysis.v1',
        'protocol_sha256':sha((OUT/'protocol.json').read_bytes()),
        'extraction_sha256':sha((OUT/'extraction.json').read_bytes()),
        'source_sha256':sha(Path(__file__).read_bytes()),'numpy':np.__version__,'scipy':scipy.__version__,
        'elapsed_seconds':time.monotonic()-start,'tables':groups,'heldout_summaries':pooled,
        'analysis_boundary':'Historical fitted associations under recorded search policies. No new independent points, arithmetic theorem or general success probability is certified.'})
    print(json.dumps({'rows':len(rows),'folds':len(folds),'predictions':len(predictions),'elapsed_seconds':time.monotonic()-start}))

def audit():
    """A second, CSV-based accounting path; no rank-certificate recalculation."""
    rows=[json.loads(l) for l in (OUT/'fibres.jsonl').read_text().splitlines()]
    index={r['id']:r for r in rows};byrow=defaultdict(list)
    with (OUT/'gain_events.csv').open() as f:events=list(csv.DictReader(f))
    for e in events:byrow[e['id']].append(e)
    checked=0
    for row in rows:
        if not row['trajectory_complete']:continue
        rank=row['generic_rank'];last=-1;deep=False
        for e in byrow[row['id']]:
            call=int(e['call']);assert last<call<=row['calls'];last=call
            assert int(e['rank_from_earlier_calls'])==rank
            deep |= rank>=23;rank+=int(e['directions'])
        assert rank==row['final_rank'] and deep==bool(row['deep_amplifier'])
        assert bool(byrow[row['id']])==bool(row['seed']);checked+=1
    vectors=[json.loads(l) for l in (OUT/'prime_features.jsonl').read_text().splitlines()]
    for vector in vectors:assert sum(vector['score_units'])==index[vector['id']]['all_frozen_features']['score_units']
    summary=read(OUT/'summary.json')
    assert summary['source_sha256']==sha(Path(__file__).read_bytes())
    for group in summary['tables']:
        selected=[r for r in rows if r['cohort']==group['cohort'] and
            (group['family']=='ALL' or r['family']==group['family']) and
            ('arm' not in group or r.get('arm')==group['arm'])]
        assert len(selected)==group['rows']
        assert sum(r.get('seed') or 0 for r in selected)==group['seed'][0]
        assert sum(r.get('deep_amplifier') or 0 for r in selected)==group['deep'][0]
    import numpy as np
    x=np.arange(-20,20,dtype=float)[:,None];y=(x[:,0]>=0).astype(float)
    lp,_=logistic(x,y,x);assert all(np.diff(lp)>0)
    tp,fit=tree(x,y,x,['x'],np.array([0]),np.array([1]))
    assert tp[0]<.1 and tp[-1]>.9 and fit['feature']=='x'
    assert metrics([0,1],[.1,.9])['auc']==1
    names=['protocol.json','audit_corrections.json','extraction.json','fibres.jsonl','fibres.csv',
           'gain_events.csv','prime_features.jsonl','summary.json','heldout_folds.json',
           'heldout_predictions.csv','cohort_summary.csv','lattices.json']
    write(OUT/'accounting_audit.json',{'status':'PASS_RETAINED_DATA_ACCOUNTING','rows':len(rows),
        'complete_trajectories_checked':checked,'acquiring_events':len(events),
        'prime_vectors_checked':len(vectors),'summary_groups_checked':len(summary['tables']),
        'scope':'Separate CSV/JSON aggregation verifies event conservation, same-call merging, deep-event definition, score sums and exact counts. Tiny analytic model sanity checks pass. No point search or elliptic arithmetic replay; no statistical success theorem.',
        'bindings':{name:sha((OUT/name).read_bytes()) for name in names}})
    print(json.dumps({'status':'PASS_RETAINED_DATA_ACCOUNTING','rows':len(rows),'trajectories':checked}))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('mode',choices=['extract','enrich','analyze','audit'])
    args=parser.parse_args()
    if args.mode=='extract':extract()
    elif args.mode=='enrich':enrich_features()
    elif args.mode=='analyze':analyze()
    else:audit()
