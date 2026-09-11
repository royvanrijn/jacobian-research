#!/usr/bin/env python3
"""Freeze and replay the completed broad campaign's requested rank >22 cohort."""
import argparse
from collections import Counter
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

import broad_rank_runtime as io
import certify_compact_r17_candidates as independent
from memory_rank_certificate import checked_rank
from parent_foundry_worker import specialize
from research_curve_refresh import ROOT, MANIFEST
from v3_warm_support import require

RUN = ROOT/'artifacts/local/elliptic-curves/broad-rank-v1'
SNAPSHOT = ROOT/'artifacts/generated-results/elliptic-curves/broad_rank_ledger_snapshot_v1.json'
REPLAY = SNAPSHOT.with_name('broad_rank_ledger_replay_v1.json')
CLAIM = 'EC-BROAD-RANK-LEDGER-20260912'


def digest(value):
    return hashlib.sha256((json.dumps(value, sort_keys=True, indent=2)+'\n').encode()).hexdigest()


def invariant(model):
    v = independent.weierstrass_invariants(tuple(map(F, model)))
    require(v['discriminant'] != 0, 'singular model')
    return v['c4']**3/v['discriminant']


def validate_record(row, parent):
    packet, receipt = row['packet'], row['certificate_replay']
    expected = 'PASS_CERTIFIED_SEARCH' if row['backend'] == 'native' else 'PASS_CERTIFIED_PARENT_EVALUATION'
    require(row['source_status'] == expected, 'uncertified source endpoint')
    require(receipt['status'] == 'PASS_TWO_FINITE_IMPLEMENTATIONS', 'missing source replay')
    require(digest(packet) == row['packet_sha256'] == receipt['packet_sha256'], 'packet binding differs')
    require(row['rank_lower_bound'] == packet['rank_lower_bound'] == len(packet['points']) > 22,
            'rank or publication threshold differs')
    require(parent['family'] == row['family'], 'parent family differs')
    model, generic = specialize(parent, row['parameter'])
    require(tuple(map(F, packet['curve'])) == model, 'specialized equation differs')
    require(tuple(tuple(map(F, p)) for p in packet['points'][:len(generic)]) == generic,
            'generic prefix differs')


def freeze():
    require(not SNAPSHOT.exists(), 'preserve snapshot; use another version')
    plan, rt = io.guard(RUN)
    complete = io.read(RUN/'COMPLETE.json')
    require(complete['status'] == 'COMPLETE_BOUNDED_CAMPAIGN'
            and complete['queue_sha256'] == io.sha(RUN/'queue.json'), 'completed queue required')
    database = ROOT/'elliptic-curves/data/research_curves/database.json'
    baseline = [{k:r[k] for k in ('id','ainvs','rank_lower_bound')} for r in io.read(database)['curves']]
    parents = {p['id']:io.read(rt/p['path']) for p in plan['parents']}
    records = []
    for row in io.read(RUN/'queue.json')['rows']:
        if row['intake_status'] != 'ACCEPTED':
            continue
        case = rt/'broad-cases'/row['id']
        state = io.read(case/'state.json')
        if (state.get('rank') or 0) <= 22:
            continue
        require(state['status'] in ('CERTIFIED', 'POLICY_EXHAUSTED'), 'uncertified final state')
        job = case/f'batch-{state["round"]:03d}'
        require(state == io.read(job/'broad-state.json') and io.sha(job/'broad-state.json') ==
                io.read(job/'broad-state-seal.json')['sha256'], 'final state changed')
        seal, result = io.completed(job)
        require(state['source_seal_sha256'] == io.sha(job/'seal.json'), 'state source differs')
        record = {'id':'broad-'+row['parent_id']+'-'+row['id'][2:],
            'source_case_id':row['id'], 'family':row['family'], 'parameter':row['parameter'],
            'backend':row['backend'], 'parent_id':row['parent_id'], 'rank_lower_bound':state['rank'],
            'source_status':result['status'], 'packet':io.read(job/'packet.json'),
            'packet_sha256':result['packet_sha256'],
            'certificate_replay':io.read(job/('packet-verified.json' if row['backend']=='native' else 'verified.json')),
            'final_state':state, 'source_bindings':{str(p.relative_to(ROOT)):io.sha(p) for p in
                (job/'result.json',job/'packet.json',job/'seal.json',job/'broad-state.json')}}
        validate_record(record, parents[row['parent_id']])
        records.append(record)
    require(len(records) == 30, 'expected requested thirty rank >22 fibres')
    io.write(SNAPSHOT, {'schema':'broad-rank.ledger-snapshot.v1', 'records':records,
        'parents':parents, 'baseline':baseline, 'baseline_database_sha256':io.sha(database),
        'run_bindings':{str((RUN/n).relative_to(ROOT)):io.sha(RUN/n) for n in
                        ('plan.json','manifest.json','queue.json','COMPLETE.json','COMPLETION_REVIEW.json')},
        'threshold_exclusive':22, 'rank_counts':dict(Counter(str(r['rank_lower_bound']) for r in records)),
        'boundary':'Certified subgroup lower bounds and exact Q-nonisomorphism relative to the pinned ledger. No exact rank, worldwide novelty, minimal-model metric or conductor claim.'})
    print('BROAD_LEDGER_FROZEN', len(records), flush=True)


def check():
    # A new in-memory finite cache for the second existing replay path avoids
    # trusting previous disk facts or publishing per-point arithmetic caches.
    import research_runtime.finite_reduction as finite
    from research_runtime.memory_store import MemoryFactStore
    data = io.read(SNAPSHOT)
    require(data['schema'] == 'broad-rank.ledger-snapshot.v1', 'unknown snapshot')
    baseline = {}
    for row in data['baseline']:
        baseline.setdefault(invariant(row['ainvs']), []).append(row)
    seen = {}; verified = []
    for row in data['records']:
        validate_record(row, data['parents'][row['parent_id']])
        packet = row['packet']; model = tuple(map(F, packet['curve']))
        points = tuple(tuple(map(F, p)) for p in packet['points']); proof = packet['proof']
        primes = [s['prime'] for s in proof['signatures']]; torsion = proof['no_rational_2_torsion_prime']
        first = checked_rank(model, points, primes, torsion)
        finite._default = finite.ReductionCache(MemoryFactStore())
        second = independent.checked_rank(model, points, primes, torsion)
        require(json.loads(json.dumps(first)) == proof and second['rank_lower_bound'] == row['rank_lower_bound'],
                'fresh rank replay differs')
        j = invariant(model)
        require(not any(independent.isomorphic(model, p['ainvs']) for p in baseline.get(j, [])),
                'candidate duplicates prior ledger')
        require(not any(independent.isomorphic(model, p) for p in seen.get(j, [])), 'duplicate cohort curve')
        seen.setdefault(j, []).append(model)
        verified.append({'id':row['id'],'rank_lower_bound':row['rank_lower_bound'],'packet_sha256':row['packet_sha256']})
        print('BROAD_LEDGER_RANK_REPLAY_PASS', row['id'], row['rank_lower_bound'], flush=True)
    require(len(verified) == 30 and len({r['id'] for r in verified}) == 30, 'cohort identity count differs')
    require(data['rank_counts'] == dict(Counter(str(r['rank_lower_bound']) for r in verified)), 'counts differ')
    io.write(REPLAY, {'status':'PASS_BROAD_LEDGER_FRESH_RANK_REPLAY', 'snapshot_sha256':io.sha(SNAPSHOT),
        'checker_sha256':io.sha(Path(__file__)), 'records':verified,
        'scope':'Both existing finite rank replay paths, all point memberships and generic prefixes, torsion exclusion, and pairwise/baseline Q-nonisomorphism. No point search or conductor calculation.'})
    print('PASS_BROAD_LEDGER_REPLAY', len(verified), flush=True)


def index():
    data, replay = io.read(SNAPSHOT), io.read(REPLAY)
    require(replay['snapshot_sha256'] == io.sha(SNAPSHOT) and replay['checker_sha256'] == io.sha(Path(__file__)),
            'fresh publication replay required')
    require(replay['status'] == 'PASS_BROAD_LEDGER_FRESH_RANK_REPLAY', 'publication replay did not pass')
    manifest = io.read(MANIFEST); rel = str(SNAPSHOT.relative_to(ROOT))
    claim = next(r for r in io.read(ROOT/'MATH_STATUS.json')['entries'] if r['id']==CLAIM)
    require(claim['state']=='proved' and rel in claim['software_lock'], 'source-bound ledger claim required')
    entries = {r['id']:r for r in manifest['curves']}
    for row in data['records']:
        entry = {k:row[k] for k in ('id','family','parameter','rank_lower_bound')}
        entry.update(kind='cohort', source=rel, record=row['id'], claim=CLAIM)
        require(row['id'] not in entries or entries[row['id']] == entry, 'existing ledger identity differs')
        entries[row['id']] = entry
    manifest['curves'] = list(entries.values()); manifest['sources'][rel] = io.sha(SNAPSHOT)
    io.write(MANIFEST, manifest, False)
    print('BROAD_LEDGER_INDEXED',len(data['records']),flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('freeze','check','index'))
    args = parser.parse_args()
    {'freeze':freeze,'check':check,'index':index}[args.action]()
