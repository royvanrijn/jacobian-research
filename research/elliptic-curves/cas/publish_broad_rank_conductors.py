#!/usr/bin/env python3
"""Freeze, freshly replay, and index the first broad-rank conductor cutoff."""
import argparse
from fractions import Fraction as F
from pathlib import Path
import tempfile

from publish_broad_rank_results import SNAPSHOT as CURVES, digest
from research_curve_refresh import MANIFEST, ROOT
import run_broad_rank_conductors as queue
from v3_warm_support import atomic, read, require, sha

SNAPSHOT = queue.OUT.parent/'broad_rank_conductor_snapshot_v1.json'
REPLAY = SNAPSHOT.with_name('broad_rank_conductor_replay_v1.json')
CLAIM = 'EC-BROAD-RANK-CONDUCTORS-V1-20260912'
SELF = Path(__file__).resolve()


def validate(identifier, entry, curves, protocol):
    source, result = entry['input'], entry['result']
    certificate = result['certificate']
    require(identifier in protocol['ids'] and identifier in curves, 'curve outside requested cohort')
    require(identifier == source['id'] == result['id'] == certificate['id'], 'conductor identity differs')
    require(result['status'] == 'PASS_INDEPENDENT_CONDUCTOR_REPLAY'
            and result['protocol_sha256'] == digest(protocol), 'unsealed conductor endpoint')
    require(certificate['status'] in ('EXACT', 'UNKNOWN'), 'invalid conductor status')
    require(digest(source) == entry['input_sha256'] == certificate['input_sha256']
            == protocol['inputs'][entry['input_origin']], 'conductor input binding differs')
    require(digest(result) == entry['origin_sha256'], 'conductor result binding differs')
    require(tuple(map(F, source['original_curve'])) == tuple(map(F, curves[identifier]['packet']['curve'])),
            'conductor source differs from published model')
    require(source['rank_lower_bound'] == certificate['rank_lower_bound']
            == curves[identifier]['rank_lower_bound'], 'conductor rank label differs')


def freeze():
    require(not SNAPSHOT.exists(), 'preserve conductor cutoff; choose a new version')
    protocol = queue.guard()
    curves = {r['id']:r for r in read(CURVES)['records']}
    # Take a single path cutoff. Later completions stay in the live queue.
    paths = sorted(queue.OUT.glob('broad-*.json'))
    records = {}
    for path in paths:
        result = read(path)
        if result['status'] != 'PASS_INDEPENDENT_CONDUCTOR_REPLAY':
            continue
        identifier = result['id']; case = queue.D/'cases'/identifier
        require(read(case/'result.json') == result, 'generated and local results differ')
        for name, expected in result['bindings'].items():
            require(sha(ROOT/name) == expected, 'conductor source binding changed: '+name)
        require(read(case/'replay.supervisor.json')['outcome'] == 'completed', 'replay did not complete')
        require(read(case/'conductor.json') == result['certificate'], 'embedded certificate differs')
        entry = dict(input=read(case/'input.json'), result=result,
                     origin=str(path.relative_to(ROOT)), origin_sha256=sha(path),
                     input_origin=str((case/'input.json').relative_to(ROOT)), input_sha256=sha(case/'input.json'))
        validate(identifier, entry, curves, protocol)
        records[identifier] = entry
    require(records, 'no completed conductor certificates')
    exact = sum(e['result']['certificate']['status'] == 'EXACT' for e in records.values())
    atomic(SNAPSHOT, dict(schema='broad-rank.conductor-snapshot.v1', records=records,
        record_count=len(records), exact_count=exact, unknown_count=len(records)-exact,
        curve_snapshot_sha256=sha(CURVES), protocol=protocol,
        boundary='Fixed completed-path cutoff; later queue results excluded. Exact only after complete local arithmetic and independent replay. Partial conductors remain UNKNOWN.'), immutable=True)
    print('BROAD_CONDUCTOR_CUTOFF_FROZEN', len(records), exact, 'exact', flush=True)


def check():
    from inventory_conductor_worker_v2 import run as replay_conductor
    data = read(SNAPSHOT)
    require(data['schema'] == 'broad-rank.conductor-snapshot.v1'
            and data['curve_snapshot_sha256'] == sha(CURVES), 'conductor snapshot binding differs')
    curves = {r['id']:r for r in read(CURVES)['records']}
    with tempfile.TemporaryDirectory(prefix='broad-conductor-publication-') as directory:
        folder = Path(directory)
        for identifier, entry in sorted(data['records'].items()):
            validate(identifier, entry, curves, data['protocol'])
            atomic(folder/'input.json', entry['input'])
            atomic(folder/'certificate.json', entry['result']['certificate'])
            replay_conductor(folder/'input.json', folder/'certificate.json', True)
            print('BROAD_CONDUCTOR_REPLAY_PASS', identifier, flush=True)
    exact = sum(e['result']['certificate']['status'] == 'EXACT' for e in data['records'].values())
    require(data['record_count'] == len(data['records']) and data['exact_count'] == exact
            and data['unknown_count'] == len(data['records'])-exact, 'conductor counts differ')
    atomic(REPLAY, dict(status='PASS_BROAD_CONDUCTOR_PUBLICATION_REPLAY',
        snapshot_sha256=sha(SNAPSHOT), checker_sha256=sha(SELF),
        record_count=len(data['records']), exact_count=exact,
        scope='Fresh exact transport, prime-certificate and Sage generic Tate/PARI local conductor replay, with no factor discovery.'), immutable=True)
    print('PASS_BROAD_CONDUCTOR_PUBLICATION_REPLAY',len(data['records']),exact,'exact',flush=True)


def index():
    data, replay = read(SNAPSHOT), read(REPLAY)
    require(replay['status'] == 'PASS_BROAD_CONDUCTOR_PUBLICATION_REPLAY'
            and replay['snapshot_sha256'] == sha(SNAPSHOT) and replay['checker_sha256'] == sha(SELF),
            'fresh conductor publication replay required')
    manifest = read(MANIFEST); rel = str(SNAPSHOT.relative_to(ROOT))
    claim = next(r for r in read(ROOT/'MATH_STATUS.json')['entries'] if r['id'] == CLAIM)
    require(claim['state'] == 'proved' and rel in claim['software_lock'], 'source-bound conductor claim required')
    entries = {r['id']:r for r in manifest['conductors']}
    for identifier in data['records']:
        entry = dict(id=identifier,source=rel,claim=CLAIM,kind='checkpoint')
        require(identifier not in entries or entries[identifier] == entry, 'existing conductor publication differs')
        entries[identifier] = entry
    manifest['conductors'] = list(entries.values()); manifest['sources'][rel] = sha(SNAPSHOT)
    atomic(MANIFEST, manifest)
    print('BROAD_CONDUCTORS_INDEXED', len(data['records']), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('freeze','check','index'))
    args = parser.parse_args()
    {'freeze':freeze,'check':check,'index':index}[args.action]()
