#!/usr/bin/env python3
"""Freeze, replay, and index the completed conductor pass for foundry-v3 V2."""
import argparse
from pathlib import Path
import tempfile

import certify_compact_r17_candidates as cert
import refresh_foundry_curve_ledger_v2 as ledger
from research_curve_refresh import MANIFEST, ROOT
from v3_warm_support import atomic, read, require, sha

ART = ROOT/'artifacts/generated-results/elliptic-curves'
CURVES = ledger.CURVES
SOURCE = ART/'foundry_v3_conductors_v1'
LOCAL = ROOT/'artifacts/local/elliptic-curves/foundry-v3-conductors-v1'
SNAPSHOT = ART/'foundry_v3_conductor_snapshot_v2.json'
CLAIM = 'EC-FOUNDRY-V3-CONDUCTORS-V2-20260912'
SELF = Path(__file__).resolve()


def freeze():
    require(not SNAPSHOT.exists(), 'preserve immutable snapshot; choose a new version')
    state = read(LOCAL/'state.json')
    require(state['status'] == 'COMPLETE_SOURCE_STOPPED' and not state['active'],
            'completed stopped-source conductor queue required')
    curves = read(CURVES)
    require(curves['status'] == 'FROZEN_FOUNDRY_LEDGER_SELECTION_V2', 'invalid curve snapshot')
    records = {}
    for identifier, entry in curves['records'].items():
        source_id = entry['result']['id']
        result_path = SOURCE/f'{source_id}.json'
        input_path = LOCAL/'cases'/source_id/'input.json'
        require(result_path.exists() and input_path.exists(), 'selected endpoint lacks conductor result')
        result, packet = read(result_path), read(input_path)
        require(result['status'] == 'PASS_INDEPENDENT_CONDUCTOR_REPLAY', 'independent conductor replay required')
        require(result['certificate']['id'] == source_id == packet['id'], 'conductor source identity differs')
        require(cert.isomorphic(result['certificate']['curve'], entry['result']['packet']['curve']),
                'conductor curve differs from selected foundry curve')
        records[identifier] = {
            'input': packet,
            'result': result,
            'origin': str(result_path.relative_to(ROOT)),
            'origin_sha256': sha(result_path),
            'input_origin': str(input_path.relative_to(ROOT)),
            'input_sha256': sha(input_path),
            'source_foundry_id': source_id,
        }
    exact = sum(entry['result']['certificate']['status'] == 'EXACT' for entry in records.values())
    require(len(records) == 94 and exact == 71, 'completed selected cohort counts differ')
    atomic(SNAPSHOT, {
        'status': 'PASS_COMPLETE_FOUNDRY_V3_CONDUCTOR_CUTOFF',
        'records': records,
        'record_count': len(records),
        'exact_count': exact,
        'unknown_count': len(records)-exact,
        'curve_snapshot_sha256': sha(CURVES),
        'queue_state_sha256': sha(LOCAL/'state.json'),
        'selection': 'All independently replayed conductor endpoints for the fixed 94-curve foundry-v3 V2 selection after its stopped source queue completed.',
        'boundary': 'Exact only when the remaining cofactor is one. Partial endpoints remain UNKNOWN bounds. No exact-rank or conductor-record claim.',
    }, immutable=True)
    print('FOUNDRY_V3_CONDUCTORS_V2_FROZEN', len(records), 'records;', exact, 'exact', flush=True)


def check():
    from inventory_conductor_worker_v2 import run as replay_conductor

    curves, data = read(CURVES), read(SNAPSHOT)
    require(data['status'] == 'PASS_COMPLETE_FOUNDRY_V3_CONDUCTOR_CUTOFF', 'invalid conductor snapshot')
    require(data['curve_snapshot_sha256'] == sha(CURVES), 'curve snapshot binding differs')
    require(data['record_count'] == len(data['records']) == 94, 'record count differs')
    with tempfile.TemporaryDirectory(prefix='foundry-v3-conductor-v2-replay-') as directory:
        folder = Path(directory)
        for identifier, entry in sorted(data['records'].items()):
            require(identifier in curves['records'], 'conductor curve outside fixed cohort')
            require(sha(ROOT/entry['origin']) == entry['origin_sha256'], 'conductor result changed')
            require(read(ROOT/entry['origin']) == entry['result'], 'embedded conductor result differs')
            require(sha(ROOT/entry['input_origin']) == entry['input_sha256'], 'conductor input changed')
            require(read(ROOT/entry['input_origin']) == entry['input'], 'embedded conductor input differs')
            atomic(folder/'input.json', entry['input'])
            atomic(folder/'certificate.json', entry['result']['certificate'])
            replay_conductor(folder/'input.json', folder/'certificate.json', True)
    exact = sum(entry['result']['certificate']['status'] == 'EXACT' for entry in data['records'].values())
    require(data['exact_count'] == exact == 71 and data['unknown_count'] == 23, 'status counts differ')
    print('PASS_FOUNDRY_V3_CONDUCTORS_V2_REPLAY', data['record_count'], 'records;', exact, 'exact', flush=True)


def index():
    data = read(SNAPSHOT)
    manifest = read(MANIFEST)
    rel = str(SNAPSHOT.relative_to(ROOT))
    claim = next(row for row in read(ROOT/'MATH_STATUS.json')['entries'] if row['id'] == CLAIM)
    require(claim['state'] == 'proved' and rel in claim['software_lock'], 'proved source-bound claim required')
    entries = {row['id']: row for row in manifest['conductors']}
    for identifier in data['records']:
        entries[identifier] = {'id': identifier, 'source': rel, 'claim': CLAIM, 'kind': 'checkpoint'}
    manifest['conductors'] = list(entries.values())
    manifest['sources'][rel] = sha(SNAPSHOT)
    atomic(MANIFEST, manifest)
    print('FOUNDRY_V3_CONDUCTORS_V2_INDEXED', len(data['records']), '71 exact; 23 partial', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('freeze', 'check', 'index'))
    args = parser.parse_args()
    {'freeze': freeze, 'check': check, 'index': index}[args.action]()
