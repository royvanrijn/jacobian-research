#!/usr/bin/env python3
"""Publish the completed broad conductor pass using the preserved V1 replay."""
import argparse
from pathlib import Path

import publish_broad_rank_conductors as v1
from v3_warm_support import atomic, read, require, sha

PREVIOUS = v1.SNAPSHOT
SNAPSHOT = PREVIOUS.with_name('broad_rank_conductor_snapshot_v2.json')
REPLAY = PREVIOUS.with_name('broad_rank_conductor_replay_v2.json')
CLAIM = 'EC-BROAD-RANK-CONDUCTORS-V2-20260912'
SELF = Path(__file__).resolve()
COMMON_SHA256 = '5d1cc7fe4665e7217683cbf402f8525816e2711fd25e2cb7cb610338fb2a1be8'


def configure():
    require(sha(Path(v1.__file__)) == COMMON_SHA256, 'preserved publication checker changed')
    v1.SNAPSHOT, v1.REPLAY, v1.CLAIM, v1.SELF = SNAPSHOT, REPLAY, CLAIM, SELF


def complete_snapshot():
    data = read(SNAPSHOT)
    require(set(data['records']) == set(data['protocol']['ids'])
            and len(data['records']) == data['record_count'] == 30,
            'complete thirty-curve cohort required')
    require(data['exact_count'] == 26 and data['unknown_count'] == 4, 'completed counts differ')
    for identifier, entry in read(PREVIOUS)['records'].items():
        require(data['records'][identifier] == entry, 'previous published certificate changed')
    for entry in data['records'].values():
        if entry['result']['certificate']['status'] == 'UNKNOWN':
            require(entry['result']['build_outcome'] == 'strict_wall_timeout',
                    'unexpected unresolved endpoint')
    return data


def freeze():
    protocol = v1.queue.guard()
    state = read(v1.queue.D/'state.json')
    require(state['status'] == 'COMPLETE' and not state['active'], 'completed conductor queue required')
    require(len(protocol['ids']) == 30, 'unexpected queue roster')
    v1.freeze()
    complete_snapshot()


def check():
    complete_snapshot()
    v1.check()


def index():
    data, replay = complete_snapshot(), read(REPLAY)
    require(replay['status'] == 'PASS_BROAD_CONDUCTOR_PUBLICATION_REPLAY'
            and replay['snapshot_sha256'] == sha(SNAPSHOT) and replay['checker_sha256'] == sha(SELF),
            'fresh completed-pass replay required')
    manifest = read(v1.MANIFEST)
    rel, previous_rel = (str(p.relative_to(v1.ROOT)) for p in (SNAPSHOT, PREVIOUS))
    require(manifest['sources'][previous_rel] == sha(PREVIOUS), 'earlier snapshot binding changed')
    claim = next(r for r in read(v1.ROOT/'MATH_STATUS.json')['entries'] if r['id'] == CLAIM)
    require(claim['state'] == 'proved' and rel in claim['software_lock'], 'proved source-bound claim required')
    entries = {r['id']:r for r in manifest['conductors']}
    previous_ids = set(read(PREVIOUS)['records'])
    for identifier in data['records']:
        entry = dict(id=identifier,source=rel,claim=CLAIM,kind='checkpoint')
        old = entries.get(identifier)
        require(old is None or old == entry or
                (identifier in previous_ids and old == dict(id=identifier,source=previous_rel,
                 claim='EC-BROAD-RANK-CONDUCTORS-V1-20260912',kind='checkpoint')),
                'unexpected existing conductor publication')
        entries[identifier] = entry
    manifest['conductors'] = list(entries.values())
    manifest['sources'][rel] = sha(SNAPSHOT)
    atomic(v1.MANIFEST, manifest)
    print('BROAD_CONDUCTORS_V2_INDEXED',len(data['records']),'26 exact; 4 partial',flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=('freeze','check','index'))
    args = parser.parse_args()
    configure()
    {'freeze':freeze,'check':check,'index':index}[args.action]()
