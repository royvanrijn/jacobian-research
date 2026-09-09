#!/usr/bin/env python3
"""Read-only post-execution audit; never selects centres or launches searches.

Use Sage -python with --proofs for fresh native-prefix/independence checks.
--complete rejects partial panels or resource-unresolved cases.
"""
import argparse
from fractions import Fraction
import json

import run_r17_60_panel as panel
from v3_warm_support import read, require, sha


def bindings(base, records):
    for name, digest in records.items():
        require(sha(base / name) == digest, f'evidence changed: {base / name}')


def case_audit(row, proofs=False):
    folder = panel.D / 'cases' / row['id']
    result = read(folder / 'result.json')
    require(result['status'] == 'PASS_INDEPENDENT_R17_60_CASE', 'case unresolved: ' + row['id'])
    for key in ('id', 'family', 'parameter', 'stratum'):
        require(result[key] == row[key], 'roster identity differs: ' + key)
    require(result['protocol_sha256'] == sha(panel.D / 'protocol.json'), 'case protocol differs')
    bindings(panel.ROOT, result['bindings'])
    seed_folder = folder / 'seed-search'
    seed = read(seed_folder / 'terminal.json')
    sp = read(seed_folder / 'protocol.json')
    replay = read(seed_folder / 'verified.json')
    require(replay['status'] == 'PASS_INDEPENDENT_FIRST_SEED_REPLAY' and
            replay['terminal_sha256'] == sha(seed_folder / 'terminal.json'), 'seed replay differs')
    require(seed['protocol_sha256'] == sha(seed_folder / 'protocol.json'), 'seed protocol changed')
    require(0 <= seed['charts'] <= sp['max_point_invocations'] <= 98, 'seed allowance exceeded')
    require(seed['rank_lower_bound'] in (17, 18), 'first-seed endpoint differs')
    require(result['first_M18_found'] == (seed['rank_lower_bound'] == 18), 'seed incidence differs')
    require(len(list(seed_folder.glob('chart-*.json'))) == seed['charts'], 'seed chart count differs')
    packet = read(folder / 'seed-reconciled.json')
    bindings(seed_folder, packet['bindings'])
    require(packet['point_searches'] == 0, 'seed reconciliation searched')
    require(packet['rank_lower_bound'] == seed['rank_lower_bound'] + len(packet['reconciliation_gains']),
            'seed cloud rank accounting differs')
    timeline = []
    if seed['rank_lower_bound'] == 18:
        timeline.append({'phase': 'seed', 'call': seed['charts'], 'before': 17, 'after': 18})
    for gain in packet['reconciliation_gains']:
        require(1 <= gain['seed_call'] <= seed['charts'], 'invalid seed-cloud call')
        timeline.append({'phase': 'seed-cloud', 'call': gain['seed_call'],
                         'before': gain['after'] - 1, 'after': gain['after']})
    used = 0
    segments = sorted(folder.glob('complement-[0-9][0-9]'))
    for index, run in enumerate(segments):
        require(run.name == f'complement-{index:02d}', 'nonconsecutive segments')
        policy = read(run / 'protocol.json')
        require(policy['max_charts'] == 100 - used, 'remaining budget differs')
        bindings(panel.ROOT, policy['inputs'])
        bindings(panel.ROOT, policy['sources'])
        terminal = read(run / 'terminal.json')
        require(terminal['protocol_sha256'] == sha(run / 'protocol.json'), 'terminal policy changed')
        replay = read(run / 'verified.json')
        require(replay['status'] == 'PASS_INDEPENDENT_PRODUCTIVE_V3_REPLAY' and
                replay['terminal_sha256'] == sha(run / 'terminal.json'), 'complement replay differs')
        require(0 <= terminal['charts'] <= 100 - used, 'complement allowance exceeded')
        before = packet['rank_lower_bound']
        offsets = {}
        subtotal = 0
        for stage in terminal['stages']:
            require(stage['before'] == before, 'subgroup chain differs')
            epoch = stage['epoch']
            offsets[epoch] = subtotal
            wd = run / f'epoch-{epoch:02d}'
            require(len(list(wd.glob('chart-*.json'))) == stage['charts'], 'epoch call count differs')
            epoch_replay = read(wd / 'reference-verified.json')
            require(epoch_replay['status'] == 'PASS_FULL_PRODUCTIVE_EPOCH_REPLAY' and
                    epoch_replay['stage_sha256'] == sha(wd / 'stage.json') and read(wd / 'stage.json') == stage,
                    'full epoch replay missing')
            subtotal += stage['charts']
            before = stage['after']
            if stage['after'] > stage['before']:
                timeline.append({'phase': 'complement', 'call': used + subtotal,
                                 'before': stage['before'], 'after': stage['after']})
        require(subtotal == terminal['charts'] and before == terminal['rank_lower_bound'],
                'terminal count/rank differs')
        packet = read(folder / f'reconciled-{index:02d}.json')
        bindings(run, packet['bindings'])
        require(packet['point_searches'] == 0 and packet['rank_lower_bound'] == before + len(packet['gains']),
                'cloud count/rank differs')
        for gain in packet['gains']:
            epoch, chart = gain['chart'].split('/')
            call = offsets[int(epoch.split('-')[1])] + int(chart.removesuffix('.json').split('-')[1]) + 1
            require(1 <= call <= subtotal, 'invalid cloud call')
            timeline.append({'phase': 'complement-cloud', 'call': used + call,
                             'before': gain['after'] - 1, 'after': gain['after']})
        used += subtotal
        if index + 1 < len(segments):
            require(terminal['stop_reason'] == 'ADDITIONAL_FINITE_RANK_REQUIRES_RECONCILIATION' and
                    packet['rank_lower_bound'] > terminal['rank_lower_bound'], 'unallowed suffix')
    require(result['packet'] == packet and result['gain_timeline'] == timeline, 'reported packet/timeline differs')
    require(result['rank_lower_bound'] == packet['rank_lower_bound'], 'reported rank differs')
    require(result['seed_calls'] == seed['charts'] and result['complement_calls'] == used and
            result['total_calls'] == seed['charts'] + used, 'reported budget differs')
    added = packet['rank_lower_bound'] - 17
    require(result['added_directions'] == added, 'added rank differs')
    ratio = str(Fraction(seed['charts'] + used, added)) if added else None
    require(result['calls_per_added_direction'] == ratio, 'calls per gain differs')
    calls = {g['call'] for g in timeline if g['phase'].startswith('complement')}
    last = max(calls, default=None)
    require(result['gaining_complement_calls'] == len(calls) and result['last_complement_gain_call'] == last,
            'reported gain timing differs')
    require(result['complement_calls_since_last_gain'] == used - (last or 0), 'trailing stall differs')
    if proofs:
        from r17_60_arithmetic import native_check
        native_check(packet, row)
    return {k: result[k] for k in ('id', 'rank_lower_bound', 'total_calls', 'last_complement_gain_call')}


def audit(complete=False, proofs=False):
    panel.check_protocol()
    require(read(panel.ART / 'r17_60_panel_protocol_v1.json') == read(panel.D / 'protocol.json'),
            'published protocol differs')
    roster = read(panel.D / 'roster.json')
    require(roster == read(panel.ART / 'r17_60_panel_roster_v1.json'), 'published roster differs')
    selection = panel.selection.select(read(panel.D / 'selection-input.json'))
    require(selection == {k: roster[k] for k in ('rows', 'skipped', 'stratum_pool_counts')}, 'selection replay differs')
    gate = read(panel.D / 'preflight.json')
    require(gate['protocol_sha256'] == sha(panel.D / 'protocol.json') and gate['point_searches'] == 0 and
            gate['generic_verified_sha256'] == sha(panel.D / 'generic/verified.json'), 'preflight differs')
    generic = read(panel.D / 'generic/prepared.json')
    require(len(generic['records']) == 60, 'generic inputs incomplete')
    for record in generic['records']:
        require(sha(panel.ROOT / record['path']) == record['sha256'], 'generic input changed')
    results, pending, unresolved = [], [], []
    for row in roster['rows']:
        path = panel.D / 'cases' / row['id'] / 'result.json'
        if not path.exists():
            pending.append(row['id'])
        elif read(path)['status'] != 'PASS_INDEPENDENT_R17_60_CASE':
            unresolved.append(row['id'])
        else:
            results.append(case_audit(row, proofs))
    if complete:
        require(not pending and not unresolved and len(results) == 60, 'sixty independently verified outcomes required')
        report = read(panel.ART / 'r17_60_panel_results_v1.json')
        require(report['completed'] == report['independently_verified'] == 60, 'final report incomplete')
        require(report['results'] == [read(panel.D / 'cases' / row['id'] / 'result.json') for row in roster['rows']],
                'exported results differ')
    return {'status': 'PASS_COMPLETE_PANEL_AUDIT' if complete else 'PASS_PARTIAL_PANEL_AUDIT',
            'fresh_sage_proofs': proofs, 'verified': results, 'pending': pending, 'unresolved': unresolved,
            'protocol_sha256': sha(panel.D / 'protocol.json'), 'point_searches': 0}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--complete', action='store_true')
    parser.add_argument('--proofs', action='store_true')
    args = parser.parse_args()
    print(json.dumps(audit(args.complete, args.proofs), indent=2))
