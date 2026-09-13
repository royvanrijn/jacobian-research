"""Two-stage adapter to the retained, independently replayed V3 machinery.

Stage one stops searching after the first gaining cloud. Stage two has a fixed
call budget. Neither stage changes parent selection or replaces an address.
"""
import argparse
from pathlib import Path

from parent_foundry_worker import (ROOT, bank, generic_seed, save, verify)
from v3_warm_support import atomic, read, require, sha


def evaluate(job):
    request = read(job/'request.json')
    parent_path = ROOT/request['parent']
    require(sha(parent_path) == request['parent_sha256'], 'parent input changed')
    parent = read(parent_path)
    if request.get('packet'):
        path = ROOT/request['packet']
        require(sha(path) == request['packet_sha256'], 'sealed seed changed')
        packet = read(path)
        verify(packet, parent, request['parameter'])
    else:
        packet = generic_seed(parent, request['parameter'], job)
    if packet is None:
        save(job/'result.json', {'status': 'UNRESOLVED_GENERIC_SPECIALIZATION',
             'rank_lower_bound': None, 'calls': 0, 'binary_outcome_complete': False,
             'request_sha256': sha(job/'request.json')})
        return
    initial = packet['rank_lower_bound']
    used, segments, events = 0, [], []
    import run_complement_seed_v3 as runner
    from parent_foundry_seed import freeze
    import pari_pointed_backend as backend
    from high_rank_foundry_point_calls import install
    from reconcile_verified_v3_cloud import run as reconcile
    from high_rank_foundry_certificate import compact_packet
    if request['mapper'] == 'factor_free':
        runner.MAPPERS = {'factor_free': runner.MAPPERS['factor_free']}
    else:
        require(request['mapper'] == 'dual', 'undeclared mapping policy')
    install(job, backend)
    for segment in range(8):
        if used >= request['allowance'] or packet['rank_lower_bound'] >= 32:
            break
        if request['phase'] == 'seed' and packet['rank_lower_bound'] > 16:
            break
        prep, dest = job/f'bank-{segment:02d}', job/f'search-{segment:02d}'
        bank(parent, packet, prep, request['bank_index']*8+segment)
        runner.PREP = runner.BANK = prep
        p = freeze(dest, prep)
        p.update(max_charts=request['allowance']-used,
                 target_rank=17 if request['phase'] == 'seed' else 32,
                 max_epochs=1 if request['phase'] == 'seed' else 33-packet['rank_lower_bound'],
                 mapping_order=list(runner.MAPPERS))
        p['inputs'][str((job/'request.json').relative_to(ROOT))] = sha(job/'request.json')
        p['sources'][str(Path(__file__).relative_to(ROOT))] = sha(Path(__file__))
        atomic(dest/'protocol.json', p)
        runner.search(dest)
        runner.replay(dest)
        reconciled = job/f'reconciled-{segment:02d}.json'
        reconcile(dest, reconciled)
        terminal, following = read(dest/'terminal.json'), read(reconciled)
        n = terminal['charts']
        require(0 <= n <= request['allowance']-used, 'call allowance exceeded')
        if request['phase'] == 'seed':
            require(len(terminal['stages']) <= 1, 'stage one must not amplify')
        offset = used
        offsets = {}
        for stage in terminal['stages']:
            offsets[stage['epoch']] = offset
            offset += stage['charts']
            if stage['after'] > stage['before']:
                events.append({'before': stage['before'], 'after': stage['after'],
                               'call': offset, 'phase': 'online'})
        for gain in following['gains']:
            chart = Path(gain['chart'])
            epoch = int(chart.parts[0].split('-')[1])
            index = int(chart.stem.split('-')[1])
            events.append({'before': gain['after']-1, 'after': gain['after'],
                           'call': offsets[epoch]+index+1, 'phase': 'cloud-reconciliation'})
        used += n
        segments.append({'search': str(dest.relative_to(ROOT)), 'calls': n,
            'before': packet['rank_lower_bound'], 'after': following['rank_lower_bound'],
            'point_timeouts': sum(s['censored'] for s in terminal['stages']),
            'map_timeouts': sum(s['censored_maps'] for s in terminal['stages']),
            'stop': terminal['stop_reason']})
        packet = following
        packet.update(family=parent['family'], parameter=request['parameter'])
        require(packet['rank_lower_bound'] >= max([terminal['rank_lower_bound']]+[
            r for s in terminal['stages'] for r in s['odd_ranks'].values()]),
            'stronger finite cloud not reconciled')
        if n == 0:
            break
    packet = compact_packet(packet)
    verify(packet, parent, request['parameter'])
    packet['status'] = 'PASS_TWO_FINITE_IMPLEMENTATIONS'
    save(job/'packet.json', packet)
    censored = any(s['point_timeouts'] or s['map_timeouts'] for s in segments)
    gained = packet['rank_lower_bound'] > initial
    complete = (request['phase'] == 'seed' and gained) or (used == request['allowance'] and not censored) or packet['rank_lower_bound'] >= 32
    save(job/'verified.json', {'status': 'PASS_TWO_FINITE_IMPLEMENTATIONS',
        'packet_sha256': sha(job/'packet.json'), 'request_sha256': sha(job/'request.json')})
    save(job/'result.json', {'status': 'PASS_CERTIFIED_EXPOSURE',
        'family': parent['family'], 'parameter': request['parameter'], 'phase': request['phase'],
        'initial_rank': initial, 'rank_lower_bound': packet['rank_lower_bound'],
        'calls': used, 'allowance': request['allowance'], 'segments': segments,
        'gain_timeline': events, 'binary_outcome_complete': complete,
        'point_or_map_censoring': censored, 'packet': str((job/'packet.json').relative_to(ROOT)),
        'packet_sha256': sha(job/'packet.json'), 'request_sha256': sha(job/'request.json'),
        'boundary': 'First gaining cloud is reconciled without further searches in stage one. '
                    'No seed under exposure is not a rank upper bound. No automatic retry.'})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--job', type=Path, required=True)
    parser.add_argument('--verify-only', action='store_true')
    args = parser.parse_args()
    job = args.job.resolve()
    if args.verify_only:
        request = read(job/'request.json')
        require(sha(ROOT/request['parent']) == request['parent_sha256'], 'parent changed')
        require(read(job/'verified.json')['packet_sha256'] == sha(job/'packet.json'), 'packet changed')
        verify(read(job/'packet.json'), read(ROOT/request['parent']), request['parameter'])
    else:
        evaluate(job)
