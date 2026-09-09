#!/usr/bin/env python3
"""Preserve a verified short-pass exposure cursor; does not launch continuation."""
import argparse
import hashlib
import json
from pathlib import Path


def record(folder):
    bindings = {}

    def read(path):
        raw = path.read_bytes()
        bindings[str(path.resolve())] = hashlib.sha256(raw).hexdigest()
        return json.loads(raw)

    protocol, terminal, verified = [read(folder/n) for n in
                                    ('protocol.json', 'terminal.json', 'verified.json')]
    if (protocol['schema'] not in ('prepared-parent-complement-boxes.v1', 'prepared-parent-complement-cached.v1') or
            terminal['protocol_sha256'] != bindings[str((folder/'protocol.json').resolve())] or
            verified['status'] != 'PASS_INDEPENDENT_PRODUCTIVE_V3_REPLAY' or
            verified['terminal_sha256'] != bindings[str((folder/'terminal.json').resolve())] or
            verified['rank_lower_bound'] != terminal['rank_lower_bound'] or
            verified['charts'] != terminal['charts']):
        raise ValueError('a sealed, independently verified short pass is required')
    for suffix in ('-supervision', '-replay-supervision'):
        supervisor = read(folder.with_name(folder.name+suffix)/'supervisor.json')
        if supervisor['outcome'] != 'completed' or supervisor['returncode'] != 0:
            raise ValueError('successful supervised completion required')
    if terminal['stop_reason'] != 'CHART_BUDGET_EXHAUSTED':
        raise ValueError('only budget-truncated passes enter this suffix queue')
    stage = terminal['stages'][-1]
    epoch = folder/f"epoch-{stage['epoch']:02d}"
    selection_path = epoch/'landscape/selection.json'
    selection = read(selection_path)
    if stage['selection_sha256'] != bindings[str(selection_path.resolve())]:
        raise ValueError('selection differs from verified terminal')
    last_path = epoch/f"chart-{stage['charts']-1:04d}.json"
    last = read(last_path)
    if stage['last_chart_sha256'] != bindings[str(last_path.resolve())]:
        raise ValueError('last receipt differs from verified terminal')
    output = dict(status='QUEUED_NOT_LAUNCHED', rank_lower_bound=terminal['rank_lower_bound'],
                  epoch=stage['epoch'], bindings=bindings,
                  censored_maps_in_final_epoch=stage['censored_maps'],
                  censored_maps_require_separate_review=True,
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  claim_boundary='Exact exposure cursor only. Use the terminal basis and proof. '
                  'Recheck projective coordinate duplication before executing a remaining map. '
                  'This suffix cursor does not resolve earlier censored maps. No suffix coverage, new rank or automatic execution is asserted.')
    if stage['after'] > stage['before']:
        output['next_action'] = 'REBUILD_FROM_CERTIFIED_GAIN'
    else:
        if selection['basis'] != terminal['points']:
            raise ValueError('unchanged epoch basis differs from final points')
        order = protocol['mapping_order']
        policy_index = order.index(last['policy'])+1
        centre = last['centre_index']
        if policy_index == len(order):
            centre += 1
            policy_index = 0
        if centre == len(selection['centres']):
            output['next_action'] = 'NO_UNVISITED_CENTRE_SUFFIX'
        elif 0 <= centre < len(selection['centres']):
            output.update(next_action='CONSIDER_REMAINING_SUFFIX',
                          next_centre_index=centre, next_policy_to_consider=order[policy_index])
        else:
            raise ValueError('invalid centre cursor')
    return output


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--folder', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    result = record(args.folder.resolve())
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write('\n')
    print(result['status'], result['next_action'], result.get('next_centre_index'),
          result.get('next_policy_to_consider'))
