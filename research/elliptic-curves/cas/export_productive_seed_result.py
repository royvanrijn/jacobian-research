#!/usr/bin/env python3
"""Export sealed productive-seed results; this is not another arithmetic replay."""
import argparse
import hashlib
import json
from pathlib import Path
from research_runtime.store import checkpoint

ROOT = Path(__file__).resolve().parents[2]


def export(folder, output):
    def read(p):
        return json.loads(p.read_text())

    def sha(p):
        return hashlib.sha256(p.read_bytes()).hexdigest()

    protocol, terminal, verified = [read(folder/n) for n in
        ('protocol.json', 'terminal.json', 'verified.json')]
    if (verified['status'] != 'PASS_INDEPENDENT_PRODUCTIVE_V3_REPLAY'
            or verified['terminal_sha256'] != sha(folder/'terminal.json')
            or verified['rank_lower_bound'] != terminal['rank_lower_bound']
            or verified['charts'] != terminal['charts']):
        raise ArithmeticError('sealed independent replay required')
    paths = [folder/n for n in ('protocol.json', 'terminal.json', 'verified.json')]
    for suffix in ('-supervision', '-replay-supervision'):
        p = folder.with_name(folder.name+suffix)/'supervisor.json'
        if read(p)['outcome'] != 'completed':
            raise ArithmeticError('successful supervised completion required')
        paths.append(p)
    for stage in sorted(folder.glob('epoch-*/stage.json')):
        paths.extend([stage, stage.parent/'reference-verified.json',
                      stage.parent/'landscape/timings.json'])
        gain = stage.parent/'gain.json'
        if gain.exists():
            paths.append(gain)
    for category in ('inputs', 'sources'):
        if any(sha(ROOT/n) != h for n, h in protocol[category].items()):
            raise ArithmeticError('frozen binding changed')
    result = {'status': verified['status'], 'initial_rank': protocol['initial_rank'],
        'final_lower_bound': terminal['rank_lower_bound'], 'charts': terminal['charts'],
        'bindings': {str(p.relative_to(ROOT)): sha(p) for p in paths},
        'records': {str(p.relative_to(ROOT)): read(p) for p in paths},
        'scope': 'Verbatim sealed evidence export, not an additional arithmetic replay. '
                 'Certified subgroup lower bound; no exact rank or conductor record claim.'}
    if output.exists():
        if read(output) != result:
            raise ArithmeticError('preserve existing export')
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        checkpoint(output, result)
    print('EXPORTED', result['initial_rank'], result['final_lower_bound'], result['charts'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--folder', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    export(args.folder.resolve(), args.output.resolve())
