#!/usr/bin/env python3
"""Read-only cost attribution for completed future V3 runs; no point searches."""
import argparse
import hashlib
import json
from pathlib import Path


def audit(folders):
    bindings = {}

    def read(path):
        raw = path.read_bytes()
        bindings[str(path.resolve())] = hashlib.sha256(raw).hexdigest()
        return json.loads(raw)

    rows = []
    for folder in folders:
        terminal = read(folder / 'terminal.json')
        supervisor = read(folder.with_name(folder.name + '-supervision') / 'supervisor.json')
        if supervisor['outcome'] != 'completed' or supervisor['returncode'] != 0:
            raise ValueError(f'incomplete search: {folder}')
        epochs = []
        for stage in terminal['stages']:
            epoch = folder / f"epoch-{stage['epoch']:02d}"
            timing = read(epoch / 'landscape/timings.json')
            previous = hashlib.sha256((epoch / 'landscape/selection.json').read_bytes()).hexdigest()
            backend = cpu = 0.0
            statuses = {}
            for i in range(stage['charts']):
                path = epoch / f'chart-{i:04d}.json'
                chart = read(path)
                if chart['index'] != i or chart['previous_sha256'] != previous:
                    raise ValueError(f'chart chain differs: {path}')
                previous = bindings[str(path.resolve())]
                search = chart['search']
                backend += search['wall_seconds']
                cpu += (search['search_cpu_ms'] or 0) / 1000
                statuses[search['status']] = statuses.get(search['status'], 0) + 1
            if previous != stage['last_chart_sha256']:
                raise ValueError(f'terminal chain differs: {folder}')
            epochs.append(dict(epoch=stage['epoch'], before=stage['before'],
                               after=stage['after'], charts=stage['charts'],
                               backend_wall_seconds=backend, backend_cpu_seconds=cpu,
                               landscape_wall_seconds=timing['landscape_wall_seconds'],
                               statuses=statuses))
        backend = sum(e['backend_wall_seconds'] for e in epochs)
        landscape = sum(e['landscape_wall_seconds'] for e in epochs)
        wall = supervisor['wall_seconds']
        if backend + landscape > wall + 0.01:
            raise ValueError('component times exceed sequential search wall time')
        verified_path = folder / 'verified.json'
        verified = read(verified_path) if verified_path.exists() else None
        if verified and verified['terminal_sha256'] != bindings[str((folder/'terminal.json').resolve())]:
            raise ValueError('verification binds another terminal')
        rows.append(dict(run=folder.name, epochs=epochs, search_wall_seconds=wall,
                         backend_wall_seconds=backend, backend_fraction=backend / wall,
                         landscape_wall_seconds=landscape,
                         other_wall_seconds=wall - backend - landscape,
                         replay_status=verified['status'] if verified else 'PENDING'))
    return dict(status='PASS_RETAINED_COST_ACCOUNTING', runs=rows, bindings=bindings,
                source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                point_searches=0,
                claim_boundary='Recorded sequential timings only, not a controlled speedup '
                'benchmark or arithmetic replay. Other time is not further attributed. '
                'CPU time is not subtracted from wall time. No rank upper bound.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, action='append', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('preserve earlier cost reports')
    report = audit(args.run)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write('\n')
    for row in report['runs']:
        print(row['run'], 'wall', round(row['search_wall_seconds'], 3),
              'backend%', round(100 * row['backend_fraction'], 2),
              'landscape', round(row['landscape_wall_seconds'], 3),
              'other', round(row['other_wall_seconds'], 3), row['replay_status'])
