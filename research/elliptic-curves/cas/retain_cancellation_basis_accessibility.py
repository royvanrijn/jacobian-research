#!/usr/bin/env python3
"""Retain/replay the finite causal audit, without dispatching point searches."""
import argparse
from copy import deepcopy
import json
from pathlib import Path
import tarfile
import tempfile
import time

import cancellation_basis_accessibility_v3 as audit
from cancellation_basis_accessibility import ROOT, LOCAL, PRIOR, read, sha, need, new_write
from finite_cancellation_corpus import digest, write

OUT, RAW = audit.OUT, audit.RAW


def regress():
    from verify_cancellation_basis_accessibility_v3 import verify
    start = time.process_time(); plan = audit.base.guard(); summary = read(OUT/'summary.json')
    example = summary['visible_examples'][0]; row = next(r for r in plan['cases'] if r['id'] == example['case'])
    folder = RAW/row['id']; maps = read(folder/'maps.json'); scan = read(folder/'scan.json'); rows = []
    for kind in ('removed_new_generator', 'changed_minimum_coordinate'):
        m, s = deepcopy(maps), deepcopy(scan)
        if kind == 'removed_new_generator':
            entry = next(e for e in m['entries'] if e['id'] == example['best_new']['entry'])
            entry['word'][18] = 0
        else:
            target = next(t for t in s['targets'] if t['endpoint_index'] == example['endpoint_index'])
            value = target['minima'][example['best_new']['entry']]
            value['coordinate'][0] = str(int(value['coordinate'][0])+1)
        with tempfile.TemporaryDirectory(prefix='basis-accessibility-tamper-') as temp:
            dest = Path(temp); write(dest/'maps.json', m); s['maps_sha256'] = sha(dest/'maps.json'); write(dest/'scan.json', s)
            try:
                verify(row, plan, dest)
            except ArithmeticError as error:
                rows.append({'control': kind, 'status': 'REJECTED', 'reason': str(error)})
            else:
                raise ArithmeticError('tampered causal witness accepted: '+kind)
    result = {'status': 'PASS_NEGATIVE_REPLAYS', 'rows': rows, 'cpu_seconds': time.process_time()-start,
        'source_sha256': sha(Path(__file__)), 'point_search_calls': 0,
        'boundary': 'Original positive case already independently checked. Temporary corruptions bypass only the outer map hash, so the arithmetic/finite-dictionary checks must reject them.'}
    new_write(OUT/'regressions.json', result); print(json.dumps(result), flush=True)


def replay(output):
    from verify_cancellation_basis_accessibility_v3 import verify
    plan = audit.base.guard(); start = time.process_time(); rows = []
    for row in plan['cases']:
        result = verify(row, plan, RAW/row['id'])
        rows.append({'case': row['id'], 'status': result['status'], 'counts': result['counts']})
        print(json.dumps(rows[-1]), flush=True)
    new_write(Path(output), {'status': 'PASS_STANDALONE_REPLAY', 'rows': rows,
        'cpu_seconds': time.process_time()-start, 'point_search_calls': 0,
        'protocol_sha256': sha(OUT/'protocol.json')})


def pack():
    plan = audit.base.guard(); need(read(OUT/'regressions.json')['status'] == 'PASS_NEGATIVE_REPLAYS', 'negative controls missing')
    paths = {ROOT/name for name in plan['source_sha256']}
    paths.add(Path(__file__).resolve())
    paths.update([PRIOR/'protocol.json', PRIOR/'endpoint-packets.json'])
    for version in ('v1', 'v2', 'v3'):
        generated = ROOT/f'artifacts/generated-results/elliptic-curves/cancellation_basis_accessibility_{version}'
        local = LOCAL/f'cancellation-basis-accessibility-{version}'
        paths.update(p for p in generated.rglob('*') if p.is_file())
        paths.update(p for p in local.rglob('*') if p.is_file())
    for row in plan['cases']:
        paths.update(ROOT/item['path'] for item in row['inputs'].values())
        paths.update(ROOT/item['path'] for item in row['executed_calls'])
    bundle = OUT/'replay-bundle.tar.gz'; need(not bundle.exists(), 'preserve existing bundle')
    entries = [{'path': str(p.relative_to(ROOT)), 'sha256': sha(p), 'bytes': p.stat().st_size} for p in sorted(paths)]
    with tarfile.open(bundle, 'w:gz') as archive:
        for entry in entries:
            archive.add(ROOT/entry['path'], arcname=entry['path'], recursive=False)
    with tarfile.open(bundle, 'r:gz') as archive:
        need(len(archive.getmembers()) == len(entries), 'bundle member count differs')
        for entry in entries:
            data = archive.extractfile(entry['path']).read()
            need(digest(data) == entry['sha256'] and len(data) == entry['bytes'], 'bundle bytes differ')
    new_write(OUT/'replay-manifest.json', {'status': 'PASS_BYTE_COMPLETE_AUDIT_BUNDLE', 'entries': entries,
        'files': len(entries), 'uncompressed_bytes': sum(e['bytes'] for e in entries), 'bundle_sha256': sha(bundle),
        'boundary': 'Includes original and failed protocols, exact input packets, all finite audit maps/scans/replays, referenced executed calls and frozen source files. Byte verification adds no new mathematical replay.'})
    new_write(OUT/'completion.json', {'status': 'RETAINED_FINITE_CAUSAL_AUDIT',
        'summary_sha256': sha(OUT/'summary.json'), 'regressions_sha256': sha(OUT/'regressions.json'),
        'manifest_sha256': sha(OUT/'replay-manifest.json'), 'bundle_sha256': sha(bundle),
        'metered_successor_worker_cpu_seconds': read(OUT/'summary.json')['charged_cpu_seconds'],
        'retained_failed_worker_cpu_seconds': plan['predecessor']['charged_cpu_seconds']+plan['retained_scan']['failed_worker_cpu_seconds'],
        'supplement_internal_cpu_seconds': plan['retained_scan']['supplement_cpu_seconds'],
        'regression_internal_cpu_seconds': read(OUT/'regressions.json')['cpu_seconds'],
        'cost_boundary': 'Every isolated audit worker includes interpreter/imports and independent verification. Supplemental/manual smoke interpreter costs and packaging are outside these diagnostic component meters; no complete development-programme cost or performance advantage is claimed.',
        'point_search_calls': 0})
    print(json.dumps({'status': 'PACKED', 'files': len(entries), 'bytes': bundle.stat().st_size}), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('command', choices=['regress', 'replay', 'pack'])
    parser.add_argument('--output'); args = parser.parse_args()
    if args.command == 'replay':
        need(args.output is not None, 'supply a new --output path for standalone replay'); replay(args.output)
    else:
        globals()[args.command]()
