#!/usr/bin/env python3
"""Retain and inspect the fixed early-exposure experiment without point search.

This post-protocol supplement never changes the frozen policy or its gates.
The optional replay command rechecks one retained arm, not an entire campaign.
"""
import argparse
from collections import Counter
import json
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile
import time

from cancellation_basis_early import CAS, OUT, RAW, ROOT, guard, need, new_write, read, sha
from cancellation_scheduler_fresh import source_closure
from finite_cancellation_corpus import canonical, digest


def endpoints():
    guard()
    preflight = {r['id']: r for r in read(OUT/'preflight.json')['rows']}
    rows = []
    for oracle in read(OUT/'oracle.json'):
        path = oracle['source']; split = path.index('.json')+5
        source = ROOT/path[:split]
        need(sha(source) == preflight[oracle['id']]['source_sha256'], 'endpoint source changed')
        packet = read(source)
        for part in path[split+1:].split('/'):
            if part:
                packet = packet[int(part)] if isinstance(packet, list) else packet[part]
        proof = packet.get('rank_certificate', packet.get('proof'))
        need(digest(canonical(proof)) == oracle['proof_sha256'], 'endpoint proof changed')
        rows.append({'case': oracle['id'], 'source': path, 'source_sha256': sha(source),
            'packet': {'curve': packet.get('curve', packet.get('ainvs', packet.get('model'))),
                       'points': packet['points'], 'proof': proof},
            'independent_preflight_rank': preflight[oracle['id']]['endpoint']})
    new_write(OUT/'endpoint-packets.json', {'status': 'EXTRACTED_ALREADY_VERIFIED_ENDPOINT_INPUTS',
        'rows': rows, 'oracle_sha256': sha(OUT/'oracle.json'), 'preflight_sha256': sha(OUT/'preflight.json'),
        'source_sha256': sha(Path(__file__)), 'point_search_calls': 0,
        'boundary': 'Byte-bound extraction of the already verified eligibility packets. Their original models, points and proofs are preserved. This adds no arithmetic replay or eligibility search.'})
    print(json.dumps({'status': 'RETAINED_ENDPOINT_PACKETS', 'cases': len(rows)}), flush=True)


def bindings():
    """Read every retained input needed for a single-arm arithmetic replay."""
    plan, inputs = guard(); supervision = read(OUT/'supervision.json')
    need(supervision['status'] == 'COMPLETE' and len(supervision['records']) == 48, 'incomplete fixed roster')
    need(len({(r['case'], r['arm']) for r in supervision['records']}) == 48, 'duplicate arm')
    need({(r['case'], r['arm']) for r in supervision['records']} ==
         {(c['id'], arm) for c in inputs for arm in plan['arms']}, 'roster mismatch')
    counts = Counter(); rows = []; pairs = {}
    for receipt in supervision['records']:
        dest = RAW/'arms'/receipt['case']/receipt['arm']; result = read(dest/'result.json')
        need(receipt['status'] == 'COMPLETE' and sha(dest/'result.json') == receipt['result_sha256'], 'result changed')
        need(result['protocol_sha256'] == receipt['protocol_sha256'] == sha(OUT/'protocol.json'), 'protocol mismatch')
        need(sha(dest/'independent-verification.json') == result['independent_verification_sha256'], 'proof changed')
        need(read(dest/'independent-verification.json')['status'] == 'PASS', 'unverified arm')
        need(sha(dest/'events.json') == result['events_sha256'], 'events changed')
        need(sha(dest/'final-rank.json') == result['final_rank_sha256'], 'final packet changed')
        read(dest/'initial-rank.json')
        local = Counter(); banks = {}; first = []; gain = 0
        events = read(dest/'events.json')
        for event in events:
            if event['kind'] == 'prepared':
                need(sha(dest/event['file']) == event['sha256'], 'prepared map changed')
            elif event['kind'] == 'bank_ready':
                folder = dest/f'epoch-{event["epoch"]:02d}'
                need(sha(folder/'verification.json') == event['verification_sha256'], 'bank receipt changed')
                verified = read(folder/'verification.json')
                need(verified['status'] == 'PASS_INDEPENDENT_BANK', 'unverified bank')
                for filename, field in (('bank.json', 'bank_sha256'), ('seed.json', 'seed_sha256'),
                        ('producer/selection.json', 'producer_sha256'), ('reference/selection.json', 'reference_sha256')):
                    need(sha(folder/filename) == verified[field], 'bank input changed')
                banks[event['epoch']] = read(folder/'bank.json')
            elif event['kind'] == 'exact_completed_duplicate':
                local['exact_completed_duplicates_skipped'] += 1
        initial = result['initial_rank']
        original = {tuple(c['point']) for c in banks.get(0, {}).get('centres', [])}
        for record in result['calls']:
            need(sha(dest/record['file']) == record['sha256'], 'call changed')
            call = read(dest/record['file']); bank = banks[call['epoch']]
            need(call['word'] == bank['centres'][call['centre']]['representative'], 'anchor binding differs')
            counts['point_call_'+call['search']['status']] += 1
            if gain == 0:
                first.append({'word': call['word'], 'mapping': call['mapping'],
                    'status': call['search']['status'], 'points': call['search']['finite_curve_points']})
            gain += call['gain']
            if call['rank_file']:
                need(sha(dest/call['rank_file']) == call['rank_sha256'], 'gain certificate changed')
            if call['epoch']:
                uses_new = any(call['word'][initial:])
                local['post_refresh_calls'] += 1
                local['post_refresh_calls_using_added_generator'] += uses_new
                local['post_refresh_directions_using_added_generator'] += call['gain']*uses_new
                local['post_refresh_directions'] += call['gain']
            witness = call['tail_witness']
            if witness:
                need(witness['new_height'] <= plan['height'] < witness['old_height'], 'tail bound differs')
                local['tail_coordinate_witness_calls'] += 1
                local['tail_coordinate_witness_complete_calls'] += call['search']['status'] == 'bounded_search_complete'
        need(gain == result['new_directions'], 'gain accounting differs')
        local['verified_rebuilds'] = max(0, len(banks)-1)
        local['exported_anchors_using_added_generator'] = sum(any(c['representative'][initial:])
            for epoch, bank in banks.items() if epoch for c in bank['centres'])
        local['exported_new_anchors'] = sum(tuple(c['point']) not in original
            for epoch, bank in banks.items() if epoch for c in bank['centres'])
        pairs[receipt['case'], receipt['arm']] = {'initial_bank': banks.get(0), 'first_cloud_prefix': first}
        rows.append({'case': receipt['case'], 'arm': receipt['arm'], 'counts': dict(local)})
    paired = []
    for case in inputs:
        a = pairs[case['id'], 'basis_refresh']; b = pairs[case['id'], 'fixed_bank']
        paired.append({'case': case['id'], 'same_initial_bank': a['initial_bank'] is not None and a['initial_bank'] == b['initial_bank'],
            'same_pre_first_gain_point_clouds': a['first_cloud_prefix'] == b['first_cloud_prefix']})
    packets = read(OUT/'endpoint-packets.json')
    need(packets['oracle_sha256'] == sha(OUT/'oracle.json') and packets['preflight_sha256'] == sha(OUT/'preflight.json'), 'eligibility binding differs')
    oracle = {r['id']: r for r in read(OUT/'oracle.json')}
    need({r['case'] for r in packets['rows']} == set(oracle), 'endpoint roster differs')
    for row in packets['rows']:
        need(digest(canonical(row['packet']['proof'])) == oracle[row['case']]['proof_sha256'], 'endpoint proof differs')
    totals = {}
    for arm in plan['arms']:
        total = Counter()
        for row in rows:
            if row['arm'] == arm:
                total.update(row['counts'])
        totals[arm] = dict(total)
    return {'status': 'PASS_RETAINED_BYTES_AND_INPUT_BINDINGS', 'rows': rows, 'totals': totals,
        'paired_initial_exposure': paired, 'point_call_statuses': dict(counts),
        'protocol_sha256': sha(OUT/'protocol.json'), 'supervision_sha256': sha(OUT/'supervision.json'),
        'point_search_calls': 0,
        'boundary': 'Post-protocol receipt accounting and input-binding inspection. Arithmetic replays are the separately retained independently checked worker results. Tail witnesses assert coordinates, not rational points. This supplement changes no selection policy or promotion rule.'}


def audit():
    result = bindings(); result['source_sha256'] = sha(Path(__file__))
    new_write(OUT/'exposure-audit.json', result)
    plan, _ = guard()
    closure = source_closure([Path(__file__), CAS/'pack_cancellation_scheduler.py'])
    extra = {name: value for name, value in closure.items() if name not in plan['source_sha256']}
    for name in set(closure) & set(plan['source_sha256']):
        need(closure[name] == plan['source_sha256'][name], 'frozen source changed')
    new_write(OUT/'dependency-audit.json', {'extra_source_sha256': extra,
        'boundary': 'Post-protocol retention and replay entry points. Original sealed source hashes remain unchanged.'})
    print(json.dumps({k: result[k] for k in ('status', 'totals', 'point_call_statuses')}, indent=2), flush=True)


def portable():
    manifest = read(OUT/'replay-manifest.json'); archive = OUT/'replay-bundle.tar.gz'
    need(sha(archive) == manifest['archive_sha256'], 'bundle changed')
    with tempfile.TemporaryDirectory(prefix='basis-early-empty-root-') as temp:
        root = Path(temp); seen = set()
        with tarfile.open(archive, 'r:gz') as bundle:
            for member in bundle.getmembers():
                name = member.name
                need(name not in seen and member.isfile() and not Path(name).is_absolute() and '..' not in Path(name).parts, 'unsafe bundle member')
                raw = bundle.extractfile(member).read()
                need(digest(raw) == manifest['files_sha256'].get(name), 'bundle bytes differ')
                path = root/name; path.parent.mkdir(parents=True, exist_ok=True); path.write_bytes(raw); seen.add(name)
        need(seen == set(manifest['files_sha256']), 'missing bundle member')
        code = "import sys;sys.path.insert(0,'elliptic-curves/cas');import retain_cancellation_basis_early as r;from verify_cancellation_basis_early import verify_arm;a=r.bindings();print(a['status'],len(a['rows']))"
        process = subprocess.run([sys.executable, '-c', code], cwd=root, text=True, capture_output=True, timeout=60)
        need(process.returncode == 0, 'empty-root preflight failed: '+process.stderr)
    result = {'status': 'PASS_EMPTY_ROOT_BYTE_IMPORT_AND_INPUT_PREFLIGHT', 'files': len(seen),
        'bundle_sha256': sha(archive), 'manifest_sha256': sha(OUT/'replay-manifest.json'),
        'source_sha256': sha(Path(__file__)), 'preflight': process.stdout.strip(), 'point_search_calls': 0,
        'boundary': 'All archive bytes, source seals, imports, 48 primary-arm input bindings and 24 extracted endpoint proof bindings checked in an empty root. No primary search or full arithmetic replay is repeated.'}
    new_write(OUT/'portable-check.json', result); print(json.dumps(result), flush=True)


def replay(case_id, arm, output):
    from verify_cancellation_basis_early import verify_arm
    need(case_id and arm and output, 'replay needs one --case, --arm and unused --output path')
    need(not Path(output).exists(), 'preserve earlier replay output')
    plan, inputs = guard(); case = next(c for c in inputs if c['id'] == case_id)
    need(arm in plan['arms'], 'unknown arm')
    dest = RAW/'arms'/case_id/arm; result = read(dest/'result.json'); start = time.process_time()
    verification = verify_arm(case, arm, dest, plan, result)
    new_write(Path(output), {'case': case_id, 'arm': arm, 'verification': verification,
        'cpu_seconds': time.process_time()-start, 'point_search_calls': 0,
        'protocol_sha256': sha(OUT/'protocol.json'), 'source_sha256': sha(Path(__file__))})
    print(json.dumps({'status': verification['status'], 'case': case_id, 'arm': arm}), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=['endpoints', 'audit', 'portable', 'replay'])
    parser.add_argument('--case'); parser.add_argument('--arm'); parser.add_argument('--output')
    args = parser.parse_args()
    if args.command == 'replay':
        replay(args.case, args.arm, args.output)
    else:
        globals()[args.command]()
