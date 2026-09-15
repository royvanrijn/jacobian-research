#!/usr/bin/env python3
"""Post-protocol receipt diagnostic: did the search reach enlarged-basis anchors?

No new searches, fitting, or selection changes. The distinction between a
verified new bank and actually exposing its new directions is descriptive.
"""
from collections import Counter
import json
from pathlib import Path

from finite_cancellation_corpus import ROOT, canonical, digest, write
from cancellation_basis_amplification import OUT, RAW, guard


def read(path): return json.loads(path.read_text())
def sha(path): return digest(path.read_bytes())
def need(condition, message):
    if not condition:
        raise ArithmeticError(message)


def audit():
    from cancellation_basis_finish import effective
    plan, inputs = guard(); cases = {c['id']: c for c in inputs}
    supervision = read(OUT/'supervision.json')
    need(supervision['status'] == 'COMPLETE' and len(supervision['records']) == 2*plan['cases'], 'incomplete fixed cohort')
    # Extract only the already checked endpoint packets, preserving proof
    # bytes and provenance. Replaying eligibility then needs no large parent
    # campaign file and never reconstructs a missing local artifact.
    preflight = {r['id']: r for r in read(OUT/'preflight.json')['rows']}
    endpoints = []
    for oracle in read(OUT/'oracle.json'):
        path = oracle['source']; split = path.index('.json')+5
        source = ROOT/path[:split]
        need(sha(source) == preflight[oracle['id']]['source_sha256'], 'preflight endpoint source changed')
        packet = read(source)
        for part in path[split+1:].split('/'):
            if part: packet = packet[int(part)] if isinstance(packet, list) else packet[part]
        proof = packet.get('rank_certificate', packet.get('proof'))
        need(digest(canonical(proof)) == oracle['proof_sha256'], 'endpoint proof changed')
        endpoints.append({'case': oracle['id'], 'source': path, 'source_sha256': sha(source),
            'packet': {'curve': packet.get('curve', packet.get('ainvs', packet.get('model'))),
                       'points': packet['points'], 'proof': proof},
            'independent_preflight_rank': preflight[oracle['id']]['endpoint']})
    write(OUT/'endpoint-packets.json', {'status': 'EXTRACTED_ALREADY_VERIFIED_ENDPOINT_INPUTS', 'rows': endpoints,
        'oracle_sha256': sha(OUT/'oracle.json'), 'preflight_sha256': sha(OUT/'preflight.json'),
        'boundary': 'Byte-bound extraction of previously independently checked endpoint packets; no new search, eligibility selection, or arithmetic replay.'})
    rows = []; paired = {}; statuses = Counter()
    for original_receipt in supervision['records']:
        dest, result, receipt = effective(original_receipt)
        initial = result['initial_rank']; banks = {}
        for epoch in result['epochs']:
            folder = dest/f'epoch-{epoch["epoch"]:02d}'
            need(sha(folder/'verification.json') == epoch['verification_sha256'], 'bank receipt changed')
            need(sha(folder/'bank.json') == epoch['bank_sha256'], 'bank changed')
            banks[epoch['epoch']] = read(folder/'bank.json')
        initial_anchors = {tuple(c['point']) for c in banks.get(0, {}).get('centres', [])}
        counts = Counter(); gains = 0; first_prefix = []
        for record in result['calls']:
            path = dest/record['file']; need(sha(path) == record['sha256'], 'call changed')
            call = read(path); bank = banks[call['epoch']]
            need(call['word'] == bank['centres'][call['centre']]['representative'], 'wrong bank anchor')
            statuses[call['search']['status']] += 1
            if gains == 0:
                first_prefix.append({'word': call['word'], 'mapping': call['mapping'],
                                     'status': call['search']['status'], 'points': call['search']['finite_curve_points']})
            gains += call['gain']
            if call['rank_file']:
                need(sha(dest/call['rank_file']) == call['rank_sha256'], 'cloud rank packet changed')
            if call['epoch'] > 0:
                counts['post_refresh_calls'] += 1
                counts['post_refresh_directions'] += call['gain']
                uses_added = any(call['word'][initial:])
                counts['post_refresh_calls_using_added_generator'] += uses_added
                counts['post_refresh_directions_using_added_generator'] += call['gain']*uses_added
                new_anchor = tuple(bank['centres'][call['centre']]['point']) not in initial_anchors
                counts['post_refresh_calls_on_new_anchors'] += new_anchor
                counts['post_refresh_directions_on_new_anchors'] += call['gain']*new_anchor
            witness = call['tail_witness']
            if witness:
                need(witness['new_height'] <= plan['height'] < witness['old_height'], 'tail height witness differs')
                counts['tail_coordinate_witness_calls'] += 1
                counts['tail_coordinate_witness_complete_calls'] += call['search']['status'] == 'bounded_search_complete'
        counts['verified_rebuilds'] = max(0, len(banks)-1)
        counts['exported_anchors_using_added_generator'] = sum(
            any(c['representative'][initial:]) for epoch, bank in banks.items() if epoch > 0 for c in bank['centres'])
        counts['exported_new_anchors'] = sum(tuple(c['point']) not in initial_anchors
            for epoch, bank in banks.items() if epoch > 0 for c in bank['centres'])
        counts['rebuilds_exporting_added_generator_anchors'] = sum(
            any(any(c['representative'][initial:]) for c in bank['centres']) for epoch, bank in banks.items() if epoch > 0)
        counts['exact_completed_duplicates_skipped'] = sum(e['kind'] == 'exact_completed_duplicate' for e in read(dest/'events.json'))
        need(gains == result['new_directions'], 'uncapped gain count differs')
        paired[receipt['case'], receipt['arm']] = {'initial_bank': banks.get(0), 'first_prefix': first_prefix}
        rows.append({'case': receipt['case'], 'arm': receipt['arm'], **dict(counts)})
    pairs = []
    for case in inputs:
        a = paired[case['id'], 'basis_refresh']; b = paired[case['id'], 'fixed_bank']
        pairs.append({'case': case['id'], 'same_initial_bank': a['initial_bank'] is not None and a['initial_bank'] == b['initial_bank'],
                      'same_pre_first_gain_point_clouds': a['first_prefix'] == b['first_prefix']})
    totals = {arm: dict(sum((Counter({k: v for k, v in row.items() if k not in ('case', 'arm')})
                           for row in rows if row['arm'] == arm), Counter())) for arm in plan['arms']}
    # Counter addition omits zeros; keep critical zero observations explicit.
    for total in totals.values():
        for name in ('verified_rebuilds', 'post_refresh_calls', 'post_refresh_directions',
                     'post_refresh_calls_using_added_generator', 'post_refresh_directions_using_added_generator',
                     'exported_anchors_using_added_generator', 'rebuilds_exporting_added_generator_anchors',
                     'post_refresh_calls_on_new_anchors', 'post_refresh_directions_on_new_anchors', 'exported_new_anchors'):
            total.setdefault(name, 0)
    result = {'status': 'PASS_POST_PROTOCOL_MECHANISM_DIAGNOSTIC', 'rows': rows, 'totals': totals,
        'paired_initial_exposure': pairs, 'point_call_statuses': dict(statuses),
        'protocol_sha256': sha(OUT/'protocol.json'), 'supervision_sha256': sha(OUT/'supervision.json'),
        'source_sha256': sha(Path(__file__)), 'point_search_calls_added': 0,
        'boundary': 'Descriptive receipt audit added after protocol sealing. Actual nonzero coefficients in newly admitted independent generators distinguish use of the enlarged subgroup from merely rebuilding a compatible bank. This diagnostic changes no policy or promotion rule. Coordinate witnesses do not assert rational points.'}
    write(OUT/'mechanism-audit.json', result)
    supplement = [Path(__file__), Path(__file__).with_name('check_cancellation_basis_transitions.py'),
                  Path(__file__).with_name('cancellation_basis_finish.py')]
    write(OUT/'dependency-audit.json', {'extra_source_sha256': {str(p.relative_to(ROOT)): sha(p) for p in supplement},
        'boundary': 'Post-protocol descriptive receipt checker, not retroactively included in the pre-execution source seal.'})
    print(json.dumps({'status': result['status'], 'totals': totals, 'point_call_statuses': dict(statuses),
                      'same_initial_banks': sum(p['same_initial_bank'] for p in pairs),
                      'same_first_cloud_prefixes': sum(p['same_pre_first_gain_point_clouds'] for p in pairs)}, indent=2), flush=True)


if __name__ == '__main__':
    audit()
