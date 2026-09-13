#!/usr/bin/env python3
"""Read-only exposure analysis; optional publication of a labelled snapshot.

This reporter never searches, continues a fibre or changes the frozen run.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil
import time

from run_x948_seed_foundry import DEFAULT, ROOT, read, sha, seed_summary
from research_runtime.store import checkpoint


def acquired_after_23(result):
    rank = result.get('initial_rank', 16)
    by_call = {}
    for event in result.get('gain_timeline', []):
        by_call[event['call']] = max(by_call.get(event['call'], 0), event['after'])
    for call, after in sorted(by_call.items()):
        if after > rank and rank >= 23:
            return True
        rank = max(rank, after)
    return False


def amplification_summary(rows, parents):
    out = []
    for parent in parents:
        data = [r for r in rows if r['phase'] == 'amplify' and r['family'] == parent['family']]
        valid = [r for r in data if r['valid'] and r['result'].get('rank_lower_bound') is not None]
        gain = [r for r in valid if r['result']['rank_lower_bound'] > r['result']['initial_rank']]
        binary_complete = [r for r in valid if r in gain or r['result']['binary_outcome_complete']]
        entry = {'family': parent['family'], 'started': len(data), 'additional_gain': len(gain),
            'resolved_another_gain_outcomes': len(binary_complete),
            'complete_fixed_continuations': sum(r['result']['binary_outcome_complete'] for r in valid),
            'p_amp': len(gain)/len(binary_complete) if binary_complete else None,
            'acquiring_after_23': sum(acquired_after_23(r['result']) for r in valid),
            'CPU_seconds': sum(r['cost']['process_tree_cpu_seconds'] for r in data)}
        for threshold in range(23, 33):
            hits = sum(r['result']['rank_lower_bound'] >= threshold for r in valid)
            resolved = sum(r['result']['rank_lower_bound'] >= threshold or r['result']['binary_outcome_complete'] for r in valid)
            entry.update({f'reached_{threshold}': hits, f'resolved_{threshold}': resolved,
                          f'p_reach_{threshold}': hits/resolved if resolved else None,
                          f'already_in_seed_{threshold}': sum(r['result']['initial_rank'] >= threshold for r in valid)})
        out.append(entry)
    return out


def report(folder, state_bytes=None):
    state_bytes = state_bytes if state_bytes is not None else (folder/'state.json').read_bytes()
    state, plan, roster = json.loads(state_bytes), read(folder/'plan.json'), read(folder/'roster.json')
    rows, parents = state['exposures'], roster['parents']
    seed = seed_summary(rows, parents)
    amp = amplification_summary(rows, parents)
    baseline = next((r for r in seed if r['family'] == plan['baseline']), None)
    setup = state['process_tree_cpu_seconds']-sum(r['cost']['process_tree_cpu_seconds'] for r in rows)
    # Calibration is common detector preparation; production CPU is separately visible.
    calibration = sum(r['cost']['process_tree_cpu_seconds'] for r in rows if r['address_index'] is None)
    for row in seed:
        row['allocated_setup_CPU_seconds'] = (setup+calibration)/len(parents)
        row['cold_CPU_per_seed'] = (row['CPU_seconds']+row['allocated_setup_CPU_seconds'])/row['seeded'] if row['seeded'] else None
    candidates = []
    if (folder/'stage1-sealed.json').exists() and baseline:
        gate = plan['winner_rule']
        for row in seed:
            if row is baseline:
                continue
            complete = row['slots'] == baseline['slots'] == 128
            censor_ok = max(row['unresolved'], baseline['unresolved'])/128 <= gate['maximum_unresolved_fraction']
            enough = min(row['seeded'], baseline['seeded']) >= gate['minimum_seed_count']
            ratio = baseline['cold_CPU_per_seed']/row['cold_CPU_per_seed'] if enough else None
            candidates.append({'family': row['family'], 'seed_CPU_enrichment': ratio,
                'practical_seed_gate': bool(complete and censor_ok and enough and ratio >= gate['seed_per_total_CPU_ratio'])})
    matched = []
    if (folder/'stage2-roster.json').exists():
        scheduled = read(folder/'stage2-roster.json')['seeds']
        base_rows = sorted((r for r in rows if r['phase'] == 'amplify' and r['family'] == plan['baseline']), key=lambda r:r['address_index'])
        for parent in parents:
            family = parent['family']
            if family == plan['baseline']:
                continue
            n = min(len(scheduled[family]), len(scheduled[plan['baseline']]))
            arm = sorted((r for r in rows if r['phase'] == 'amplify' and r['family'] == family), key=lambda r:r['address_index'])[:n]
            base = base_rows[:n]
            complete = n >= plan['winner_rule']['minimum_matched_amp_seeds'] and len(arm) == len(base) == n and all(
                r['valid'] and r['result']['binary_outcome_complete'] for r in arm+base)
            ah = sum((r['result'].get('rank_lower_bound') or 0) >= 23 for r in arm)
            bh = sum((r['result'].get('rank_lower_bound') or 0) >= 23 for r in base)
            matched.append({'family': family, 'matched_prefix': n, 'all_exposures_complete': complete,
                'reached_23': ah, 'baseline_reached_23': bh,
                'practical_amp_gate': bool(complete and ah-bh >= 3 and ah >= 2*bh)})
    return {'status': state['status'], 'snapshot_unix': time.time(), 'plan_sha256': sha(folder/'plan.json'),
        'roster_sha256': sha(folder/'roster.json'), 'state_sha256': hashlib.sha256(state_bytes).hexdigest(),
        'seed_summary': seed, 'amplification_summary': amp, 'seed_gate': candidates, 'matched_amp_gate': matched,
        'charged_wall_seconds': state['charged_wall_seconds'], 'process_tree_CPU_seconds': state['process_tree_cpu_seconds'],
        'construction_and_nonexposure_CPU_seconds': setup, 'calibration_CPU_seconds': calibration,
        'limit_boundary': 'Ten-hour ceiling includes preparation, admission, calibration, failed attempts and verification. '
            'Completed-stage CPU totals exclude a still-running job until its meter closes.',
        'CPU_accounting_scope': 'Reaped worker process trees, stage supervision and initial snapshot preparation are metered. '
            'Controller hashing and ledger work between supervised stages are included in charged elapsed time '
            'but are not separately attributed to per-fibre CPU. Cold CPU ratios use this same scope in both arms; '
            'they are not complete host-energy or end-to-end speed claims.',
        'interpretation': 'Binary outcomes concern the frozen detector. No-seed is not rank16. '
            'Censored outcomes are separate. This is a finite fibration/chart/exposure comparison, '
            'not a causal decomposition of surface geometry. No winner launches further work automatically.'}


def export(folder, output):
    output.mkdir(parents=True, exist_ok=True)
    state_bytes = (folder/'state.json').read_bytes()
    result = report(folder, state_bytes)
    (output/'exposure-state.json').write_bytes(state_bytes)
    for name in ('plan.json', 'roster.json', 'proposal-order.json', 'calibration.json', 'stage1-sealed.json', 'stage2-roster.json'):
        if (folder/name).exists():
            shutil.copyfile(folder/name, output/name)
    roster = read(folder/'roster.json')
    for parent in roster['parents']:
        if parent['baseline']:
            continue
        directory = Path(parent['admission']).parent
        for name in ('parent.json','admission.json','verified.json'):
            dest = output/directory/name
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(folder/directory/name, dest)
    checkpoint(output/'summary.json', result)
    for name, records in [('seed-summary.csv',result['seed_summary']), ('amplification-summary.csv',result['amplification_summary'])]:
        if records:
            with (output/name).open('w', newline='') as stream:
                writer = csv.DictWriter(stream, fieldnames=list(records[0]))
                writer.writeheader()
                writer.writerows(records)
    checkpoint(output/'manifest.json', {'schema':'x948-foundry.snapshot.v1',
        'local_runtime':str(folder), 'snapshot_unix': result['snapshot_unix'],
        'files': {str(p.relative_to(output)): sha(p) for p in output.rglob('*') if p.is_file() and p.name != 'manifest.json'},
        'boundary':'Point-search status is a dated snapshot. Live process state is checked separately; admission packets are exact retained mathematics.'})
    print(output/'summary.json')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--folder',type=Path,default=DEFAULT)
    parser.add_argument('--export',type=Path)
    args = parser.parse_args()
    if args.export:
        export(args.folder.resolve(), args.export.resolve())
    else:
        print(json.dumps(report(args.folder.resolve()), indent=2))
