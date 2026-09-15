#!/usr/bin/env python3
"""Retain interrupted bank work and certify its completed point-search prefix.

This supplement never restarts a point worker, completes an interrupted CVP,
changes a primary receipt, or promotes a failed protocol. It charges a separate
isolated final verification to the interrupted arm and reports all frozen cases.
"""
import argparse
from collections import Counter
from copy import deepcopy
import json
from pathlib import Path
import resource
import shutil
import signal
import subprocess
import tempfile
import time

from cancellation_basis_amplification import OUT, RAW, guard
from finite_cancellation_corpus import ROOT, canonical, digest, write


def read(path): return json.loads(path.read_text())
def sha(path): return digest(path.read_bytes())
def need(condition, message):
    if not condition: raise ArithmeticError(message)


def replay_interrupted(case_id, arm):
    from cancellation_basis_epoch import native_rank
    from verify_cancellation_basis_amplification import verify_arm
    plan, cases = guard(); case = next(c for c in cases if c['id'] == case_id)
    dest = RAW/'arms'/case_id/arm; recipe = read(dest/'recovery-protocol.json')
    need(recipe['source_sha256'] == sha(Path(__file__)), 'supplementary source changed')
    need(recipe['events_sha256'] == sha(dest/'events.json'), 'interrupted events changed')
    resource.setrlimit(resource.RLIMIT_CPU, (recipe['hard_cpu_seconds'],)*2)
    events = read(dest/'events.json')
    need(events and events[-1]['kind'] == 'bank_start', 'failure is not an interrupted bank with a sealed point prefix')
    need('AlarmInterrupt' in (dest/'worker.log').read_text(), 'unexpected failure type')
    need(not (dest/'result.json').exists(), 'do not replace a primary result')
    prefix = events[:-1]; calls = [{k: e[k] for k in ('file', 'sha256', 'epoch', 'gain', 'status', 'tail_witness')}
                                  for e in prefix if e['kind'] == 'call']
    epochs = [{k: e[k] for k in ('epoch', 'rank', 'verification_sha256', 'bank_sha256', 'cpu_seconds')}
              for e in prefix if e['kind'] == 'bank_ready']
    packet = case['seed']; first = 0; gains = 0; last_cpu = None
    for record in calls:
        need(sha(dest/record['file']) == record['sha256'], 'point receipt changed')
        row = read(dest/record['file'])
        if row['gain']:
            need(sha(dest/row['rank_file']) == row['rank_sha256'], 'rank receipt changed')
            packet = read(dest/row['rank_file'])['packet']
            if not gains: first = row['gain']; last_cpu = row['cpu_after_certification']
        gains += row['gain']
    need(len(packet['points']) == len(case['seed']['points'])+gains, 'prefix gain count differs')
    result = {'case': case_id, 'family': case['family'], 'stratum': case['stratum'], 'arm': arm,
        'initial_rank': len(case['seed']['points']), 'rank_lower_bound': len(packet['points']),
        'new_directions': gains, 'first_cloud_gain': first, 'first_gain_cpu_seconds': last_cpu,
        'later_cloud_directions': gains-first, 'success': gains >= plan['target_directions'],
        'calls': calls, 'epochs': epochs, 'unknowns': sum(e['kind'] in ('bank_unknown', 'preparation_unknown') for e in prefix),
        'status': 'CPU_CAP_BEFORE_REFRESH' if calls else 'CPU_CAP_BEFORE_BANK',
        'protocol_sha256': sha(OUT/'protocol.json')}
    # Verify a literal executed prefix ending just before the failed bank.
    # The original log retains the attempted bank and its complete CPU cost.
    # This proof view cannot assert that the interrupted bank was verified.
    with tempfile.TemporaryDirectory(prefix='basis-interrupted-prefix-') as tmp:
        view = Path(tmp)/'arm'; shutil.copytree(dest, view)
        write(view/'events.json', prefix)
        write(view/'final-rank.json', {'packet': packet, 'independent_rank': native_rank(packet)})
        result.update(events_sha256=sha(view/'events.json'), final_rank_sha256=sha(view/'final-rank.json'))
        verification = verify_arm(case, arm, view, plan, result)
        write(dest/'point-prefix-proof-view.json', {'events': prefix, 'result': result,
            'omitted_incomplete_work': events[-1],
            'boundary': 'Executed point prefix only. The following bank attempt remains interrupted, unverified and fully charged; no point call followed it.'})
        write(dest/'interrupted-final-rank.json', read(view/'final-rank.json'))
    verification.update(interrupted_bank=events[-1], incomplete_banks=1,
                        original_events_sha256=sha(dest/'events.json'), point_search_calls_added=0,
                        supplementary_source_sha256=sha(Path(__file__)))
    write(dest/'interrupted-verification.json', verification)
    result.update(status='CERTIFIED_POINT_PREFIX_BANK_REPLAY_INTERRUPTED', success=False,
        unknowns=result['unknowns']+1, events_sha256=sha(dest/'events.json'),
        final_rank_sha256=sha(dest/'interrupted-final-rank.json'),
        independent_verification_sha256=sha(dest/'interrupted-verification.json'),
        components={'bank': sum(e['cpu_seconds'] for e in epochs)},
        incomplete_bank_cpu_seconds='INCLUDED_IN_OUTER_CPU_NOT_SEPARATELY_MEASURED',
        boundary='The sealed point prefix is independently certified. Failed rational-CVP work is retained, charged, and never used for search. The original worker failure remains a failed gate; this is not a successful bank rebuild or timing retry.')
    write(dest/'interrupted-result.json', result)
    print(json.dumps({'case': case_id, 'arm': arm, 'status': result['status'], 'new_directions': gains, 'calls': len(calls)}), flush=True)


def recover():
    plan, _ = guard(); supervision = read(OUT/'supervision.json')
    need(len(supervision['records']) == 2*plan['cases'], 'finish all still-unstarted arms first')
    for receipt in supervision['records']:
        if receipt['status'] == 'COMPLETE': continue
        dest = RAW/'arms'/receipt['case']/receipt['arm']
        if (dest/'recovery-supervisor.json').exists():
            need(read(dest/'recovery-supervisor.json')['status'] == 'COMPLETE_POINT_PREFIX_REPLAY', 'retain failed supplementary replay')
            continue
        need(not (dest/'recovery-protocol.json').exists(), 'preserve unreceipted recovery')
        budget = min(30, int(plan['hard_process_cpu_seconds']-receipt['charged_cpu_seconds']))
        need(budget > 0, 'no verification budget remains')
        write(dest/'recovery-protocol.json', {'status': 'FROZEN_POINT_PREFIX_VERIFICATION_ONLY',
            'source_sha256': sha(Path(__file__)), 'events_sha256': sha(dest/'events.json'),
            'original_supervisor_sha256': sha(dest/'supervisor.json'), 'protocol_sha256': sha(OUT/'protocol.json'),
            'hard_cpu_seconds': budget, 'wall_seconds': 60, 'point_search_calls': 0,
            'boundary': 'No search or CVP retry. The unverified bank remains failed; verify only the previously sealed point prefix and charge this entire child process to that arm.'})
        command = ['sage', '-python', str(Path(__file__).resolve()), 'replay-interrupted', '--case', receipt['case'], '--arm', receipt['arm']]
        before = resource.getrusage(resource.RUSAGE_CHILDREN); tick = time.monotonic()
        with (dest/'recovery.log').open('w') as log:
            process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
            try: code = process.wait(timeout=60)
            except subprocess.TimeoutExpired:
                import os
                os.killpg(process.pid, signal.SIGKILL); code = process.wait()
        after = resource.getrusage(resource.RUSAGE_CHILDREN)
        charged = after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime
        success = code == 0 and (dest/'interrupted-result.json').exists()
        record = {'status': 'COMPLETE_POINT_PREFIX_REPLAY' if success else 'UNKNOWN_FAILED_PREFIX_REPLAY',
            'returncode': code, 'charged_cpu_seconds': charged, 'wall_seconds': time.monotonic()-tick,
            'original_supervisor_sha256': sha(dest/'supervisor.json'), 'point_search_calls': 0,
            'result_sha256': sha(dest/'interrupted-result.json') if success else None}
        write(dest/'recovery-supervisor.json', record)
        need(success, 'point-prefix replay failed; do not rerun')
        print(json.dumps({'case': receipt['case'], 'arm': receipt['arm'], **record}), flush=True)


def effective(receipt):
    """Read original outcomes plus explicitly charged interrupted-prefix proofs."""
    dest = RAW/'arms'/receipt['case']/receipt['arm']; accounting = dict(receipt)
    if receipt['status'] == 'COMPLETE':
        filename = 'result.json'; verification = 'independent-verification.json'
    else:
        recovery = read(dest/'recovery-supervisor.json')
        need(recovery['status'] == 'COMPLETE_POINT_PREFIX_REPLAY', 'uncertified interrupted prefix')
        need(sha(dest/'supervisor.json') == recovery['original_supervisor_sha256'], 'failed receipt changed')
        filename = 'interrupted-result.json'; verification = 'interrupted-verification.json'
        accounting.update(result_sha256=recovery['result_sha256'],
            charged_cpu_seconds=receipt['charged_cpu_seconds']+recovery['charged_cpu_seconds'],
            supplementary_verification_cpu_seconds=recovery['charged_cpu_seconds'])
    accounting['result_file'] = filename; accounting['verification_file'] = verification
    result = read(dest/filename)
    need(sha(dest/filename) == accounting['result_sha256'], 'result receipt differs')
    need(sha(dest/verification) == result['independent_verification_sha256'], 'verification bytes differ')
    proof = read(dest/verification)
    need(proof['status'] == 'PASS' and proof['rank']['rank'] == result['rank_lower_bound'], 'uncertified point count')
    need(sha(dest/'events.json') == result['events_sha256'], 'events differ')
    for row in result['calls']: need(sha(dest/row['file']) == row['sha256'], 'call changed')
    return dest, result, accounting


def report():
    from cancellation_cloud_programme import compare
    plan, _ = guard(); supervision = read(OUT/'supervision.json')
    need(len(supervision['records']) == 2*plan['cases'], 'missing frozen arm')
    rows = []; receipts = []
    for original in supervision['records']:
        dest, r, receipt = effective(original); receipts.append(receipt)
        rows.append({**{k: r[k] for k in ('case', 'family', 'stratum', 'arm', 'success', 'new_directions',
            'later_cloud_directions', 'first_cloud_gain', 'initial_rank', 'rank_lower_bound', 'status', 'unknowns')},
            'calls': len(r['calls']), 'epochs': len(r['epochs']), 'cpu_seconds': receipt['charged_cpu_seconds'],
            'verified_bank_cpu_seconds': r['components']['bank'],
            'primary_process_status': original['status'],
            'supplementary_verification_cpu_seconds': receipt.get('supplementary_verification_cpu_seconds', 0),
            'tail_witness_calls': sum(c['tail_witness'] for c in r['calls'])})
    def comparison(later):
        converted = [{**r, 'arm': 'factor_free' if r['arm'] == 'fixed_bank' else 'adaptive_cloud',
            'new_directions': r['later_cloud_directions'] if later else r['new_directions']} for r in rows]
        totals, result = compare(converted)
        return {('fixed_bank' if k == 'factor_free' else 'basis_refresh'): v for k, v in totals.items()}, result
    totals, overall = comparison(False); later_totals, later = comparison(True)
    unknowns = sum(r['unknowns'] for r in rows); failures = sum(r['status'] != 'COMPLETE' for r in supervision['records'])
    gate = later['gate_passed']
    if unknowns or failures or totals['basis_refresh']['directions'] < totals['fixed_bank']['directions']:
        gate = False
    shared = read(OUT/'preflight.json')['cpu_seconds']; a, b = later_totals['basis_refresh'], later_totals['fixed_bank']
    conservative = (a['directions']/(a['cpu_seconds']+shared))/(b['directions']/b['cpu_seconds']) if b['directions'] else None
    summary = {'status': 'COMPLETE_FIXED_COHORT_WITH_RETAINED_INTERRUPTION', 'rows': rows, 'totals': totals,
        'all_direction_comparison': overall, 'later_cloud_totals': later_totals, 'later_cloud_comparison': later,
        'promotion_gate': gate, 'preparation_unknowns': unknowns, 'interrupted_primary_arms': failures,
        'shared_preflight_cpu_seconds': shared, 'conservative_later_rate_ratio': conservative,
        'protocol_sha256': sha(OUT/'protocol.json'), 'rank32': 'UNKNOWN', 'fresh_fibres_run': 0,
        'supplementary_source_sha256': sha(Path(__file__)),
        'boundary': 'All frozen arms attempted once. Interrupted bank work stays failed and fully charged. Independently replayed captured point prefixes contribute only certified gains; no interrupted bank is used, completed, or retried. Supplementary verification CPU is added to its arm. The protocol fails promotion irrespective of rate estimates.'}
    write(OUT/'completion.json', {'status': summary['status'], 'accounting': receipts,
        'original_supervision_sha256': sha(OUT/'supervision.json'), 'protocol_sha256': sha(OUT/'protocol.json'),
        'boundary': summary['boundary']})
    write(OUT/'summary.json', summary)
    write(OUT/'audit.json', {'status': 'PASS_CERTIFIED_OUTCOMES_WITH_RETAINED_INTERRUPTION',
        'arms': len(rows), 'interrupted_primary_arms': failures, 'new_directions': sum(r['new_directions'] for r in rows),
        'summary_sha256': sha(OUT/'summary.json'), 'completion_sha256': sha(OUT/'completion.json'),
        'point_search_calls_added': 0, 'supplementary_source_sha256': sha(Path(__file__))})
    print(json.dumps({k: summary[k] for k in ('totals', 'later_cloud_totals', 'later_cloud_comparison', 'promotion_gate', 'interrupted_primary_arms')}, indent=2), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=['recover', 'replay-interrupted', 'report'])
    parser.add_argument('--case'); parser.add_argument('--arm'); args = parser.parse_args()
    replay_interrupted(args.case, args.arm) if args.command == 'replay-interrupted' else globals()[args.command]()
