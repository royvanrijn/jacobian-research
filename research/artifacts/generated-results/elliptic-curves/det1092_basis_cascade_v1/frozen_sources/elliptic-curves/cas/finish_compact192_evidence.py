#!/usr/bin/env python3
"""Finish inventory and isolated evidence after the frozen192 cohort closes."""
import argparse
import sys
import time
from collections import Counter
from pathlib import Path
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run, Limits

ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / 'elliptic-curves/cas'
ART = ROOT / 'artifacts/generated-results/elliptic-curves'
LOCAL = ROOT / 'artifacts/local/elliptic-curves'
D = LOCAL / 'compact192-evidence-finalization-v1'


def prepare():
    if (D/'protocol.json').exists():
        raise FileExistsError('preserve compact192 evidence finalization')
    names = ['finish_compact192_evidence.py', 'export_new_high_rank_curve_index_v13.py',
             'replay_inventory_v13_memory.py', 'package_compact192_results.py',
             'verify_compact192_bundle.py', 'report_compact192_portable.py',
             'report_r17_13_scaling_geometry_v2.py']
    paths = [*(CAS/n for n in names),
             LOCAL/'compact192-r17-pari-v1/post-batch/protocol.json',
             ART/'new_high_rank_curve_index_v12.json',
             ART/'endpoint_point_trial_evidence_v1.json']
    checkpoint(D/'protocol.json', {
        'schema': 'elliptic-curves.compact192-evidence-finalization.v1',
        'sources': {str(p.relative_to(ROOT)): cert.hashed(p) for p in paths},
        'wait_deadline_seconds': 120000, 'inventory_build_seconds': 600,
        'inventory_replay_seconds': 1200, 'package_seconds': 1800,
        'isolated_replay_seconds': 60000, 'report_seconds': 180,
        'rss_bytes': 4294967296,
        'scope': 'Wait for the existing fixed192 finalizer to pass every point, geometry, catalogue and odd-cloud stage. Any failure, censoring or stronger odd-prime bound stops this process for review. Export only exactly certified catalogue-unmatched curves of lower bound at least22, preserving all V12 IDs; independently replay the full V13 inventory and CSV. Render the completed finite-experiment note, package immutable evidence and perform exactly585 isolated replay stages, with one verifier at a time. No point search, new candidate, adaptive wave, automatic retry or mathematical-status promotion occurs here.'})


def write_note():
    trial = cert.read(ART/'compact192_experiment_v1.json')
    inventory = cert.read(ART/'new_high_rank_curve_index_v13.json')
    old = cert.read(ART/'new_high_rank_curve_index_v12.json')
    counts = dict(sorted(Counter(r['rank_lower_bound'] for r in inventory['curves']).items(), reverse=True))
    text = [
        '# Completed fixed192 compact-parameter point experiment', '',
        f"All192 candidates and their exact point histories, geometry and full-cloud proofs pass. The fixed8640-chart cap gives {trial['attempted_point_boxes']} attempted charts and {trial['completed_point_boxes']} completed point boxes. There are {trial['catalogue_unmatched_candidates_at_least22']} catalogue-unmatched candidates with certified lower bound at least22.", '',
        f"The [V13 inventory](../../artifacts/generated-results/elliptic-curves/new_high_rank_curve_index_v13.json) contains {len(inventory['curves'])} distinct curves, an increase of {len(inventory['curves'])-len(old['curves'])} over V12. Its exact point proofs and equation CSV pass independent local replay. Lower-bound counts are {counts}. Mathematical-status promotion and completion of the separate portable replay must be read from MATH_STATUS.json and the replay report.", '',
        'The [experiment certificate](../../artifacts/generated-results/elliptic-curves/compact192_experiment_v1.json) binds the unchanged selection, exact point proofs and pinned catalogue comparisons. The [inventory replay](../../artifacts/generated-results/elliptic-curves/new_high_rank_curve_index_v13_memory_replay_v1.json) rechecks every retained independent point certificate and the equation CSV.', '',
        '## Fixed selection and exposure', '',
        'The roster uses32 previously unsearched candidates from each of six R17 families, selected from6144 already scored H4096 rows. Exact rational-isomorphism exclusions cover528 previous measured equations and21 endpoints before selection. Family order is fixed; S1 through32749, good-prime count, denominator and signed numerator determine ordering. Validation32771through65521 does not enter ties. No new parameter or prime-score scan is part of this trial.', '',
        'All192 map files precede point execution. Each curve starts with17 generic sections and uses all43 or49 exact generic maximum-parity classes, height125000 and ten seconds per chart. Each can stop at a provisional28-point subgroup pending exact replay. There is no adaptive wave or refill. Catalogue comparison starts only after the complete batch and streaming point proofs pass.', '',
        '## Strongest cohort certificates', '',
        '| Candidate | Parameter | Certified lower bound | Pinned catalogue matches |',
        '| --- | --- | ---: | --- |',
    ]
    for row in trial['rows']:
        if row['rank_lower_bound'] >= 26:
            text.append(f"| {row['id']} | {row['parameter']} | {row['rank_lower_bound']} | {row['catalogue_matches']} |")
    text += ['', '## Scope and replay', '',
             'The catalogue comparison uses593 pinned ICARM equations. Absence from that snapshot is not universal novelty. These point certificates prove lower bounds, not exact whole-curve ranks, rank upper bounds, saturation, point absence or optimality of the selector.', '',
             'The [first26 minimal-model export](../../artifacts/generated-results/elliptic-curves/compact192_first26_candidate_v1.json) is retained with its original pre-comparison status. The completed cohort certificate supplies the subsequent catalogue result; the earlier artifact is not rewritten.', '',
             'The [evidence supplement](../../artifacts/generated-results/elliptic-curves/compact192_evidence_v1.json) supports585 isolated checks, including the [complete13-scaling classification](R17_INTEGRAL_13_PARAMETER_CHARTS_2026-09-06.md). The [portable replay report](../../artifacts/generated-results/elliptic-curves/compact192_portable_replay_v1.json), when present, records whether all585 passed. No new point search runs during replay.', '']
    (ROOT/'elliptic-curves/notes/COMPACT192_UNSEARCHED_TRIAL_2026-09-06.md').write_text('\n'.join(text))


def launch():
    p = cert.read(D/'protocol.json')
    if any(cert.hashed(ROOT/n) != h for n, h in p['sources'].items()):
        raise ArithmeticError('frozen evidence sources changed')
    out = D/'ledger.json'
    if out.exists():
        raise FileExistsError('preserve single evidence finalization')
    ledger = {'status': 'WAITING_FOR_FIXED_JOBS', 'rows': []}
    checkpoint(out, ledger)
    deadline = time.monotonic() + p['wait_deadline_seconds']
    while True:
        upstream = cert.read(LOCAL/'compact192-r17-pari-v1/post-batch/ledger.json')
        state = upstream['status']
        if any(row['status'] != 'PASS' for row in upstream['rows']):
            raise ArithmeticError('upstream failed or censored stage requires review')
        for name, allowed in [('ledger.json', ('RUNNING', 'COMPLETE_FIXED_BATCH_ATTEMPTS')),
                              ('stream-verification-v1/ledger.json', ('RUNNING', 'PASS'))]:
            if cert.read(LOCAL/'compact192-r17-pari-v1'/name)['status'] not in allowed:
                raise ArithmeticError('upstream point or proof failure requires review')
        if state == 'PASS':
            break
        if state not in ('WAITING_FOR_FIXED_JOBS', 'RUNNING'):
            raise ArithmeticError('upstream failure or stronger odd bound requires review')
        if time.monotonic() > deadline:
            raise TimeoutError('fixed upstream wait deadline')
        time.sleep(5)
    ledger['status'] = 'RUNNING'
    checkpoint(out, ledger)
    jobs = [
        ('inventory-build', 'export_new_high_rank_curve_index_v13.py', [], p['inventory_build_seconds']),
        ('inventory-replay', 'replay_inventory_v13_memory.py', ['--output', str(ART/'new_high_rank_curve_index_v13_memory_replay_v1.json')], p['inventory_replay_seconds']),
        ('package', 'package_compact192_results.py', [], p['package_seconds']),
        ('isolated-replay', 'verify_compact192_bundle.py', [], p['isolated_replay_seconds']),
        ('portable-report', 'report_compact192_portable.py', [], p['report_seconds']),
    ]
    for name, script, args, seconds in jobs:
        if name == 'package':
            write_note()
        path = D/(name+'.supervisor.json')
        if path.exists():
            raise FileExistsError('preserve one evidence stage')
        s = run([sys.executable, str(CAS/script), *args], limits=Limits(seconds, p['rss_bytes']),
                log_path=D/(name+'.log'), checkpoint_path=path, cwd=ROOT)
        ok = s['outcome'] == 'completed' and s['returncode'] == 0
        ledger['rows'].append({'name': name, 'status': 'PASS' if ok else 'FAILED_OR_CENSORED', 'supervision': s})
        checkpoint(out, ledger)
        print('COMPACT192 EVIDENCE', name, ledger['rows'][-1]['status'], flush=True)
        if not ok:
            raise ArithmeticError('evidence stage failed or censored; no retry')
    ledger['status'] = 'PASS'
    checkpoint(out, ledger)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('stage', choices=['prepare', 'launch'])
    stage = parser.parse_args().stage
    try:
        globals()[stage]()
    except Exception as error:
        if stage == 'launch' and (D/'ledger.json').exists():
            record = cert.read(D/'ledger.json')
            record.update(status='FAILED_OR_CENSORED', reason=str(error))
            checkpoint(D/'ledger.json', record)
        raise
