#!/usr/bin/env python3
"""Freeze and supervise exactly two bounded equation-only relation workers."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import json
import os
import signal
import subprocess
import time

import historical_external_arithmetic as historical
import arithmetic_frozen_snapshot as frozen
import two_class_relation_core as core
import wide_arithmetic_profile_core as storage

ROOT = historical.ROOT
OUT = ROOT/'artifacts/generated-results/elliptic-curves/two_class_relation_pilot_v1'
WORKER = Path(__file__).with_name('two_class_relation_worker.sage')
read, sha, save, require = historical.read, historical.sha, historical.save, historical.require


def prepare(out):
    require(not (out/'plan.json').exists(), 'pilot already frozen')
    old = historical.DEFAULT_OUTPUT
    history, source_resolution = frozen.verify()
    census = Path(history['census'])
    protected = {str(p): historical.tree_hashes(p) for p in (old, census)}
    prospective = read(old/'prospective_profiles.json')['rows']
    external = read(old/'historical_profiles.json')['rows']
    selected = [(next(r for r in prospective if r['family']=='11952' and r['t']=='921/653'), census),
                (next(r for r in external if r['family']=='11952' and r['t']=='110314/102227'), old)]
    protocol = {'factor_base_bound': 20000, 'candidate_budget': 16384,
                'candidate_prefixes': [128, 1024, 8192, 16384],
                'height_reduction_max_steps': 128, 'wall_seconds_per_field': 180, 'max_jobs': 2,
                'factor_base_rule': 'all prime ideals above rational p<=20000 plus all denominator-prime ideals',
                'relation_rule': 'all canonical (p) relations, then primitive (a,b), b>0, ordered by max(|a|,b)',
                'comparison': 'fixed attempted-candidate prefixes; report achieved relation counts, not an assumed N',
                'classes_or_points_supplied': False, 'new_point_searches': False,
                'large_prime_variant': False, 'full_class_group': False,
                'interpretation': 'factor-base image bound only; full g and generation UNKNOWN'}
    runtime = read(census/'base-runtime.json')
    sources = {historical.worker.rel(p): sha(p) for p in (Path(__file__), WORKER, Path(core.__file__), Path(storage.__file__), Path(frozen.__file__))}
    rows = []
    out.mkdir(parents=True, exist_ok=True)
    (out/'inputs').mkdir()
    for profile, directory in selected:
        key = profile['curve_key']
        b, l = (read(directory/mode/f'{key}.json') for mode in ('base', 'local'))
        require(b['status']==l['status']=='PASS' and b['rational_2torsion_rank']==0, 'target arithmetic not certified')
        a1, a2, a3, a4, a6 = map(int, b['minimal_ainvs'])
        cubic = [16*(a3*a3+4*a6), 8*(2*a4+a1*a3), a1*a1+4*a2, 1]
        inp = {'schema': 'elliptic-curves.two-class-relation-input.v1', 'curve_key': key,
               'integral_monic_cubic_ascending': list(map(str, cubic)),
               'field_discriminant': l['field_discriminant'], 'field_signature': l['field_signature'],
               'certified_bad_rational_primes': [v[0] for v in l['minimal_discriminant_factorization']],
               'protocol': protocol}
        save(out/'inputs'/f'{key}.json', inp)
        (out/key).mkdir()
        rows.append({k: profile[k] for k in ('curve_key', 'family', 't', 'final_rank_lower_bound', 'cohort', 'selection_mode')})
        rows[-1].update(bk_local_term=l['bk_local_term'],
                       forced_g_lower_from_known_rank=profile['final_rank_lower_bound']-l['bk_local_term'])
    require(all(historical.tree_hashes(Path(p))==hashes for p, hashes in protected.items()), 'old evidence changed')
    save(out/'plan.json', {'schema': 'elliptic-curves.two-class-relation-pilot-plan.v1',
         'protocol': protocol, 'rows': rows, 'sources': sources, 'runtime': runtime,
         'protected_trees': protected,
         'original_inventory_resolution': source_resolution[frozen.INVENTORY_SOURCE],
         'input_hashes': {r['curve_key']: sha(out/'inputs'/f"{r['curve_key']}.json") for r in rows},
         'selection': 'user-directed retrospective two-field engineering pilot, not a prospective validation panel',
         'reference': 'https://arxiv.org/abs/1606.07178 (sections 4-5)',
         'method_difference': 'finite exact index-form height descent and smooth-norm enumeration; no claim to implement KSW production Julia reduction/NFS sieving'})
    print('MOD2_PREPARE|PASS|fields=2|180_seconds_each|16384_candidates_each|no_BNF_or_points')


def check(out, complete=False):
    plan = read(out/'plan.json')
    for path, digest in plan['sources'].items():
        require(sha(ROOT/path)==digest, f'pilot source changed: {path}')
    for path, hashes in plan['protected_trees'].items():
        require(historical.tree_hashes(Path(path))==hashes, f'protected evidence changed: {path}')
    require(sha(Path(plan['runtime']['sage']))==plan['runtime']['sage_sha256'], 'Sage changed')
    for r in plan['rows']:
        key = r['curve_key']
        require(sha(out/'inputs'/f'{key}.json')==plan['input_hashes'][key], 'worker input changed')
        if complete:
            s = read(out/key/'supervisor.json')
            require(s['input_sha256']==plan['input_hashes'][key], 'supervisor input mismatch')
            require(s['wall_limit_seconds']==plan['protocol']['wall_seconds_per_field'], 'worker budget mismatch')
            if s['status']=='PASS_BOUNDED_WORKER':
                result = read(out/key/'result.json')
                require(result['input_sha256']==plan['input_hashes'][key], 'result input mismatch')
                require(result['global_class_2rank_estimate'] is None and result['global_class_2rank_upper_bound'] is None, 'unsupported class-rank claim')
                require(result['status']=='PASS_BOUNDED_RELATION_PROTOTYPE_NOT_CLASS_RANK', 'invalid worker status')
                require(result['candidate_prefixes'][-1]['attempted']==plan['protocol']['candidate_budget'], 'incomplete attempt count')
            else:
                require(s['status'].startswith('UNKNOWN'), 'invalid failure status')
    return plan


def run_one(out, row, plan):
    key = row['curve_key']
    dest = out/key
    if (dest/'supervisor.json').exists():
        return read(dest/'supervisor.json')['status']
    require(not list(dest.iterdir()), 'unresolved worker output: do not overwrite partial evidence')
    command = [plan['runtime']['sage'], '-python', str(WORKER), '--input', str(out/'inputs'/f'{key}.json'), '--output', str(dest)]
    limit = plan['protocol']['wall_seconds_per_field']
    start = time.monotonic()
    proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, start_new_session=True)
    try:
        output, _ = proc.communicate(timeout=limit)
        status = 'PASS_BOUNDED_WORKER' if proc.returncode==0 and (dest/'result.json').exists() else 'UNKNOWN_WORKER_FAILURE'
    except subprocess.TimeoutExpired:
        os.killpg(proc.pid, signal.SIGKILL)
        output, _ = proc.communicate()
        status = 'UNKNOWN_TIMEOUT'
    (dest/'worker.log').write_text(output)
    save(dest/'supervisor.json', {'status': status, 'command': command, 'returncode': proc.returncode,
         'input_sha256': plan['input_hashes'][key], 'wall_limit_seconds': limit,
         'wall_seconds': round(time.monotonic()-start, 6), 'no_new_point_searches': True, 'no_full_class_group': True})
    return status


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('command', choices=['prepare', 'run', 'check'])
    ap.add_argument('--output', type=Path, default=OUT)
    args = ap.parse_args()
    out = args.output.resolve()
    if args.command=='prepare':
        prepare(out)
    elif args.command=='run':
        plan = check(out)
        with ThreadPoolExecutor(max_workers=plan['protocol']['max_jobs']) as pool:
            jobs = {pool.submit(run_one, out, r, plan): r for r in plan['rows']}
            for job in as_completed(jobs):
                print(f"MOD2_WORKER|t={jobs[job]['t']}|status={job.result()}", flush=True)
        check(out, True)
    else:
        check(out, True)
        print('MOD2_CHECK|PASS|old_evidence_unchanged=PASS')


if __name__=='__main__':
    main()
