#!/usr/bin/env sage -python
"""Frozen, equally charged basis-amplification control; no residue refitting."""
import argparse
from fractions import Fraction as F
import gzip
import json
from pathlib import Path
import resource
import signal
import time

import cancellation_scheduler_prepare as preparation
import cancellation_scheduler_round4 as supervisor
from cancellation_scheduler_fresh import source_closure
from finite_cancellation_corpus import ROOT, LOCAL, OUT as CORPUS, canonical, digest, write

CAS = Path(__file__).resolve().parent
OUT = ROOT/'artifacts/generated-results/elliptic-curves/cancellation_basis_amplification_v1'
RAW = LOCAL/'cancellation-basis-amplification-v1'
PRIOR = ROOT/'artifacts/generated-results/elliptic-curves/cancellation_cloud_scheduler_v1'


def read(path): return json.loads(path.read_text())
def sha(path): return digest(path.read_bytes())
def need(condition, message):
    if not condition:
        raise ArithmeticError(message)


def parent_subset(landscape):
    generic = landscape['generic_rank']; rows = []
    for entry in landscape['anchors'][:16]:
        row = entry['anchor']; word = row['representative']
        need(not any(word[generic:]), 'parent anchor uses an exceptional direction')
        rows.append({'mask': row['orbit'], 'word': word[:generic], 'norm': int(2*row['shell'])})
    need(len(rows) == 16 and len({r['mask'] for r in rows}) == 16, 'retained parent subset differs')
    return {'status': 'COMPLETE_FROZEN_COMPLEMENT_PARENT_SUBSET', 'dimension': generic,
            'generic_gram_scale': 2, 'shells': sorted({r['norm'] for r in rows}), 'rows': rows,
            'boundary': 'First16 generic anchors of the retained compatible landscape, unchanged words and shell labels. No new generic census or productivity claim.'}


def prepare():
    design = {
        'status': 'DESIGN_BEFORE_BASIS_AMPLIFICATION_OUTCOMES',
        'arms': ['fixed_bank', 'basis_refresh'], 'cases': 24,
        'target_directions': 2, 'initial_rank': 18, 'maximum_centres': 48,
        'height': 125000, 'maximum_calls': 288, 'maximum_extension_dimension': 4,
        'working_cpu_seconds': 40, 'point_wall_seconds': 5, 'bank_wall_seconds': 30,
        'hard_process_cpu_seconds': 100, 'arm_wall_seconds': 150, 'prime_bound': 500,
        'bank_policy': {'anchors_per_shell': 16, 'canonical_per_shell': 25, 'exact_cvp_node_limit': 2000000},
        'selection': 'Four new whole-j CPU groups from each of six retained R17 families. Exclude every preceding scheduler/control j group using the cloud roster. Require retained rank18 seed and endpoint at least20. Select lexical first compatible landscape per j, then SHA256(cancellation-basis-amplification-v1/ + j) within family. Endpoint eligibility constructs known-direction controls only. No outcomes choose the roster, no replacement.',
        'hypothesis': 'After a complete cloud is reconciled and its full admitted basis independently certified, replacing the actual search subgroup and verifying new compatible anchors improves later-direction yield per complete CPU.',
        'ablation': 'Both arms build the same initial bank from the same16 retained generic anchors and rank18 seed, charge production plus rational-CVP replay, and use the same factor-free/prime-neighbour box order. fixed_bank keeps that search basis while accumulating certified output. basis_refresh rebuilds immediately after a positive cloud unless its two-direction target or CPU cap has been reached. There is no residue-weight, radial, hazard or cloud-yield fit.',
        'exploration': 'For each anchor execute factor-free H125000 then every distinct retained neighbour (at most2) at H125000 before opening another anchor. A zero empirical radius mass never skips a job. Only exactly equivalent completed boxes on the same actual anchor are skipped, including infinity. A finite exact boundary probe records coordinates beyond old H125000 and outside completed same-anchor boxes; failure to find a probe witness remains UNKNOWN and never prunes. No optimality or exhaustive-neighbour claim.',
        'subgroup': 'Consider every returned point with the same finite places through500 and independently Sage-certify every enlarged basis before reuse. Finite columns in the span stay UNKNOWN. Rebuild using the complete admitted cloud, never just enough points for the stopping target. Reset all epoch-local anchors/maps/job indices; retain only exact completed-box identities across epochs.',
        'cost': '40 CPU seconds per arm for interpreter, imports, initial bank, every rebuild, producer, rational reference, anchor identities, all maps, point calls, complete-cloud classification and each independent rank certificate. Atomic operations can overrun this launch cap. Final transcript/policy/rank replay is mandatory, separately measured, and charged in complete isolated child-process CPU. Both arms share the same100 CPU/150 wall hard limits. Retained seed/parent construction is the common input interface. Shared input preflight CPU is reported separately and charged to the candidate in conservative sensitivity. No elapsed-time proxy, free rebuild, or outcome-dependent retry.',
        'primary_gate': 'Candidate must match or exceed baseline two-direction target completions and uncapped total directions. Later-cloud directions (all gains after the first positive full cloud) per COMPLETE CPU must improve by at least10 percent, with central97.5-percent paired family bootstrap interval wholly above1;10000 draws, seed20260914. Include zero-gain and censored cases in all denominators. Zero baseline later yield or any undefined bootstrap sample makes the gate UNKNOWN. Any infrastructure or preparation unknown blocks promotion. Report uncapped all-direction rate separately.',
        'continuation': 'Only this fixed new control is authorized. No fresh-fibre continuation, repeat of failed controls, enlarged allowance, or automatic campaign follows. A passed gate would support proposing a separately frozen transfer experiment.',
        'boundary': 'Masked retained-corpus controls, not externally pristine population data, new rank discoveries, exact ranks, saturation, or rank32. Initial rank18 limits extension cost near the one-case integration evidence. The prior rank20 fixed-bank cloud failure is unchanged.'}
    preparation.new_write(OUT/'design.json', design)
    old = read(PRIOR/'roster.json'); excluded = set(old['excluded_training_j_groups'])
    eligible = json.loads(gzip.decompress((preparation.V2/'eligible_banks.json.gz').read_bytes()))
    banks = {}
    for row in sorted(eligible, key=lambda r: r['source_landscape']):
        if row['j_group'] not in excluded and row['initial_rank'] == 18 and row['endpoint_rank'] >= 20:
            banks.setdefault(row['j_group'], row)
    chosen = []
    for family in preparation.FAMILIES:
        rows = sorted((r for r in banks.values() if r['family'] == family),
                      key=lambda r: digest(('cancellation-basis-amplification-v1/'+r['j_group']).encode()))
        need(len(rows) >= 4, 'insufficient family; no refill'); chosen.extend(dict(r) for r in rows[:4])
    corpus = {r['id']: r for r in json.loads(gzip.decompress((CORPUS/'corpus.json.gz').read_bytes()))}
    inputs = []; oracle = []
    for row in chosen:
        need(sha(ROOT/row['source_seed']) == row['seed_sha256'] and sha(ROOT/row['source_landscape']) == row['landscape_sha256'], 'retained input changed')
        seed = read(ROOT/row['source_seed']); landscape = read(ROOT/row['source_landscape'])
        need(seed['points'] == landscape['basis'] and landscape['generic_rank'] == 17, 'seed/landscape interface differs')
        ident = digest(canonical(['basis-amplification-v1', row['j_group'], row['source_landscape']]))[:20]; row['id'] = ident
        inputs.append({'id': ident, 'j_group': row['j_group'], 'family': row['family'], 'stratum': 'M18',
                       'seed': {k: seed[k] for k in ('curve', 'points', 'proof')},
                       'parent_bank': parent_subset(landscape)})
        endpoint = corpus[row['corpus_id']]
        oracle.append({'id': ident, 'curve': endpoint['curve'], 'points': endpoint['generic_points']+endpoint['targets'],
                       'source': endpoint['source'], 'proof_sha256': endpoint['proof_sha256']})
    preparation.new_write(OUT/'inputs.json', inputs); preparation.new_write(OUT/'oracle.json', oracle)
    preparation.new_write(OUT/'roster.json', {'chosen': chosen, 'prior_cpu_j_groups': sorted(excluded),
        'excluded_training_j_groups': sorted(excluded | {c['j_group'] for c in inputs}),
        'previous_roster_sha256': sha(PRIOR/'roster.json'), 'eligible_banks_sha256': sha(preparation.V2/'eligible_banks.json.gz'),
        'inputs_sha256': sha(OUT/'inputs.json'), 'oracle_sha256': sha(OUT/'oracle.json')})
    print(json.dumps({'status': 'FROZEN_BASIS_ABLATION_ROSTER', 'cases': len(inputs), 'arms': design['arms']}), flush=True)


def preflight():
    preparation.OUT = OUT; preparation.preflight()
    need(all(r['seed']['rank'] == 18 and r['endpoint']['rank'] >= 20 for r in read(OUT/'preflight.json')['rows']), 'two-direction eligibility not independently certified')


def smoke():
    """Retained development cloud -> full certified basis -> compatible bank.

    No point search. Neither its directions nor timing enters the holdout.
    """
    from cancellation_basis_epoch import rebuild, native_rank
    from cancellation_scheduler_cpu import cpu
    from finite_cancellation_features import alarm
    from future_point_admission import FinitePointAdmission
    from memory_rank_certificate import checked_rank
    from pointed_quartic_search import PointedQuarticSearch
    from pari_pointed_backend import replay
    from v3_warm_engine import certified_state
    need(not (OUT/'development-smoke.json').exists() and not (RAW/'development/start.json').exists(), 'preserve development attempt')
    row = next(r for r in read(PRIOR/'cloud-fit.json')['records'] if r['directions'] >= 2)
    old = next(c for c in read(preparation.V2/'inputs.json') if c['id'] == row['case'])
    eligible = json.loads(gzip.decompress((preparation.V2/'eligible_banks.json.gz').read_bytes()))
    source = next(r for r in eligible if r['j_group'] == old['j_group'] and r['initial_rank'] == len(old['seed']['points']))
    landscape = read(ROOT/source['source_landscape']); need(landscape['basis'] == old['seed']['points'], 'development source basis differs')
    case = {**old, 'parent_bank': parent_subset(landscape)}; plan = read(OUT/'design.json'); dest = RAW/'development'
    preparation.new_write(dest/'start.json', {'status': 'DEVELOPMENT_ONLY_NO_POINT_SEARCH', 'source': row, 'case': case})
    resource.setrlimit(resource.RLIMIT_CPU, (60, 60)); signal.signal(signal.SIGALRM, alarm); signal.alarm(90)
    tick = cpu(); seed = case['seed']; curve = tuple(map(F, seed['curve'])); basis = [tuple(map(F, p)) for p in seed['points']]
    state = certified_state(curve, basis, seed['proof']); call = read(ROOT/row['source'])
    need(sha(ROOT/row['source']) == row['sha256'], 'retained cloud changed')
    search = PointedQuarticSearch(state=state, centre={'coefficients': call['centre']}, coordinate_policy=call['mapping']['coordinate_policy'])
    admission = FinitePointAdmission(curve, basis, prime_bound=plan['prime_bound'])
    for point in replay(search, call['mapping'], call['search']): admission.consider(point)
    proof = checked_rank(curve, admission.points, admission.primes, seed['proof']['no_rational_2_torsion_prime'])
    packet = {'curve': seed['curve'], 'points': [list(map(str, p)) for p in admission.points], 'proof': proof}
    independent = native_rank(packet); need(independent['rank'] == len(basis)+row['directions'], 'full development cloud differs')
    bank, receipt = rebuild(packet, case, plan, dest/'rebuilt')
    need(bank['seed']['points'] == packet['points'] and receipt['rank'] == independent['rank'], 'rebuild dropped a direction')
    signal.alarm(0)
    preparation.new_write(OUT/'development-smoke.json', {'status': 'PASS_COMPLETE_CLOUD_REBUILD', 'initial_rank': len(basis),
        'enlarged_rank': len(packet['points']), 'cloud_directions': row['directions'], 'independent_rank': independent,
        'bank_verification_sha256': sha(dest/'rebuilt/verification.json'), 'cpu_seconds': cpu()-tick, 'point_search_calls': 0,
        'boundary': 'Previously tested development cloud, independent full-rank reconciliation and new compatible bank; no performance outcome.'})
    print(json.dumps(read(OUT/'development-smoke.json')), flush=True)


def seal():
    from pari_pointed_backend import sources
    from sage.all import pari
    from sage.env import SAGE_VERSION
    import sys
    plan = read(OUT/'design.json')
    need(read(OUT/'preflight.json')['status'] == 'PASS', 'input preflight missing')
    need(read(OUT/'development-smoke.json')['status'] == 'PASS_COMPLETE_CLOUD_REBUILD', 'development integration missing')
    development = RAW/'development-worker-v2'
    need(read(development/'independent-verification.json')['status'] == 'PASS', 'complete runner integration missing')
    preparation.new_write(OUT/'development-worker.json', {
        'status': 'PASS_DEVELOPMENT_WORKER_ONLY',
        'result_sha256': sha(development/'result.json'),
        'verification_sha256': sha(development/'independent-verification.json'),
        'retained_configuration_failure_sha256': sha(RAW/'development-worker/failure.json'),
        'boundary': 'One previously tested group; two directions in its first positive cloud. A prior development recipe omitted the GP hash and failed before any point call; corrected recipe retained separately. Neither attempt enters holdout outcomes.'})
    inputs = read(OUT/'inputs.json'); roster = read(OUT/'roster.json')
    need(len(inputs) == plan['cases'] and len({c['j_group'] for c in inputs}) == plan['cases'], 'whole-j roster differs')
    need(not {c['j_group'] for c in inputs} & set(roster['prior_cpu_j_groups']), 'previous controls reused')
    source_names = [Path(__file__), CAS/'cancellation_basis_epoch.py', CAS/'cancellation_basis_exploration.py',
                    CAS/'verify_cancellation_basis_amplification.py', CAS.parent/'tests/test_cancellation_basis_exploration.py']
    closure = source_closure(source_names + [ROOT/p for p in sources()])
    preparation.new_write(OUT/'environment.json', {'sage': SAGE_VERSION, 'pari': str(pari('version()')), 'python': sys.version,
        'threads': {'OPENBLAS_NUM_THREADS': '1', 'OMP_NUM_THREADS': '1'}})
    plan.update(status='FROZEN_BASIS_AMPLIFICATION_CONTROL', source_sha256=closure,
                input_sha256={n: sha(OUT/n) for n in ('design.json', 'inputs.json', 'oracle.json', 'roster.json', 'preflight.json', 'development-smoke.json', 'development-worker.json', 'environment.json')},
                gp_sha256=sha(Path('/usr/bin/gp')))
    preparation.new_write(OUT/'protocol.json', plan)
    print(json.dumps({'status': plan['status'], 'sources': len(closure), 'protocol_sha256': sha(OUT/'protocol.json')}), flush=True)


def guard():
    plan = read(OUT/'protocol.json')
    for name, h in plan['source_sha256'].items(): need(sha(ROOT/name) == h, 'sealed source changed: '+name)
    for name, h in plan['input_sha256'].items(): need(sha(OUT/name) == h, 'sealed input changed: '+name)
    need(sha(Path('/usr/bin/gp')) == plan['gp_sha256'], 'GP executable changed')
    return plan, read(OUT/'inputs.json')


def worker(case_id, arm):
    from cancellation_basis_epoch import run as execute
    plan, inputs = guard(); case = next(c for c in inputs if c['id'] == case_id)
    execute(case, arm, plan, RAW/'arms'/case_id/arm, sha(OUT/'protocol.json'))


def run():
    guard(); supervisor.OUT = OUT; supervisor.RAW = RAW; supervisor.ENTRY = Path(__file__).resolve(); supervisor.run()


def report():
    from cancellation_cloud_programme import compare
    from verify_cancellation_basis_amplification import check
    check(); plan, _ = guard(); rows = []
    for receipt in read(OUT/'supervision.json')['records']:
        dest = RAW/'arms'/receipt['case']/receipt['arm']; r = read(dest/'result.json')
        rows.append({**{k: r[k] for k in ('case', 'family', 'stratum', 'arm', 'success', 'new_directions',
                     'later_cloud_directions', 'first_cloud_gain', 'initial_rank', 'rank_lower_bound', 'status', 'unknowns')},
                     'calls': len(r['calls']), 'epochs': len(r['epochs']), 'cpu_seconds': receipt['charged_cpu_seconds'],
                     'bank_cpu_seconds': r['components']['bank'], 'final_replay_cpu_seconds': r['final_replay_cpu_seconds'],
                     'tail_witness_calls': sum(c['tail_witness'] for c in r['calls'])})
    def comparison(later):
        renamed = [{**r, 'arm': 'factor_free' if r['arm'] == 'fixed_bank' else 'adaptive_cloud',
                    'new_directions': r['later_cloud_directions'] if later else r['new_directions']} for r in rows]
        totals, outcome = compare(renamed)
        return {('fixed_bank' if k == 'factor_free' else 'basis_refresh'): v for k, v in totals.items()}, outcome
    totals, overall = comparison(False); later_totals, later = comparison(True)
    candidate, baseline = totals['basis_refresh'], totals['fixed_bank']
    unknowns = sum(r['unknowns'] for r in rows)
    gate = later['gate_passed']
    if candidate['directions'] < baseline['directions'] or candidate['target_completions'] < baseline['target_completions'] or unknowns: gate = False
    shared = read(OUT/'preflight.json')['cpu_seconds']
    ca, cb = later_totals['basis_refresh'], later_totals['fixed_bank']
    conservative = (ca['directions']/(ca['cpu_seconds']+shared))/(cb['directions']/cb['cpu_seconds']) if cb['directions'] else None
    summary = {'status': 'COMPLETE_BASIS_AMPLIFICATION_CONTROL', 'rows': rows, 'totals': totals,
        'all_direction_comparison': overall, 'later_cloud_totals': later_totals, 'later_cloud_comparison': later,
        'promotion_gate': gate, 'preparation_unknowns': unknowns, 'shared_preflight_cpu_seconds': shared,
        'conservative_later_rate_ratio': conservative, 'development_cpu_seconds': read(OUT/'development-smoke.json')['cpu_seconds'],
        'protocol_sha256': sha(OUT/'protocol.json'), 'rank32': 'UNKNOWN', 'fresh_fibres_run': 0,
        'boundary': 'Same unfitted box policy, paired whole-j controls and complete CPU accounting. Later yield counts every gain after the first full positive cloud. Known control directions are not new rank discoveries. No fresh continuation or automatic escalation.'}
    write(OUT/'summary.json', summary)
    print(json.dumps({k: summary[k] for k in ('totals', 'later_cloud_totals', 'later_cloud_comparison', 'promotion_gate', 'preparation_unknowns')}, indent=2), flush=True)


def pack():
    import pack_cancellation_scheduler as bundle
    bundle.OUT = OUT; bundle.RAW = RAW; bundle.pack()


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('command', choices=['prepare', 'preflight', 'smoke', 'seal', 'worker', 'run', 'report', 'pack'])
    p.add_argument('--case'); p.add_argument('--arm'); args = p.parse_args()
    worker(args.case, args.arm) if args.command == 'worker' else globals()[args.command]()
