#!/usr/bin/env sage -python
"""Two fixed validation blocks for latest-generator-first amplification."""
import argparse
from collections import Counter
import gzip
import json
from pathlib import Path

import cancellation_scheduler_prepare as preparation
import cancellation_scheduler_round4 as supervisor
from cancellation_basis_amplification import parent_subset
from cancellation_scheduler_fresh import source_closure
from finite_cancellation_corpus import ROOT, LOCAL, OUT as CORPUS, canonical, digest, write

CAS = Path(__file__).resolve().parent
OUT = ROOT/'artifacts/generated-results/elliptic-curves/cancellation_basis_early_v1'
RAW = LOCAL/'cancellation-basis-early-v1'
PRIOR = ROOT/'artifacts/generated-results/elliptic-curves/cancellation_basis_amplification_v1'
MECHANISM = ROOT/'artifacts/generated-results/elliptic-curves/cancellation_basis_accessibility_v3'


def read(path): return json.loads(path.read_text())
def sha(path): return digest(path.read_bytes())
def need(condition, message):
    if not condition: raise ArithmeticError(message)
def new_write(path, data): preparation.new_write(path, data)


def prepare():
    prior = read(PRIOR/'design.json'); diagnostic = read(MECHANISM/'summary.json')
    need(diagnostic['counts']['strict_coefficient_deletion_targets'] == 4, 'exact causal prerequisite missing')
    plan = {**prior, 'status': 'DESIGN_BEFORE_EARLY_EXPOSURE_OUTCOMES', 'target_directions': 3,
        'validation_blocks': 2, 'cases_per_block': 12,
        'selection': '24 previously untested whole-j groups with retained M18 seed and independently certified endpoint at least21. Exclude the entire prior basis-ablation roster and every earlier scheduler CPU group. Lexical first compatible bank per j, then SHA256(cancellation-basis-early-v1/ + j) within each of six families. First four groups per family, first two assigned block0 and next two block1. Both blocks are fixed before any new outcome; no substitutions or interim policy changes.',
        'hypothesis': 'A newly admitted direction makes some subsequent directions accessible through anchors requiring that generator. Visiting those compatible anchors early, instead of allowing old anchors to consume the budget, improves repeatable later-direction yield per complete CPU.',
        'ablation': 'Both arms start from the same M18 seed and reconstruct the same initial16-parent-subset bank with independent rational-CVP and native anchor checks. Both retain and certify every full cloud. fixed_bank keeps its search basis and bank. basis_refresh rebuilds from the complete enlarged subgroup after each positive cloud until three added directions or the unchanged40-working-CPU allowance. The candidate stably visits anchors using the latest admitted block first, then anchors using earlier admitted generators, then original-subgroup anchors. All other map, box, admission and cost choices match. No missing point, residue fit, Schur diagnostic or radius model enters selection.',
        'early_order': 'At epoch0 use the original exported order. At each later epoch, fresh_start is the rank of the preceding search bank. Stable partition the new exported bank into nonzero coefficients at indices fresh_start onward; otherwise nonzero coefficients at indices18 onward; otherwise zero at indices18 onward. Execute every selected anchor factor-free then at-most-two retained neighbours, subject to the same gain/CPU stops. The stored bank itself is unchanged.',
        'primary_gate': 'Across all24 pairs, later-cloud directions per complete CPU must improve at least10 percent with central97.5-percent paired family bootstrap lower bound above1,10000 draws seed20260914. Candidate must also match or exceed uncapped total directions and three-direction target completions. Both separately frozen12-pair blocks must have later-direction/CPU ratio strictly above1 and at least one actual later gain on an anchor requiring a newly admitted generator. Undefined reference yield or bootstrap samples, any preparation/infrastructure unknown, or any absent independent proof blocks promotion. Report both blocks even if the first fails; never change policy or select the second after outcomes.',
        'mechanism_gate': 'After both blocks finish, audit every later candidate gain on a new-generator anchor. Bind its exact pre-call subgroup, required generator coefficients, new point, square and coordinate map. Prepare old-bank and single-generator-deletion maps before evaluating the discovered point and its negative. Report finite exclusive witnesses and exact counterfactual failures without claiming an entire coset is inaccessible. Reuse numerical rank_growth.py diagnostics only after execution. Performance promotion also requires a strict finite generator-deletion visibility witness in each block; missing postmortem work is UNKNOWN, not success.',
        'continuation': 'Only this fixed24-pair control and its bounded post-execution mechanism audit run here. The40 working CPU,100 hard CPU and150 wall seconds per arm remain unchanged. No failed-primary search is retried, no automatic fresh-fibre search or larger cap follows. An actual candidate-only development integration on one already tested control is separately labelled and never changes the primary roster or counts.',
        'development': 'One previously tested whole-j case6f5afe31891f7dc1938b checks actual early-generator bank transitions and complete final replay with target3 at the same limits. This is a software integration control selected from the completed accessibility audit, not a new validation result. Its original failed-comparison arm remains unchanged. Development timing is outside per-arm performance meters and reported separately.',
        'prior_protocol_sha256': sha(PRIOR/'protocol.json'), 'causal_summary_sha256': sha(MECHANISM/'summary.json'),
        'boundary': 'Development-designed policy on disjoint retained-corpus known controls, in two fixed blocks. No externally pristine population claim, new rank discovery, exact rank, rank32 or state-of-the-art comparison. First-recovery gains cannot pass the later-gain gate. Numerical diagnostics and bounded misses prove no arithmetic exclusion.'}
    new_write(OUT/'design.json', plan)
    previous = read(PRIOR/'roster.json'); excluded = set(previous['excluded_training_j_groups'])
    eligible = json.loads(gzip.decompress((preparation.V2/'eligible_banks.json.gz').read_bytes())); banks = {}
    for row in sorted(eligible, key=lambda r: r['source_landscape']):
        if row['j_group'] not in excluded and row['initial_rank'] == 18 and row['endpoint_rank'] >= 21:
            banks.setdefault(row['j_group'], row)
    chosen = []
    for family in preparation.FAMILIES:
        rows = sorted((r for r in banks.values() if r['family'] == family), key=lambda r: digest(('cancellation-basis-early-v1/'+r['j_group']).encode()))
        need(len(rows) >= 4, 'insufficient family; no refill')
        chosen.extend({**row, 'validation_block': i//2} for i, row in enumerate(rows[:4]))
    chosen.sort(key=lambda r: (r['validation_block'], preparation.FAMILIES.index(r['family'])))
    corpus = {c['id']: c for c in json.loads(gzip.decompress((CORPUS/'corpus.json.gz').read_bytes()))}
    inputs = []; oracle = []
    for row in chosen:
        need(sha(ROOT/row['source_seed']) == row['seed_sha256'] and sha(ROOT/row['source_landscape']) == row['landscape_sha256'], 'retained source changed')
        seed = read(ROOT/row['source_seed']); landscape = read(ROOT/row['source_landscape'])
        need(seed['points'] == landscape['basis'] and landscape['generic_rank'] == 17, 'incompatible seed')
        ident = digest(canonical(['basis-early-v1', row['j_group'], row['source_landscape']]))[:20]; row['id'] = ident
        inputs.append({'id': ident, 'j_group': row['j_group'], 'family': row['family'], 'stratum': 'M18',
            'validation_block': row['validation_block'], 'seed': {k: seed[k] for k in ('curve','points','proof')},
            'parent_bank': parent_subset(landscape)})
        endpoint = corpus[row['corpus_id']]
        oracle.append({'id': ident, 'curve': endpoint['curve'], 'points': endpoint['generic_points']+endpoint['targets'],
            'source': endpoint['source'], 'proof_sha256': endpoint['proof_sha256']})
    new_write(OUT/'inputs.json', inputs); new_write(OUT/'oracle.json', oracle)
    new_write(OUT/'roster.json', {'chosen': chosen, 'prior_cpu_j_groups': sorted(excluded),
        'excluded_training_j_groups': sorted(excluded | {c['j_group'] for c in inputs}),
        'previous_roster_sha256': sha(PRIOR/'roster.json'), 'eligible_banks_sha256': sha(preparation.V2/'eligible_banks.json.gz'),
        'inputs_sha256': sha(OUT/'inputs.json'), 'oracle_sha256': sha(OUT/'oracle.json')})
    print(json.dumps({'status': 'FROZEN_EARLY_EXPOSURE_ROSTER', 'cases': len(inputs), 'blocks': Counter(c['validation_block'] for c in inputs)}), flush=True)


def preflight():
    preparation.OUT = OUT; preparation.preflight()
    need(all(r['seed']['rank'] == 18 and r['endpoint']['rank'] >= 21 for r in read(OUT/'preflight.json')['rows']), 'three-direction eligibility missing')


def development():
    from cancellation_basis_early_epoch import run
    plan = read(OUT/'design.json'); plan['gp_sha256'] = sha(Path('/usr/bin/gp'))
    case = next(c for c in read(PRIOR/'inputs.json') if c['id'] == '6f5afe31891f7dc1938b')
    need(case['j_group'] not in {c['j_group'] for c in read(OUT/'inputs.json')}, 'development overlaps validation')
    recipe = {'case': case, 'plan': plan, 'status': 'DEVELOPMENT_INTEGRATION_ONLY',
              'source_sha256': source_closure([Path(__file__), CAS/'cancellation_basis_early_epoch.py', CAS/'verify_cancellation_basis_early.py'])}
    new_write(RAW/'development-recipe.json', recipe)
    run(case, 'basis_refresh', plan, RAW/'development', 'DEVELOPMENT:'+sha(RAW/'development-recipe.json'))


def seal():
    from pari_pointed_backend import sources
    from sage.all import pari
    from sage.env import SAGE_VERSION
    import sys
    plan = read(OUT/'design.json'); inputs = read(OUT/'inputs.json'); pre = read(OUT/'preflight.json')
    result = read(RAW/'development/result.json'); verification = read(RAW/'development/independent-verification.json')
    need(pre['status'] == 'PASS' and len(pre['rows']) == 24, 'input preflight missing')
    need(verification['status'] == 'PASS' and len(result['epochs']) >= 2, 'actual early-bank integration missing')
    need(read(OUT/'order-regressions.json')['status'] == 'PASS', 'order/alarm regressions missing')
    need(len({c['j_group'] for c in inputs}) == 24 and Counter(c['validation_block'] for c in inputs) == {0:12,1:12}, 'invalid validation split')
    new_write(OUT/'development.json', {'status': 'PASS_DEVELOPMENT_ONLY', 'result_sha256': sha(RAW/'development/result.json'),
        'verification_sha256': sha(RAW/'development/independent-verification.json'),
        'cpu_seconds': result['total_internal_cpu_seconds'], 'directions': result['new_directions'],
        'epochs': len(result['epochs']), 'boundary': plan['development']})
    new_write(OUT/'environment.json', {'sage': SAGE_VERSION, 'pari': str(pari('version()')), 'python': sys.version,
        'threads': {'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'}})
    names = [Path(__file__), CAS/'cancellation_basis_early_epoch.py', CAS/'cancellation_basis_early_order.py',
             CAS/'verify_cancellation_basis_early.py', CAS/'audit_cancellation_basis_early.py',
             CAS/'check_cancellation_basis_early.py']+[ROOT/n for n in sources()]
    plan.update(status='FROZEN_TWO_BLOCK_EARLY_GENERATOR_CONTROL', source_sha256=source_closure(names),
        gp_sha256=sha(Path('/usr/bin/gp')),
        input_sha256={n: sha(OUT/n) for n in ('design.json','inputs.json','oracle.json','roster.json','preflight.json','development.json','environment.json','order-regressions.json')})
    new_write(OUT/'protocol.json', plan)
    print(json.dumps({'status': plan['status'], 'protocol_sha256': sha(OUT/'protocol.json'), 'sources': len(plan['source_sha256'])}), flush=True)


def guard():
    plan = read(OUT/'protocol.json')
    for n, h in plan['source_sha256'].items(): need(sha(ROOT/n) == h, 'sealed source changed: '+n)
    for n, h in plan['input_sha256'].items(): need(sha(OUT/n) == h, 'sealed input changed: '+n)
    need(sha(Path('/usr/bin/gp')) == plan['gp_sha256'], 'GP executable changed')
    return plan, read(OUT/'inputs.json')


def worker(case_id, arm):
    from cancellation_basis_early_epoch import run
    plan, inputs = guard(); case = next(c for c in inputs if c['id'] == case_id)
    run(case, arm, plan, RAW/'arms'/case_id/arm, sha(OUT/'protocol.json'))


def run():
    guard(); supervisor.OUT = OUT; supervisor.RAW = RAW; supervisor.ENTRY = Path(__file__).resolve(); supervisor.run()


def report():
    from cancellation_cloud_programme import compare
    from verify_cancellation_basis_early import check
    check(); plan, inputs = guard(); cases = {c['id']:c for c in inputs}; rows = []
    for receipt in read(OUT/'supervision.json')['records']:
        r = read(RAW/'arms'/receipt['case']/receipt['arm']/'result.json')
        rows.append({**{k:r[k] for k in ('case','family','stratum','arm','success','new_directions','later_cloud_directions','first_cloud_gain','unknowns','status')},
            'validation_block': cases[r['case']]['validation_block'], 'calls': len(r['calls']), 'epochs': len(r['epochs']),
            'cpu_seconds': receipt['charged_cpu_seconds'], 'bank_cpu_seconds':r['components']['bank']})
    def comparison(selected, later):
        renamed = [{**r,'arm':'factor_free' if r['arm']=='fixed_bank' else 'adaptive_cloud',
            'new_directions':r['later_cloud_directions'] if later else r['new_directions']} for r in selected]
        totals, outcome = compare(renamed)
        return {'totals': {('fixed_bank' if k=='factor_free' else 'basis_refresh'):v for k,v in totals.items()}, 'comparison':outcome}
    all_gain = comparison(rows, False); later = comparison(rows, True)
    blocks = {str(b): {'all':comparison([r for r in rows if r['validation_block']==b],False),
                      'later':comparison([r for r in rows if r['validation_block']==b],True)} for b in (0,1)}
    unknowns = sum(r['unknowns'] for r in rows); a,b=all_gain['totals']['basis_refresh'],all_gain['totals']['fixed_bank']
    performance = bool(later['comparison']['gate_passed'] and a['directions']>=b['directions'] and a['target_completions']>=b['target_completions'] and not unknowns)
    repeat = all(block['later']['comparison']['uncapped_direction_rate_ratio'] is not None and block['later']['comparison']['uncapped_direction_rate_ratio']>1 for block in blocks.values())
    mechanism = read(OUT/'mechanism.json') if (OUT/'mechanism.json').exists() else None
    causal = bool(mechanism and mechanism['status']=='PASS' and all(mechanism['blocks'][str(i)]['strict_deletion_witnesses']>0 for i in (0,1)))
    shared = read(OUT/'preflight.json')['cpu_seconds']; la,lb=later['totals']['basis_refresh'],later['totals']['fixed_bank']
    conservative = (la['directions']/(la['cpu_seconds']+shared))/(lb['directions']/lb['cpu_seconds']) if lb['directions'] else None
    conservative_audit = (la['directions']/(la['cpu_seconds']+shared+mechanism['cpu_seconds']))/(lb['directions']/lb['cpu_seconds']) if lb['directions'] and mechanism else None
    summary = {'status':'COMPLETE_FIXED_TWO_BLOCK_EARLY_EXPOSURE_CONTROL','rows':rows,'all':all_gain,'later':later,'blocks':blocks,
        'aggregate_performance_gate':performance,'repeatability_rate_gate':repeat,'actual_mechanism_gate':causal,
        'promotion_gate':performance and repeat and causal,'preparation_unknowns':unknowns,
        'shared_preflight_cpu_seconds':shared,'conservative_later_rate_ratio':conservative,
        'mechanism_cpu_seconds':None if mechanism is None else mechanism['cpu_seconds'],
        'conservative_ratio_including_mechanism':conservative_audit,
        'development':read(OUT/'development.json'),'protocol_sha256':sha(OUT/'protocol.json'),
        'fresh_fibres_run':0,'rank32':'UNKNOWN','boundary':plan['boundary']}
    write(OUT/'summary.json', summary)
    print(json.dumps({k:summary[k] for k in ('all','later','aggregate_performance_gate','repeatability_rate_gate','actual_mechanism_gate','promotion_gate','preparation_unknowns')},indent=2),flush=True)


def pack():
    import pack_cancellation_scheduler as p
    p.OUT=OUT; p.RAW=RAW; p.pack()


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('command',choices=['prepare','preflight','development','seal','worker','run','report','pack']);p.add_argument('--case');p.add_argument('--arm');a=p.parse_args()
    worker(a.case,a.arm) if a.command=='worker' else globals()[a.command]()
