#!/usr/bin/env sage -python
"""Audit retained full-cloud receipts and the failed, frozen promotion gate.

Default execution checks integrity and accounting, not new arithmetic. The
original independent arithmetic replay ran inside every charged worker.
"""
import argparse
from collections import Counter
import json
from pathlib import Path

from cancellation_cloud_training import OUT, RAW
from finite_cancellation_corpus import ROOT, canonical, digest, write


def read(path): return json.loads(path.read_text())
def sha(path): return digest(path.read_bytes())
def need(value, message):
    if not value: raise ArithmeticError(message)


def audit(before_pack=False):
    plan = read(OUT/'protocol.json'); ph = sha(OUT/'protocol.json')
    supervision = read(OUT/'supervision.json'); summary = read(OUT/'summary.json')
    cases = {c['id']: c for c in read(OUT/'inputs.json')}
    need(len(cases) == 24 and supervision['status'] == 'COMPLETE', 'cohort incomplete')
    need(len(supervision['records']) == 48, 'missing arm receipt')
    need(ph == supervision['protocol_sha256'] == summary['protocol_sha256'], 'protocol differs')
    for name, h in plan['source_sha256'].items():
        need(sha(ROOT/name) == h, 'frozen source differs: '+name)
    for name, h in plan['input_sha256'].items():
        need(sha(OUT/name) == h, 'frozen input differs: '+name)
    for name, h in read(OUT/'dependency-audit.json')['extra_source_sha256'].items():
        need(sha(ROOT/name) == h, 'supplementary checker differs: '+name)
    test = {c['j_group'] for c in cases.values()}
    roster = read(OUT/'roster.json'); fit = read(OUT/'fit.json'); cloud = read(OUT/'cloud-fit.json')
    need(len(test) == 24 and not test & set(roster['prior_cpu_j_groups']), 'CPU holdout reused')
    need(not test & (set(fit['training_j_groups']) | set(cloud['development_j_groups'])), 'fit leakage')
    need(read(OUT/'cloud-training-replay.json')['cloud_fit_sha256'] == sha(OUT/'cloud-fit.json'), 'cloud fit replay differs')
    totals = {a: dict(cases=0, target_completions=0, directions=0, cpu_seconds=0., calls=0) for a in plan['arms']}
    statuses = Counter(); endings = {a: Counter() for a in plan['arms']}; seen = set(); unknowns = 0
    records = []
    for receipt in supervision['records']:
        ident, arm = key = receipt['case'], receipt['arm']
        need(key not in seen and ident in cases and arm in totals, 'unexpected or duplicate arm')
        seen.add(key); dest = RAW/'arms'/ident/arm; result = read(dest/'result.json')
        proof = read(dest/'independent-verification.json'); packet = read(dest/'rank-input.json')
        case = cases[ident]; seed = case['seed']
        need(receipt == read(dest/'supervisor.json') and receipt['status'] == 'COMPLETE', 'supervisor differs')
        need(receipt['result_sha256'] == sha(dest/'result.json'), 'result bytes differ')
        need(result['case'] == ident and result['arm'] == arm and result['protocol_sha256'] == ph, 'result identity differs')
        need(not result['development_only'] and result['target_directions'] == 2, 'control scope differs')
        need(result['events_sha256'] == sha(dest/'events.json'), 'event bytes differ')
        need(result['independent_verification_sha256'] == sha(dest/'independent-verification.json') and proof['status'] == 'PASS', 'rank receipt differs')
        n, g = result['initial_rank'], result['new_directions']
        need(n == len(seed['points']) == 20 and g == result['candidate_directions'] == proof['new_directions'], 'direction count differs')
        need(n+g == result['rank_lower_bound'] == proof['rank']['rank'] == len(packet['points']), 'rank accounting differs')
        need(packet['curve'] == seed['curve'] and packet['points'][:n] == seed['points'], 'starting subgroup differs')
        need(result['success'] == (g >= 2), 'target accounting differs')
        new = []
        for call in result['calls']:
            path = dest/call['file']; row = read(path)
            need(sha(path) == call['sha256'], 'point call bytes differ')
            need(row['centre'] == case['centres'][row['decision']['index']], 'fixed bank differs')
            need(call['new_directions'] == len(row['new_points']), 'cloud accounting differs')
            need(call['status'] == row['search']['status'], 'call status differs')
            new.extend(row['new_points']); statuses[call['status']] += 1
        need(packet['points'][n:] == new and len(new) == g, 'full cloud certificate differs')
        total = totals[arm]; total['cases'] += 1; total['target_completions'] += result['success']
        total['directions'] += g; total['cpu_seconds'] += receipt['charged_cpu_seconds']; total['calls'] += len(result['calls'])
        unknowns += proof['preparation_unknowns']; endings[arm][result['status']] += 1
        records.append(dict(case=ident, arm=arm, rank=n+g, directions=g,
            result_sha256=sha(dest/'result.json'), packet_sha256=sha(dest/'rank-input.json'),
            independent_verification_sha256=sha(dest/'independent-verification.json')))
    for arm, total in totals.items():
        for name, value in total.items():
            need(abs(value-summary['totals'][arm][name]) < 1e-9, 'reported total differs')
    need(dict(statuses) == summary['point_call_statuses'] and unknowns == summary['preparation_unknowns'], 'status accounting differs')
    from cancellation_cloud_programme import compare
    replay_totals, comparison = compare(summary['rows'])
    need(replay_totals == summary['totals'] and comparison == summary['comparison'], 'frozen statistical gate differs')
    need(comparison['gate_passed'] is False and summary['fresh_fibre_gate'] is False, 'failed gate was promoted')
    boundary = read(OUT/'subgroup-boundary.json')
    need(boundary['status'] == 'PASS_SAME_CANCELLATION_DIFFERENT_SUBGROUP_LABEL', 'subgroup boundary missing')
    need([r['rank'] for r in boundary['rank_chain']] == [18, 19, 20], 'subgroup witness scope differs')
    for name, h in boundary['source_sha256'].items():
        need(sha(ROOT/name) == h, 'subgroup witness source differs')
    result = dict(status='PASS_RETAINED_FULL_CLOUD_RECEIPTS', arms=len(seen), totals=totals,
        point_statuses=dict(statuses), endings={a:dict(v) for a,v in endings.items()},
        preparation_unknowns=unknowns, records=records,
        bindings_sha256={n:sha(OUT/n) for n in ('protocol.json','summary.json','cloud-fit.json','cloud-training-replay.json','subgroup-boundary.json','dependency-audit.json')},
        checker_sha256=sha(Path(__file__)), rank32='UNKNOWN',
        boundary='Integrity/accounting audit of independent arithmetic receipts, all charged in the original 48 workers. The two-direction policy failed its frozen gate; no fresh two-direction trial is authorized by this result. The exact subgroup example rules out universally correct subgroup-blind independence labels, not useful average cancellation prediction.')
    target = OUT/'retained-audit.json'
    if target.exists(): need(read(target) == result, 'retained audit differs')
    else:
        need(before_pack, 'missing pre-pack audit'); write(target, result)
    if not before_pack:
        import pack_cancellation_scheduler as bundle
        bundle.OUT = OUT; bundle.RAW = RAW; bundle.verify()
    print(json.dumps({k:result[k] for k in ('status','arms','totals','endings','point_statuses')}, indent=2), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--before-pack', action='store_true'); p.add_argument('--arithmetic', action='store_true')
    args = p.parse_args()
    if args.arithmetic:
        from verify_cancellation_cloud import main
        main()
    audit(args.before_pack)
