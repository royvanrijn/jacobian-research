#!/usr/bin/env python3
"""Audit the frozen 64-case V4 screen's records; never call PARI or Sage.

This checks coverage, identities, censoring and summaries of retained results.
It does not independently prove the recorded rank bounds.
"""

from collections import Counter
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PRODUCER = 'elkies-k3/scripts/screen_r17_norm12_11952_v4_base_jacobian_ranks.py'
SHORTLIST = 'artifacts/generated-results/elkies-k3-r17-norm12-11952-v4-pair-shortlist-64-v1.json'
SCREEN = 'artifacts/generated-results/elkies-k3-r17-norm12-11952-v4-base-rank-screen-v1.json'
CAMPAIGN = 'artifacts/generated-results/elkies-k3-r17-all17-product-toric-frobenius-campaign-v1.json'
SINGLE = 'artifacts/generated-results/elkies-k3-r17-product-19bad-083ad-rank-zero-v1.json'
SWEEP = 'artifacts/generated-results/elkies-k3-r17-product-regulator-sweep-v1.json'
HASHES = {
    PRODUCER: 'd139fd59dbf66745d90896d3fab7a3fcd71c82467dd4e09e80b6a04b14186e20',
    SHORTLIST: '4739e0b24b00e276269b228d7c85a01743caa1dcbc62ef58ea59ca3152d116a0',
    SCREEN: 'ed41526fca5724f529aaf73d3f695fe2a781eaf4ea6b6fb792a2af5676759f69',
    CAMPAIGN: '9b9467f6c1a754f41f9feeed6be0ae8c13d275e4203f9cd358df08705ef7318c',
    SINGLE: '8c02b033f4ae85ba3b490863756054328cca3ec35ebc9794ebe358e4058498f0',
    SWEEP: 'a7a9f98e3186deb3b437cb80b9c3cfb68f7a237a182a9a21e56c4e3aaf48201f',
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def validate(payload, shortlist):
    require(payload['schema'] == 'elkies-k3.r17-norm12-11952-v4-base-rank-screen.v1', 'screen schema')
    require(payload['status'] == 'PASS_BOUNDED_EXACT_BASE_JACOBIAN_RANK_INTERVAL_SCREEN', 'screen status')
    require(shortlist['status'] == 'PASS_EXACT_BOUNDED_RATIONAL_V4_PAIR_SHORTLIST', 'shortlist status')
    inputs = {name: HASHES[name] for name in (PRODUCER, SHORTLIST)}
    require(payload['inputs'] == inputs, 'input identities')
    require(payload['limits'] == dict(shortlist_prefix=64, timeout_seconds_per_pair=10.0,
                                     concurrent_workers=2, pari_effort=0), 'frozen screen limits')
    pairs, rows = shortlist['pairs'], payload['results']
    require(len(pairs) == len(rows) == 64, '64 distinct cases required')
    require([p['shortlist_rank'] for p in pairs] == list(range(1, 65)), 'shortlist coverage')
    require(len({p['pair_key'] for p in pairs}) == 64, 'duplicate shortlist identity')
    require([(r['shortlist_rank'], r['pair_key']) for r in rows] ==
            [(p['shortlist_rank'], p['pair_key']) for p in pairs], 'screen coverage or identity')
    completed, censored = [], []
    for row in rows:
        require(row['run_key'] == dict(pair_key=row['pair_key'], pari_effort=0,
                                      script_sha256=HASHES[PRODUCER],
                                      shortlist_sha256=HASHES[SHORTLIST], timeout_seconds=10.0),
                'case provenance or execution mode')
        if row['status'] == 'completed':
            require(row['rank_status'] == 'EXACT_INTERVAL', 'completed interval status')
            lower, upper = row['rank_lower_bound'], row['rank_upper_bound']
            require(type(lower) is int and type(upper) is int and 0 <= lower <= upper,
                    'invalid rank interval')
            require(row['pari_effort'] == 0, 'case effort')
            completed.append(row)
        else:
            require(row['status'] == 'timeout' and row['rank_status'] == 'UNKNOWN',
                    'censored case must stay UNKNOWN')
            require(not {'rank_lower_bound', 'rank_upper_bound', 'sha_information',
                         'independent_points_found'} & row.keys(), 'timeout has a rank result')
            censored.append(row)
    require([(r['shortlist_rank'], r['pair_key']) for r in censored] == [
        (3, 'alternate-orbit-1463f:alternate-orbit-083ad'),
        (30, 'alternate-orbit-19b4e:alternate-orbit-146dc')], 'frozen timeout identities')
    intervals = Counter((r['rank_lower_bound'], r['rank_upper_bound']) for r in completed)
    require(intervals == {(1, 1): 17, (0, 2): 27, (0, 4): 8, (1, 3): 10}, 'interval histogram')
    summary = dict(completed=len(completed), timeouts=len(censored), errors=0,
                   exact_rank_count=sum(r['rank_lower_bound'] == r['rank_upper_bound'] for r in completed),
                   minimum_completed_upper_bound=min(r['rank_upper_bound'] for r in completed))
    require(payload['summary'] == summary, 'summary differs from the cases')
    ranked = sorted(completed, key=lambda r: (r['rank_upper_bound'], r['rank_lower_bound'], r['shortlist_rank']))
    require(payload['completed_pairs_ranked_by_upper_then_lower_bound'] ==
            [r['pair_key'] for r in ranked], 'ranked list differs from the cases')
    return summary


def validate_product_closures(screen, campaign, single, sweep):
    selected = {r['pair_key'] for r in screen['results'] if r['status'] == 'completed'
                and r['rank_lower_bound'] == r['rank_upper_bound'] == 1}
    require(len(selected) == 17, 'seventeen selected base Jacobians required')
    require(len(campaign['targets']) == 17 and
            {r['pair_key'] for r in campaign['targets']} == selected, 'product target coverage')
    geometric = {r['pair_key'] for r in campaign['targets'] if r['geometric_rank_zero']}
    require(len(geometric) == campaign['rank_zero_count'] == 12, 'twelve geometric closures required')
    for row in campaign['targets']:
        require(row['best_geometric_mw_rank_upper_bound'] == (0 if row['geometric_rank_zero'] else 2),
                'geometric rank boundary')
    require(single['status'] == 'PROVED_ARITHMETIC_PRODUCT_TWIST_RANK_ZERO'
            and type(single['rank_over_QQ_u']) is int and single['rank_over_QQ_u'] == 0,
            'single-target arithmetic closure')
    arithmetic = {single['pair_key']} | {r['pair_key'] for r in sweep['targets']}
    require(len(sweep['targets']) == sweep['rank_zero_count'] == 4 and len(arithmetic) == 5,
            'five distinct arithmetic closures required')
    require(not geometric & arithmetic and selected == geometric | arithmetic,
            'arithmetic closures do not cover the five geometric survivors')
    for row in sweep['targets']:
        require(row['rank_over_QQ_u'] == dict(lower=0, upper=0, status='PROVED')
                and row['section_solving_eligible'] is False, 'arithmetic closure or section queue')
        require(row['rank_over_QQbar_u'] == dict(lower=0, upper=2, status='UNKNOWN'),
                'arithmetic rank zero must not become geometric rank zero')
    require(sweep['regulator_compatible_count'] == 0 and sweep['explicit_section_solving_queue'] == [],
            'retired section-solving queue reopened')


def audit():
    for name, digest in HASHES.items():
        require(hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == digest,
                f'frozen source or output changed: {name}')
    summary = validate(json.loads((ROOT / SCREEN).read_text()), json.loads((ROOT / SHORTLIST).read_text()))
    validate_product_closures(*(json.loads((ROOT / path).read_text())
                                for path in (SCREEN, CAMPAIGN, SINGLE, SWEEP)))
    print(f"PASS V4 screen records: {summary['completed']} completed intervals, "
          f"{summary['exact_rank_count']} exact ranks, {summary['timeouts']} timeouts remain UNKNOWN; no rank replay")
    print('PASS selected product records: 12 geometric + 5 arithmetic rank-zero closures; section queue empty; no cohomology replay')


if __name__ == '__main__':
    audit()
