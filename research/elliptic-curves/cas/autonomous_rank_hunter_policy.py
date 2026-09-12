#!/usr/bin/env python3
"""Pure scheduling policy for the autonomous elliptic-curve rank hunter."""
from __future__ import annotations

import math

FAMILIES = ('074d9', '07ca9', '08234', '08f72', '103b2', '11952')
TARGET_RANK = 32


def parameter_height(row):
    return max(abs(int(row['numerator'])), int(row['denominator']))


def height_stratum(row):
    return 'low' if parameter_height(row) <= 1024 else 'high'


def stale_call_limit(rank):
    """No-gain calls tolerated on one fixed parent bank before representation pivot."""
    if rank <= 24:
        return 100
    if rank == 25:
        return 125
    if rank == 26:
        return 150
    if rank == 27:
        return 200
    if rank == 28:
        return 300
    if rank == 29:
        return 450
    if rank == 30:
        return 600
    return 750


def max_bank_count(rank, total_gains, late_gain=False):
    """Maximum distinct parent/representation banks justified by observed cascade evidence."""
    count = 1
    if total_gains >= 3 or late_gain:
        count = 2
    if rank >= 26 and total_gains >= 5:
        count = 3
    if rank >= 28:
        count = max(count, 4)
    if rank >= 30:
        count = max(count, 5)
    return count


def _recent_gain_count(gain_calls, calls, window):
    cutoff = max(0, int(calls) - int(window))
    return sum(int(c) > cutoff for c in gain_calls)


def continuation_decision(summary):
    """Return an autonomous exploit/pivot/retire decision from certified trajectory data.

    Expected fields are plain JSON values. ``complement_calls``/``gain_calls``
    are curve-global; ``bank_calls``/``bank_last_gain_call`` are local to the
    current parent bank. Rank is deliberately only a soft term in the score.
    """
    rank = int(summary.get('rank_lower_bound', 17))
    if rank >= TARGET_RANK:
        return {'action': 'target', 'priority': math.inf, 'reason': 'certified target reached'}

    calls = int(summary.get('complement_calls', 0))
    gain_calls = [int(x) for x in summary.get('gain_calls', [])]
    total_gains = int(summary.get('complement_directions', len(gain_calls)))
    bank_calls = int(summary.get('bank_calls', 0))
    bank_last = summary.get('bank_last_gain_call')
    bank_last = int(bank_last) if bank_last is not None else None
    bank_stale = bank_calls - bank_last if bank_last is not None else bank_calls
    generation = int(summary.get('bank_generation', 0))
    has_suffix = bool(summary.get('has_suffix', False))
    late_gain = bool(summary.get('late_gain', False))
    seed_cloud = int(summary.get('seed_cloud_directions', 0))

    recent100 = _recent_gain_count(gain_calls, calls, 100)
    recent200 = _recent_gain_count(gain_calls, calls, 200)
    limit = stale_call_limit(rank)

    rank_bonus = max(0, rank - 23) * 1.6
    momentum = 9.0 * recent100 + 3.5 * max(0, recent200 - recent100)
    recency = max(0.0, 8.0 * (1.0 - bank_stale / max(1, limit)))
    efficiency = 0.0 if calls == 0 else min(8.0, 80.0 * total_gains / calls)
    cloud_bonus = min(5.0, 1.25 * seed_cloud)
    late_bonus = 6.0 if late_gain else 0.0
    stale_penalty = 7.0 * max(0.0, bank_stale / max(1, limit) - 0.5)
    priority = rank_bonus + momentum + recency + efficiency + cloud_bonus + late_bonus - stale_penalty

    if has_suffix and bank_stale < limit and (recent200 or late_gain or rank >= 28):
        return {'action': 'continue_bank', 'priority': priority,
                'reason': f'productive/recent bank; {bank_stale} no-gain calls < stale limit {limit}'}

    allowed = max_bank_count(rank, total_gains, late_gain)
    if generation + 1 < allowed:
        pivot_bonus = 4.0 if bank_stale >= limit else 1.0
        return {'action': 'new_bank', 'priority': priority + pivot_bonus,
                'reason': f'representation pivot {generation + 1}/{allowed - 1}; current bank stale/exhausted'}

    return {'action': 'retire', 'priority': priority,
            'reason': f'no justified continuation: bank stale={bank_stale}, limit={limit}, banks={generation + 1}/{allowed}'}


def exploration_slots(workers, exploit_summaries):
    """Guarantee breadth; only a genuinely hot >=29 cascade may use three of four slots."""
    workers = max(1, int(workers))
    if workers == 1:
        return 1
    hottest = None
    for summary in exploit_summaries:
        decision = continuation_decision(summary)
        if decision['action'] in ('continue_bank', 'new_bank'):
            key = (decision['priority'], int(summary.get('rank_lower_bound', 17)))
            if hottest is None or key > hottest[0]:
                hottest = (key, summary)
    if workers >= 4 and hottest:
        summary = hottest[1]
        calls = int(summary.get('complement_calls', 0))
        recent = _recent_gain_count(summary.get('gain_calls', []), calls, 200)
        if int(summary.get('rank_lower_bound', 17)) >= 29 and recent:
            return 1
    return max(1, (workers + 1) // 2)


def lane_for_dispatch(index):
    """Deterministic breadth schedule; score is a signal, never an exclusion gate."""
    lanes = (
        'low_score', 'low_spread', 'high_score', 'low_hash',
        'low_score', 'high_spread', 'low_spread', 'high_hash',
        'low_score', 'low_hash', 'high_score', 'low_spread',
    )
    return lanes[int(index) % len(lanes)]
