"""Deterministic scheduling policy. These numbers allocate work, never rank."""
from __future__ import annotations
import hashlib
import math

FAMILIES = ('074d9', '07ca9', '08234', '08f72', '103b2', '11952')
DEFAULTS = {
    'workers': 4, 'batch_calls': 100, 'fresh_share': 0.60,
    'daily_point_calls': 60000, 'daily_worker_seconds': 172800,
    'batch_growth_step': 25, 'batch_growth_every': 100, 'max_batch_calls': 300,
    'phase_seconds': 1800, 'job_seconds': 10800, 'rss_bytes': 3*1024**3,
    'min_free_gib': 20, 'min_free_inodes': 100000,
    'height': 125000, 'point_seconds': 10, 'map_seconds': 5,
    'target_rank': 32, 'max_crash_retries': 2, 'max_consecutive_failures': 4,
    'conductor_seconds': 30, 'conductor_share': 0.02,
    'revival_every': 40, 'max_revivals': 2,
}


def effective_batch_calls(config, fresh_completed):
    """Increase bounded per-job work in fixed steps as fresh fibres accumulate."""
    base = int(config['batch_calls'])
    step = int(config.get('batch_growth_step', 0))
    every = max(1, int(config.get('batch_growth_every', 100)))
    ceiling = int(config.get('max_batch_calls', base))
    return min(ceiling, base + step * (int(fresh_completed) // every))


def hashed(*parts):
    return hashlib.sha256('/'.join(map(str, parts)).encode()).hexdigest()


def stale_limit(rank):
    return 900 if rank >= 30 else 600 if rank >= 28 else 300 if rank >= 25 else 200


def utility(c):
    rank, stale, calls = c['rank'], c['stale_calls'], c['total_calls']
    recent = c.get('last_batch_gain', 0)
    speed = (rank - 17) / max(25, calls)
    near = {27: 1.5, 28: 3, 29: 5, 30: 8, 31: 12}.get(rank, 0)
    return (1 + near + 5*recent + 30*speed + min(c.get('gaining_batches', 0), 5)) * math.exp(-stale/150)


def next_bank(c):
    return c.get('bank_index', 0) + (1 if c.get('last_batch_gain', 0) == 0 or c['stale_calls'] >= 100 else 0)


def apply_result(c, result):
    if result['rank_lower_bound'] < c['rank']:
        raise ArithmeticError('a certified subgroup cannot decrease')
    if result['calls'] < 0:
        raise ArithmeticError('negative point accounting')
    out = dict(c)
    gain = result['rank_lower_bound'] - c['rank']
    offset = c['total_calls']
    events = [{**e, 'call': offset+e['call']} for e in result['gain_timeline']]
    if any(e['call'] < offset or e['call'] > offset+result['calls'] for e in events):
        raise ArithmeticError('gain lies outside the dispatched batch')
    expected=c['rank']
    for e in events:
        if e['before']!=expected or e['after']<=e['before']:
            raise ArithmeticError('broken certified gain trajectory')
        expected=e['after']
    if expected!=result['rank_lower_bound']:
        raise ArithmeticError('gain trajectory does not establish final lower bound')
    new_last = max((e['call'] for e in events), default=None)
    if gain and new_last is None:
        raise ArithmeticError('gain without originating point call')
    out.update(rank=result['rank_lower_bound'], total_calls=offset+result['calls'],
               stale_calls=(offset+result['calls']-new_last if new_last is not None
                            else c['stale_calls']+result['calls']),
               last_batch_gain=gain, gaining_batches=c.get('gaining_batches',0)+bool(gain),
               batches=c.get('batches',0)+1,
               packet=result['packet_path'], packet_sha256=result['packet_sha256'],
               head=result.get('head'), bank_index=result.get('bank_index',c.get('bank_index',0)),
               history=c.get('history',[])+events)
    out['state'] = ('TARGET_CERTIFIED' if out['rank'] >= 32 else
                    'COOLED' if out['rank'] == 17 or out['stale_calls'] >= stale_limit(out['rank']) else 'READY')
    return out


def choose_lane(fresh_seconds, exploit_seconds, config, *, has_exploit=True):
    total = fresh_seconds + exploit_seconds
    return 'fresh' if not has_exploit or total == 0 or fresh_seconds < config['fresh_share']*total else 'exploit'


def parent_masks(family, maximum_masks, bank_index, count=16):
    """Disjoint maximum blocks, then full-space samples without within-cycle repeats."""
    maximum = sorted(set(maximum_masks))
    blocks = math.ceil(len(maximum)/count)
    if bank_index < blocks:
        return maximum[bank_index*count:(bank_index+1)*count], 'maximum'
    seed = int(hashed('foundry-parents-v1',family),16)
    mul = (seed & 131071) | 1
    add = (seed >> 17) & 131071
    offset = (bank_index-blocks)*128
    masks = []
    for i in range(offset, offset+128):
        mask = (mul*i+add) % 131072
        if mask and mask not in maximum:
            masks.append(mask)
    return masks, 'sampled_full_space'
