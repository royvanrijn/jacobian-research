"""Durable receipts and adapters to existing, independently replayed search workers."""
from __future__ import annotations
from fractions import Fraction as F
import hashlib
import json
import os
from pathlib import Path
import tempfile

from broad_rank_policy import require


def read(path):
    return json.loads(Path(path).read_text())


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''): h.update(block)
    return h.hexdigest()


def write(path, obj, immutable=True):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    data = (json.dumps(obj, sort_keys=True, indent=2)+'\n').encode()
    if immutable and path.exists():
        require(path.read_bytes() == data, 'immutable receipt differs: '+str(path))
        return
    fd, tmp = tempfile.mkstemp(prefix=path.name+'.', suffix='.tmp', dir=path.parent)
    try:
        with os.fdopen(fd, 'wb') as f:
            f.write(data); f.flush(); os.fsync(f.fileno())
        os.replace(tmp, path)
        directory = os.open(path.parent, os.O_RDONLY)
        try: os.fsync(directory)
        finally: os.close(directory)
    finally:
        if os.path.exists(tmp): os.unlink(tmp)


def token(pid):
    try:
        # Fields start after ')'; field 22 is process start time. Reject zombies.
        fields = Path(f'/proc/{int(pid)}/stat').read_text().rsplit(')',1)[1].split()
        return fields[19] if fields[0] != 'Z' else None
    except (OSError, ValueError, IndexError):
        return None


def alive(lease):
    return bool(lease and lease.get('token') and token(lease.get('pid')) == lease['token'])


def guard(folder):
    folder = Path(folder); manifest = read(folder/'manifest.json')
    require(sha(folder/'plan.json') == manifest['plan_sha256'], 'plan changed')
    plan = read(folder/'plan.json'); root = Path(plan['root'])
    for name, digest in manifest['files'].items():
        require(sha(root/name) == digest, 'frozen input/source changed: '+name)
    for name, digest in manifest['executables'].items():
        require(sha(name) == digest, 'executable changed: '+name)
    return plan, root


def completed(job):
    """A checksum receipt is not a mathematical proof; backend proof replay is separate."""
    job = Path(job)
    if not (job/'seal.json').exists(): return None
    seal = read(job/'seal.json')
    require(sha(job/'dispatch.json') == seal['dispatch_sha256'], 'dispatch changed')
    for name, digest in seal['files'].items():
        require(sha(job/name) == digest, 'completed evidence changed: '+str(job/name))
    if not seal['completed']: return seal, None
    require('result.json' in seal['files'], 'completed task has no sealed result')
    return seal, read(job/'result.json')


def normalize_parent(source, gram, family):
    if 'A_coefficients_low_to_high' in source:
        result = {k:source[k] for k in ('A_coefficients_low_to_high','B_coefficients_low_to_high','sections')}
    else:
        def polynomial(record):
            den = list(map(F, record['denominator']))
            require(den and den[0] and not any(den[1:]), 'nonpolynomial parent invariant')
            return [str(F(c)/den[0]) for c in record['numerator']]
        inv = [polynomial(a) for a in source['a_invariants']]
        require(not any(F(c) for a in inv[:3] for c in a), 'expected short parent')
        def rat(r):
            return {'numerator_coefficients_low_to_high':list(map(str,r['numerator'])),
                    'denominator_coefficients_low_to_high':list(map(str,r['denominator']))}
        result = {'A_coefficients_low_to_high':inv[3], 'B_coefficients_low_to_high':inv[4],
                  'sections':[{'basis_index':i, 'X':rat(p[0]), 'Y':rat(p[1])}
                              for i,p in enumerate(source['basis_weierstrass_coordinates'])]}
    require(len(result['sections']) == 17 and len(gram) == 17 and all(len(r)==17 for r in gram), 'not an M17 parent')
    result.update(family=family, generic_rank_lower_bound=17, generic_height_gram=gram)
    return result


def absorb(previous, result, row, spec, job, root):
    """Only backend results tied to replayed packet, request and exact model advance rank."""
    native = row['backend'] == 'native'
    expected = 'PASS_CERTIFIED_SEARCH' if native else 'PASS_CERTIFIED_PARENT_EVALUATION'
    require(result['status'] == expected, 'not a certified search result')
    require(result['request_sha256'] == sha(job/'request.json'), 'backend request changed')
    require(result['family'] == row['family'] and F(result['parameter']) == F(row['parameter']), 'wrong fibre')
    packet = read(job/'packet.json'); proof = read(job/('packet-verified.json' if native else 'verified.json'))
    digest = sha(job/'packet.json')
    require(digest == result['packet_sha256'] == proof['packet_sha256'], 'packet/replay binding differs')
    require(proof['status'] == 'PASS_TWO_FINITE_IMPLEMENTATIONS', 'independent replay missing')
    rank = result['rank_lower_bound']; old_rank = previous['rank'] if previous else 17
    require(type(rank) is int and rank == packet['rank_lower_bound'] == len(packet['points']) and rank >= old_rank,
            'rank/prefix regression')
    require(list(map(F,packet['curve'])) == list(map(F,row['model'])), 'packet equation changed')
    if previous:
        old = read(root/previous['packet'])
        require(packet['points'][:len(old['points'])] == old['points'], 'continuation lost certified points')
    calls = result['calls']; maximum = spec['allowance']+(98 if spec['kind']=='fresh' or spec['revival'] else 0)
    if not native: maximum = spec['generic_allowance']
    require(type(calls) is int and 0 <= calls <= maximum, 'point-call budget violated')
    events = result.get('gain_timeline', [])
    require(all(type(e['call']) is int and 0 <= e['call'] <= calls for e in events), 'gain origin outside batch')
    # Existing full-cloud events need not be chronological by discovery rank;
    # they still have to account for a continuous sequence of certified ranks.
    cursor = old_rank
    for event in sorted(events, key=lambda e:(e['before'],e['after'])):
        require(event['before'] == cursor and event['after'] > cursor, 'broken gain accounting')
        cursor = event['after']
    require(cursor == rank, 'unaccounted rank gain')
    if native:
        censored = result['point_timeouts'] or result['map_timeouts'] or result['unresolved_cloud']
    else:
        censored = any(s['point_timeouts'] or s['map_timeouts'] for s in result['segments'])
    status = 'CENSORED_WITH_CERTIFIED_LOWER_BOUND' if censored else 'CERTIFIED'
    if not native and not result['exposure_complete'] and not censored:
        status = 'POLICY_EXHAUSTED'
    if calls == 0 and rank < 32: status = 'UNKNOWN_NO_PROGRESS'
    old_calls = previous['calls'] if previous else 0
    last = max((e['call'] for e in events), default=None)
    stale = calls-last if last is not None else (previous['stale_calls'] if previous else 0)+calls
    return {'status':status, 'rank':rank, 'calls':old_calls+calls, 'stale_calls':stale,
            'round':previous['round']+1 if previous else 0,
            'packet':str((job/'packet.json').relative_to(root)), 'packet_sha256':digest,
            'head':result.get('head'), 'rescue_used':bool(spec['revival'] or (previous and previous.get('rescue_used'))),
            'history':(previous['history'] if previous else [])+
                      [{**e,'call':old_calls+e['call']} for e in events]}
