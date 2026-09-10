"""Scheduling and comparable exposure rules; no heuristic is a rank exclusion."""
from fractions import Fraction
from functools import lru_cache
import hashlib
import math

PANEL_SIZE=12
BASE_CALLS=64


@lru_cache(maxsize=32768)
def parameter(family,index):
    """Distinct, target-free rational addresses, in expanding fixed height bands.

    Each 12-fibre block has eight small and four larger-height addresses. A
    prefix is frozen before outcomes; singular/uncertified slots are retained.
    """
    block,slot=divmod(index,PANEL_SIZE)
    band=(2+8*block,9+8*block) if slot<8 else (17+64*block,65+64*block)
    pool=[Fraction(n,d) for d in range(1,band[1]) for n in range(-band[1]+1,band[1])
          if n and band[0]<=max(abs(n),d)<band[1] and math.gcd(n,d)==1]
    pool.sort(key=lambda q:hashlib.sha256(f'parent-panel-v1/{family}/{block}/{q}'.encode()).digest())
    used={parameter(family,i) for i in range(index)}
    return next(str(q) for q in pool if str(q) not in used)


def tails(rows):
    """Detection proportions, with unresolved slots exposed, not true rank laws."""
    n=len(rows)
    good=[r for r in rows if r.get('rank_lower_bound') is not None]
    unresolved=sum(not r.get('exposure_complete',False) for r in rows)
    result={'slots':n,'certified_slots':len(good),'incomplete_exposures':unresolved}
    for threshold in (3,5,8):
        hits=sum(r['certified_jump_lower_bound']>=threshold for r in good)
        unseen=sum(not r.get('exposure_complete',False) and
                   (r.get('rank_lower_bound') is None or r['certified_jump_lower_bound']<threshold) for r in rows)
        result[str(threshold)]={'certified_hits':hits,'fraction':f'{hits}/{n}' if n else None,
            'incomplete_slot_upper_fraction':f'{hits+unseen}/{n}' if n else None}
    return result


def promote(rows):
    """Promote the parent after a complete scheduled panel, including its nulls."""
    if len(rows)<PANEL_SIZE:return False
    gains=[r.get('certified_jump_lower_bound') or 0 for r in rows]
    return max(gains,default=0)>=5 or sum(g>=1 for g in gains)>=3


def utility(fibre):
    rank=fibre['rank'];jump=rank-fibre['generic_rank'];calls=fibre.get('calls',0)
    stale=fibre.get('stale_calls',0);recent=fibre.get('last_gain',0)
    near={27:2,28:4,29:7,30:11,31:18}.get(rank,0)
    return (1+near+4*recent+40*jump/max(32,calls))*math.exp(-stale/96)


def continuation(fibre):
    if fibre.get('rank') is None or fibre['rank']>=32:return False
    if fibre['rank']<=fibre['generic_rank']:return False
    limit=256 if fibre['rank']>=28 else 128
    return fibre.get('stale_calls',0)<limit and fibre.get('empty_batches',0)<2


def exploit_allowance(fibre):
    return min(256,BASE_CALLS+32*fibre.get('batches',0))


def guardian_should_restart(stop_requested,returncode):
    # Clean bounded runs and crashes both renew. Only explicit stop ends search.
    return not stop_requested
