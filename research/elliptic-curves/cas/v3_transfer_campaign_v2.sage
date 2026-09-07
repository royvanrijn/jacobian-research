#!/usr/bin/env sage-python
"""V2 wrapper for the failed first 11952 transfer campaign.

Preserves v3_transfer_campaign.sage and its failed v1 artifacts.  The only
behavioral change is a local guard against Python object-id reuse in
research_runtime.search_state's per-cache memoization during preparation.
All new artifacts live under v3-transfer-11952-v2.
"""
from importlib.machinery import SourceFileLoader
from pathlib import Path
import sys

SELF = Path(__file__).resolve()
CAS = SELF.parent
base = SourceFileLoader('v3_transfer_campaign_v1_preserved', str(CAS/'v3_transfer_campaign.sage')).load_module()
base.D = base.ROOT/'artifacts/local/elliptic-curves/v3-transfer-11952-v2'
base.__file__ = str(SELF)  # supervised children must re-enter this wrapper

# Bind both this wrapper and the preserved v1 implementation into every new
# protocol.  No v1 artifact or source is rewritten.
_original_sources = base.driver_sources
def driver_sources():
    rows = dict(_original_sources())
    old = CAS/'v3_transfer_campaign.sage'
    rows[str(old.relative_to(base.ROOT))] = base.sha(old)
    rows[str(SELF.relative_to(base.ROOT))] = base.sha(SELF)
    return rows
base.driver_sources = driver_sources

# search_state historically keyed memoized states by id(cache), without
# retaining the cache object.  A dead cache's id can therefore be reused in a
# long preparation process.  Keep every explicitly supplied cache alive for
# this campaign and purge a stale numeric-id namespace before first use.
import research_runtime.search_state as search_state
_raw_state = search_state.raw_state
_live_cache_ids = {}
_cache_refs = []
def safe_raw_state(curve, points, *, cache=None, prime_bound=1000):
    if cache is not None:
        cid = id(cache)
        known = _live_cache_ids.get(cid)
        if known is not cache:
            for key in [k for k in search_state._states if k[0] == cid]:
                search_state._states.pop(key, None)
                search_state._state_locks.pop(key, None)
            _live_cache_ids[cid] = cache
            _cache_refs.append(cache)
    return _raw_state(curve, points, cache=cache, prime_bound=prime_bound)
search_state.raw_state = safe_raw_state

if __name__ == '__main__':
    base.main()
