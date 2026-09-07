#!/usr/bin/env sage-python
"""V3 transfer wrapper with deterministic seed replay and isolated supervision.

This wrapper preserves the failed v1/v2 transfer attempts. New work remains in
v3-transfer-11952-v3. Compared with the earlier wrapper, seed intake is now
certificate-driven: after bind_job has validated the frozen seed/proof files,
raw_state reconstructs the exact seed basis using the exact finite-reduction
primes and no-rational-2-torsion witness from seed-proof.json. It no longer asks
the greedy incremental admission path to rediscover an already-certified basis.
"""
from importlib.machinery import SourceFileLoader
from pathlib import Path
from fractions import Fraction as F
import os
import sys

SELF = Path(__file__).resolve()
CAS = SELF.parent
v2 = SourceFileLoader('v3_transfer_campaign_v2_preserved', str(CAS/'v3_transfer_campaign_v2.sage')).load_module()
base = v2.base
base.D = base.ROOT/'artifacts/local/elliptic-curves/v3-transfer-11952-v3'
base.__file__ = str(SELF)

# Bind v1, v2 and this wrapper into every newly prepared protocol.
_prev_sources = base.driver_sources
def driver_sources():
    rows = dict(_prev_sources())
    for p in (CAS/'v3_transfer_campaign_v2.sage', SELF):
        rows[str(p.relative_to(base.ROOT))] = base.sha(p)
    return rows
base.driver_sources = driver_sources

# Preparation supervision must belong to the active campaign, not v1.
_original_supervise = base.supervise
def supervise(action, case=None, v3_replay=None):
    if action != 'prepare-worker':
        return _original_supervise(action, case=case, v3_replay=v3_replay)
    from research_runtime.supervisor import Limits, run
    logdir = base.D.parent/(base.D.name+'-preparation')
    logdir.mkdir(parents=True, exist_ok=True)
    stem = 'prepare'
    base.require(not (logdir/(stem+'.supervisor.json')).exists(),
                 'Preserve previous supervised preparation attempt. Use repair-resume; do not overwrite.')
    command = [sys.executable, str(SELF), action]
    if v3_replay:
        command += ['--v3-replay', str(v3_replay.resolve())]
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1',
               MKL_NUM_THREADS='1', V3_TRANSFER_SUPERVISED='1')
    report = run(command,
        limits=Limits(base.LIMITS['prepare_wall_seconds'], base.LIMITS['rss_bytes']),
        cwd=base.ROOT, env=env, log_path=logdir/(stem+'.log'),
        checkpoint_path=logdir/(stem+'.supervisor.json'))
    base.require(report['outcome'] == 'completed',
                 f"prepare-worker stopped: {report['outcome']}; inspect {logdir/'prepare.log'}")
base.supervise = supervise

# Deterministic theorem-backed seed replay. bind_job registers only a seed whose
# exact files and protocol hashes have already passed the preserved v1 checks.
import research_runtime.search_state as search_state
from research_runtime.arithmetic import ArithmeticContext, CurveModel, rationals
from research_runtime.mw_state import MWState
_fallback_raw_state = search_state.raw_state
_certified_seeds = {}

def _normal_curve(curve):
    curve = tuple(map(F, curve))
    return (F(0),F(0),F(0),*curve) if len(curve)==2 else curve

def _seed_key(curve, points, prime_bound):
    return (_normal_curve(curve), tuple(rationals(p) for p in points), int(prime_bound))

def certified_raw_state(curve, points, *, cache=None, prime_bound=1000):
    key = _seed_key(curve, points, prime_bound)
    proof = _certified_seeds.get(key)
    if proof is None:
        return _fallback_raw_state(curve, points, cache=cache, prime_bound=prime_bound)
    cache = cache or search_state.reduction_cache()
    model = CurveModel(key[0])
    context = ArithmeticContext.for_search(model)
    state = MWState.empty(context, cache=cache, primes=proof['primes'],
                          no_two_torsion_prime=proof['torsion_prime'])
    for index, point in enumerate(key[1]):
        before = state.rank
        state = state.adjoin(point, cache=cache, extra_primes=())
        base.require(state.rank == before + 1,
            f'Frozen seed proof failed at column {index}; exact certificate primes do not certify this basis prefix')
    base.require(tuple(state.basis) == key[1], 'Certificate-driven seed replay changed exact basis/order')
    return state
search_state.raw_state = certified_raw_state

_original_bind_job = base.bind_job
def bind_job(case):
    engine, folder, policy = _original_bind_job(case)
    seed = base.read_json(folder/'seed-input.json')
    proof = base.read_json(folder/'seed-proof.json')
    primes = tuple(int(row['prime']) for row in proof['signatures'])
    base.require(bool(primes) and len(set(primes)) == len(primes), 'Malformed frozen seed certificate prime list')
    torsion = int(proof['no_rational_2_torsion_prime'])
    key = _seed_key(seed['curve'], seed['points'], 1000)
    _certified_seeds[key] = {'primes':primes, 'torsion_prime':torsion}
    # Replay immediately at the trust boundary. This is intentionally stronger
    # than merely checking the JSON digest and catches a stale/mismatched proof
    # before any landscape or chart is built.
    from research_runtime.memory_store import MemoryFactStore
    from research_runtime.quotient_only_reduction import QuotientOnlyReductionCache as Cache
    replay = certified_raw_state(key[0], key[1], cache=Cache(MemoryFactStore()), prime_bound=1000)
    base.require(replay.rank == policy['initial_rank'] and tuple(replay.basis) == key[1],
                 'Frozen seed certificate does not replay to declared initial rank')
    return engine, folder, policy
base.bind_job = bind_job

if __name__ == '__main__':
    base.main()
