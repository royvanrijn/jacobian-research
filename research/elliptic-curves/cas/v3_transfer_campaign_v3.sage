#!/usr/bin/env sage-python
"""V3 transfer wrapper with deterministic seed replay and isolated supervision.

This wrapper preserves the failed v1/v2 transfer attempts. New work remains in
v3-transfer-11952-v3. Seed intake is certificate-driven and case-bound: after
bind_job validates the frozen seed/proof files, the worker reconstructs exactly
that subgroup from the certificate primes and no-rational-2-torsion witness.
The first raw_state call for the frozen curve+point list must consume that exact
certified state; it may not silently fall back to greedy admission.
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

_prev_sources = base.driver_sources
def driver_sources():
    rows = dict(_prev_sources())
    for p in (CAS/'v3_transfer_campaign_v2.sage', SELF):
        rows[str(p.relative_to(base.ROOT))] = base.sha(p)
    return rows
base.driver_sources = driver_sources

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

import research_runtime.search_state as search_state
from research_runtime.arithmetic import ArithmeticContext, CurveModel, rationals
from research_runtime.mw_state import MWState
_fallback_raw_state = search_state.raw_state
_active_seed = None

def _normal_curve(curve):
    curve = tuple(map(F, curve))
    return (F(0),F(0),F(0),*curve) if len(curve)==2 else curve

def _seed_identity(curve, points):
    return (_normal_curve(curve), tuple(rationals(p) for p in points))

def _build_certified_state(identity, proof, cache):
    model = CurveModel(identity[0])
    context = ArithmeticContext.for_search(model)
    state = MWState.empty(context, cache=cache, primes=proof['primes'],
                          no_two_torsion_prime=proof['torsion_prime'])
    for index, point in enumerate(identity[1]):
        before = state.rank
        state = state.adjoin(point, cache=cache, extra_primes=())
        base.require(state.rank == before + 1,
            f'Frozen seed proof failed at column {index}; certificate does not certify this basis prefix')
    base.require(tuple(state.basis) == identity[1], 'Certificate-driven seed replay changed exact basis/order')
    return state

def certified_raw_state(curve, points, *, cache=None, prime_bound=1000):
    global _active_seed
    identity = _seed_identity(curve, points)
    active = _active_seed
    if active is not None and identity == active['identity']:
        cache = cache or search_state.reduction_cache()
        state = _build_certified_state(identity, active['proof'], cache)
        base.require(state.rank == active['initial_rank'],
                     'Frozen seed certificate replay rank differs from declared initial rank')
        active['consumed'] += 1
        return state
    return _fallback_raw_state(curve, points, cache=cache, prime_bound=prime_bound)
search_state.raw_state = certified_raw_state

_original_bind_job = base.bind_job
def bind_job(case):
    global _active_seed
    engine, folder, policy = _original_bind_job(case)
    seed = base.read_json(folder/'seed-input.json')
    proof = base.read_json(folder/'seed-proof.json')
    primes = tuple(int(row['prime']) for row in proof['signatures'])
    base.require(bool(primes) and len(set(primes)) == len(primes), 'Malformed frozen seed certificate prime list')
    torsion = int(proof['no_rational_2_torsion_prime'])
    identity = _seed_identity(seed['curve'], seed['points'])
    frozen = {'primes':primes, 'torsion_prime':torsion}
    from research_runtime.memory_store import MemoryFactStore
    from research_runtime.quotient_only_reduction import QuotientOnlyReductionCache as Cache
    replay = _build_certified_state(identity, frozen, Cache(MemoryFactStore()))
    base.require(replay.rank == policy['initial_rank'] and tuple(replay.basis) == identity[1],
                 'Frozen seed certificate does not replay to declared initial rank')
    _active_seed = {'case':case, 'identity':identity, 'proof':frozen,
                    'initial_rank':policy['initial_rank'], 'consumed':0}
    return engine, folder, policy
base.bind_job = bind_job

if __name__ == '__main__':
    base.main()
