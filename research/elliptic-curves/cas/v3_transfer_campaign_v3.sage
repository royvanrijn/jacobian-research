#!/usr/bin/env sage-python
"""V3 transfer wrapper with deterministic seed replay and isolated supervision.

This wrapper preserves failed transfer attempts while keeping new work in the
v3-transfer-11952-v3 namespace. Seed intake is certificate-driven: bind_job
validates the frozen seed/proof files, replays the exact finite-reduction
certificate, then binds that certified MWState to the exact (model,basis) tuple
returned by engine.v1.seed(None). The worker dispatch is patched at the actual
run_case boundary so the already-certified seed state is used directly; only
later enlarged subgroups use the normal raw_state path.
"""
from importlib.machinery import SourceFileLoader
from pathlib import Path
from fractions import Fraction as F
import inspect
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

def _normal_points(points):
    return tuple(rationals(p) for p in points)

def _first_diff(left, right):
    n=min(len(left),len(right))
    for i in range(n):
        if left[i] != right[i]:
            return i,left[i],right[i]
    if len(left)!=len(right):
        return n,('<missing>' if n>=len(left) else left[n]),('<missing>' if n>=len(right) else right[n])
    return None

def certified_raw_state(curve, points, *, cache=None, prime_bound=1000):
    global _active_seed
    c = _normal_curve(curve)
    pts = _normal_points(points)
    if _active_seed is not None and int(prime_bound)==1000:
        expected_curve, expected_points, state = _active_seed
        if c == expected_curve and pts == expected_points:
            return state
        if c == expected_curve and len(pts) == len(expected_points):
            diff = _first_diff(pts, expected_points)
            raise ValueError('V3 initial seed identity mismatch: '+repr(diff))
    return _fallback_raw_state(curve, points, cache=cache, prime_bound=prime_bound)
search_state.raw_state = certified_raw_state

_original_bind_job = base.bind_job
def bind_job(case):
    global _active_seed
    engine, folder, policy = _original_bind_job(case)
    seed = base.read_json(folder/'seed-input.json')
    proof = base.read_json(folder/'seed-proof.json')
    frozen_curve = _normal_curve(seed['curve'])
    frozen_points = _normal_points(seed['points'])

    engine_curve, engine_points = engine.v1.seed(None)
    engine_curve = _normal_curve(engine_curve)
    engine_points = _normal_points(engine_points)
    base.require(engine_curve == frozen_curve, 'Frozen engine seed curve differs from seed-input.json')
    diff = _first_diff(engine_points, frozen_points)
    base.require(diff is None, 'Frozen engine seed points differ from seed-input.json: '+repr(diff))

    primes = tuple(int(row['prime']) for row in proof['signatures'])
    base.require(bool(primes) and len(set(primes)) == len(primes), 'Malformed frozen seed certificate prime list')
    torsion = int(proof['no_rational_2_torsion_prime'])

    from research_runtime.memory_store import MemoryFactStore
    from research_runtime.quotient_only_reduction import QuotientOnlyReductionCache as Cache
    cache = Cache(MemoryFactStore())
    model = CurveModel(frozen_curve)
    context = ArithmeticContext.for_search(model)
    state = MWState.empty(context, cache=cache, primes=primes, no_two_torsion_prime=torsion)
    for index, point in enumerate(frozen_points):
        before = state.rank
        state = state.adjoin(point, cache=cache, extra_primes=())
        base.require(state.rank == before + 1,
            f'Frozen seed certificate failed at column {index}; certificate primes do not certify this marked prefix')
    base.require(state.rank == policy['initial_rank'],
                 f'Frozen seed certificate rank {state.rank} != declared {policy["initial_rank"]}')
    diff = _first_diff(tuple(state.basis), frozen_points)
    base.require(diff is None, 'Certificate-driven marked basis changed: '+repr(diff))

    _active_seed = (engine_curve, engine_points, state)
    return engine, folder, policy
base.bind_job = bind_job


def _initial_worker_state(model, basis, policy):
    if _active_seed is None:
        raise ValueError('V3 certified seed state was not installed by bind_job')
    expected_curve, expected_points, state = _active_seed
    actual_curve = _normal_curve(model)
    actual_points = _normal_points(basis)
    if actual_curve != expected_curve:
        raise ValueError('V3 worker seed curve differs after bind_job')
    diff = _first_diff(actual_points, expected_points)
    if diff is not None:
        raise ValueError('V3 worker seed points differ after bind_job: '+repr(diff))
    if state.rank != policy['initial_rank']:
        raise ValueError(f'V3 certified state rank {state.rank} != declared {policy["initial_rank"]}')
    diff = _first_diff(tuple(state.basis), actual_points)
    if diff is not None:
        raise ValueError('V3 certified state basis differs from worker seed: '+repr(diff))
    return state

# Patch the actual dispatch target, not merely search_state.raw_state. Keep the
# original run_case body unchanged except for the redundant startup rebuild and
# generic seed equality assertion. This makes the trust boundary explicit while
# preserving every subsequent landscape/chart/replay operation verbatim.
_run_case_source = inspect.getsource(base.run_case)
_old_startup = "    state = raw_state(model,basis,cache=cache,prime_bound=1000)\n    require(state.rank == p['initial_rank'] and tuple(state.basis) == basis, 'Seed rank replay failed')\n"
_new_startup = "    state = _initial_worker_state(model,basis,p)\n"
base.require(_run_case_source.count(_old_startup) == 1,
             'V3 run_case patch anchor changed; review upstream worker before execution')
_run_case_source = _run_case_source.replace(_old_startup, _new_startup)
_worker_globals = dict(base.__dict__)
_worker_globals['bind_job'] = bind_job
_worker_globals['_initial_worker_state'] = _initial_worker_state
exec(compile(_run_case_source, str(SELF)+':patched-run-case', 'exec'), _worker_globals)
base.run_case = _worker_globals['run_case']

if __name__ == '__main__':
    base.main()
