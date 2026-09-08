"""Raw-search MWState admission, without a minimal model or number field."""
from dataclasses import replace
from fractions import Fraction
import threading

from .arithmetic import ArithmeticContext, CurveModel, rationals
from .finite_reduction import ReductionCache
from .mw_state import MWState, PointObservation
from .store import default_store

_cache = None
_states = {}
_state_locks = {}


def reduction_cache():
    global _cache
    if _cache is None:
        _cache = ReductionCache(default_store())
    return _cache


def _build_raw_state(curve, points, *, cache=None, prime_bound=1000):
    """Preserve input generators; certify their independent subset incrementally.

    Dependence of an arbitrary input list is permitted. Its cardinality is
    never called a rank. The exact generator list is retained as observations;
    independent ``state.basis`` is a separate reusable invariant.
    """
    from mod2_reduction_independence import _is_prime
    cache = cache or reduction_cache()
    curve = tuple(curve)
    if len(curve) == 2:
        curve = (0,0,0,*curve)
    model = CurveModel(curve)
    points = tuple(rationals((p['x'],p['y']) if isinstance(p,dict) else p) for p in points)
    if any(not model.contains(p) for p in points):
        raise ValueError('input subgroup contains an off-curve point')
    key = (id(cache),model.key,points,prime_bound)
    if key in _states:
        return _states[key]
    context = ArithmeticContext.for_search(model)
    primes = tuple(p for p in range(3,prime_bound+1) if _is_prime(p))
    polynomial = tuple(map(Fraction,model.two_division_polynomial))
    torsion_prime = next((p for p in primes if all(c.denominator%p for c in polynomial)
        and not any(sum(c.numerator*pow(c.denominator,-1,p)*x**i for i,c in enumerate(polynomial))%p==0
                    for x in range(p))),None)
    # Reuse a previously admitted prefix when an adaptive caller adds a point.
    previous = next((_states[(id(cache),model.key,points[:n],prime_bound)] for n in range(len(points)-1,0,-1)
                    if (id(cache),model.key,points[:n],prime_bound) in _states),None)
    initial_count = len([r for r in previous.observations if r.status.startswith('INPUT_GENERATOR:')]) if previous else 0
    state = previous or MWState.empty(context,cache=cache,primes=(),no_two_torsion_prime=torsion_prime)
    for i, point in enumerate(points[initial_count:],initial_count):
        state = state.adjoin(point,cache=cache,extra_primes=primes)
        state = replace(state,observations=(*state.observations,PointObservation(point,f'INPUT_GENERATOR:{i}',None)))
    _states[key] = state
    return state


def raw_state(curve, points, *, cache=None, prime_bound=1000):
    cache=cache or reduction_cache()
    curve=tuple(curve)
    points=tuple(rationals((p['x'],p['y']) if isinstance(p,dict) else p) for p in points)
    key=(id(cache),rationals(curve),points,prime_bound)
    lock=_state_locks.setdefault(key,threading.RLock())
    with lock:
        return _build_raw_state(curve,points,cache=cache,prime_bound=prime_bound)


def input_generators(state):
    rows = [r for r in state.observations if r.status.startswith('INPUT_GENERATOR:')]
    return tuple(r.point for r in sorted(rows,key=lambda r:int(r.status.split(':')[1]))) if rows else state.basis
