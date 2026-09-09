"""Per-curve finite column admission for newly frozen search campaigns.

Prepare complete finite quotients once; retain their ordered bit coordinates
across charts and epochs. A failed test is UNKNOWN, never rational dependence.
Every admitted basis must still receive a standalone exact rank certificate.
"""
from fractions import Fraction as F
from research_runtime.arithmetic import CurveModel
from research_runtime.finite_reduction import ReductionCache
from research_runtime.memory_store import MemoryFactStore
from mod2_reduction_independence import _primes_up_to


class FinitePointAdmission:
    def __init__(self, model, basis, *, prime_bound=1000):
        self.model = CurveModel(tuple(map(F, model)))
        if any(map(F, self.model.coefficients[:3])):
            raise ValueError('short model required')
        cache = ReductionCache(MemoryFactStore())
        self.places, self.pivots, self.points, self.columns = [], {}, [], {}
        self.known_keys = set()
        offset = 0
        for prime in _primes_up_to(prime_bound):
            if prime == 2:
                continue
            try:
                record, table, _ = cache.quotient(self.model.coefficients, prime)
            except ValueError:
                continue
            self.places.append((prime, table, offset))
            offset += record['dimension']
        self.dimension = offset
        for p in basis:
            result = self.consider(p)
            if result['status'] != 'INDEPENDENT_FINITE_COLUMN':
                raise ArithmeticError('initial seed not independent in declared finite places')

    @property
    def primes(self):
        return [prime for prime, _, _ in self.places]

    def consider(self, point):
        point = tuple(map(F, point))
        if not self.model.contains(point):
            raise ValueError('point does not lie on the exact curve')
        key = point[0], abs(point[1])
        if key in self.known_keys:
            return {'status': 'KNOWN_POINT_UP_TO_SIGN', 'rank': len(self.points)}
        if key not in self.columns:
            col = 0
            for prime, table, offset in self.places:
                reduced = None if any(x.denominator % prime == 0 for x in point) else tuple(
                    x.numerator*pow(x.denominator, -1, prime) % prime for x in point)
                col |= table[reduced] << offset
            self.columns[key] = col
        value = self.columns[key]
        while value:
            bit = (value & -value).bit_length() - 1
            if bit not in self.pivots:
                self.pivots[bit] = value
                self.points.append(point)
                self.known_keys.add(key)
                return {'status': 'INDEPENDENT_FINITE_COLUMN', 'rank': len(self.points)}
            value ^= self.pivots[bit]
        return {'status': 'UNKNOWN_FINITE_COLUMN_IN_SPAN', 'rank': len(self.points)}
