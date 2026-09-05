"""Optional bounded point cache; durable witnesses are complete finite quotients.

Point masks are recomputed from exact rational coordinates and a replayed coset
table. They need no separate immutable fact. Existing ReductionCache and frozen
campaigns retain their original behavior. Use this in a newly frozen worker.
"""
from collections import OrderedDict
from .arithmetic import CurveModel, rationals
from .finite_reduction import ReductionCache, short_presentation


class QuotientOnlyReductionCache(ReductionCache):
    def __init__(self, store=None, *, point_cache_limit=8192):
        if type(point_cache_limit) is not int or point_cache_limit < 0:
            raise ValueError('point cache limit must be a nonnegative integer')
        super().__init__(store)
        self.point_cache_limit = point_cache_limit
        self._points = OrderedDict()

    def point_signature(self, coefficients, point, prime):
        model = CurveModel(tuple(coefficients))
        point = rationals(point)
        identity = (model.key, point, prime)
        if identity in self._points:
            self._points.move_to_end(identity)
            return self._points[identity]
        if not model.contains(point):
            raise ValueError('point does not lie on the exact curve')
        row, table, _ = self.quotient(model.coefficients, prime)
        affine = short_presentation(model, point)
        reduced = (None if any(q.denominator % prime == 0 for q in affine) else
                   tuple(q.numerator * pow(q.denominator, -1, prime) % prime for q in affine))
        result = table[reduced], row['dimension']
        self.point_evaluations += 1
        if self.point_cache_limit:
            self._points[identity] = result
            if len(self._points) > self.point_cache_limit:
                self._points.popitem(last=False)
        return result
