"""Optional prepared exact models/points; existing frozen caches remain untouched.

The inherited quotient builder still verifies every complete finite coset table.
Only deterministic rational model validation and short-coordinate transport are
reused across primes. No rank assertion or ambiguous/dependent status is cached.
"""
from collections import OrderedDict
from .arithmetic import CurveModel, rationals
from .finite_reduction import short_presentation
from .quotient_only_reduction import QuotientOnlyReductionCache


class PreparedQuotientReductionCache(QuotientOnlyReductionCache):
    def __init__(self, store=None, *, point_cache_limit=8192,
                 prepared_point_limit=8192, prepared_model_limit=64):
        super().__init__(store, point_cache_limit=point_cache_limit)
        if any(type(v) is not int or v < 1 for v in (prepared_point_limit, prepared_model_limit)):
            raise ValueError('prepared cache limits must be positive integers')
        self.prepared_point_limit = prepared_point_limit
        self.prepared_model_limit = prepared_model_limit
        self._prepared_models = OrderedDict()
        self._prepared_points = OrderedDict()
        self._verified_quotient_aliases = {}

    def _model(self, coefficients):
        raw = tuple(coefficients)
        # Only exact strings take the fast cache path. Typed numeric equality
        # must not let a previously valid input hide an invalid bool/float/etc.
        if not all(type(v) is str for v in raw):
            model = CurveModel(raw)
            raw = model.coefficients
        if raw not in self._prepared_models:
            model = CurveModel(raw)
            self._prepared_models[raw] = (model, model.key)
            if len(self._prepared_models) > self.prepared_model_limit:
                self._prepared_models.popitem(last=False)
        else:
            self._prepared_models.move_to_end(raw)
        return self._prepared_models[raw]

    def quotient(self, coefficients, prime):
        if type(prime) is not int:
            return super().quotient(coefficients, prime)
        model, key = self._model(coefficients)
        identity = (model.coefficients, prime)
        if identity not in self._verified_quotient_aliases:
            self._verified_quotient_aliases[identity] = super().quotient(model.coefficients, prime)
        return self._verified_quotient_aliases[identity]

    def point_signature(self, coefficients, point, prime):
        if type(prime) is not int:
            return super().point_signature(coefficients, point, prime)
        model, key = self._model(coefficients)
        raw = tuple(point)
        if not all(type(v) is str for v in raw):
            raw = rationals(raw)
        prepared_key = (key, raw)
        if prepared_key not in self._prepared_points:
            normalized = rationals(raw)
            if not model.contains(normalized):
                raise ValueError('point does not lie on the exact curve')
            self._prepared_points[prepared_key] = (normalized, short_presentation(model, normalized))
            if len(self._prepared_points) > self.prepared_point_limit:
                self._prepared_points.popitem(last=False)
        else:
            self._prepared_points.move_to_end(prepared_key)
        normalized, affine = self._prepared_points[prepared_key]
        identity = (key, normalized, prime)
        if identity in self._points:
            self._points.move_to_end(identity)
            return self._points[identity]
        row, table, _ = self.quotient(model.coefficients, prime)
        reduced = (None if any(q.denominator % prime == 0 for q in affine) else
                   tuple(q.numerator * pow(q.denominator, -1, prime) % prime for q in affine))
        result = table[reduced], row['dimension']
        self.point_evaluations += 1
        if self.point_cache_limit:
            self._points[identity] = result
            if len(self._points) > self.point_cache_limit:
                self._points.popitem(last=False)
        return result
