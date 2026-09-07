#!/usr/bin/env python3
"""Rank-agnostic pointed quartic search over Q.

The preserved v1 module supplies exact arithmetic and the single GMP worker.
Its old search entry points are regression controls, never called here.
Centre selection and certification of quotient gains belong to the caller.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as Q
from hashlib import sha256
from math import gcd, isfinite, isqrt, lcm
from pathlib import Path
import json
import os
import tempfile
import time

import half_lattice_pointed_sieve as base

ROOT = base.ROOT
BACKEND_NAME = "pointed_quartic_search_v1"
QuarticSearchResult = base.QuarticSearchResult
compiled_worker = base.compiled_worker


def rational(value):
    if isinstance(value, float):
        raise ValueError("supply exact rationals, not floating point coordinates")
    return Q(value)


def point(value):
    values = (value["x"], value["y"]) if isinstance(value, dict) else value
    if values is None or len(values) != 2:
        raise ValueError("expected a finite rational point")
    return tuple(map(rational, values))


def point_record(value):
    return dict(zip(("x", "y"), map(str, value)))


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def sources():
    path = Path(__file__).resolve()
    return {**base.provenance(), str(path.relative_to(ROOT)): sha256(path.read_bytes()).hexdigest()}


provenance = sources


def multiply(m, n):
    a, b, c, d = m
    e, f, g, h = n
    return a*e+b*g, a*f+b*h, c*e+d*g, c*f+d*h


def inverse(m):
    a, b, c, d = m
    determinant = a*d-b*c
    if not determinant:
        raise ValueError("singular PGL2 matrix")
    return tuple(Q(x)/determinant for x in (d, -b, -c, a))


def normalize(values):
    """Clear denominators/remove only verified square content; never factor."""
    values = tuple(map(Q, values))
    denominator = lcm(*(v.denominator for v in values))
    integers = tuple(int(v*denominator**2) for v in values)
    content = gcd(*integers)
    if not content:
        raise ValueError("zero quartic")
    root = 1
    for prime in base.SMALL_PRIMES:
        exponent = 0
        while content % prime == 0:
            content //= prime
            exponent += 1
        root *= prime**(exponent//2)
    residual = isqrt(content)
    if residual**2 == content:
        root *= residual
    scale = Q(denominator, root)
    coefficients = tuple(base.divide(v, root**2) for v in integers)
    if coefficients != tuple(v*scale**2 for v in values):
        raise ArithmeticError("quartic square scaling failed")
    return coefficients, scale


@dataclass(frozen=True)
class CoordinatePolicy:
    """z -> (a*z+b)/(c*z+d) after the selected denominator/Gauss chart.

    kind='metric', weight=16 is the calibrated MW16 metric; its sensitivity
    on another family is not asserted. 'gauss' is weight one. 'raw' uses the
    original rational chord slope. Matrix entries are exact rationals.
    """
    kind: str = "metric"
    weight: Q = Q(16)
    matrix: tuple = (1, 0, 0, 1)

    def __post_init__(self):
        if self.kind not in ("metric", "gauss", "raw"):
            raise ValueError("coordinate kind must be metric, gauss, or raw; PARI reduction is regression-only")
        weight = rational(self.weight)
        if weight <= 0:
            raise ValueError("metric weight must be positive")
        matrix = tuple(map(rational, self.matrix))
        if len(matrix) != 4:
            raise ValueError("PGL2 matrix must have four entries")
        inverse(matrix)
        # Eliminate scalar ambiguity without identifying different boxes.
        first = next(v for v in matrix if v)
        matrix = tuple(v/first for v in matrix)
        object.__setattr__(self, "weight", weight if self.kind == "metric" else Q(1))
        object.__setattr__(self, "matrix", matrix)

    @classmethod
    def parse(cls, value=None):
        if value is None:
            return cls()
        if isinstance(value, cls):
            return value
        if isinstance(value, dict):
            return cls(**value)
        parts = value.split(":")
        kind = parts.pop(0)
        weight = rational(parts.pop(0)) if kind == "metric" and parts else Q(16 if kind == "metric" else 1)
        matrix = (1, 0, 0, 1)
        if parts:
            p, q, r = map(rational, parts.pop(0).split(","))
            if p <= 0 or q <= 0:
                raise ValueError("rational slope scales must be positive")
            matrix = (p, r, 0, q)
        if parts:
            raise ValueError("unconsumed coordinate specification")
        return cls(kind, weight, matrix)

    def record(self):
        return {"kind": self.kind, "weight": str(self.weight), "matrix": list(map(str, self.matrix))}

    def horizontal(self, chart):
        if self.kind == "raw":
            den, k, u = chart.denominator, chart.shift, chart.curve_scale
            m = inverse(multiply((Q(den)/u, Q(k, den)/u, 0, 1), chart.matrix))
        else:
            den, k = chart.denominator, chart.shift
            aa = chart.model[3]*chart.curve_scale**4
            a = chart.base_point[0]*chart.curve_scale**2*den**2
            weight = (abs(a)+den**2*(isqrt(abs(int(aa)))+1))*self.weight
            m = multiply(inverse(chart.matrix), base.gauss_matrix(den, k, weight))
        m = multiply(m, self.matrix)
        coefficients, scale = normalize(base.binary_transform(chart.coefficients, m))
        return m, coefficients, scale


class PointedQuarticSearch:
    """Accept curve + subgroup + centre + coordinate policy, at any rank.

    curve: [A,B] or [a1,a2,a3,a4,a6] over Q.
    subgroup: finite points (no independence/saturation assumption).
    centre: {'coefficients': [...]} and/or {'point': [x,y]}; if both are
    supplied their equality is checked. An empty subgroup is allowed with
    an explicit point. Every output point is on the original input model.
    """
    def __init__(self, curve, subgroup, centre, coordinate_policy=None):
        curve = tuple(map(rational, curve))
        if len(curve) == 2:
            curve = (Q(0), Q(0), Q(0), *curve)
        if len(curve) != 5:
            raise ValueError("expected two or five Weierstrass coefficients")
        self.curve = curve
        a1, a2, a3, a4, a6 = curve
        c2, c4, c6 = a2+a1*a1/4, a4+a1*a3/2, a6+a3*a3/4
        self.x_shift = c2/3
        self.short_model = (Q(0), Q(0), Q(0), c4-c2*c2/3,
                            c6-c2*c4/3+2*c2**3/27)
        if 4*self.short_model[3]**3+27*self.short_model[4]**2 == 0:
            raise ValueError("singular elliptic curve")
        self.subgroup = tuple(point(p) for p in subgroup)
        if any(not self.on_curve(p) for p in self.subgroup):
            raise ValueError("subgroup point is off the input curve")
        if not isinstance(centre, dict) or not centre or set(centre)-{"coefficients", "point"}:
            raise ValueError("centre must specify coefficients and/or a point")
        self.centre_spec = {}
        computed = None
        if "coefficients" in centre:
            coefficients = tuple(map(rational, centre["coefficients"]))
            if len(coefficients) != len(self.subgroup) or any(v.denominator != 1 for v in coefficients):
                raise ValueError("centre needs one integral coefficient per subgroup generator")
            coefficients = tuple(map(int, coefficients))
            computed = base.linear_combination(self.short_model,
                tuple(self.to_short(p) for p in self.subgroup), coefficients) if coefficients else None
            if computed is None:
                raise ValueError("pointed centre is at infinity")
            computed = self.from_short(computed)
            self.centre_spec["coefficients"] = list(coefficients)
        if "point" in centre:
            explicit = point(centre["point"])
            if not self.on_curve(explicit):
                raise ValueError("centre point is off the input curve")
            if computed is not None and computed != explicit:
                raise ValueError("centre point differs from its subgroup combination")
            computed = explicit
            self.centre_spec["point"] = point_record(explicit)
        self.centre = computed
        self.policy = CoordinatePolicy.parse(coordinate_policy)
        self.chart = base.make_chart(self.short_model, self.to_short(self.centre))
        self.matrix, self.coefficients, self.ordinate_scale = self.policy.horizontal(self.chart)

    def on_curve(self, p):
        x, y = p
        a1, a2, a3, a4, a6 = self.curve
        return y*y+a1*x*y+a3*y == x**3+a2*x*x+a4*x+a6

    def to_short(self, p):
        x, y = p
        return x+self.x_shift, y+(self.curve[0]*x+self.curve[2])/2

    def from_short(self, p):
        x = p[0]-self.x_shift
        return x, p[1]-(self.curve[0]*x+self.curve[2])/2

    def map_hit(self, n, den, root):
        a, b, c, d = self.matrix
        p = self.chart.map_point(a*n+b*den, c*n+d*den, Q(root)/self.ordinate_scale)
        if p is None:
            # The two known pointed endpoints O and centre are recorded
            # separately; they are never counted as search discoveries.
            return None
        p = self.from_short(p)
        if not self.on_curve(p):
            raise ArithmeticError("recovered point is off the original curve")
        return p

    def input_record(self):
        return {"curve": list(map(str, self.curve)),
                "subgroup": [point_record(p) for p in self.subgroup],
                "centre": self.centre_spec, "coordinate_policy": self.policy.record()}

    def chart_record(self):
        return {"input": self.input_record(), "base_point": point_record(self.centre),
                "short_model": list(map(str, self.short_model)),
                "short_model_x_shift": str(self.x_shift), "pointed_chart": self.chart.record(),
                "horizontal_matrix": list(map(str, self.matrix)),
                "ordinate_scale": str(self.ordinate_scale), "coefficients": list(map(str, self.coefficients)),
                "maximum_coefficient_bits": max(abs(v).bit_length() for v in self.coefficients),
                "known_pointed_endpoints": [None, point_record(self.centre)]}

    def verify_record(self, record):
        """Exact replay of maps, hits, and coverage markers, without resieving.

        As with all saved worker transcripts, completeness trusts the pinned
        worker execution; this method is not a second point enumeration.
        """
        if record["backend"] != BACKEND_NAME:
            raise ArithmeticError("checkpoint backend mismatch")
        if record["source_hashes"] != sources() or any(record[k] for k in
                ("hyperellminimalmodel_called", "hyperellred_called", "hyperellratpoints_called")):
            raise ArithmeticError("checkpoint implementation mismatch")
        for key, value in self.chart_record().items():
            if record[key] != value:
                raise ArithmeticError("checkpoint chart mismatch: "+key)
        h, first, last = (record[k] for k in ("height_bound", "denominator_start", "denominator_end"))
        done = record["completed_denominator"]
        if not 1 <= first <= last <= h <= 1_000_000 or not first-1 <= done <= last:
            raise ArithmeticError("invalid checkpoint coverage interval")
        expected_status = "bounded_search_complete" if done == last else "bounded_search_timeout"
        if record["status"] != expected_status or record["integer_pairs_covered"] != (done-first+1)*(2*h+1):
            raise ArithmeticError("checkpoint completion census mismatch")
        found, finite, infinity, seen = set(), 0, [], set()
        for raw in record["primitive_square_hits"]:
            n, den, root = map(int, raw)
            if (n, den) in seen:
                raise ArithmeticError("duplicate projective square hit")
            seen.add((n, den))
            if den == 0:
                infinity.append((n, den, root))
            else:
                if not -h <= n <= h or not first <= den <= done or gcd(n, den) != 1:
                    raise ArithmeticError("checkpoint hit outside completed primitive box")
                finite += 1
            if root < 0 or root*root != sum(v*n**i*den**(4-i) for i, v in enumerate(self.coefficients)):
                raise ArithmeticError("checkpoint square identity failed")
            for signed in {root, -root}:
                p = self.map_hit(n, den, signed)
                if p is not None:
                    found.add(p)
        lead = self.coefficients[4]
        expected_infinity = [(1, 0, isqrt(lead))] if lead >= 0 and isqrt(lead)**2 == lead else []
        if infinity != expected_infinity or not record["infinity_checked"] or finite != record["nonnegative_square_hits"]:
            raise ArithmeticError("checkpoint projective hit census mismatch")
        if record["finite_curve_points"] != [point_record(p) for p in sorted(found)]:
            raise ArithmeticError("checkpoint mapped points mismatch")
        return QuarticSearchResult(record, tuple(sorted(found)))

    def search(self, height, seconds, *, denominator_start=1, denominator_end=None, checkpoint_dir=None):
        last = height if denominator_end is None else denominator_end
        if any(isinstance(v, bool) or int(v) != v for v in (height, denominator_start, last)):
            raise ValueError("height and interval endpoints must be integers")
        height, first, last = map(int, (height, denominator_start, last))
        if not 1 <= first <= last <= height <= 1_000_000 or not isfinite(seconds) or seconds <= 0:
            raise ValueError("invalid bounded search budget")
        bounds = {"height_bound": height, "denominator_start": first, "denominator_end": last,
                  "timeout_seconds": seconds}
        key = {"sources": sources(), "input": self.input_record(), "bounds": bounds}
        path = None
        if checkpoint_dir is not None:
            path = Path(checkpoint_dir)/(sha256(canonical(key).encode()).hexdigest()+".json")
            if path.exists():
                saved = json.loads(path.read_text())
                record = saved["record"]
                if saved["input"] != key or sha256(canonical(record).encode()).hexdigest() != saved["record_sha256"]:
                    raise ArithmeticError("pointed checkpoint integrity failed")
                if any(record[k] != v for k, v in bounds.items()):
                    raise ArithmeticError("pointed checkpoint budget changed")
                outcome = self.verify_record(record)
                if record["status"] == "bounded_search_complete":
                    return outcome
        started = time.monotonic()
        worker, hits = base.search_box(self.coefficients, height, seconds, first, last)
        lead = self.coefficients[4]
        if lead >= 0 and isqrt(lead)**2 == lead:
            hits += ((1, 0, isqrt(lead)),)
        found = set()
        for n, den, root in hits:
            for signed in {root, -root}:
                p = self.map_hit(n, den, signed)
                if p is not None:
                    found.add(p)
        record = {"backend": BACKEND_NAME, "source_hashes": sources(), **self.chart_record(), **bounds,
                  "integer_pairs_covered": 0, "nonnegative_square_hits": 0, **worker,
                  "infinity_checked": True, "primitive_square_hits": [list(map(str, p)) for p in hits],
                  "finite_curve_points": [point_record(p) for p in sorted(found)],
                  "wall_seconds": time.monotonic()-started,
                  "hyperellminimalmodel_called": False, "hyperellred_called": False,
                  "hyperellratpoints_called": False,
                  "claim_boundary": "Finite box coverage and exact points only; no quotient rank or rank upper bound is inferred."}
        outcome = self.verify_record(record)
        if path is not None:
            path.parent.mkdir(parents=True, exist_ok=True)
            document = {"input": key, "record": record,
                        "record_sha256": sha256(canonical(record).encode()).hexdigest()}
            with tempfile.NamedTemporaryFile(mode="w", dir=path.parent, suffix=".tmp", delete=False) as stream:
                temporary = Path(stream.name)
                stream.write(canonical(document)+"\n")
            try:
                os.replace(temporary, path)
            finally:
                temporary.unlink(missing_ok=True)
        return outcome


def run_quartic_search(*, mask, representative, short_model, generic_points,
                       height_bound, timeout_seconds, stack_bytes=0,
                       coordinate_policy=None, checkpoint_dir=None):
    """Compatibility adapter for half-lattice callers; no rank constants."""
    search = PointedQuarticSearch(short_model, generic_points,
        {"coefficients": representative}, coordinate_policy)
    outcome = search.search(height_bound, timeout_seconds,
        checkpoint_dir=checkpoint_dir or ROOT/"artifacts/local/elliptic-curves/pointed-quartic-search/charts")
    record = {**outcome.record, "mask": int(mask), "hex": hex(int(mask)),
              "representative": list(map(int, representative))}
    return QuarticSearchResult(record, outcome.curve_points)


class CheckpointedBackend:
    def __init__(self, directory, coordinate_policy=None):
        self.directory = Path(directory)
        self.coordinate_policy = CoordinatePolicy.parse(coordinate_policy)
        self.sources = sources()

    def run_quartic_search(self, **kwargs):
        return run_quartic_search(**kwargs, coordinate_policy=self.coordinate_policy, checkpoint_dir=self.directory)


def checkpoint(directory, *, model, points, representative, mask, specification, height, seconds):
    """Sensitivity-runner adapter; coordinate exploration uses the same service."""
    result = run_quartic_search(mask=mask, representative=representative,
        short_model=model, generic_points=points, height_bound=height, timeout_seconds=seconds,
        coordinate_policy=specification, checkpoint_dir=directory)
    return {**result.record, "specification": specification}
