#!/usr/bin/env python3
"""Factorization-free pointed-quartic charts and exact modular slope search.

For an integral E and Q=(a/d^2,b/d^3), put k=-b/a (mod d^2).
Then f(d*s+k/d)/d^2 is integral. Its five coefficients are constructed by
checked divisions. Exact Gauss reduction balances the slope lattice; no
hyperelliptic minimization/reduction, real roots, or large factorization occurs.
The C++ worker sieves a declared projective box before GMP square tests.
Every hit is independently checked and transported over QQ in Python.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as Q
from functools import lru_cache
from hashlib import sha256
from math import comb, gcd, isqrt, lcm
from pathlib import Path
import os
import json
import subprocess
import time

from alternate_quartic_covers import alternate_cover, point_on_short_curve, short_add

ROOT = Path(__file__).resolve().parents[2]
WORKER = Path(__file__).with_name("pointed_quartic_sieve.cpp")
BACKEND_NAME = "pointed_denominator_lattice_modular_sieve_v1"
SMALL_PRIMES = tuple(p for p in range(2,301) if all(p % q for q in range(2,isqrt(p)+1)))


def integer_root(n, k):
    if n < 0 or k < 1:
        raise ValueError("expected nonnegative integer and positive root degree")
    if n < 2:
        return n
    x = 1 << ((n.bit_length()+k-1)//k)
    while True:
        y = ((k-1)*x+n//x**(k-1))//k
        if y >= x:
            return x
        x = y


def divide(n, d):
    quotient, remainder = divmod(n, d)
    if remainder:
        raise ArithmeticError("pointed integral-chart divisibility failed")
    return quotient


def denominator_scale(n, degree):
    scale = 1
    for p in SMALL_PRIMES:
        v = 0
        while n % p == 0:
            n //= p
            v += 1
        scale *= p**((v+degree-1)//degree)
    r = integer_root(n, degree)
    # A failed perfect-power test is harmless: the remaining denominator is
    # itself a valid clearing scale. Never factor it or guess its valuations.
    return scale*(r if r**degree == n else n)


@lru_cache(maxsize=128)
def integral_short_scale(model):
    if len(model) != 5 or any(model[:3]):
        raise ValueError("expected a short Weierstrass model")
    a, b = map(Q, model[3:])
    if 4*a**3+27*b**2 == 0:
        raise ValueError("singular elliptic curve")
    # Combine valuations before taking a root: cancellation in A alone must
    # not hide the common scale still visible in B (or conversely).
    d = denominator_scale(lcm(a.denominator**3,b.denominator**2),12)
    aa, bb = a*d**4, b*d**6
    if aa.denominator != 1 or bb.denominator != 1:
        raise ArithmeticError("curve denominator clearing failed")
    aa, bb = int(aa), int(bb)
    content = gcd(aa**3, bb**2)
    n = 1
    for p in SMALL_PRIMES:
        v = 0
        while content % p == 0:
            content //= p
            v += 1
        n *= p**(v//12)
    r = integer_root(content, 12)
    if r**12 == content:
        n *= r
    integral = (Q(0), Q(0), Q(0), Q(divide(aa, n**4)), Q(divide(bb, n**6)))
    return Q(d, n), integral


def binary_transform(coefficients, matrix):
    """Coefficients of F(a*s+b,c*s+d), in ascending powers of s."""
    a, b, c, d = matrix
    out = [0]*5
    for i, value in enumerate(coefficients):
        for j in range(i+1):
            for k in range(5-i):
                out[j+k] += value*comb(i,j)*a**j*b**(i-j)*comb(4-i,k)*c**k*d**(4-i-k)
    return tuple(out)


def invariants(f):
    e, d, c, b, a = f
    return 12*a*e-3*b*d+c*c, 72*a*c*e+9*b*c*d-27*a*d*d-27*b*b*e-2*c**3


def gauss_matrix(d, k, weight_squared):
    """Reduce <(d^2,0),(k,1)> with exact metric diag(1,weight_squared)."""
    left, right = (d*d, 0), (k, 1)
    u, v = (1, 0), (0, 1)
    def dot(a,b):
        return a[0]*b[0]+weight_squared*a[1]*b[1]
    while True:
        if dot(right,right) < dot(left,left):
            left,right,u,v = right,left,v,u
        norm, product = dot(left,left), dot(left,right)
        if 2*abs(product) <= norm:
            break
        multiple = (2*product+norm)//(2*norm)
        right = tuple(a-multiple*b for a,b in zip(right,left))
        v = tuple(a-multiple*b for a,b in zip(v,u))
    answer = (u[0],v[0],u[1],v[1])
    if abs(answer[0]*answer[3]-answer[1]*answer[2]) != 1:
        raise ArithmeticError("slope-lattice basis change is not unimodular")
    return answer


@dataclass(frozen=True)
class PointedChart:
    model: tuple
    base_point: tuple
    curve_scale: Q
    denominator: int
    shift: int
    matrix: tuple
    coefficients: tuple

    def map_point(self, numerator, denominator, root):
        """Map [n:d:root] with root^2=F(n,d); d=0 is allowed."""
        a,b,c,d = self.matrix
        upper, lower = a*numerator+b*denominator, c*numerator+d*denominator
        if lower == 0:
            # These are the two rational points over the original slope infinity.
            return None
        slope = Q(self.denominator**2*upper+self.shift*lower,
                  self.denominator*lower)/self.curve_scale
        ordinate = Q(self.denominator*root, lower**2)/self.curve_scale**2
        return alternate_cover(self.model, self.base_point).cover_point_to_curve((slope,ordinate))

    def record(self):
        return {
            "curve_coordinate_scale": str(self.curve_scale),
            "point_denominator_root": str(self.denominator),
            "shift_mod_denominator_squared": str(self.shift),
            "unimodular_horizontal_matrix": [str(v) for v in self.matrix],
            "integral_quartic_coefficients_ascending": [str(v) for v in self.coefficients],
            "maximum_coefficient_bits": max(abs(v).bit_length() for v in self.coefficients),
            "invariants_checked": True,
        }


def make_chart(model, base_point):
    model, base_point = tuple(map(Q,model)), tuple(map(Q,base_point))
    if not point_on_short_curve(model,base_point):
        raise ValueError("base point is not on the supplied curve")
    scale, integral = integral_short_scale(model)
    x, y = base_point[0]*scale**2, base_point[1]*scale**3
    d = isqrt(x.denominator)
    if d*d != x.denominator or (y*d**3).denominator != 1:
        raise ArithmeticError("integral curve point lacks square/cube denominators")
    a, b, aa = int(x*d*d), int(y*d**3), int(integral[3])
    k = (-b*pow(a,-1,d*d)) % (d*d) if d > 1 else 0
    if 2*k > d*d:
        k -= d*d
    f = (divide(k**4-6*a*k*k-8*b*k-3*a*a-4*aa*d**4,d**6),
         divide(4*(k**3-3*a*k-2*b),d**4),
         divide(6*(k*k-a),d*d), 4*k, d*d)
    matrix = gauss_matrix(d,k,abs(a)+d*d*(isqrt(abs(aa))+1))
    f = binary_transform(f,matrix)
    if invariants(f) != (-48*integral[3], -1728*integral[4]):
        raise ArithmeticError("pointed quartic invariants changed")
    return PointedChart(model,base_point,scale,d,k,matrix,f)


@lru_cache(maxsize=1)
def compiled_worker():
    digest = sha256(WORKER.read_bytes()).hexdigest()
    binary = ROOT / f"artifacts/local/elliptic-curves/pointed-sieve-build/sieve-{digest}"
    if not binary.is_file():
        binary.parent.mkdir(parents=True,exist_ok=True)
        temporary = binary.with_name(binary.name+f".{os.getpid()}.tmp")
        try:
            subprocess.run(["g++","-O3","-std=c++17",str(WORKER),"-lgmpxx","-lgmp","-o",str(temporary)],
                           check=True,capture_output=True,text=True,timeout=60)
            temporary.replace(binary)
        finally:
            temporary.unlink(missing_ok=True)
    return binary


def search_box(coefficients, height, seconds, first=1, last=None):
    last = height if last is None else last
    if len(coefficients) != 5 or any(int(v) != v for v in coefficients):
        raise ValueError("expected five integral quartic coefficients")
    if not 1 <= first <= last <= height <= 1_000_000 or seconds <= 0:
        raise ValueError("invalid bounded slope box")
    binary = compiled_worker()
    data = f"{height} {first} {last} {seconds}\n"+'\n'.join(map(str,coefficients))+'\n'
    try:
        process = subprocess.run([str(binary)],input=data,text=True,capture_output=True,
                                 check=True,timeout=seconds+5)
    except subprocess.TimeoutExpired:
        # An external interruption discards this chunk; never infer coverage
        # from a partial stdout stream. Replay it from the last saved chunk.
        return {"status":"bounded_search_timeout", "completed_denominator":first-1}, ()
    points, primes, done = [], None, None
    for line in process.stdout.splitlines():
        fields = line.split()
        if fields[0] == "PRIMES":
            primes = list(map(int,fields[1:]))
        elif fields[0] == "POINT" and len(fields) == 4:
            point = tuple(map(int,fields[1:]))
            n,d,r = point
            if not -height <= n <= height or not first <= d <= last or gcd(n,d) != 1 or r < 0:
                raise ArithmeticError("worker returned a point outside its primitive box")
            if r*r != sum(f*n**i*d**(4-i) for i,f in enumerate(coefficients)):
                raise ArithmeticError("worker square test failed independent replay")
            points.append(point)
        elif fields[0] == "DONE" and len(fields) == 7:
            done = fields[1:]
        else:
            raise ArithmeticError("unrecognized sieve worker output")
    if primes is None or done is None:
        raise ArithmeticError("sieve worker omitted its completion markers")
    completed, words, modular, exact, count = map(int,done[:5])
    if not first-1 <= completed <= last or count != len(points) or any(p[1] > completed for p in points):
        raise ArithmeticError("sieve worker completion census failed")
    return {
        "status":"bounded_search_complete" if completed == last else "bounded_search_timeout",
        "height_bound":height, "denominator_start":first, "denominator_end":last,
        "completed_denominator":completed,
        "integer_pairs_covered":(completed-first+1)*(2*height+1),
        "word_sieve_survivors":words, "all_prime_survivors":modular,
        "exact_square_tests":exact, "nonnegative_square_hits":count,
        "sieve_primes":primes, "worker_seconds":float(done[5]),
    }, tuple(points)


@dataclass(frozen=True)
class QuarticSearchResult:
    record: dict
    curve_points: tuple


def linear_combination_python(model, points, coefficients):
    if len(points) != len(coefficients):
        raise ValueError("point/coefficient length mismatch")
    answer = None
    for point, coefficient in zip(points,coefficients):
        if int(coefficient) != coefficient:
            raise ValueError("nonintegral lattice coefficient")
        n = int(coefficient)
        point = tuple(map(Q,point))
        if n < 0:
            point = (point[0],-point[1])
        n = abs(n)
        while n:
            if n & 1:
                answer = short_add(model,answer,point)
            n >>= 1
            if n:
                point = short_add(model,point,point)
    return answer


def linear_combination(model, points, coefficients):
    if len(points) != len(coefficients) or not points:
        raise ValueError("nonempty point/coefficient lists of equal length are required")
    if len(model) != 5 or any(model[:3]) or any(int(c) != c for c in coefficients):
        raise ValueError("expected a short model and integral lattice coefficients")
    if any(not point_on_short_curve(model,p) for p in points):
        raise ValueError("a supplied generic point is off curve")
    program = f"{len(points)} {model[3]} {model[4]}\n"+'\n'.join(
        f"{p[0]} {p[1]} {int(c)}" for p,c in zip(points,coefficients))+'\n'
    process = subprocess.run([str(compiled_worker()),"--base-point"],input=program,
                             text=True,capture_output=True,check=True,timeout=60)
    fields = process.stdout.split()
    if fields == ["INFINITY"]:
        return None
    if len(fields) != 3 or fields[0] != "BASE":
        raise ArithmeticError("GMP group worker omitted its point")
    point = (Q(fields[1]),Q(fields[2]))
    if not point_on_short_curve(model,point):
        raise ArithmeticError("GMP lattice combination left the original curve")
    return point


def run_quartic_search(*, mask, representative, short_model, generic_points,
                       height_bound, timeout_seconds, stack_bytes):
    started = time.monotonic()
    base_point = linear_combination(short_model,generic_points,representative)
    if base_point is None:
        raise ArithmeticError("nonzero chart class produced infinity")
    chart = make_chart(short_model,base_point)
    record, hits = search_box(chart.coefficients,height_bound,timeout_seconds)
    curve_points = set()
    # Include the transformed slope infinity explicitly. Its image need not
    # be one of the two old points over raw slope infinity.
    lead = chart.coefficients[4]
    infinity_hits = ()
    if lead >= 0 and isqrt(lead)**2 == lead:
        infinity_hits = ((1,0,isqrt(lead)),)
    for n,d,r in hits+infinity_hits:
        for root in sorted({r,-r}):
            point = chart.map_point(n,d,root)
            if point is not None:
                curve_points.add(point)
    record.update({
        "backend":BACKEND_NAME, "mask":int(mask),
        "representative":list(map(int,representative)),
        "base_point":{"x":str(base_point[0]),"y":str(base_point[1])},
        "chart":chart.record(), "timeout_seconds":timeout_seconds,
        "wall_seconds":time.monotonic()-started,
        "hyperellminimalmodel_called":False, "hyperellred_called":False,
        "height_coordinate":"s in the recorded exact denominator/Gauss chart",
        "infinity_checked":True,
        "primitive_square_hits":[list(map(str,p)) for p in hits+infinity_hits],
        "finite_curve_points":[{"x":str(p[0]),"y":str(p[1])} for p in sorted(curve_points)],
    })
    return QuarticSearchResult(record,tuple(sorted(curve_points)))


def provenance():
    return {str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()
            for p in (Path(__file__).resolve(),WORKER,Path(__file__).with_name("alternate_quartic_covers.py"))}


class CheckpointedBackend:
    """Content-addressed, atomic chart checkpoints; incomplete charts replay."""
    def __init__(self, directory):
        self.directory = Path(directory)
        self.directory.mkdir(parents=True,exist_ok=True)
        self.sources = provenance()

    def run_quartic_search(self, **kwargs):
        key = {"sources":self.sources,"arguments":kwargs}
        canonical = json.dumps(key,sort_keys=True,default=str,separators=(",",":"))
        path = self.directory / (sha256(canonical.encode()).hexdigest()+".json")
        if path.is_file():
            cached = json.loads(path.read_text())
            record = cached["record"]
            if cached["input"] != json.loads(canonical):
                raise ArithmeticError("chart checkpoint input mismatch")
            encoded = json.dumps(record,sort_keys=True,separators=(",",":"))
            if sha256(encoded.encode()).hexdigest() != cached["record_sha256"]:
                raise ArithmeticError("chart checkpoint checksum mismatch")
            if record["status"] == "bounded_search_complete":
                points = tuple((Q(p["x"]),Q(p["y"])) for p in record["finite_curve_points"])
                if any(not point_on_short_curve(kwargs["short_model"],p) for p in points):
                    raise ArithmeticError("cached point is off curve")
                return QuarticSearchResult(record,points)
        outcome = run_quartic_search(**kwargs)
        encoded = json.dumps(outcome.record,sort_keys=True,separators=(",",":"))
        document = {"input":json.loads(canonical),"record":outcome.record,
                    "record_sha256":sha256(encoded.encode()).hexdigest()}
        temporary = path.with_suffix(f".{os.getpid()}.tmp")
        temporary.write_text(json.dumps(document,indent=2,sort_keys=True)+'\n')
        temporary.replace(path)
        return outcome
