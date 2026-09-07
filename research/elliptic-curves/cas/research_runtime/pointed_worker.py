"""The active GMP worker adapter, with shared resource ownership and logs."""
from functools import lru_cache
from hashlib import sha256
from math import gcd
from pathlib import Path
import os
import subprocess

import half_lattice_pointed_sieve as base
from .store import default_store
from .supervisor import captured_run, WorkerFailure


@lru_cache(maxsize=1)
def compiled_worker():
    content=sha256(base.WORKER.read_bytes()).hexdigest()
    binary=base.ROOT/f'artifacts/local/elliptic-curves/pointed-sieve-build/sieve-{content}'
    with default_store().lock({'compiled-pointed-worker':content}):
        if not binary.is_file():
            binary.parent.mkdir(parents=True,exist_ok=True)
            temporary=binary.with_name(binary.name+f'.{os.getpid()}.tmp')
            try:
                captured_run(['g++','-O3','-std=c++17',str(base.WORKER),'-lgmpxx','-lgmp','-o',str(temporary)],
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
        process = captured_run([str(binary)],input=data,text=True,capture_output=True,
                                 check=True,timeout=seconds+5)
    except (subprocess.TimeoutExpired, WorkerFailure) as error:
        # An external interruption discards this chunk; never infer coverage
        # from a partial stdout stream. Replay it from the last saved chunk.
        return {"status":"bounded_search_timeout", "completed_denominator":first-1,
                "supervisor_failure":getattr(error,"record",{"outcome":"strict_wall_timeout"})}, ()
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
        "supervisor":process.supervision,
    }, tuple(points)

