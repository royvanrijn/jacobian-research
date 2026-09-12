"""Small shared I/O and binary algebra for the frozen transfer commissioning."""
from fractions import Fraction as F
import hashlib
import json
import os
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT/'elliptic-curves/cas'
sys.path.insert(0, str(CAS))


def read(path):
    return json.loads(Path(path).read_text())


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def clean(value):
    if isinstance(value, dict): return {str(k): clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)): return [clean(x) for x in value]
    if value is None or isinstance(value, (str, int, float, bool)): return value
    return str(value)


def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(clean(value), indent=2, sort_keys=True)+'\n')
    os.replace(temporary, path)


def indices(mask):
    while mask:
        bit = mask & -mask
        yield bit.bit_length()-1
        mask ^= bit


def xor(mask, values):
    out = 0
    for i in indices(mask): out ^= values[i]
    return out


def pack(bits):
    return sum(int(x) << i for i, x in enumerate(bits))


def reduce(value, basis):
    while value:
        pivot = value.bit_length()-1
        if pivot not in basis: return value
        value ^= basis[pivot]
    return 0


def rank(rows):
    basis = {}
    for row in rows:
        r = reduce(row, basis)
        if r: basis[r.bit_length()-1] = r
    return len(basis)


def kernel(rows):
    # Same provenance-preserving elimination as early_relation_pool.elimination.
    basis, result = {}, []
    for i, row in enumerate(rows):
        word = 1 << i
        while row:
            p = row.bit_length()-1
            if p not in basis:
                basis[p] = (row, word)
                break
            r, w = basis[p]
            row ^= r
            word ^= w
        if not row: result.append(word)
    return result


def seed(packet):
    # Integral short model; this scaling is a square in the Kummer map.
    d = F(packet['parameter']).denominator
    model = [F(0)]*3 + [F(packet['short_model'][0])*d**8,
                       F(packet['short_model'][1])*d**12]
    points = [[F(x)*d**4, F(y)*d**6] for x, y in packet['generic_points']]
    assert all(v.denominator == 1 for v in model)
    return model, points


def admission(model, points, bound=1000):
    from future_point_admission import FinitePointAdmission
    from memory_rank_certificate import checked_rank
    from mod2_reduction_independence import find_two_torsion_certificate_prime
    a = FinitePointAdmission(tuple(model), tuple(map(tuple, points)), prime_bound=bound)
    proof = checked_rank(tuple(model), tuple(map(tuple, points)), a.primes,
                         find_two_torsion_certificate_prime(tuple(model)))
    return a, proof


def guard(folder):
    """Application read boundary, not OS isolation; all code is snapshotted."""
    folder = Path(folder).resolve()
    reads = set()
    def audit(event, args):
        if event != 'open' or not isinstance(args[0], (str, bytes, os.PathLike)): return
        path = Path(os.fsdecode(args[0])).resolve()
        # No research data except this arm's packet and its own outputs.
        code = path.is_relative_to(ROOT/'elliptic-curves') and path.suffix in {'.py','.pyc','.sage','.so','.cpp','.h'}
        if 'artifacts' in path.parts and not path.is_relative_to(folder) and not code:
            raise PermissionError('transfer data allowlist rejected '+str(path))
        if path.is_relative_to(folder): reads.add(str(path.relative_to(folder)))
    sys.addaudithook(audit)
    return reads


def progress(folder, stage, **details):
    write(Path(folder)/'progress.json', dict(stage=stage, updated_unix=time.time(), **details))
    print(stage, json.dumps(clean(details)), flush=True)
