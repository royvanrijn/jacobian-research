#!/usr/bin/env python3
"""Equation-only, transactional determinant1092 intake; no Sage or point backend.

SQLite commits the compressed block, deduplication index, cursor and retention
heaps together. A killed transaction can be replayed without losing addresses.
Use the controller for supervised execution; importing this module does no IO.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
from fractions import Fraction as F
import hashlib
import heapq
import json
import math
from pathlib import Path
import sqlite3
import shutil
import subprocess
import sys
import zlib

from v3_warm_support import atomic, bindings, read, require, sha

CAS = Path(__file__).resolve().parent
ROOT = CAS.parents[1]
ART = ROOT/'artifacts/generated-results/elliptic-curves'
PARENT = ART/'det1092_reduced_parameter_chart_v1/reduced-parent.json'
CHART = ART/'det1092_reduced_parameter_chart_v1/generic-proof.json'
COVER = ART/'det1092_orbit8044_rank18_base_change_v2.json'
CALIBRATION = ROOT/'artifacts/local/elliptic-curves/adaptive-visibility-cascade-v3'
ORBITS = ART/'curve302_parent_degree2_multisection_orbits_v1.tsv'
PRIMES = [p for p in range(5, 98) if all(p % d for d in range(2, math.isqrt(p)+1))]


def packed(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':')).encode()


def digest(value):
    return hashlib.sha256(value).hexdigest()


def homogeneous(co, a, b):
    value, power = co[-1], b
    for c in reversed(co[:-1]):
        value = value*a+c*power
        power *= b
    return value


def source_bindings():
    # Deliberately overbind the local Python/Sage dependency tree. An unrelated
    # source change stops a resume instead of silently changing arithmetic.
    paths = sorted(p for p in CAS.rglob('*') if p.suffix in ('.py', '.sage', '.cpp'))
    paths += sorted((ROOT/'elliptic-curves/ecsearch').rglob('*.py'))
    return {str(p.relative_to(ROOT)): sha(p) for p in paths}


def freeze(folder, profile='production'):
    require(profile in ('smoke', 'production'), 'unknown profile')
    require(not (folder/'protocol.json').exists(), 'already frozen; use resume')
    smoke = profile == 'smoke'
    original = read(CALIBRATION/'protocol.json')
    runtime = subprocess.run([shutil.which('sage') or 'sage','-python','-c',
        "import json,sys,sage.version; from sage.all import pari; print(json.dumps(dict(sage=sage.version.version,pari=str(pari.version()),python=sys.version)))"],
        capture_output=True,text=True,timeout=30,check=True)
    policy = dict(original)
    # Retain numerical settings, replace inherited oracle paths and metadata.
    for key in ('inputs', 'oracle_boundary', 'scope'):
        policy.pop(key, None)
    policy.update(initial_rank=18, target_rank=32, calibration_only=False)
    inputs = {str(p.relative_to(ROOT)): sha(p) for p in
              (PARENT, CHART, COVER, ORBITS, CALIBRATION/'protocol.json')}
    p = dict(schema='det1092-funnel-run.v1', profile=profile,
             population=256 if smoke else 10_000_000,
             maximum_draws=2560 if smoke else 50_000_000,
             domain='det1092-funnel-'+profile+'-v1', shell_bits=[12, 16, 20, 24, 28, 32],
             block_size=64 if smoke else 4096, height_bands=list(range(7, 18)),
             retain_strong_per_band=8 if smoke else 810,
             retain_control_per_band=2 if smoke else 90,
             seed_strong_per_band=1 if smoke else 8,
             seed_control_per_band=0 if smoke else 1,
             arithmetic_cap=100 if smoke else 10000, seed_cap=2 if smoke else 100,
             primes=PRIMES, inputs=inputs, sources=source_bindings(),
             unknown_features=['strict', 'selmer', 'kummer_extra_direction', 'xi'],
             score='sum round(1e12*(2-ap)*log(p)/(p+1-ap)) at frozen good primes; UNKNOWN features have no weight',
             seed=dict(max_charts=1 if smoke else 49, parity_samples=64 if smoke else 2048,
                       height=100 if smoke else 125000,
                       seconds_per_chart=2 if smoke else 10,
                       centres='fixed SHA parities, 384-bit rounded metric, largest computed norms; no adaptive seed wave',
                       proof_prime_bound=1000, torsion_prime_bound=200),
             v3=policy,
             resources=dict(workers=1, rss_bytes=3*1024**3, storage_bytes=8*1024**3,
                            intake_seconds=120 if smoke else 14400,
                            seed_seconds=120 if smoke else 1200,
                            amplifier_seconds=7200, replay_seconds=7200,
                            campaign_seconds=600 if smoke else 172800),
             software=dict(python=sys.version, sqlite=sqlite3.sqlite_version,
                           zlib=zlib.ZLIB_VERSION),
             worker_software=json.loads(runtime.stdout),
             scope='Finite prospective intake and first-seed construction. Unknown features are not exclusions. Conic controls have a separate label. No automatic refill or new rank claim from scores.')
    atomic(folder/'protocol.json', p, immutable=True)
    return p


def protocol(folder):
    p = read(folder/'protocol.json')
    bindings(ROOT, p['inputs']); bindings(ROOT, p['sources'])
    require(p['software'] == dict(python=sys.version, sqlite=sqlite3.sqlite_version,
                                 zlib=zlib.ZLIB_VERSION), 'intake runtime changed')
    require(sha('/usr/bin/gp') == p['v3']['gp_sha256'], 'frozen GP changed')
    return p


def draw(p, counter):
    bits = p['shell_bits'][counter % len(p['shell_bits'])]
    height = 1 << bits
    h = hashlib.sha256((p['domain']+':'+str(counter)).encode()).digest()
    a = int.from_bytes(h[:8], 'big') % (2*height+1)-height
    b = int.from_bytes(h[8:16], 'big') % height+1
    g = math.gcd(a, b); a //= g; b //= g
    if not a or max(abs(a), b) <= height//2:
        return None
    return a, b


class Arithmetic:
    def __init__(self, p):
        parent = read(PARENT)
        require(all(v['denominator'] == ['1'] for v in parent['a_invariants']), 'nonpolynomial model')
        require(all(not v['numerator'] for v in parent['a_invariants'][:3]), 'nonshort parent')
        self.A, self.B = [list(map(int, v['numerator'])) for v in parent['a_invariants'][3:]]
        require((len(self.A), len(self.B)) == (9, 13), 'unexpected model degrees')
        chart = read(CHART)
        require(chart['export_sha256'] == sha(PARENT), 'chart/parent binding changed')
        self.matrix = list(map(int, chart['parameter_matrix']))
        aa, bb, cc, dd = self.matrix
        q0, q1, q2 = map(int, read(COVER)['curve_over_Q']['q_coefficients'])
        self.conic = [q0*dd*dd+q1*bb*dd+q2*bb*bb,
                      2*q0*cc*dd+q1*(aa*dd+bb*cc)+2*q2*aa*bb,
                      q0*cc*cc+q1*aa*cc+q2*aa*aa]
        self.primes = p['primes']
        self.tables = []
        self.fallback = {}
        self.weights = {prime: {trace:round(10**12*(2-trace)*math.log(prime)/(prime+1-trace))
                               for trace in range(-math.isqrt(4*prime),math.isqrt(4*prime)+1)}
                        for prime in self.primes}
        for prime in self.primes:
            table = []
            for residue in range(prime+1):
                a, b = (residue, 1) if residue < prime else (1, 0)
                A, B = homogeneous(self.A, a, b) % prime, homogeneous(self.B, a, b) % prime
                table.append(self.local(A, B, prime))
            self.tables.append(table)

    @staticmethod
    def local(A, B, p):
        if (4*A**3+27*B*B) % p == 0:
            return [None, None, int(A % p == B % p == 0)]
        trace, roots = 0, 0
        for x in range(p):
            v = (x*x*x+A*x+B) % p
            roots += v == 0
            if v:
                trace -= 1 if pow(v, (p-1)//2, p) == 1 else -1
        return [trace, roots, 0]

    def model(self, a, b):
        return homogeneous(self.A, a, b), homogeneous(self.B, a, b)

    def record(self, counter, a, b):
        A, B = self.model(a, b)
        disc = 4*A**3+27*B*B
        if not disc:
            return [counter, a, b, 'SINGULAR']
        num = 6912*A**3; g = math.gcd(num, disc)
        nb, db = abs(num//g).bit_length(), abs(disc//g).bit_length()
        local, score = [], 0
        for prime, table in zip(self.primes, self.tables):
            r = a*pow(b, -1, prime) % prime if b % prime else prime
            trace, roots, ambiguous = table[r]
            scale = 0
            if ambiguous:
                av, bv = A, B
                while av % prime**4 == bv % prime**6 == 0:
                    av //= prime**4; bv //= prime**6; scale += 1
                key = (prime, av % prime, bv % prime)
                if key not in self.fallback:
                    self.fallback[key] = self.local(av % prime, bv % prime, prime)
                trace, roots, _ = self.fallback[key]
            local.append([trace, roots, scale])
            if trace is not None:
                score += self.weights[prime][trace]
        q = homogeneous(self.conic, a, b)
        pole = self.matrix[2]*a+self.matrix[3]*b == 0
        split = 'CHART_POLE' if pole else ('SPLIT' if q >= 0 and math.isqrt(q)**2 == q else 'NONSPLIT')
        return [counter, a, b, 'SMOOTH', nb, db, max(nb, db)//64, score, local, split]

    def candidate(self, row, role):
        i, a, b = row[:3]; A, B = self.model(a, b)
        aa, bb, cc, dd = self.matrix
        den = cc*a+dd*b
        return dict(id=f'funnel-{i:09d}', parameter=str(F(a, b)),
                    original_parameter=str(F(aa*a+bb*b, den)) if den else None,
                    family='det1092-reduced', model=['0', '0', '0', str(A), str(B)],
                    j_numerator_bits=row[4], j_denominator_bits=row[5], height_bin=row[6],
                    score_units=row[7], local_fingerprint=row[8], conic_splitting=row[9],
                    role=role, strict='UNKNOWN', selmer='UNKNOWN', kummer_extra_direction='UNKNOWN', xi='UNKNOWN')


def retain(state, row, p):
    if row[3] != 'SMOOTH' or row[6] not in p['height_bands']:
        return
    band = str(row[6]); counter, a, b = row[:3]
    control = int(digest((p['domain']+':control:'+str(counter)).encode()), 16)
    for role, key, size in [('strong', [row[7], -b, -a, -counter], p['retain_strong_per_band']),
                            ('control', [-control, -counter], p['retain_control_per_band'])]:
        heap = state['heaps'].setdefault(band+':'+role, [])
        value = [key, row]
        if len(heap) < size:
            heapq.heappush(heap, value)
        elif key > heap[0][0]:
            heapq.heapreplace(heap, value)


def initial_state():
    return dict(cursor=0, accepted=0, blocks=0, chain='0'*64, heaps={},
                height_counts={}, rejections={'shell_or_zero': 0, 'duplicate': 0}, singular=0)


@contextmanager
def database(folder):
    folder.mkdir(parents=True, exist_ok=True)
    db = sqlite3.connect(folder/'intake.sqlite', timeout=0)
    try:
        db.execute('PRAGMA journal_mode=WAL'); db.execute('PRAGMA synchronous=FULL')
        db.execute('CREATE TABLE IF NOT EXISTS meta (key TEXT PRIMARY KEY, value TEXT NOT NULL)')
        db.execute('CREATE TABLE IF NOT EXISTS seen (a INTEGER, b INTEGER, PRIMARY KEY(a,b)) WITHOUT ROWID')
        db.execute('CREATE TABLE IF NOT EXISTS blocks (id INTEGER PRIMARY KEY, sha256 TEXT, payload BLOB NOT NULL)')
        db.commit()
        yield db
    finally:
        db.close()


def scan(folder, p, arithmetic, max_blocks=None, failpoint=None):
    with database(folder) as db:
        binding = digest(packed(p))
        old = db.execute("SELECT value FROM meta WHERE key='protocol'").fetchone()
        require(old is None or old[0] == binding, 'checkpoint protocol mismatch')
        if old is None:
            with db:
                db.execute('INSERT INTO meta VALUES (?,?)', ('protocol', binding))
                db.execute('INSERT INTO meta VALUES (?,?)', ('state', packed(initial_state()).decode()))
        blocks_done = 0
        while True:
            db.execute('BEGIN IMMEDIATE')
            state = json.loads(db.execute("SELECT value FROM meta WHERE key='state'").fetchone()[0])
            if state['accepted'] >= p['population'] or state['cursor'] >= p['maximum_draws'] or (
                    max_blocks is not None and blocks_done >= max_blocks):
                db.rollback(); break
            rows = []; start = state['cursor']
            while len(rows) < p['block_size'] and state['accepted'] < p['population'] and state['cursor'] < p['maximum_draws']:
                i = state['cursor']; state['cursor'] += 1
                pair = draw(p, i)
                if pair is None:
                    state['rejections']['shell_or_zero'] += 1; continue
                a, b = pair
                if not db.execute('INSERT OR IGNORE INTO seen VALUES (?,?)', (a,b)).rowcount:
                    state['rejections']['duplicate'] += 1; continue
                row = arithmetic.record(i, a, b); rows.append(row); state['accepted'] += 1
                if row[3] == 'SMOOTH':
                    key = str(row[6]); state['height_counts'][key] = state['height_counts'].get(key, 0)+1
                else:
                    state['singular'] += 1
                retain(state, row, p)
            payload = zlib.compress(packed(dict(start=start, end=state['cursor'], rows=rows)), 6)
            h = digest(payload)
            db.execute('INSERT INTO blocks VALUES (?,?,?)', (state['blocks'], h, payload))
            state['blocks'] += 1; state['chain'] = digest((state['chain']+h).encode())
            db.execute("UPDATE meta SET value=? WHERE key='state'", (packed(state).decode(),))
            if failpoint:
                failpoint('before_commit')
            db.commit(); blocks_done += 1
            if failpoint:
                failpoint('after_commit')
            atomic(folder/'progress.json', {k:v for k,v in state.items() if k != 'heaps'})
            require(sum(f.stat().st_size for f in folder.glob('intake.sqlite*')) <= p['resources']['storage_bytes'], 'intake storage cap reached')
            print('INTAKE', state['accepted'], '/', p['population'], 'draws', state['cursor'], flush=True)
        return state


def selection(state, p, arithmetic):
    require(state['accepted'] == p['population'], 'intake incomplete; no selection')
    candidates, seeds, used, seed_used = [], [], set(), set()
    for band in p['height_bands']:
        for role in ('strong', 'control'):
            heap = state['heaps'].get(str(band)+':'+role, [])
            rows = [r for _,r in sorted(heap, reverse=True)]
            count = 0
            for row in rows:
                if row[0] not in used and len(candidates) < p['arithmetic_cap']:
                    candidates.append(arithmetic.candidate(row, role)); used.add(row[0])
                if row[0] in used and row[0] not in seed_used and count < p['seed_'+role+'_per_band'] and len(seeds) < p['seed_cap']:
                    seeds.append(arithmetic.candidate(row, role)); seed_used.add(row[0]); count += 1
    return dict(status='SEALED_EQUATION_ONLY_SELECTION', protocol_digest=digest(packed(p)),
                block_chain=state['chain'], parameters=state['accepted'], draws=state['cursor'],
                height_counts=state['height_counts'], arithmetic_candidates=candidates, seed_inputs=seeds,
                scope='Local trace scheduling and conic splitting only; strict/Selmer/xi UNKNOWN. No points or ranks used.')


def audit(folder, p, arithmetic):
    """Full draw/rejection/block/retention replay, without modifying production DB."""
    import tempfile
    with database(folder) as db, tempfile.TemporaryDirectory(prefix='funnel-audit-') as tmp:
        saved = json.loads(db.execute("SELECT value FROM meta WHERE key='state'").fetchone()[0])
        # Bounded disk-backed replay, including the complete duplicate index.
        rebuilt = scan(Path(tmp), p, arithmetic, max_blocks=saved['blocks'])
        require(rebuilt == saved, 'intake state differs under replay')
        with database(Path(tmp)) as other:
            left = db.execute('SELECT id,sha256,payload FROM blocks ORDER BY id')
            right = other.execute('SELECT id,sha256,payload FROM blocks ORDER BY id')
            for l, r in zip(left, right, strict=True):
                require(l == r and digest(l[2]) == l[1], 'block replay differs')
            require(db.execute('SELECT count(*) FROM seen').fetchone()[0] == saved['accepted'], 'dedup index count differs')
            for left_pair,right_pair in zip(db.execute('SELECT a,b FROM seen ORDER BY a,b'),
                                            other.execute('SELECT a,b FROM seen ORDER BY a,b'),strict=True):
                require(left_pair == right_pair, 'dedup index content differs')
        if (folder/'selection.json').exists():
            require(read(folder/'selection.json') == selection(rebuilt, p, arithmetic), 'selection replay differs')
    result = dict(status='PASS_FULL_INTAKE_REPLAY', parameters=saved['accepted'],
                  block_chain=saved['chain'], protocol_digest=digest(packed(p)))
    atomic(folder/'intake-replay.json', result, immutable=True)
    return result


def queue(rank):
    require(type(rank) is int and rank >= 18, 'queue needs a certified lower bound >=18')
    for threshold, name in [(32,'TARGET_LOWER_BOUND'), (30,'AGGRESSIVE'), (28,'DEEP'), (20,'INTERESTING'), (18,'RETAIN')]:
        if rank >= threshold:
            return name


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['freeze','scan','audit'])
    parser.add_argument('--directory', type=Path, required=True)
    parser.add_argument('--profile', choices=['smoke','production'], default='smoke')
    args = parser.parse_args(); folder = args.directory.resolve()
    if args.action == 'freeze':
        freeze(folder, args.profile); return
    p = protocol(folder); arithmetic = Arithmetic(p)
    tables = read(folder/'tables-replay.json')
    require(tables['status'] == 'PASS_INDEPENDENT_PARI_TABLES' and
            tables['table_digest'] == digest(packed(arithmetic.tables)), 'independent table gate missing or changed')
    if args.action == 'audit':
        audit(folder, p, arithmetic); return
    state = scan(folder, p, arithmetic)
    require(state['accepted'] == p['population'], 'finite draw cap reached; preserve partial intake')
    atomic(folder/'selection.json', selection(state, p, arithmetic), immutable=True)


if __name__ == '__main__':
    main()
