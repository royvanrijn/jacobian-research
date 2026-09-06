#!/usr/bin/env python3
"""Real-cache regressions and one untouched higher-slice corrected-score gate."""
import argparse
import math
from hashlib import sha256
from pathlib import Path
import certify_compact_r17_candidates as cert
import scan_broad_mw16_higher_annuli as broad
import encode_mw16_joint_caches as base
import encode_mw16_local_corrections as local
import benchmark_11952_annulus_cache_v3 as reader
import build_mw16_local_score_corrections as models
from research_runtime.store import checkpoint
from research_runtime.supervisor import capture, Limits

ROOT = broad.ROOT
CAS = broad.CAS
D = ROOT/'artifacts/local/elliptic-curves/corrected-mw16-annulus-benchmark-v1'
OUT = broad.ART/'corrected_mw16_annulus_benchmark_v1.json'
CPP = CAS/'newfamily/scan_corrected_mw16_annulus.cpp'
TEST = ROOT/'elliptic-curves/tests/test_corrected_mw16_annulus.py'
BINARY = D/'scanner'


def sources():
    paths = [Path(__file__).resolve(), CPP, TEST, CAS/'newfamily/scan_joint_cache_annulus.cpp',
             base.OUT, local.OUT, models.OUT, Path(models.__file__), reader.BINARY, reader.OUT,
             broad.D/'protocol.json', broad.D/'replay.json', broad.narrow.D/'protocol.json',
             CAS/'research_runtime/supervisor.py', CAS/'research_runtime/store.py']
    return {str(p.relative_to(ROOT)): cert.hashed(p) for p in paths}


def prepare():
    if (D/'protocol.json').exists():
        raise FileExistsError('preserve one corrected-score benchmark')
    if (cert.read(broad.D/'controller/ledger.json')['status'] != 'PASS'
            or cert.read(local.OUT)['status'] != 'PASS' or cert.read(base.OUT)['status'] != 'PASS'):
        raise ArithmeticError('complete base/local cache and broad scan gates required')
    rows = cert.read(broad.narrow.D/'protocol.json')['rows']
    selected = next(r for r in rows if (r['family'], r['band'], r['sign']) == ('a1-fibration-01', 2, 1))
    excluded = {selected['shard'], selected['excluded_previous_shard']}
    excluded.update(r['shard'] for r in cert.read(broad.D/'protocol.json')['rows']
                    if (r['family'], r['band'], r['sign']) == ('a1-fibration-01', 2, 1))
    choices = [(sha256(f'corrected-mw16-first-new-slice-v1|{q}'.encode()).hexdigest(), q)
               for q in range(selected['shard'] % 2, selected['shards'], 2) if q not in excluded]
    choice, shard = min(choices)
    real = {**selected, 'id': 'corrected-first-new-slice', 'shard': shard,
            'excluded_previous_shards': sorted(excluded), 'choice_sha256': choice}
    real['primitive_population'] = (broad.old.counts.population(real['N'], real['M'], shard, real['shards'])
                                   - broad.old.counts.population(real['inner'], real['inner'], shard, real['shards']))
    cases = [{**r, 'N': 351, 'M': 2*r['shards'], 'inner': r['shards']} for r in rows]
    bindings = {}
    for family in sorted({r['family'] for r in rows}):
        for label in ('short', 'extended'):
            path = base.cache(family, label)
            bindings[str(path.relative_to(ROOT))] = cert.hashed(path)
        path = local.D/(family+'.bin')
        bindings[str(path.relative_to(ROOT))] = cert.hashed(path)
    checkpoint(D/'protocol.json', {'sources': sources(), 'cache_bindings': bindings,
        'cases': cases, 'real_case': real, 'keep': 4096, 'seconds_per_call': 120,
        'cost_gate_seconds': 45, 'outer_seconds': 1800, 'rss_bytes': 1610612736,
        'scope': 'Three signed/annular/ordering test groups and twenty complete real-cache '
                 'signed frames across all five families and both higher-band denominator moduli '
                 'precede exactly one new positive family01 slice in16384<H<=65536. The slice '
                 'excludes the previous outer residue, narrow-trial residue and all sixteen broad '
                 'residues. Score all3510 primes plus the exact restored5/13 terms before retaining4096. '
                 'Independent retained-cache components and exact full-model local scaling check all '
                 'returned scores. One worker,120 seconds per call,45-second full-slice cost gate, '
                 '1800 seconds total. No retry, automatic broader scan, point search, catalogue, '
                 'known-record input or change to either existing experiment.'})


def protocol():
    p = cert.read(D/'protocol.json')
    if p['sources'] != sources() or any(cert.hashed(ROOT/n) != h for n, h in p['cache_bindings'].items()):
        raise ArithmeticError('frozen corrected benchmark input changed')
    return p


def execute(name, command, check=False):
    path = D/(name+'.json')
    if not check:
        if path.exists():
            raise FileExistsError('preserve corrected benchmark invocation')
        c = capture(command, limits=Limits(120, 1610612736), log_path=D/(name+'.log'), separate_stderr=True, check=False)
        checkpoint(path, {'command': command, 'stdout': c.stdout, 'stderr': c.stderr, 'supervision': c.supervision})
    raw = cert.read(path)
    if (raw['command'] != command or raw['supervision']['command'] != command
            or raw['supervision']['outcome'] != 'completed' or raw['supervision']['returncode'] != 0):
        raise ArithmeticError('corrected benchmark failed/censored')
    return raw


def references(rows, family, name, check):
    path = D/(name+'-candidates.txt')
    text = 'R17-CANDIDATES-V1 '+str(len(rows))+'\n'+''.join(f'{n} {d}\n' for n, d, *_ in rows)
    if check:
        if path.read_text() != text:
            raise ArithmeticError('reference roster differs')
    else:
        if path.exists():
            raise FileExistsError('preserve reference roster')
        path.write_text(text)
    sums = []
    for label, count in [('short', 562), ('extended', 2948)]:
        raw = execute(name+'-'+label, [str(reader.BINARY), str(base.cache(family, label)), str(path), '262144'], check)
        lines = raw['stdout'].splitlines()
        if raw['stderr'] or len(lines) != len(rows)+1 or lines[-1] != f'S {len(rows)} {count}':
            raise ArithmeticError('reference cache frame differs')
        values = []
        for i, line in enumerate(lines[:-1]):
            v = line.split()
            if len(v) != 4 or v[:2] != ['R', str(i)]:
                raise ArithmeticError('reference index differs')
            values.append(tuple(map(int, v[2:])))
        sums.append(values)
    f = next(f for f in cert.read(models.spec.ATLAS)['families'] if f['fibration_id'] == family)
    output = []
    trace_cache = {}
    for row, first, second in zip(rows, *sums):
        n, d = row[:2]
        A, B = models.coefficients(f, n, d)
        units, good = first[0]+second[0], first[1]+second[1]
        for q in (5, 13):
            a, b, exponent = A, B, 0
            while a % q**4 == 0 and b % q**6 == 0:
                a //= q**4; b //= q**6; exponent += 1
            if exponent:
                key = (a % q, b % q, q)
                if key not in trace_cache:
                    trace_cache[key] = models.diagnostic.trace(*key)
                ap = trace_cache[key]
                if ap is not None:
                    units += round((2-ap)/(q+1-ap)*math.log(q)*10**12); good += 1
        output.append([n, d, units, good])
    return output


def scan(c, keep, name, check):
    args = [c['sign'], c['N'], c['M'], keep, c['shard'], c['shards'], c['inner'], 3510]
    raw = execute(name, [str(BINARY), str(base.cache(c['family'], 'short')), str(base.cache(c['family'], 'extended')),
                         *map(str, args), str(local.D/(c['family']+'.bin'))], check)
    count = (broad.old.counts.population(c['N'], c['M'], c['shard'], c['shards'])
             - broad.old.counts.population(min(c['N'], c['inner']), min(c['M'], c['inner']), c['shard'], c['shards']))
    lines = raw['stdout'].splitlines()
    if (raw['stderr'] or len(lines) != min(count, keep)+3
            or lines[:2] != ['CORRECTED_MW16_ANNULUS_V1', 'P '+' '.join(map(str, args))]
            or lines[-1] != f'S {count} {min(count, keep)}'):
        raise ArithmeticError('corrected scanner frame differs')
    rows = []
    for line in lines[2:-1]:
        v = line.split()
        if len(v) != 5 or v[0] != 'C':
            raise ArithmeticError('corrected row framing differs')
        n, d, units, good = map(int, v[1:])
        if (n*c['sign'] <= 0 or math.gcd(n, d) != 1 or not 1 <= abs(n) <= c['N']
                or not 1 <= d <= c['M'] or max(abs(n), d) <= c['inner']
                or (d-1) % c['shards'] != c['shard'] or not 0 <= good <= 3510):
            raise ArithmeticError('corrected candidate outside primitive frame')
        rows.append([n, d, units, good])
    if len({tuple(r[:2]) for r in rows}) != len(rows) or rows != sorted(rows, key=lambda r: (-r[2], -r[3], r[1], abs(r[0]))):
        raise ArithmeticError('corrected score ordering or uniqueness differs')
    return rows, raw, count


def run(check=False):
    p = protocol()
    execute('compile', ['g++', '-O3', '-std=c++17', str(CPP), '-o', str(BINARY)], check)
    test = execute('tests', ['/usr/bin/python3', str(TEST)], check)
    if 'Ran 3 tests' not in test['stderr'] or not test['stderr'].rstrip().endswith('OK'):
        raise ArithmeticError('all three corrected scanner regression groups required')
    frames = []
    for i, c in enumerate(p['cases']):
        name = f'frame-{i:02}'
        bare = [[c['sign']*n, d] for d in range(c['shard']+1, c['M']+1, c['shards'])
                for n in range(1, c['N']+1) if math.gcd(n, d) == 1 and max(n, d) > c['inner']]
        expected = references(bare, c['family'], name+'-reference', check)
        expected.sort(key=lambda r: (-r[2], -r[3], r[1], abs(r[0])))
        for keep in (1000, 7):
            actual, raw, count = scan(c, keep, name+'-keep'+str(keep), check)
            if actual != expected[:keep] or count != len(expected):
                raise ArithmeticError('complete real-cache corrected frame/top7 differs')
        frames.append({'id': c['id'], 'primitive_addresses': len(expected), 'status': 'PASS'})
    rows, raw, count = scan(p['real_case'], p['keep'], 'new-full-slice', check)
    if count != p['real_case']['primitive_population'] or references(rows, p['real_case']['family'], 'new-full-reference', check) != rows:
        raise ArithmeticError('new full-slice count or local-corrected scores differ')
    result = {'schema': 'elliptic-curves.corrected-mw16-annulus-benchmark.v1', 'status': 'PASS',
              'sources': sources(), 'protocol_sha256': cert.hashed(D/'protocol.json'),
              'binary_sha256': cert.hashed(BINARY), 'frames': frames, 'primitive_addresses': count,
              'retained_rows': rows, 'wall_seconds': raw['supervision']['wall_seconds'],
              'cost_gate_passed': raw['supervision']['wall_seconds'] <= p['cost_gate_seconds'],
              'claim_boundary': p['scope']}
    if check:
        if cert.read(OUT) != result:
            raise ArithmeticError('corrected benchmark replay differs')
    else:
        if OUT.exists():
            raise FileExistsError('preserve corrected benchmark proof')
        checkpoint(OUT, result)
    print('CORRECTED MW16 FULL-SLICE BENCHMARK', count, result['wall_seconds'], result['cost_gate_passed'], flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('stage', choices=['prepare', 'run', 'check'])
    args = parser.parse_args()
    prepare() if args.stage == 'prepare' else run(args.stage == 'check')
