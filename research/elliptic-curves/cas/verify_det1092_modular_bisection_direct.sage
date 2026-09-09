#!/usr/bin/env sage-python
"""Independent replay of the fixed 15-orbit, three-prime incidence panel.

Manual elliptic addition in reversed word order; no EllipticCurve, nullspace
solver, producer imports, rational point additions, or additional exposures.
Written maximal-pole lifting lemma is a separate mathematical dependency.
"""
import csv
import hashlib
import json
import signal
from collections import Counter
from pathlib import Path
from sage.all import QQ, ZZ, GF, PolynomialRing, matrix, vector

signal.alarm(25)
ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / 'artifacts/generated-results/elliptic-curves'
OUT = ART / 'det1092_modular_bisection_direct_v1'
PRE = ART / 'det1092_modular_bisection_exclusion_v1'
REG = ART / 'det1092_rational_bisection_index_v1'
PARENT = ART / 'curve302_recovered_mw17_parent_v1.json'
TABLE = ART / 'curve302_parent_degree2_multisection_orbits_v1.tsv'
CENSUS = ART / 'curve302_parent_degree2_multisection_lattice_v1.json'
PRIMES = [149, 151, 157]
REGRESSIONS = [8044, 47755, 103186]


def read(path):
    return json.loads(path.read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def binding(record):
    for name, digest in record['inputs'].items():
        assert sha(ROOT / name) == digest, name


protocol, summary = read(OUT / 'protocol.json'), read(OUT / 'summary.json')
binding(protocol)
binding(summary)
binding(read(PRE / 'protocol.json'))
binding(read(PRE / 'summary.json'))
assert protocol['prime_pool'] == PRIMES
assert protocol['limits'] == {
    'complete_atlas_runs': 0, 'new_orbits': 12, 'point_searches': 0,
    'prime_trials': 45, 'rational_RR_solves': 0,
    'rational_point_additions': 0, 'regression_orbits': 3, 'seconds': 25,
}
assert read(CENSUS)['orbits_tsv_sha256'] == sha(TABLE)
with TABLE.open() as stream:
    rows = [row for row in csv.DictReader(stream, delimiter='\t')
            if row['category'] == 'rational']
assert len(rows) == 40917
words = {int(row['orbit_mask']): list(map(int, row['parent_MW17_w'].split()))
         for row in rows}
assert len(words) == 40917
selected = sorted(words)[:12]
assert selected == [61, 87, 103, 107, 109, 110, 111, 115, 117, 121, 122, 123]
masks = sorted(set(selected + REGRESSIONS))
assert len(masks) == 15
parent = read(PARENT)
G = matrix(ZZ, parent['generic_height_gram'])
records = {m: read(OUT / ('orbit-%d.json' % m)) for m in masks}
prior = {m: read(PRE / ('orbit-%d.json' % m)) for m in masks}
traces = {m: read(PRE / ('trace-%d.json' % m)) for m in masks}
for m in masks:
    w = vector(ZZ, words[m])
    assert w * G * w == 10
    assert records[m]['mask'] == m and records[m]['word'] == words[m]
    assert [r['prime'] for r in records[m]['trials']] == PRIMES
    assert prior[m]['trace_sha256'] == sha(PRE / ('trace-%d.json' % m))

# Construct short coefficients by invariant formulas, without EllipticCurve.
RQ = PolynomialRing(QQ, 't')
FQ = RQ.fraction_field()


def decq(record):
    return FQ(RQ(record['numerator'])) / RQ(record['denominator'])


a1, a2, a3, a4, a6 = [decq(r) for r in parent['a_invariants']]
b2 = a1*a1 + 4*a2
b4 = a1*a3 + 2*a4
b6 = a3*a3 + 4*a6
AQ = RQ(-(b2*b2 - 24*b4)/48)
BQ = RQ(-(-b2**3 + 36*b2*b4 - 216*b6)/864)
assert AQ.degree() <= 8 and BQ.degree() <= 12
assert -16*(4*AQ(0)**3 + 27*BQ(0)**2)
assert len(parent['basis_weierstrass_coordinates']) == 17
for pair in parent['basis_weierstrass_coordinates']:
    x, y = map(decq, pair)
    assert y*y + a1*x*y + a3*y == x**3 + a2*x*x + a4*x + a6

# Saved characteristic-zero traces are only a regression, not the source
# of the modular word calculation or the lifting argument.
for m in masks:
    tr = traces[m]
    assert tr['word'] == words[m] and tr['mask'] == m
    x, y = map(decq, tr['short_trace'])
    h, nx, ny = [RQ(tr[key]) for key in ['pole', 'nx', 'ny']]
    assert h.is_monic() and h.degree() <= 3
    assert nx.degree() <= 10 and ny.degree() <= 15
    assert x == nx/h**2 and y == ny/h**3
    assert y*y == x**3 + AQ*x + BQ

counts = Counter()
certificates = []
comparisons = 0
old_squareclass_checks = 0

for ip, p in enumerate(PRIMES):
    fp = GF(p)
    R = PolynomialRing(fp, 't')
    K = R.fraction_field()
    t = R.gen()
    A, B = R(AQ), R(BQ)
    assert -16*(4*A**3 + 27*B**2)  # Good reduction for the Gauss valuation.

    def dec(record):
        numerator = R([fp(QQ(c)) for c in record['numerator']])
        denominator = R([fp(QQ(c)) for c in record['denominator']])
        assert denominator
        return K(numerator)/denominator

    ap1, ap3, bp2 = K(a1), K(a3), K(b2)
    basis = []
    for pair in parent['basis_weierstrass_coordinates']:
        xx, yy = map(dec, pair)
        X, Y = xx + bp2/12, yy + (ap1*xx + ap3)/2
        assert Y*Y == X**3 + A*X + B
        basis.append((X, Y))

    def add(P, Q):
        if P is None:
            return Q
        if Q is None:
            return P
        x, y = P
        u, v = Q
        if x == u:
            if y == -v:
                return None
            assert y == v and y
            slope = (3*x*x + A)/(2*y)
        else:
            slope = (v-y)/(u-x)
        z = slope*slope - x - u
        return (z, -y + slope*(x-z))

    def mul(n, P):
        if n < 0:
            return mul(-n, (P[0], -P[1]))
        result = None
        while n:
            if n & 1:
                result = add(result, P)
            P = add(P, P)
            n //= 2
        return result

    for m in masks:
        row = records[m]['trials'][ip]
        trace = None
        for n, P in reversed(list(zip(words[m], basis))):
            trace = add(trace, mul(-n, P))
        status = row['status']
        counts[status] += 1
        if trace is None or trace[0].denominator().degree() != 6:
            assert row == {'prime': p, 'status': 'UNKNOWN_MAXIMAL_POLE_GATE'}
            continue
        X, Y = trace
        assert Y*Y == X**3 + A*X + B
        h, nx, ny = [R(row[key]) for key in ['pole', 'nx', 'ny']]
        assert h.is_monic() and h.degree() == 3
        assert h*h == X.denominator().monic()
        assert nx.degree() <= 10 and ny.degree() <= 15
        assert X == nx/h**2 and Y == ny/h**3
        assert all(R(RQ(traces[m][key])) == R(row[key]) for key in ['pole', 'nx', 'ny'])
        if not h(0) or not -16*(4*A(0)**3 + 27*B(0)**2):
            assert status == 'UNKNOWN_BAD_TARGET_OR_TRACE'
            continue
        columns = ([t**i*h**3 for i in range(10)]
                   + [t**i*nx*h for i in range(6)]
                   + [t**i*ny for i in range(4)])
        assert all(c.degree() <= 18 for c in columns)
        mat = matrix(fp, [[c[i] for c in columns] for i in range(19)])
        if status == 'UNKNOWN_RR_RANK_DROP':
            assert mat.rank() < 19
            continue
        kernel = vector(fp, row['kernel'])
        assert len(kernel) == 20 and kernel and mat*kernel == 0
        pivots = row['pivot_columns']
        assert len(pivots) == len(set(pivots)) == 19
        assert all(0 <= j < 20 for j in pivots)
        det = mat.matrix_from_columns(pivots).det()
        assert det and int(det) == row['pivot_minor_determinant']
        g0, g1, g2 = [R(list(v)) for v in [kernel[:10], kernel[10:16], kernel[16:]]]
        assert g0*h**3 + g1*nx*h + g2*ny == 0
        if not g2(0):
            assert status == 'UNKNOWN_AFFINE_CHART'
            continue
        # Polynomial division, not the producer's synthetic-division formulas.
        RX = PolynomialRing(fp, 'X')
        u = RX.gen()
        eliminated = (g0(0) + g1(0)*u)**2 - g2(0)**2*(u**3 + A(0)*u + B(0))
        residual, remainder = eliminated.quo_rem(u - nx(0)/h(0)**2)
        assert remainder == 0 and residual.degree() == 2
        assert list(map(int, residual.list())) == row['residual_coefficients']
        disc = residual.discriminant()
        assert int(disc) == row['discriminant']
        roots = [int(v) for v in fp if residual(v) == 0]
        expected = 'EXCLUDED_NONSQUARE' if not roots else 'UNKNOWN_SPLIT_OR_REPEATED_REDUCTION'
        assert status == expected
        if status == 'EXCLUDED_NONSQUARE':
            assert pow(int(disc), (p-1)//2, p) == p-1
            certificates.append({'mask': m, 'prime': p, 'discriminant': int(disc),
                                 'residual_coefficients': row['residual_coefficients']})
        old = prior[m]['trials'][ip]
        for key in ['status', 'kernel', 'pivot_columns', 'pivot_minor_determinant',
                    'residual_coefficients', 'discriminant']:
            assert old[key] == row[key]
        comparisons += 1
        if m in REGRESSIONS and disc:
            q0 = QQ(read(REG / ('orbit-%d.json' % m))['q'][0])
            valuation = q0.valuation(p)
            assert valuation % 2 == 0
            unit = int(fp(q0/QQ(p)**valuation))
            assert pow(unit, (p-1)//2, p) == pow(int(disc), (p-1)//2, p)
            old_squareclass_checks += 1

assert sum(counts.values()) == 45
assert dict(counts) == summary['status_counts']
excluded = sorted({row['mask'] for row in certificates})
for m in masks:
    assert records[m]['excluded'] == (m in excluded)
assert summary['new_masks'] == selected and summary['regression_masks'] == REGRESSIONS
assert summary['primes'] == PRIMES
assert summary['new_excluded'] == [m for m in selected if m in excluded]
assert summary['new_unknown'] == [m for m in selected if m not in excluded]
assert summary['regression_excluded'] == [m for m in REGRESSIONS if m in excluded]
report = {
    'status': 'PASS_INDEPENDENT_FINITE_FIELD_BISECTION_EXCLUSIONS',
    'classification': 'verified application of the written maximal-pole lifting and unit-minor lemma',
    'new_excluded': summary['new_excluded'], 'new_unknown': summary['new_unknown'],
    'regression_excluded': summary['regression_excluded'], 'prime_trials': 45,
    'status_counts': dict(sorted(counts.items())),
    'certificates': sorted(certificates, key=lambda row: (row['mask'], row['prime'])),
    'rational_trace_modular_RR_regressions': comparisons,
    'old_exact_cover_squareclass_checks': old_squareclass_checks,
    'independent_methods': ['manual reverse-order elliptic addition', 'saved nonzero minor and kernel validation',
                            'polynomial division', 'enumeration of every residue for residual root exclusion'],
    'boundary': 'Twelve new orbits and three regressions only. No rational incidence or rank upper bound follows from survival. No full-atlas run, extra prime, exceptional point, or point search.',
    'inputs': {str(path.relative_to(ROOT)): sha(path) for path in [
        OUT/'protocol.json', OUT/'summary.json', PRE/'protocol.json', PRE/'summary.json',
        *[REG/('orbit-%d.json' % m) for m in REGRESSIONS], Path(__file__)]},
}
destination = OUT / 'independent-replay.json'
if destination.exists():
    assert read(destination) == report
else:
    with destination.open('x') as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write('\n')
print(report['status'], 'new excluded', report['new_excluded'],
      'unknown', report['new_unknown'], 'all 45 exposures retained', flush=True)
