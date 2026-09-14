#!/usr/bin/env sage -python
"""Freeze generic equations, then test a fixed small Mestre branch gate.

No section enumeration, auxiliary-function search, or specialization input.
Missing finite witnesses leave the corresponding assertion UNKNOWN.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import resource
import time

from sage.all import EllipticCurve, GF, PolynomialRing, QQ, prime_range, version

ROOT = Path(__file__).resolve().parents[2]
SOURCES = {
    'published-r17': 'artifacts/generated-results/elkies-k3-r17-mestre-shared-twist-v1/input.json',
    'alternate-q80': 'artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json',
    'curve302-parent': 'artifacts/generated-results/elliptic-curves/curve302_recovered_mw17_parent_v1.json',
    'x1092-class1': 'artifacts/generated-results/elliptic-curves/x1092_class1_realization_compact_parent_v1.json',
}
R = PolynomialRing(QQ, 't')


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def write_new(path, obj):
    with path.open('x') as stream:
        json.dump(obj, stream, indent=2, sort_keys=True)
        stream.write('\n')


def freeze(directory):
    parents = []
    for name, source in SOURCES.items():
        data = json.loads((ROOT / source).read_text())
        if name == 'published-r17':
            A, B = R(data['A']), R(data['B'])
        elif name == 'alternate-q80':
            model = data['weierstrass_model']
            A = R(model['A_coefficients_low_to_high'])
            B = R(model['B_coefficients_low_to_high'])
        else:
            ainvs = [R(a['numerator']) / R(a['denominator']) for a in data['a_invariants']]
            E = EllipticCurve(R.fraction_field(), ainvs)
            A, B = R(-E.c4()/48), R(-E.c6()/864)
        parents.append({'name': name, 'source': source,
                        'source_sha256': digest(ROOT/source),
                        'A': [str(c) for c in A.list()],
                        'B': [str(c) for c in B.list()]})
    directory.mkdir(parents=True, exist_ok=True)
    write_new(directory/'input.json', {
        'schema': 'mestre-parent-branch-input-v1', 'parents': parents,
        'selection': 'Four retained generic equations only; u=2 control. No exceptional points or specialization data are projected.',
        'prime_min': 5, 'prime_max_exclusive': 2000,
        'cpu_seconds': 40, 'address_space_bytes': 4294967296,
        'u': '2', 'sage_version': version(),
        'producer_sha256': digest(Path(__file__)),
    })


def run(directory):
    packet = json.loads((directory/'input.json').read_text())
    assert packet['producer_sha256'] == digest(Path(__file__))
    resource.setrlimit(resource.RLIMIT_CPU, (packet['cpu_seconds'], packet['cpu_seconds']))
    resource.setrlimit(resource.RLIMIT_AS, (packet['address_space_bytes'], packet['address_space_bytes']))
    start = time.process_time()
    rows = []
    # Persist each completed parent; an interrupted run has no final certificate.
    for parent in packet['parents']:
        A, B = R(parent['A']), R(parent['B'])
        core = 9261*B**2 + 400*A**3
        D = -5*A*B*core
        delta = 4*A**3 + 27*B**2
        row = {'name': parent['name'], 'fields': {'A': {}, 'B': {}},
               'degrees': [int(f.degree()) for f in (A, B, core, D, delta)],
               'smooth_disjoint_prime': None}
        for p in prime_range(packet['prime_min'], packet['prime_max_exclusive']):
            p = int(p)
            F = PolynomialRing(GF(p), 't')
            try:
                a, b, d, e = [F(f) for f in (A, B, D, delta)]
            except (ValueError, ZeroDivisionError):
                continue
            for name, f, original in [('A', a, A), ('B', b, B)]:
                entry = row['fields'][name]
                if f.degree() != original.degree():
                    continue
                if 'irreducible_prime' not in entry and f.is_irreducible():
                    entry['irreducible_prime'] = p
                congruence = p % 12 == 11 if name == 'A' else p % 4 == 3
                if 'local_place' not in entry and congruence:
                    roots = sorted(int(x) for x in f.roots(multiplicities=False) if f.derivative()(x))
                    if roots:
                        entry['local_place'] = {'prime': p, 'root': roots[0]}
            if (row['smooth_disjoint_prime'] is None and d.degree() == 44
                    and e.degree() == 24 and d.is_squarefree()
                    and e.is_squarefree() and d.gcd(e).degree() == 0):
                row['smooth_disjoint_prime'] = p
            if row['smooth_disjoint_prime'] and all(len(v) == 2 for v in row['fields'].values()):
                break
        row['status'] = ('PASS' if row['smooth_disjoint_prime'] and
                         all(len(v) == 2 for v in row['fields'].values()) else 'UNKNOWN')
        rows.append(row)
        write_new(directory/(parent['name']+'-checkpoint.json'), row)
    write_new(directory/'result.json', {
        'schema': 'mestre-parent-branch-result-v1',
        'input_sha256': digest(directory/'input.json'), 'parents': rows,
        'status': 'PASS' if all(r['status'] == 'PASS' for r in rows) else 'UNKNOWN',
        'component_cpu_seconds': round(time.process_time()-start, 6),
        'positive_low_genus_cover': False,
    })
    print(json.dumps(rows, sort_keys=True))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['freeze', 'run'])
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    (freeze if args.action == 'freeze' else run)(args.output)
