#!/usr/bin/env python3
"""Bounded, point-free local valuation-defect audit of fixed quartet carriers."""
import argparse
from pathlib import Path
import subprocess
import retrospective as r

HERE = Path(__file__).resolve().parent
INPUT = r.OUT / 'rank_jump_soluble_quartet_compression_inputs_v1.json'
SUPPORT = r.OUT / 'rank_jump_collision_prime_lift_v1.json'
PROTOCOL = HERE / 'COLLISION_DEFECT_PROTOCOL.json'
OUTPUT = r.OUT / 'rank_jump_collision_defect_v1.json'
WORK = r.ROOT / 'artifacts/local/rank-jump-collision-defect-v1'


def evaluate(q, t):
    return q[0] + q[1]*t + q[2]*t*t


def worker(index):
    from sage.all import GF, PolynomialRing, ZZ, prod, version
    case = r.read(INPUT)['cases'][index]
    source = r.read(SUPPORT)['rows'][index]['result']
    forms = case['covers']
    qs = [list(map(int, c['form'])) for c in forms]
    rows = []
    for pstr in source['collision_primes']:
        p = int(pstr)
        ring = PolynomialRing(GF(p), 'x')
        pols = [ring(q) for q in qs]
        pair_rows = []
        collision_roots = set()
        for i in range(4):
            for j in range(i):
                g = pols[j].gcd(pols[i]).monic()
                roots = sorted(int(a) for a in g.roots(multiplicities=False))
                infinity = qs[j][2] % p == qs[i][2] % p == 0
                collision_roots.update(roots)
                pair_rows.append({'indices': [j, i], 'gcd': [int(a) for a in g.list()],
                                  'finite_roots': list(map(str, roots)), 'infinity': infinity})
        clusters = []
        for root in sorted(collision_roots):
            mask = sum(1 << i for i, q in enumerate(qs) if evaluate(q, root) % p == 0)
            clusters.append({'chart': 'finite', 'root': str(root), 'mask': mask})
        infinity_mask = sum(1 << i for i, q in enumerate(qs) if q[2] % p == 0)
        if infinity_mask.bit_count() >= 2:
            clusters.append({'chart': 'infinity', 'root': '0', 'mask': infinity_mask})
        possible_masks = sorted({0} | {m for c in clusters for m in range(1, 16)
                                      if m & c['mask'] == m and m.bit_count() % 2 == 0})
        row = {'prime': pstr, 'pairs': pair_rows, 'collision_clusters': clusters,
               'necessary_parity_masks': possible_masks, 'prunable': not clusters}
        incident = [(pair['indices'], int(e)) for pair in source['pair_resultants']
                    for pp, e in pair['prime_factorization'] if pp == pstr]
        ordinary = p >= 5 and len(incident) == 1 and incident[0][1] == 1
        ordinary = ordinary and len(clusters) == 1 and clusters[0]['chart'] == 'finite'
        if ordinary:
            root = int(clusters[0]['root']); pair = incident[0][0]
            ordinary = all((q[1] + 2*q[2]*root) % p for q in (qs[i] for i in pair))
        row['ordinary_simple_pair_collision'] = bool(ordinary)
        if ordinary:
            i, j = pair
            a = [(evaluate(qs[k], root)//p) % p for k in pair]
            b = [(qs[k][1] + 2*qs[k][2]*root) % p for k in pair]
            assert (a[0]*b[1]-a[1]*b[0]) % p != 0
            other = prod(evaluate(qs[k], root) % p for k in range(4) if k not in pair) % p
            leading = int(other*b[0]*b[1] % p)
            chi = 1 if pow(leading, (p-1)//2, p) == 1 else -1
            row['odd_pair_residue_count'] = str((p-2-chi)//2)
            row['odd_pair_witness'] = {'status': 'UNKNOWN'}
            for s in range(min(p, 64)):
                residue = root+p*s
                values = [evaluate(q, residue) for q in qs]
                units = [(v//p if k in pair else v) % p for k, v in enumerate(values)]
                unit = int(prod(units) % p)
                if unit and pow(unit, (p-1)//2, p) == 1:
                    row['odd_pair_witness'] = {'status': 'CERTIFIED', 'residue_mod_p_squared': str(residue),
                        's': s, 'mask': (1 << i) | (1 << j), 'product_unit': str(unit)}
                    break
            row['zero_parity_witness'] = {'status': 'UNKNOWN'}
            for residue in range(min(p, 64)):
                values = [evaluate(q, residue) % p for q in qs]
                unit = int(prod(values) % p)
                if unit and pow(unit, (p-1)//2, p) == 1:
                    row['zero_parity_witness'] = {'status': 'CERTIFIED', 'residue_mod_p': str(residue),
                        'mask': 0, 'product_unit': str(unit)}
                    break
        rows.append(row)
    masks = [x['odd_pair_witness']['mask'] for x in rows if x.get('odd_pair_witness', {}).get('status') == 'CERTIFIED']
    result = {'id': case['id'], 'rows': rows, 'retained_support': [x['prime'] for x in rows if not x['prunable']],
              'pruned_primes': [x['prime'] for x in rows if x['prunable']],
              'certified_pair_masks': sorted(set(masks)), 'native_mask_span_dimension': r.rank(masks),
              'separate_prime_defect_witnesses': len(masks), 'sage_version': version()}
    r.write_new(WORK / str(index) / 'result.json', result)
    print(case['id'], 'pruned', result['pruned_primes'], 'pair masks', result['certified_pair_masks'],
          'span', r.rank(masks), 'prime witnesses', len(masks), flush=True)


def run():
    rows = []
    for index in range(3):
        wd = WORK / str(index); wd.mkdir(parents=True, exist_ok=True)
        log = wd / 'worker.log'; execution = wd / 'execution.json'
        if not log.exists():
            with log.open('x') as out:
                try:
                    proc = subprocess.run(['/home/royvanrijn/.local/bin/sage', '-python', str(Path(__file__).resolve()),
                                           'worker', '--case', str(index)], stdout=out, stderr=out, timeout=60)
                    state = {'status': 'COMPLETE' if proc.returncode == 0 else 'FAILED', 'returncode': proc.returncode}
                except subprocess.TimeoutExpired:
                    state = {'status': 'TIMEOUT'}
            r.write_new(execution, state)
        result = wd / 'result.json'
        rows.append({'execution': r.read(execution), 'log': log.read_text(),
                     'result': r.read(result) if result.exists() else {'status': 'UNKNOWN'}})
        print(rows[-1]['execution'], rows[-1]['log'], flush=True)
    r.write_new(OUTPUT, {'schema': 'rank-jump.collision-defect.v1', 'layer': 'solubility', 'rows': rows,
        'bindings': {str(p.relative_to(r.ROOT)): r.digest(p.read_bytes())
                     for p in (INPUT, SUPPORT, PROTOCOL, Path(__file__), HERE / 'retrospective.py')},
        'boundary': r.read(PROTOCOL)['failure_semantics']})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('mode', choices=['run', 'worker'])
    parser.add_argument('--case', type=int); args = parser.parse_args()
    if args.mode == 'run': run()
    else: worker(args.case)
