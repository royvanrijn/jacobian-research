#!/usr/bin/env python3
"""Small immutable binary score corrections, with byte-for-byte replay."""
import argparse
import struct
from pathlib import Path
import build_mw16_local_score_corrections as local
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint

ROOT = local.ROOT
D = ROOT/'artifacts/local/elliptic-curves/mw16-local-binary-corrections-v1'
OUT = local.reduction.first.support.ART/'mw16_local_binary_corrections_v1.json'


def encoded(tables):
    data = b'MW16LC01'+struct.pack('<I', 2)
    if [(t['prime'], t['modulus']) for t in tables] != [(5, 125), (13, 169)]:
        raise ArithmeticError('fixed local ring ordering required')
    for table in tables:
        rows = table['entries']
        if [r['index'] for r in rows] != list(range(table['modulus']+table['modulus']//table['prime'])):
            raise ArithmeticError('complete local ring frame required')
        data += struct.pack('<III', table['prime'], table['modulus'], len(rows))
        data += struct.pack('<'+'q'*len(rows), *[r['correction_units'] for r in rows])
        data += bytes(int(r['restored_good']) for r in rows)
    return data+b'ENDLC001'


def run(check=False):
    source = cert.read(local.OUT)
    if source['status'] != 'PASS' or source['projective_ring_entries'] != 1660:
        raise ArithmeticError('complete exact local tables required')
    if any(cert.hashed(ROOT/n) != h for n, h in source['sources'].items()):
        raise ArithmeticError('local proof inputs changed')
    rows = []
    for family in sorted({t['family'] for t in source['tables']}):
        data = encoded([t for t in source['tables'] if t['family'] == family])
        path = D/(family+'.bin')
        if not check:
            if path.exists():
                raise FileExistsError('preserve local binary cache')
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(data)
        if path.read_bytes() != data:
            raise ArithmeticError('local binary byte replay differs')
        rows.append({'family': family, 'path': str(path.relative_to(ROOT)),
                     'bytes': len(data), 'sha256': cert.hashed(path)})
    paths = [Path(__file__).resolve(), local.OUT]
    result = {'schema': 'elliptic-curves.mw16-local-binary-corrections.v1', 'status': 'PASS',
              'sources': {str(p.relative_to(ROOT)): cert.hashed(p) for p in paths}, 'rows': rows,
              'claim_boundary': 'Exact binary encoding of five complete local correction tables. '
                  'No change to an existing scanner, selection or point experiment.'}
    if check:
        if cert.read(OUT) != result:
            raise ArithmeticError('local encoding certificate differs')
    else:
        if OUT.exists():
            raise FileExistsError('preserve local binary certificate')
        checkpoint(OUT, result)
    print('FIVE LOCAL BINARY CACHES PASS', sum(r['bytes'] for r in rows), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--check', action='store_true')
    run(parser.parse_args().check)
