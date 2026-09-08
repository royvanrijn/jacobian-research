#!/usr/bin/env python3
"""Read-only complete address/retention audit of the frozen million-row intake."""
import gzip
import hashlib
import heapq
import json
import math
from collections import Counter
from fractions import Fraction as F
import numpy as np
import det1092_record_scale_selection as c
from research_runtime.store import checkpoint

OUT = c.ART / 'det1092_record_scale_population_audit_v1.json'


def main():
    assert not OUT.exists()
    p = c.protocol()
    pop = c.cert.read(c.D / 'population-result.json')
    tab = c.cert.read(c.D / 'tables.json')
    assert pop['status'] == tab['status'] == 'PASS'
    hmax = p['parameter_height_max']
    seen, heaps, counts = set(), {k: [] for k in p['height_bins']}, Counter()
    draw = index = scaled = 0
    tables = []
    for t in tab['tables']:
        q = t['prime']
        weights = [0 if a is None else round(10**12*(2-a)*math.log(q)/(q+1-a))
                   for a in t['traces']]
        tables.append((q, np.array([0]+[pow(v, -1, q) for v in range(1, q)], dtype=np.int64),
                       np.array(weights, dtype=np.int64), np.array(t['ambiguous'], dtype=bool)))
    ac, bc = c.coefficients()
    for block in pop['blocks']:
        path = c.ROOT / block['path']
        assert c.cert.hashed(path) == block['sha256']
        rows = json.loads(gzip.decompress(path.read_bytes()))
        assert len(rows) == block['rows']
        models, admitted = [], []
        for r in rows:
            while True:
                assert draw < p['maximum_draws']
                digest = hashlib.sha256((p['sample_domain']+':'+str(draw)).encode()).digest()
                m = int.from_bytes(digest[:8], 'big') % (2*hmax+1)-hmax
                n = int.from_bytes(digest[8:16], 'big') % hmax+1
                g = math.gcd(m, n)
                m, n = m//g, n//g
                source_draw = draw
                draw += 1
                if m and max(abs(m), n) >= p['parameter_height_min'] and (m, n) not in seen:
                    seen.add((m, n))
                    break
            assert r[:4] == [index, source_draw, m, n]
            index += 1
            a = c.homogeneous(ac, m, n)
            b = c.homogeneous(bc, m, n)
            j = F(1728*4*a**3, 4*a**3+27*b*b)
            assert r[4:6] == [abs(j.numerator).bit_length(), j.denominator.bit_length()]
            band = r[4]//p['height_bin_width']
            counts[str(band)] += 1
            if band in heaps:
                admitted.append(r)
                models.append((a, b))
            else:
                assert r[6] is None
        nums = np.array([r[2] for r in admitted], dtype=np.int64)
        dens = np.array([r[3] for r in admitted], dtype=np.int64)
        scores = np.zeros(len(admitted), dtype=np.int64)
        for q, inverses, weights, ambiguous in tables:
            residues = np.where(dens % q == 0, q, (nums % q)*inverses[dens % q] % q)
            scores += weights[residues]
            for k in np.flatnonzero(ambiguous[residues]):
                a, b = models[k]
                if a % q**4 == 0 and b % q**6 == 0:
                    # Fail closed rather than reuse the original fallback unnoticed.
                    scaled += 1
                    raise AssertionError('An actual local scaling requires an independent scalar fallback audit')
        for r, score in zip(admitted, scores):
            # Operation order in rounded heuristic weights is part of the frozen protocol.
            assert int(score) == r[6], (r[0], int(score), r[6])
            key = (int(score), -r[3], -r[2], -r[0])
            heap = heaps[r[4]//p['height_bin_width']]
            if len(heap) < p['retain_per_bin']:
                heapq.heappush(heap, key)
            elif key > heap[0]:
                heapq.heapreplace(heap, key)
        if index % 131072 == 0:
            print('AUDITED', index, flush=True)
    assert index == pop['parameters'] == p['population']
    assert draw == pop['draws'] and dict(counts) == pop['height_bin_counts']
    for band, heap in heaps.items():
        ids = ['scale-'+str(-r[3]).zfill(7) for r in sorted(heap, reverse=True)]
        retained = [r for r in pop['retained'] if r['height_bin'] == band]
        assert ids == [r['id'] for r in retained]
        # Independent direct polynomial evaluation on all retained equations.
        for r in retained:
            t = F(r['parameter']);m,n = t.numerator,t.denominator
            values = [sum(v*m**i*n**(len(co)-1-i) for i,v in enumerate(co)) for co in [ac,bc]]
            assert list(map(str, values)) == r['model'][3:]
    checkpoint(OUT, dict(schema='elliptic-curves.det1092-record-scale-population-audit.v1',
        status='PASS', parameters=index, draws=draw, retained=2048,
        height_bin_counts=dict(counts), actual_local_scalings=scaled,
        protocol_sha256=c.cert.hashed(c.D/'protocol.json'),
        population_sha256=c.cert.hashed(c.D/'population-result.json'),
        tables_sha256=c.cert.hashed(c.D/'tables.json'),
        checker_sha256=c.cert.hashed(__import__('pathlib').Path(__file__)),
        scope='Reconstruct every SHA draw including rejection/deduplication, exact j heights and all admitted table scores; independently rebuild both retention heaps and directly evaluate retained polynomials. Trace tables were independently checked in the original run. This does not rerun ellcard, point searches, or introduce new addresses. Actual local scaling fails closed for a separate audit. No rank inference.'))
    print('PASS complete population and first-stage retention audit', flush=True)


if __name__ == '__main__':
    main()
