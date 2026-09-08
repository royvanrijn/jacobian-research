#!/usr/bin/env python3
"""Immutable intake package; replay selection and equation checks without new scoring."""
import argparse
import gzip
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path
import det1092_record_scale_selection as c
from research_runtime.store import checkpoint

OUT = c.ART / 'det1092_record_scale_intake_v1'
NAMES = ['protocol.json', 'tables.json', 'population-result.json',
         'selection-result.json', 'ledger.json']


def audit():
    p = c.protocol()
    pop = c.cert.read(c.D / 'population-result.json')
    selected = c.cert.read(c.D / 'selection-result.json')
    ledger = c.cert.read(c.D / 'ledger.json')
    table = c.cert.read(c.D / 'tables.json')
    assert all(d['status'] == 'PASS' for d in [pop, selected, ledger, table])
    for d in [pop, selected, table]:
        assert d['protocol_sha256'] == c.cert.hashed(c.D / 'protocol.json')
    assert selected['population_sha256'] == c.cert.hashed(c.D / 'population-result.json')
    assert pop['table_sha256'] == c.cert.hashed(c.D / 'tables.json')
    assert pop['parameters'] == p['population'] == sum(b['rows'] for b in pop['blocks'])
    assert sum(pop['height_bin_counts'].values()) == pop['parameters']
    for block in pop['blocks']:
        assert c.cert.hashed(c.ROOT / block['path']) == block['sha256']
    assert len(selected['rows']) == 2048
    retained = {r['id']: r for r in pop['retained']}
    ac, bc = c.coefficients()
    raw_files = {}
    for r in selected['rows']:
        assert all(r[k] == v for k, v in retained[r['id']].items())
        t = F(r['parameter'])
        a, b = [c.homogeneous(co, t.numerator, t.denominator) for co in [ac, bc]]
        assert r['model'] == ['0', '0', '0', str(a), str(b)]
        j = F(6912*a**3, 4*a**3+27*b*b)
        assert r['j_numerator_bits'] == abs(j.numerator).bit_length()
        assert r['j_denominator_bits'] == j.denominator.bit_length()
        path = c.D / 'extended-scores' / r['id'] / 'raw.json'
        assert c.cert.hashed(path) == r['raw_sha256']
        raw_files[str(path.relative_to(c.ROOT))] = r['raw_sha256']
    expected = []
    for band in p['height_bins']:
        rows = sorted((r for r in selected['rows'] if r['height_bin'] == band),
                      key=lambda r: (-r['score_units'], F(r['parameter']).denominator,
                                     F(r['parameter']).numerator, r['id']))
        n = len(rows)
        assert n == 1024 and selected['orders'][str(band)] == [r['id'] for r in rows]
        expected += [dict(r, stratum='strong') for r in rows[:16]]
        expected += [dict(r, stratum='moderate') for r in rows[n//3:n//3+4]]
        lower = sorted(rows[2*n//3:], key=lambda r: hashlib.sha256(
            ('det1092-record-scale-lower-v1:'+r['id']).encode()).hexdigest())
        expected += [dict(r, stratum='lower_fixed') for r in lower[:4]]
    assert selected['selected'] == expected and len(expected) == 48
    js = []
    for r in expected:
        a, b = map(int, r['model'][3:])
        js.append(F(6912*a**3, 4*a**3+27*b*b))
    assert len(set(js)) == 48
    stages = [s['supervision'] for s in ledger['stages']]
    assert all(s['outcome'] == 'completed' and s['returncode'] == 0 for s in stages)
    return dict(schema='elliptic-curves.det1092-record-scale-intake.v1', status='PASS',
                parameters=pop['parameters'], draws=pop['draws'],
                admitted=sum(pop['height_bin_counts'].get(str(k), 0) for k in p['height_bins']),
                retained=2048, selected=48, height_bin_counts=pop['height_bin_counts'],
                selected_j_numerator_bits=[min(r['j_numerator_bits'] for r in expected),
                                           max(r['j_numerator_bits'] for r in expected)],
                stage_seconds=sum(s['wall_seconds'] for s in stages),
                raw_score_files=raw_files,
                input_hashes={n: c.cert.hashed(c.D / n) for n in NAMES},
                checker_sha256=c.cert.hashed(Path(__file__)),
                scope='Immutable completed intake only. This checker replays selection, exact models/heights, source and raw-score hashes; it does not recompute all character sums or all million scores. Original scoring independently checked projective tables and retained scalar traces. No point outcome or novelty claim. Validation primes remain separate.')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--check', action='store_true')
    args = ap.parse_args()
    result = audit()
    if args.check:
        assert c.cert.read(OUT / 'manifest.json') == result
        for name in NAMES:
            assert gzip.decompress((OUT / (name+'.gz')).read_bytes()) == (c.D / name).read_bytes()
    else:
        assert not OUT.exists()
        OUT.mkdir()
        for name in NAMES:
            (OUT / (name+'.gz')).write_bytes(gzip.compress((c.D / name).read_bytes(), mtime=0))
        checkpoint(OUT / 'manifest.json', result)
    print('PASS intake:', result['parameters'], 'parameters;', result['selected'], 'fixed fibres')


if __name__ == '__main__':
    main()
