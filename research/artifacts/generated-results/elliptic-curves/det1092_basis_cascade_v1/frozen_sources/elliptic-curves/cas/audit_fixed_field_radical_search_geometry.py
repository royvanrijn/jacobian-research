#!/usr/bin/env sage
"""Exact, offline audit of the six recorded PointsQI height boxes.

Default: verify a pinned certificate using rational Sturm sequences and
integer inequalities. --write-certificate uses algebraic root isolation only
to propose intervals; the same rational verifier must accept the result.
This excludes a bounded coordinate box, never a rational torsor globally.
"""
from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

from sage.all import AA, QQ, ZZ, PolynomialRing, RealIntervalField, matrix
from sage.version import version as sage_version
import run_fixed_field_radical_covers as models

ROOT = models.ROOT
OUTPUT = ROOT / 'artifacts/generated-results/elliptic-curves/fixed_field_radical_search_geometry_v1.json'
R = PolynomialRing(QQ, 'x')
BOUND = ZZ(10)**7


def sturm(f):
    seq = [f, f.derivative()]
    while seq[-1]:
        seq.append(-seq[-2].quo_rem(seq[-1])[1])
    return seq[:-1]


def variations(values):
    signs = [1 if v > 0 else -1 for v in values if v]
    return sum(a != b for a, b in zip(signs, signs[1:]))


def root_count(f, lo=None, hi=None):
    """Sturm count on (lo, hi), with None denoting the relevant infinity."""
    seq = sturm(f)
    if lo is not None:
        assert f(lo), 'root at left endpoint'
    if hi is not None:
        assert f(hi), 'root at right endpoint'
    left = [p(lo) for p in seq] if lo is not None else [
        p.leading_coefficient() * (-1)**p.degree() for p in seq]
    right = [p(hi) for p in seq] if hi is not None else [
        p.leading_coefficient() for p in seq]
    return variations(left) - variations(right)


def geometry(record):
    row = record['cover']
    assert not any(map(QQ, row['reduced_quartic_model'][1]))
    f = R(row['reduced_quartic_model'][0])
    assert f.degree() == 4 and all(c in ZZ for c in f)
    assert f.leading_coefficient() > BOUND**2
    assert not ZZ(f.leading_coefficient()).is_square(), 'infinity needs separate treatment'
    S = matrix(QQ, 4, record['quadric_model']['variable_transform'])
    assert all(c in ZZ for c in S.list()) and abs(S.det()) == 1
    assert list(S.column(3)) == [1, 0, 0, 0], 'old ordinate must be new first coordinate'
    z_norm = ZZ(sum(abs(c) for c in S.column(2)))
    D = (BOUND * z_norm).isqrt()
    return f, int(D), int(z_norm)


def propose(record):
    f, D, z_norm = geometry(record)
    row = record['cover']
    result = {'anchor_mask': row['mask'],
              'translated_by_Q': row['translated_by_universal_point'],
              'old_Z_column_l1_norm': z_norm, 'parameter_denominator_bound': D}
    g = f - BOUND**2
    if root_count(g) == 0:
        # These six fixed inputs have a much stronger uniform real bound.
        powers = [k for k in range(28, 34)
                  if f(0) > 10**k and root_count(f - 10**k) == 0]
        assert powers
        result['uniform_f_lower_bound'] = str(10**max(powers))
        result['root_bands'] = []
    else:
        assert root_count(f) == root_count(g) == 4
        roots = [f.roots(AA, multiplicities=False), g.roots(AA, multiplicities=False)]
        dyadic = ZZ(2)**128
        bands = []
        for a, b in zip(*roots):
            intervals = [RealIntervalField(256)(r) for r in (a, b)]
            lo = min(v.lower().exact_rational() for v in intervals)
            hi = max(v.upper().exact_rational() for v in intervals)
            lo = (lo * dyadic).floor() / dyadic
            hi = (hi * dyadic).ceil() / dyadic
            bands.append([str(lo), str(hi)])
        result['uniform_f_lower_bound'] = None
        result['root_bands'] = bands
    return result


def verify_row(record, certificate):
    f, D, z_norm = geometry(record)
    row = record['cover']
    assert certificate['anchor_mask'] == row['mask']
    assert certificate['translated_by_Q'] == row['translated_by_universal_point']
    assert certificate['old_Z_column_l1_norm'] == z_norm
    assert certificate['parameter_denominator_bound'] == D
    g = f - BOUND**2
    lower = certificate['uniform_f_lower_bound']
    if lower is not None:
        lower = ZZ(lower)
        assert lower >= BOUND**2
        assert f(0) > lower and root_count(f - lower) == 0
        assert certificate['root_bands'] == []
        return 0
    bands = [tuple(map(QQ, band)) for band in certificate['root_bands']]
    assert len(bands) == 4 and root_count(f) == root_count(g) == 4
    assert f.gcd(f.derivative()).degree() == g.gcd(g.derivative()).degree() == 0
    for i, (lo, hi) in enumerate(bands):
        assert lo < hi
        assert root_count(f, lo, hi) == root_count(g, lo, hi) == 1
        if i:
            previous = bands[i-1][1]
            assert previous < lo
            q = (previous + lo) / 2
            assert f(q) < 0 or g(q) > 0, 'a feasible complementary interval was omitted'
        # Exact integer arithmetic: no rational with denominator <= D lies
        # in this closed band. No square tests or point enumeration occur.
        for denominator in range(1, D+1):
            assert (lo * denominator).ceil() > (hi * denominator).floor()
    # Both tails have g > 0 by its positive leading coefficient and the
    # absence of any further real roots. Thus 0 <= f <= B^2 is covered.
    return len(bands) * D


def inputs():
    summary = json.loads(models.SUMMARY.read_text())
    assert models.sha(models.EVIDENCE) == summary['evidence_sha256']
    for path, digest in summary['source_hashes'].items():
        assert models.sha(ROOT / path) == digest
    evidence = json.loads(gzip.decompress(models.EVIDENCE.read_bytes()))
    assert evidence['source_sha256'] == models.sha(models.base.SOURCE)
    assert evidence['pairing_summary_sha256'] == models.sha(models.CT)
    _, _, _, _, E = models.base.context(models.base.SOURCE, -1)
    records = evidence['covers']
    assert len(records) == 6
    assert {(r['cover']['mask'], r['cover']['translated_by_universal_point'])
            for r in records} == {(m, bool(t)) for m in models.approved_masks() for t in (0, 1)}
    for record in records:
        models.verify_quadric_model(record['cover'], record['quadric_model'], E)
    return records


def verify(records, certificate):
    assert certificate['schema'] == 'elliptic-curves.fixed-field-radical-search-geometry.v1'
    assert certificate['status'] == 'SIX_NOMINAL_HEIGHT_BOXES_EXACTLY_EMPTY'
    assert certificate['naive_projective_height_bound'] == BOUND
    assert certificate['source_evidence_sha256'] == models.sha(models.EVIDENCE)
    assert certificate['checker_sha256'] == models.sha(__file__)
    assert certificate['point_or_sha'] == {str(m): 'UNKNOWN' for m in models.MASKS}
    assert len(certificate['models']) == len(records)
    checked = sum(verify_row(r, c) for r, c in zip(records, certificate['models']))
    assert certificate['rational_denominator_inequalities'] == checked
    return checked


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-certificate', action='store_true')
    args = parser.parse_args()
    records = inputs()
    if args.write_certificate:
        rows = [propose(record) for record in records]
        certificate = {
            'schema': 'elliptic-curves.fixed-field-radical-search-geometry.v1',
            'status': 'SIX_NOMINAL_HEIGHT_BOXES_EXACTLY_EMPTY',
            'source_evidence': str(models.EVIDENCE.relative_to(ROOT)),
            'source_evidence_sha256': models.sha(models.EVIDENCE),
            'checker_sha256': models.sha(__file__), 'sage_version': sage_version,
            'naive_projective_height_bound': int(BOUND),
            'models': rows,
            'rational_denominator_inequalities': sum(verify_row(r, c) for r, c in zip(records, rows)),
            'point_or_sha': {str(m): 'UNKNOWN' for m in models.MASKS},
            'claim_boundary': [
                'Only primitive integral projective coordinates of height <= 10000000 are excluded.',
                'PointsQI used ExactBound=false and may also explore points outside this box.',
                'No global rational-point, Sha, full Selmer or rank conclusion is obtained.']}
        verify(records, certificate)
        models.base.save(OUTPUT, certificate)
    certificate = json.loads(OUTPUT.read_text())
    count = verify(records, certificate)
    print(f'PASS_SIX_EMPTY_NOMINAL_BOXES: {count} exact denominator inequalities; all three classes UNKNOWN')


if __name__ == '__main__':
    main()
