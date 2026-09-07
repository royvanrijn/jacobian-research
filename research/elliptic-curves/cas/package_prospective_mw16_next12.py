#!/usr/bin/env python3
"""Portable evidence for the fixed next12 batch and its small-conductor curve."""
import argparse
from hashlib import sha256
from pathlib import Path
import zipfile
import certify_compact_r17_candidates as cert
from package_recorded_mod2_audit import dependencies

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT/'artifacts/generated-results/elliptic-curves'
LOCAL = ROOT/'artifacts/local/elliptic-curves'

def package(output):
    archive = output.with_suffix('.zip')
    if output.exists() or archive.exists():
        raise FileExistsError('preserve next12 evidence')
    directory = LOCAL/'prospective-mw16-f5-next12-v1'
    ledger = cert.read(directory/'verification/ledger.json')
    if len(ledger['rows']) != 12 or any(s['status'] != 'PASS' for r in ledger['rows'] for s in r['rows']):
        raise ArithmeticError('twelve chart and full-cloud replays required')
    if any(r['status'] != 'PASS' for r in cert.read(directory/'verification/certificate-results.json')['rows']):
        raise ArithmeticError('batch point certificates not replayed')
    if any(r['outcome'] != 'completed' or r['returncode'] != 0 for r in cert.read(LOCAL/'small-conductor-python-proof-v1/result.json')['rows']):
        raise ArithmeticError('conductor Python proof not replayed')
    paths = {Path(__file__).resolve(), ROOT/'elliptic-curves/notes/NEW_SMALL_CONDUCTOR_CURVE_2026-09-05.md'}
    paths.update(p for p in directory.rglob('*') if p.is_file())
    names = ['prospective_mw16_next12_results_v1.json', 'new_high_rank_curve_index_v1.json',
             'new_high_rank_curve_index_v2.json', 'next12_rank22_exact_conductor_v1.json',
             'next12_conductor_lucas_primes_v1.json', 'small_conductor_rank22_proof_v1.json']
    names += [f'prospective_mw16_next12_recorded_mod2_{i:02}_v1.json' for i in range(12)]
    names += list(cert.read(ART/'new_high_rank_curve_index_v2.json')['source_certificate_hashes'])
    for name in set(names):
        p = ART/name; paths.add(p)
        for key in ('sources', 'source_hashes'):
            paths.update(ROOT/s for s in cert.read(p).get(key, {}))
    paths.add(ART/'new_high_rank_curve_index_v2.csv')
    for folder in ('next12-new-curve-conductor-bounds-v1', 'next12-rank22-exact-conductor-v1',
                   'next12-conductor-lucas-v1', 'small-conductor-python-proof-v1',
                   'next12-current-catalogue-v1', 'new-high-rank-index-v2'):
        paths.update(p for p in (LOCAL/folder).rglob('*') if p.is_file())
    paths.update(ROOT/s for s in cert.read(directory/'point-protocol.json')['sources'])
    parent = LOCAL/'prospective-mw16-h1024-v1'
    for name in ('protocol.json', 'point-protocol.json', 'point-ledger.json'):
        paths.add(parent/name)
    # This batch reuses its parent's ranking. Retain that family’s complete
    # selection inputs; prior point campaigns remain in their own bundles.
    family = parent/'a1-fibration-05'
    paths.update(p for p in family.iterdir() if p.is_file())
    paths.update(p for p in (family/'trace-tables').rglob('*') if p.is_file())
    paths.update(ROOT/'elliptic-curves/cas'/name for name in (
        'run_prospective_mw16_next12.py', 'replay_prospective_mw16_search.py',
        'audit_recorded_point_mod2_rank_v2.py', 'export_new_high_rank_curve_index_v2.py',
        'select_prospective_mw16_atlas.sage', 'search_prospective_mw16_atlas.sage',
        'newfamily/scan_rational_nagao_tables.cpp', 'certify_prospective_mw16_next12_results.py',
        'build_next12_conductor_prime_certificates.sage', 'certify_small_conductor_curve.py',
        'audit_next12_conductor_bounds.sage', 'certify_next12_rank22_conductor.sage'))
    paths = dependencies(paths)
    missing = [str(p.relative_to(ROOT)) for p in paths if not p.is_file()]
    if missing:
        raise FileNotFoundError(missing)
    rows = []
    with zipfile.ZipFile(archive, 'x', compression=zipfile.ZIP_DEFLATED) as z:
        for p in sorted(paths):
            raw = p.read_bytes(); name = str(p.relative_to(ROOT))
            info = zipfile.ZipInfo(name, date_time=(2026, 9, 5, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED; info.external_attr = 0o644 << 16
            z.writestr(info, raw)
            rows.append({'path': name, 'bytes': len(raw), 'sha256': sha256(raw).hexdigest()})
    with zipfile.ZipFile(archive) as z:
        if z.testzip() is not None or any(sha256(z.read(r['path'])).hexdigest() != r['sha256'] for r in rows):
            raise ArithmeticError('archive integrity failed')
    cert.write(output, {'schema': 'elliptic-curves.prospective-mw16-next12-evidence.v1',
        'builder_sha256': cert.hashed(Path(__file__).resolve()), 'archive': str(archive.relative_to(ROOT)),
        'archive_sha256': cert.hashed(archive), 'archive_bytes': archive.stat().st_size, 'files': rows,
        'curves': 12, 'charts_replayed': 516, 'new_inventory_curves': 4,
        'complete_cloud_point_occurrences': 5940,
        'claim_boundary': 'Fixed next12 batch, exact lower-bound proofs, 36-curve index and independent rank22/global-minimal/exact-conductor proof. Includes all original and continued batch transcripts and parent family05 selection inputs. Subsequent small-conductor follow-up is separate. Catalogue586 placement is relative to recorded conductors, not universal novelty or an exact-rank claim.'})
    print('PACKAGED NEXT12', len(rows), 'files;', archive.stat().st_size, 'bytes', flush=True)

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    package(p.parse_args().output.resolve())
