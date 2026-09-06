#!/usr/bin/env python3
"""Package the six fresh censuses, balanced R17 search and new rank26 proof."""
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
        raise FileExistsError('preserve R17 evidence')
    directory = LOCAL/'compact-six-r17-h4096-v1'
    ledger = cert.read(directory/'verification/ledger.json')
    if len(ledger['rows']) != 24 or any(r['status'] != 'PASS' for r in ledger['rows']):
        raise ArithmeticError('all24 chart and cloud checks required')
    for path in (directory/'certification/check.supervisor.json',
                 LOCAL/'new-rank26-publication-v1/curve-v2/check.supervisor.json',
                 LOCAL/'new-rank26-publication-v1/index/check.supervisor.json'):
        row = cert.read(path)
        if row['outcome'] != 'completed' or row['returncode'] != 0:
            raise ArithmeticError('standalone proof replay missing')
    paths = {Path(__file__).resolve(), ROOT/'elliptic-curves/notes/NEW_COMPACT_RANK26_CURVE_2026-09-05.md'}
    for folder in (directory, LOCAL/'compact-r17-fresh-generic-census-v1',
                   LOCAL/'new-rank26-publication-v1', LOCAL/'next12-current-catalogue-v1'):
        paths.update(p for p in folder.rglob('*') if p.is_file())
    names = ['compact_r17_wide_results_v1.json', 'new_compact_rank26_proof_v1.json',
             'compact_r17_wide_rank26_initial_proof_v1.json', 'new_high_rank_curve_index_v3.json',
             'new_high_rank_curve_index_v2.json', 'new_high_rank_curve_index_v1.json',
             'compact_r17_wide_selection_audit_v1.json']
    names += list(cert.read(ART/'new_high_rank_curve_index_v3.json')['source_certificate_hashes'])
    families = cert.read(directory/'protocol.json')['families']
    names += [f'compact_r17_wide_recorded_mod2_{f}_{i:02}_v1.json' for f in families for i in range(4)]
    names += [f'compact_r17_generic_census_{f}_replay_v1.json' for f in families]
    paths.update(ART/name for name in names)
    paths.update(ART/name for name in ('new_high_rank_curve_index_v3.csv', 'new_compact_rank26_curve.sage'))
    parent = LOCAL/'compact-six-r17-h1024-v2'
    paths.add(parent/'trace-table-check.json')
    for f in families:
        paths.update(parent/f/f'tables-{sign}.txt' for sign in (-1, 1))
        paths.update(p for p in (parent/f/'trace-tables').rglob('*') if p.is_file())
    paths.add(LOCAL/'compact-r17-wide-v1/scanner')
    paths.update(ROOT/'elliptic-curves/cas'/name for name in (
        'census_compact_r17_generic_parities.sage', 'replay_compact_r17_generic_census.py',
        'select_compact_r17_wide.py', 'search_compact_r17_wide.sage', 'run_compact_r17_wide.py',
        'replay_compact_r17_wide.py', 'verify_compact_r17_wide_batch.py',
        'audit_compact_r17_wide_selection.py', 'certify_compact_r17_wide_results.py',
        'certify_new_compact_rank26.py', 'certify_new_compact_rank26_v2.py',
        'export_new_compact_rank26_sage.py', 'export_new_high_rank_curve_index_v3.py',
        'newfamily/scan_rational_nagao_tables.cpp'))
    # Include all directly source-hashed dependencies, including protocols
    # which were frozen before workers started. Follow the closure once each.
    seen = set()
    while True:
        todo = [p for p in paths-seen if p.suffix == '.json']
        if not todo:
            break
        for p in todo:
            seen.add(p)
            data = cert.read(p)
            if isinstance(data, dict):
                for key in ('sources', 'source_hashes'):
                    value = data.get(key, {})
                    if isinstance(value, dict):
                        paths.update(ROOT/s for s in value)
    paths = dependencies(paths)
    missing = [str(p) for p in paths if not p.is_file()]
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
    cert.write(output, {'schema': 'elliptic-curves.compact-r17-wide-evidence.v1',
        'builder_sha256': cert.hashed(Path(__file__).resolve()), 'archive': str(archive.relative_to(ROOT)),
        'archive_sha256': cert.hashed(archive), 'archive_bytes': archive.stat().st_size, 'files': rows,
        'generic_parities_checked': 786432, 'scored_addresses': 122400468,
        'fixed_point_attempts': 24, 'charts_replayed': 1032, 'complete_cloud_point_occurrences': 6845,
        'new_high_rank_curves': 3, 'inventory_curves': 39,
        'claim_boundary': 'Fresh finite generic parity witnesses, bounded selection and all24 point measurements; exact rank lower bounds and a new globally minimal rank-at-least26 curve. No exact rank, conductor or universal novelty. The separate adaptive follow-up is not included.'})
    print('PACKAGED R17 H4096', len(rows), 'files;', archive.stat().st_size, 'bytes', flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, default=ART/'compact_r17_wide_evidence_v1.json')
    package(p.parse_args().output.resolve())
