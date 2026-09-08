#!/usr/bin/env python3
"""Portable rank26 follow-up, latest-seven incidence, orbit and conductor diagnostics."""
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
        raise FileExistsError('preserve follow-up evidence')
    directory = LOCAL/'compact-r17-new26-followup-v1'
    rows = cert.read(directory/'verification/result.json')['rows']
    if len(rows) != 3 or any(r['status'] != 'PASS' for r in rows):
        raise ArithmeticError('all retained chart and cloud replays must pass')
    result = cert.read(directory/'candidate-00/result.json')
    cloud = cert.read(ART/'compact_r17_new26_followup_recorded_mod2_v1.json')
    if cloud['input_sha256'] != cert.hashed(directory/'candidate-00/result.json'):
        raise ArithmeticError('terminal cloud input differs')
    terminal = cert.read(directory/'terminal.json')
    if terminal['result_sha256'] != cert.hashed(directory/'candidate-00/result.json'):
        raise ArithmeticError('terminal follow-up changed')
    paths = {Path(__file__).resolve(), ROOT/'elliptic-curves/notes/NEW_RANK26_FOLLOWUP_2026-09-05.md',
             ROOT/'elliptic-curves/tests/test_pointed_orbit_compression.py'}
    for folder in ('compact-r17-new26-followup-v1', 'latest7-cross-family-incidence-v1',
                   'new-rank26-conductor-bound-v1', 'pointed-orbit-compression-audit-v1',
                   'compact-r17-omitted-generic-class-gap-v1'):
        paths.update(p for p in (LOCAL/folder).rglob('*') if p.is_file())
    paths.update(p for p in (LOCAL/'compact-r17-wide-portable-verification-v1').iterdir() if p.is_file())
    for name in ('compact_r17_new26_followup_recorded_mod2_v1.json', 'new_compact_rank26_proof_v1.json',
                 'compact_r17_wide_results_v1.json', 'new_high_rank_curve_index_v3.json',
                 'latest7_cross_family_j_incidence_v1.json', 'latest7_cross_family_j_incidence_replay_v1.json',
                 'pointed_orbit_compression_audit_v1.json'):
        p = ART/name; paths.add(p)
        paths.update(ROOT/s for s in cert.read(p).get('sources', {}))
    protocol = cert.read(directory/'protocol.json')
    paths.update(ROOT/s for s in protocol['sources'])
    paths.update((ROOT/protocol['input_path'], ROOT/protocol['generic_census_path']))
    for name in ('followup_compact_r17_new26.sage', 'replay_compact_r17_new26_followup.py',
                 'verify_compact_r17_new26_followup.py', 'audit_latest7_cross_family_incidence.sage',
                 'replay_latest7_cross_family_incidence.py', 'audit_new_rank26_conductor_bound.sage',
                 'audit_pointed_orbit_compression.py'):
        paths.add(ROOT/'elliptic-curves/cas'/name)
    seen = set()
    while True:
        todo = [p for p in paths-seen if p.suffix == '.json']
        if not todo: break
        for p in todo:
            seen.add(p); data = cert.read(p)
            if isinstance(data, dict):
                for key in ('sources', 'source_hashes'):
                    value = data.get(key, {})
                    if isinstance(value, dict): paths.update(ROOT/q for q in value)
    paths = dependencies(paths)
    members = []
    with zipfile.ZipFile(archive, 'x', compression=zipfile.ZIP_DEFLATED) as z:
        for p in sorted(paths):
            raw = p.read_bytes(); name = str(p.relative_to(ROOT))
            info = zipfile.ZipInfo(name, date_time=(2026, 9, 5, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED; info.external_attr = 0o644 << 16
            z.writestr(info, raw)
            members.append({'path': name, 'bytes': len(raw), 'sha256': sha256(raw).hexdigest()})
    with zipfile.ZipFile(archive) as z:
        if z.testzip() is not None or any(sha256(z.read(r['path'])).hexdigest() != r['sha256'] for r in members):
            raise ArithmeticError('archive member integrity failed')
    cert.write(output, {'schema': 'elliptic-curves.new-rank26-followup-evidence.v1',
        'builder_sha256': cert.hashed(Path(__file__).resolve()), 'archive': str(archive.relative_to(ROOT)),
        'archive_sha256': cert.hashed(archive), 'archive_bytes': archive.stat().st_size, 'files': members,
        'charts_replayed': len(result['charts']), 'planned_charts': protocol['charts'],
        'cloud_points_up_to_sign': len(cloud['points']), 'rank_lower_bound': cloud['rank_lower_bound'],
        'supervision_outcome': terminal['supervision']['outcome'],
        'claim_boundary': 'All retained adaptive charts and their complete point cloud replay; incomplete search coverage remains censored and supplies no upper rank. Includes exact latest-seven family incidence, a span-preserving orbit compression audit and a partial conductor bound. The completed initial24-candidate experiment has a separate bundle.'})
    print('PACKAGED NEW RANK26 FOLLOWUP', len(members), 'files;', archive.stat().st_size, 'bytes', flush=True)

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    package(p.parse_args().output.resolve())
