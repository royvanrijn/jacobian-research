#!/usr/bin/env python3
"""Retain the censored follow-up, complete replays and failed descent diagnostic."""
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
    directory = LOCAL/'prospective-mw16-small-conductor-followup-v2'
    rows = cert.read(directory/'verification/result.json')['rows']
    if len(rows) != 2 or any(r['status'] != 'PASS' for r in rows):
        raise ArithmeticError('all retained chart and cloud replays must pass')
    result = cert.read(directory/'candidate-00/result.json')
    cloud = cert.read(ART/'small_conductor_followup_recorded_mod2_v1.json')
    if cloud['input_sha256'] != cert.hashed(directory/'candidate-00/result.json'):
        raise ArithmeticError('terminal cloud input differs')
    if cert.read(directory/'supervision-result.json')['outcome'] != 'strict_wall_timeout':
        raise ArithmeticError('fixed censored run differs')
    if cert.read(LOCAL/'small-conductor-descent-v1/supervision-result.json')['outcome'] != 'strict_wall_timeout':
        raise ArithmeticError('fixed descent outcome differs')
    paths = {Path(__file__).resolve(), ROOT/'elliptic-curves/notes/SMALL_CONDUCTOR_FOLLOWUP_2026-09-05.md',
             ROOT/'elliptic-curves/notes/NEW_SMALL_CONDUCTOR_CURVE_2026-09-05.md',
             ROOT/'elliptic-curves/tests/test_cached_observation_state.py',
             ART/'new_rank22_small_conductor_curve.sage',
             LOCAL/'small-conductor-python-proof-v1/exported-curve-replay.log'}
    for folder in ('prospective-mw16-small-conductor-followup-v1', 'prospective-mw16-small-conductor-followup-v2',
                   'small-conductor-root-number-v1', 'small-conductor-descent-v1',
                   'retained-admission-profile-v1', 'cached-observation-validation-v1',
                   'refreshed-inventory-replay-v1', 'next12-current-catalogue-v1'):
        paths.update(p for p in (LOCAL/folder).rglob('*') if p.is_file())
    paths.update(p for p in (LOCAL/'next12-portable-verification-v1').iterdir() if p.is_file())
    for name in ('small_conductor_followup_recorded_mod2_v1.json', 'cached_observation_state_validation_v1.json',
                 'refreshed_new_curve_inventory_replay_v1.json', 'small_conductor_rank22_proof_v1.json',
                 'prospective_mw16_next12_results_v1.json', 'new_high_rank_curve_index_v2.json'):
        p = ART/name; paths.add(p)
        paths.update(ROOT/s for s in cert.read(p).get('sources', {}))
    protocol = cert.read(directory/'protocol.json')
    paths.update(ROOT/s for s in protocol['sources'])
    paths.update((ROOT/protocol['input_path'], ROOT/protocol['generic_census_path']))
    for name in ('followup_small_conductor_mw16.sage', 'followup_small_conductor_mw16_v2.sage',
                 'replay_prospective_mw16_followup.py', 'replay_prospective_mw16_followup_v2.py',
                 'audit_small_conductor_root_number.sage', 'audit_small_conductor_descent.sage',
                 'profile_retained_mwstate_admission.py', 'audit_cached_observation_state.py',
                 'replay_refreshed_curve_inventory.py'):
        paths.add(ROOT/'elliptic-curves/cas'/name)
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
    cert.write(output, {'schema': 'elliptic-curves.small-conductor-followup-evidence.v1',
        'builder_sha256': cert.hashed(Path(__file__).resolve()), 'archive': str(archive.relative_to(ROOT)),
        'archive_sha256': cert.hashed(archive), 'archive_bytes': archive.stat().st_size, 'files': members,
        'charts_replayed': len(result['charts']), 'planned_charts': protocol['charts'],
        'cloud_points_up_to_sign': len(cloud['points']), 'rank_lower_bound': cloud['rank_lower_bound'],
        'claim_boundary': 'The127-chart retained prefix and7753-point cloud replay exactly, still at lower bound22. The301-chart plan is censored; no box is exhausted. Root number+1 is analytic information, and the descent initialization timeout gives no rank upper bound. The optional exact membership cache preserves state records. Earlier completed next12 experiments have their separate evidence bundle.'})
    print('PACKAGED SMALL CONDUCTOR FOLLOWUP', len(members), 'files;', archive.stat().st_size, 'bytes', flush=True)

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    package(p.parse_args().output.resolve())
