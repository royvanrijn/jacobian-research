#!/usr/bin/env python3
"""Retain the final301-chart new-MW16 followup and the fixed specialized parity audit."""
import argparse
from hashlib import sha256
from pathlib import Path
import shutil
import zipfile
import certify_compact_r17_candidates as cert
from package_recorded_mod2_audit import dependencies

ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT/'artifacts/local/elliptic-curves'
ART = ROOT/'artifacts/generated-results/elliptic-curves'

def package(output):
    archive = output.with_suffix('.zip')
    if output.exists() or archive.exists():
        raise FileExistsError('preserve prior diagnostics bundle')
    directory = LOCAL/'prospective-mw16-wide-rank25-followup-v1'
    final = directory/'candidate-00/result.json'
    data = cert.read(final)
    verification = cert.read(directory/'verification/results.json')
    if data['status'] != 'COMPLETE_DECLARED_PILOT' or len(data['charts']) != 301 or any(r['status'] != 'PASS' for r in verification['rows']) or len(verification['rows']) != 3:
        raise ArithmeticError('followup or full verification incomplete')
    sample = LOCAL/'specialized-rank25-parity-sample-v1'
    if cert.read(sample/'ledger.json')['status'] != 'COMPLETE_DECLARED_ATTEMPTS' or cert.read(sample/'replay-result.json')['status'] != 'PASS':
        raise ArithmeticError('specialized sample audit incomplete')
    paths = {Path(__file__).resolve(), final}
    for pattern in ('*.json', '*.log', 'continuation-*/*.json', 'continuation-*/*.log', 'verification/*'):
        paths.update(p for p in directory.glob(pattern) if p.is_file())
    paths.update(p for p in sample.iterdir() if p.is_file())
    paths.update(ROOT/name for name in cert.read(directory/'protocol.json')['sources'])
    paths.update(ROOT/name for name in cert.read(sample/'protocol.json')['sources'])
    paths.update(ART/name for name in ('prospective_mw16_wide_followup_mod2_v1.json', 'specialized_rank25_parity_sample_replay_v1.json'))
    for name in ('replay_prospective_mw16_followup.py', 'audit_recorded_point_mod2_rank_v2.py', 'replay_specialized_rank25_parities.py', 'package_recorded_mod2_audit.py'):
        paths.add(ROOT/'elliptic-curves/cas'/name)
    paths = dependencies(paths)
    rows = []
    with zipfile.ZipFile(archive, 'x', compression=zipfile.ZIP_DEFLATED) as z:
        for p in sorted(paths):
            name = str(p.relative_to(ROOT))
            info = zipfile.ZipInfo(name, date_time=(2026,9,5,0,0,0))
            info.compress_type = zipfile.ZIP_DEFLATED; info.external_attr = 0o644 << 16
            h = sha256()
            with p.open('rb') as source, z.open(info, 'w') as target:
                for chunk in iter(lambda:source.read(1 << 20), b''):
                    target.write(chunk); h.update(chunk)
            rows.append({'path':name, 'bytes':p.stat().st_size, 'sha256':h.hexdigest()})
    with zipfile.ZipFile(archive) as z:
        if z.testzip() is not None:
            raise ArithmeticError('archive integrity failed')
        for row in rows:
            h = sha256()
            with z.open(row['path']) as source:
                for chunk in iter(lambda:source.read(1 << 20), b''):
                    h.update(chunk)
            if h.hexdigest() != row['sha256']:
                raise ArithmeticError('member hash differs')
    proof = cert.read(ART/'prospective_mw16_wide_followup_mod2_v1.json')
    cert.write(output, {'schema':'elliptic-curves.rank25-followup-diagnostics.v2',
        'builder_sha256':cert.hashed(Path(__file__).resolve()), 'archive':str(archive.relative_to(ROOT)),
        'archive_sha256':cert.hashed(archive), 'archive_bytes':archive.stat().st_size, 'files':rows,
        'family':data['family'], 'parameter':data['parameter'], 'charts_replayed':301,
        'rank_lower_bound':proof['rank_lower_bound'], 'recorded_points_up_to_sign':len(proof['points']),
        'specialized_parity_representatives_checked':6144,
        'claim_boundary':'Final301-chart transcript, all continuation limits and exact replays, complete-cloud mod2 proof, and three fixed specialized parity samples. Previous217/290-chart snapshots remain preserved by hash locally; the final superset is bundled. No rank upper bound, new rank gain, or CVP optimality.'})
    print('PACKAGED RANK25 FOLLOWUP DIAGNOSTICS', len(rows), 'files;', archive.stat().st_size, 'bytes', flush=True)

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    package(p.parse_args().output.resolve())
