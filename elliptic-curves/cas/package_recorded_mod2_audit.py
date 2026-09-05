#!/usr/bin/env python3
"""Portable point-cloud lower-bound proofs; original search snapshots stay separate."""
import argparse
import ast
from hashlib import sha256
from pathlib import Path
import zipfile
import certify_compact_r17_candidates as cert

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT/'artifacts/generated-results/elliptic-curves'
LOCAL = ROOT/'artifacts/local/elliptic-curves'

def dependencies(paths):
    paths = set(paths)
    paths.update((ROOT/'elliptic-curves/cas/research_runtime').glob('*.py'))
    todo = list(paths)
    while todo:
        path = todo.pop()
        if path.suffix != '.py':
            continue
        for node in ast.walk(ast.parse(path.read_text())):
            names = [n.name for n in node.names] if isinstance(node, ast.Import) else ([node.module] if isinstance(node, ast.ImportFrom) and node.module else [])
            for name in names:
                for base in (path.parent, ROOT/'elliptic-curves/cas', ROOT/'elliptic-curves', ROOT):
                    relative = Path(*name.split('.'))
                    for candidate in (base/relative.with_suffix('.py'), base/relative/'__init__.py'):
                        if candidate.is_file() and candidate not in paths:
                            paths.add(candidate); todo.append(candidate)
    return paths

def package(output):
    archive = output.with_suffix('.zip')
    if output.exists() or archive.exists():
        raise FileExistsError('preserve prior audit bundle')
    paths = {Path(__file__).resolve()}
    cases = []
    for dirname, rank_key in [('recorded-high-rank-mod2-recertification-v2', 'rank_lower_bound'),
                              ('recorded-low-rank-mod2-recertification-v1', 'original_rank_lower_bound')]:
        directory = LOCAL/dirname
        ledger = cert.read(directory/'ledger.json')
        if ledger['status'] != 'COMPLETE_DECLARED_AUDIT' or any(r['status'] != 'PASS' for r in ledger['rows']):
            raise ArithmeticError('cohort or certificate replay unfinished')
        paths.update(p for p in directory.iterdir() if p.is_file())
        for row in ledger['rows']:
            proof = ROOT/row['output']
            if cert.hashed(proof) != row['output_sha256']:
                raise ArithmeticError('proof changed after replay')
            data = cert.read(proof)
            paths.add(proof)
            paths.update(ROOT/p for p in data['sources'])
            cases.append({'input':row['input'], 'input_sha256':row['sha256'],
                'input_status':row['input_status'], 'certificate':row['output'],
                'point_count':row['point_count'], 'old_lower_bound':row[rank_key],
                'lower_bound':row['new_rank_lower_bound']})
    for name in ('curve542_initial_visibility_v1.json', 'curve542_translated_visibility_v1.json',
                 'curve542_translated_visibility_replay_v1.json', 'quotient_only_cache_validation_v1.json'):
        paths.add(ART/name)
    for p in (ART/'curve542_initial_visibility_v1.json', ART/'quotient_only_cache_validation_v1.json'):
        paths.update(ROOT/name for name in cert.read(p)['sources'])
    paths.update(ROOT/'elliptic-curves/cas'/name for name in ('run_recorded_point_mod2_audit.py',
        'audit_recorded_point_mod2_rank.py', 'audit_recorded_point_mod2_rank_v2.py',
        'replay_curve542_translated_visibility.py', 'audit_curve542_translated_visibility.sage'))
    paths.update(ART/f'recorded_rank25_mod2_recertification_{i}_v1.json' for i in range(4))
    for dirname in ('recorded-rank25-mod2-recertification-v1', 'curve542-translated-visibility-replay-v1', 'quotient-only-cache-validation-v1'):
        paths.update(p for p in (LOCAL/dirname).iterdir() if p.is_file())
    paths = dependencies(paths)
    rows = []
    with zipfile.ZipFile(archive, 'x', compression=zipfile.ZIP_DEFLATED) as z:
        for p in sorted(paths):
            raw = p.read_bytes(); name = str(p.relative_to(ROOT))
            info = zipfile.ZipInfo(name, date_time=(2026,9,5,0,0,0))
            info.compress_type = zipfile.ZIP_DEFLATED; info.external_attr = 0o644 << 16
            z.writestr(info, raw)
            rows.append({'path':name, 'bytes':len(raw), 'sha256':sha256(raw).hexdigest()})
    with zipfile.ZipFile(archive) as z:
        if z.testzip() is not None:
            raise ArithmeticError('archive integrity failed')
        for row in rows:
            if sha256(z.read(row['path'])).hexdigest() != row['sha256']:
                raise ArithmeticError('archive member hash changed')
    cert.write(output, {'schema':'elliptic-curves.recorded-mod2-audit-evidence.v1',
        'builder_sha256':cert.hashed(Path(__file__).resolve()),
        'archive':str(archive.relative_to(ROOT)), 'archive_sha256':cert.hashed(archive),
        'archive_bytes':archive.stat().st_size, 'files':rows, 'cases':cases,
        'transcripts':len(cases), 'point_occurrences_up_to_sign':sum(r['point_count'] for r in cases),
        'upgrades':[r for r in cases if r['lower_bound'] > r['old_lower_bound']],
        'claim_boundary':'Standalone rational point and exact finite quotient lower bounds, with all successful build/replay logs. Original multi-gigabyte search snapshots are hash-bound and remain in their earlier evidence bundles/local archives; they are unnecessary to replay these rank proofs. Replaying the original full-chart oracle translation script additionally uses the wider-MW16 input bundle. No upper ranks or universal novelty.'})
    print('PACKAGED MOD2 ADMISSION AUDIT', len(cases), 'transcripts;', len(rows), 'files;', archive.stat().st_size, 'bytes', flush=True)

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    package(p.parse_args().output.resolve())
