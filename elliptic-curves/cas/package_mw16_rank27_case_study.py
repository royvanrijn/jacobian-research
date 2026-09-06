#!/usr/bin/env python3
"""Package and replay only this case study, without historical large archives."""
import argparse
import hashlib
import json
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path

import certify_mw16_rank27_translation_classes as translation
import followup_mw16_rank27_visibility_v2 as trial
import study_mw16_rank27_visibility as visibility
from research_runtime.store import checkpoint

ROOT, ART = trial.ROOT, trial.ART
ARCHIVE = ART / 'mw16_rank27_case_study_evidence_v1.zip'
MANIFEST = ART / 'mw16_rank27_case_study_evidence_v1.json'
REPLAY = ART / 'mw16_rank27_case_study_portable_replay_v1.json'
NOTE = ROOT / 'elliptic-curves/notes/MW16_RANK27_VISIBILITY_AND_TRANSLATION_CLASSES_2026-09-06.md'
TEST = ROOT / 'elliptic-curves/tests/test_mw16_rank27_visibility_study.py'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def package():
    if ARCHIVE.exists() or MANIFEST.exists():
        raise FileExistsError('preserve case study evidence')
    # Exercise the replay paths so lazily imported local modules are included.
    if visibility.expected() != trial.cert.read(visibility.OUTPUT):
        raise ArithmeticError('visibility replay differs')
    if translation.expected() != trial.cert.read(translation.OUT):
        raise ArithmeticError('translation replay differs')
    trial.check()
    files = {Path(__file__).resolve(), NOTE, TEST, visibility.INPUT, visibility.OUTPUT,
             translation.OUT, trial.MOD2, trial.ODD, trial.REPORT,
             ROOT / 'elliptic-curves/cas/followup_mw16_rank27_visibility.py'}
    files.update(ROOT / name for name in trial.sources())
    for version in (1, 2):
        folder = trial.LOCAL / f'mw16-rank27-visibility-followup-v{version}'
        files.update(p for p in folder.rglob('*') if p.is_file())
    for module in list(sys.modules.values()):
        filename = getattr(module, '__file__', None)
        if filename:
            p = Path(filename).resolve()
            if p.is_relative_to(ROOT) and p.suffix == '.py':
                files.add(p)
    rows = []
    with zipfile.ZipFile(ARCHIVE, 'x', compression=zipfile.ZIP_DEFLATED, compresslevel=9) as z:
        for path in sorted(files):
            name = str(path.relative_to(ROOT))
            content = path.read_bytes()
            info = zipfile.ZipInfo(name, date_time=(2026, 9, 6, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o100644 << 16
            z.writestr(info, content)
            rows.append({'path': name, 'bytes': len(content), 'sha256': sha(path)})
    checkpoint(MANIFEST, {'schema': 'elliptic-curves.mw16-rank27-case-study-evidence.v1',
        'archive': str(ARCHIVE.relative_to(ROOT)), 'archive_sha256': sha(ARCHIVE),
        'archive_bytes': ARCHIVE.stat().st_size, 'files': rows,
        'scope': 'Self-contained compact geometry, exact finite-reduction witnesses, '
                 'new twelve-chart trial transcripts and failed preflight evidence. '
                 'The original 43/301 historical full admission transcripts remain in '
                 'their existing archives; this bundle replays their retained geometric '
                 'and point witnesses, not their entire historical worker executions.'})
    print('PACKAGED', len(rows), 'files;', ARCHIVE.stat().st_size, 'bytes', flush=True)


def replay():
    if REPLAY.exists():
        raise FileExistsError('preserve portable replay')
    manifest = trial.cert.read(MANIFEST)
    if sha(ARCHIVE) != manifest['archive_sha256']:
        raise ArithmeticError('archive binding differs')
    commands = [
        ['elliptic-curves/cas/study_mw16_rank27_visibility.py', 'check'],
        ['elliptic-curves/cas/certify_mw16_rank27_translation_classes.py', '--check'],
        ['elliptic-curves/cas/followup_mw16_rank27_visibility_v2.py', 'check'],
        ['-m', 'unittest', 'discover', '-s', 'elliptic-curves/tests',
         '-p', 'test_mw16_rank27_visibility_study.py']]
    result = {'schema': 'elliptic-curves.mw16-rank27-case-study-portable-replay.v1',
              'status': 'RUNNING', 'manifest_sha256': sha(MANIFEST),
              'archive_sha256': sha(ARCHIVE), 'python': sys.version,
              'point_search_performed': False, 'stages': []}
    checkpoint(REPLAY, result)
    with tempfile.TemporaryDirectory(prefix='mw16-rank27-study-') as tmp:
        base = Path(tmp)
        with zipfile.ZipFile(ARCHIVE) as z:
            if set(z.namelist()) != {r['path'] for r in manifest['files']}:
                raise ArithmeticError('archive members differ')
            for row in manifest['files']:
                name = row['path']
                p = base / name
                if not p.resolve().is_relative_to(base):
                    raise ArithmeticError('unsafe archive path')
                content = z.read(name)
                if len(content) != row['bytes'] or hashlib.sha256(content).hexdigest() != row['sha256']:
                    raise ArithmeticError('archive member binding differs')
                p.parent.mkdir(parents=True, exist_ok=True)
                p.write_bytes(content)
        for args in commands:
            completed = subprocess.run([sys.executable, *args], cwd=base, capture_output=True,
                                       text=True, timeout=60)
            result['stages'].append({'command': args, 'returncode': completed.returncode,
                                     'stdout': completed.stdout, 'stderr': completed.stderr})
            if completed.returncode:
                result['status'] = 'FAILED'
                checkpoint(REPLAY, result)
                raise ArithmeticError('isolated replay failed: ' + completed.stderr)
            checkpoint(REPLAY, result)
            print('ISOLATED PASS', args[0], flush=True)
    result['status'] = 'PASS'
    checkpoint(REPLAY, result)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage', choices=['package', 'replay'])
    globals()[parser.parse_args().stage]()
