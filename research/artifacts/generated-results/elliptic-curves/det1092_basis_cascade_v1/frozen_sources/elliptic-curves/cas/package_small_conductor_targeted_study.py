#!/usr/bin/env python3
"""Package and independently replay the fixed small-conductor study."""
import ast
import subprocess
import sys
import zipfile
from pathlib import Path

import study_small_conductor_target as study
from research_runtime.store import checkpoint
from research_runtime.supervisor import Limits, run

ROOT, D, ART = study.ROOT, study.D, study.ART
ZIP = ART / 'small_conductor_targeted_study_evidence_v3.zip'
MANIFEST = ART / 'small_conductor_targeted_study_evidence_v3.json'
REPORT = ART / 'small_conductor_targeted_study_portable_replay_v3.json'


def package():
    if ZIP.exists() or MANIFEST.exists():
        raise FileExistsError('preserve evidence')
    protocol = study.protocol()
    if study.cert.read(D / 'ledger.json')['status'] != 'PASS':
        raise ArithmeticError('terminal replay required')
    files = {ROOT / p for p in (*protocol['sources'], *protocol['inputs'])}
    files.update(p for p in D.rglob('*') if p.is_file())
    failed = D.with_name('small-conductor-targeted-study-v1')
    files.update(p for p in failed.rglob('*') if p.is_file())
    files.update((ROOT / 'elliptic-curves/ecsearch').rglob('*.py'))
    files.update((ROOT / 'elliptic-curves/cas/research_runtime').rglob('*.py'))
    for module in list(sys.modules.values()):
        name = getattr(module, '__file__', None)
        if name:
            path = Path(name).resolve()
            if path.is_relative_to(ROOT) and path.is_file():
                files.add(path)
    files.update([study.OUT, Path(__file__).resolve(),
                  ROOT / 'elliptic-curves/cas/audit_small_conductor_targeted_cloud.py',
                  ART / 'small_conductor_targeted_cloud_modl_v2.json',
                  ROOT / study.cert.read(study.PROOF)['catalogue']['path']])
    # Include static local imports recursively, plus all explicitly pinned source
    # files; the extracted replay detects any missing runtime dependency.
    pending = [p for p in files if p.suffix == '.py']
    seen = set()
    cas = ROOT / 'elliptic-curves/cas'
    while pending:
        p = pending.pop()
        if p in seen:
            continue
        seen.add(p)
        tree = ast.parse(p.read_text())
        for node in ast.walk(tree):
            names = []
            base = cas
            if isinstance(node, ast.Import):
                names = [a.name for a in node.names]
            elif isinstance(node, ast.ImportFrom):
                if node.level:
                    base = p.parent
                    for _ in range(node.level - 1):
                        base = base.parent
                if node.module:
                    names = [node.module]
                else:
                    names = [a.name for a in node.names]
            for name in names:
                stem = base.joinpath(*name.split('.'))
                for candidate in [stem.with_suffix('.py'), stem / '__init__.py']:
                    if candidate.is_file() and candidate not in files:
                        files.add(candidate)
                        pending.append(candidate)
                for parent in stem.parents:
                    if parent == cas.parent:
                        break
                    init = parent / '__init__.py'
                    if init.is_file() and init not in files:
                        files.add(init)
                        pending.append(init)
    entries = [{'path': str(p.relative_to(ROOT)), 'bytes': p.stat().st_size,
                'sha256': study.cert.hashed(p)} for p in sorted(files)]
    with zipfile.ZipFile(ZIP, 'x', zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for row in entries:
            archive.write(ROOT / row['path'], row['path'])
    checkpoint(MANIFEST, {'schema': 'elliptic-curves.small-conductor-study-evidence.v3',
                         'zip_sha256': study.cert.hashed(ZIP), 'zip_bytes': ZIP.stat().st_size,
                         'python_version': sys.version,
                         'pari_version': subprocess.run(['/usr/bin/gp', '--version-short'], capture_output=True, text=True, check=True, timeout=5).stdout.strip(),
                         'superseded_packaging_attempt': 'V2 omitted the ecsearch import root and failed before replay. Its archive and failed report are retained. V3 adds import dependencies without changing study inputs or results.',
                         'files': entries, 'requirements': 'Python3 and the hash-pinned /usr/bin/gp; no Sage, factorization or point search needed for replay.'})
    print('PACKAGED', len(entries), 'files', ZIP.stat().st_size, 'bytes', flush=True)


def replay():
    if REPORT.exists():
        raise FileExistsError('preserve isolated replay report')
    manifest = study.cert.read(MANIFEST)
    if study.cert.hashed(ZIP) != manifest['zip_sha256']:
        raise ArithmeticError('archive hash differs')
    destination = D.with_name('small-conductor-targeted-isolated-v3')
    if destination.exists():
        raise FileExistsError('preserve extracted replay')
    destination.mkdir()
    with zipfile.ZipFile(ZIP) as archive:
        for name in archive.namelist():
            if Path(name).is_absolute() or '..' in Path(name).parts:
                raise ArithmeticError('unsafe archive member')
        archive.extractall(destination)
    for row in manifest['files']:
        if study.cert.hashed(destination / row['path']) != row['sha256']:
            raise ArithmeticError('extracted hash differs')
    result = {'schema': 'elliptic-curves.small-conductor-study-isolated.v3',
              'status': 'RUNNING', 'zip_sha256': manifest['zip_sha256'], 'stages': []}
    for stage in ['replay', 'check', 'odd-cloud']:
        command = [sys.executable, 'elliptic-curves/cas/study_small_conductor_target.py', stage]
        if stage == 'odd-cloud':
            command = [sys.executable, 'elliptic-curves/cas/audit_small_conductor_targeted_cloud.py', '--check']
        outcome = run(command,
                      limits=Limits(300, 1610612736), cwd=destination,
                      log_path=destination / (stage + '.log'),
                      checkpoint_path=destination / (stage + '.supervisor.json'))
        result['stages'].append({'stage': stage, 'supervision': outcome})
        if outcome['outcome'] != 'completed' or outcome['returncode'] != 0:
            result['status'] = 'FAILED_OR_CENSORED'
            checkpoint(REPORT, result)
            return
    result['status'] = 'PASS'
    checkpoint(REPORT, result)
    print('ISOLATED REPLAY PASS: exact chart histories and arithmetic summary', flush=True)


if __name__ == '__main__':
    {'package': package, 'replay': replay}[sys.argv[1]]()
