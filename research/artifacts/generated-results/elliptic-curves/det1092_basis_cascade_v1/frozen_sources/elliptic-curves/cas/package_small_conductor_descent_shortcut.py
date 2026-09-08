#!/usr/bin/env python3
"""Retain and replay the fixed norm pilot without rerunning slow descent arms."""
import ast
from pathlib import Path
import sys
import zipfile
import prepare_small_conductor_norm_form as forms
from research_runtime.store import checkpoint
from research_runtime.supervisor import Limits, run

ROOT, ART, cert = forms.ROOT, forms.ART, forms.target.original.cert
STEM = 'small_conductor_descent_shortcut'
ZIP = ART / (STEM + '_evidence_v1.zip')
MANIFEST = ART / (STEM + '_evidence_v1.json')
REPORT = ART / (STEM + '_portable_replay_v1.json')
CHECKER = 'elliptic-curves/cas/check_small_conductor_descent_shortcut.sage'


def package():
    if ZIP.exists() or MANIFEST.exists():
        raise FileExistsError('preserve evidence')
    summary = ART / (STEM + '_v1.json')
    if cert.read(summary)['status'] != 'PASS':
        raise ArithmeticError('aggregate proof required')
    files = {summary, ROOT / CHECKER, Path(__file__).resolve(), forms.target.PROOF,
             ROOT / cert.read(forms.target.PROOF)['catalogue']['path'],
             ROOT / 'elliptic-curves/notes/SMALL_CONDUCTOR_DESCENT_SHORTCUT_2026-09-06.md'}
    for directory in [forms.D, forms.D.with_name('small-conductor-norm-pilot-v1')]:
        files.update(p for p in directory.rglob('*') if p.is_file())
    for directory in [ROOT / 'elliptic-curves/ecsearch', ROOT / 'elliptic-curves/cas/research_runtime']:
        files.update(directory.rglob('*.py'))
    # Follow explicit source bindings and local imports to a fixed point.
    pending = list(files)
    seen = set()
    cas = ROOT / 'elliptic-curves/cas'
    while pending:
        path = pending.pop()
        if path in seen:
            continue
        seen.add(path)
        candidates = []
        if path.suffix == '.json':
            data = cert.read(path)
            if isinstance(data, dict):
                for key in ['sources', 'inputs']:
                    bound = data.get(key, {})
                    if isinstance(bound, dict):
                        for name in bound:
                            candidate = ROOT / name
                            if candidate.is_relative_to(ROOT) and candidate.is_file():
                                candidates.append(candidate)
        if path.suffix in ['.py', '.sage']:
            for node in ast.walk(ast.parse(path.read_text())):
                names = []
                base = cas
                if isinstance(node, ast.Import):
                    names = [a.name for a in node.names]
                elif isinstance(node, ast.ImportFrom):
                    if node.level:
                        base = path.parent
                        for _ in range(node.level-1):
                            base = base.parent
                    names = [node.module] if node.module else [a.name for a in node.names]
                for name in names:
                    stem = base.joinpath(*name.split('.'))
                    candidates.extend([stem.with_suffix('.py'), stem / '__init__.py'])
        for candidate in candidates:
            if candidate.is_file() and candidate not in files:
                files.add(candidate)
                pending.append(candidate)
    entries = [{'path': str(p.relative_to(ROOT)), 'sha256': cert.hashed(p),
                'bytes': p.stat().st_size} for p in sorted(files)]
    with zipfile.ZipFile(ZIP, 'x', zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for row in entries:
            archive.write(ROOT / row['path'], row['path'])
    checkpoint(MANIFEST, {'schema': 'elliptic-curves.small-conductor-descent-evidence.v1',
                         'files': entries, 'zip_sha256': cert.hashed(ZIP),
                         'zip_bytes': ZIP.stat().st_size,
                         'requirements': 'Sage with PARI2.17.3, plus hash-pinned /usr/bin/gp for root numbers. Timed class jobs are retained, not rerun.'})
    print('PACKAGED', len(entries), 'files', ZIP.stat().st_size, 'bytes', flush=True)


def replay():
    if REPORT.exists():
        raise FileExistsError('preserve replay report')
    manifest = cert.read(MANIFEST)
    if cert.hashed(ZIP) != manifest['zip_sha256']:
        raise ArithmeticError('archive hash differs')
    destination = forms.D.with_name('small-conductor-descent-isolated-v1')
    destination.mkdir(exist_ok=False)
    with zipfile.ZipFile(ZIP) as archive:
        if any(Path(name).is_absolute() or '..' in Path(name).parts for name in archive.namelist()):
            raise ArithmeticError('unsafe archive member')
        archive.extractall(destination)
    for row in manifest['files']:
        if cert.hashed(destination / row['path']) != row['sha256']:
            raise ArithmeticError('extracted hash differs')
    command = ['/home/royvanrijn/.local/bin/sage', '-python', CHECKER, '--check']
    outcome = run(command, limits=Limits(60, 1610612736), cwd=destination,
                  log_path=destination / 'replay.log',
                  checkpoint_path=destination / 'replay.supervisor.json')
    passed = outcome['outcome'] == 'completed' and outcome['returncode'] == 0
    checkpoint(REPORT, {'schema': 'elliptic-curves.small-conductor-descent-isolated.v1',
                        'status': 'PASS' if passed else 'FAILED_OR_CENSORED',
                        'zip_sha256': manifest['zip_sha256'], 'supervision': outcome})
    print('ISOLATED REPLAY', 'PASS' if passed else 'FAILED_OR_CENSORED', flush=True)
    if not passed:
        raise SystemExit(1)


if __name__ == '__main__':
    {'package': package, 'replay': replay}[sys.argv[1]]()
