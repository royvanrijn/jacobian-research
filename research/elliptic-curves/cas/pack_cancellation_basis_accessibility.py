#!/usr/bin/env python3
"""Close the audit bundle's import dependency and check an empty-root preflight.

The first byte-checked archive remains unchanged. It omitted the supplemental
finish module imported by the audit entry point; no arithmetic result changes.
"""
import json
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile

from finite_cancellation_corpus import ROOT, digest
from cancellation_scheduler_fresh import source_closure
from cancellation_basis_accessibility import read, sha, need, new_write

OUT = ROOT/'artifacts/generated-results/elliptic-curves/cancellation_basis_accessibility_v3'


def main():
    old = read(OUT/'replay-manifest.json')
    paths = {ROOT/r['path'] for r in old['entries']}
    source = Path(__file__).resolve(); paths.add(source)
    closure = source_closure([ROOT/'elliptic-curves/cas'/name for name in (
        'cancellation_basis_accessibility_v3.py', 'verify_cancellation_basis_accessibility_v3.py',
        'retain_cancellation_basis_accessibility.py')])
    additions = sorted(set(closure)-{r['path'] for r in old['entries']})
    paths.update(ROOT/name for name in additions)
    paths.update([OUT/'replay-manifest.json', OUT/'completion.json'])
    entries = [{'path': str(p.relative_to(ROOT)), 'sha256': sha(p), 'bytes': p.stat().st_size} for p in sorted(paths)]
    bundle = OUT/'portable-replay.tar.gz'; need(not bundle.exists(), 'preserve previous portable bundle')
    with tarfile.open(bundle, 'w:gz') as archive:
        for item in entries:
            archive.add(ROOT/item['path'], arcname=item['path'], recursive=False)
    with tempfile.TemporaryDirectory(prefix='basis-accessibility-empty-root-') as temp:
        root = Path(temp)
        with tarfile.open(bundle, 'r:gz') as archive:
            for member in archive.getmembers():
                need(member.isfile() and not Path(member.name).is_absolute() and '..' not in Path(member.name).parts,
                     'non-file or escaping member')
                data = archive.extractfile(member).read(); dest = root/member.name; dest.parent.mkdir(parents=True, exist_ok=True)
                dest.write_bytes(data)
        for item in entries:
            need(sha(root/item['path']) == item['sha256'], 'portable bytes differ')
        code = "import sys;sys.path.insert(0,'elliptic-curves/cas');import cancellation_basis_accessibility_v3 as a;from verify_cancellation_basis_accessibility_v3 import verify;p=a.base.guard();[a.loaded_inputs(r) for r in p['cases']];print('PASS_EMPTY_ROOT_IMPORT_AND_ALL_INPUT_BINDINGS')"
        check = subprocess.run([sys.executable, '-c', code], cwd=root, text=True, capture_output=True)
        need(check.returncode == 0, 'empty-root preflight failed: '+check.stderr)
    result = {'status': 'PASS_PORTABLE_BYTE_IMPORT_AND_INPUT_PREFLIGHT', 'entries': entries, 'files': len(entries),
        'bundle_sha256': sha(bundle), 'previous_manifest_sha256': sha(OUT/'replay-manifest.json'),
        'added_import_dependencies': additions, 'source_sha256': sha(source), 'preflight': check.stdout.strip(),
        'point_search_calls': 0,
        'boundary': 'All archive bytes, standalone imports, frozen source hashes and eighteen case input bindings checked in an empty root. This supplements the original worker arithmetic replays; it is not another full arithmetic replay.'}
    new_write(OUT/'portable-manifest.json', result)
    print(json.dumps({k: v for k, v in result.items() if k != 'entries'}), flush=True)


if __name__ == '__main__': main()
