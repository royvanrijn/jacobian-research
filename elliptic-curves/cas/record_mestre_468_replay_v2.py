#!/usr/bin/env python3
"""Bind the completed fresh-directory replay to its exact portable inputs."""
import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / 'artifacts/generated-results/elliptic-curves'
LOCAL = ROOT / 'artifacts/local/elliptic-curves/mestre-468-standalone-v2'
OUT = ART / 'mestre_468_portable_replay_v2.json'


def compute():
    supervisor = json.loads((LOCAL / 'supervisor.json').read_text())
    transcript = (LOCAL / 'replay.log').read_bytes()
    assert supervisor['outcome'] == 'completed' and supervisor['returncode'] == 0
    assert supervisor['failure_reason'] is None
    assert hashlib.sha256(transcript).hexdigest() == supervisor['log_sha256']
    assert transcript.decode().splitlines() == [
        f'PASS u{u} 714-height-chain / counts / rank / saturation / discriminant witnesses'
        for u in (11, 13, 17, 19, 23, 29)
    ] + ['PASS6 INDEPENDENT DET468 PARENT REPLAYS']
    files = {}
    for source in (ROOT / 'elliptic-curves/cas/verify_mestre_468_bundle_v2.sage',
                   ART / 'mestre_468_replay_bundle_v1.json'):
        content = source.read_bytes()
        assert content == (LOCAL / source.name).read_bytes()
        files[str(source.relative_to(ROOT))] = {
            'sha256': hashlib.sha256(content).hexdigest(), 'bytes': len(content)}
    return {
        'schema': 'elliptic-curves.mestre-468-portable-replay.v2',
        'status': 'PASS', 'parents': [11, 13, 17, 19, 23, 29],
        'portable_files': files,
        'command': 'sage -python verify_mestre_468_bundle_v2.sage --input mestre_468_replay_bundle_v1.json',
        'supervisor': supervisor, 'transcript': transcript.decode(),
        'checks': '78 direct rational generic halvings, 714 heights total, exact generic relations and base involutions, finite saturation witnesses, every saved finite fibre count, Picard bounds and determinant arithmetic',
        'scope': 'Completed independent arithmetic replay using the mathematical implications in MESTRE_DETERMINANT468_PARENTS_2026-09-07.md. No formal proof-assistant verification or new point search.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--check', action='store_true')
    result = compute()
    if parser.parse_args().check:
        assert result == json.loads(OUT.read_text())
    else:
        if OUT.exists():
            raise FileExistsError('preserve completed replay record')
        OUT.write_text(json.dumps(result, indent=2) + '\n')
    print('PASS6 portable files match the completed independent replay')
