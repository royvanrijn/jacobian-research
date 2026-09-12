#!/usr/bin/env python3
"""Verify programme relocation and preservation without executing research code."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

from research_programmes import is_active, relocate_reference

ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / 'archive/non-elliptic'


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit() -> None:
    manifest = json.loads((ARCHIVE / 'MANIFEST.json').read_text())
    for record in manifest['files'] + manifest['copies'] + manifest.get('archival_navigation', []):
        path = ROOT / record['path']
        assert path.is_file() and digest(path) == record['sha256'], f'archive source changed: {path}'
    for record in manifest['shared_links']:
        path = ROOT / record['path']
        assert path.is_symlink() and str(path.readlink()) == record['target'] and path.exists(), (
            f'archival replay link changed or broken: {path}'
        )
    for path in manifest['root_moves']:
        assert not (ROOT / path).exists(), f'archived project still occupies active tree: {path}'
    old = json.loads((ARCHIVE / 'MATH_STATUS.json').read_text())
    current = json.loads((ROOT / 'MATH_STATUS.json').read_text())
    by_id = {entry['id']: entry for entry in current['entries']}
    archived = set(manifest['archived_claim_ids'])
    for entry in old['entries']:
        assert entry['id'] in by_id, f"lost historical claim: {entry['id']}"
        if entry['id'] not in archived:
            assert is_active(by_id[entry['id']]), f"EC/K3 claim accidentally archived: {entry['id']}"
            continue
        expected = {**entry, 'programme_status': 'archived'}
        for field in ['canonical_source', 'checker']:
            if expected[field] is not None:
                expected[field] = relocate_reference(expected[field])
        for field in ['software_lock', 'consumers']:
            expected[field] = [relocate_reference(p) for p in expected[field]]
        assert by_id[entry['id']] == expected, f"archiving changed mathematical metadata: {entry['id']}"
    for record in manifest.get('link_relocations', []):
        snapshot = ROOT / record['before_snapshot']
        assert digest(snapshot) == record['before_sha256'], f'link relocation snapshot changed: {snapshot}'
    old_legacy = json.loads((ARCHIVE / 'knowledge/legacy_work_review.json').read_text())
    legacy = json.loads((ROOT / 'knowledge/legacy_work_review.json').read_text())
    assert [(i['id'], i['text']) for i in legacy['items']] == [(i['id'], i['text']) for i in old_legacy['items']], (
        'archiving lost or changed an inherited obligation'
    )
    old_items = {item['id']: item for item in old_legacy['items']}
    for item in legacy['items']:
        if not is_active(item):
            assert ('resolution' in item) == ('resolution' in old_items[item['id']]), (
                'archiving must not supply a completion record'
            )
    print(f"PASS programme archive: {len(manifest['files'])} byte-identical moved files, "
          f"{len(archived)} claims preserved, shared replay inputs and all inherited obligations retained")


if __name__ == '__main__':
    audit()
