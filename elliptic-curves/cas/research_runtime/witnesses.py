"""Cheap retained-witness replay, separate from optional full regeneration."""

from copy import deepcopy
from hashlib import sha256
from pathlib import Path
import os


SOURCE_ARCHIVE="archive/elliptic-curves/runtime-source-snapshots-2026-09-05"


def retained_source(root,name,expected):
    """Resolve generation-time source bytes without pretending they are current.

    Source snapshots are provenance only. The active replayer must independently
    check the mathematical witness; retrieving old code does not verify a claim.
    """
    root=Path(root).resolve();path=(root/name).resolve()
    if not path.is_relative_to(root):raise ValueError("source path outside repository")
    if path.exists() and sha256(path.read_bytes()).hexdigest()==expected:return path
    snapshot=root/SOURCE_ARCHIVE/expected/Path(name).name
    if not snapshot.exists() and os.environ.get('EC_RUNTIME_SOURCE_ARCHIVE'):
        if len(expected)!=64 or any(c not in '0123456789abcdef' for c in expected):raise ValueError('invalid source digest')
        snapshot=Path(os.environ['EC_RUNTIME_SOURCE_ARCHIVE'])/expected/Path(name).name
    if not snapshot.exists() or sha256(snapshot.read_bytes()).hexdigest()!=expected:
        raise ValueError(f"missing generation-time source witness: {name}")
    return snapshot


def compare_replay(stored,replayed,*,root,source_paths=()):
    """Compare exact outputs while checking both generations' source provenance.

    Only explicitly named implementation-source hashes may differ. All data,
    mathematical outputs, algorithms' recorded options and non-source inputs
    must match exactly. This function does not regenerate discoveries or mutate
    the stored certificate. Full regeneration may publish a new provenance set.
    """
    before,after=deepcopy(stored),deepcopy(replayed)
    for group in ("inputs","source_hashes"):
        if group not in before or group not in after:continue
        for name in source_paths:
            old=before[group].get(name);new=after[group].get(name)
            if old==new:continue
            if old is None or new is None:raise ValueError("source provenance shape changed")
            retained_source(root,name,old)
            active=Path(root)/name
            if not active.exists() or sha256(active.read_bytes()).hexdigest()!=new:
                raise ValueError("active replayer source hash is incorrect")
            # Normalize only the comparison copies. The returned/written
            # provenance remains the actual generation's provenance.
            before[group][name]=after[group][name]="VERIFIED_GENERATION_SOURCE"
    if before!=after:raise ArithmeticError("retained mathematical witness does not replay exactly")
