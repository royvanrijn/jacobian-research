#!/usr/bin/env python3
"""Replay frozen626 controls independently of the moving current catalogue."""
from pathlib import Path
import prepare_recent_icarm_controls as controls
import refresh_icarm_local_database as database

ROOT = Path(__file__).resolve().parents[2]
PINNED = ROOT/'elliptic-curves/data/icarm_snapshots/626_manifest.json'

if __name__ == '__main__':
    controls.CURRENT = PINNED
    database.CURRENT = PINNED
    controls.run(check=True)
