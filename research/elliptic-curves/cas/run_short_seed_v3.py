#!/usr/bin/env python3
"""Frozen 100-invocation adaptive V3 first pass on one prepared native seed.

Uses the existing parent runner unchanged. The total budget spans all epochs;
a certified gain triggers an immediate rebuild while budget remains. Terminal
receipts preserve exposure for a separately authorized deeper continuation.
"""
import argparse
import fcntl
from pathlib import Path

import run_parent_seed_v3 as parent

original_freeze = parent.freeze


def freeze(folder):
    protocol = original_freeze(folder)
    protocol.update(schema='prepared-parent-short-pass.v1', max_charts=100,
                    search_wall_limit_seconds=1800,
                    replay_wall_limit_seconds=1800)
    protocol['sources'][str(Path(__file__).resolve().relative_to(parent.ROOT))] = parent.sha(Path(__file__))
    protocol['scope'] += (' First-pass scheduling experiment: at most100 actual '
                          'invocations in total, retaining the original centre/map order. '
                          'Budget exhaustion leaves the untested suffix UNKNOWN; '
                          'it is not policy exhaustion or an exact rank.')
    parent.checkpoint(folder/'protocol.json', protocol)
    return protocol


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('search', 'replay'))
    parser.add_argument('--folder', required=True, type=Path)
    parser.add_argument('--preparation', type=Path)
    args = parser.parse_args()
    parent.PREP = args.preparation.resolve() if args.preparation else None
    parent.BANK = parent.PREP
    parent.freeze = freeze
    parent.require(args.folder.exists() or parent.PREP is not None, 'new pass needs preparation')
    if args.folder.exists():
        protocol = parent.read(args.folder/'protocol.json')
        parent.require(protocol['schema'] == 'prepared-parent-short-pass.v1' and
                       protocol['max_charts'] == 100, 'not a frozen short pass')
    with (args.folder.parent/(args.folder.name+'.lock')).open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        (parent.search if args.mode == 'search' else parent.replay)(args.folder)
