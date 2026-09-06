#!/usr/bin/env python3
"""Bounded portable verification of a principal-relation ledger for Cl(K)/2."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil
import sys
import time
from research_runtime.store import checkpoint
from research_runtime.supervisor import Limits, run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--check', action='store_true', help='replay an existing certificate')
    parser.add_argument('--wall-seconds', type=int, default=120)
    parser.add_argument('--rss-mib', type=int, default=1536)
    parser.add_argument('--run-dir', type=Path)
    parser.add_argument('--worker', action='store_true', help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.wall_seconds <= 0 or args.rss_mib <= 0:
        parser.error('positive resource limits required')
    if not args.worker:
        sage = shutil.which('sage')
        if not sage:
            raise RuntimeError('Sage is required')
        directory = args.run_dir or Path('artifacts/local/class-span-runs')/str(time.time_ns())
        directory.mkdir(parents=True, exist_ok=False)
        command = [sage, '-python', str(Path(__file__).resolve()), str(args.input.resolve()),
                   '--output', str(args.output.resolve()), '--worker']
        if args.check:
            command.append('--check')
        result = run(command, limits=Limits(args.wall_seconds, args.rss_mib*1024**2),
                     cwd=Path.cwd(), log_path=directory/'replay.log',
                     checkpoint_path=directory/'supervisor.json')
        print((directory/'replay.log').read_text(), end='')
        if result['outcome'] != 'completed' or result['returncode'] != 0:
            raise SystemExit('No new certificate: '+result['outcome'])
        return
    import class_span_grh as engine
    from sage.all import pari
    from sage.version import version
    document = json.loads(args.input.read_text())
    certificate = engine.verify_document(document)
    certificate['software'] = {'sage':version, 'pari':str(pari('version()')),
        'engine_sha256':sha256(Path(engine.__file__).read_bytes()).hexdigest(),
        'checker_sha256':sha256(Path(__file__).read_bytes()).hexdigest()}
    if args.check:
        if json.loads(args.output.read_text()) != certificate:
            raise ArithmeticError('class-span replay differs')
    else:
        if args.output.exists():
            raise FileExistsError('preserve existing certificate')
        checkpoint(args.output, certificate)
    print(certificate['status'], 'class-2-rank upper bound:',
          certificate['class_two_rank_upper_bound_under_grh'])


if __name__ == '__main__':
    main()
