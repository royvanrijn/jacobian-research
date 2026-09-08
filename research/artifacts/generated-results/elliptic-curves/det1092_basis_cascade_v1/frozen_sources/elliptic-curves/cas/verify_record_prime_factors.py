#!/usr/bin/env python3
"""Rigorous, bounded primality supplement for the pinned ICARM 302/356 replays.

The historical producers used SymPy isprime, which above 2^64 is a probable
prime test. Preserve those source-bound artifacts and supplement their large
factor claims with PARI isprime (not ispseudoprime). No factor search is run.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from math import isfinite
from pathlib import Path
import subprocess

import verify_icarm_curve302_rank31 as curve302
import verify_icarm_curve356_rank29 as curve356


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elliptic-curves/record_prime_factor_proofs_20260904.json"


def prove_primes(values, *, timeout: float = 15.0):
    """Prove every supplied factor prime, or raise on failure/timeout."""
    if not isfinite(timeout) or timeout <= 0:
        raise ValueError("a finite positive wall limit is required")
    values = sorted(set(values))
    if not values or any(type(value) is not int or value < 2 for value in values):
        raise ValueError("supply integer factors at least two")
    program = 'print(version());\n' + ''.join(
        f'print("{value}|",isprime({value}));\n' for value in values
    ) + 'quit\n'
    result = subprocess.run(
        ['gp', '-q', '-f'], input=program, capture_output=True, text=True,
        timeout=timeout, check=True,
    )
    lines = result.stdout.strip().splitlines()
    if '***' in result.stderr or len(lines) != len(values) + 1:
        raise ArithmeticError("PARI primality proof failed or returned incomplete output")
    if lines[1:] != [f'{value}|1' for value in values]:
        raise ArithmeticError("a supplied factor is not proved prime")
    return lines[0], [str(value) for value in values]


def build_certificate(*, timeout: float = 15.0):
    records = {}
    factors = set()
    for label, module in [('302', curve302), ('356', curve356)]:
        tables = {name: getattr(module, name) for name in (
            'C4_FACTORIZATION', 'DISCRIMINANT_FACTORIZATION', 'CONDUCTOR_FACTORIZATION'
        )}
        c4 = module.weierstrass_invariants()[4]
        expected = [c4, abs(module.PUBLIC_DISCRIMINANT), module.PUBLIC_CONDUCTOR]
        for table, value in zip(tables.values(), expected):
            if module.factor_product(table) != value:
                raise ArithmeticError(f'curve {label}: factor product mismatch')
            factors.update(prime for prime, _ in table)
        source = Path(module.__file__).resolve()
        records[label] = {
            'source': str(source.relative_to(ROOT)),
            'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
            'factorizations': {name: [[str(p), e] for p, e in table]
                               for name, table in tables.items()},
        }
    version, proved = prove_primes(factors, timeout=timeout)
    return {
        'schema': 'elliptic-curves.record-prime-factor-proofs.v1',
        'method': 'PARI isprime: unconditional primality test',
        'pari_version': version,
        'records': records,
        'proved_primes': proved,
        'claim_boundary': 'Proves the listed factors prime and their products; the separate record checkers establish point independence and local reduction.',
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--timeout', type=float, default=15.0)
    parser.add_argument('--output', type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    certificate = build_certificate(timeout=args.timeout)
    if args.check:
        pinned = json.loads(args.output.read_text())
        # Software version is provenance; all factors, proofs, products and
        # source hashes are still replayed and compared exactly.
        pinned.pop('pari_version')
        current = dict(certificate)
        current.pop('pari_version')
        if current != pinned:
            raise SystemExit('record prime-factor certificate differs')
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(certificate, indent=2, sort_keys=True)+'\n')
    print(f"PASS rigorous record-factor proofs: {len(certificate['proved_primes'])} primes")


if __name__ == '__main__':
    main()
