"""Cheap projective trace lookup for the certified determinant1092 parent.

Input parameters are supplied by the caller; this module generates no search
population. Singular table entries remain missing, never rank exclusions.
The score is a finite-prime scheduling heuristic, not a rank estimate.
"""
import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import log
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DEFAULT = ROOT/'artifacts/generated-results/elliptic-curves/curve302_parent_trace_tables_v1.json'


def load_tables(path=DEFAULT):
    data = json.loads(Path(path).read_text())
    if data['schema'] != 'curve302.projective-traces.v1':
        raise ValueError('unknown trace table schema')
    parent = ROOT/data['parent_path']
    if sha256(parent.read_bytes()).hexdigest() != data['parent_sha256']:
        raise ValueError('parent changed: rebuild and reverify trace tables')
    for row in data['tables']:
        if len(row['traces']) != row['prime']+1:
            raise ValueError('incomplete projective table')
    return data


def traces_at(data, parameter):
    infinity = parameter == 'infinity'
    q = None if infinity else Fraction(parameter)
    result = []
    for row in data['tables']:
        p = row['prime']
        denominator = 0 if infinity else q.denominator % p
        residue = p if not denominator else q.numerator*pow(denominator,-1,p) % p
        result.append(row['traces'][residue])
    return result


def score(data, parameter):
    traces = traces_at(data, parameter)
    bands = {}
    for band in ['discovery','held_out']:
        entries = [(row['prime'],a) for row,a in zip(data['tables'],traces) if row['band']==band]
        bands[band] = {'score':sum(-a*log(p)/p for p,a in entries if a is not None),
                       'used_primes':sum(a is not None for p,a in entries),
                       'missing_primes':[p for p,a in entries if a is None]}
    return {'parameter':str(parameter),'bands':bands,'rank_bound':None}


def legacy_text(data):
    """Export the same trace weights to the existing C++ rational scanner.

    Integer score units are 10^-12 times -a_p*log(p)/p. A missing trace has
    good=0; the placeholder zero is not a mathematical trace claim.
    """
    lines=['RATIONAL_NAGAO_LOCAL_TABLE_V1','F CURVE302_DET1092_MW17 8 12']
    for band,label in [('discovery','D'),('held_out','H')]:
        selected=[row for row in data['tables'] if row['band']==band]
        lines.append(f'B {label} {len(selected)}')
        for row in selected:
            p=row['prime'];lines.append(f'P {p}')
            for a in row['traces']:
                lines.append('0 0 0' if a is None else f'1 {a} {round(-a*log(p)/p*10**12)}')
    return '\n'.join(lines+['END',''])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tables',type=Path,default=DEFAULT)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--parameters',type=Path,help='one rational parameter or infinity per line')
    group.add_argument('--export-cpp',type=Path,help='export tables for the existing C++ rational scanner')
    args = parser.parse_args()
    data = load_tables(args.tables)
    if args.export_cpp:
        if args.export_cpp.exists():
            raise FileExistsError('use a fresh export path')
        args.export_cpp.write_text(legacy_text(data))
        raise SystemExit(0)
    for line in args.parameters.read_text().splitlines():
        if line.strip() and not line.lstrip().startswith('#'):
            print(json.dumps(score(data,line.strip()),sort_keys=True))
