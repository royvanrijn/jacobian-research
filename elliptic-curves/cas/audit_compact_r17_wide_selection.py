#!/usr/bin/env python3
"""Recompute all24 selected full-prime scores and independent cardinality samples."""
import argparse
from fractions import Fraction as F
from math import log
from pathlib import Path
import certify_compact_r17_candidates as cert
import compact_atlas_specialization as spec
from mod2_reduction_independence import _is_prime
from research_runtime.store import digest
ROOT = Path(__file__).resolve().parents[2]
PARENT = ROOT/'artifacts/local/elliptic-curves/compact-six-r17-h1024-v2'
DIRECTORY = ROOT/'artifacts/local/elliptic-curves/compact-six-r17-h4096-v1'

def audit(output):
    if output.exists():
        raise FileExistsError('preserve selection audit')
    protocol = cert.read(DIRECTORY/'protocol.json')
    primes = [p for p in range(5, 4094) if _is_prime(p)]
    if len(primes) != 562:
        raise ArithmeticError('prime universe differs')
    phi = list(range(4097))
    for p in range(2, 4097):
        if phi[p] == p:
            for j in range(p, 4097, p):
                phi[j] -= phi[j]//p
    count = 4*sum(phi[1:])-2
    if count != 20400078:
        raise ArithmeticError('independent primitive address count differs')
    rows, tables = [], {}
    for family in cert.read(spec.ATLAS)['families']:
        label = family['family']
        population_path = DIRECTORY/label/'population.json'
        population = cert.read(population_path)
        if population['protocol_hash'] != digest(protocol) or population['candidate_count'] != count:
            raise ArithmeticError('population binding differs')
        model = {key: family[key]+['0']*(n-len(family[key])) for key, n in (
            ('A_coefficients_low_to_high', 9), ('B_coefficients_low_to_high', 13))}
        candidates = population['finalists']
        totals = [0]*4; good_counts = [0]*4; checks = []
        for p in primes:
            path = PARENT/label/'trace-tables'/f'{p}.json'
            table = cert.read(path)
            if table['input'] != {'family': label, 'model_hash': digest(model), 'prime': p}:
                raise ArithmeticError('trace model binding differs')
            tables[str(path.relative_to(ROOT))] = cert.hashed(path)
            for i, row in enumerate(candidates):
                n, d = row['numerator'], row['denominator']
                t = n*pow(d, -1, p) % p if d % p else p
                ap, good = table['traces'][t], table['good'][t]
                if good:
                    totals[i] += round((2-ap)/(p+1-ap)*log(p)*10**12)
                    good_counts[i] += 1
                if p in (5, 7, 11, 13, 997, 4093):
                    values = []
                    for key in model:
                        cs = [F(q) for q in model[key]]
                        residues = [q.numerator*pow(q.denominator, -1, p) % p for q in cs]
                        value = residues[-1] if t == p else sum(c*pow(t, j, p) for j, c in enumerate(residues)) % p
                        values.append(value)
                    a, b = values
                    actual_good = (4*a**3+27*b*b) % p != 0
                    # Independent scalar Euler-criterion character sum.
                    trace = 0
                    for x in range(p):
                        value = (x*x*x+a*x+b) % p
                        if value:
                            trace -= 1 if pow(value, (p-1)//2, p) == 1 else -1
                    if actual_good != good or trace != ap:
                        raise ArithmeticError('independent selected-fibre cardinality differs')
                    checks.append({'parameter': row['parameter'], 'prime': p, 'trace': trace, 'good': good})
        if any(total != row['score_units'] or good != row['good_primes'] for total, good, row in zip(totals, good_counts, candidates)):
            raise ArithmeticError('selected full-prime score differs')
        rows.append({'family': label, 'population_sha256': cert.hashed(population_path),
                     'finalists': candidates, 'independent_cardinality_checks': checks})
    cert.write(output, {'schema': 'elliptic-curves.wide-r17-selection-audit.v1', 'status': 'PASS',
        'sources': {str(p.relative_to(ROOT)): cert.hashed(p) for p in (Path(__file__).resolve(), spec.ATLAS, DIRECTORY/'protocol.json')},
        'trace_table_hashes': tables, 'families': rows, 'primes_per_score': 562,
        'primitive_parameters_per_family': count, 'total_parameters': 6*count,
        'claim_boundary': 'All24 chosen scores are recomputed from all declared trace residues;144 sampled selected-fibre traces are independently counted, including primes997 and4093. The full population count is independently recomputed. This does not independently rerun every scanner score or certify rank.'})
    print('AUDITED WIDE R17 SELECTION:24 full scores,144 scalar cardinalities,122400468 primitive addresses', flush=True)

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    audit(p.parse_args().output.resolve())
