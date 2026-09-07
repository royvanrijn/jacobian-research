#!/usr/bin/env sage-python
"""Independent five-family polynomial resultants and bounded prime support."""
import sys
from pathlib import Path
from sage.all import ZZ, PolynomialRing, prime_range
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'elliptic-curves/cas'))
import certify_compact_r17_candidates as cert
import audit_mw16_scaling_prime_support as audit

p = cert.read(audit.D/'protocol.json')
d = cert.read(audit.OUT)
if (p['sources'] != audit.sources() or d['sources'] != audit.sources()
        or d['protocol_sha256'] != cert.hashed(audit.D/'protocol.json')):
    raise ArithmeticError('support proof binding differs')
families = cert.read(audit.INPUT)['families']
R = PolynomialRing(ZZ, 't')
if len(families) != 5 or len(d['rows']) != 5:
    raise ArithmeticError('complete five-family roster required')
for f, row in zip(families, d['rows']):
    if f['fibration_id'] != row['family']:
        raise ArithmeticError('family ordering differs')
    A = R([ZZ(v) for v in f['A_coefficients_low_to_high']])
    B = R([ZZ(v) for v in f['B_coefficients_low_to_high']])
    if A.degree() != 8 or B.degree() != 12:
        raise ArithmeticError('full homogeneous degrees required')
    result = A.resultant(B)
    if not result or result != ZZ(row['homogeneous_resultant']):
        raise ArithmeticError('independent polynomial resultant differs')
    remaining = abs(result)
    factors = []
    for q in prime_range(2, p['maximum_trial_prime']+1):
        e = remaining.valuation(q)
        if e:
            factors.append([int(q), int(e)])
            remaining //= q**e
    if (factors != row['trial_prime_factors'] or str(remaining) != row['remaining_cofactor']
            or row['complete_prime_support'] != (remaining == 1)
            or row['score_range_support'] != [q for q, e in factors if q >= 5]):
        raise ArithmeticError('independent bounded prime stripping differs')
print('INDEPENDENT FIVE MW16 RESULTANTS AND PRIME SUPPORT PASS', flush=True)
