#!/usr/bin/env sage-python
"""Independent Sage polynomial-resultant verification of the prime-support audit."""
import sys
from pathlib import Path
from sage.all import ZZ, PolynomialRing, prime_range
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'elliptic-curves/cas'))
import certify_compact_r17_candidates as cert
import audit_r17_scaling_prime_support as audit

p=cert.read(audit.D/'protocol.json');d=cert.read(audit.OUT)
if p['sources']!=audit.sources() or d['sources']!=audit.sources() or d['protocol_sha256']!=cert.hashed(audit.D/'protocol.json'):raise ArithmeticError('exact source binding differs')
families=cert.read(audit.INPUT)['families'];R=PolynomialRing(ZZ,'t')
if len(families)!=6 or len(d['rows'])!=6:raise ArithmeticError('complete six-family roster required')
for f,row in zip(families,d['rows']):
    if f['family']!=row['family']:raise ArithmeticError('family order differs')
    A=R([ZZ(v) for v in f['A_coefficients_low_to_high']]);B=R([ZZ(v) for v in f['B_coefficients_low_to_high']])
    if A.degree()!=8 or B.degree()!=12:raise ArithmeticError('full degrees required for this independent resultant witness')
    result=A.resultant(B)
    if result!=ZZ(row['homogeneous_resultant']) or not result:raise ArithmeticError('independent polynomial resultant differs')
    remainder=abs(result);factors=[]
    for prime in prime_range(2,998):
        exponent=remainder.valuation(prime)
        if exponent:factors.append([int(prime),int(exponent)]);remainder//=prime**exponent
    if factors!=row['trial_prime_factors'] or str(remainder)!=row['remaining_cofactor'] or row['complete_prime_support']!=(remainder==1):raise ArithmeticError('prime-power stripping differs')
    if row['support_primes_at_least5']!=[q for q,e in factors if q>=5]:raise ArithmeticError('known support differs')
print('INDEPENDENT SIX POLYNOMIAL RESULTANTS AND BOUNDED PRIME SUPPORT PASS',flush=True)
