#!/usr/bin/env sage-python
"""Two fixed maps versus the existing orbit8044 conic, no parameter search.

Squarefree degree2d pullback yields genus d-1 and excludes rational lifting
of the new P1 parameter map through the old conic. Modular certificates are
replayed algebraically by Bezout. The fixed12-prime limit is not extended.
"""
import hashlib,json
from pathlib import Path
from sage.all import QQ,GF,PolynomialRing
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_pencil_multiples_v2'
CONIC=ART/'det1092_orbit8044_rank18_base_change_v2.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
protocol=json.loads((OUT/'protocol.json').read_text())
q=list(map(QQ,json.loads(CONIC.read_text())['curve_over_Q']['q_coefficients']))
R=PolynomialRing(QQ,'z');rows=[]
for n in [2,3]:
    path=OUT/('multiple%d.json'%n);data=json.loads(path.read_text());d=data['parameter_map_degree']
    N,D=R(data['t_of_z']['numerator']),R(data['t_of_z']['denominator'])
    witness=None;attempts=[]
    for p in protocol['ramification_certificate_primes']:
        try:
            Rp=PolynomialRing(GF(p),'z');np,dp=Rp(N),Rp(D)
            qp=list(map(GF(p),q))
        except (ZeroDivisionError,ValueError):
            attempts.append({'prime':p,'status':'denominator_bad'});continue
        f=qp[0]*dp*dp+qp[1]*np*dp+qp[2]*np*np
        if f.degree()!=2*d:
            attempts.append({'prime':p,'status':'degree_drop'});continue
        g,a,b=f.xgcd(f.derivative()*dp)
        assert a*f+b*f.derivative()*dp==g
        attempts.append({'prime':p,'gcd_degree':int(g.degree())})
        if g==1:
            witness={'prime':p,'polynomial_mod_p':list(map(int,f.list())),
                     'degree':2*d,'genus':d-1,'squarefree_and_coprime_to_denominator':True}
            break
    rows.append({'n':n,'status':'NO_RATIONAL_MAP_LIFT_TO_ORBIT8044' if witness else 'UNKNOWN',
                 'certificate':witness,'attempts':attempts,'map_sha256':sha(path)})
result={'classification':'new verified obstruction',
        'cases':rows,'conic_sha256':sha(CONIC),'checker_sha256':sha(Path(__file__)),
        'boundary':'No rational-function lift through the old conic. No rank19 construction or assertion about individual rational intersections.'}
path=OUT/'conic-lift-obstruction.json'
if path.exists():assert json.loads(path.read_text())==result
else:
    with path.open('x') as stream:json.dump(result,stream,indent=2,sort_keys=True);stream.write('\n')
print(json.dumps({'classification':result['classification'],
    'cases':[{'n':r['n'],'status':r['status'],'genus':None if r['certificate'] is None else r['certificate']['genus'],
              'prime':None if r['certificate'] is None else r['certificate']['prime']} for r in rows]},sort_keys=True))
