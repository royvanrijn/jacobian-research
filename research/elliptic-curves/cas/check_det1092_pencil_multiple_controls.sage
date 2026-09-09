#!/usr/bin/env sage-python
"""Root-free modular certificates for the unchanged nine-address panel.

Only evaluate the frozen degree20 and58 maps. No point search, root hunt over
Q, new address, or construction retuning. All noncertified cases stay UNKNOWN.
The default replays; --write preserves the first deterministic report.
"""
import argparse
import hashlib
import json
from pathlib import Path
from sage.all import QQ, ZZ, GF, PolynomialRing, gcd, lcm

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT/'artifacts/generated-results/elliptic-curves'
OUT = ART/'det1092_pencil_multiples_v2'
ROSTER = ART/'det1092_rr_generic_point_controls_v2/protocol.json'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def report():
    protocol = json.loads((OUT/'protocol.json').read_text())
    primes = protocol['ramification_certificate_primes']
    cases = json.loads(ROSTER.read_text())['cases']
    R = PolynomialRing(QQ, 'z')
    data = []
    for n in protocol['multiples']:
        path = OUT/('multiple%d.json'%n)
        row = json.loads(path.read_text())
        tmap = row['t_of_z']
        nn,dd = R(tmap['numerator']), R(tmap['denominator'])
        degree = row['parameter_map_degree']
        for case in cases:
            tau = QQ(case['parameter'])
            polynomial = nn-tau*dd
            integral = polynomial*lcm(c.denominator() for c in polynomial.list())
            integral /= gcd([ZZ(c) for c in integral.list()])
            integral = R(integral)
            attempts=[]; witness=None
            # A degree drop means infinity itself is a rational preimage.
            if integral.degree()==degree:
                for p in primes:
                    Rp = PolynomialRing(GF(p),'z'); fp=Rp(integral)
                    values = [int(fp(i)) for i in range(p)]
                    projective_roots = sum(v==0 for v in values)+int(fp.degree()<degree)
                    attempts.append({'prime': p,'projective_roots': projective_roots})
                    if projective_roots==0:
                        witness={'prime':p,'primitive_polynomial_mod_p':list(map(int,fp.list())),
                                 'affine_values':values,'infinity_value':int(fp[degree])}
                        break
            data.append({'n':n,'label':case['label'],'parameter':str(tau),
                'status':'NO_RATIONAL_PREIMAGE' if witness else 'UNKNOWN',
                'certificate':witness,'attempts':attempts})
    return {'classification':'verified equation-only incidence comparison',
        'selection':'Unchanged nine addresses and two already constructed maps; unchanged12-prime certificate pool.',
        'cases':data, 'point_searches':0,'new_addresses':0,
        'inputs':{str(p.relative_to(ROOT)):sha(p) for p in
            [ROSTER, OUT/'protocol.json',OUT/'multiple2.json',OUT/'multiple3.json',Path(__file__)]},
        'boundary':'Absence of rational preimages excludes only these fixed covers. It is not a rank upper bound or an exclusion of seeds from other covers.'}


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true')
    args=parser.parse_args(); result=report()
    target=OUT/'controls.json'
    if args.write:
        with target.open('x') as stream:
            json.dump(result,stream,indent=2,sort_keys=True);stream.write('\n')
    else:
        assert result==json.loads(target.read_text())
    for row in result['cases']:
        print(row['n'],row['label'],row['status'],
              None if not row['certificate'] else row['certificate']['prime'])
