#!/usr/bin/env python3
"""Executable equations and independent point witnesses for twelve constructed fibres."""
import argparse
from pathlib import Path
import certify_compact_r17_candidates as cert
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';INPUT=ART/'native_rank3_carrier_subgroups_v1.json';OUT=ART/'new_native_carrier_rank19_curves.sage'
def expected():
    d=cert.read(INPUT);ledger=cert.read(ROOT/'artifacts/local/elliptic-curves/native-rank3-carrier-subgroups-v1/ledger.json')
    if d['status']!='PASS' or ledger['status']!='PASS' or len(d['rows'])!=12 or any(r['rank_lower_bound']!=19 for r in d['rows']):raise ArithmeticError('all twelve rank19 proofs and replay required')
    if any(cert.hashed(ROOT/n)!=h for n,h in d['sources'].items()):raise ArithmeticError('point certificate source differs')
    lines=['# Twelve constructed rank-at-least19 fibres; no exact-rank or record claim.','# Proof: native_rank3_carrier_subgroups_v1.json; models are displayed integral equations.','# No global minimality or universal novelty assertion.','from sage.all import QQ, EllipticCurve','curves = {}']
    for i,r in enumerate(d['rows']):
        name=f'native-carrier-19-{i+1:02}'
        lines += [f"E = EllipticCurve(QQ, {r['curve']!r})",f"points = [E([QQ(x), QQ(y)]) for x,y in {r['independent_points']!r}]",'assert len(points) == 19',f"curves[{name!r}] = (E, points)",f"print({name!r}, 'compact08234 parameter', {r['parameter']!r}, 'rank >= 19')"]
    return '\n'.join(lines)+'\n'
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();s=expected()
    if a.check:
        if OUT.read_text()!=s:raise ArithmeticError('Sage export differs')
    else:
        if OUT.exists():raise FileExistsError('preserve executable witness export')
        OUT.write_text(s)
    print('TWELVE RANK19 SAGE EXPORT PASS',flush=True)
